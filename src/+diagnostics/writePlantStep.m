function writePlantStep(registry, built, kt, dt_s)
%WRITEPLANTSTEP  Export accepted canopy, hydraulic, and root-uptake outputs.
%
%   All values come from the accepted outer-step transaction. Root-sink
%   totals use the same finite-element trapezoid as the Richards equation:
%
%       q_root = sum_e dz_e (Srt_e,left + Srt_e,right) / 2.   [cm s^-1]
%
%   Multiplication by dt_s and 10 converts the accepted rate to millimetres
%   per outer step. No diagnostic here feeds back into model physics.
%   Editor: Jizixing Zhu; Date: 2026-07-16.
    if ~isfield(built.workspace, 'scope_step')
        return;
    end
    scope = built.workspace.scope_step;
    writeCanopyState(registry, built, scope, kt);
    writeHydraulics(registry, built, scope, kt);
    writeRootUptake(registry, built, scope, kt, dt_s);
end

function writeCanopyState(registry, built, scope, kt)
% Export LAI, VPD, and accepted canopy-water state.
    if isfield(scope, 'canopy') && isfield(scope.canopy, 'LAI')
        registry.write('leaf_area_index', scope.canopy.LAI, kt);
    end
    if isfield(built.forcing, 'derived') ...
            && isfield(built.forcing.derived, 'es') ...
            && isfield(built.forcing.derived, 'ea') ...
            && numel(built.forcing.derived.es) >= kt ...
            && numel(built.forcing.derived.ea) >= kt
        vpd = max(built.forcing.derived.es(kt) ...
            - built.forcing.derived.ea(kt), 0) / 1000;
        registry.write('vapor_pressure_deficit', vpd, kt);
    end
    if ~isfield(built, 'diagnostics') ...
            || ~isfield(built.diagnostics, 'canopy_water')
        return;
    end
    canopy = built.diagnostics.canopy_water;
    mappings = { ...
        'h2ocan', 'h2ocan_history'; ...
        'h2ocan_sun', 'h2ocan_sun_history'; ...
        'h2ocan_sha', 'h2ocan_sha_history'; ...
        'fwet', 'fwet_history'; ...
        'fwet_sun', 'fwet_sun_history'; ...
        'fwet_sha', 'fwet_sha_history'; ...
        'lE_canopy_intercept', 'lE_CI_history'; ...
        'canopy_water_residual', 'canopy_water_residual'};
    for i = 1:size(mappings, 1)
        sourceName = mappings{i, 2};
        if isfield(canopy, sourceName) && size(canopy.(sourceName), 1) >= kt
            registry.write(mappings{i, 1}, canopy.(sourceName)(kt, 1), kt);
        end
    end
end

function writeHydraulics(registry, built, scope, kt)
% Export algebraic plant-water potentials and conductances.
    if ~isfield(scope, 'hydraulic') || ~isstruct(scope.hydraulic)
        return;
    end
    hydraulic = scope.hydraulic;
    mappings = { ...
        'leaf_potential', 'psil'; ...
        'root_potential', 'psir'; ...
        'soil_potential_weighted', 'psis'; ...
        'gsr_node_sum', 'gsr_total'; ...
        'k_plant_sun_diff', 'Kplant_diff_sun'; ...
        'k_plant_sha_diff', 'Kplant_diff_sha'; ...
        'k_plant_canopy_diff', 'Kplant_diff_flux_weighted'; ...
        'k_plant_sun_app', 'Kplant_app_sun'; ...
        'k_plant_sha_app', 'Kplant_app_sha'; ...
        'k_plant_canopy_app', 'Kplant_app_total'; ...
        'k_soil_leaf', 'Kwhole_app_total'; ...
        'k_soil_root', 'Ksoilroot_app_total'};
    writeScalarMappings(registry, hydraulic, mappings, kt);

    if isfield(built, 'diagnostics') ...
            && isfield(built.diagnostics, 'v2_hydraulic')
        dynamic = built.diagnostics.v2_hydraulic;
        dynamicMappings = { ...
            'w_root_v2', 'W_root'; 'w_stem_v2', 'W_stem'; ...
            'w_leaf_v2', 'W_leaf'; 'psi_leaf_v2', 'psi_leaf_v2'; ...
            'psi_root_v2', 'psi_root_v2'; 'psi_stem_v2', 'psi_stem_v2'; ...
            'rwc_stem_v2', 'RWC_stem'; 'rwc_leaf_v2', 'RWC_leaf'};
        writeScalarMappings(registry, dynamic, dynamicMappings, kt);
    end
end

function writeRootUptake(registry, built, scope, kt, dt_s)
% Export endpoint root sink, hydraulic redistribution, and column closure.
    if ~isfield(built.workspace, 'flux') ...
            || ~isfield(built.workspace.flux, 'Srt')
        return;
    end
    rootSink = built.workspace.flux.Srt;
    if ~isequal(size(rootSink), [built.domain.NL, 2])
        return;
    end
    registry.write('water_uptake_profile', sum(rootSink, 2), kt);

    hydraulic = struct();
    if isfield(scope, 'hydraulic') && isstruct(scope.hydraulic)
        hydraulic = scope.hydraulic;
    end
    if isfield(hydraulic, 'RWU_pos_node') ...
            && numel(hydraulic.RWU_pos_node) >= built.domain.NL
        registry.write('rwu_pos_profile', ...
            hydraulic.RWU_pos_node(1:built.domain.NL) * 100, kt);
    end
    if isfield(hydraulic, 'RWU_neg_node') ...
            && numel(hydraulic.RWU_neg_node) >= built.domain.NL
        registry.write('rwu_neg_profile', ...
            hydraulic.RWU_neg_node(1:built.domain.NL) * 100, kt);
    end
    if isfield(hydraulic, 'rwu_diag') && isstruct(hydraulic.rwu_diag)
        liftMappings = { ...
            'hydraulic_lift_intensity', 'qflx_phs_neg_m_s', -100 * 36000; ...
            'hydraulic_lift_pos', 'F_pos_m_s', 100 * 36000; ...
            'hydraulic_lift_net', 'F_net_m_s', 100 * 36000};
        for i = 1:size(liftMappings, 1)
            sourceName = liftMappings{i, 2};
            if isfield(hydraulic.rwu_diag, sourceName)
                registry.write(liftMappings{i, 1}, ...
                    hydraulic.rwu_diag.(sourceName) * liftMappings{i, 3}, kt);
            end
        end
    end

    dz = built.domain.dz(:);
    netRate = sum(0.5 * (rootSink(:, 1) + rootSink(:, 2)) .* dz);
    positiveRate = sum(0.5 * (max(rootSink(:, 1), 0) ...
        + max(rootSink(:, 2), 0)) .* dz);
    negativeRate = sum(0.5 * (min(rootSink(:, 1), 0) ...
        + min(rootSink(:, 2), 0)) .* dz);
    netMm = netRate * dt_s * 10;
    registry.write('root_sink_srt_net', netMm, kt);
    registry.write('root_sink_srt_pos', positiveRate * dt_s * 10, kt);
    registry.write('root_sink_srt_neg', negativeRate * dt_s * 10, kt);
    registry.write('root_sink_minus_transpiration', ...
        netMm - currentTranspirationMm(scope, dt_s), kt);
end

function value = currentTranspirationMm(scope, dt_s)
% Convert accepted canopy transpiration from centimetres per second to mm.
    value = 0;
    if isfield(scope, 'hydraulic') ...
            && isfield(scope.hydraulic, 'Tp_t')
        value = scope.hydraulic.Tp_t * dt_s * 10;
    end
end

function writeScalarMappings(registry, source, mappings, kt)
% Apply a two-column output/source mapping for finite scalar diagnostics.
    for i = 1:size(mappings, 1)
        sourceName = mappings{i, 2};
        if isfield(source, sourceName) ...
                && isscalar(source.(sourceName)) ...
                && isreal(source.(sourceName))
            registry.write(mappings{i, 1}, source.(sourceName), kt);
        end
    end
end
