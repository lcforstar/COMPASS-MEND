function response = estimateDynamicSurfaceAdmittance( ...
        built, kt, dt_s, heatRequest)
%ESTIMATEDYNAMICSURFACEADMITTANCE  Linearize the upper 20 cm heat block at z=0.
%
%   Inputs:
%     built       [struct] Accepted pre-soil model state.
%     kt          [-] One-based forcing index.
%     dt_s        [s] Atmospheric outer-step duration.
%     heatRequest [struct] Materialized heat request; its bottom policy is
%                  retained while its top is replaced by diagnostic fluxes.
%
%   Returns:
%     response.K_W_m2_K                  [W m-2 K-1] upper-soil one-step
%                                        thermal admittance.
%     response.zero_flux_temperature_C   [degC] predicted z=0 temperature
%                                        under zero top heat flux.
%     response.temperature_per_unit_flux_K_m2_W [K m2 W-1] reciprocal K.
%
%   Method:
%     The current water, ice, gas, and thermal coefficients are frozen at the
%     accepted state. The existing implicit heat matrix is assembled once.
%     A zero-flux solve supplies the affine intercept, and a +1 W m-2 top-RHS
%     perturbation supplies dT_top/dG. Its reciprocal is the Dirichlet-to-
%     Neumann admittance. The active block extends from z=0 to the first grid
%     node at or below 0.20 m; that node is held at its accepted temperature.
%     This fixed hourly approximation prevents a tens-of-metres soil column
%     from increasing the two predictor-solve costs. If the whole column is
%     shallower than 0.20 m, its configured bottom condition is retained.
%
%   The predictor freezes water-carried boundary heat at zero. The formal
%   coupled h/U/P_g solve remains authoritative and uses the matched fixed
%   heat flux produced from this response.
%
%   Reference:
%     soil.thermal.solveThermalTrial
%     surface.energy.matchDynamicSoilHeatResponse
% Editor: Jizixing Zhu; Date: 2026-07-20.

    validateInputs(built, kt, dt_s, heatRequest);
    predictor = built;
    predictor.workspace.step.trial.h = built.Y.h;
    predictor.workspace.step.trial.T = built.Y.T;
    predictor.workspace.step.trial.U = built.Y.U;
    if isfield(built.Y, 'p_gas')
        predictor.workspace.step.trial.p_gas = built.Y.p_gas;
    end
    refreshOptions = struct( ...
        'L_f', built.workspace.step.L_f, ...
        'dt_s', dt_s, 'KT', kt, 'NoTime_KT', kt, ...
        'DSTMAX', built.workspace.surface.DSTMAX, ...
        'refreshThermalBoundary', false);

    [acceptedView, acceptedEquations] = ...
        soil.refreshSoilCoefficients(predictor, refreshOptions);
    energyAccepted = acceptedEquations.heat.energyStorage;
    [predictor, equations] = ...
        soil.refreshSoilCoefficients(acceptedView, refreshOptions);
    rootSink = zeros(built.domain.NL, 2);
    if isfield(built.workspace, 'flux') ...
            && isfield(built.workspace.flux, 'Srt') ...
            && isequal(size(built.workspace.flux.Srt), size(rootSink))
        rootSink = built.workspace.flux.Srt;
    end
    matrices = soil.assembleEquationMatrices( ...
        equations, built.domain.dz, rootSink);

    request = heatRequest;
    request.top.kind = 'specified_heat_flux';
    request.top.requested_flux_into_domain_J_cm2_s = 0;
    request.top.value_source = 'dynamic_admittance_zero_flux_probe';
    waterFlux = struct('QMB', 0, 'Throughfall', 0);
    surfaceWaterFlux = struct('EVAP', 0);
    [~, ~, ~, system] = soil.thermal.solveThermalTrial( ...
        predictor, equations.heat, matrices.heat, energyAccepted, ...
        dt_s, kt, waterFlux, surfaceWaterFlux, true, request);

    matrix = system.matrix;
    [matrix, rhsZero, activeRange, responseDepth_m] = ...
        truncateAtInteractionDepth( ...
            matrix, system.rhs, built.Y.T(:), built.domain.dz, 0.20);
    temperatureZero = solveSystem( ...
        matrix, rhsZero, activeRange);
    rhsUnit = rhsZero;
    rhsUnit(end) = rhsUnit(end) + 1e-4; % 1 W m-2 = 1e-4 J cm-2 s-1
    temperatureUnit = solveSystem(matrix, rhsUnit, activeRange);
    temperaturePerFlux = temperatureUnit(end) - temperatureZero(end);
    assert(isfinite(temperaturePerFlux) && temperaturePerFlux > 0, ...
        'estimateDynamicSurfaceAdmittance:nonpositiveResponse', ...
        'A positive inward top heat flux must warm the top soil node.');

    response = struct( ...
        'schema_version', 1, ...
        'sign_convention', 'flux_into_soil_positive', ...
        'outer_step_s', dt_s, ...
        'K_W_m2_K', 1 / temperaturePerFlux, ...
        'zero_flux_temperature_C', temperatureZero(end), ...
        'temperature_per_unit_flux_K_m2_W', temperaturePerFlux, ...
        'interaction_depth_m', responseDepth_m, ...
        'active_node_count', activeRange(2) - activeRange(1) + 1, ...
        'column_depth_m', sum(built.domain.dz) / 100);
end

function temperature = solveSystem(matrix, rhs, activeRange)
% Solve one retained tridiagonal heat system.
% Editor: Jizixing Zhu; Date: 2026-07-20.
    temperature = solver.solveTridiagonal( ...
        matrix.C5_a, matrix.C5(:, 1), matrix.C5(:, 2), rhs, activeRange);
end

function [matrix, rhs, activeRange, actualDepth_m] = ...
        truncateAtInteractionDepth(matrix, rhs, acceptedTemperature_C, ...
        dz_cm, requestedDepth_m)
% Fix the lowest node of the upper-soil predictor at its accepted temperature.
% Editor: Jizixing Zhu; Date: 2026-07-20.
    nodeDepth_m = [0; cumsum(flip(dz_cm(:))) / 100];
    NN = numel(nodeDepth_m);
    nFromTop = find(nodeDepth_m >= requestedDepth_m, 1, 'first');
    if isempty(nFromTop)
        activeRange = [1, NN];
        actualDepth_m = nodeDepth_m(end);
        return;
    end
    first = NN - nFromTop + 1;
    activeRange = [first, NN];
    actualDepth_m = nodeDepth_m(nFromTop);
    if first == 1
        return;
    end

    fixedTemperature_C = acceptedTemperature_C(first);
    rhs(first) = fixedTemperature_C;
    matrix.C5(first, 1) = 1;
    matrix.C5(first, 2) = 0;
    rhs(first + 1) = rhs(first + 1) ...
        - matrix.C5_a(first) * fixedTemperature_C;
    matrix.C5_a(first) = 0;
end

function validateInputs(built, kt, dt_s, request)
% Reject incomplete states instead of manufacturing a thermal response.
% Editor: Jizixing Zhu; Date: 2026-07-20.
    assert(isstruct(built) && isfield(built, 'Y') ...
            && isfield(built.Y, 'T') && all(isfinite(built.Y.T(:))), ...
        'estimateDynamicSurfaceAdmittance:badState', ...
        'Accepted soil temperature must be finite.');
    assert(isscalar(kt) && isfinite(kt) && kt >= 1 ...
            && isscalar(dt_s) && isfinite(dt_s) && dt_s > 0, ...
        'estimateDynamicSurfaceAdmittance:badTime', ...
        'kt and dt_s must be positive finite scalars.');
    assert(isstruct(request) && isfield(request, 'enabled') ...
            && request.enabled && isfield(request, 'top') ...
            && isfield(request, 'bottom'), ...
        'estimateDynamicSurfaceAdmittance:badRequest', ...
        'An enabled materialized heat request is required.');
end
