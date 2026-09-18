function ctx = prepareResidualContext(built, dt_s, kt)
%PREPARERESIDUALCONTEXT Freeze accepted anchors for one exact evaluation.
%
%   The returned context contains only accepted-time storage anchors, the
%   materialized boundary requests, forcing-side water data, and numerical
%   scales.  It performs no nonlinear solve and never mutates BUILT.
%
%   Editor: Jizixing Zhu; Date: 2026-08-14.
%   Change: 2026-08-15 Jizixing Zhu -- expose the accepted internal-energy
%           capacity so the main Picard driver reuses this single accepted-
%           state refresh instead of building a second anchor bank.

    assert(isstruct(built) && isscalar(dt_s) && isfinite(dt_s) ...
            && dt_s > 0 && isscalar(kt) && isfinite(kt), ...
        'prepareResidualContext:badInput', ...
        'built, a positive dt_s, and a finite kt are required.');

    airEnabled = enabled(built, 'lSoilairefc');
    heatEnabled = enabled(built, 'lThmrlefc');
    refreshOpts = struct('L_f', localField(built.workspace.step, 'L_f', 0), ...
        'dt_s', dt_s, 'KT', kt, 'NoTime_KT', kt, ...
        'DSTMAX', built.workspace.surface.DSTMAX, ...
        'refreshThermalBoundary', true);

    accepted = built;
    accepted.workspace.step.trial.h = built.Y.h;
    accepted.workspace.step.trial.T = built.Y.T;
    if isfield(built.Y, 'U')
        accepted.workspace.step.trial.U = built.Y.U;
    end
    if isfield(built.Y, 'p_gas')
        accepted.workspace.step.trial.p_gas = built.Y.p_gas;
    end
    [accepted, acceptedEquations] = soil.refreshSoilCoefficients( ...
        accepted, refreshOpts);
    thetaEqN = waterEquivalentEndpoints(accepted.workspace.derived);

    ctx = struct();
    ctx.schema_version = 1;
    ctx.dt_s = dt_s;
    ctx.kt = kt;
    ctx.refresh_opts = refreshOpts;
    ctx.air_enabled = airEnabled;
    ctx.heat_enabled = heatEnabled;
    ctx.theta_equivalent_accepted = thetaEqN;
    ctx.dry_air_storage_accepted = [];
    ctx.energy_storage_accepted = [];
    ctx.energy_capacity_accepted = [];
    if airEnabled
        ctx.dry_air_storage_accepted = ...
            acceptedEquations.dry_air.dryAirStorage;
        acceptedAir = soil.gas.evaluateDryAirAdmissibility( ...
            accepted.workspace.derived.RHODA, ...
            ctx.dry_air_storage_accepted);
        assert(isequal(size(ctx.dry_air_storage_accepted), ...
                size(thetaEqN)) && acceptedAir.admissible, ...
            'prepareResidualContext:badAcceptedDryAir', ...
            ['Accepted dry-air density/storage must be finite, nonnegative, ', ...
             'and storage must match the endpoint water array.']);
    end
    if heatEnabled
        ctx.energy_storage_accepted = ...
            acceptedEquations.heat.energyStorage;
        assert(isequal(size(ctx.energy_storage_accepted), size(thetaEqN)) ...
                && all(isfinite(ctx.energy_storage_accepted(:))), ...
            'prepareResidualContext:badAcceptedEnergy', ...
            'Accepted internal energy must be a finite endpoint array.');
        [~, ctx.energy_capacity_accepted] = ...
            soil.thermal.assembleNodalInternalEnergy( ...
                acceptedEquations.heat.energyStorage, ...
                acceptedEquations.heat.CTh, ...
                acceptedEquations.heat.CTT, ...
                acceptedEquations.heat.CTa, built.domain.dz);
    end
    ctx.water_root_sink = built.workspace.flux.Srt;
    ctx.water_boundary_flux = built.workspace.flux;
    ctx.water_boundary_request = [];
    ctx.heat_boundary_request = [];
    ctx.air_boundary_request = [];
    if isfield(built.workspace, 'boundary') ...
            && isstruct(built.workspace.boundary)
        ctx = readRequests(ctx, built.workspace.boundary);
    end
    if isfield(built.workspace.step, 'boundary_request') ...
            && isstruct(built.workspace.step.boundary_request)
        ctx = readRequests(ctx, built.workspace.step.boundary_request);
    end
    ctx.active_richards_first_node = ...
        solver.activeRichardsFirstNode(built);
    ctx.scales = struct('water_cm_s', 1e-6, ...
        'dry_air_g_cm2_s', 1e-12, 'energy_J_cm2_s', 1e-2);
    ctx.scale_id = 'compass.c2.residual_scales.v1';
end

function thetaEq = waterEquivalentEndpoints(derived)
    NL = size(derived.Theta_UU, 1);
    rhoVapor = derived.RHOV(:);
    vaporEndpoints = [rhoVapor(1:NL), rhoVapor(2:NL + 1)];
    rhoLiquid = equations.getPhysicalConstant('RHOL');
    thetaEq = derived.Theta_UU + derived.Theta_V ...
        .* vaporEndpoints / rhoLiquid;
end

function ctx = readRequests(ctx, source)
    if isfield(source, 'water')
        ctx.water_boundary_request = source.water;
    end
    if isfield(source, 'heat')
        ctx.heat_boundary_request = source.heat;
    end
    if isfield(source, 'air')
        ctx.air_boundary_request = source.air;
    end
end

function tf = enabled(built, name)
    tf = isfield(built.control.switches, name) ...
        && built.control.switches.(name) == 1;
end

function value = localField(container, name, fallback)
    value = fallback;
    if isfield(container, name) && ~isempty(container.(name))
        value = container.(name);
    end
end
