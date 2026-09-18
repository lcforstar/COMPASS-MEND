function [model, report] = runWaterOnlyExperimentStep(model, kt, input)
%RUNWATERONLYEXPERIMENTSTEP  Advance one standalone Richards experiment step.
%
%   [model, report] = solver.runWaterOnlyExperimentStep(model, kt, input)
%   bypasses canopy, energy-balance, and plant production while reusing the
%   same advanceSoilEquationOuterStep scheduler and attemptSolveSoilEquations kernel
%   as the natural switched-soil chain. It is intended for water-only soil-column boundary
%   experiments and benchmarks, not as a second Richards implementation.
%
%   input is a transient plain-data struct with optional fields:
%     dt_s             [s] outer duration (default configured Delt_t)
%     t_start_s        [s] explicit start (default runtime.step.TIME or zero)
%     water_root_sink  [1/s, NL x 2] explicit internal sink; defaults to zero
%                      because the standalone component set has no plant
%     solve_ctl        solver controls; defaults to fixed HYDRUS convergence controls
%     solve_fn         transient test/diagnostic injection only
%
%   Laboratory head and flux values are read only from the site's configured
%   water-boundary CSV. The public flux sign is positive into the soil domain:
%   a negative top fixed flux is prescribed upward removal, such as a controlled
%   evaporation experiment. Atmospheric top boundaries belong to the natural
%   field driver and are rejected here instead of receiving fake zero forcing.
%
%   The final soil state and clock commit only when every tile is accepted.
%   report.attempts retains immutable request and actual boundary results for
%   both accepted and rejected trials; these audit records are not model state.
%
%   Editor: Jizixing Zhu; Date: 2026-07-11.
%   Change: 2026-07-11 Jizixing Zhu -- make the laboratory driver CSV-only, reject atmospheric top boundaries and the removed boundary_flux input, and expose the laboratory context in its report.
%   Change: 2026-07-12 Jizixing Zhu -- return canonical actual-flux water balance and
%           accepted h=0 interface diagnostics for laboratory columns.
%   Change: 2026-07-12 Jizixing Zhu -- preserve the configured freeze latent heat in
%           accepted-storage diagnostics for isothermal freeze-thaw experiments.
%   Change: 2026-07-12 Jizixing Zhu -- allow isothermal freeze-thaw Richards
%           experiments while continuing to reject heat, air, and nitrogen solves.
%   Change: 2026-07-17 Jizixing Zhu -- retain accepted finite-ponding
%           complementarity diagnostics in the standalone output workspace.

    if nargin < 3 || isempty(input), input = struct(); end
    assert(isstruct(input) && isscalar(input), ...
        'runWaterOnlyExperimentStep:badInput', 'input must be a scalar struct.');
    assertWaterOnly(model);
    assertLaboratoryBoundary(model.built);
    if isfield(input, 'boundary_flux')
        error('runWaterOnlyExperimentStep:boundaryFluxRemoved', ...
            ['Laboratory boundary values come only from water_boundary.csv. ' ...
             'Remove input.boundary_flux and select a prescribed top boundary mode.']);
    end

    built = model.built;
    built.workspace.step.L_f = experimentLatentHeat(built.control.switches);
    model.built = built;
    dt_s = readField(input, 'dt_s', built.control.time.Delt_t);
    tStart_s = readField(input, 't_start_s', readRuntimeTime(model.runtime));
    assert(isnumeric(dt_s) && isscalar(dt_s) && isfinite(dt_s) && dt_s > 0, ...
        'runWaterOnlyExperimentStep:badDt', 'input.dt_s must be finite and >0 [s].');
    assert(isnumeric(tStart_s) && isscalar(tStart_s) && isfinite(tStart_s) && tStart_s >= 0, ...
        'runWaterOnlyExperimentStep:badTime', 'input.t_start_s must be finite and >=0 [s].');

    currentFlux = struct('Throughfall', 0, 'EVAP', 0);   % [cm/s] inert legacy payload; prescribed requests override both fields
    solveCtl = readField(input, 'solve_ctl', defaultSolveControls(built));
    solveCtl.water_boundary_flux = currentFlux;
    solveCtl.water_root_sink = readField( ...
        input, 'water_root_sink', zeros(size(built.workspace.flux.Srt)));
    solveFn = readField(input, 'solve_fn', @solver.attemptSolveSoilEquations);
    storageBeforeMm = acceptedStorageMm(built, dt_s, kt);

    outerSnapshot = solver.snapshotSubstepState(model.built);
    model.built = resetTrialFromAccepted(model.built);
    waterAccum = solver.initWaterFluxAccum();
    model.built.workspace.flux.water_hourly = solver.finalizeWaterFluxAccum(waterAccum);
    model.built.workspace.flux.surface_boundary_hourly = struct('has_data', false);

    [model, advance] = solver.advanceSoilEquationOuterStep( ...
        model, kt, dt_s, built.control.time.Delt_t_min, ...
        built.control.numerics.substep, solveCtl, currentFlux, tStart_s, solveFn);

    if advance.accepted
        model.built = solver.carrySubstepState(model.built);
        model.built = surface.commitWaterCandidate( ...
            model.built, model.built.workspace.step.trial.surface_water);
        model.built.workspace.flux.water_hourly = ...
            solver.finalizeWaterFluxAccum(advance.water_flux_accum);
        if advance.surface_boundary_summary.has_data
            model.built.workspace.flux.surface_boundary_hourly = ...
                advance.surface_boundary_summary;
        end
        model.built.workspace.flux.Throughfall = currentFlux.Throughfall;
        model.built.workspace.flux.EVAP = currentFlux.EVAP;
        model.runtime.step.TIME = tStart_s + dt_s;
        storageAfterMm = acceptedStorageMm(model.built, dt_s, kt);
        rootSinkCmS = solver.integrateRootSinkCmS( ...
            solveCtl.water_root_sink, model.built.domain.dz);
        waterBalance = solver.computeWaterBalanceResidual( ...
            model.built.workspace.flux.water_hourly, storageBeforeMm, ...
            storageAfterMm, rootSinkCmS, dt_s);
        model.built.workspace.flux.water_balance = waterBalance;
    else
        model.built = solver.restoreSubstepState(model.built, outerSnapshot);
        waterBalance = emptyWaterBalance();
    end
    model.runtime.substep = advance.substep;

    zeroHead = soil.water.computeZeroHeadDiagnostics( ...
        model.built.Y.h, model.built.domain.z_node);
    report = struct( ...
        'context', 'laboratory_column', ...
        'kt', kt, ...
        'accepted', advance.accepted, ...
        't_start_s', tStart_s, ...
        't_end_s', tStart_s + advance.elapsed_s, ...
        'dt_s', dt_s, ...
        'elapsed_s', advance.elapsed_s, ...
        'residual', advance.residual, ...
        'iters', advance.iters, ...
        'water_flux', solver.finalizeWaterFluxAccum(advance.water_flux_accum), ...
        'water_balance', waterBalance, ...
        'zero_head', zeroHead, ...
        'attempts', advance.attempts);
end

% -------------------------------------------------------------------------
function mm = acceptedStorageMm(built, dt_s, kt)
% Clone-refresh accepted state with the exact trial time context before integration.
    probe = built;
    probe.workspace.step.trial.h = built.Y.h;
    probe.workspace.step.trial.T = built.Y.T;
    if isfield(built.Y, 'p_gas') && isfield(probe.workspace.step.trial, 'p_gas')
        probe.workspace.step.trial.p_gas = built.Y.p_gas;
    end
    probe = soil.refreshSoilCoefficients(probe, struct( ...
        'L_f', built.workspace.step.L_f, 'dt_s', dt_s, 'KT', kt));
    mm = solver.soilWaterEquivalentStorageMm( ...
        probe.workspace.derived.Theta_UU, probe.workspace.derived.Theta_V, ...
        probe.workspace.derived.RHOV, built.domain.dz);
end

% -------------------------------------------------------------------------
function wb = emptyWaterBalance()
% Return missing diagnostics for a fully rejected outer experiment step.
% Editor: Jizixing Zhu; Date: 2026-07-17.
    wb = solver.emptyWaterBalance();
end

% -------------------------------------------------------------------------
function assertLaboratoryBoundary(built)
% Reject natural-surface policy in a driver that has no atmosphere or canopy.
    assert(isfield(built, 'boundary') && isfield(built.boundary, 'water') ...
        && isfield(built.boundary.water, 'top'), ...
        'runWaterOnlyExperimentStep:missingWaterBoundary', ...
        'built.boundary.water.top is required for a laboratory experiment.');
    if strcmp(built.boundary.water.top.kind, 'atmospheric')
        error('runWaterOnlyExperimentStep:atmosphericTop', ...
            ['The laboratory driver cannot use an atmospheric top boundary. ' ...
             'Select constant_head, variable_head, constant_flux, variable_flux, or zero_flux.']);
    end
end

% -------------------------------------------------------------------------
function built = resetTrialFromAccepted(built)
% Start an outer experiment step from accepted prognostic solve fields.
    fields = {'h','T','p_gas','h_frez'};
    for i = 1:numel(fields)
        name = fields{i};
        if isfield(built.Y, name) && isfield(built.workspace.step.trial, name)
            built.workspace.step.trial.(name) = built.Y.(name);
        end
    end
end

% -------------------------------------------------------------------------
function solveCtl = defaultSolveControls(built)
% Minimal controls consumed by the unified solver in its water-only mode.
% Editor: Jizixing Zhu; Date: 2026-07-12.
% Change: 2026-07-12 Jizixing Zhu -- remove the obsolete water-form selector;
% mixed storage is now fixed behavior inside the Richards solver.
% Change: 2026-07-12 Jizixing Zhu -- remove the obsolete raw-head and
% water-convergence-policy controls; HYDRUS convergence is fixed behavior.
    nm = built.control.numerics;
    solveCtl = struct('NIT', nm.NIT, ...
        'tol_theta', readField(nm, 'tol_theta', 1e-4), ...
        'tol_head_saturated_cm', readField( ...
            nm, 'tol_head_saturated_cm', 0.1), ...
        'tol_water_mass_mm', readField(nm, 'tol_water_mass_mm', 1e-4));
end

% -------------------------------------------------------------------------
function assertWaterOnly(model)
% Fail instead of silently running an experiment under an unsupported switch set.
    assert(isstruct(model) && isfield(model, 'built') && isfield(model, 'runtime'), ...
        'runWaterOnlyExperimentStep:badModel', ...
        'model must be initialized and runtime-prepared.');
    sw = model.built.control.switches;
    nm = model.built.control.numerics;
    ok = readField(sw, 'lThmrlefc', 1) == 0 ...
        && readField(sw, 'lSoilairefc', 1) == 0 ...
        && readField(sw, 'lNefc', 1) == 0 ...
        && strcmpi(readField(nm, 'MatricKernel', 'array'), 'scalar');
    assert(ok, 'runWaterOnlyExperimentStep:notWaterOnly', ...
        'Standalone experiment driver requires the scalar water-only switch set.');
end

% -------------------------------------------------------------------------
function L_f = experimentLatentHeat(switches)
% Match the natural driver while holding temperature prescribed in this driver.
    if readField(switches, 'lFreezing', 0) == 1
        L_f = equations.getPhysicalConstant('lambda_fus');
    else
        L_f = 0;
    end
    override = str2double(io.getEnv('LFUSION_OVERRIDE'));
    if isfinite(override), L_f = override; end
end

% -------------------------------------------------------------------------
function t_s = readRuntimeTime(runtime)
% Read committed elapsed time without advancing it.
    t_s = 0;
    if isfield(runtime, 'step') && isfield(runtime.step, 'TIME') ...
            && ~isempty(runtime.step.TIME)
        t_s = runtime.step.TIME;
    end
end

% -------------------------------------------------------------------------
function value = readField(s, name, dflt)
% Read one optional field without adding stateful defaults to the model.
    value = dflt;
    if isfield(s, name) && ~isempty(s.(name)), value = s.(name); end
end
