function [model, candidates] = stagePostEbal( ...
        model, kt, dt_outer_s, lai_state_next, snowTransactionOverride)
%STAGEPOSTEBAL  Finalize the surface-energy pass and downstream candidates.
%
%   Inputs:
%     model          [inout] Model shell after stageCanopyEnergy.
%     kt             [-] One-based outer-step index.
%     dt_outer_s     [s] Actual outer-step duration.
%     lai_state_next [in] Pending LAI state from stageCanopyInputs.
%
%   Returns:
%     model          [out] Model with Srt and diagnostic radiation staged.
%     candidates     [struct] Pending hydraulic, LAI-stress, LUNA, surface,
%                    canopy-water, snow, and water-flux values.
%
%   The order is load-bearing: guard the provisional SCOPE surface; advance
%   canopy water and snow from accepted state; optionally replace the SCOPE
%   snow surface and repeat RTMo/ebal once; only then run plant hydraulics,
%   water-stress accumulation, root uptake, LUNA, and SIF.  Thus all downstream
%   consumers read the final ebal pass and snow storage is never double-stepped.
%   Editor: Jizixing Zhu; Date: 2026-07-16.
%   Change: 2026-07-17 Jizixing Zhu -- stage the ClimaLand snow transaction
%     and its canonical diagnostics after the guarded surface-energy solve.
%   Change: 2026-07-17 Jizixing Zhu -- finalize the snow/SCOPE interface
%     before hydraulic and fluorescence consumers, with at most one bounded
%     correction when snow cover appears, disappears, or changes abruptly.
%   Change: 2026-07-17 Jizixing Zhu -- replace the obsolete pre-snow consumer
%     order in the interface documentation with the implemented final order.
%   Change: 2026-07-18 Jizixing Zhu -- expose the finalized SCOPE force-restore
%     history as an accept-only outer-step candidate.
%   Change: 2026-07-19 Jizixing Zhu -- reuse a reconciled terminal-meltout snow
%     transaction during the bounded post-soil surface event restage.

    if nargin < 5
        snowTransactionOverride = struct();
    end

    % Finalize the surface transaction first so every downstream consumer
    % reads one corresponding ebal/snow state.  Each pass restarts snow from
    % accepted Y.snow; the predictor/correction never advances storage twice.
    scope_step = model.built.workspace.scope_step;
    scope_step = vegetation.energy.guardSurfaceTemperature( ...
        scope_step, model.built.Y.surface.Tss, kt);
    model.built.workspace.scope_step = scope_step;
    [model, water_snow] = vegetation.integrated.advanceCanopyWaterSnow( ...
        model, dt_outer_s, snowTransactionOverride);
    if isempty(fieldnames(snowTransactionOverride))
        [correctSurface, correctionDiagnostic] = ...
            needsSnowSurfaceCorrection(model.built, water_snow.snow);
    else
        correctSurface = false;
        correctionDiagnostic = struct( ...
            'requested', false, ...
            'reason', 'pre_reconciled_terminal_event');
    end
    if correctSurface
        model = vegetation.integrated.stageCanopyEnergy( ...
            model, kt, dt_outer_s, water_snow.snow, false);
        scope_step = model.built.workspace.scope_step;
        scope_step = vegetation.energy.guardSurfaceTemperature( ...
            scope_step, model.built.Y.surface.Tss, kt);
        model.built.workspace.scope_step = scope_step;
        [model, water_snow] = vegetation.integrated.advanceCanopyWaterSnow( ...
            model, dt_outer_s);
    end
    model.built.workspace.snow.surface_correction_applied = correctSurface;
    model.built.workspace.snow.surface_correction = correctionDiagnostic;

    % Plant hydraulics, root uptake, LUNA, and SIF now consume only the final
    % ebal pass rather than a provisional pre-correction surface.
    [plant_water, hydraulic_diagnostics] = ...
        vegetation.hydraulics.stepDynamicStorage(model, dt_outer_s);
    lai_state_next = ...
        vegetation.integrated.phenology.accumulateLaiWaterStress( ...
        model, lai_state_next);
    model = vegetation.root.stageWaterUptake(model);
    luna_climate_entry = ...
        vegetation.integrated.LUNA.buildClimateEntry(model, kt);
    model = surface.radiation.runStepSIF(model);

    scope_step = model.built.workspace.scope_step;
    acceptedSurfaceState = model.built.Y.surface;
    acceptedSurfaceState.Ts = scope_step.thermal.Ts;
    acceptedSurfaceState.Tss = scope_step.thermal.Tss;
    scopeSoil = model.built.Y.scope_soil;
    if isfield(scope_step, 'Soil') && isfield(scope_step.Soil, 'Tsold') ...
            && ~isempty(scope_step.Soil.Tsold)
        scopeSoil.Tsold = scope_step.Soil.Tsold;
    end
    candidates = struct( ...
        'plant_hydraulic_W', plant_water, ...
        'hydraulic_diagnostics', hydraulic_diagnostics, ...
        'lai_state', lai_state_next, ...
        'luna_climate_entry', luna_climate_entry, ...
        'surface', acceptedSurfaceState, ...
        'scope_soil', scopeSoil, ...
        'canopy_water', water_snow.canopy_water, ...
        'snow', water_snow.snow, ...
        'water_flux', water_snow.water_flux, ...
        'canopy_water_diagnostics', water_snow.history_row, ...
        'snow_diagnostics', water_snow.snow_diagnostics, ...
        'Srt', model.built.workspace.scope_step.Srt, ...
        'surface_soil_water', water_snow.surface_soil_water, ...
        'root_soil_exchange', ...
        model.built.workspace.scope_step.root_soil_exchange);
end

function [needsCorrection, diagnostic] = ...
        needsSnowSurfaceCorrection(built, candidateState)
% Request one correction only for a snow-regime flip or large cover change.
% Editor: Jizixing Zhu; Date: 2026-07-17.
    surfaceState = built.Y.snow;
    if isfield(built.workspace, 'snow') ...
            && isfield(built.workspace.snow, 'surface_state')
        surfaceState = built.workspace.snow.surface_state;
    end
    diagnosedSurface = surface.snow.diagnoseState( ...
        surfaceState, built.parameters.snow);
    candidate = surface.snow.diagnoseState( ...
        candidateState, built.parameters.snow);
    tolerance = readNumerics( ...
        built.control.numerics, 'tol_energy_balance_J_cm2', 1e-2);
    surfaceActive = diagnosedSurface.cover_fraction > 0 ...
        && abs(diagnosedSurface.energy_J_cm2) > tolerance;
    candidateActive = candidate.cover_fraction > 0 ...
        && abs(candidate.energy_J_cm2) > tolerance;
    regimeFlip = xor(surfaceActive, candidateActive);
    coverChange = abs(candidate.cover_fraction ...
        - diagnosedSurface.cover_fraction);
    % Current precipitation is already present in surfaceState.  The second
    % solve is reserved for meltout/freeze onset or a material geometry jump,
    % avoiding a costly outer ebal loop during ordinary snow-covered hours.
    needsCorrection = regimeFlip || coverChange > 0.15;
    diagnostic = struct( ...
        'surface_cover_fraction', diagnosedSurface.cover_fraction, ...
        'candidate_cover_fraction', candidate.cover_fraction, ...
        'cover_fraction_change', coverChange, ...
        'surface_thermal_active', surfaceActive, ...
        'candidate_thermal_active', candidateActive, ...
        'regime_flip', regimeFlip, ...
        'requested', needsCorrection);
end

function value = readNumerics(numerics, name, fallback)
% Read one optional numerical control without mutating configuration.
% Editor: Jizixing Zhu; Date: 2026-07-17.
    value = fallback;
    if isfield(numerics, name) && ~isempty(numerics.(name))
        value = numerics.(name);
    end
end
