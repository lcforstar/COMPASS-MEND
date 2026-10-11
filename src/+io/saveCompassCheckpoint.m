function file = saveCompassCheckpoint(model, kt)
%SAVECOMPASSCHECKPOINT  Write a resume checkpoint (built.Y + runtime prognostic memory).
%
%   file = io.saveCompassCheckpoint(model, kt) serializes the minimal state needed
%   to continue a run at kt+1 into runtime.output.checkpoint_dir. What is stored is
%   deliberately scoped:
%     - built.Y                          prognostic soil/plant state (the core); this now includes the whole
%                                          Y.plant_acclimation (LUNA state/last_* + leaf_resp + the 10-day
%                                          climate buffer as a serialized struct), increment 6a + 6b
%     - runtime.step/substep/calendar    time / adaptive-step / wall-clock counters
%     - runtime.groundwater              as-is (no handles)
%   NOT stored (rebuilt after load): workspace/aux (re-derived from Y via
%   init.rebuildSoilAuxiliaryFromAccepted on resume; never initializeSoilAux,
%   which would mutate Y.h / rewrite Y.U) and the diagnostics Registry
%   (re-filled from the resume KT).
%
%   Inputs:
%     model [in] struct  a runCOMPASS model after prepareRuntime + beginRun.
%     kt    [in] [-]      the just-completed outer step this checkpoint restarts from.
%
%   Returns:
%     file  [char]  full path of the written .mat (contains a single `ckpt` struct).
%
%   Editor: Jizixing Zhu; Date: 2026-07-01.
%   Change: 2026-09-29 Wenjie Liu -- retain actual DIN observation parameters
%     so sample scoring cannot silently use baseline adsorption capacity.
%   Change: 2026-07-05 Jizixing Zhu -- stop serializing runtime.snow; SWE now lives in Y.snow (prognostic) and rides ckpt.Y whole-store (state-ownership migration).
%   Change: 2026-07-05 Jizixing Zhu -- doc: Y.surface.Tss now rides ckpt.Y (increment 2); previously runtime.thermal was never serialized (Tss lost on resume) -- gap now fixed.
%   Change: 2026-07-08 Jizixing Zhu -- doc: Y.surface.Ts also rides ckpt.Y as the accepted SCOPE/ebal thermal.Ts carry.
%   Change: 2026-07-05 Jizixing Zhu -- delete stateOnlyCanopyWater + ckpt.runtime.canopy_water; canopy_water state now lives in Y.canopy_water and rides ckpt.Y whole-store (increment 3 stage 2).
%   Change: 2026-07-05 Jizixing Zhu -- stop serializing runtime.phenology/lai_state; now Y.vegetation rides ckpt.Y whole-store (increment 5).
%   Change: 2026-07-05 Jizixing Zhu -- delete stateOnlyLeafResp + ckpt.runtime.leaf_resp; reduce snapshotLuna to climate_buffer only. LUNA state + last_* + leaf_resp.Ta_history now ride ckpt.Y (Y.plant_acclimation); increment 6a.
%   Change: 2026-07-06 Jizixing Zhu -- delete snapshotLuna + ckpt.runtime.luna; the climate buffer now rides ckpt.Y (Y.plant_acclimation.luna_climate_buffer); runtime.luna gone; increment 6b.
%   Change: 2026-07-19 Jizixing Zhu -- doc: resume rebuild uses rebuildSoilAuxiliaryFromAccepted, not initializeSoilAux.
%   Change: 2026-08-21 Jizixing Zhu -- checkpoint schema 2 records a strict
%     DOM/MEMC/litter compatibility signature.

    if ~isstruct(model) || ~isfield(model, 'runtime') || ~isfield(model.runtime, 'output')
        error('saveCompassCheckpoint:badModel', 'model.runtime.output is required.');
    end
    dir = model.runtime.output.checkpoint_dir;
    if ~isfolder(dir)
        error('saveCompassCheckpoint:noDir', ...
            'checkpoint_dir does not exist: %s (run beginRun first).', dir);
    end

    ckpt = buildCheckpoint(model, kt);
    file = fullfile(dir, sprintf('checkpoint_KT%05d.mat', kt));
    save(file, 'ckpt', '-v7.3');
end

% ---------------------------------------------------------------------------
function ckpt = buildCheckpoint(model, kt)
% Assemble checkpoint schema 2 with optional-biogeochemistry compatibility.
% Editor: Jizixing Zhu; Date: 2026-09-15.
% Change: 2026-09-15 Jizixing Zhu -- retain accepted root-water history used
% by the next outer-step predictor (bottom-up layer endpoints, s^-1).
% Change: 2026-09-30 Wenjie Liu -- retain the two accepted surface fluxes
% consumed by the next abrupt-wetting predictor; they are causal memory,
% not reconstructible transport coefficients or diagnostic-only records.
    rt = model.runtime;
    ckpt = struct();
    ckpt.schema_version = 2;
    ckpt.kt = kt;
    ckpt.created_at = datetime('now');

    ckpt.Y = model.built.Y;
    ckpt.compatibility = io.checkpointCompatibilitySignature(model.built);
    if isfield(model.built.parameters,'soil_n')
        p=model.built.parameters.soil_n;
        ckpt.din_observation_contract=struct('schema_version',1, ...
            'dz_cm',model.built.domain.dz(:), ...
            'RHO_bulk',p.RHO_bulk(:),'RHOKG',p.RHOKG(:), ...
            'observable','dissolved_plus_equilibrium_adsorbed');
    end
    rootWater = model.built.workspace.flux.Srt;
    assert(all(isfinite(rootWater(:))), ...
        'saveCompassCheckpoint:rootWaterHistory', 'Root-water history must be finite.');
    ckpt.coupling_memory = struct('completed_kt',kt, ...
        'root_water_element_endpoint_rate_1_s',rootWater);
    flux=model.built.workspace.flux;
    for name={'Throughfall','EVAP'}
        value=flux.(name{1});
        assert(isnumeric(value)&&isscalar(value)&&isreal(value)&&isfinite(value), ...
            'saveCompassCheckpoint:surfaceWaterHistory', ...
            'Accepted surface-water history %s must be a finite real scalar [cm/s].',name{1});
    end
    ckpt.coupling_memory.surface_water_history=struct( ...
        'throughfall_cm_s',flux.Throughfall,'evap_cm_s',flux.EVAP);

    ckpt.runtime = struct();
    ckpt.runtime.step        = rt.step;
    ckpt.runtime.substep     = rt.substep;
    ckpt.runtime.calendar    = rt.calendar;
    % canopy_water state now lives in Y.canopy_water (prognostic) -> rides ckpt.Y whole-store; histories
    % are built.diagnostics (not checkpointed, rebuilt on resume); no runtime.canopy_water serialize (2026-07-05).
    % snow.SWE now lives in Y.snow (prognostic) -> rides ckpt.Y whole-store; no runtime.snow serialize (2026-07-05)
    % surface.Ts/Tss likewise live in Y.surface -> ride ckpt.Y. runtime.thermal was NEVER serialized here
    % (thermal fields lost on resume, re-init NaN); the Y.surface owner fixes that gap.
    % phenology + lai_state now live in Y.vegetation (prognostic) -> ride ckpt.Y whole-store; no runtime
    % serialize (2026-07-05 increment 5).
    % LUNA state (+ last_*), leaf_resp.Ta_history, AND the LUNA climate buffer (serialized struct) now all live
    % in Y.plant_acclimation -> ride ckpt.Y whole-store; runtime owns no LUNA field (2026-07-06 increment 6b).
    ckpt.runtime.groundwater = rt.groundwater;
    ckpt.runtime.management  = rt.management;                % fertilizer once-only state (serializable)
    % Change: 2026-09-29 Wenjie Liu -- persist optional replay seam memory
    % alongside its authoritative accepted runtime.step clock.
    if isfield(rt,'biogeochemistry_cycle')
        assert(rt.step.KT==kt,'saveCompassCheckpoint:cycleClock', ...
            'Replay context must be saved at its actual accepted KT.');
        init.validateBiogeochemistryCycleContext(rt.biogeochemistry_cycle,model.built,kt);
        ckpt.runtime.biogeochemistry_cycle=rt.biogeochemistry_cycle;
    end
    if isfield(rt,'biogeochemistry_spinup_identity')
        ckpt.runtime.biogeochemistry_spinup_identity=rt.biogeochemistry_spinup_identity;
        ckpt.runtime.biogeochemistry_spinup_cycle_steps=rt.biogeochemistry_spinup_cycle_steps;
    end
end
