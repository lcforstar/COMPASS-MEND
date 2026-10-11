function model = beginRun(model)
%BEGINRUN  Pre-loop one-time runtime actions (layer 3 of runCOMPASS).
%
%   model = beginRun(model) performs the side-effecting actions that must happen
%   ONCE after init.prepareRuntime but before the KT loop, and that are NOT pure
%   initialization: restart-state loading, output-directory creation, and
%   groundwater-coupling pre-setup.
%
%   Deliberately kept OUT of init.prepareRuntime: prepareRuntime is a pure
%   allocator (no prognostic-state writes, no file I/O, repeatable). The actions
%   here write prognostic state and the filesystem, so they live in their own
%   layer that runs exactly once before stepping.
%
%   Run flow (runCOMPASS):
%       1) init.initialize  2) init.prepareRuntime  3) beginRun
%       4) for kt: runOneStep  5) finalizeRun
%
%   Inputs/Returns:
%     model [inout] struct  from init.prepareRuntime. May be updated: resume
%                           overwrites built.Y/workspace and KT/TIME; output
%                           directories are created on disk.
%
%   Status: all one-time actions implemented -- 3.1 resume load/apply,
%   3.2 output-dir creation, 3.3 groundwater setup + Dunnian forcing adapter. See
%   docs/init_refactor_progress.md (五层执行流).
%
%   Editor: Jizixing Zhu; Date: 2026-06-29.
%   Change: 2026-06-29 Jizixing Zhu — implement 3.2 output-directory creation (ensureOutputDirs); 3.1/3.3 still pending.
%   Change: 2026-07-01 Jizixing Zhu — validate runtime.output before the one-time actions.
%   Change: 2026-07-01 Jizixing Zhu — implement 3.1 resume load/apply (return-value, no assignin/evalin); 3.3 groundwater still pending.
%   Change: 2026-07-01 Jizixing Zhu — implement 3.3 groundwater main (gw_Dep/index/tempBotm/NBCh_bottom); Dunnian precip deferred to a forcing adapter.
%   Change: 2026-07-01 Jizixing Zhu — wire the Dunnian forcing adapter (init.applyGroundwaterForcingAdjustment) after applyGroundwaterSetup; 3.3 now complete.
%   Change: 2026-07-20 Jizixing Zhu -- materialize missing soil-N solver history after a legacy checkpoint resume.
%   Change: 2026-07-20 Jizixing Zhu -- seed explicitly absent legacy PlantCN pools from
%     the current initial-state schema and mark the non-parity compatibility path.
%   Change: 2026-08-21 Jizixing Zhu -- validate biogeochemistry checkpoint
%     compatibility and guard blank DOM/DON production cold starts.
%   Change: 2026-09-12 Jizixing Zhu -- allow finite nonnegative blank DOC/DON
%     production cold starts and record the transient-start provenance.

    validateBeginRun(model);

    % 3.1 resume: if requested, load the checkpoint ONCE here (not per step) and
    %   rebuild the live state from it (Y + workspace/aux + runtime memory + resume
    %   point). Return-value load — replaces legacy assignin/evalin restore.
    if isfield(model.runtime, 'resume') && model.runtime.resume.requested
        loaded = io.loadCompassCheckpoint(model.runtime.resume.file);
        io.validateCheckpointCompatibility(model.built,loaded);
        model = applyResumeCheckpoint(model, loaded);
    end
    model = validateDomColdStart(model);

    % 3.2 output directories: prepareRuntime only RECORDED the paths (pure
    %   allocator, no disk). Create them once here so finalizeRun / checkpoint
    %   writers can assume the directories exist.
    ensureOutputDirs(model.runtime.output);

    % 3.3 groundwater pre-loop setup -- only when coupling is on (changwu = 0 -> no-op).
    %   applyGroundwaterSetup sets water-table depth, bottom-layer index, bottom
    %   temperature and the Dirichlet bottom water BC. The Dunnian forcing adapter
    %   (init.applyGroundwaterForcingAdjustment) then applies saturation-excess
    %   runoff to daily raw.Precip and REBUILDS hourly+derived+gas so the change
    %   reaches the loop -- reproducing legacy COMPASS.m:325-326, which ran before
    %   the in-loop daily->hourly downscale. It runs AFTER applyGroundwaterSetup
    %   because it consumes the gw_Dep that setup computes (not a per-step runOneStep
    %   correction). See docs/init_refactor_progress.md.
    if isfield(model.runtime, 'groundwater') && model.runtime.groundwater.GroundwaterCoupling == 1
        model = applyGroundwaterSetup(model);
        model.built.forcing = init.applyGroundwaterForcingAdjustment( ...
            model.built.forcing, model.runtime.groundwater.gw_Dep, model.built);
    end
end

function model = validateDomColdStart(model)
% Blank DOC/DON profiles are legal cold-start states with transient provenance.
% Editor: Jizixing Zhu; Date: 2026-09-12.
    sw=model.built.control.switches;
    if ~isfield(sw,'lSoilDOM')||sw.lSoilDOM~=1, return; end
    resumed=isfield(model.runtime,'resume')&&model.runtime.resume.requested;
    spinup=isfield(model.runtime,'biogeochemistry_spinup') ...
        &&logical(model.runtime.biogeochemistry_spinup);
    blank=isfield(model.built.Y,'soil_dom') ...
        &&isfield(model.built.Y.soil_dom,'initialized_from_blank') ...
        &&model.built.Y.soil_dom.initialized_from_blank;
    if blank&&~resumed&&~spinup
        dom=model.built.Y.soil_dom;
        values=[dom.doc_total_g_cm3(:);dom.don_total_g_cm3(:)];
        assert(all(isfinite(values))&&all(values>=0), ...
            'beginRun:invalidBlankDom', ...
            'Blank-start DOC/DON stores must be finite and nonnegative.');
        model.runtime.biogeochemistry_cold_start=struct( ...
            'enabled',true,'doc_blank',logical(dom.initial_doc_blank), ...
            'don_blank',logical(dom.initial_don_blank), ...
            'interpretation','valid_zero_inventory_with_initial_transient');
        warning('beginRun:blankDomColdStart', ...
            ['Starting production with zero DOC/DON inventory. Early output ', ...
             'contains a biogeochemical initialization transient.']);
    else
        model.runtime.biogeochemistry_cold_start=struct( ...
            'enabled',false,'doc_blank',false,'don_blank',false, ...
            'interpretation','resumed_or_initialized_inventory');
    end
end

% ---------------------------------------------------------------------------
function validateBeginRun(model)
% Guard that init.prepareRuntime ran and produced the output-path record beginRun
% consumes (runtime.output with the four writer directories).
% Editor: Jizixing Zhu; Date: 2026-07-01.
    if ~isstruct(model) || ~isfield(model, 'runtime') || ~isstruct(model.runtime)
        error('beginRun:missingRuntime', ...
            'beginRun requires init.prepareRuntime(model) before the run loop.');
    end
    if ~isfield(model.runtime, 'output') || ~isstruct(model.runtime.output)
        error('beginRun:missingOutput', ...
            'model.runtime.output is required (produced by init.prepareRuntime).');
    end
    req = {'output_dir', 'checkpoint_dir', 'nc_dir', 'csv_dir'};
    for i = 1:numel(req)
        if ~isfield(model.runtime.output, req{i})
            error('beginRun:badOutput', 'model.runtime.output.%s is missing.', req{i});
        end
    end
end

% ---------------------------------------------------------------------------
function model = applyResumeCheckpoint(model, ckpt)
% Rebuild the live run state from a checkpoint (io.loadCompassCheckpoint output):
% restore accepted Y, RE-DERIVE workspace/aux from it (allocateWorkspace warm-starts
% the trial from the resumed Y, then rebuildSoilAuxiliaryFromAccepted recomputes the
% disposable derived cache WITHOUT one-time init mutations), restore runtime
% prognostic memory, and set the resume point so the loop continues at ckpt.kt+1.
% The diagnostics Registry is NOT restored (it is re-filled from the resume KT),
% matching the checkpoint schema.
%
% CRITICAL: do NOT call init.initializeSoilAux on resume. That path runs
% soil.water.initializeFreezingHead (mutates Y.h on frozen nodes) and always
% rewrites Y.U from the energy kernel. Both break restart equivalence with a
% continuous accepted state.
% Editor: Jizixing Zhu; Date: 2026-07-01.
% Change: 2026-07-05 Jizixing Zhu -- snow.SWE restored via Y whole-restore + materializeCheckpointStateDefaults (pre-move ckpt compat); drop the runtime.snow restore (state-ownership migration).
% Change: 2026-07-05 Jizixing Zhu -- drop the runtime.canopy_water.state_* restore (rides ckpt.Y); pass nCell to materializeCheckpointStateDefaults (increment 3 stage 2).
% Change: 2026-07-05 Jizixing Zhu -- drop runtime.phenology/lai_state restore; capture fresh-init Y.vegetation before Y=ckpt.Y and overlay legacy ckpt.runtime.{phenology,lai_state} (increment 5).
% Change: 2026-07-05 Jizixing Zhu -- drop runtime.leaf_resp + runtime.luna.{enabled,steps_per_day,state} restore; capture fresh-init Y.plant_acclimation + overlay legacy ckpt.runtime.leaf_resp/luna.state/last_*; KEEP the runtime.luna.climate_buffer handle restore (buffer moves to Y in 6b); increment 6a.
% Change: 2026-07-06 Jizixing Zhu -- drop the runtime.luna.climate_buffer restore (runtime.luna deleted); two-layer buffer compat -- overlay old ckpt.runtime.luna.climate_buffer into Y.plant_acclimation.luna_climate_buffer in BOTH the container-missing branch (pre-6a) and a new subfield-missing branch (6a-era); increment 6b.
% Change: 2026-07-19 Jizixing Zhu -- resume aux rebuild via rebuildSoilAuxiliaryFromAccepted; stop initializeSoilAux so frozen Y.h and finite Y.U stay checkpoint-exact.
% Change: 2026-09-15 Jizixing Zhu -- restore accepted root-water history after
% workspace allocation; explicitly flag legacy checkpoints without it.
    % --- accepted state + workspace/aux rebuilt from the resumed Y ---
    % capture the fresh-init Y.vegetation default before the whole-Y restore drops it (increment 5): its
    % default is a complex nested struct (PhenologyState + lai_state template), so use the fresh build, not
    % a hand-assembled default.
    freshVeg = [];
    if isfield(model.built.Y, 'vegetation'), freshVeg = model.built.Y.vegetation; end
    % capture the fresh-init Y.plant_acclimation default too (increment 6a): its default carries the
    % fixed luna_state field set with last_* = NaN + the tS-length Ta_history, so use the fresh build.
    freshPa = [];
    if isfield(model.built.Y, 'plant_acclimation'), freshPa = model.built.Y.plant_acclimation; end
    % A checkpoint written before the PlantCN ownership migration can contain
    % an empty Y.plant_cn even when the resumed configuration now enables the
    % plant C/N cycle.  No historical pool trajectory can be reconstructed
    % from that empty payload.  Keep the only defensible compatibility action
    % explicit: use the configured initial pool vector and mark the resume as
    % non-parity, rather than silently failing later in the LAI/CN consumer.
    freshPlantCn = [];
    if isfield(model.built.Y, 'plant_cn'), freshPlantCn = model.built.Y.plant_cn; end
    model.built.Y = ckpt.Y;
    nCell = size(model.built.Y.h, 2);   % nCell = 2nd dim of the resumed matric-head profile (getDomainDim not in beginRun scope)
    model.built.Y = materializeCheckpointStateDefaults(model.built.Y, ckpt.runtime, nCell);   % re-materialize moved state (pre-move ckpt compat)
    % vegetation (moved runtime.{phenology,lai_state} -> Y.vegetation): old ckpts lack Y.vegetation but carry
    % ckpt.runtime.{phenology,lai_state}; start from the fresh-init default, overlay the legacy fields.
    if ~isfield(model.built.Y, 'vegetation')
        veg = freshVeg;
        if isempty(veg), veg = struct(); end
        if isfield(ckpt.runtime, 'phenology'), veg.phenology = ckpt.runtime.phenology; end
        if isfield(ckpt.runtime, 'lai_state'), veg.lai_state = ckpt.runtime.lai_state; end
        model.built.Y.vegetation = veg;
    end
    plantCnEnabled = isfield(model.built.control, 'switches') ...
        && isfield(model.built.control.switches, 'lPlantCN') ...
        && model.built.control.switches.lPlantCN == 1;
    plantCnMissing = ~isfield(model.built.Y, 'plant_cn') ...
        || ~isstruct(model.built.Y.plant_cn) ...
        || ~isfield(model.built.Y.plant_cn, 'C') ...
        || ~isfield(model.built.Y.plant_cn, 'N');
    if plantCnEnabled && plantCnMissing
        assert(~isempty(freshPlantCn) && isfield(freshPlantCn, 'C') ...
                && isfield(freshPlantCn, 'N'), ...
            'beginRun:legacyPlantCnStateUnavailable', ...
            'The current configuration enables PlantCN but has no initial pool schema.');
        model.built.Y.plant_cn = freshPlantCn;
        model.runtime.resume.plant_cn_compatibility = ...
            'seeded_from_current_initial_pools_not_continuous_parity';
    else
        model.runtime.resume.plant_cn_compatibility = 'checkpoint_exact';
    end
    % plant_acclimation (moved runtime.leaf_resp + runtime.luna.state/last_* -> Y.plant_acclimation, increment
    % 6a): old ckpts lack Y.plant_acclimation but carry ckpt.runtime.leaf_resp.Ta_history + ckpt.runtime.luna.state
    % (+ maybe some last_*). Start from the fresh-init default (last_* = NaN), overlay the legacy fields that
    % exist; anything absent keeps the fresh default (esp. last_* = NaN).
    if ~isfield(model.built.Y, 'plant_acclimation')
        pa = freshPa;
        if isempty(pa), pa = struct(); end
        if isfield(ckpt.runtime, 'leaf_resp') && isfield(ckpt.runtime.leaf_resp, 'Ta_history')
            pa.leaf_resp.Ta_history = ckpt.runtime.leaf_resp.Ta_history;
        end
        if isfield(ckpt.runtime, 'luna') && isfield(ckpt.runtime.luna, 'state')
            base = ckpt.runtime.luna.state;
            for f = {'Vcmax25_z','Jmax25_z','PNlc_z','enzs_z','first_day'}
                if isfield(base, f{1}), pa.luna_state.(f{1}) = base.(f{1}); end
            end
        end
        if isfield(ckpt.runtime, 'luna')
            for f = {'last_Vcmax25','last_Jmax25','last_FNCa','last_PNlc','last_PNcb','last_Cab_nlc'}
                if isfield(ckpt.runtime.luna, f{1}), pa.luna_state.(f{1}) = ckpt.runtime.luna.(f{1}); end
            end
        end
        % 6b: a PRE-6a ckpt stored the serialized LUNA buffer under ckpt.runtime.luna.climate_buffer; overlay it
        % onto pa (else pa keeps freshPa's empty buffer struct).
        if isfield(ckpt.runtime, 'luna') && isfield(ckpt.runtime.luna, 'climate_buffer')
            pa.luna_climate_buffer = ckpt.runtime.luna.climate_buffer;
        end
        model.built.Y.plant_acclimation = pa;
    end
    % 6b buffer subfield compat: a 6a-era ckpt HAS Y.plant_acclimation (luna_state + leaf_resp) but LACKS the
    % luna_climate_buffer subfield (the buffer was stored under ckpt.runtime.luna.climate_buffer). The container
    % branch above did not run for it, so overlay the buffer here; else fall back to the fresh empty buffer.
    if ~isfield(model.built.Y.plant_acclimation, 'luna_climate_buffer')
        if isfield(ckpt.runtime, 'luna') && isfield(ckpt.runtime.luna, 'climate_buffer')
            model.built.Y.plant_acclimation.luna_climate_buffer = ckpt.runtime.luna.climate_buffer;
        elseif ~isempty(freshPa) && isfield(freshPa, 'luna_climate_buffer')
            model.built.Y.plant_acclimation.luna_climate_buffer = freshPa.luna_climate_buffer;
        end
    end
    % Pure aux/cache rebuild from accepted Y. allocateWorkspace is inside
    % rebuildSoilAuxiliaryFromAccepted; prognostic Y is preserved (dt_s=0 B1,
    % no Tbh fake-KT refresh). Pass completed ckpt.kt for audit/trace only.
    model.built = init.rebuildSoilAuxiliaryFromAccepted(model.built, ckpt.kt);
    if isfield(ckpt,'coupling_memory')
        memory=ckpt.coupling_memory;
        assert(isstruct(memory)&&isfield(memory,'completed_kt') ...
            &&memory.completed_kt==ckpt.kt ...
            &&isfield(memory,'root_water_element_endpoint_rate_1_s'), ...
            'beginRun:rootWaterHistory','Invalid checkpoint root-water history.');
        rootWater=memory.root_water_element_endpoint_rate_1_s;
        assert(isnumeric(rootWater)&&isequal(size(rootWater),size(model.built.workspace.flux.Srt)) ...
            &&all(isfinite(rootWater(:))), ...
            'beginRun:rootWaterHistory','Checkpoint root-water history shape/values are invalid.');
        model.built.workspace.flux.Srt=rootWater;
        model.runtime.resume.root_water_history='checkpoint_exact';
        if isfield(memory,'surface_water_history')
            history=memory.surface_water_history;
            assert(isstruct(history)&&isscalar(history) ...
                &&all(isfield(history,{'throughfall_cm_s','evap_cm_s'})), ...
                'beginRun:surfaceWaterHistory','Invalid checkpoint surface-water history.');
            for name={'throughfall_cm_s','evap_cm_s'}
                value=history.(name{1});
                assert(isnumeric(value)&&isscalar(value)&&isreal(value)&&isfinite(value), ...
                    'beginRun:surfaceWaterHistory', ...
                    'Checkpoint %s must be a finite real scalar [cm/s].',name{1});
            end
            % Change: 2026-09-30 Wenjie Liu -- hasAbruptWetting consumes
            % these previous accepted rates before the current soil solve.
            model.built.workspace.flux.Throughfall=history.throughfall_cm_s;
            model.built.workspace.flux.EVAP=history.evap_cm_s;
            model.runtime.resume.surface_water_history='checkpoint_exact';
        else
            model.runtime.resume.surface_water_history='legacy_missing_zero_seed_not_exact';
            warning('beginRun:missingSurfaceWaterHistory', ...
                ['Checkpoint lacks previous accepted surface-water fluxes; using zero. ', ...
                 'Exact continuous/restart equivalence is not guaranteed.']);
        end
    else
        model.runtime.resume.root_water_history='legacy_missing_zero_seed_not_exact';
        model.runtime.resume.surface_water_history='legacy_missing_zero_seed_not_exact';
        warning('beginRun:missingRootWaterHistory', ...
            ['Checkpoint lacks previous accepted root-water history; using zero. ', ...
             'Exact continuous/restart equivalence is not guaranteed.']);
    end
    % Checkpoints written before the soil-N ownership migration have no
    % two-time-level N memory.  Materialize a marked compatibility seed from
    % their accepted endpoint; modern checkpoints pass through unchanged.
    [model.built.Y, ~] = soil.nitrogen.materializeAcceptedState( ...
        model.built.Y, model.built.workspace.derived, model.built.domain);

    % --- runtime prognostic memory ---
    rt = ckpt.runtime;
    model.runtime.step     = rt.step;
    model.runtime.substep  = rt.substep;
    model.runtime.calendar = rt.calendar;
    % canopy_water / snow / surface / vegetation / plant_acclimation (incl. the LUNA climate buffer) all ride
    % ckpt.Y now; runtime owns no LUNA field (increment 6b). Old-ckpt buffer overlay is in the capture-fresh block above.
    model.runtime.groundwater = rt.groundwater;
    if isfield(rt, 'management')          % fertilizer once-only state (older ckpts lack it)
        model.runtime.management = rt.management;
    end

    % Change: 2026-09-29 Wenjie Liu -- restore explicit replay seam history;
    % ordinary checkpoints without this optional contract retain their path.
    if isfield(rt,'biogeochemistry_cycle')
        context=rt.biogeochemistry_cycle;
        init.validateBiogeochemistryCycleContext(context,model.built,ckpt.kt);
        assert(isequaln(context.identity,init.buildBiogeochemistryCycleIdentity(model)), ...
            'beginRun:cycleIdentity','Checkpoint replay model/source/forcing identity differs.');
        assert(context.forcing_cycle_steps==model.runtime.run_control.n_steps, ...
            'beginRun:cycleWindow','Checkpoint replay window differs from the requested window.');
        model.runtime.biogeochemistry_cycle=context;
    end
    if isfield(rt,'biogeochemistry_spinup_identity')
        assert(isequaln(rt.biogeochemistry_spinup_identity, ...
            init.buildBiogeochemistryCycleIdentity(model)), ...
            'beginRun:spinupIdentity','Checkpoint spin-up model/source/forcing identity differs.');
        assert(isfield(rt,'biogeochemistry_spinup_cycle_steps') ...
            &&rt.biogeochemistry_spinup_cycle_steps==model.runtime.run_control.n_steps, ...
            'beginRun:spinupWindow','Checkpoint spin-up window differs from the requested window.');
        model.runtime.biogeochemistry_spinup_identity=rt.biogeochemistry_spinup_identity;
        model.runtime.biogeochemistry_spinup_cycle_steps=rt.biogeochemistry_spinup_cycle_steps;
    end

    % --- resume point (runCOMPASS starts the loop at next_kt) ---
    model.runtime.run_control.last_completed_kt = ckpt.kt;
    model.runtime.run_control.next_kt = ckpt.kt + 1;
    model.runtime.resume.loaded = true;
end

% ---------------------------------------------------------------------------
function Y = materializeCheckpointStateDefaults(Y, rt, nCell)
% Re-materialize prognostic state that a whole-Y restore may drop when resuming from a PRE-move
% checkpoint (a checkpoint saved before a field's ownership moved runtime.* -> Y.*): map the legacy
% runtime.* value if present, else the init default. One clause per moved field -- THE growth point
% every state-ownership migration increment extends. nCell sizes the array-valued defaults (regional
% nCell>1 must stay 1 x nCell). New checkpoints already carry Y.<field>, so each guard is a no-op.
% Editor: Jizixing Zhu; Date: 2026-07-05.
% Change: 2026-07-05 Jizixing Zhu -- add surface.Tss clause (increment 2; defensive rt.thermal map, else NaN).
% Change: 2026-07-05 Jizixing Zhu -- add nCell arg + canopy_water clause (increment 3; per-field remap ckpt.runtime.canopy_water.state_*, else zeros(1,nCell)).
% Change: 2026-07-08 Jizixing Zhu -- add surface.Ts compatibility default [NaN;NaN], with defensive rt.thermal.Ts remap when present.
% Change: 2026-07-17 Jizixing Zhu -- materialize liquid-water and energy
%   defaults for checkpoints from the former scalar-SWE snow scheme.
    % snow.SWE (moved runtime.snow -> Y.snow, 2026-07-05)
    if ~isfield(Y, 'snow') || ~isfield(Y.snow, 'SWE')
        swe0 = 0;
        if isfield(rt, 'snow') && isfield(rt.snow, 'SWE'), swe0 = rt.snow.SWE; end   % legacy path map
        Y.snow = struct('SWE', swe0);
    end
    if ~isfield(Y.snow, 'liquid_water')
        Y.snow.liquid_water = 0;
    end
    if ~isfield(Y.snow, 'energy')
        phys = equations.physicalConstants();
        Y.snow.energy = -phys.lambda_fus * 1e-3 * max(Y.snow.SWE, 0);
    end
    % surface.Ts/Tss (moved runtime.thermal -> Y.surface, 2026-07-05/2026-07-08).
    % runtime.thermal was never serialized by the new checkpoint writer, so the rt.thermal
    % maps are defensive only; old checkpoints usually default both fields to NaN.
    if ~isfield(Y, 'surface') || ~isstruct(Y.surface), Y.surface = struct(); end
    if ~isfield(Y.surface, 'Ts')
        Ts0 = NaN(2, 1);
        if isfield(rt, 'thermal') && isfield(rt.thermal, 'Ts') && numel(rt.thermal.Ts) == 2
            Ts0 = rt.thermal.Ts(:);
        end
        Y.surface.Ts = Ts0;
    end
    if ~isfield(Y.surface, 'Tss')
        tss0 = NaN;
        if isfield(rt, 'thermal') && isfield(rt.thermal, 'Tss'), tss0 = rt.thermal.Tss; end
        Y.surface.Tss = tss0;
    end
    % canopy_water state (moved runtime.canopy_water -> Y.canopy_water, 2026-07-05). Old ckpts stored
    % state under ckpt.runtime.canopy_water.state_* (stateOnlyCanopyWater); per-field remap, else zeros(1,nCell).
    if ~isfield(Y, 'canopy_water') || ~isstruct(Y.canopy_water), Y.canopy_water = struct(); end
    rc = struct();
    if isfield(rt, 'canopy_water') && isstruct(rt.canopy_water), rc = rt.canopy_water; end
    cwFields = {'state_h2ocan', 'state_h2ocan_sun', 'state_h2ocan_sha'};
    for i = 1:numel(cwFields)
        f = cwFields{i};
        if ~isfield(Y.canopy_water, f)
            if isfield(rc, f), Y.canopy_water.(f) = rc.(f); else, Y.canopy_water.(f) = zeros(1, nCell); end
        end
    end
    % Y.U (prognostic internal energy, added 2026-07-14).
    % Old checkpoints lack this field.  Recompute from the restored
    % (h,T,P_g) using the unified kernel rather than defaulting to NaN.
    if ~isfield(Y, 'U') || all(isnan(Y.U(:)))
        Y.U = nan(size(Y.T));
    end
end

% ---------------------------------------------------------------------------
function model = applyGroundwaterSetup(model)
% Groundwater-coupling pre-loop setup (beginRun 3.3), run once when
% runtime.groundwater.GroundwaterCoupling == 1: compute the water-table depth,
% locate the bottom (water-table) layer, seed the bottom temperature from the
% initial/resumed soil profile, and switch the bottom water BC to Dirichlet.
% Dunnian runoff (precip modification) is deferred — see the 3.3 note in beginRun.
% Editor: Jizixing Zhu; Date: 2026-07-01.
    gw = model.runtime.groundwater;
    domain = model.built.domain;

    gwDep = groundwater.calculateGroundWaterDepth(gw.topLevel, gw.headBotmLayer, domain.tot_depth);
    if gwDep > domain.tot_depth
        error('beginRun:groundwaterBelowColumn', ...
            ['water table depth (%.1f cm) is below the soil column (%.1f cm); ', ...
             'enlarge total soil thickness.'], gwDep, domain.tot_depth);
    end

    grid = struct('NL', domain.NL, 'NN', domain.NN);
    [idx, idxR] = groundwater.calculateIndexBottomLayer(gw.soilLayerThickness, gwDep, grid);

    gw.gw_Dep = gwDep;
    gw.indexBotmLayer   = idx;    % STEMMUS bottom-to-top (1 = deepest)
    gw.indexBotmLayer_R = idxR;   % MODFLOW top-to-bottom
    if ~isfinite(gw.tempBotm)
        gw.tempBotm = model.built.Y.T(idx);   % Y.T is bottom-up; idx is bottom-to-top
    end
    gw.status = "enabled";
    model.runtime.groundwater = gw;

    % The external water table supplies an accepted lower-column boundary
    % state as well as the next Richards boundary row. Initialise the
    % MODFLOW-saturated segment hydrostatically before the first trial: these
    % nodes are not unconstrained Picard unknowns and leaving an old deep-dry
    % profile there creates a fictitious one-step multi-kilometre head jump.
    fixedHead = gw.headBotmLayer - gw.topLevel ...
        + gw.soilLayerThickness(idxR);
    h = model.built.Y.h(:);
    h(idx) = fixedHead;
    for node = idx:-1:2
        h(node - 1) = h(node) + domain.dz(node - 1);
    end
    model.built.Y.h = reshape(h, size(model.built.Y.h));
    model.built.workspace.step = solver.warmStartTrial( ...
        model.built.workspace.step, model.built.Y);
    model.built = init.rebuildSoilAuxiliaryFromAccepted( ...
        model.built, 0);
    model.built = synchronizeExternalGroundwaterAcceptedState(model.built);

    % Dirichlet bottom water BC fed by the water table (new-structure field).
    model.built.boundary.NBCh_bottom = 1;
end

% ---------------------------------------------------------------------------
function built = synchronizeExternalGroundwaterAcceptedState(built)
% Synchronize accepted soil storage after imposing an external water table.
%
% A MODFLOW water table is an initial/boundary condition, not a Richards
% increment over COMPASS's first outer step.  Once its hydrostatic head has
% replaced the inherited profile, every accepted water and energy storage used
% in the next time derivative must represent that same (h,T,P_g) state.
% Leaving Y.Theta_L or Y.U at the former dry profile makes the first Richards
% residual interpret the externally prescribed saturated segment as an
% unphysical one-hour source of water and internal energy.
%
% Editor: Jizixing Zhu; Date: 2026-07-20.
    d = built.workspace.derived;
    assert(isfield(d, 'Theta_LL') && isfield(d, 'Theta_UU') ...
            && isfield(d, 'Theta_II'), ...
        'beginRun:missingGroundwaterStorage', ...
        'Rebuilt soil auxiliary state lacks accepted phase storage.');
    built.Y.Theta_L = d.Theta_LL;
    built.Y.Theta_U = d.Theta_UU;
    built.Y.Theta_I = d.Theta_II;
    built.Y.theta = [d.Theta_LL(:,1); d.Theta_LL(end,2)];

    L_f = built.workspace.step.L_f;
    dt_s = built.control.time.Delt_t;
    [estate, ederiv] = soil.thermal.refreshInternalEnergyState( ...
        built, built.Y.h(:), built.Y.T(:), built.Y.p_gas(:), ...
        dt_s, L_f, 0);
    U_endpoint = soil.thermal.evaluateInternalEnergyEndpoints( ...
        estate, ederiv, struct('L_f', L_f, 'validate', true));
    built = init.setPrognosticInternalEnergy(built, U_endpoint);
    built.workspace.step = solver.warmStartTrial( ...
        built.workspace.step, built.Y);
end

% ---------------------------------------------------------------------------
function ensureOutputDirs(output)
% Create the run output directories once (idempotent: skip existing).
% Editor: Jizixing Zhu; Date: 2026-06-29.
    dirs = {output.output_dir, output.checkpoint_dir, output.nc_dir, output.csv_dir};
    for i = 1:numel(dirs)
        d = dirs{i};
        if ~isempty(d) && ~isfolder(d)
            [ok, msg] = mkdir(d);
            if ~ok
                error('beginRun:mkdirFailed', ...
                    'could not create output directory %s: %s', d, msg);
            end
        end
    end
end
