function model = prepareRuntime(model)
%PREPARERUNTIME  Build pre-loop runtime services for an initialized model.
%
%   model = init.prepareRuntime(model) adds model.runtime to a model returned by
%   init.initialize. This stage owns run services that are needed before the KT
%   loop starts but are not physical input/state construction:
%       diagnostics registry, checkpoint/resume policy, Picard diagnostics,
%       the LUNA 10-day climate buffer HANDLE (state/config moved out -- see below),
%       groundwater status, soil-mode status and SCOPE runtime scaffolding (runtime.scope;
%       the static optics cache is built.p.scope_static, an init cache, 2026-07-05 increment 4).
%       (canopy-water is no longer a runtime container -- state -> Y.canopy_water,
%       histories -> built.diagnostics.canopy_water, config -> params.veg.canopy, 2026-07-05.)
%       (phenology + lai_state are no longer runtime -- Y.vegetation, 2026-07-05 increment 5.)
%       (LUNA state/last_* + leaf_resp + the 10-day climate buffer are no longer runtime -- they are
%       Y.plant_acclimation.{luna_state,leaf_resp,luna_climate_buffer}; LUNA/leaf_resp config ->
%       built.p.plant_acclimation; runtime owns no LUNA field, 2026-07-05 increment 6a + 6b.)
%
%   Inputs:
%     model  [inout] struct  init.initialize output. Uses model.built but does
%                            not modify built.Y, built.parameters, built.forcing
%                            or built.workspace.
%
%   Returns:
%     model.runtime output/diagnostics/checkpoint/resume/groundwater/
%                   picard/soil_modes/scope containers, plus
%                   pre-loop state slots substep/calendar/debug/run_control
%                   (phenology + lai_state are no longer runtime -- they are Y.vegetation, 2026-07-05 increment 5)
%                   (leaf_resp + the WHOLE luna container are no longer runtime -- LUNA state/last_* + the
%                    climate buffer are Y.plant_acclimation.{luna_state,luna_climate_buffer}, config in
%                    built.p.plant_acclimation; runtime owns no LUNA field, 2026-07-05 increment 6a + 6b)
%                   (mutable counters runOneStep advances; strategy/params stay in
%                   control.numerics/time).
%
%   Closure:
%     No fluxes are introduced and no prognostic state is advanced here. This is
%     a pre-loop service allocator, matching the old COMPASS.m blocks before the
%     first KT while preserving the new Y/workspace separation.
% Editor: Jizixing Zhu; Date: 2026-06-27.
% Change: 2026-06-28 Jizixing Zhu — env reads go through io.getEnv (COMPASS_* canonical, STEMMUS_* fallback).
% Change: 2026-06-29 Jizixing Zhu — add pre-loop runtime slots substep/calendar/debug/run_control/leaf_resp.
% Change: 2026-06-29 Jizixing Zhu — reorder local functions to match the assembly-line call order; helpers grouped at end.
% Change: 2026-07-05 Jizixing Zhu — remove runtime.snow (SWE now prognostic Y.snow, buildInitialState); state-ownership migration.
% Change: 2026-07-05 Jizixing Zhu — remove runtime.thermal (Tss now prognostic Y.surface, buildInitialState); state-ownership migration increment 2.
% Change: 2026-07-08 Jizixing Zhu -- Y.surface now also owns Ts [2x1], the accepted SCOPE/ebal thermal.Ts carry.
% Change: 2026-07-05 Jizixing Zhu — prepareCanopyWater stripped to state-only (histories -> built.diagnostics, config -> params.veg.canopy); delete now-unused readEnvDouble local (-> shared io.getEnvDouble); increment 3 stage 1.
% Change: 2026-07-05 Jizixing Zhu — delete prepareCanopyWater + runtime.canopy_water (state now prognostic Y.canopy_water, buildInitialState); state-ownership migration increment 3 stage 2.
% Change: 2026-07-05 Jizixing Zhu — delete prepareScopeStatic + runtime.scope_static (now built.p.scope_static via init.buildScopeStatic in buildState); prepareScopeRuntime reads built.p; runtime stays a pure allocator; increment 4.
% Change: 2026-07-05 Jizixing Zhu — delete preparePhenology + prepareLaiState + runtime.phenology/lai_state (now prognostic Y.vegetation, buildInitialState); increment 5.
% Change: 2026-07-05 Jizixing Zhu -- delete prepareLeafResp + runtime.leaf_resp (Y.plant_acclimation.leaf_resp), reduce prepareLuna to climate_buffer only (state -> Y.plant_acclimation.luna_state, config -> built.p.plant_acclimation); increment 6a.
% Change: 2026-07-06 Jizixing Zhu -- delete prepareLuna + runtime.luna entirely; the climate buffer is now Y.plant_acclimation.luna_climate_buffer (serialized struct); runtime owns no LUNA field; increment 6b.
% Change: 2026-07-11 Jizixing Zhu -- dispatch laboratory water sites to the minimal water-only runtime allocator.
% Change: 2026-07-13 Jizixing Zhu -- apply the canonical DIAG_INTERVAL
%         environment override to the new runtime checkpoint service.
% Change: 2026-07-20 Jizixing Zhu -- materialize the accepted t0 clock for a first-step failure checkpoint.

    if nargin < 1 || ~isstruct(model) || ~isfield(model, 'built')
        error('prepareRuntime:badModel', ...
            'prepareRuntime expects the struct returned by init.initialize.');
    end

    built = model.built;
    if isfield(built, 'control') && isfield(built.control, 'simulation') ...
            && isfield(built.control.simulation, 'is_laboratory') ...
            && built.control.simulation.is_laboratory
        model = init.prepareWaterOnlyRuntime(model);
        return;
    end
    tS = built.control.time.tS;
    NL = built.domain.NL;
    NN = built.domain.NN;
    nCell = getDomainDim(built.domain, 'nCell', 1);

    runtime = struct();

    runtime.output      = prepareOutput(model);
    runtime.diagnostics = prepareDiagnostics(built);
    runtime.checkpoint  = prepareCheckpoint(model, runtime.output);
    runtime.resume      = prepareResume();
    runtime.groundwater = prepareGroundwater(built);
    % phenology + lai_state moved to prognostic Y.vegetation (buildInitialState) 2026-07-05 -- runtime no longer owns them.
    % canopy_water state moved to prognostic Y.canopy_water (buildInitialState) 2026-07-05 -- runtime no longer owns it.
    % snow.SWE moved to prognostic Y.snow (buildInitialState) 2026-07-05 -- runtime no longer owns it.
    % surface temperatures Ts/Tss moved to prognostic Y.surface -- runtime no longer owns them.
    runtime.picard      = preparePicard(tS, built);
    runtime.soil_modes  = prepareSoilModes(built);
    % scope_static moved to the init cache built.p.scope_static (buildState) 2026-07-05 -- runtime no longer owns it.
    runtime.scope       = init.prepareScopeRuntime(built, built.p.scope_static);

    % pre-loop runtime state slots the loop owns (strategy/params live in
    % control.numerics/time; these are the mutable counters runOneStep advances).
    runtime.substep     = prepareSubstep(built);
    runtime.calendar    = prepareCalendar(built);
    runtime.step        = prepareInitialStep(built);
    runtime.debug       = prepareDebug();
    runtime.run_control = prepareRunControl(built);
    runtime.management  = prepareManagement(built);

    runtime.dimensions = struct('tS', tS, 'NL', NL, 'NN', NN, 'nCell', nCell);
    runtime.status = "prepared";
    runtime.created_at = datetime("now");

    model.runtime = runtime;
    model.status = "runtime_prepared";
end

% ---------------------------------------------------------------------------
function s = prepareInitialStep(built)
%PREPAREINITIALSTEP  Materialize the accepted t0 clock for restart safety.
%
%   A checkpoint may be requested after the first outer step is rejected.
%   That checkpoint represents the initialized accepted state t0, so it needs
%   the same runtime.step schema as an ordinary accepted endpoint rather than
%   an absent field or a fabricated KT=1 result.
%
%   Editor: Jizixing Zhu; Date: 2026-07-20.
    time = built.control.time;
    start = time.start_time;
    s = struct('KT', 0, 'TIME', 0, 'time_step', 0, ...
        'elapsed_seconds', 0, 'clock_time', start, 'NoTime', 0, ...
        'forcing_index', 0, 'day_hour', 0, ...
        'current_day', day(start), 'current_hour', hour(start));
end

% ===========================================================================
%  Service builders — listed in the same order the assembly line calls them.
% ===========================================================================

% ---------------------------------------------------------------------------
function o = prepareOutput(model)
% Resolve the run output location and common writer subdirectories.
% Editor: Jizixing Zhu; Date: 2026-06-28.
    [outputDir, envName] = io.getEnv('OUTPUT');
    source = "env:" + string(envName);
    if isempty(outputDir)
        [outputDir, source] = configuredOutputDir(model);
    end

    outputCfg = struct();
    if isfield(model.raw.config, 'output')
        outputCfg = model.raw.config.output;
    end

    o = struct();
    o.output_dir = outputDir;
    o.source = source;
    o.checkpoint_dir = fullfile(outputDir, 'checkpoints');
    o.nc_dir = fullfile(outputDir, 'nc');
    o.csv_dir = fullfile(outputDir, 'csv');
    o.legacy_csv = getLogicalField(outputCfg, 'legacy_csv', false);
    o.periodic_partial = getLogicalField(outputCfg, 'periodic_partial', false);
end

% ---------------------------------------------------------------------------
function d = prepareDiagnostics(built)
% Create the central diagnostics registry and register all known output fields.
% Editor: Jizixing Zhu; Date: 2026-06-27.
% Change: 2026-08-20 Codex -- support a sparse driver-generation diagnostic
% profile; this changes output allocation only, never model physics.
    cfg = init.makeLegacyRuntimeCfg(built);
    profile = lower(strtrim(io.getEnv('DIAGNOSTIC_PROFILE')));
    if isempty(profile)
        profile = 'full';
    end
    sparseProfile = strcmp(profile, 'driver');
    reg = diagnostics.Registry(built.control.time.tS, built.domain.NL, ...
        built.domain.NN, ~sparseProfile);
    diagnostics.registerAllDiagnostics(reg, cfg, profile);

    d = struct();
    d.registry = reg;
    d.profile = profile;
    d.online_gs_diag_enabled = strcmpi(strtrim(io.getEnv('ONLINE_GS_DIAG')), '1');
    d.n_variables = numel(reg.list());
    d.n_groups = numel(reg.groupNames());
end

% ---------------------------------------------------------------------------
function c = prepareCheckpoint(model, output)
% Collect restart/partial-output policy without writing files yet.
% Editor: Jizixing Zhu; Date: 2026-06-27.
    built = model.built;

    checkpointInterval = built.control.numerics.checkpointInterval;
    if checkpointInterval <= 0
        checkpointInterval = 24;
    end

    diagnosticInterval = built.control.numerics.diagnosticInterval;
    diagnosticIntervalEnv = str2double(io.getEnv('DIAG_INTERVAL'));
    if isfinite(diagnosticIntervalEnv) && diagnosticIntervalEnv > 0
        diagnosticInterval = diagnosticIntervalEnv;
    end
    if diagnosticInterval <= 0
        diagnosticInterval = 360;
    end

    partialInterval = checkpointInterval;
    partialEnv = str2double(io.getEnv('PARTIAL_OUTPUT_INTERVAL'));
    if isfinite(partialEnv) && partialEnv > 0
        partialInterval = partialEnv;
    end

    periodicPartial = false;
    if isfield(model.raw.config, 'output') && isfield(model.raw.config.output, 'periodic_partial')
        periodicPartial = logical(model.raw.config.output.periodic_partial);
    end

    c = struct();
    c.output_dir = output.output_dir;
    c.checkpoint_dir = output.checkpoint_dir;
    c.checkpoint_interval = checkpointInterval;
    c.diagnostic_interval = diagnosticInterval;
    c.partial_output_interval = partialInterval;
    c.write_periodic_partial = periodicPartial;
    c.in_optimizer = ~isempty(io.getEnv('WORKER_ROOT'));
end

% ---------------------------------------------------------------------------
function r = prepareResume()
% Record resume request status. Loading is kept outside this pure allocator.
% Editor: Jizixing Zhu; Date: 2026-06-27.
    resumeFile = strtrim(io.getEnv('RESUME_FILE'));
    r = struct();
    r.requested = ~isempty(resumeFile);
    r.file = resumeFile;
    r.loaded = false;
    if r.requested && ~isfile(resumeFile)
        error('prepareRuntime:missingResumeFile', ...
            'STEMMUS_RESUME_FILE does not exist: %s', resumeFile);
    end
end

% ---------------------------------------------------------------------------
function g = prepareGroundwater(built)
% Build the disabled-by-default groundwater runtime container.
% Editor: Jizixing Zhu; Date: 2026-06-27.
    NL = built.domain.NL;
    NN = built.domain.NN;
    nCell = getDomainDim(built.domain, 'nCell', 1);
    topDownDz = flip(built.domain.dz(:));

    soilLayerThickness = zeros(NN, 1);
    for i = 2:NL
        soilLayerThickness(i) = soilLayerThickness(i - 1) + topDownDz(i - 1);
    end
    soilLayerThickness(NN) = built.domain.tot_depth;

    g = struct();
    g.GroundwaterCoupling = 0;
    g.headBotmLayer = 1750.0;
    g.topLevel = 2000.0;
    g.tempBotm = NaN;
    g.indexBotmLayer = 40;
    if nCell > 1
        soilLayerThickness = repmat(soilLayerThickness, 1, nCell);
    end
    g.soilLayerThickness = soilLayerThickness;
    g.status = "disabled";
end

% ---------------------------------------------------------------------------
function p = preparePicard(tS, built)
% Allocate compact Picard diagnostics.
% Editor: Jizixing Zhu; Date: 2026-06-27.
    cfg = init.makeLegacyRuntimeCfg(built);
    [saveHistory, summary, history] = solver.initPicardDiagnostics(tS, cfg);
    p = struct();
    p.save_history = saveHistory;
    p.summary = summary;
    p.history = history;
end

% ---------------------------------------------------------------------------
function s = prepareSoilModes(built)
% Record advanced soil-water mode status; old structure adapters are pending.
% Editor: Jizixing Zhu; Date: 2026-06-27.
    sw = built.control.switches;
    s = struct();
    s.hysteresis.enabled = isfield(sw, 'lHystrs') && logical(sw.lHystrs);
    if s.hysteresis.enabled
        s.hysteresis.status = "initialized";
    else
        s.hysteresis.status = "disabled";
    end
    s.dual_porosity.enabled = isfield(sw, 'lDualPor') && logical(sw.lDualPor);
    if s.dual_porosity.enabled
        s.dual_porosity.status = "candidate_transaction";
    else
        s.dual_porosity.status = "disabled";
    end
end


% ---------------------------------------------------------------------------
function s = prepareSubstep(built)
% Adaptive-substep runtime state slots; the strategy lives in control.numerics.substep.
% Editor: Jizixing Zhu; Date: 2026-06-29.
% Change: 2026-06-29 Jizixing Zhu — add last_iters slot (solver.decideSubstep grow gate).
    dt = built.control.time.Delt_t;
    s = struct();
    s.Delt_t        = dt;    % [s] current step length (starts at the outer dt)
    s.Delt_old      = dt;    % [s] previous step length (temporal-derivative term)
    s.retry_count   = 0;     % [-] substep retries on the current step
    s.success_count = 0;     % [-] consecutive accepted steps (-> grow dt)
    s.hold_count    = 0;     % [-] steps left to hold small dt after a reject
    s.last_iters    = 0;     % [-] solver iterations of the last accepted step (grow gate)
end

% ---------------------------------------------------------------------------
function c = prepareCalendar(built)
% Wall-clock trackers the loop advances (year reset + daily phenology hooks).
% Editor: Jizixing Zhu; Date: 2026-06-29.
    dt = built.control.time.Delt_t;
    c = struct();
    c.current_year  = year(built.control.time.start_time);  % year tracker (max-temp reset)
    c.current_date  = NaT;                                  % current calendar day (set in loop)
    c.steps_per_day = max(1, round(86400 / dt));            % [-] outer steps per day
end

% ---------------------------------------------------------------------------
function d = prepareDebug()
% Optional diagnostic toggles read from the environment.
% Editor: Jizixing Zhu; Date: 2026-06-29.
    d = struct();
    d.spike_diag_enabled = strcmpi(strtrim(io.getEnv('SPIKE_DIAG')), '1');
end

% ---------------------------------------------------------------------------
function r = prepareRunControl(built)
% Run-window control: loop bounds + optional short-window override (smoke / N-step).
% Editor: Jizixing Zhu; Date: 2026-06-29.
% Change: 2026-07-01 Jizixing Zhu — add next_kt/last_completed_kt (loop start; resume overrides).
    tS = built.control.time.tS;
    r = struct();
    r.n_steps = tS;          % [-] outer steps to run (default: full window)
    r.max_kt  = NaN;         % [-] requested MAX_KT override (NaN = none)
    r.next_kt = 1;           % [-] loop start step (resume sets this; runOneStep advances it)
    r.last_completed_kt = 0; % [-] last accepted step (0 = none yet)
    maxKt = str2double(io.getEnv('MAX_KT'));
    if isfinite(maxKt) && maxKt > 0 && maxKt < tS
        r.max_kt  = maxKt;
        r.n_steps = maxKt;
    end
end

% ---------------------------------------------------------------------------
function m = prepareManagement(built)
% Serializable management-event state (checkpoint/resume safe). Fertilizer once-only
% flags live HERE (not in a closure) so a resume does not re-fire past events; the
% matcher soil.nitrogen.fertilizer.matchFertilizerEvent reads/updates them each step.
% Editor: Jizixing Zhu; Date: 2026-07-02.
% Change: 2026-08-16 Codex -- allocate accepted-only organic-fertilizer flags.
    nFert = 0;
    if isfield(built, 'sources') && isfield(built.sources, 'fertilization_schedule') ...
            && isfield(built.sources.fertilization_schedule, 'n_events')
        nFert = built.sources.fertilization_schedule.n_events;
    end
    m = struct();
    m.fertilizer = struct('applied', false(1, nFert));
    nOrganic = 0;
    if isfield(built,'sources') ...
            && isfield(built.sources,'organic_fertilization_schedule')
        nOrganic = built.sources.organic_fertilization_schedule.n_events;
    end
    m.organic_fertilizer = struct('applied',false(1,nOrganic));
end

% ===========================================================================
%  Shared helpers (used by several builders above).
% ===========================================================================

% ---------------------------------------------------------------------------
function [outputDir, source] = configuredOutputDir(model)
% Read [output].output_dir with a deterministic site-relative fallback.
% Editor: Jizixing Zhu; Date: 2026-06-28.
    source = "default:data/output/output_initialized";
    dataDir = fileparts(fileparts(model.raw.siteDir));
    outputDir = fullfile(dataDir, 'output', 'output_initialized');
    if ~isfield(model.raw.config, 'output')
        return;
    end

    outputCfg = model.raw.config.output;
    for name = {'output_dir','output_path','path'}
        if isfield(outputCfg, name{1}) && ~isempty(outputCfg.(name{1}))
            rawPath = char(outputCfg.(name{1}));
            if isAbsolutePath(rawPath)
                outputDir = rawPath;
            else
                outputDir = fullfile(model.raw.siteDir, rawPath);
            end
            source = "site.toml:output." + string(name{1});
            return;
        end
    end
end

% ---------------------------------------------------------------------------
function tf = isAbsolutePath(pathValue)
% True for Windows drive/UNC paths and POSIX absolute paths.
% Editor: Jizixing Zhu; Date: 2026-06-28.
    pathValue = char(pathValue);
    tf = ~isempty(regexp(pathValue, '^[A-Za-z]:[\\/]', 'once')) || ...
        startsWith(pathValue, '\\') || startsWith(pathValue, '/');
end

% ---------------------------------------------------------------------------
function tf = getLogicalField(s, name, defaultValue)
% Read one optional logical output setting.
% Editor: Jizixing Zhu; Date: 2026-06-28.
    tf = defaultValue;
    if isfield(s, name) && ~isempty(s.(name))
        tf = logical(s.(name));
    end
end

% ---------------------------------------------------------------------------
% (makeLegacyRuntimeCfg moved to src/+init/makeLegacyRuntimeCfg.m — package
%  function, single source shared by prepareRuntime + runOneStep partial output.)

% ---------------------------------------------------------------------------
function value = getDomainDim(domain, name, fallback)
% Read a domain dimension with fallback for single-column legacy paths.
% Editor: Jizixing Zhu; Date: 2026-06-27.
    value = fallback;
    if isfield(domain, name) && ~isempty(domain.(name))
        value = domain.(name);
    end
end
