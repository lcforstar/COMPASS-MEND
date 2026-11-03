function src = buildSources(rawInput)
%BUILDSOURCES  Collect source/sink management events into built.sources.
%
%   src = init.buildSources(rawConfig) gathers the management source/sink terms
%   (fertilizer events from [fertilization]; irrigation etc. would be added
%   here) into a clean struct. These are management inputs, not boundary
%   conditions, so they live under built.sources rather than built.boundary.
%   Pure config translation -- no physics.
%
%   Inputs:
%     rawConfig  struct  raw.config (uses optional .fertilization)
%
%   Returns src with:
%     src.fertilization           struct  raw event view: events_time {cellstr},
%                                 events_amount [kgN/ha], events_type {cellstr},
%                                 NH4_ratio/NO3_ratio/urea_ratio, depth_start/end [cm]
%     src.fertilization_schedule  struct  normalized schedule consumed by the run:
%                                 DEPTHSTART/DEPTHEND, t_fert [datetime], fert_kg_ha,
%                                 NH4_ratio/NO3_ratio/urea_ratio, pool_entry {cellstr},
%                                 depth_cm, n_events. pool_entry names the target state
%                                 pool per event (nh4_granule/nh4_liquid/no3_liquid/urea;
%                                 defaulted from events_type). Plain data (no closure);
%                                 the once-only applied-state lives in
%                                 runtime.management.fertilizer, and matching is
%                                 soil.nitrogen.fertilizer.matchFertilizerEvent.
%
%   Reference: legacy BoundarySettings (fertilizer field names/units).
%
%   Editor: Jizixing Zhu; Date: 2026-06-26.
%   Change: 2026-08-16 Codex -- parse AOM fertilizer only while the soil
%     organic cycle is enabled; configured events remain inert when disabled.
%   Change: 2026-07-02 Jizixing Zhu — add src.fertilization_schedule (moved from buildBoundary bc.Fertilizer): fertilizer is a management source, not a boundary condition.
%   Change: 2026-07-02 Jizixing Zhu — add per-event pool_entry (config-first, else defaulted from events_type) routing solid NH4 (nh4_granule) vs dissolved NH4/NO3/urea.
%   Change: 2026-07-02 Jizixing Zhu — build-time fail-fast: an NH4-bearing event must route to an NH4 pool (validatePoolFormConsistency).
% Change: 2026-10-30 Wenjie Liu -- append paired external root-residue
%   events to the existing once-only organic-management schedule.

    if isfield(rawInput,'config'), rawConfig=rawInput.config;
    else, rawConfig=rawInput; end
    src = struct();

    if isfield(rawConfig, 'fertilization')
        f = rawConfig.fertilization;
        src.fertilization = struct( ...
            'events_time',   {cellstr(string(f.events_time))}, ...
            'events_amount', f.events_amount(:)', ...
            'events_type',   {cellstr(string(f.events_type))}, ...
            'NH4_ratio',     f.NH4_ratio(:)', ...
            'NO3_ratio',     f.NO3_ratio(:)', ...
            'urea_ratio',    f.urea_ratio(:)', ...
            'depth_start',   f.depth_start, ...
            'depth_end',     f.depth_end);
    else
        src.fertilization = struct('events_time', {{}}, 'events_amount', [], ...
            'events_type', {{}}, 'NH4_ratio', [], 'NO3_ratio', [], 'urea_ratio', [], ...
            'depth_start', [], 'depth_end', []);
    end

    src.fertilization_schedule = buildFertilizerSchedule(rawConfig);
    organicEnabled = false;
    if isfield(rawConfig, 'switches')
        switches = init.buildSwitches(rawConfig);
        organicEnabled = logical(switches.lSoilOrganic);
    end
    src.organic_fertilization_schedule = ...
        soil.organic.buildFertilizationSchedule(rawConfig, organicEnabled);
    src.organic_fertilization_schedule = ...
        soil.organic.appendStrawReturnSchedule( ...
        src.organic_fertilization_schedule,rawConfig,organicEnabled);
    src.organic_fertilization_schedule = ...
        soil.organic.appendRootResidueSchedule( ...
        src.organic_fertilization_schedule,rawConfig,organicEnabled);
end

% ---------------------------------------------------------------------------
function F = buildFertilizerSchedule(rawConfig)
% Normalized fertilizer event schedule from raw.config.fertilization (site.toml
% [fertilization]): datetime t_fert + per-event amount/form ratios + injection depth
% band, as plain data (no closure). n_events = 0 when there is no [fertilization]
% section or no events (fertilizer then a no-op). Consumed by runOneStep.applyFertStep
% via soil.nitrogen.fertilizer.matchFertilizerEvent; once-only state is in
% runtime.management.fertilizer (checkpoint/resume safe).
% Editor: Jizixing Zhu; Date: 2026-07-02.
    F = struct('DEPTHSTART', 0, 'DEPTHEND', 0, 't_fert', datetime.empty, ...
               'fert_kg_ha', [], 'NH4_ratio', [], 'NO3_ratio', [], 'urea_ratio', [], ...
               'pool_entry', {{}}, 'depth_cm', 0, 'n_events', 0);
    if ~isfield(rawConfig, 'fertilization') || ~isstruct(rawConfig.fertilization)
        return;
    end
    f = rawConfig.fertilization;
    if ~isfield(f, 'events_time') || isempty(f.events_time)
        return;
    end
    F.DEPTHSTART = f.depth_start;
    F.DEPTHEND   = f.depth_end;
    F.t_fert     = datetime(string(f.events_time), 'InputFormat', 'yyyy-MM-dd''T''HH:mm:ss');
    F.fert_kg_ha = f.events_amount(:).';
    F.NH4_ratio  = f.NH4_ratio(:).';
    F.NO3_ratio  = f.NO3_ratio(:).';
    F.urea_ratio = f.urea_ratio(:).';
    F.depth_cm   = F.DEPTHEND - F.DEPTHSTART;
    F.n_events   = numel(F.t_fert);
    % pool_entry: which model state pool each event enters. Explicit config wins;
    % otherwise defaulted from events_type. This decides SOLID NH4 (nh4_granule, slow
    % release) vs dissolved NH4 (nh4_liquid), and keeps nitrate/urea from being
    % mis-routed to the granule pool.
    if isfield(f, 'pool_entry') && ~isempty(f.pool_entry)
        F.pool_entry = cellstr(string(f.pool_entry));
    else
        F.pool_entry = defaultPoolEntry(f.events_type);
    end
    validatePoolEntry(F.pool_entry, F.n_events);
    validatePoolFormConsistency(F);
end

% ---------------------------------------------------------------------------
function pe = defaultPoolEntry(events_type)
% Map events_type -> pool_entry (nh4_granule / no3_liquid / urea). Mixed or custom
% types have no unambiguous single pool and MUST set pool_entry explicitly.
% Editor: Jizixing Zhu; Date: 2026-07-02.
    types = cellstr(string(events_type));
    pe = cell(1, numel(types));
    for i = 1:numel(types)
        switch lower(strtrim(types{i}))
            case 'ammonium',               pe{i} = 'nh4_granule';   % solid NH4 -> granule
            case 'nitrate',                pe{i} = 'no3_liquid';    % NO3 fast-dissolve
            case {'urea', 'lime_nitrogen'}, pe{i} = 'urea';
            otherwise
                error('buildSources:fertPoolEntry', ...
                    ['fertilization event %d type "%s" has no default pool_entry; ', ...
                     'set [fertilization].pool_entry explicitly.'], i, types{i});
        end
    end
end

% ---------------------------------------------------------------------------
function validatePoolEntry(pe, nEvents)
% pool_entry length must match events, and each value must be a known pool.
% Editor: Jizixing Zhu; Date: 2026-07-02.
    if numel(pe) ~= nEvents
        error('buildSources:fertPoolEntryLen', ...
            'pool_entry has %d entries but there are %d fertilization events.', numel(pe), nEvents);
    end
    known = {'nh4_granule', 'nh4_liquid', 'no3_liquid', 'urea'};
    for i = 1:numel(pe)
        if ~ismember(lower(strtrim(pe{i})), known)
            error('buildSources:fertPoolEntryValue', ...
                'pool_entry "%s" (event %d) not in {nh4_granule, nh4_liquid, no3_liquid, urea}.', pe{i}, i);
        end
    end
end

% ---------------------------------------------------------------------------
function validatePoolFormConsistency(F)
% Fail-fast (residual guard): an event that carries an NH4 fraction (NH4_ratio>0) must
% route it to an NH4 pool (nh4_granule/nh4_liquid). Otherwise the granule-vs-dissolved
% split is undefined and matchFertilizerEvent could not allocate the NH4 mass. NO3/urea
% fractions have fixed pools, so they need no such check.
% Editor: Jizixing Zhu; Date: 2026-07-02.
    nh4Pools = {'nh4_granule', 'nh4_liquid'};
    for i = 1:F.n_events
        if F.NH4_ratio(i) > 0 && ~ismember(lower(strtrim(F.pool_entry{i})), nh4Pools)
            error('buildSources:fertPoolFormMismatch', ...
                ['fertilization event %d carries NH4 (NH4_ratio=%.3g) but pool_entry="%s" ', ...
                 'is not an NH4 pool; use nh4_granule or nh4_liquid.'], i, F.NH4_ratio(i), F.pool_entry{i});
        end
    end
end
