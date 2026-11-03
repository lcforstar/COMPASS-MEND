function raw = loadSite(siteDir, plan, config)
%LOADSITE  Stage-1 assembler: read a self-contained site folder into raw containers.
%
%   raw = io.loadSite(siteDir) reads the legacy full site package. The optional
%   raw = io.loadSite(siteDir, plan, config) form reads only inputs enabled by
%   plan.inputs; CONFIG is the result of io.loadSiteConfig and avoids parsing
%   site.toml twice.
%   single struct of RAW, format-agnostic data at the files' native granularity.
%   This is stage 1 (pure read): it performs NO grid construction, NO
%   interpolation, NO van-Genuchten inversion, NO config-name -> code-name
%   mapping and NO unit conversion -- all of that is stage 2 (init.buildState).
%   Keeping read separate from transform lets gate A compare raw containers
%   (read fidelity) independently of the interpolated grid state.
%
%   Inputs:
%     siteDir  [char]  path to a site folder, e.g. data/sites/changwu
%
%   Returns (fields absent if the file is missing, except config which is required):
%     raw.config   struct  parsed site.toml (nested: site/time/grid/switches/...)
%     raw.soil     table   per-depth soil properties        (soil.csv; horizon
%                          boundaries or samples selected in site.toml)
%     raw.initial  table   per-depth soil initial conditions (soil_initial_profile.csv)
%     raw.organic_initial table organic-pool fractions (organic_initial_pools.csv)
%     raw.plant    table   plant C/N pool spin-up state      (plant_initial_pools.csv)
%     raw.root     table   root-length-density profile       (root.csv)
%     raw.forcing  table   meteorological time series         (forcing.csv)
%     raw.water_boundary table optional prescribed water BC   (configured CSV)
%     raw.water_boundary_path char resolved prescribed-BC file path
%     raw.vegetation table optional measured canopy series    (canopy_measurements.csv)
%     raw.siteDir  char    the source folder (provenance)
%
%   Files are read by the format adapters io.adapters.readToml / readCsvTable.
%
%   Editor: Jizixing Zhu; Date: 2026-06-26.
%   Change: 2026-07-11 Jizixing Zhu -- read optional water_boundary.csv as raw prescribed BC data.
%   Change: 2026-07-11 Jizixing Zhu -- resolve boundary.water_boundary_file explicitly and expose its path as provenance.
%   Change: 2026-07-11 Jizixing Zhu -- add plan-aware reads so inactive laboratory component files are never touched.
%   Change: 2026-07-12 Jizixing Zhu -- document the two public soil.csv depth semantics.
%   Change: 2026-07-16 Jizixing Zhu -- replace ambiguous legacy input filenames with descriptive site-interface names.
%   Change: 2026-08-16 Codex -- read organic initialization only for an active cycle.
%   Change: 2026-08-28 Codex -- keep organic management in site.toml beside
%     mineral fertilization; no separate management CSV is read.

if ~isfolder(siteDir)
    error('loadSite:noDir', 'Site folder not found: %s', siteDir);
end

    raw = struct();
    raw.siteDir = siteDir;

    if nargin < 2 || isempty(plan)
        plan = fullInputPlan();
    end
    if nargin < 3 || isempty(config)
        config = io.loadSiteConfig(siteDir);
    end
    raw.config = config;
    raw.input_plan = plan;

    % --- CSV tables: one read each, into one raw field (absent file -> []).
    %     buildState/checkSite decide which are mandatory. ---
    csvTables = { ...
        'soil',    'soil.csv'; ...
        'initial', 'soil_initial_profile.csv'; ...
        'organic_initial', 'organic_initial_pools.csv'; ...
        'plant',   'plant_initial_pools.csv'; ...
        'root',           'root.csv'; ...
        'forcing',        'forcing.csv'};
    for k = 1:size(csvTables, 1)
        field = csvTables{k, 1};
        path = fullfile(siteDir, csvTables{k, 2});
        if inputEnabled(plan, field) && exist(path, 'file') == 2
            raw.(csvTables{k, 1}) = io.adapters.readCsvTable(path);
        else
            raw.(csvTables{k, 1}) = [];
        end
    end

    waterBoundaryFile = 'water_boundary.csv';
    if isfield(raw.config, 'boundary') ...
            && isfield(raw.config.boundary, 'water_boundary_file')
        waterBoundaryFile = readConfiguredPath( ...
            raw.config.boundary.water_boundary_file, ...
            'boundary.water_boundary_file');
    end
    raw.water_boundary_path = resolveSitePath(siteDir, waterBoundaryFile);
    if inputEnabled(plan, 'water_boundary') && exist(raw.water_boundary_path, 'file') == 2
        raw.water_boundary = io.adapters.readCsvTable(raw.water_boundary_path);
    else
        raw.water_boundary = [];
    end

    raw.vegetation = [];
    if inputEnabled(plan, 'vegetation')
        raw.vegetation = readOptionalVegetationTable(raw.config, siteDir);
    end
end

% ---------------------------------------------------------------------------
function plan = fullInputPlan()
% Backward-compatible direct-load policy for callers without an init plan.
% Editor: Jizixing Zhu; Date: 2026-07-11.
    plan.inputs = struct('soil', true, 'initial', true, 'plant', true, ...
        'root', true, 'forcing', true, 'water_boundary', true, ...
        'vegetation', true, 'scope', true, 'organic_initial', false);
end

% ---------------------------------------------------------------------------
function tf = inputEnabled(plan, name)
% Read one input-plan flag; absent flags default to the legacy full read.
% Editor: Jizixing Zhu; Date: 2026-07-11.
    tf = true;
    if isfield(plan, 'inputs') && isfield(plan.inputs, name)
        tf = logical(plan.inputs.(name));
    end
end

% ---------------------------------------------------------------------------
function value = readConfiguredPath(raw, label)
% Read one non-empty scalar text path from site.toml.
% Editor: Jizixing Zhu; Date: 2026-07-11.
    if ~(ischar(raw) || (isstring(raw) && isscalar(raw)))
        error('loadSite:badConfiguredPath', '%s must be a scalar text path.', label);
    end
    value = strtrim(char(raw));
    if isempty(value)
        error('loadSite:badConfiguredPath', '%s must not be empty.', label);
    end
end

% ---------------------------------------------------------------------------
function tbl = readOptionalVegetationTable(config, siteDir)
% Read optional measured vegetation structure data for LAI/hc/SAI data modes.
% Editor: Jizixing Zhu; Date: 2026-06-28.
    tbl = [];
    explicitFile = '';
    if isfield(config, 'vegetation') && isfield(config.vegetation, 'vegetation_data_file')
        explicitFile = strtrim(char(config.vegetation.vegetation_data_file));
    end

    if ~isempty(explicitFile)
        tbl = io.adapters.readTableFile(resolveSitePath(siteDir, explicitFile));
        return;
    end

    candidates = {'canopy_measurements.csv'};
    for i = 1:numel(candidates)
        path = fullfile(siteDir, candidates{i});
        if exist(path, 'file') == 2
            tbl = io.adapters.readTableFile(path);
            return;
        end
    end
end

% ---------------------------------------------------------------------------
function path = resolveSitePath(siteDir, value)
% Resolve a site-relative path unless the input is already absolute.
% Editor: Jizixing Zhu; Date: 2026-06-28.
    value = char(value);
    if ~isempty(regexp(value, '^[A-Za-z]:[\\/]', 'once')) || ...
            startsWith(value, '\\') || startsWith(value, '/')
        path = value;
    else
        path = fullfile(siteDir, value);
    end
end
