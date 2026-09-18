function model = initializeRegional(siteDir, regionalInput, options)
%INITIALIZEREGIONAL  Build a regional model from a space struct or NetCDF static file.
%
%   model = init.initializeRegional(siteDir, regionalInput, options) is the
%   regional/global initialization entry point. It reuses the validated site
%   package as a column template, expands it onto the supplied horizontal cells,
%   optionally applies NetCDF static overrides, and allocates a matching
%   multi-cell workspace.
%
%   Inputs:
%     siteDir        [char] path to a site folder, e.g. 'data/sites/changwu'
%     regionalInput  [struct|char] either:
%       (1) horizontal cell descriptor:
%       .lat      [nCell x 1] latitude  [degree_north]
%       .lon      [nCell x 1] longitude [degree_east]
%       .mask     [nCell x 1] logical land mask (optional, default true)
%       .cell_id  [nCell x 1] cell identifiers (optional, default 1:nCell)
%       (2) NetCDF static file readable by init.readRegionalStatic.
%     options       [struct] optional NetCDF block selection:
%       .active_index [nBlock x 1] full-grid linear land indices
%       .max_cells    [scalar]     first N active land cells
%
%   Returns model with model.built.domain.kind = 'regional' and nCell-aware
%   Y/parameters/forcing/workspace shapes. This function is solver-neutral:
%   it prepares initial arrays and the pre-step soil auxiliary cache, but never
%   advances one time step.
% Editor: Jizixing Zhu; Date: 2026-06-27.

    if nargin < 2
        error('initializeRegional:badSpace', ...
            'initializeRegional requires a space struct or NetCDF static path.');
    end
    if nargin < 3 || isempty(options)
        options = struct();
    end

    regional = struct();
    if ischar(regionalInput) || isstring(regionalInput)
        regional = init.readRegionalStatic(regionalInput, options);
        space = regional.space;
    elseif isstruct(regionalInput) && isfield(regionalInput, 'space') ...
            && isfield(regionalInput, 'overrides')
        regional = regionalInput;
        space = regional.space;
    elseif isstruct(regionalInput)
        space = regionalInput;
    else
        error('initializeRegional:badSpace', ...
            'regionalInput must be a space struct, regional struct, or NetCDF static path.');
    end

    raw    = io.loadSite(siteDir);
    report = io.checkSite(raw);
    built  = init.buildState(raw);
    built  = init.expandBuiltToRegional(built, space);
    built  = init.applyRegionalStatic(built, regional);
    built  = init.allocateWorkspace(built);
    built  = init.initializeSoilAux(built);

    model              = struct();
    model.siteDir      = raw.siteDir;
    model.raw          = raw;
    model.check_report = report;
    model.built        = built;
    model.status       = "initialized_regional";
    model.created_at   = datetime("now");
end
