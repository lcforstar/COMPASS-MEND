function model = run_site(siteSelector, options)
%RUN_SITE  Select and run one configured COMPASS site from the repository root.
%
%   model = run_site("changwu") runs data/sites/changwu as a point-field task.
%   model = run_site("lab_water") runs data/sites/lab_water as a standalone
%   laboratory Richards task. SITESELECTOR may also be a site-folder path.
%   The [simulation] block in site.toml selects the calculation path; this
%   launcher does not infer context from which files happen to exist.
%
%   options are passed to runCOMPASS or runWaterOnlySite. For either context,
%   options.stop_after_initialization=true loads and validates the selected
%   site without advancing time.
%
%   Editor: Jizixing Zhu; Date: 2026-07-11.

    if nargin < 1 || isempty(siteSelector), siteSelector = "changwu"; end
    if nargin < 2 || isempty(options), options = struct(); end
    rootDir = fileparts(mfilename('fullpath'));
    addpath(fullfile(rootDir, 'src'));

    candidate = char(string(siteSelector));
    if isfolder(candidate)
        siteDir = candidate;
    else
        siteDir = fullfile(rootDir, 'data', 'sites', candidate);
    end
    config = io.loadSiteConfig(siteDir);
    plan = init.resolveSimulationPlan(config);

    if strcmp(plan.spatial_mode, 'regional')
        if ~isfield(options, 'regional_input') || isempty(options.regional_input)
            error('run_site:missingRegionalInput', ...
                'A regional site requires options.regional_input.');
        end
        model = runCOMPASS(siteDir, options);
    elseif plan.is_laboratory
        model = runWaterOnlySite(siteDir, options);
    else
        model = runCOMPASS(siteDir, options);
    end
end
