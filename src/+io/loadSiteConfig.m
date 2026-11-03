function config = loadSiteConfig(siteDir)
% Editor: Wenjie Liu; Date: 2026-10-17.
% Change: 2026-10-17 Wenjie Liu -- default and validate the optional surface-resistance multiplier.
%LOADSITECONFIG  Read only the site-level TOML configuration.
%
%   config = io.loadSiteConfig(siteDir) resolves site.toml inside one site
%   package and parses it without reading any CSV or model-library input. This
%   first pass lets initialization decide which component inputs are active
%   before touching optional files.
%
%   Editor: Jizixing Zhu; Date: 2026-07-11.

    if ~isfolder(siteDir)
        error('loadSiteConfig:noDir', 'Site folder not found: %s', siteDir);
    end
    path = fullfile(siteDir, 'site.toml');
    if exist(path, 'file') ~= 2
        error('loadSiteConfig:noConfig', 'Required site.toml not found in %s', siteDir);
    end
    config = io.adapters.readToml(path);
    if ~isfield(config,'scope'),config.scope=struct();end
    if ~isfield(config.scope,'soil_optical'),config.scope.soil_optical=struct();end
    if ~isfield(config.scope.soil_optical,'rss_scale')
        config.scope.soil_optical.rss_scale=1.0;
    end
    validateattributes(config.scope.soil_optical.rss_scale,{'numeric'}, ...
        {'real','scalar','positive','finite'},mfilename,'scope.soil_optical.rss_scale');
end
