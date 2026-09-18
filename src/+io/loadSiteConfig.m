function config = loadSiteConfig(siteDir)
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
end
