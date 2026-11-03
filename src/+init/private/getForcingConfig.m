function fcfg = getForcingConfig(built)
%GETFORCINGCONFIG  Return forcing-owned static driver configuration.
%
%   Inputs:
%     built  struct  initialized build container with built.forcing.config
%
%   Returns:
%     fcfg   struct  wind/radiation/background-atmosphere scalar settings
%
%   The new initialization path has no built.env fallback; [environment] must be
%   split by init.buildEnv before buildForcing is called. Kept in +init/private so
%   buildForcing, buildHourlyForcing and buildDerivedForcing share one copy.
%
%   Editor: Jizixing Zhu; Date: 2026-06-26.
%   Change: 2026-07-01 Jizixing Zhu — moved from buildForcing nested to +init/private (shared by the split public F2/F3 builders).
    if isfield(built, 'forcing') && isfield(built.forcing, 'config')
        fcfg = built.forcing.config;
    else
        error('buildForcing:missingForcingConfig', ...
            'built.forcing.config is required; split [environment] with init.buildEnv first.');
    end
end
