function scopeState = buildScopeState(built)
%BUILDSCOPESTATE  Prognostic canopy-air + SCOPE soil-heat initial state (built.Y.*).
%
%   scopeState = init.buildScopeState(built) builds the cross-step SCOPE/canopy
%   states the layer-4 canopy energy balance carries between outer steps, so they
%   are accepted state in built.Y (checkpoint/resume safe), NOT per-step scratch:
%
%     scopeState.canopy_air  struct  Ta_layer/ea_layer/CO2_layer [nlayers x 1] --
%                            the canopy-air profile OLD layers used as the old time
%                            level in the implicit canopy-air scalar balances
%                            (ebal_clm45_profile). Initialised to the atmospheric
%                            value broadcast over the canopy layers, which is exactly
%                            ebal's own first-step guess; the EB solve overwrites them
%                            each accepted step (profiles.*_layer -> next step's OLD).
%     scopeState.scope_soil  struct  Tsold [12 x 2] -- the SCOPE force-restore
%                            soil-heat temperature history (12 = force-restore history
%                            depth, 2 = sunlit/shaded soil faces). DORMANT under the
%                            default soil_heat_method = 2; only exercised when
%                            soil_heat_method < 2.
%
%   Inputs:
%     built  struct  needs parameters.scope.canopy_rt.nlayers + forcing.hourly.Ta/Ca
%                    + forcing.derived.ea (built before Y in init.buildState).
%
%   Reference: ebal_clm45_profile.m:137-139 (Ta/ea/CO2_layer = atm * ones(nl,1)),
%   run_stemmus_initialize.m:402-404 (Tsold = Ta * ones(12,2)).
%   Editor: Jizixing Zhu; Date: 2026-07-02.
%   Change: 2026-07-18 Jizixing Zhu -- initialize force-restore skin-temperature
%           history from the supplied top-soil temperature instead of air.
    nl = built.parameters.scope.canopy_rt.nlayers;   % [-] canopy layer count (30)
    h  = built.forcing.hourly;
    d  = built.forcing.derived;

    Ta = firstFinite(h.Ta, 20);              % [degC]      first-step air temperature
    ea = firstFinite(d.ea, 1000) / 100;      % [hPa]       first-step vapour pressure (Pa -> hPa, SCOPE unit)
    Ca = firstFinite(h.Ca, 400);             % [umol/mol]  first-step CO2 mixing ratio

    scopeState = struct();
    scopeState.canopy_air = struct( ...
        'Ta_layer',  Ta * ones(nl, 1), ...    % [degC]
        'ea_layer',  ea * ones(nl, 1), ...    % [hPa]
        'CO2_layer', Ca * ones(nl, 1));       % [umol/mol]
    initialSkinTemperature = firstFinite(built.Y.T(end), Ta); % [degC]
    scopeState.scope_soil = struct( ...
        'Tsold', initialSkinTemperature * ones(12, 2)); % [degC, 12 x 2]
end

% ---------------------------------------------------------------------------
function v = firstFinite(x, fallback)
% First finite element of a forcing series, or a fallback when absent/non-finite.
% Editor: Jizixing Zhu; Date: 2026-07-02.
    v = fallback;
    if ~isempty(x) && isfinite(x(1))
        v = x(1);
    end
end
