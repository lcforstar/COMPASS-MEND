function rbs = calcSoilBoundaryResistanceKustas(Ts_C, Ta_C, u_soil_m_s, params)
%CALCSOILBOUNDARYRESISTANCEKUSTAS  Soil boundary-layer resistance for TSEB-style ebal.
%
%   Inputs:
%     Ts_C        [degC]       Soil surface temperature. Vectors are averaged to one
%                              bulk soil-surface temperature because COMPASS carries one rbs.
%     Ta_C        [degC]       Air temperature.
%     u_soil_m_s  [m s-1]      Wind speed at the soil boundary layer.
%     params      struct       Optional: b [-], c [m s-1 K^(-1/3)],
%                              rbs_min [s m-1], rbs_max [s m-1].
%
%   Output:
%     rbs         [s m-1]      Soil boundary-layer resistance.
%
%   Method:
%     R_s = 1 / (c * max(Ts - Ta, 0)^(1/3) + b * u_soil).
%     A temperature difference in degC has the same increment as K, so deltaT_K
%     is computed directly from Ts_C - Ta_C and labelled as [K].
%     u_soil is passed from resistances_clm45. That module already applies
%     COMPASS wind floors: forcing wind U >= 1 m s-1 at MeteoData input,
%     ustar >= 0.001 m s-1, canopy-top wind uh >= 0.01 m s-1, and then
%     computes the attenuated near-soil wind uz0.
%
%   Reference:
%     Kustas and Norman 1999 Agric. For. Meteorol. 94:13-29,
%     DOI:10.1016/S0168-1923(99)00005-2. Formula follows pyTSEB v2.1
%     resistances.calc_R_S_Kustas. Defaults follow pyTSEB source constants:
%     KN_b = 0.012 [-] (commented there as Kustas et al. 1999) and
%     KN_c = 0.0038 [m s-1 K^(-1/3)] (commented there as Kustas et al. 2016).

    if nargin < 4 || isempty(params)
        params = struct();
    end

    b = get_param(params, 'b', 0.012);           % [-]
    c = get_param(params, 'c', 0.0038);          % [m s-1 K^(-1/3)]
    rbs_min = get_param(params, 'rbs_min', 0);   % [s m-1]
    rbs_max = get_param(params, 'rbs_max', Inf); % [s m-1]

    Ts_eff_C = finite_mean(Ts_C);
    Ta_eff_C = finite_mean(Ta_C);
    if ~isfinite(Ta_eff_C)
        Ta_eff_C = 0;
    end
    if ~isfinite(Ts_eff_C)
        Ts_eff_C = Ta_eff_C;
    end

    deltaT_K = max(Ts_eff_C - Ta_eff_C, 0);                  % [K]
    u_soil_eff = max(finite_mean(u_soil_m_s), 0);             % [m s-1]

    denom_m_s = c .* deltaT_K.^(1 / 3) + b .* u_soil_eff;     % [m s-1]
    denom_m_s = max(denom_m_s, realmin('double'));
    rbs = 1 ./ denom_m_s;                                     % [s m-1]
    rbs = min(max(rbs, rbs_min), rbs_max);

end

function v = get_param(params, name, defaultValue)
    if isstruct(params) && isfield(params, name)
        candidate = params.(name);
    else
        candidate = [];
    end
    if isnumeric(candidate) && isscalar(candidate) && isfinite(candidate)
        v = candidate;
    else
        v = defaultValue;
    end
end

function m = finite_mean(x)
    x = double(x(:));
    x = x(isfinite(x));
    if isempty(x)
        m = NaN;
    else
        m = mean(x);
    end
end
