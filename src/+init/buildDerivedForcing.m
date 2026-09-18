function d = buildDerivedForcing(hourly, built)
%BUILDDERIVEDFORCING  Compute hourly thermodynamic derivatives from meteorology.
%
%   Inputs:
%     hourly  struct  hourly driver fields from init.buildHourlyForcing
%     built   struct  domain.site.elevation and forcing.config
%
%   Returns d with saturation vapour pressure, actual vapour pressure, humidity
%   ratios, vapour density and (virtual) potential temperature -- thermodynamic
%   derivatives only. The resolved drivers Pg/Ca/Oa live in built.forcing.hourly
%   (buildHourlyForcing); Esat reads the already-resolved hourly.Pg.
%
%   Public (F3) so a groundwater/Dunnian forcing adapter can rebuild the derived
%   layer after hourly drivers change (init.applyGroundwaterForcingAdjustment).
%
%   Reference: faithful port of legacy MeteoData.compute_derived
%   (legacy/src/+core/MeteoData.m:495-557), reusing equations.Calc_Esat_Flatau and
%   src-level fixed physical constants.
%
%   Editor: Jizixing Zhu; Date: 2026-06-26.
%   Change: 2026-06-27 Jizixing Zhu -- derived holds only thermodynamic
%     derivatives; Pg/Ca/Oa moved to buildHourly (resolved drivers). Esat reads
%     the already-resolved hourly.Pg.
%   Change: 2026-07-01 Jizixing Zhu -- extracted verbatim from buildForcing nested
%     buildDerived to a public init.buildDerivedForcing (rebuildable by the adapter).
    fcfg    = getForcingConfig(built);
    z_meas  = fcfg.Z_measure;                % [m] reference height (potential temp)

    phys    = equations.physicalConstants();
    Rv      = phys.R_wv_SI;   % water-vapour gas constant [J/(kg K)]
    Gamma_d = phys.Gamma_d;   % dry adiabatic lapse rate [K/m]
    eps_wv  = phys.epsilon;   % MW_wv/MW_da ~ 0.622

    Ta = hourly.Ta;
    Pg = hourly.Pg;           % resolved driver (measured | barometric) from buildHourlyForcing

    % saturation vapour pressure + derivative + saturation specific humidity (Flatau)
    [es_Pa, des_dT, qs_sat, ~] = equations.Calc_Esat_Flatau(Ta, Pg);
    d.es    = es_Pa;        % [Pa] saturation vapour pressure
    d.delta = des_dT;       % [Pa/K] des/dT (Penman-Monteith Delta)
    d.qs    = qs_sat;       % [kg/kg] saturation specific humidity

    ea_Pa  = es_Pa .* hourly.RH;          % [Pa] actual vapour pressure
    d.ea   = ea_Pa;
    d.qa   = (eps_wv .* ea_Pa) ./ (Pg - (1 - eps_wv) .* ea_Pa);   % [kg/kg]
    d.RHOV_A = ea_Pa ./ (Rv .* (Ta + phys.C_to_K));               % vapour density

    TaK       = Ta + phys.C_to_K;
    d.Theta   = TaK + Gamma_d * z_meas;                % [K] potential temperature
    d.Theta_v = d.Theta .* (1 + 0.61 .* d.qa);         % [K] virtual potential temp
    % Pg/Ca/Oa are resolved DRIVERS (built.forcing.hourly), not re-stored here.
end
