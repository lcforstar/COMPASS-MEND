function phase = partitionPrecipitation(total_cm_s, air_temperature_C, relative_humidity, parameters)
%PARTITIONPRECIPITATION  Split total precipitation into rain and snow rates.
%
%   The ClimaLand bulk SnowModel consumes prescribed P_liq and P_snow. When
%   only total precipitation is available, COMPASS applies the Jennings et al.
%   (2018) logistic phase partition used by ClimaLand forcing adapters. Call
%   after canopy interception so the rates sum exactly to ground throughfall.
%
%   Inputs:
%     total_cm_s           [in]  [cm/s]  Total water rate (>=0); usually
%                                  interception.throughfall_cms, not sky P.
%     air_temperature_C    [in]  [degC]  Air temperature (scalar).
%     relative_humidity    [in]  [-] or [%]  RH as fraction [0,1] or percent
%                                  [0,100]; normalized to [0,1] inside.
%     parameters           [in]  struct  Snow params from init.buildSnowParameters:
%                                  .precipitation_phase_model  (must be
%                                    'climaland_jennings2018')
%                                  .phase_alpha                [-] default -10.04
%                                  .phase_beta_C_inverse       [1/degC] default 1.41
%                                  .phase_gamma_rh_fraction_inverse  [1/RH_frac]
%                                    default 9.0 (= 0.09 per RH-percent)
%
%   Returns:
%     phase [out] struct
%       .liquid_cm_s               [cm/s]  Rain rate
%       .snow_cm_s                 [cm/s]  Snow rate
%       .snow_fraction             [-]     f_snow in (0,1)
%       .relative_humidity_fraction [-]    RH used in the formula [0,1]
%     Identity: liquid_cm_s + snow_cm_s == total_cm_s (within float).
%
%   Method:
%     (1) Validate scalars; map RH percent -> fraction and clamp to [0,1].
%     (2) Require precipitation_phase_model == climaland_jennings2018
%         (prescribed separate P_liq/P_snow drivers are not wired here).
%     (3) Logistic snow fraction (warmer/moister -> smaller f_snow):
%           x = alpha + beta*Ta + gamma*RH
%           f_snow = 1 / (1 + exp(x))
%           P_snow = f_snow * P_tot;  P_rain = P_tot - P_snow
%
%   Reference:
%     Jennings et al. (2018) rain/snow partition (logistic form).
%     ClimaLand ext/neural_snow/DataTools.jl::snowsplit (RH-fraction gamma=9).
%     Defaults: init.buildSnowParameters (phase_alpha/beta/gamma).
%     Caller: vegetation.integrated.advanceCanopyWaterSnow (Step 2a).
% Editor: Jizixing Zhu; Date: 2026-07-17.
% Change: 2026-07-17 Jizixing Zhu — add Inputs/Method/Reference and step comments.
% Change: 2026-07-18 Jizixing Zhu -- when the snow process is disabled,
%   preserve total precipitation as liquid without evaluating phase physics.

    % --- Step 1: validate drivers ---
    assert(isscalar(total_cm_s) && isfinite(total_cm_s) && total_cm_s >= 0, ...
        'partitionPrecipitation:badTotal', ...
        'Total precipitation must be finite and nonnegative [cm/s].');
    assert(isscalar(air_temperature_C) && isfinite(air_temperature_C), ...
        'partitionPrecipitation:badTemperature', ...
        'Air temperature must be a finite scalar [degC].');
    assert(isscalar(relative_humidity) && isfinite(relative_humidity), ...
        'partitionPrecipitation:badHumidity', ...
        'Relative humidity must be finite.');

    % --- Step 2: RH as fraction [0,1] (accept percent if > 1) ---
    rh = double(relative_humidity);
    if rh > 1 + 1e-12
        rh = rh / 100;
    end
    rh = min(max(rh, 0), 1);

    if isfield(parameters, 'enabled') && ~logical(parameters.enabled)
        phase = struct( ...
            'liquid_cm_s', double(total_cm_s), ...
            'snow_cm_s', 0, ...
            'snow_fraction', 0, ...
            'relative_humidity_fraction', rh);
        return;
    end

    % --- Step 3: only Jennings adapter is implemented ---
    model = parameters.precipitation_phase_model;
    assert(strcmp(model, 'climaland_jennings2018'), ...
        'partitionPrecipitation:prescribedNeedsDrivers', ...
        ['The prescribed phase model requires explicit liquid and snow ', ...
         'precipitation drivers.']);

    % --- Step 4: logistic f_snow; rain + snow = total ---
    % x rises with Ta and RH => exp(x) large => f_snow small (more rain).
    exponent = parameters.phase_alpha ...
        + parameters.phase_beta_C_inverse * air_temperature_C ...
        + parameters.phase_gamma_rh_fraction_inverse * rh;
    snowFraction = 1 / (1 + exp(exponent));
    snowRate = total_cm_s * snowFraction;          % [cm/s]
    rainRate = total_cm_s - snowRate;              % [cm/s], exact complement

    phase = struct( ...
        'liquid_cm_s', rainRate, ...
        'snow_cm_s', snowRate, ...
        'snow_fraction', snowFraction, ...
        'relative_humidity_fraction', rh);
end
