function U_target = convertWindHeight(U_obs, forcingConfig)
%CONVERTWINDHEIGHT  Convert observed wind speed to the model reference height.
%
%   U_target = convertWindHeight(U_obs, forcingConfig) applies a neutral
%   logarithmic wind profile to convert wind speed from Wind_Obs_Height to
%   Z_measure. It is the initialization-layer replacement for the legacy
%   legacy MeteoData.convert_wind_height static method.
%
%   Inputs:
%     U_obs          [m/s] wind speed at Wind_Obs_Height
%     forcingConfig  struct with fields Wind_Height_Convert, Wind_Obs_Height,
%                    Z_measure, Wind_z0_grass and Wind_d_grass
%
%   Returns:
%     U_target       [m/s] wind speed at Z_measure; unchanged when conversion is
%                    disabled or the geometry is invalid.
%
%   Method:
%     U_target / U_obs = log((z_target - d) / z0) / log((z_obs - d) / z0)
%
%   Editor: Jizixing Zhu; Date: 2026-06-26.

    U_target = U_obs;

    doConv  = getConfigField(forcingConfig, 'Wind_Height_Convert', 0);
    if doConv ~= 1
        return;
    end

    z_obs    = getConfigField(forcingConfig, 'Wind_Obs_Height', 2.0);
    z_target = getConfigField(forcingConfig, 'Z_measure', 10.0);
    z0       = getConfigField(forcingConfig, 'Wind_z0_grass', 0.03);
    d0       = getConfigField(forcingConfig, 'Wind_d_grass', 0.0);

    if z_obs <= d0 + z0
        warning('convertWindHeight:invalidGeometry', ...
            'Wind_Obs_Height=%.2f must be greater than Wind_d_grass+Wind_z0_grass=%.3f; skipping conversion.', ...
            z_obs, d0 + z0);
        return;
    end
    if z_target <= d0 + z0
        warning('convertWindHeight:invalidGeometry', ...
            'Z_measure=%.2f must be greater than Wind_d_grass+Wind_z0_grass=%.3f; skipping conversion.', ...
            z_target, d0 + z0);
        return;
    end

    if abs(z_target - z_obs) < 1e-6
        return;
    end

    ratio = log((z_target - d0) / z0) / log((z_obs - d0) / z0);
    U_target = U_obs * ratio;

    persistent logged
    if isempty(logged)
        logged = false;
    end
    if ~logged
        fprintf('[buildForcing] Wind height conversion: U(%.1fm) -> U(%.1fm), ratio=%.4f (z0=%.3fm, d=%.3fm)\n', ...
            z_obs, z_target, ratio, z0, d0);
        logged = true;
    end
end

function value = getConfigField(s, name, defaultValue)
%GETCONFIGFIELD  Read a scalar config field with a finite fallback.
%
%   Inputs:
%     s             struct config container
%     name          char field name
%     defaultValue  scalar fallback value
%
%   Returns:
%     value         field value, or defaultValue when absent/invalid
%
%   Editor: Jizixing Zhu; Date: 2026-06-26.

    if isstruct(s) && isfield(s, name)
        value = s.(name);
    else
        value = defaultValue;
        return;
    end

    if isempty(value) || any(~isfinite(value(:)))
        value = defaultValue;
    end
end
