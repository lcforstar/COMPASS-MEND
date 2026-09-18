function groups = buildEnv(rawConfig)
%BUILDENV  Split [environment] scalars into their physical owner containers.
%
%   groups = init.buildEnv(rawConfig) maps friendly site.toml [environment]
%   names into grouped scalar settings. It is a pure config translation layer:
%   no forcing time series, no interpolation, no physics.
%
%   Physical ownership:
%     groups.forcing       Z_measure/Wind_*/Rin_unit_is_watt + Ca/Oa
%     groups.boundary      Tav/At for the bottom soil-temperature climatology
%     groups.scope         dry/wet soil thermal-infrared reflectance endpoints
%
%   Inputs:
%     rawConfig  struct  raw.config (uses .site and .environment)
%
%   Returns groups with fields:
%     forcing     struct  driver/background-atmosphere scalars
%     boundary    struct  thermal-bottom-boundary climatology scalars
%     scope       struct  soil optical thermal-reflectance endpoints
%
%   Reference: legacy EnvironmentSettings (field names).
%
%   Editor: Jizixing Zhu; Date: 2026-06-26.

    e    = rawConfig.environment;
    groups = struct();

    forcing = struct();
    forcing.Rin_unit_is_watt    = double(logical(e.rin_unit_is_watt));
    forcing.Z_measure           = e.z_measure;
    forcing.Wind_Height_Convert = double(logical(e.wind_convert));
    forcing.Wind_Obs_Height     = e.wind_obs_height;
    forcing.Wind_z0_grass       = e.wind_z0_grass;
    forcing.Wind_d_grass        = e.wind_d_grass;
    forcing.Ca                  = e.Ca;
    if isfield(e, 'Oa'); forcing.Oa = e.Oa; else; forcing.Oa = NaN; end

    boundary = struct();
    boundary.thermal_bottom = struct('Tav', e.Tav, 'At', e.At);

    scope = struct();
    scope.soil_optical = struct( ...
        'Rho_Thermal_Soil_Dry', e.rho_thermal_soil_dry, ...
        'Rho_Thermal_Soil_Wet', e.rho_thermal_soil_wet);

    groups.forcing    = forcing;
    groups.boundary   = boundary;
    groups.scope      = scope;
end
