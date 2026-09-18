function response = computeGroundHeatResponse(built, snowState, parameters)
%COMPUTEGROUNDHEATRESPONSE  Linear snow-to-soil conduction at the surface node.
%
%   COMPASS places the upper soil node at z=0, the physical snow-soil
%   interface.  The bulk snow temperature represents the snow-layer mean, so
%   the conductive flux into soil is
%
%     G_si(T_s) = C_si (T_snow - T_s),
%     C_si      = k_snow / (0.5 max(d_snow,d_min)).
%
%   The returned derivative is with respect to the soil surface temperature.
%   This affine response is inserted directly in the U-primary top row; it is
%   not frozen from the accepted soil temperature and it does not add a second
%   artificial soil skin resistance.
%
%   Reference: CLM5 SoilTemperatureMod combines snow and soil conduction in
%   one implicit thermal system.  This reduced one-snow-node response is the
%   matching Schur form for COMPASS's explicit z=0 soil node.
%
%   Editor: Jizixing Zhu; Date: 2026-07-19.

    d = surface.snow.diagnoseState(snowState, parameters);
    assert(~isempty(built.Y.T) && isfinite(built.Y.T(end)), ...
        'computeGroundHeatResponse:missingSurfaceTemperature', ...
        'The accepted z=0 soil-surface temperature is required.');
    referenceTemperature = double(built.Y.T(end));
    if isfield(built.workspace, 'surface') ...
            && isfield(built.workspace.surface, ...
                'soil_surface_temperature_predictor_C') ...
            && isfinite(built.workspace.surface. ...
                soil_surface_temperature_predictor_C)
        referenceTemperature = double(built.workspace.surface. ...
            soil_surface_temperature_predictor_C);
    end
    if d.cover_fraction <= 0 || d.ice_water_equivalent_cm <= 0
        conductance = 0;
        flux = 0;
    else
        exchangeDepth = max(d.depth_m, parameters.minimum_exchange_depth_m);
        conductance = 2 * d.thermal_conductivity_W_m_K ...
            / max(exchangeDepth, eps);
        flux = conductance * (d.temperature_C - referenceTemperature);
    end
    response = struct( ...
        'reference_soil_surface_temperature_C', referenceTemperature, ...
        'snow_temperature_C', d.temperature_C, ...
        'conductance_W_m2_K', conductance, ...
        'reference_flux_into_soil_W_m2', flux, ...
        'dflux_dsoil_surface_temperature_W_m2_K', -conductance, ...
        'cover_fraction', d.cover_fraction, ...
        'snow_state', d);
end
