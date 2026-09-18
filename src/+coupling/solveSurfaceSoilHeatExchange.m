function exchange = solveSurfaceSoilHeatExchange(surface_response, soil_response)
%SOLVESURFACESOILHEATEXCHANGE Solve the frozen surface--soil heat interface.
%
%   exchange = coupling.solveSurfaceSoilHeatExchange(surface_response,
%   soil_response) intersects two affine maps at one outer-step interface:
%
%       q_atm(Ti)  = q_ref + a (Ti - T_ref),
%       q_soil(Ti) = K (Ti - T_zero).
%
%   Inputs:
%     surface_response [struct] Frozen post-ebal atmospheric response, with
%                      inward-positive flux into the soil [W m-2].
%     soil_response    [struct] Upper-soil implicit thermal response [W m-2 K-1].
%
%   Returns:
%     exchange [struct] One interface temperature and one common energy flux.
%
%   Physical contract:
%     Ti is an algebraic interface temperature, not a second soil prognostic
%     node and not the diagnostic radiation skin temperature.  The returned
%     flux is the only requested surface-to-soil energy boundary for the
%     following coupled h/U/P_g solve.  `surface.energy.match...` owns the
%     algebra so legacy and new callers use identical mathematics.
%
%   Editor: Jizixing Zhu; Date: 2026-07-20.
    match = surface.energy.matchDynamicSoilHeatResponse( ...
        surface_response, soil_response);
    exchange = struct( ...
        'schema_version', 1, ...
        'producer', 'coupling.surface_soil_heat', ...
        'receiver', 'soil', ...
        'sign_convention', 'flux_into_soil_positive', ...
        'interface_temperature_C', match.interface_temperature_C, ...
        'flux_into_soil_W_m2', match.matched_flux_into_soil_W_m2, ...
        'surface_flux_W_m2', match.surface_flux_at_match_W_m2, ...
        'linear_residual_W_m2', match.linear_closure_residual_W_m2, ...
        'surface_response', surface_response, ...
        'soil_response', soil_response, ...
        'match', match);
end
