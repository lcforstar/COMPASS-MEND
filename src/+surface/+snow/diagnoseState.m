function d = diagnoseState(state, parameters)
%DIAGNOSESTATE  Diagnose bulk-snow thermodynamics and geometry from S/S_l/U.
%
%   State follows ClimaLand's bulk SnowModel, expressed in COMPASS units:
%     SWE          [cm]       total water equivalent S
%     liquid_water [cm]       retained liquid water S_l
%     energy       [J cm-2]   total snow internal energy U
%
%   The zero of energy is liquid water at 0 degC. Ice therefore contributes
%   negative latent energy. Air heat capacity is neglected, as in ClimaLand.
%
%   Source model: ClimaLand SnowModel and snow_parameterizations.jl,
%   Apache License 2.0, Copyright 2022 California Institute of Technology.
%   MATLAB adaptation changes names, units, and discrete coupling.
%
%   Editor: Jizixing Zhu; Date: 2026-07-17.

    S = max(finiteScalar(state, 'SWE', 0), 0);
    Sl = min(max(finiteScalar(state, 'liquid_water', 0), 0), S);
    Ucm = finiteScalar(state, 'energy', 0);
    Sm = S / 100;
    Slm = Sl / 100;
    U = Ucm * 1e4;
    rhoL = parameters.rho_liquid_kg_m3;
    latent = parameters.latent_heat_fusion_J_kg;
    cpIce = parameters.specific_heat_ice_J_kg_K;
    cpLiquid = parameters.specific_heat_liquid_J_kg_K;

    if S <= 1e-14
        ql = 0;
        temperature = 0;
    else
        ql = Slm / Sm;
        cp = ql * cpLiquid + (1 - ql) * cpIce;
        temperatureK = parameters.reference_temperature_K ...
            + (U + rhoL * latent * Sm * (1 - ql)) ...
            / max(rhoL * cp * (Sm + parameters.regularization_swe_m), eps);
    end
    if S <= 1e-14, temperatureK = parameters.reference_temperature_K; end
    density = parameters.minimum_density_kg_m3 * (1 - ql) ...
        + rhoL * ql;
    depth = Sm * rhoL / max(density, eps); % [m]
    x = density / parameters.rho_ice_kg_m3;
    kAir = parameters.thermal_conductivity_air_W_m_K;
    conductivity = kAir ...
        + (0.07 * x + 0.93 * x^2) ...
        * (parameters.thermal_conductivity_ice_W_m_K - kAir);
    scaledDepth = depth / parameters.cover_depth_scale_m;
    cover = min(parameters.cover_fraction_beta * scaledDepth ...
        / (scaledDepth + 1), 1);
    if S <= 1e-14, cover = 0; end

    d = struct( ...
        'SWE_cm', S, 'liquid_water_cm', Sl, ...
        'ice_water_equivalent_cm', S - Sl, ...
        'energy_J_cm2', Ucm, 'liquid_fraction', ql, ...
        'temperature_K', temperatureK, ...
        'temperature_C', temperatureK - parameters.freezing_temperature_K, ...
        'density_kg_m3', density, 'depth_m', depth, ...
        'thermal_conductivity_W_m_K', conductivity, ...
        'cover_fraction', cover, ...
        'albedo', parameters.albedo, ...
        'emissivity', parameters.emissivity);
end

function value = finiteScalar(state, name, fallback)
% Read a finite scalar state field.
% Editor: Jizixing Zhu; Date: 2026-07-17.
    value = fallback;
    if isstruct(state) && isfield(state, name) && ~isempty(state.(name))
        value = double(state.(name));
    end
    assert(isscalar(value) && isfinite(value), ...
        'diagnoseSnowState:badState', 'snow.%s must be finite scalar.', name);
end
