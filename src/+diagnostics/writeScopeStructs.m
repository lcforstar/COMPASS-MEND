function writeScopeStructs(reg, fluxes, thermal, profiles, rad, KT, cfg)
%WRITESCOPESTRUCTS  Flatten the four SCOPE diagnostic structs (fluxes/thermal/
%profiles/rad) into the diagnostics registry.
%
% Replaces the legacy helpers.append_numeric_snapshot pattern that flattened
% these structs into Scope_output.mat. Each registered name in the radiation /
% energy_balance / carbon_fluxes / canopy_aerodynamics / canopy_profiles
% groups has a SCOPE source field (e.g. fluxes.Rntot) — this function
% materialises that mapping at write time. Missing source fields are silently
% skipped so the call works for SCOPE variants that produce a subset.
%
% Inputs:
%   reg       [in,inout]  diagnostics.Registry
%   fluxes    [in]        SCOPE fluxes struct (may be empty)
%   thermal   [in]        SCOPE thermal struct (may be empty)
%   profiles  [in]        SCOPE profiles struct (canopy-layer fields)
%   rad       [in]        SCOPE rad struct
%   KT        [in] [-]    Outer time step index
%   cfg       [in]        config/switch struct; respiration is emitted only
%                         when the full soil-N + PlantCN chain is active
%
% Author: Zhu Jizixing, 2026-05-15
% Change: 2026-07-18 Jizixing Zhu -- convert mol PAR sources to the public
%   umol contract and keep ecosystem respiration missing when C/N is off.
% Change: 2026-07-18 Jizixing Zhu -- export explicit bare/snow radiative,
%   turbulent, vapor, conductive, and temperature partitions.
% Change: 2026-07-19 Jizixing Zhu -- write the physical soil-top response
%   separately from the diagnostic SCOPE force-restore heat-flux predictor.
% Change: 2026-07-20 Jizixing Zhu -- materialize the RTMo spectral integrals that
%   are actually produced on accepted steps; retired Lout/Loutt/Eout scalar
%   placeholders are never used as physical output sources.

    % ----- radiation (12) -----
    tryWrite(reg, 'net_radiation_total',       fluxes, 'Rntot',     KT);
    tryWrite(reg, 'net_radiation_canopy',      fluxes, 'Rnctot',    KT);
    tryWrite(reg, 'net_radiation_soil',        fluxes, 'Rnstot',    KT);
    % RTMo provides integrated hemispherical fluxes as Eouto (optical),
    % Eoutt (thermal atmospheric band), and Eoutte (thermal emission).
    % Lout/Loutt/Eout are legacy preallocated scalar placeholders and stay
    % NaN in the production RTMo path, so they cannot be public sources.
    tryWriteSum(reg, 'longwave_outgoing', rad, {'Eoutt', 'Eoutte'}, KT);
    tryWrite(reg, 'longwave_thermal_emission', rad, 'Eoutte', KT);
    tryWrite(reg, 'shortwave_outgoing_hemis',  rad, 'Eouto', KT);
    tryWrite(reg, 'shortwave_outgoing_obs',    rad,    'Eouto',     KT);
    tryWriteScaled(reg, 'par_at_canopy_top',   rad,    'PAR',      1e6, KT);
    tryWriteScaled(reg, 'apar_canopy',         fluxes, 'aPAR',     1e6, KT);
    tryWriteScaled(reg, 'apar_cab',            fluxes, 'aPAR_Cab', 1e6, KT);
    tryWrite(reg, 'apar_canopy_wm2',           fluxes, 'aPAR_Wm2',  KT);
    tryWrite(reg, 'fluorescence_outgoing',     rad,    'Eoutf',     KT);

    % ----- energy_balance (13 own + 1 cross conservation) -----
    tryWrite(reg, 'land_surface_upward_latent_heat_flux', fluxes, 'lEtot', KT);
    tryWrite(reg, 'land_surface_upward_sensible_heat_flux', fluxes, 'Htot', KT);
    tryWrite(reg, 'soil_top_downward_energy_flux', fluxes, ...
        'soil_top_downward_energy_flux_W_m2', KT);
    tryWrite(reg, 'bare_ground_external_energy_flux', fluxes, ...
        'bare_ground_external_energy_flux_W_m2', KT);
    tryWrite(reg, ...
        'bare_ground_external_energy_flux_temperature_sensitivity', fluxes, ...
        'bare_ground_external_energy_flux_temperature_derivative_W_m2_K', KT);
    tryWrite(reg, 'scope_empirical_bare_ground_heat_flux', fluxes, ...
        'scope_empirical_bare_ground_heat_flux_W_m2', KT);
    tryWrite(reg, 'soil_top_energy_flux_reference', fluxes, ...
        'soil_top_energy_flux_reference_W_m2', KT);
    tryWrite(reg, 'latent_heat_canopy',        fluxes, 'lEctot',          KT);
    tryWrite(reg, 'latent_heat_canopy_wet',    fluxes, 'lE_wet_canopy',   KT);
    tryWrite(reg, 'latent_heat_canopy_dry',    fluxes, 'lE_dry_canopy',   KT);
    tryWrite(reg, 'ground_upward_latent_heat_flux', fluxes, ...
        'ground_upward_latent_heat_flux_W_m2', KT);
    tryWrite(reg, 'sensible_heat_canopy',      fluxes, 'Hctot',           KT);
    tryWrite(reg, 'ground_upward_sensible_heat_flux', fluxes, ...
        'ground_upward_sensible_heat_flux_W_m2', KT);
    tryWrite(reg, 'bare_ground_net_radiation_flux_contribution', fluxes, ...
        'bare_ground_net_radiation_contribution_W_m2', KT);
    tryWrite(reg, 'snow_net_radiation_flux_contribution', fluxes, ...
        'snow_net_radiation_contribution_W_m2', KT);
    tryWrite(reg, 'bare_ground_upward_latent_heat_flux_contribution', fluxes, ...
        'bare_ground_upward_latent_heat_flux_contribution_W_m2', KT);
    tryWrite(reg, 'snow_upward_latent_heat_flux_contribution', fluxes, ...
        'snow_upward_latent_heat_flux_contribution_W_m2', KT);
    tryWrite(reg, 'bare_ground_upward_sensible_heat_flux_contribution', fluxes, ...
        'bare_ground_upward_sensible_heat_flux_contribution_W_m2', KT);
    tryWrite(reg, 'snow_upward_sensible_heat_flux_contribution', fluxes, ...
        'snow_upward_sensible_heat_flux_contribution_W_m2', KT);
    tryWrite(reg, 'bare_ground_downward_soil_heat_flux_contribution', fluxes, ...
        'bare_ground_downward_soil_heat_flux_contribution_W_m2', KT);
    tryWrite(reg, 'snow_downward_soil_energy_flux_contribution', fluxes, ...
        'snow_downward_soil_energy_flux_contribution_W_m2', KT);
    tryWrite(reg, 'bare_ground_upward_water_vapor_flux', fluxes, ...
        'bare_ground_upward_water_vapor_flux_kg_m2_s', KT);
    tryWrite(reg, 'snow_upward_water_vapor_flux', fluxes, ...
        'snow_upward_water_vapor_flux_kg_m2_s', KT);
    tryWrite(reg, 'bare_ground_skin_temperature', thermal, ...
        'bare_ground_skin_temperature_C', KT);
    tryWrite(reg, 'snow_surface_temperature', thermal, ...
        'snow_surface_temperature_C', KT);
    tryWrite(reg, 'ground_radiative_temperature', thermal, ...
        'ground_radiative_temperature_C', KT);
    tryWrite(reg, 'latent_heat_wet_sun',       fluxes, 'lE_wet_sun',      KT);
    tryWrite(reg, 'latent_heat_wet_sha',       fluxes, 'lE_wet_sha',      KT);
    tryWrite(reg, 'latent_heat_dry_sun',       fluxes, 'lE_dry_sun',      KT);
    tryWrite(reg, 'latent_heat_dry_sha',       fluxes, 'lE_dry_sha',      KT);
    % energy_balance_residual is computed downstream; not in raw SCOPE structs

    % ----- carbon_fluxes (6) -----
    tryWrite(reg, 'gross_primary_productivity', fluxes, 'GPP',     KT);
    tryWrite(reg, 'net_ecosystem_exchange',     fluxes, 'NEE',     KT);
    tryWrite(reg, 'assimilation_canopy',        fluxes, 'Actot',   KT);
    tryWrite(reg, 'assimilation_total',         fluxes, 'Atot',    KT);
    if switchEnabled(cfg, 'lPlantCN') && switchEnabled(cfg, 'lNefc')
        tryWrite(reg, 'respiration_ecosystem',  fluxes, 'Resp',    KT);
        tryWrite(reg, 'respiration_plant', fluxes, 'Ra_plant', KT);
        tryWrite(reg, 'respiration_aom', fluxes, 'Rh_aom', KT);
        tryWrite(reg, 'respiration_mend', fluxes, 'Rh_mend', KT);
        tryWrite(reg, 'respiration_heterotrophic', fluxes, 'Rh', KT);
        tryWrite(reg, 'respiration_ecosystem_carbon', fluxes, 'Reco', KT);
    end
    tryWrite(reg, 'leaf_dark_respiration',      fluxes, 'Rd_leaf', KT);  % [Phase C M2-B.1a 2026-05-19]
    tryWrite(reg, 'stomatal_conductance_sun',   fluxes, 'gs_sun',  KT);  % [gs output 2026-06-13] gs=1/rcw, set in vegetation.energy.ebal / ebalProfile
    tryWrite(reg, 'stomatal_conductance_sha',   fluxes, 'gs_sha',  KT);  % [gs output 2026-06-13]

    % ----- canopy_aerodynamics (18) -----
    tryWrite(reg, 'air_temperature',                thermal, 'Ta',       KT);
    tryWrite(reg, 'soil_surface_temperature_avg',   thermal, 'Tsave',    KT);
    tryWrite(reg, 'canopy_temperature_avg',         thermal, 'Tcave',    KT);
    tryWrite(reg, 'air_density',                    thermal, 'rhoa',     KT);
    tryWrite(reg, 'resistance_aa',                  thermal, 'raa',      KT);
    tryWrite(reg, 'resistance_ah',                  thermal, 'rah',      KT);
    tryWrite(reg, 'resistance_aw',                  thermal, 'raw',      KT);
    tryWrite(reg, 'resistance_canopy_bulk',         thermal, 'rawc',     KT);
    tryWrite(reg, 'resistance_soil_aero',           thermal, 'raws',     KT);
    tryWrite(reg, 'resistance_soil_surface',        thermal, 'rss',      KT);
    tryWrite(reg, 'resistance_leaf_boundary',       thermal, 'rbc',      KT);
    tryWrite(reg, 'friction_velocity',              thermal, 'ustar',    KT);
    tryWrite(reg, 'roughness_length_momentum',      thermal, 'zo',       KT);
    tryWrite(reg, 'roughness_length_heat',          thermal, 'z0h',      KT);
    tryWrite(reg, 'displacement_height',            thermal, 'd',        KT);
    tryWrite(reg, 'stem_area_index',                thermal, 'SAI',      KT);
    tryWrite(reg, 'plant_area_index',               thermal, 'L_plus_S', KT);
    tryWrite(reg, 'soil_surface_saturation',        thermal, 'Sat_surf', KT);

    % ----- canopy_profiles (14, dims time × canopy_layer) -----
    tryWrite(reg, 'air_temperature_layer',          profiles, 'Ta_layer',     KT);
    tryWrite(reg, 'vapor_pressure_layer',           profiles, 'ea_layer',     KT);
    tryWrite(reg, 'co2_layer',                      profiles, 'CO2_layer',    KT);
    tryWrite(reg, 'eddy_diffusivity_layer',         profiles, 'Kh_layer',     KT);
    tryWrite(reg, 'leaf_temperature_layer',         profiles, 'Tleaf_layer',  KT);
    tryWrite(reg, 'canopy_temperature_layer',       profiles, 'Tch',          KT);
    tryWrite(reg, 'lai_layer',                      profiles, 'LAI_layer',    KT);
    tryWrite(reg, 'layer_weight',                   profiles, 'layer_weight', KT);
    tryWrite(reg, 'net_photosynthesis_layer',       profiles, 'Pn1d',         KT);
    tryWrite(reg, 'net_photosynthesis_cab_layer',   profiles, 'Pn1d_Cab',     KT);
    tryWrite(reg, 'sensible_heat_layer',            profiles, 'Hc1d',         KT);
    tryWrite(reg, 'latent_heat_layer',              profiles, 'lEc1d',        KT);
    tryWrite(reg, 'assimilation_layer',             profiles, 'A1d',          KT);
    tryWrite(reg, 'canopy_avg_temperature_3d',      profiles, 'Tchave',       KT);
end


function tryWriteScaled(reg, regName, srcStruct, srcField, scale, KT)
% Convert a numeric source to its public output unit before registry storage.
    if isempty(srcStruct) || ~isstruct(srcStruct) || ~isfield(srcStruct, srcField)
        return
    end
    value = srcStruct.(srcField);
    if isempty(value) || ~(isnumeric(value) || islogical(value)) || ~isreal(value)
        return
    end
    reg.write(regName, double(value) .* scale, KT);
end

function tryWriteSum(reg, regName, srcStruct, srcFields, KT)
% Write a scalar sum only when every accepted physical component is finite.
% Editor: Jizixing Zhu; Date: 2026-07-20.
    if isempty(srcStruct) || ~isstruct(srcStruct)
        return
    end
    values = nan(size(srcFields));
    for i = 1:numel(srcFields)
        name = srcFields{i};
        if ~isfield(srcStruct, name) || ~isscalar(srcStruct.(name)) ...
                || ~isfinite(srcStruct.(name))
            return
        end
        values(i) = double(srcStruct.(name));
    end
    reg.write(regName, sum(values), KT);
end


function tf = switchEnabled(cfg, name)
% Read one logical switch from either a normalized struct or ModelSettings.
    tf = false;
    if isstruct(cfg) && isfield(cfg, name) && ~isempty(cfg.(name))
        tf = logical(cfg.(name));
    elseif isobject(cfg) && isprop(cfg, name) && ~isempty(cfg.(name))
        tf = logical(cfg.(name));
    end
end


function tryWrite(reg, regName, srcStruct, srcField, KT)
% Silently skip if SCOPE source struct does not contain the field (some SCOPE
% variants emit a reduced subset). Keeps main loop oblivious to SCOPE branch.
    if isempty(srcStruct) || ~isstruct(srcStruct) || ~isfield(srcStruct, srcField)
        return
    end
    value = srcStruct.(srcField);
    if isempty(value) || ~(isnumeric(value) || islogical(value)) || ~isreal(value)
        return
    end
    reg.write(regName, double(value), KT);
end
