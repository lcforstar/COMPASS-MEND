function flushToCsv(reg, output_dir)
%FLUSHTOCSV  Compatibility writer: dump each registered variable as a single CSV.
%
% Mirrors the legacy `Sim_output_<Name>.csv` layout so old plot scripts
% (PlotResults/*.m) keep working during the migration period. Schedule for
% removal once all plots read NetCDF.
%
% File naming: data/output_<tag>/Sim_output_<VariableName>.csv
%   - Column 1 = KT index (1-based), like the original writeCompassOutputs.m
%   - Columns 2:end = data
%   - First row = header (NaN-padded for backward compat with readmatrix)
%
% Inputs:
%   reg         [in]  diagnostics.Registry
%   output_dir  [in]  char
%
% Author: Zhu Jizixing, 2026-05-15
% Editor: Jizixing Zhu; Date: 2026-07-17.
% Change: 2026-07-17 Jizixing Zhu -- map finite-ponding complementarity
%         diagnostics to stable compatibility CSV names.
% Change: 2026-07-17 Jizixing Zhu -- replace legacy infiltration CSV aliases
%         with canonical boundary rate and amount names.

    if ~exist(output_dir, 'dir'); mkdir(output_dir); end

    % [2026-05-16] CSV 按 group 子目录: output_<tag>/csv/<group>/Sim_output_<X>.csv
    csv_root = fullfile(output_dir, 'csv');
    if ~exist(csv_root, 'dir'); mkdir(csv_root); end

    names = reg.list();
    n_done = 0;
    for k = 1:numel(names)
        name = names{k};
        md   = reg.getMeta(name);
        buf  = reg.get(name);

        % Skip if buffer all NaN (variable was never written)
        if all(isnan(buf(:)))
            continue;
        end

        % Use legacy variable casing for back-compat (snake_case → CamelCase known mapping)
        legacy_name = legacyCaseFor(name);

        % 子目录 = group; conservation 跨 group, 用 'conservation' 即可
        grp = md.group;
        grp_dir = fullfile(csv_root, grp);
        if ~exist(grp_dir, 'dir'); mkdir(grp_dir); end
        fpath = fullfile(grp_dir, sprintf('Sim_output_%s.csv', legacy_name));
        % Build CSV: header NaN + 1..N indices, then KT,value rows
        n_kt = size(buf, 1);
        n_cols = size(buf, 2);
        out = nan(n_kt + 1, n_cols + 1);
        out(1, 2:end) = 1:n_cols;             % header column indices
        out(2:end, 1) = (1:n_kt)';
        out(2:end, 2:end) = buf;
        % NaN → empty in csv (legacy writers used printable form)
        out(isnan(out)) = md.fill_value;
        writematrix(out, fpath);
        n_done = n_done + 1;
    end
    fprintf('[diag] wrote %d legacy CSVs to %s\n', n_done, output_dir);
end

% Change: 2026-09-02 Jizixing Zhu -- map the crop xsmr debt-writeoff ledger.


function legacy = legacyCaseFor(name)
% Map new snake_case to old CamelCase for back-compat. Anything not mapped
% stays as-is (snake_case is acceptable for new variables).
    persistent map
    if isempty(map)
        map = containers.Map();
        % soil_water
        map('soil_moisture')                = 'SoilMoisture';
        map('soil_matric_head')             = 'SoilMatricHead';
        map('soil_water_capacity')          = 'SoilWaterCapacity';
        map('hydraulic_conductivity')       = 'HydraulicConductivity';
        map('effective_saturation_vg')      = 'EffectiveSaturation_VG';
        map('water_storage_column')         = 'WaterStorage_mm';
        map('water_storage_change')         = 'WaterStorageChange_mm';
        % soil_thermal
        map('soil_temperature')             = 'SoilTemperature';
        map('bottom_temperature')           = 'BottomTemperature';
        map('surface_heat_flux')            = 'SurfaceHeatFlux';
        map('bottom_heat_flux')             = 'BottomHeatFlux';
        % soil_freezing
        map('soil_ice')                     = 'SoilIce';
        map('soil_unfrozen_water')          = 'SoilUnfrozenWater';
        map('ice_ratio')                    = 'IceRatio';
        map('frozen_matric_head')           = 'FrozenMatricHead';
        map('frozen_hydraulic_conductivity')= 'FrozenHydraulicConductivity';
        % soil_air
        map('air_pressure')                 = 'AirPressure';
        % soil_nitrogen
        map('nh4_liquid')                   = 'AmmoniumConc';
        map('no3_liquid')                   = 'NitrateConc';
        map('urea_liquid')                  = 'UreaConc';
        map('nh4_total')                    = 'AmmoniumTotalConc';
        map('no3_total')                    = 'NitrateTotalConc';
        map('nitrification_flux')           = 'NitrificationFlux';
        map('denitrification_flux')         = 'DenitrificationFlux';
        map('mineralization_flux')          = 'MineralizationFlux';
        map('urea_hydrolysis_flux')         = 'UreaHydrolysisFlux';
        map('ammonia_volatilization')       = 'AmmoniaVolatilization';
        map('root_sink_nh4')                = 'RootSinkNH4';
        map('root_sink_no3')                = 'RootSinkNO3';
        map('bot_flux_nh4')                 = 'BotFluxNH4';
        map('bot_flux_no3')                 = 'BotFluxNO3';
        map('n2o_nitrif_flux')              = 'N2O_NitrifFlux';
        map('nh4_storage')                  = 'NH4_storage';
        map('no3_storage')                  = 'NO3_storage';
        map('nh4_cum_residual')             = 'NH4_cum_residual';
        map('no3_cum_residual')             = 'NO3_cum_residual';
        map('fert_nh4_added')               = 'FertNH4_added';
        map('fert_no3_added')               = 'FertNO3_added';
        map('fert_urea_added')              = 'FertUrea_added';
        % M2-B.2 demand-driven uptake diagnostics (2026-05-19)
        map('root_demand_nh4')              = 'RootDemandNH4';
        map('root_demand_no3')              = 'RootDemandNO3';
        map('nitrogen_supply_limit_factor') = 'NSupplyLimitFactor';
        % plant
        map('leaf_area_index')              = 'LeafAreaIndex';
        map('leaf_potential')               = 'Leafpotential';
        map('root_potential')               = 'Rootpotential';
        map('soil_potential_weighted')      = 'Soilpotential';
        map('fwet')                         = 'fwet';
        map('fwet_sun')                     = 'fwet_sun';
        map('fwet_sha')                     = 'fwet_sha';
        map('h2ocan')                       = 'h2ocan_mm';
        map('h2ocan_sun')                   = 'h2ocan_sun_mm';
        map('h2ocan_sha')                   = 'h2ocan_sha_mm';
        map('canopy_water_residual')        = 'canopy_water_residual_mm';
        map('lE_canopy_intercept')          = 'lE_CanopyIntercept_W';
        map('vapor_pressure_deficit')       = 'VaporPressureDeficit';
        % plant_hydraulics (k_*/c_*/w_*/q_* unified naming, 2026-05-15 v2)
        map('k_plant_sun_app')              = 'KplantAppSun';
        map('k_plant_sha_app')              = 'KplantAppSha';
        map('k_plant_canopy_app')           = 'KplantAppTotal';
        map('k_plant_sun_diff')             = 'KplantDiffSun';
        map('k_plant_sha_diff')             = 'KplantDiffSha';
        map('k_plant_canopy_diff')          = 'KplantDiffFluxWeighted';
        map('k_soil_root')                  = 'KsoilrootAppTotal';
        map('k_soil_leaf')                  = 'KwholeAppTotal';
        map('gsr_node_sum')                 = 'GsrTotal';
        % surface_water
        map('precipitation_flux')           = 'PrecipitationFlux';
        map('precipitation_amount')         = 'PrecipitationAmount_mm';
        map('surface_liquid_water_supply_flux') = 'SurfaceLiquidWaterSupplyFlux';
        map('surface_liquid_water_supply_amount') = 'SurfaceLiquidWaterSupplyAmount_mm';
        map('wet_canopy_evaporation_flux')  = 'WetCanopyEvaporationFlux';
        map('wet_canopy_evaporation_amount') = 'WetCanopyEvaporationAmount_mm';
        map('soil_top_water_flux')          = 'SoilTopWaterFlux';
        map('soil_top_water_amount')        = 'SoilTopWaterAmount_mm';
        map('soil_bottom_water_flux')       = 'SoilBottomWaterFlux';
        map('soil_bottom_water_amount')     = 'SoilBottomWaterAmount_mm';
        map('surface_runoff_flux')          = 'SurfaceRunoffFlux';
        map('surface_runoff_amount')        = 'SurfaceRunoffAmount_mm';
        map('surface_water_storage')         = 'SurfaceWaterStorage';
        map('surface_water_storage_change')  = 'SurfaceWaterStorageChange_mm';
        map('surface_pressure_head')        = 'SurfacePressureHead';
        map('surface_water_available_flux') = 'SurfaceWaterAvailableFlux';
        map('surface_boundary_mode')        = 'SurfaceBoundaryMode';
        map('surface_mass_residual_max')    = 'SurfaceMassResidualMax';
        map('surface_lower_complementarity_residual_max') = ...
            'SurfaceLowerComplementarityResidualMax';
        map('surface_upper_complementarity_residual_max') = ...
            'SurfaceUpperComplementarityResidualMax';
        % evapotranspiration
        map('soil_evaporation_flux')        = 'SoilEvaporationFlux';
        map('soil_evaporation_amount')      = 'SoilEvaporationAmount_mm';
        map('soil_evaporation_demand_flux') = 'SoilEvaporationDemandFlux';
        map('soil_evaporation_demand_amount') = 'SoilEvaporationDemandAmount_mm';
        map('transpiration_flux')           = 'TranspirationFlux';
        map('transpiration_amount')         = 'TranspirationAmount_mm';
        % root_uptake
        map('water_uptake_profile')         = 'WaterUptake';
        map('rwu_pos_profile')              = 'RWU_pos';
        map('rwu_neg_profile')              = 'RWU_neg';
        map('hydraulic_lift_pos')           = 'HydraulicLift_F_pos_mmh';
        map('hydraulic_lift_net')           = 'HydraulicLift_F_net_mmh';
        map('hydraulic_lift_intensity')     = 'HydraulicLift_intensity';
        map('root_sink_srt_pos')            = 'RootSinkSrtPos_mm';
        map('root_sink_srt_neg')            = 'RootSinkSrtNeg_mm';
        map('root_sink_srt_net')            = 'RootSinkSrtNet_mm';
        map('root_sink_minus_transpiration')= 'RootSinkMinusTranspiration_mm';
        % dual_porosity
        map('mobile_water')                 = 'MobileWater';
        map('immobile_water')               = 'ImmobileWater';
        map('sink_immobile')                = 'SinkIm';
        % conservation
        map('water_balance_residual')       = 'WaterBalanceResidual_mm';
        map('surface_soil_liquid_input')    = 'SurfaceSoilLiquidInput_mm';
        map('surface_soil_evaporation_loss') = 'SurfaceSoilEvaporationLoss_mm';
        map('surface_soil_runoff_loss')     = 'SurfaceSoilRunoffLoss_mm';
        map('surface_soil_bottom_loss')     = 'SurfaceSoilBottomLoss_mm';
        map('surface_soil_root_sink')       = 'SurfaceSoilRootSink_mm';
        map('surface_soil_storage_change')  = 'SurfaceSoilStorageChange_mm';
        map('surface_soil_water_balance_residual') = ...
            'SurfaceSoilWaterBalanceResidual_mm';
        % numerics
        map('step_failed')                  = 'StepFailed';

        % --- 2026-05-15 新加: SCOPE 迁移字段 (旧名来自 Scope_output.<oldname>) ---
        % radiation
        map('net_radiation_total')          = 'Rntot';
        map('net_radiation_canopy')         = 'Rnctot';
        map('net_radiation_soil')           = 'Rnstot';
        map('longwave_outgoing')            = 'Lout';
        map('longwave_thermal_emission')    = 'Loutt';
        map('shortwave_outgoing_hemis')     = 'Eout';
        map('shortwave_outgoing_obs')       = 'Eouto';
        map('par_at_canopy_top')            = 'PAR';
        map('apar_canopy')                  = 'aPAR';
        map('apar_cab')                     = 'aPAR_Cab';
        map('apar_canopy_wm2')              = 'aPAR_Wm2';
        map('fluorescence_outgoing')        = 'Eoutf';
        % energy_balance
        map('land_surface_upward_latent_heat_flux') = 'land_surface_upward_latent_heat_flux';
        map('land_surface_upward_sensible_heat_flux') = 'land_surface_upward_sensible_heat_flux';
        map('soil_top_downward_energy_flux') = 'soil_top_downward_energy_flux';
        map('latent_heat_canopy')           = 'lEctot';
        map('latent_heat_canopy_wet')       = 'lE_wet_canopy';
        map('latent_heat_canopy_dry')       = 'lE_dry_canopy';
        map('ground_upward_latent_heat_flux') = 'ground_upward_latent_heat_flux';
        map('sensible_heat_canopy')         = 'Hctot';
        map('ground_upward_sensible_heat_flux') = 'ground_upward_sensible_heat_flux';
        map('latent_heat_wet_sun')          = 'lE_wet_sun';
        map('latent_heat_wet_sha')          = 'lE_wet_sha';
        map('latent_heat_dry_sun')          = 'lE_dry_sun';
        map('latent_heat_dry_sha')          = 'lE_dry_sha';
        map('energy_balance_residual')      = 'EnergyBalanceResidual_W';
        % carbon_fluxes
        map('gross_primary_productivity')   = 'GPP';
        map('net_ecosystem_exchange')       = 'NEE';
        map('assimilation_canopy')          = 'Actot';
        map('assimilation_total')           = 'Atot';
        map('respiration_ecosystem')        = 'Resp';
        map('leaf_dark_respiration')        = 'Rd_leaf';
        % plant_nitrogen (Phase C M2-B.1b 2026-05-19)
        map('plant_n_labile')               = 'PlantN_labile';
        map('plant_n_leaf')                 = 'PlantN_leaf';
        map('plant_n_root_fine')            = 'PlantN_root_fine';
        map('plant_n_livestem')             = 'PlantN_livestem';   % [13c] live wood N
        map('plant_n_deadstem')             = 'PlantN_deadstem';   % [13c] dead wood N
        map('plant_n_fruit')                = 'PlantN_fruit';
        map('net_primary_productivity')     = 'NPP';
        map('autotrophic_respiration')      = 'R_a_plant';
        map('maintenance_respiration_leaf') = 'R_m_leaf';
        map('maintenance_respiration_root') = 'R_m_root';
        map('r_m_leaf_actual')              = 'R_m_leaf_actual';
        map('r_m_root_actual')              = 'R_m_root_actual';
        map('r_m_stem_actual')              = 'R_m_stem_actual';
        map('r_m_livecroot_actual')         = 'R_m_livecroot_actual';
        map('r_m_fruit_actual')             = 'R_m_fruit_actual';
        map('r_g_actual')                   = 'R_g_actual';
        map('ra_above_actual')              = 'Ra_above_actual';
        map('ra_below_actual')              = 'Ra_below_actual';
        map('resp_theta_actual')            = 'resp_theta_actual';
        map('maintenance_respiration_stem') = 'R_m_stem';
        map('growth_respiration')           = 'R_g';
        map('plant_n_growth_factor')        = 'f_N_growth';
        map('plant_n_alloc_leaf')           = 'PlantN_alloc_leaf';
        map('plant_n_alloc_root')           = 'PlantN_alloc_root';
        map('plant_n_alloc_stem')           = 'PlantN_alloc_stem';
        map('plant_n_alloc_fruit')          = 'PlantN_alloc_fruit';
        map('plant_c_alloc_leaf')           = 'PlantC_alloc_leaf';
        map('plant_c_alloc_root')           = 'PlantC_alloc_root';
        map('plant_c_alloc_stem')           = 'PlantC_alloc_stem';
        map('plant_c_alloc_fruit')          = 'PlantC_alloc_fruit';
        % LUNA (M2-B.3 2026-05-20)
        map('vcmax25_luna')                 = 'Vcmax25_LUNA';
        map('jmax25_luna')                  = 'Jmax25_LUNA';
        map('luna_fnca')                    = 'LUNA_FNCa';
        map('luna_pnlc')                    = 'LUNA_PNlc';
        map('luna_pncb')                    = 'LUNA_PNcb';
        % M2-B.4 turnover/harvest (2026-05-21)
        map('plant_n_storage')              = 'PlantN_storage';
        map('n_harvest_today')              = 'N_harvest_today';
        map('n_leaf_litter_today')          = 'N_leaf_litter_today';
        map('n_root_litter_today')          = 'N_root_litter_today';
        map('n_spring_release_today')       = 'N_spring_release_today';
        map('cum_n_harvest')                = 'cum_N_harvest';
        map('cum_n_leaf_litter')            = 'cum_N_leaf_litter';
        map('cum_n_root_litter')            = 'cum_N_root_litter';
        map('cum_n_prune_litter')           = 'cum_N_prune_litter';
        map('cum_n_mortality')              = 'cum_N_mortality';
        map('frphu')                        = 'frPHU';
        % M2-B.8 Plant C 6 pools + cum_C + residual (2026-05-21)
        map('plant_c_labile')               = 'PlantC_labile';
        map('plant_c_leaf')                 = 'PlantC_leaf';
        map('plant_c_root_fine')            = 'PlantC_root_fine';
        map('plant_c_livestem')             = 'PlantC_livestem';   % [13c] live wood C
        map('plant_c_deadstem')             = 'PlantC_deadstem';   % [13c] dead wood C
        map('plant_c_fruit')                = 'PlantC_fruit';
        map('plant_c_storage')              = 'PlantC_storage';
        map('cum_c_harvest')                = 'cum_C_harvest';
        map('cum_c_leaf_litter')            = 'cum_C_leaf_litter';
        map('cum_c_root_litter')            = 'cum_C_root_litter';
        map('cum_c_prune_litter')           = 'cum_C_prune_litter';
        map('cum_c_mortality')              = 'cum_C_mortality';
        map('cum_c_xsmr_writeoff')          = 'cum_C_xsmr_writeoff';
        map('cum_gpp_plant')                = 'cum_GPP_plant';
        map('plant_c_residual')             = 'PlantC_residual';
        map('plant_n_residual')             = 'PlantN_residual';
        % M2-B.8b NSC 联动
        map('r_a_potential')                = 'R_a_potential';
        map('r_a_actual')                   = 'R_a_actual';
        map('storage_drawdown_event')       = 'StorageDrawdownEvent';
        map('nsc_total')                    = 'NSC_total';
        map('nsc_net')                      = 'NSC_net';
        map('c_storage_drawdown')           = 'C_storage_drawdown';
        map('cum_r_a_potential')            = 'cum_R_a_potential';
        map('cum_storage_drawdown_count')   = 'cum_StorageDrawdownCount';
        % canopy_aerodynamics
        map('air_temperature')              = 'Ta';
        map('soil_surface_temperature_avg') = 'Tsave';
        map('canopy_temperature_avg')       = 'Tcave';
        map('air_density')                  = 'rhoa';
        map('resistance_aa')                = 'raa';
        map('resistance_ah')                = 'rah';
        map('resistance_aw')                = 'raw';
        map('friction_velocity')            = 'ustar';
        map('roughness_length_momentum')    = 'zo';
        map('displacement_height')          = 'd';
        % soil_oxygen
        map('o2_concentration')             = 'C_O2';
        map('aerobic_fraction')             = 'f_aerobic';
        map('anaerobic_fraction')           = 'f_anaerobic';
        % soil_nitrogen 新名修正
        map('n2o_from_nitrification')       = 'Cvnit_N2o';
        map('n2o_from_denitrification')     = 'Cvdenit_N2o';
        map('n2_from_denitrification')      = 'Cvdenit_N2';
        map('top_flux_nh4')                 = 'TopFluxNH4';
        map('top_flux_no3')                 = 'TopFluxNO3';
    end
    if map.isKey(name)
        legacy = map(name);
    else
        legacy = name;  % new variables keep snake_case
    end
end
