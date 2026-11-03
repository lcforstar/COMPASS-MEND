function ts = load_compass_state(nc_root, depth_node_cm)
%LOAD_COMPASS_STATE  Read full-year COMPASS output and pack into time-series state struct.
%
%   Returns ts struct with the following time series (each NT×1 unless noted):
%     ts.depth_cm         [cm, NN×1]    Node depths, top-down (from NC depth_node)
%     ts.theta_prof       [m³/m³, NN×NT]  Soil liquid water content profile
%     ts.T_prof_K         [K, NN×NT]    Soil temperature profile
%     ts.Ts_top_K         [K, NT×1]     Surface skin temperature (canopy avg)
%     ts.Tc_K             [K, NT×1]     Canopy temperature (LAI-weighted average)
%     ts.LAI              [m²/m², NT×1] Leaf area index
%     ts.SAI              [m²/m², NT×1] Stem area index (Phase 2 VOD_stem proxy)
%     ts.h2ocan           [mm, NT×1]    Canopy intercepted water
%     ts.psi_leaf_MPa     [MPa, NT×1]   Leaf water potential (steady-state diag)
%     ts.datetime         [datetime, NT×1]  Hourly timestamps, starts at 2021-01-01 00:00
%
%   Inputs:
%     nc_root             [string]  Path to output run folder containing nc/*.nc
%                                   e.g. 'data/output_v40_phaseD_year2021_v3'
%     depth_node_cm       [cm, NN×1] Node depths (top-down) — pass through from caller
%                                   so we don't read it twice
%
%   Implementation:
%     COMPASS NC stores profile variables as (node, time) — we transpose to
%     (NN, NT) so column t is the snapshot at hour t. Single-dimension
%     time series stay column vectors.
%
%   Unit conversions inside this loader:
%     soil_temperature       °C  → K   (+273.15)
%     canopy_temperature_avg °C  → K   (+273.15)
%     h2ocan                 mm  → mm  (no change; downstream divides by 1000)
%     leaf_potential         MPa → MPa (output already in MPa per COMPASS convention)
%
%   References:
%     COMPASS metadata.yaml (groups soil_water, soil_thermal, plant_canopy,
%     canopy_aerodynamics) — output_v33 schema, 8760 hourly steps per year.

% Resolve files
sw  = fullfile(nc_root, 'nc', 'soil_water.nc');
st  = fullfile(nc_root, 'nc', 'soil_thermal.nc');
pc  = fullfile(nc_root, 'nc', 'plant_canopy.nc');
aero = fullfile(nc_root, 'nc', 'canopy_aerodynamics.nc');

% ─── Soil profiles ─────────────────────────────────────────────────────────
% NetCDF stores (node, time), but ncread returns column-major Matlab orientation
% which collapses to (NT, NN). Transpose so column t is the profile at hour t.
ts.depth_cm   = depth_node_cm(:);
theta_raw     = ncread(sw, 'soil_moisture');
Tsoil_raw     = ncread(st, 'soil_temperature');
NN_in         = numel(ts.depth_cm);
if size(theta_raw, 1) ~= NN_in
    theta_raw = theta_raw.';
    Tsoil_raw = Tsoil_raw.';
end
ts.theta_prof = theta_raw;                              % (NN, NT)
ts.T_prof_K   = Tsoil_raw + 273.15;

% ─── Surface + canopy + air temperatures ─────────────────────────────────
%   Tsave = soil_surface_temperature_avg is the COMPASS-computed surface
%           skin T (includes shaded/sunlit average), more accurate as
%           microwave Tg than the topmost soil node (which is the bulk
%           top layer temperature).
%   Tcave = canopy_temperature_avg is LAI-weighted average over 30 layers
%           (canopy_temperature_layer × layer_weight), equivalent to the
%           γ-attenuated effective canopy emission temperature for L-band.
%   Ta    = air_temperature is the reference-height (Z_measure) Ta from
%           forcing, used by Pellarin atmospheric model (compute_atmosphere).
ts.Ts_top_K = ncread(aero, 'soil_surface_temperature_avg') + 273.15;
ts.Tc_K     = ncread(aero, 'canopy_temperature_avg')       + 273.15;
ts.Tair_K   = ncread(aero, 'air_temperature')              + 273.15;

% ─── Canopy structure + water ──────────────────────────────────────────────
ts.LAI    = ncread(pc, 'leaf_area_index');
ts.LAI_max = max(ts.LAI);     % annual max LAI for dormancy factor in compute_VOD
ts.SAI    = ncread(aero, 'stem_area_index');
ts.h2ocan = ncread(pc, 'h2ocan');                       % [mm]
ts.psi_leaf_MPa = ncread(pc, 'leaf_potential');         % [MPa, v1 steady-state]

% ─── Prognostic v2 plant water storage (optional, present if lProgHydraulic active) ──
%   These come from updateOrganStorageFromPsi (organ storage schemes).

%   If NC file lacks them (older runs), fields stay empty and downstream
%   compute_VOD will fall back to empirical (LAI/AGB) path.
ph = fullfile(nc_root, 'nc', 'plant_hydraulics.nc');
ts.W_root_m = [];
ts.W_stem_m = [];
ts.W_leaf_m = [];
ts.rwc_stem = [];           % [-] stem relative water content, drives L-VOD seasonality
ts.psi_leaf_v2_MPa = [];
if exist(ph, 'file')
    try
        info_ph = ncinfo(ph);
        ph_vars = {info_ph.Variables.Name};
        if any(strcmp(ph_vars, 'w_stem_v2'))
            ts.W_root_m = ncread(ph, 'w_root_v2');           % [m]
            ts.W_stem_m = ncread(ph, 'w_stem_v2');           % [m]
            ts.W_leaf_m = ncread(ph, 'w_leaf_v2');           % [m]
            ts.psi_leaf_v2_MPa = ncread(ph, 'psi_leaf_v2');  % [MPa, v2 prognostic]
        end
        % Stem relative water content (W/W_max): summer-low/winter-high in apple
        % orchard (transpiration draws down xylem water in growing season). This
        % is the physically correct driver of L-band VOD seasonality, matching
        % VODCA-L observed phase. Used by compute_VOD to modulate stem VWC.
        if any(strcmp(ph_vars, 'rwc_stem_v2'))
            ts.rwc_stem = ncread(ph, 'rwc_stem_v2');         % [-]
        end
        % Leaf relative water content (W/W_max), symmetric with rwc_stem. Lets
        % compute_VOD drive leaf VWC as M_leaf = C_leaf*WDR*rwc_leaf (carbon-driven
        % saturated water modulated by capacitance RWC), replacing the under-scaled
        % W_leaf_m*1000 path. Same rwc_leaf also feeds fluspect Cw downstream.
        if any(strcmp(ph_vars, 'rwc_leaf_v2'))
            ts.rwc_leaf = ncread(ph, 'rwc_leaf_v2');         % [-]
        end
    catch
        % NC has no v2 fields — silent fallback
    end
end

% ─── Carbon pools (optional, present if lNefc=1 = N module active) ───────
%   COMPASS plant_nitrogen.nc writes C pools when CN cycling is on [g C/m²].
%   [13c] stem split into livestem (live wood, respires) + deadstem (heartwood);
%   older runs wrote a single plant_c_stem_wood. We read whichever is present.
%   For microwave biomass driving, convert to dry biomass [kg/m²]:
%     dry_biomass = C × 2 / 1000      (C content ≈ 50% of dry biomass)
%   If pools absent (lNefc=0 or older runs), fields stay empty and
%   compute_VOD falls back to AGB_static constant via iBiomass=2 (auto).
ts.plant_c_stem_kgm2     = [];      % [kg/m²] dry biomass equivalent (live+dead total)
ts.plant_c_livestem_kgm2 = [];      % [kg/m²] [13c] live wood only
ts.plant_c_deadstem_kgm2 = [];      % [kg/m²] [13c] dead wood only
ts.plant_c_leaf_kgm2 = [];
ts.plant_c_root_kgm2 = [];
ts.plant_c_fruit_kgm2 = [];
ts.has_carbon_pools  = false;
ts.has_dual_stem     = false;       % [13c] true if live/dead split available
pn_file = fullfile(nc_root, 'nc', 'plant_nitrogen.nc');
if exist(pn_file, 'file')
    try
        info_pn = ncinfo(pn_file);
        pn_vars = {info_pn.Variables.Name};
        if any(strcmp(pn_vars, 'plant_c_livestem'))
            % [13c] dual stem pool
            ts.plant_c_livestem_kgm2 = ncread(pn_file, 'plant_c_livestem') * 2 / 1000;
            ts.plant_c_deadstem_kgm2 = ncread(pn_file, 'plant_c_deadstem') * 2 / 1000;
            ts.plant_c_stem_kgm2     = ts.plant_c_livestem_kgm2 + ts.plant_c_deadstem_kgm2;
            ts.plant_c_leaf_kgm2  = ncread(pn_file, 'plant_c_leaf')      * 2 / 1000;
            ts.plant_c_root_kgm2  = ncread(pn_file, 'plant_c_root_fine') * 2 / 1000;
            ts.plant_c_fruit_kgm2 = ncread(pn_file, 'plant_c_fruit')     * 2 / 1000;
            ts.has_carbon_pools   = true;
            ts.has_dual_stem      = true;
        elseif any(strcmp(pn_vars, 'plant_c_stem_wood'))
            % legacy single stem pool (pre-13c runs)
            ts.plant_c_stem_kgm2  = ncread(pn_file, 'plant_c_stem_wood') * 2 / 1000;
            ts.plant_c_leaf_kgm2  = ncread(pn_file, 'plant_c_leaf')      * 2 / 1000;
            ts.plant_c_root_kgm2  = ncread(pn_file, 'plant_c_root_fine') * 2 / 1000;
            ts.plant_c_fruit_kgm2 = ncread(pn_file, 'plant_c_fruit')     * 2 / 1000;
            ts.has_carbon_pools   = true;
        end
    catch
        % silent fallback
    end
end

% ─── Hourly time vector ────────────────────────────────────────────────────
NT = numel(ts.LAI);
ts.datetime = datetime(2021,1,1,0,0,0) + hours(0:NT-1)';
ts.NT = NT;
ts.NN = numel(ts.depth_cm);

% Sanity check
assert(size(ts.theta_prof, 1) == ts.NN, ...
       'theta_prof first dim = %d ≠ NN = %d', size(ts.theta_prof,1), ts.NN);
assert(size(ts.theta_prof, 2) == NT, ...
       'theta_prof second dim = %d ≠ NT = %d', size(ts.theta_prof,2), NT);

end
