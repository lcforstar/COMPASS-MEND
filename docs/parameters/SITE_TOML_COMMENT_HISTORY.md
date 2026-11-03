# site.toml历史整行注释（归档）

当前文件只保留必要说明；以下原注释不代表当前活动配置。

```text
# =============================================================================
# Change: 2026-10-31 Wenjie Liu -- adopt accepted K1 mineral0-20cm annual
#   baseline: FAST10 and three N-origin cohorts enabled; prescribed root
#   residue events archived and inactive in this baseline. Numeric layered
#   kinetics and B_ref initial pools are in the corresponding CSV files.
# Change: 2026-10-23 Wenjie Liu -- 纯地下水热氮碳模式；关闭植物碳氮和还田，保留根系水氮汇。
# Change: 2026-10-30 Wenjie Liu -- add prescribed harvest root-residue C/N;
#   PlantCN and automatic litter stay off, and Jackson live roots stay unchanged.
#  COMPASS Site Configuration
#
#  Author:  Wenjie Liu
#  Contact: lcforstar@gmail.com
#  Editor: Wenjie Liu; Date: 2026-10-17.
#  Change: 2026-10-23 Wenjie Liu -- disable whole-plant C/N and its automatic
#    litter inputs; retain independent passive root N uptake and soil MEND.
#  Change: 2026-10-22 Wenjie Liu -- remove inactive apple/EPIC site fields
#    and unused legacy AOM fractions for the observed-LAI winter-wheat run.
#  Change: 2026-10-20 Wenjie Liu -- enable five core MEND C/N switches for Stage 2;
#    retain disabled microbial transport/tracer and the existing organic-input calendar.
#  Change: 2026-10-17 Wenjie Liu -- disable pH response and uncalibrated cell transport;
#    plant/manure/straw organic inputs now pass through their AOM buffers.
#
#  Function:
#    One self-contained configuration for a single COMPASS site run. It stores
#    scalar settings and process switches; tabular and array inputs remain in
#    the site folder. io.loadSite reads the following files when their owning
#    component is enabled by the simulation plan:
#      soil.csv            soil hydraulic, thermal, nitrogen carbon and microbial parameters
#      soil_initial_profile.csv  initial theta, T, gas pressure, bulk density, NH4, NO3, SOC
#      organic_initial_pools.csv native POC/MOC/MB/IOM fractions by depth
#      plant_initial_pools.csv   initial prognostic plant carbon and nitrogen pools
#      root.csv            measured root-length-density profile
#      forcing.csv         meteorological forcing time series
#      water_boundary.csv  prescribed water-head or water-flux time series
#      canopy_measurements.csv   measured LAI and canopy-structure time series
#
#    SCOPE also reads the shared NetCDF libraries selected in [model_library]:
#      data/model_library_nc/leaf_optics/       leaf optical constants
#      data/model_library_nc/atmosphere/        atmospheric transfer spectra
#      data/model_library_nc/soil_spectrum/     soil reflectance spectra
#
#    Format: TOML. on/off flags are true/false (loader maps to 1/0); units are in
#    # comments. No absolute paths; everything resolves relative to the site folder.
#    site-config-schema.json validates this file at load time and in VS Code;
#    the workspace .taplo.toml associates the schema with every site.toml.
# =============================================================================
# Spatial assembly: "point" uses one vertical column; "regional" assembles many columns.
# Point driver: "field" uses meteorology, canopy, plants, and the natural surface;
# "laboratory" uses only soil, initial state, and prescribed water boundaries.
# soil.csv depth semantics. "sample_linear" (default) treats depth_cm as a
# measurement location and interpolates primitive parameters to element centers.
# "horizon_bounds" treats depth_cm as each constant horizon's lower boundary.
# MEND samples use the same coordinates, with independently configurable mapping.
# --- wind height conversion (log-profile to z_measure) ---
# --- soil thermal optical properties ---
# --- atmospheric CO2 + soil-T climatology ---
# Change: 2026-10-18 Wenjie Liu -- empirical mineralization verification; disable the seven soil-C/DOM/microbe/litter/tracer switches.
# --- core equations: what the model solves ---
# --- soil water retention / freezing curves ---
# --- heat: methods + couplings ---
# --- soil water: transport ---
# --- biogeochemistry (C/N cycle settings) ---
# Change: 2026-10-23 Wenjie Liu -- user-selected unit kinetic pH factors;
# actual soil pH and NH4/NH3 chemical equilibrium remain separate inputs.
# ClimaLand v1.0 bulk SnowModel: prognostic total SWE, retained liquid water,
# and internal energy. Total post-canopy precipitation is split with the
# Jennings et al. (2018) temperature-relative-humidity relation used by
# ClimaLand's neural-snow forcing tools; it is not a hard 0-degC threshold.
# Native AOM/soil-carbon fixed structure. The soil_carbon_model switch is the
# sole master; organic N follows it. These rates are site-global because the
# current AOM equations do not support independent values by soil layer.
  # CoLM MOD_Initialize tau_l1/tau_l2_l3/tau_cwd; days = 365/year.
  # Ordered metabolic/cellulose/lignin/CWD. k_litter_1d is legacy-only below.
  # CoLM MOD_Const_PFT lf_* and fr_*; C and N apportioned once on entry.
# Provisional Changwu calibration candidate: faster AOM solubilization is
# required to avoid routing nearly all incorporated manure C into immobile POC.
# Frozen CoLM environmental-response constants:
# MOD_Initialize.F90 (Q10/frozen Q10/smp limits/rf_decomp) and
# MOD_BGC_Soil_BiogeochemDecompCascadeBGC.F90 (depth e-folding).
# PlantCN material quality is retained by the native litter adapter.
# Change: 2026-10-23 Wenjie Liu -- restore source MEND independent D/Q and kinetic sorption.
# Transport moves D/DON only; Q/QON stay in their source reaction owner.
# soil.csv DOM Kd is unused in this mode; its inorganic-N parameters are unchanged.
# Dispersivity remains depth resolved; microbial C/N uses the source hourly map.
# Pool-specific initialization/receiver C:N [g C g-1 N]. These values create
# the first N guess only; prognostic N pools are integrated independently and
# are never overwritten from C:N during the run.
# Change: 2026-10-31 Wenjie Liu -- Qmax/Kba/Kdes are solely owned by soil.csv; no duplicate scalar defaults.
# Change: 2026-10-27 Wenjie Liu -- move six uniform MEND values to site ownership and four DIN/dormancy values to soil.csv.
# The 2026-10-12 ten-column ownership note below is historical and superseded by this migration.
# Change: 2026-10-10 Wenjie Liu -- select native MEND actual-uptake respiration allocation.
# MEND commit 92323c7 MM/Arrhenius source values. Most kinetic entries below
# follow the final Pinitial assignment in MEND_namelist.nml; the three
# microbial DIN uptake entries remain its initial Pinitial(28:30) prior, not
# the final optimized values. Rates use MEND's hourly clock. Water,
# inorganic-N, and winter-wheat PFT parameters retain their COMPASS values.
# Change: 2026-10-12 Wenjie Liu -- MEND 8+2 reaction candidates are owned
# exclusively by soil.csv. Values are initially constant across its 18 rows;
# future depth changes belong in that table. Legacy all-scalar sites remain
# supported, but this site must not duplicate these ten keys in TOML.
# POM2/MOM sharing ratio is now read from the soil.csv profile.
# Change: 2026-10-29 Wenjie Liu -- restore microbial minimum C:N from unnatural 2.0 to 5.5
# to enable physiological overflow mineralization (deamination) under N excess.
# Source Fraction_type1=[0.07,0.37,0.56] for plant POM1/POM2/DOM.
# Change: 2026-10-10 Wenjie Liu -- direct plant/manure entry; N follows the
# source C:N closure while conserving the actual supplied organic N.
# User-selected trial prior: composted manure water-soluble C ~=10% of total C
# (Morvan & Nicolardot 2009, doi:10.1007/s00374-009-0355-1); remaining C in 7:37.
# This is an input mapping hypothesis, not measured Changwu pool fractions.
# Numerical and initialization controls. Porous-medium attachment,
# detachment, straining and dispersivity profiles are in soil.csv.
# Select all soil boundary conditions here. The current combination is suited
# to a field orchard with natural surface forcing and a deep groundwater table:
#   water: atmospheric infiltration/evaporation + free drainage
#   heat: atmospheric surface coupling + free-drainage liquid energy export
#   dry air: atmospheric surface pressure + hydrostatic bottom pressure
# --- Water boundary ---------------------------------------------------------
# Top choices:
#   atmospheric   Natural rainfall, infiltration, and evaporation (field default)
#   constant_head / variable_head   Known ponding head or laboratory head
#   constant_flux / variable_flux   Known infiltration, evaporation, or irrigation
#   zero_flux     Sealed surface with no water exchange
# Bottom choices:
#   free_drainage Gravity drainage above a deep water table (field default)
#   constant_head / variable_head   Known groundwater head
#   constant_flux / variable_flux   Known recharge or drainage flux
#   zero_flux     Impermeable lower boundary
#
# For prescribed heads or fluxes, edit time_s, top_value, and bottom_value in
# water_boundary.csv. Head uses [cm]; flux uses [cm/s]. A constant mode requires
# one repeated value. Flux is positive into the soil: positive top flux is
# downward infiltration, negative top flux is evaporation, positive bottom flux
# is upward recharge, and negative bottom flux is downward drainage.
# variable_head uses linear interpolation. variable_flux is piecewise constant
# and integrated exactly. time_s starts at simulation time zero and must cover
# the complete simulation period.
# Change: 2026-09-29 Wenjie Liu -- mobile DOC/DON leave with liquid drainage.
# Explicit bottom phase closure: zero vapor exchange, accepted QMB is liquid.
# Numbered organic boundary types: top 1 concentration, 2 inflow/nonvolatile;
# bottom 1 concentration, 2 inflow concentration/outflow zero gradient, 3 outlet.
# All *_mg_L values are mg element/L WATER, not soil mg/kg. No DIN keys changed.
# Mobile active/dormant cells and intracellular cargo follow bottom liquid drainage.
# Vectors are [MB, BD, ENZN]; C/N must be supplied together. No external cells.
# Attached/retained phases have no direct boundary exchange.
# Maximum liquid depression storage [cm]. Zero means no retained pond: once
# the atmospheric flux would make h_surface positive, h_surface=0 and the
# excess becomes runoff. Use a measured/calibrated microtopographic depth to
# retain ponding; do not substitute saturated conductivity for this capacity.
# --- Heat boundary ----------------------------------------------------------
# Top choices:
#   surface_energy_response   Experimental natural-surface boundary. Insert the
#                             temperature-linearized bare/snow external-energy
#                             response into the U-primary top row; SCOPE's
#                             empirical force-restore G remains diagnostic.
#   surface_temperature       Prescribe the frozen SCOPE Tss [degC]
#                             (production operator-split contract)
#   atmospheric              Match the frozen surface-energy tangent to the
#                             upper-20-cm dynamic soil heat response, then
#                             prescribe the unique matched G for this hour
#   specified_heat_flux       Use heat_top_value [J cm-2 s-1]
# Bottom choices:
#   seasonal_bottom_temperature  Use the calculated seasonal Tbh (field default)
#   fixed_temperature            Use heat_bottom_value [degC]
#   specified_heat_flux          Use heat_bottom_value [J cm-2 s-1]
#   free_drainage_internal_energy  Let bottom drainage carry liquid sensible
#                                  energy using the coupled bottom T and QMB;
#                                  heat_bottom_value is inactive
# Heat flux is positive into the soil. heat_top_value is inactive for the
# selected atmospheric mode. Tss remains a surface diagnostic and
# is not imposed on the prognostic soil-temperature node.
# --- Dry-air boundary -------------------------------------------------------
# Top choices:
#   atmospheric_pressure  Follow atmospheric pressure forcing (field default)
#   fixed_pressure        Use dry_air_top_value [Pa]
#   specified_mass_flux   Use dry_air_top_value [g cm-2 s-1]
# Bottom choices:
#   hydrostatic_pressure  Apply the depth-dependent hydrostatic pressure (default)
#   zero_mass_flux        Seal the bottom against dry-air exchange
#   specified_mass_flux   Use dry_air_bottom_value [g cm-2 s-1]
# Dry-air mass flux is positive into the soil and is not a pressure-change rate.
# --- Nitrogen boundary (legacy integer adapter) -----------------------------
# 1=fixed concentration, 2=specified solute flux, 3=natural/outflow boundary.
# Concentration uses [mol/m3]. The current setup has zero top input and solute
# transport with bottom drainage. Urea enters through fertilization and
# hydrolysis, not through a boundary-transport condition.
# one entry per event; arrays align by index
# which model state pool each event enters. Optional: when omitted it defaults from
# events_type (ammonium->nh4_granule, nitrate->no3_liquid, urea->urea).
#   nh4_granule = solid/granular NH4, slow release via the C_NH4_granule pool
#   nh4_liquid  = dissolved NH4 input, straight into C_NH4_new
#   no3_liquid  = nitrate (fast-dissolve, no granule pool) -> C_NO3_new
#   urea        = urea -> C_Urea_new, later hydrolysis
# Annual composted farmyard manure: 75 Mg ha-1 yr-1.
# User-supplied experimental composition (2026-09-16):
# organic C 17.68 g/kg, total N 1.97 g/kg, alkali-hydrolysable N 91 mg/kg.
# Applying composition to the same 75000 kg material basis gives:
# C = 1326 kg/ha; total N = 147.75 kg/ha; C:total-N = 17.68/1.97.
# The 2.3:1 manure:soil preparation ratio is not applied a second time:
# reported composition and application rate are treated as the applied material.
# Alkali-hydrolysable N (6.825 kg/ha) is part of total N, not extra N or TAN.
# TAN fraction below remains an unverified prior, pending measurement.
# User-selected composted-manure TAN/total-N prior: 5% (2026-09-16); replace when
# available. TAN enters the existing COMPASS NH4/NH3 pathway; organic C/N
# enters MEND directly using the manure-specific fractions above.
# Annual cattle manure is incorporated through the full cultivated layer.
# stageStep distributes event mass conservatively by exact layer overlap.
# Change: 2026-10-10 Wenjie Liu -- all event mass is placed uniformly by
# depth over the specified 0-20 cm band, following mineral-fertilizer placement.
# Surface events instead enter the top soil cell. No implicit surface remainder.
# No independent aboveground straw-return treatment is applied in the Changwu
# NM experiment. Prescribed roots above are separate from automatic PlantCN
# litter, which remains disabled. No unsupported stubble mass is added.
# --- aerodynamic resistance ---
# --- leaf inclination distribution (Verhoef LIDF, |LIDFa|+|LIDFb|<=1; feeds RTMo/RTMt lidf) ---
# --- foliage clumping for canopy radiation only; physical LAI remains unchanged ---
# Winter-wheat canopy clumping is selected from the PFT table.
# --- canopy layer profiles ---
# --- canopy-air profile feedback (ebal_clm45_profile research branch) ---
# --- LAI and canopy height ---
# --- vertical LAI profile shape ---
# --- canopy interception ---
# --- plant hydraulics ---
# hydraulic_scheme:
#   steady_spac      = v1 steady SPAC only (no organ water storage W)
#   storage_from_psi = diagnose W from psi via retention (VOD; current default)
#   storage_dynamic  = prognostic dW/dt path (same wrapper as storage_from_psi today)
# Legacy aliases still accepted: equilibrium -> steady_spac, dynamic -> storage_from_psi
# --- organ storage (active when storage_from_psi or storage_dynamic) ---
# current default: calSPAC psi + linear PV curve -> W (VWC for VOD).
# PV curve per organ:  water_stored = water_capacity * (retention_slope * psi + 1)
#   at psi=0 -> water_stored = capacity (full);  at psi = -1/slope -> water_stored = 0 (empty).
# --- tree structure for gpmax fallback ---
# tree_DBH   = 21.3   # [cm]          diameter at breast height; only used to estimate gpmax when hydraulic_gpmax <= 0
# tree_CAH   = 10.2   # [m2]          crown area per tree; only used to estimate gpmax when hydraulic_gpmax <= 0
# tree_Ksmax = 5      # [kg/m/s/MPa]  max stem hydraulic conductance; only used to estimate gpmax when hydraulic_gpmax <= 0
# --- root system ---
# root-soil hydraulic resistance
# root nitrogen uptake solver (DAISY single-root + Nacry membrane M-M)
# --- phenology (GDD model, GA+fmincon calibrated, RMSE=0.348) ---
# --- plant C/N coupling (C-dominant FlexibleCN; species Vmax/Km/CN in PFT table) ---
# --- photosynthesis ---
# --- stomatal conductance (Medlyn + optional heat term) ---
# Medlyn gs [mol/m2/s] and slope
# Optional heat-induced additive gs (KEEP OFF until calibrated)
# docs/plant/stomata_heat_induced_gs_design.md
# gs magnitudes [mol/m2/s]  (only used when heat_gs_enable=true)
# linear gate ramps: 0→1 between (lo, hi); g_heat = max * φ_T * φ_w * φ_D * φ_A
# =============================================================================
#  SCOPE radiative-transfer settings  (was input_data.xlsx 'options' sheet)
#  Fully merged here; no separate SCOPE Excel. File references (spectral/soil/
#  atmosphere) now resolve to data/model_library_nc/, not per-site.
# =============================================================================
# =============================================================================
#  SCOPE leaf & soil optics  (was input_data.xlsx 'inputdata' sheet)
#  !! values below are SCOPE TEMPLATE DEFAULTS carried over from the Yangling
#     example input_data.xlsx; NOT yet Changwu-calibrated. TODO: calibrate.
# =============================================================================
# BRDF observer angles (scope.calc_directional=true). NetCDF var "data" = N x 2 [tto_deg, psi_deg].
# fluorescence NPQ tuning (scope.calc_fluor)
# BSM spectral soil model (used only when scope.soilspectrum=1)
# --- iteration limits ---
# --- convergence tolerances ---
# --- frozen surface / soil-interface scalar closure ---
# --- numerical kernels ---
# --- discretization and stabilization ---
# --- runtime diagnostics and restart ---
# Editor: Wenjie Liu; Change: 2026-10-17 field-capacity evaporation scaling.
# Change: 2026-10-18 Wenjie Liu -- freeze sample-39 evaporation resistance.
```
