# Leaf Potential Output Inventory And Plot Audit

Data source checked:

- `data/output_v30_4q_yearfull`
- Leaf water potential target: `Sim_output_Leafpotential.csv`
- Current analysis subset: only positive transpiration hours, `Transpiration_mm > 1e-8`
- Positive transpiration sample size: `N = 6210`

## Confirmed Keep List

These variables have been confirmed for the next redraw pass:

- `HydraulicLift_F_net_mmh`: root net uptake `F_net`, unit `mm h^-1`
- `Rootpotential`: root water potential, unit `MPa`
- `Soilpotential`: representative soil water potential seen by roots, unit `MPa`
- `Transpiration_mm`: transpiration, unit `mm h^-1`
- `VaporPressureDeficit`: atmospheric vapor pressure deficit, unit `kPa`
- `fluxes_lE_dry_canopy`: dry-canopy latent heat, unit `W m^-2`
- `fluxes_lEctot`: canopy latent heat, unit `W m^-2`
- `fluxes_aPAR_Wm2`: absorbed PAR energy, unit `W m^-2`
- `fluxes_Rnctot`: canopy net radiation, unit `W m^-2`
- `rad_Eoutf`: SIF-related emitted fluorescence, unit `W m^-2`
- `fluxes_Atot`: photosynthesis / CO2 uptake, unit `umol m^-2 s^-1`
- `thermal_Tcave`: canopy mean temperature, unit `degC`
- `Evaporation_mm`: evaporation, unit `mm h^-1`
- `LeafAreaIndex`: leaf area index, unitless

## 1. Sim Output Inventory

`Sim_output_*.csv` is mostly the model state, soil water/heat outputs, water fluxes, plant hydraulic diagnostics, and surface/canopy water balance terms.

### Plant / Canopy / Atmosphere

- `Leafpotential`
- `Rootpotential`
- `Soilpotential`
- `LeafAreaIndex`
- `VaporPressureDeficit`
- `Transpiration_mm`
- `Evaporation_mm`
- `fwet`
- `fwet_sha`
- `fwet_sun`
- `h2ocan_mm`
- `h2ocan_sha_mm`
- h2ocan_sun_mm
- CanopyInterception_mm
- lE_CanopyIntercept_W`

### Soil Water State / Hydraulic State

- `SoilMoisture`
- `SoilUnfrozenWater`
- `SoilIce`
- `IceRatio`
- `SoilMatricHead`
- `SoilWaterCapacity`
- 
- `EffectiveSaturation_VG`
- 
- `HydraulicConductivity`
- `FrozenHydraulicConductivity`
- `FrozenMatricHead`
- `WaterStorage_mm`
- 
- `WaterUptake`
- `RWU_pos`
- `RWU_neg`
- ``

### Soil / Surface Heat

- `SoilTemperature`
- 
- `SurfaceHeatFlux`
- ``

### Boundary Fluxes / Water Balance

- 
- `Precipitation_mm`
- 
- ``
- ``
- `SoilTopWaterAmount_mm`
- 
- `SurfaceLiquidWaterSupplyAmount_mm`
- 
- `SurfaceWaterAvailableFlux`
- 
- `SurfaceWaterStorageChange_mm`
- ``
- `SoilBottomWaterAmount_mm`
- 
- `SurfaceRunoffAmount_mm`
- `WaterBalanceResidual_mm`
- `AirPressure`
- `StepFailed`

### Hydraulic Lift Diagnostics

- `HydraulicLift_F_net_mmh`
- `HydraulicLift_F_pos_mmh`
- `HydraulicLift_intensity`

Important check: in the current correlation table, `HydraulicLift_F_net_mmh` has essentially the same relationship as `Transpiration_mm`:

```text
Transpiration_mm          r = -0.981267
HydraulicLift_F_net_mmh   r = -0.981267
```

So this variable should not be used as independent evidence for hydraulic lift before checking whether it is being written from the same flux branch or has a naming/export issue.

## 2. Scope Output Inventory

`Scope_output.mat` contains four major blocks by naming convention: `fluxes_*`, `thermal_*`, `profiles_*`, and `rad_*`.

### Fluxes: Canopy / Soil / Total Energy And Carbon

Scalar hourly outputs, all `[8760 1]` unless noted.

- `fluxes_Rntot`: total net radiation
- `fluxes_lEtot`: total latent heat
- `fluxes_Htot`: total sensible heat
- `fluxes_Atot`: total assimilation / photosynthesis
- `fluxes_Rnctot`: canopy net radiation
- `fluxes_lEctot`: canopy latent heat
- `fluxes_lE_wet_canopy`: wet canopy latent heat
- `fluxes_lE_dry_canopy`: dry canopy latent heat
- `fluxes_Hctot`: canopy sensible heat
- `fluxes_Actot`: canopy assimilation
- `fluxes_Rnstot`: soil net radiation
- `fluxes_lEstot`: soil evaporation latent heat
- `fluxes_Hstot`: soil sensible heat
- `fluxes_Gtot`: ground / soil heat flux
- `fluxes_Resp`: respiration
- `fluxes_aPAR`: absorbed PAR
- `fluxes_aPAR_Cab`: chlorophyll-related absorbed PAR
- `fluxes_aPAR_Wm2`: absorbed PAR in W m-2
- `fluxes_aPAR_Cab_eta`: effective absorbed PAR
- `fluxes_GPP`: gross primary productivity
- `fluxes_NEE`: net ecosystem exchange
- `fluxes_Ah`: canopy-layer assimilation profile `[30 8760]`

### Thermal / Aerodynamic / Resistance Outputs

- Air and canopy state:
  - `thermal_Ta`
  - `thermal_rhoa`
  - `thermal_Ts` `[2 8760]`
  - `thermal_Tcave`
  - `thermal_Tsave`
  - `thermal_Ta_layer` `[30 8760]`
  - `thermal_ea_layer` `[30 8760]`
  - `thermal_CO2_layer` `[30 8760]`
  - `thermal_Tch` `[30 8760]`

- Resistance and aerodynamic variables:
  - `thermal_raa`
  - `thermal_rah`
  - `thermal_raw`
  - `thermal_rawc`
  - `thermal_raws`
  - `thermal_rbc`
  - `thermal_rbc_wv`
  - `thermal_rss`
  - `thermal_rbs`
  - `thermal_rbs_init`
  - `thermal_rbs_static`
  - `thermal_ustar`
  - `thermal_ustar_clm`
  - `thermal_Va`
  - `thermal_Kh`
  - `thermal_uh`
  - `thermal_Kh_layer` `[30 8760]`

- Canopy geometry / roughness:
  - `thermal_z0h`
  - `thermal_z0w`
  - `thermal_zo`
  - `thermal_d`
  - `thermal_zo_static`
  - `thermal_d_static`
  - `thermal_SAI`
  - `thermal_L_plus_S`
  - `thermal_Kb_direct`
  - `thermal_Tao`
  - `thermal_n_canopy`
  - `thermal_Sat_surf`
  - `thermal_epsilon_s`
  - `thermal_rs_thermal`
  - `thermal_rs_thermal_static`

- - 

### Radiation Outputs

- Scalar radiation outputs:
  - `rad_Lout`
  - `rad_Loutt`
  - `rad_Eoutte`
  - `rad_PAR`
  - `rad_Eouto`
  - `rad_Eout`
  - `rad_Lout_`
  - `rad_Eoutt`
  - `rad_Rnhs`
  - `rad_Rnus`
  - `rad_Etoto`
  - `rad_Lot`
  - `rad_Rnust`
  - `rad_Rnhst`
  - `rad_Eoutf`

- - 

## 3. Current Plot Audit

All current response plots use:

- x-axis: one candidate predictor
- y-axis: `Sim_output_Leafpotential`
- gray points: hourly positive-transpiration samples
- black line: median binned response curve
- blue band: interquartile range
- `r`: Pearson correlation during positive-transpiration hours

### Figure: `leaf_potential_regression_formula.png`

This is the focused regression figure.

Panel `(a)`:

- Relationship: `Ψleaf` vs `Transpiration_mm`
- Meaning: transpiration demand / plant water flow effect on leaf water potential
- Result: `r = -0.98`
- Interpretation: this is the strongest direct relation. Higher transpiration makes leaf water potential more negative.

Panel `(b)`:

- Relationship: `Ψleaf` vs `Soilpotential`
- Meaning: soil water supply-side control
- Result: `r ≈ 0.34` in positive-transpiration hours
- Interpretation: physically meaningful, but not tight alone because the same soil potential can occur under different radiation/VPD/transpiration states.

Panel `(c)`:

- Relationship: `Ψleaf` vs `VaporPressureDeficit`
- Meaning: atmospheric demand control
- Result: `r ≈ -0.60`
- Interpretation: VPD lowers leaf water potential, but much of its effect is transmitted through transpiration.

Panel `(d)`:

- Relationship: actual/model `Ψleaf` vs fitted `Ψleaf`
- Fitted formula:

```text
Ψleaf = -0.042 + 0.973 Ψsoil - 0.0085 VPD - 8.993 T
```

- Result:

```text
R2 = 0.976
RMSE = 0.110 MPa
```

- Interpretation: the three-factor empirical formula captures most of the hourly positive-transpiration variability.

### Figure: `leaf_potential_top16_predictors.png`

This is an automatic ranking figure. It is useful for screening, not final publication as-is.

Top variables are:

- `Transpiration`
- `Transpiration_mm`
- `HydraulicLift_F_net_mmh`
- `fluxes_lE_dry_canopy`
- `fluxes_lEctot`
- `fluxes_lEtot`
- `fluxes_aPAR_Cab_eta`
- `fluxes_Rnctot`
- `rad_Eoutf`
- `fluxes_aPAR_Wm2`
- `fluxes_aPAR`
- `fluxes_aPAR_Cab`
- `fluxes_GPP`
- `fluxes_Atot`
- `fluxes_Actot`
- `fluxes_Ah_max`

Audit:

- `Transpiration` and `Transpiration_mm` are duplicates in different units.
- `GPP`, `Atot`, and `Actot` are effectively duplicates for this diagnostic.
- `aPAR`, `aPAR_Wm2`, `aPAR_Cab`, and `aPAR_Cab_eta` are related radiation/absorbed-radiation variants.
- `HydraulicLift_F_net_mmh` should be checked before use because it tracks `Transpiration_mm` almost exactly here.
- Good final candidates from this figure: `Transpiration_mm`, `lE_dry_canopy` or `lEctot`, `aPAR_Wm2` or `aPAR_Cab_eta`, `Rnctot`, `GPP/Atot`.

### Figure: `leaf_potential_photosynthesis_radiation.png`

This is the radiation/photosynthesis screening figure.

Current panels include:

- `Ψleaf` vs `aPAR_Cab_eta`
- `Ψleaf` vs `Rnctot`
- `Ψleaf` vs `rad_Eoutf`
- `Ψleaf` vs `aPAR_Wm2`
- `Ψleaf` vs `aPAR`
- `Ψleaf` vs `aPAR_Cab`
- `Ψleaf` vs `GPP`
- `Ψleaf` vs `Atot`
- `Ψleaf` vs `Actot`
- `Ψleaf` vs `Rntot`
- `Ψleaf` vs `rad_Rnhc_max`
- `Ψleaf` vs `rad_PAR`
- `Ψleaf` vs `rad_Emint_min`
- `Ψleaf` vs `rad_Eoutte`
- `Ψleaf` vs `rad_Eplut_min`
- `Ψleaf` vs `rad_Eplut_sum`

Audit:

- Best physical variables for final figure:
  - `aPAR_Wm2` or `aPAR_Cab_eta`: absorbed radiation driving photosynthesis and transpiration.
  - `Rnctot`: canopy net radiation.
  - `GPP` or `Atot`: photosynthetic uptake.
  - `rad_PAR`: incoming/available PAR signal.
- Variables such as `rad_Eoutf`, `rad_Emint_min`, `rad_Eplut_min`, and `rad_Eplut_sum` are radiative-transfer internal outputs. They may correlate strongly but are harder to explain in a paper figure.

### Figure: `leaf_potential_water_energy_hydraulics.png`

This is the water/energy/hydraulic screening figure.

Current panels include:

- `Ψleaf` vs `Transpiration`
- `Ψleaf` vs `Transpiration_mm`
- `Ψleaf` vs `HydraulicLift_F_net_mmh`
- `Ψleaf` vs `lEtot`
- `Ψleaf` vs `Evaporation_mm`
- `Ψleaf` vs `Evaporation`
- `Ψleaf` vs `Htot`
- `Ψleaf` vs `VaporPressureDeficit`
- `Ψleaf` vs `Gtot`
- `Ψleaf` vs `SurfaceHeatFlux`
- `Ψleaf` vs `Rootpotential`
- `Ψleaf` vs `Soilpotential`
- `Ψleaf` vs `SoilMoisture_mean`
- `Ψleaf` vs `SoilMoisture_max`
- `Ψleaf` vs `SoilMoisture_deep`
- `Ψleaf` vs `fwet_sha`

Audit:

- Best final candidates:
  - `Transpiration_mm`: direct plant water demand.
  - `VaporPressureDeficit`: atmospheric demand.
  - `Soilpotential`: soil supply-side water potential.
  - `Rootpotential`: plant hydraulic supply node.
  - `SoilMoisture_mean` or a better root-zone soil moisture metric: water storage proxy.
  - `lEtot` or `lEctot`: energy equivalent of evapotranspiration.
- Less ideal / needs caution:
  - `HydraulicLift_F_net_mmh`: suspiciously identical to transpiration in this output.
  - `SurfaceWaterAvailableFlux`, `SoilTopWaterFlux`, `WaterStorageChange_mm`: useful boundary diagnostics but not clean explanatory variables for a leaf-potential figure.
  - `SoilMoisture_deep`: in this output it has a very narrow range and is probably not useful as a final panel.

## 4. Recommended Redraw

For the next version, avoid automatic top-16 panels and use physically selected variables.

### Recommended Main Figure: Four Panels

Use the regression formula figure:

```text
(a) Ψleaf vs Transpiration_mm
(b) Ψleaf vs VPD
(c) Ψleaf vs Soilpotential
(d) Observed/Model Ψleaf vs fitted Ψleaf
```

This directly supports:

```text
Ψleaf = a + b Ψsoil + c VPD + d T
```

### Recommended Supplementary Figure: Radiation / Photosynthesis Link

```text
(a) Ψleaf vs aPAR_Wm2 or aPAR_Cab_eta
(b) Ψleaf vs Rnctot
(c) Ψleaf vs GPP or Atot
(d) Ψleaf vs Tcave or Tleaf_layer_mean
```

This tells the chain:

```text
radiation -> photosynthesis / canopy energy -> transpiration -> leaf water potential
```

### Recommended Supplementary Figure: Soil / Hydraulic Supply Link

```text
(a) Ψleaf vs Soilpotential
(b) Ψleaf vs Rootpotential
(c) Ψleaf vs root-zone SoilMoisture
(d) Ψleaf vs WaterUptake or RWU_pos
```

Before drawing this final version, root-zone soil moisture should be computed with a meaningful depth weighting, not only `profile_mean/max/deep`.

## 5. Missing Or Worth Adding

The current output already covers most useful relationships, but the following would make the figure stronger:

- A clean root-zone soil moisture index, e.g. 0-100 cm or root-weighted soil moisture.
- A clean root-zone soil water potential index if `Soilpotential` is not already root-weighted.
- A stomatal conductance / canopy conductance output if available from SCOPE internals.
- Incoming shortwave radiation or incoming PAR forcing if not already represented by `rad_PAR`.
- A checked hydraulic redistribution variable if hydraulic lift is central to the story.
