# COMPASS Output Catalog

## Surface, snow, and soil energy contract (2026-07-20)

Changwu production uses `heat_top_kind = atmospheric` for bare, partial-snow,
and full-snow states. After ebal and snow staging, the frozen inward-positive
surface heat response is intersected with the upper-20-cm soil admittance. The
matched total-energy flux is prescribed to the coupled `h/U/P_g` solve as one
fixed Neumann boundary; ebal skin temperature is not a soil Dirichlet.

```text
r = G_soil,actual - G_surface,frozen(T_evaluation)
```

The strict target is `abs(r) <= 5 W m-2`. A freeze/thaw discontinuity first
triggers finer accepted soil tiles from the same outer-step state. If strict
closure still fails, a finite fully converged soil candidate may be accepted
only when `r` passes both the configured absolute and flux-relative fallback
caps. Such a step remains `converged = 0`, is explicitly marked bounded, and
publishes the real residual to the energy ledger. An out-of-bound candidate
rolls back the complete outer step without advancing time. Ebal is not rerun
inside these trials. SCOPE's empirical force-restore `G` remains a separate
predictor/diagnostic.

`soil_interface_closure_initial_temperature`,
`soil_interface_closure_final_temperature`, initial/final residuals, both
matched fluxes, iteration count, and convergence flag audit the outer loop.
`soil_interface_closure_acceptance_mode` codes rejected/strict/refined/bounded
as -1/1/2/3. `soil_interface_closure_accepted_bounded_residual`,
`soil_interface_closure_bounded_tolerance`,
`soil_interface_closure_substep_retry_used`, and
`soil_interface_closure_chosen_substep_cap` expose the fallback decision.
`soil_top_energy_flux_total_actual` remains the authoritative U-primary
reaction; `surface_to_soil_interface_power_residual` checks it against the
accepted frozen surface target.

One canopy RTMo/RTMt calculation uses effective ground optics and fourth-power
thermal emission. Canopy PAR/aPAR is one field and is not area-weighted a
second time. Ground radiation and turbulent exchange are split into bare and
snow tile contributions and each contribution is area weighted exactly once.
`bare_ground_skin_temperature`, `snow_surface_temperature`, and
`ground_radiative_temperature` are distinct; snow temperature is missing when
no thermodynamic snow domain exists.

`snow.nc` owns prognostic `snow_water_equivalent` [cm],
`snow_liquid_water` [cm], and `snow_internal_energy` [J m-2]. Rain/snow phase
partition, snow drainage, vapor exchange, and frozen depression-water transfer
are transactional. A terminal solid trace is removed only when its latent
energy is below the coupled energy resolution; total SWE is never truncated by
a simple depth threshold. Complete meltwater and its retained enthalpy are
transferred once through `snow_complete_meltout_*` diagnostics.

The combined surface-snow-soil ledger uses the outer-closed frozen surface map
as the boundary target. Snow-to-soil conduction is internal and must cancel when
snow and soil storage ledgers are added. Main fields are:

| Field | Meaning |
|---|---|
| `bare_ground_external_energy_flux_actual` | accepted frozen surface target after removing internal snow-to-soil transfer |
| `bare_ground_external_energy_flux_reconstructed` | independent SCOPE-component reconstruction using accepted evaporation |
| `bare_ground_external_energy_flux_reconstruction_error` | accepted minus reconstructed bare flux; exposes dry-active-set and linearization differences |
| `snow_external_energy_flux` | radiation/turbulence/precipitation energy entering snow |
| `snow_to_soil_internal_energy_flux` | internal snow-to-soil transfer, positive into soil |
| `surface_external_energy_flux` | net external input to the combined snow+soil control volume |
| `soil_top_energy_flux_total_actual` | accepted U-primary top reaction |
| `surface_to_soil_interface_power_residual` | soil reaction minus the outer-closed frozen surface target |
| `scope_ground_heat_flux_empirical` | SCOPE force-restore diagnostic only |
| `scope_ground_ebal_residual` | SCOPE algebraic diagnostic, not the coupled conservation residual |

`surface_soil_total_energy_residual` and
`surface_soil_interface_cancellation_error` are written to `conservation.nc`
in J m-2. The latter must be near machine precision. Estimated liquid and dry-
air sensible-energy components are explicitly named
`diagnosed_bare_liquid_infiltration_sensible_energy_flux` and
`diagnosed_dry_air_sensible_energy_transport`; they are not added again to the
already total-energy boundary response.

Two initialization provenance fields remain in `surface_energy.nc`:

| Field | Meaning |
|---|---|
| `initial_bottom_temperature_mismatch_before_alignment` | Input 10 m node minus same-time temperature boundary before initialization alignment, K |
| `initial_bottom_temperature_mismatch_after_alignment` | Remaining 10 m node/boundary mismatch after rebuilding U and thermal auxiliaries, K |

## Canonical boundary and conservation contract (2026-07-17)

This section is authoritative for `surface_water.nc`,
`evapotranspiration.nc`, and the water fields in `conservation.nc`; older
tables later in this historical catalog are superseded for these groups.

All boundary rates use `cm s-1`. Their paired accepted outer-step amounts use
`mm`. A flux and an amount are never stored under the same physical name.

| Rate field | Amount field | Sign and meaning |
|---|---|---|
| `precipitation_flux` | `precipitation_amount` | Atmospheric precipitation before canopy and snow processing |
| `surface_liquid_water_supply_flux` | `surface_liquid_water_supply_amount` | Liquid water reaching the ground after canopy and snow processing |
| `wet_canopy_evaporation_flux` | `wet_canopy_evaporation_amount` | Wet-canopy evaporation removed before the ground supply |
| `soil_top_water_flux` | `soil_top_water_amount` | Accepted top soil flux, positive into soil (`-QMT`) |
| `soil_bottom_water_flux` | `soil_bottom_water_amount` | Accepted bottom flux, positive into soil (`QMB`) |
| `surface_runoff_flux` | `surface_runoff_amount` | Overflow leaving the surface reservoir |
| `soil_evaporation_flux` | `soil_evaporation_amount` | Accepted realized soil-surface evaporation loss |
| `soil_evaporation_demand_flux` | `soil_evaporation_demand_amount` | Upstream SCOPE evaporation demand before the dry-head constraint |
| `soil_evaporation_deficit_flux` | `soil_evaporation_deficit_amount` | Nonnegative demand minus accepted soil evaporation |
| `transpiration_flux` | `transpiration_amount` | Plant transpiration |
| `bare_surface_dew_flux` | `bare_surface_dew_amount` | Bare condensation routed as nonnegative liquid supply |
| `canopy_dew_flux` | `canopy_dew_amount` | Condensation retained by the canopy water store |
| `snow_sublimation_flux` | `snow_sublimation_amount` | Positive snow vapor loss |
| `snow_deposition_flux` | `snow_deposition_amount` | Atmospheric vapor deposited into snow |
| `evapotranspiration_total_flux` | `evapotranspiration_total_amount` | Positive ecological water loss: soil + wet canopy + transpiration + snow |
| `net_surface_vapor_exchange_flux` | `net_surface_vapor_exchange_amount` | Positive losses minus bare/canopy dew and snow deposition; may be negative |

`surface_water.nc` also contains surface storage, storage change, pressure
head, available-water flux, boundary mode, and the three surface active-set
residual maxima. Mode codes are 0 supply-limited, 1 finite storage, 2
capacity-limited runoff, and 3 dry-head evaporation-limited. The fourth mode
has zero retained surface water and is entered only when the converged top
endpoint's liquid and total water both reach the retention model's residual
state. It preserves that state's finite pressure head rather than resetting
the soil to the legacy `hcrit`; any liquid supply releases it automatically.
It does not disable liquid or vapor transport inside the top soil element.

The soil-column ledger is
`water_balance_residual = soil_top_input + soil_bottom_input - root_sink -
soil_storage_change`. The combined surface-reservoir plus soil-column ledger
is
`liquid_supply - soil_evaporation - runoff - bottom_loss - root_sink -
(surface_storage_change + soil_storage_change)`. Its seven terms are written
as `surface_soil_*` fields in `conservation.nc`, all in `mm`.

The former public fields named potential/actual infiltration, potential top
flux, and mass flux were removed because they mixed forcing candidates,
kernel compatibility values, and accepted interface fluxes. Use
`soil_top_water_*` for the accepted soil exchange. Kernel compatibility names
remain private to the legacy adapter and never enter acceptance or public
output.

完整字段目录（**191 个 variables × 20 个物理 group**），跟 `src/+diagnostics/` 注册中心一一对应。

**v2.1 (2026-05-15)** — 在 v2 注册中心基础上把 Scope_output (123 字段) + flux_vars (17 字段) + NitrogenFlux 反硝化 + O2 模块全部并入。从"散装 `Sim_output.X = ...`"升级到 `+diagnostics/` 注册式 NetCDF 输出 (CF-1.10)。

每个字段在 `src/+diagnostics/+groups/register*.m` 里**只声明一次**，写出时按 group 分到 `data/output_<tag>/<group>.nc`。本文档是这些注册的人话索引。

> 旧文件名为 `Sim_output_naming_log.md`，2026-05-15 改名为 `output_catalog.md`（更贴 NetCDF data-archive 行话）。

---

## 设计原则（v2）

| 项 | v1 (旧) | v2 (新) |
|---|---|---|
| 注册位置 | 散在 `COMPASS.m` 103 处 | 集中在 `src/+diagnostics/+groups/register*.m` |
| 命名风格 | 驼峰混蛇形（如 `SoilMoisture` / `MobileWater`） | **snake_case + CF standard_name 双标** |
| 单位记录 | 注释或 md 文件，易漂移 | 写进 NetCDF attribute，跟数据绑定 |
| 输出格式 | 99 个 CSV | 13 个 .nc（按物理 group）+ 旧 CSV 兼容期 |
| 维度元数据 | 看 shape 猜 | NetCDF 显式 dims = `{time, node}` 等 |
| 守恒残差归类 | 跟物理字段混 | 全部进 `conservation.nc` 单独 group |
| Conventions | 无 | CF-1.10 |
| paper data archive | 不可用 | Zenodo 直接上传 |

---

## 20 个 Group 总览

| Group | 文件 | 触发条件 | 变量数 |
|---|---|---|---|
| `soil_water` | soil_water.nc | always | 7 |
| `soil_thermal` | soil_thermal.nc | `cfg.lThmrlefc=1` (默认) | 10 |
| `soil_freezing` | soil_freezing.nc | `cfg.lFreezing=1` | 5 |
| `soil_air` | soil_air.nc | `cfg.lSoilairefc=1` | 8 |
| `soil_nitrogen` | soil_nitrogen.nc | `cfg.lNefc=1` | **28** |
| `soil_oxygen` | soil_oxygen.nc | `cfg.lNefc=1` | 5 |
| `plant_canopy` | plant_canopy.nc | always | 12 |
| `plant_hydraulics` | plant_hydraulics.nc | always | 9 |
| `surface_water` | surface_water.nc | always | 20 |
| `evapotranspiration` | evapotranspiration.nc | always | 6 |
| `root_uptake` | root_uptake.nc | always | 9 |
| `dual_porosity` | dual_porosity.nc | `cfg.lDualPor=1` | 3 |
| `conservation` | conservation.nc | 跨 group 收集 | 17 |
| `numerics` | numerics.nc | always | 1 |
| **`radiation`** | radiation.nc | always | 12 |
| **`energy_balance`** | energy_balance.nc | always | 13 |
| **`carbon_fluxes`** | carbon_fluxes.nc | always | 5 |
| **`canopy_aerodynamics`** | canopy_aerodynamics.nc | always | 18 |
| **`canopy_profiles`** | canopy_profiles.nc | always | 14 |
| **`solver_internal`** | solver_internal.nc | always | 17 |

**新增 6 个 group**（2026-05-15 v2.1）：把 `Scope_output.mat` (123 字段) 和 `flux_vars` (17 字段) 分别按物理过程拆分进 NetCDF 体系，不再用单 .mat 大杂烩。

---

## 通用维度

| 维度 | 物理含义 | 大小 |
|---|---|---|
| `time` | KT 索引（外步末，CF datetime） | tS（如 8760） |
| `node` | 土壤节点（含表层、底层），node 1 = surface | NN（如 51） |
| `layer` | 土壤层（元素中点） | NL（如 50） |
| `nD` | 元素端点 dual node | 2 |

**坐标变量**：
- `time`: `units = "hours since 2021-01-01 00:00:00"`, `calendar = "standard"`
- `depth_node`: `units = "cm"`, `positive = "down"`, `axis = "Z"`，顶 → 底 = [0, 2, 5, ..., 1000]

---

## Group 1: soil_water (12 fields)

**输出**: `soil_water.nc`
**注册文件**: [src/+diagnostics/+groups/registerSoilWater.m](../src/+diagnostics/+groups/registerSoilWater.m)
**Source-of-truth**: src/+soil/+water/*

| 字段 (snake_case) | 旧 CSV 名 | 单位 | 维度 | 来源 | standard_name |
|---|---|---|---|---|---|
| `soil_moisture` | SoilMoisture | m3 m-3 | time × node | `SoilVariables.Theta_LL` | mass_content_of_water_in_soil_layer |
| `soil_matric_head` | SoilMatricHead | cm | time × node | `SoilVariables.hh` | soil_matric_potential |
| `soil_saturated_node_mask` | — | 1 | time × node | accepted `h >= 0` | — |
| `zero_head_crossing_depth` | — | cm | time × layer | all linearly interpolated `h=0` crossings | — |
| `zero_head_crossing_count` | — | 1 | time | crossing count, including perched zones | — |
| `shallowest_zero_head_depth` | — | cm | time | shallowest `h=0` crossing | depth_below_geopotential_datum |
| `saturated_profile_fraction` | — | 1 | time | integrated profile fraction with `h >= 0` | — |
| `soil_water_capacity` | SoilWaterCapacity | cm-1 | time × node | `MatricCoeff.Chh` | — |
| `hydraulic_conductivity` | HydraulicConductivity | cm s-1 | time × node | `MatricCoeff.Khh` | soil_hydraulic_conductivity_at_saturation |
| `effective_saturation_vg` | EffectiveSaturation_VG | 1 | time × node | `SoilVariables.Se` | — |
| `water_storage_column` | WaterStorage_mm | mm | time | `solver.calculateSoilWaterStorageMm` | — |
| `water_storage_change` | WaterStorageChange_mm | mm h-1 | time | computed in COMPASS.m | — |

---

## Group 2: soil_thermal (10 字段，另有 1 个残差写入 conservation)

**输出**: `soil_thermal.nc`
**注册文件**: [src/+diagnostics/+groups/registerSoilThermal.m](../src/+diagnostics/+groups/registerSoilThermal.m)

| 字段 | 旧 CSV 名 | 单位 | 维度 | 来源 | standard_name |
|---|---|---|---|---|---|
| `soil_temperature` | SoilTemperature | degC | time × node | `SoilVariables.TT` | soil_temperature |
| `bottom_temperature` | BottomTemperature | degC | time | seasonal bottom-temperature driver (active only in seasonal mode) | — |
| `surface_heat_flux` | SurfaceHeatFlux | W m-2 | time | `QET(KT)` | downward_heat_flux_in_soil |
| `bottom_heat_flux` | BottomHeatFlux | W m-2 | time | `QEB(KT)` | — |
| `soil_bottom_heat_boundary_mode` | — | 1 | time | 1=fixed T; 2=fixed flux; 3=free-drainage liquid internal energy; 4=mixed | — |
| `soil_bottom_boundary_temperature` | — | degC | time | requested bottom Dirichlet temperature; otherwise NaN | — |
| `soil_bottom_requested_heat_flux` | — | W m-2 | time | prescribed bottom flux; state-dependent drainage mode is NaN | — |
| `soil_internal_energy` | — | J m-2 | time | `workspace.flux.soil_energy_balance` | — |
| `soil_internal_energy_change` | — | J m-2 | time | 同上 | — |
| `soil_energy_top_input` | — | J m-2 | time | 同上 | — |
| `soil_energy_bottom_input` | — | J m-2 | time | 同上 | — |
| `soil_energy_boundary_net_input` | — | J m-2 | time | 同上 | — |
| `soil_energy_root_sink` | — | J m-2 | time | 同上 | — |

`soil_energy_balance_residual` 使用同一台账，但写入 `conservation.nc`。

---

## Group 3: soil_freezing (5 字段, lFreezing=1)

**输出**: `soil_freezing.nc`
**注册文件**: [src/+diagnostics/+groups/registerSoilFreezing.m](../src/+diagnostics/+groups/registerSoilFreezing.m)

| 字段 | 旧 CSV 名 | 单位 | 维度 | 来源 | standard_name |
|---|---|---|---|---|---|
| `soil_ice` | SoilIce | m3 m-3 | time × node | `SoilVariables.Theta_II` | volume_fraction_of_frozen_water_in_soil |
| `soil_unfrozen_water` | SoilUnfrozenWater | m3 m-3 | time × node | `SoilVariables.Theta_LL` (during freezing) | — |
| `ice_ratio` | IceRatio | 1 | time × node | `SoilVariables.Ratio_ice` | — |
| `frozen_matric_head` | FrozenMatricHead | cm | time × node | `SoilVariables.hh_frez` | — |
| `frozen_hydraulic_conductivity` | FrozenHydraulicConductivity | cm s-1 | time × node | `SoilVariables.KfL_h` (Lundin 1990 impedance) | — |

---

## Group 4: soil_air (8 字段，另有 1 个残差写入 conservation；lSoilairefc=1)

**输出**: `soil_air.nc`
**注册文件**: [src/+diagnostics/+groups/registerSoilAir.m](../src/+diagnostics/+groups/registerSoilAir.m)

| 字段 | 旧 CSV 名 | 单位 | 维度 | 来源 | standard_name |
|---|---|---|---|---|---|
| `air_pressure` | AirPressure | Pa | time × node | `SoilVariables.P_g` | air_pressure |
| `top_dry_air_mass_flux` | — | g cm-2 s-1 | time | `soil_boundary_hourly.dry_air` | — |
| `bottom_dry_air_mass_flux` | — | g cm-2 s-1 | time | `soil_boundary_hourly.dry_air` | — |
| `soil_dry_air_mass` | — | g m-2 | time | `soil_dry_air_mass_balance` | — |
| `soil_dry_air_mass_change` | — | g m-2 | time | 同上 | — |
| `soil_dry_air_top_input` | — | g m-2 | time | 同上 | — |
| `soil_dry_air_bottom_input` | — | g m-2 | time | 同上 | — |
| `soil_dry_air_boundary_net_input` | — | g m-2 | time | 同上 | — |

`soil_dry_air_mass_balance_residual` 使用同一台账，但写入 `conservation.nc`。

---

## Group 5: soil_nitrogen (30 字段写入 + 2 跨 conservation, lNefc=1)

**输出**: `soil_nitrogen.nc`
**注册文件**: [src/+diagnostics/+groups/registerSoilNitrogen.m](../src/+diagnostics/+groups/registerSoilNitrogen.m)

**注意**：`nh4_cum_residual` 和 `no3_cum_residual` 虽然在本文件里注册，但 `group='conservation'`，所以写入 `conservation.nc`（守恒残差跨 group 集中，见 §Group 12）。本 group 实际写入 `soil_nitrogen.nc` = 28 字段。

### 5.1 浓度场（5 字段, time × node）

| 字段 | 旧 CSV 名 | 单位 | 来源 | standard_name |
|---|---|---|---|---|
| `nh4_liquid` | AmmoniumConc | ug cm-3 | `SoilVariables.C_NH4_new` | mass_concentration_of_ammonium_in_soil_water |
| `no3_liquid` | NitrateConc | ug cm-3 | `SoilVariables.C_NO3_new` | mass_concentration_of_nitrate_in_soil_water |
| `urea_liquid` | UreaConc | ug cm-3 | `SoilVariables.C_Urea` | — |
| `nh4_total` | AmmoniumTotalConc | ug g-1 | `(θ+RHOKG) × C_NH4 / RHO_bulk` | — |
| `no3_total` | NitrateTotalConc | ug g-1 | `θ × C_NO3 / RHO_bulk` | — |

### 5.2 反应剖面（4 字段, time × node, 单位 ug cm-3 s-1）

NL×2 体积源/汇密度，可画反应深度云图，深度积分后等于 §5.3 通量。

| 字段 | 旧 CSV 名 | 物理含义 |
|---|---|---|
| `nitrification_profile` | — (新增) | NH4→NO3 反应源密度 |
| `denitrification_profile` | — | NO3→N2/N2O 汇密度 |
| `mineralization_profile` | — | SOM→NH4 源密度 |
| `root_sink_profile` | — | 根吸 N 总汇密度 (NH4+NO3) |

### 5.3 反应/边界通量（14 字段, time, 单位 ug cm-2 s-1）

剖面积分后的"深度积分通量"，给守恒账本用。分四组：反应、根系、边界、气体逸出。

**反应通量** (5)：

| 字段 | 旧 CSV 名 | 物理含义 | 来源代码 |
|---|---|---|---|
| `nitrification_flux` | NitrificationFlux | NH4 → NO3（深度积分） | nitrification.m |
| `denitrification_flux` | DenitrificationFlux | NO3 → N2/N2O | denitrification.m |
| `mineralization_flux` | MineralizationFlux | SOM → NH4 | mineralization.m |
| `urea_hydrolysis_flux` | UreaHydrolysisFlux | Urea → NH4 | ureaHydrolysis.m |
| `ammonia_volatilization` | AmmoniaVolatilization | NH3 大气逸出 | ammoniumVolatilization.m |

**根系吸收** (2)：

| 字段 | 旧 CSV 名 | 物理含义 |
|---|---|---|
| `root_sink_nh4` | RootSinkNH4 | 根系吸 NH4 |
| `root_sink_no3` | RootSinkNO3 | 根系吸 NO3 |

**边界通量** (4)：原仅 bot_flux 两件，v2.1 补 top_flux 两件（顶层 NH4/NO3 沉降+施肥地面输入）。

| 字段 | 旧 CSV 名 | 物理含义 | 来源代码 |
|---|---|---|---|
| `bot_flux_nh4` | BotFluxNH4 | 底部 NH4 通量 | calculateBoundaryFlux.m |
| `bot_flux_no3` | BotFluxNO3 | 底部 NO3 通量 | calculateBoundaryFlux.m |
| `top_flux_nh4` | — (新增) | 顶部 NH4 输入 | applyBoundaryConditions.m |
| `top_flux_no3` | — (新增) | 顶部 NO3 输入 | applyBoundaryConditions.m |

**气体逸出 N₂O/N₂** (3)：v2.1 拆分硝化与反硝化路径，**取代旧的单一 `n2o_nitrif_flux`**。

| 字段 | 旧 CSV 名 | 物理含义 |
|---|---|---|
| `n2o_from_nitrification` | (旧 N2O_NitrifFlux 拆分) | 硝化路径产 N₂O |
| `n2o_from_denitrification` | — (新增) | 反硝化路径产 N₂O |
| `n2_from_denitrification` | — (新增) | 反硝化路径产 N₂（最终产物） |

### 5.4 总储量（2 字段, time, 单位 ug cm-2）

| 字段 | 旧 CSV 名 | 含义 |
|---|---|---|
| `nh4_storage` | NH4_storage | 全剖面 NH4 总（液+吸附+冰） |
| `no3_storage` | NO3_storage | 全剖面 NO3 总（液+冰） |

### 5.5 施肥事件（3 字段, time, 单位 ug cm-2）

每次施肥那一个 KT 写入实际加入量，其他 KT = 0。

| 字段 | 旧 CSV 名 | 来源 |
|---|---|---|
| `fert_nh4_added` | FertNH4_added | `FertDiagBuffer.nh4_added` |
| `fert_no3_added` | FertNO3_added | `FertDiagBuffer.no3_added` |
| `fert_urea_added` | FertUrea_added | `FertDiagBuffer.urea_added` |

### 5.6 跨 conservation 的 2 字段（在 conservation.nc）

| 字段 | 单位 | 物理含义 |
|---|---|---|
| `nh4_cum_residual` | ug cm-2 | ΔS_NH4 − Σ(in − out)，N 守恒残差 |
| `no3_cum_residual` | ug cm-2 | 同上 NO3 |

### 5.7 DAISY MM 硝化诊断（2 字段, time × node）— 2026-05-16 新增

跟 2026-05-16 的 N 化学改造对应（线性 K_nit → DAISY 风格 Michaelis-Menten 显式源汇）：

| 字段 | 单位 | 物理含义 |
|---|---|---|
| `nit_mm_factor` | - (0..1) | M_total / (K_m + M_total)，1=底物充足无限制，~0.15=低 NH4 受限 |
| `nh4_substrate_total` | ug cm-3 soil | M_total = (θ + ρKd) · C_NH4_liq，DAISY 风格底物总量 |

源代码: `src/+soil/+nitrogen/+chemistry/nitrification.m` (重写后 Snit = V_max × eT × f_aero × ePH × mm_factor)

**子节合计**: 5.1 (5) + 5.2 (4) + 5.3 (14) + 5.4 (2) + 5.5 (3) + 5.7 (2) = **30 写入 soil_nitrogen.nc** + 5.6 (2) 跨 conservation = **32 注册**。

---

## Group 6: plant_canopy (12 字段)

**输出**: `plant_canopy.nc`
**注册文件**: [src/+diagnostics/+groups/registerPlantCanopy.m](../src/+diagnostics/+groups/registerPlantCanopy.m)

| 字段 | 旧 CSV 名 | 单位 | 维度 | 备注 |
|---|---|---|---|---|
| `leaf_area_index` | LeafAreaIndex | m2 m-2 | time | standard_name = leaf_area_index |
| `leaf_potential` | Leafpotential | MPa | time | — |
| `root_potential` | Rootpotential | MPa | time | — |
| `soil_potential_weighted` | Soilpotential | MPa | time | gsr-加权 |
| `fwet` | fwet | 1 | time | 全冠湿润比例 |
| `fwet_sun` | fwet_sun | 1 | time | 4-quadrant only |
| `fwet_sha` | fwet_sha | 1 | time | 4-quadrant only |
| `h2ocan` | h2ocan_mm | mm | time | 冠层持水 |
| `h2ocan_sun` | h2ocan_sun_mm | mm | time | — |
| `h2ocan_sha` | h2ocan_sha_mm | mm | time | — |
| `lE_canopy_intercept` | lE_CanopyIntercept_W | W m-2 | time | — |
| `vapor_pressure_deficit` | VaporPressureDeficit | kPa | time | standard_name = water_vapor_partial_pressure_deficit_in_air |

`canopy_water_residual` → 归到 `conservation.nc`

---

## Group 7: plant_hydraulics (9 字段)

**输出**: `plant_hydraulics.nc`
**注册文件**: [src/+diagnostics/+groups/registerPlantHydraulics.m](../src/+diagnostics/+groups/registerPlantHydraulics.m)

所有字段单位 `mm s-1 MPa-1`，维度 `time`。

**命名约定 v2.1**: `k_*` = 导度 / `c_*` = 容量（预留植被含水）/ `w_*` = 含水 / `q_*` = 通量。

**两种导度定义**: `_app` (apparent = Tr/Δψ 事后) vs `_diff` (differential = dT/dψ ODE 用)。

| 字段 | 旧 CSV 名 | 物理含义 |
|---|---|---|
| `k_plant_sun_app` | KplantAppSun | 阳叶视导度 = Tp_sun / (ψ_stem - ψ_leaf_sun) |
| `k_plant_sha_app` | KplantAppSha | 阴叶视导度 |
| `k_plant_canopy_app` | KplantAppTotal | 全冠视导度 (4q 组合) |
| `k_plant_sun_diff` | KplantDiffSun | 阳叶切线导度 (VulnerabilityCurve) |
| `k_plant_sha_diff` | KplantDiffSha | 阴叶切线导度 |
| `k_plant_canopy_diff` | KplantDiffFluxWeighted | 通量权重切线 |
| `k_soil_root` | KsoilrootAppTotal | 土→根视导度 |
| `k_soil_leaf` | KwholeAppTotal | 土→叶整路径视导度 |
| `gsr_node_sum` | GsrTotal | Σ gsr_j（校验：= 1/R_total） |

**预留位** (将来植被含水动态):
- `c_leaf` / `c_stem` / `c_root` — 叶/茎/根容量
- `w_leaf` / `w_stem` / `w_root` — 含水量
- `q_stem_to_leaf` 等 — 段内瞬时通量

---

## Group 8: surface_water (20 fields)

**输出**: `surface_water.nc`
**注册文件**: [src/+diagnostics/+groups/registerSurfaceWater.m](../src/+diagnostics/+groups/registerSurfaceWater.m)

每对水量都有 `model unit` (cm/timestep) + `display unit` (mm h-1) 两个字段。

| 字段 (model) | 字段 (mm h-1) | 旧 CSV 名（旧/新） |
|---|---|---|
| `precipitation_flux` | `precipitation_amount` | PrecipitationFlux / PrecipitationAmount_mm |
| `surface_liquid_water_supply_flux` | `surface_liquid_water_supply_amount` | SurfaceLiquidWaterSupplyFlux / SurfaceLiquidWaterSupplyAmount_mm |
| `wet_canopy_evaporation_flux` | `wet_canopy_evaporation_amount` | WetCanopyEvaporationFlux / WetCanopyEvaporationAmount_mm |
| `soil_top_water_flux` | `soil_top_water_amount` | SoilTopWaterFlux / SoilTopWaterAmount_mm |
| `soil_bottom_water_flux` | `soil_bottom_water_amount` | SoilBottomWaterFlux / SoilBottomWaterAmount_mm |
| `surface_runoff_flux` | `surface_runoff_amount` | SurfaceRunoffFlux / SurfaceRunoffAmount_mm |

有限积水边界另输出 8 个状态/互补诊断：

| 字段 | 单位 | 含义 |
|---|---|---|
| `surface_pressure_head` | cm | 接受子步的顶节点压力头时间均值 |
| `surface_water_storage` | cm | 接受后的液态表面积水深度 |
| `surface_water_storage_change` | mm | 接受外步的表面储水变化 |
| `surface_water_available_flux` | cm s-1 | 液态供给加步初 `DSTOR/dt` |
| `surface_boundary_mode` | 1 | 末个接受子步：0 供给、1 积水、2 满库溢流、3 残余含水量干燥限蒸发 |
| `surface_mass_residual_max` | cm | 外步内最大表面水库质量残差绝对值 |
| `surface_lower_complementarity_residual_max` | cm2 | 最大 `abs(DSTOR*(DSTOR-psi_top))` |
| `surface_upper_complementarity_residual_max` | cm2 s-1 | 最大 `abs(reaction*(DSTMAX-DSTOR))` |

Use `soil_top_water_flux/amount` as the realized QMT-derived soil-interface
exchange. Potential/actual infiltration aliases are not public outputs.

CF standard_name 用例：`precipitation_rate`、`surface_runoff_amount`。

---

## Group 9: evapotranspiration (6 fields)

**输出**: `evapotranspiration.nc`
**注册文件**: [src/+diagnostics/+groups/registerEvapotranspiration.m](../src/+diagnostics/+groups/registerEvapotranspiration.m)

| 字段 | 旧 CSV 名 | 单位 | 维度 | standard_name |
|---|---|---|---|---|
| `soil_evaporation_flux` | SoilEvaporationFlux | cm s-1 | time | water_evaporation_flux_from_soil |
| `soil_evaporation_amount` | SoilEvaporationAmount_mm | mm | time | accepted-step amount |
| `soil_evaporation_demand_flux` | SoilEvaporationDemandFlux | cm s-1 | time | upstream SCOPE demand |
| `soil_evaporation_demand_amount` | SoilEvaporationDemandAmount_mm | mm | time | requested-step amount |
| `transpiration_flux` | TranspirationFlux | cm s-1 | time | transpiration_flux |
| `transpiration_amount` | TranspirationAmount_mm | mm | time | accepted-step amount |

---

## Group 10: root_uptake (9 字段)

**输出**: `root_uptake.nc`
**注册文件**: [src/+diagnostics/+groups/registerRootUptake.m](../src/+diagnostics/+groups/registerRootUptake.m)

| 字段 | 旧 CSV 名 | 单位 | 维度 |
|---|---|---|---|
| `water_uptake_profile` | WaterUptake | cm s-1 | time × layer |
| `rwu_pos_profile` | RWU_pos | mm h-1 | time × layer |
| `rwu_neg_profile` | RWU_neg | mm h-1 | time × layer |
| `hydraulic_lift_pos` | HydraulicLift_F_pos_mmh | mm h-1 | time |
| `hydraulic_lift_net` | HydraulicLift_F_net_mmh | mm h-1 | time |
| `hydraulic_lift_intensity` | HydraulicLift_intensity | 1 | time |
| `root_sink_srt_pos` | RootSinkSrtPos_mm | mm h-1 | time |
| `root_sink_srt_neg` | RootSinkSrtNeg_mm | mm h-1 | time |
| `root_sink_srt_net` | RootSinkSrtNet_mm | mm h-1 | time |

`root_sink_minus_transpiration` → 归到 `conservation.nc`


---
## Group 11: dual_porosity (3 字段, lDualPor=1)

**输出**: `dual_porosity.nc`
**注册文件**: [src/+diagnostics/+groups/registerDualPorosity.m](../src/+diagnostics/+groups/registerDualPorosity.m)

| 字段 | 旧 CSV 名 | 单位 | 维度 |
|---|---|---|---|
| `mobile_water` | MobileWater | m3 m-3 | time × node |
| `immobile_water` | ImmobileWater | m3 m-3 | time × node |
| `sink_immobile` | SinkIm | s-1 | time × node |

---

## Group 12: conservation (17 fields, cross-group collection)

**Output**: `conservation.nc`
**Registration**: fields with `group='conservation'` are collected from multiple `register*.m` files.

| Field | Legacy CSV name | Unit | Dimension | Meaning |
|---|---|---|---|---|
| `water_balance_residual` | WaterBalanceResidual_mm | mm | time | accepted soil top input + bottom input - integrated Srt - dS |
| `nh4_cum_residual` | NH4_cum_residual | ug cm-2 | time | NH4 cumulative mass residual |
| `no3_cum_residual` | NO3_cum_residual | ug cm-2 | time | NO3 cumulative mass residual |
| `canopy_water_residual` | canopy_water_residual_mm | mm | time | Canopy water residual |
| `root_sink_minus_transpiration` | RootSinkMinusTranspiration_mm | mm h-1 | time | Root sink minus transpiration closure |
| `plant_c_residual` | PlantC_residual | g m-2 | time | Plant C conservation residual |
| `plant_n_residual` | PlantN_residual | g m-2 | time | Plant N conservation residual |
| `energy_balance_residual` | EnergyBalanceResidual_W | W m-2 | time | Rn - LE - H - G energy closure |
| `soil_energy_balance_residual` | SoilEnergyBalanceResidual_J_cm2 | J cm-2 | time | Exact accepted soil-column internal-energy residual |
| `soil_dry_air_mass_balance_residual` | SoilDryAirMassBalanceResidual_g_cm2 | g cm-2 | time | Exact accepted dry-air column mass residual |
| `surface_soil_liquid_input` | SurfaceSoilLiquidInput_mm | mm | time | Post-canopy/post-snow liquid input |
| `surface_soil_evaporation_loss` | SurfaceSoilEvaporationLoss_mm | mm | time | Soil-surface evaporation loss |
| `surface_soil_runoff_loss` | SurfaceSoilRunoffLoss_mm | mm | time | Surface overflow runoff |
| `surface_soil_bottom_loss` | SurfaceSoilBottomLoss_mm | mm | time | Signed lower-boundary loss |
| `surface_soil_root_sink` | SurfaceSoilRootSink_mm | mm | time | Signed root sink |
| `surface_soil_storage_change` | SurfaceSoilStorageChange_mm | mm | time | Surface plus soil storage change |
| `surface_soil_water_balance_residual` | SurfaceSoilWaterBalanceResidual_mm | mm | time | Combined surface-reservoir plus soil-column residual |

Opening `conservation.nc` should expose all closure diagnostics in one place.
---


## Group 13: numerics (1 字段)

**输出**: `numerics.nc`

| 字段 | 旧 CSV 名 | 单位 | 维度 | 含义 |
|---|---|---|---|---|
| `step_failed` | StepFailed | 1 | time | 0=OK, 1=Picard 不收敛 |

---

## Group 14: soil_oxygen (5 字段, lNefc=1) — 2026-05-15 注册, 2026-05-16 主循环写入修复

**输出**: `soil_oxygen.nc`
**注册文件**: [src/+diagnostics/+groups/registerSoilOxygen.m](../src/+diagnostics/+groups/registerSoilOxygen.m)

O2 是 N 模块的耦合输入（控制硝化 / 反硝化），单独 group 让用户能直接看"什么时候哪个层缺氧"。

**2026-05-16 修复 (Bug A)**: 5 字段在 v33 IO 重构时只注册不写入, 导致 yearfull_v3 NC 里 100% NaN.
COMPASS.m 主循环已补 5 行 `reg.write`, 同时 `solveNitrogen.m` 把 `R_O2_sink` 暴露给主循环
(`NitrogenVariables.R_O2_sink`). 100KT smoke 验证: 5 字段 100% 有效, aerobic_fraction 冬季均值 ~0.38.

| 字段 | 单位 | 维度 | 来源 |
|---|---|---|---|
| `o2_concentration` | mol m-3 | time × node | `NitrogenVariables.C_O2` |
| `aerobic_fraction` | 1 | time × node | `f_aerobic` (NL×2) |
| `anaerobic_fraction` | 1 | time × node | `f_anaerobic` (= 1 - f_aerobic) |
| `soil_gas_phase_fraction` | m3 m-3 | time × node | `Theta_g = POR - θ_LL - θ_II` |
| `o2_sink_total` | mol m-3 s-1 | time × node | `R_O2_sink` (硝化 + 矿化耗 O2) |

---

## Group 15: radiation (12 字段) — 2026-05-15 新增

**输出**: `radiation.nc`
**注册文件**: [src/+diagnostics/+groups/registersurface.radiation.m](../src/+diagnostics/+groups/registersurface.radiation.m)

从原 `Scope_output.rad_*` + 部分 `fluxes_*` 字段迁过来。Paper 报告的辐射收支 + PAR + SIF。

| 字段 | 旧 SCOPE 名 | 单位 |
|---|---|---|
| `net_radiation_total` | fluxes.Rntot | W m-2 |
| `net_radiation_canopy` | fluxes.Rnctot | W m-2 |
| `net_radiation_soil` | fluxes.Rnstot | W m-2 |
| `longwave_outgoing` | rad.Lout | W m-2 |
| `longwave_thermal_emission` | rad.Loutt | W m-2 |
| `shortwave_outgoing_hemis` | rad.Eout | W m-2 |
| `shortwave_outgoing_obs` | rad.Eouto | W m-2 |
| `par_at_canopy_top` | rad.PAR × 1e6 | umol m-2 s-1 |
| `apar_canopy` | fluxes.aPAR × 1e6 | umol m-2 s-1 |
| `apar_cab` | fluxes.aPAR_Cab × 1e6 | umol m-2 s-1 |
| `apar_canopy_wm2` | fluxes.aPAR_Wm2 | W m-2 |
| `fluorescence_outgoing` | rad.Eoutf | W m-2 |

---

## Group 16: energy_balance (13 字段) — 2026-05-15 新增

**输出**: `energy_balance.nc`
**注册文件**: [src/+diagnostics/+groups/registerEnergyBalance.m](../src/+diagnostics/+groups/registerEnergyBalance.m)

地表能量平衡四元组 + 冠层 / 土壤拆分 + 4-quadrant 阴阳叶细分。

| 字段 | 旧名 | 单位 | 备注 |
|---|---|---|---|
| `land_surface_upward_latent_heat_flux` | fluxes.lEtot | W m-2 | Positive upward; canopy plus ground |
| `land_surface_upward_sensible_heat_flux` | fluxes.Htot | W m-2 | Positive upward; canopy plus ground |
| `soil_top_downward_energy_flux` | fluxes.soil_top_downward_energy_flux_W_m2 | W m-2 | Positive downward into shared soil |
| `latent_heat_canopy` | fluxes.lEctot | W m-2 | LE_canopy |
| `latent_heat_canopy_wet` | fluxes.lE_wet_canopy | W m-2 | 冠层截留湿叶蒸发 |
| `latent_heat_canopy_dry` | fluxes.lE_dry_canopy | W m-2 | 冠层干叶蒸腾 |
| `ground_upward_latent_heat_flux` | fluxes.ground_upward_latent_heat_flux_W_m2 | W m-2 | Bare plus snow, positive upward |
| `sensible_heat_canopy` | fluxes.Hctot | W m-2 | H 冠 |
| `ground_upward_sensible_heat_flux` | fluxes.ground_upward_sensible_heat_flux_W_m2 | W m-2 | Bare plus snow, positive upward |
| `latent_heat_wet_sun` | lE_wet_sun | W m-2 | **4q 拆分**：阳叶截留 LE |
| `latent_heat_wet_sha` | lE_wet_sha | W m-2 | 阴叶截留 LE (含凝露负值) |
| `latent_heat_dry_sun` | lE_dry_sun | W m-2 | 阳叶干蒸 LE |
| `latent_heat_dry_sha` | lE_dry_sha | W m-2 | 阴叶干蒸 LE |

**`energy_balance_residual`** 归到 `conservation.nc`（Rn − LE − H − G 应 ≈ 0）。

---

## Group 17: carbon_fluxes (5 字段) — 2026-05-15 新增

**输出**: `carbon_fluxes.nc`
**注册文件**: [src/+diagnostics/+groups/registerCarbonFluxes.m](../src/+diagnostics/+groups/registerCarbonFluxes.m)

光合 / 呼吸 / GPP / NEE 全套。

| 字段 | 旧名 | 单位 | standard_name |
|---|---|---|---|
| `gross_primary_productivity` | fluxes.GPP | umol m-2 s-1 | gross_primary_productivity_of_biomass_expressed_as_carbon |
| `net_ecosystem_exchange` | fluxes.NEE | umol m-2 s-1 | net_primary_productivity_of_biomass_expressed_as_carbon |
| `assimilation_canopy` | fluxes.Actot | umol m-2 s-1 | — |
| `assimilation_total` | fluxes.Atot | umol m-2 s-1 | — |
| `respiration_ecosystem` | fluxes.Resp（仅 lPlantCN=1 且 lNefc=1） | umol m-2 s-1 | ecosystem_respiration_carbon_flux；链关闭时为缺失值，不把 SCOPE 占位 0 当成生态系统呼吸 |

---

## plant_nitrogen — 自养呼吸口径 (potential vs actual) — 2026-06-15 新增

**输出**: `plant_nitrogen.nc` + 同名 CSV ｜ **注册文件**: [src/+diagnostics/+groups/registerPlantNitrogen.m](../src/+diagnostics/+groups/registerPlantNitrogen.m)

植株自养呼吸有两套口径，逐小时输出（`lNefc=1`）。`θ = resp_theta_actual = R_m_final / R_m_potential`（含 nodebt cap）；器官 actual = θ × potential。**actual 口径与 `NEE` / `Reco` / 植物碳守恒闭合；potential 仅作诊断**（暴露被 NSC 节流掩盖的呼吸量）。

| 字段 | 口径 | 单位 | 说明 |
|---|---|---|---|
| `autotrophic_respiration` | **actual** | kg m-2 s-1 | 整株 Ra = R_m+R_g 节流后（进 NEE/Reco 的真实 CO2 通量） |
| `r_a_potential` / `cum_r_a_potential` | potential | g m-2 | 节流前 Ra = R_m_pot+R_g，诊断用 |
| `net_primary_productivity` / `npp_net` | **actual** | kg m-2 s-1 | GPP − R_a_actual（clamped / unclamped）；真 CUE |
| `r_m_leaf_actual` / `r_m_stem_actual` / `r_m_root_actual` / `r_m_livecroot_actual` / `r_m_fruit_actual` | **actual** | kg m-2 s-1 | θ × 对应 potential 器官维持呼吸 |
| `r_g_actual` | **actual** | kg m-2 s-1 | 生长呼吸（不节流） |
| `ra_above_actual` / `ra_below_actual` | **actual** | kg m-2 s-1 | 地上(leaf+stem+fruit Rm + R_g) / 地下(root+livecroot Rm) |
| `resp_theta_actual` | — | 1 | 有效节流因子 = actual R_m / potential R_m（解释上面 actual 怎么来；>1 可能出现在 substitution 方案如 jones） |
| `maintenance_respiration_{leaf,root,stem,livecroot,fruit}` | potential | kg m-2 s-1 | 节流前器官维持呼吸（对照 actual 看节流强度） |

---

## Group 18: canopy_aerodynamics (18 字段) — 2026-05-15 新增

**输出**: `canopy_aerodynamics.nc`
**注册文件**: [src/+diagnostics/+groups/registerCanopyAerodynamics.m](../src/+diagnostics/+groups/registerCanopyAerodynamics.m)

气动阻抗 + 表面温度 + 友速 + 粗糙度长度。从 `Scope_output.thermal_*` 拆出来。dims 全部 `time`。

### 18.1 温度与空气状态 (4 字段)

| 字段 | 旧 SCOPE 名 | 单位 | standard_name |
|---|---|---|---|
| `air_temperature` | thermal.Ta | degC | air_temperature |
| `soil_surface_temperature_avg` | thermal.Tsave | degC | — |
| `canopy_temperature_avg` | thermal.Tcave | degC | — |
| `air_density` | thermal.rhoa | kg m-3 | air_density |

### 18.2 气动阻抗 (7 字段, 单位 s m-1)

| 字段 | 旧 SCOPE 名 | 物理含义 |
|---|---|---|
| `resistance_aa` | thermal.raa | 冠上方气动阻抗 |
| `resistance_ah` | thermal.rah | 热量气动阻抗 |
| `resistance_aw` | thermal.raw | 水汽气动阻抗 |
| `resistance_canopy_bulk` | thermal.rawc | 冠层 bulk 阻抗 |
| `resistance_soil_aero` | thermal.raws | 土壤气动阻抗 |
| `resistance_soil_surface` | thermal.rss | 土壤表面阻抗 |
| `resistance_leaf_boundary` | thermal.rbc | 叶边界层阻抗 |

### 18.3 风/粗糙度 (4 字段)

| 字段 | 旧 SCOPE 名 | 单位 | standard_name |
|---|---|---|---|
| `friction_velocity` | thermal.ustar | m s-1 | wind_friction_velocity_at_sea_level |
| `roughness_length_momentum` | thermal.zo | m | — |
| `roughness_length_heat` | thermal.z0h | m | — |
| `displacement_height` | thermal.d | m | — |

### 18.4 林分结构 (3 字段)

| 字段 | 旧 SCOPE 名 | 单位 | 物理含义 |
|---|---|---|---|
| `stem_area_index` | thermal.SAI | m2 m-2 | 茎面积指数 |
| `plant_area_index` | thermal.L_plus_S | m2 m-2 | 植物面积指数 (LAI + SAI) |
| `soil_surface_saturation` | thermal.Sat_surf | 1 | ebal 土表阻抗用 |

---

## Group 19: canopy_profiles (14 字段) — 2026-05-15 新增

**输出**: `canopy_profiles.nc`
**注册文件**: [src/+diagnostics/+groups/registerCanopyProfiles.m](../src/+diagnostics/+groups/registerCanopyProfiles.m)

冠层 30-layer 垂直剖面（**新增 `canopy_layer` 维度，size = 30**）。从 `Scope_output.profiles_*` 拆出来。除 `canopy_avg_temperature_3d` 维度为 `time` 外，其他字段 dims 全部 `time × canopy_layer`。

### 19.1 冠层微气候 (4 字段)

| 字段 | 旧 SCOPE 名 | 单位 | 物理含义 |
|---|---|---|---|
| `air_temperature_layer` | profiles.Ta_layer | degC | 每层空气温度 |
| `vapor_pressure_layer` | profiles.ea_layer | hPa | 每层水汽压 |
| `co2_layer` | profiles.CO2_layer | umol mol-1 | 每层 CO₂ 混合比 |
| `eddy_diffusivity_layer` | profiles.Kh_layer | m2 s-1 | 每层湍流扩散系数 |

### 19.2 叶/冠温度 (2 字段)

| 字段 | 旧 SCOPE 名 | 单位 | 物理含义 |
|---|---|---|---|
| `leaf_temperature_layer` | profiles.Tleaf_layer | degC | 叶温（每层）|
| `canopy_temperature_layer` | profiles.Tch | degC | 冠层元素温度（每层）|

### 19.3 结构权重 (2 字段)

| 字段 | 旧 SCOPE 名 | 单位 | 物理含义 |
|---|---|---|---|
| `lai_layer` | profiles.LAI_layer | m2 m-2 | LAI 分层 |
| `layer_weight` | profiles.layer_weight | 1 | 垂直加权权重 |

### 19.4 光合 / 通量 (5 字段)

| 字段 | 旧 SCOPE 名 | 单位 | 物理含义 |
|---|---|---|---|
| `net_photosynthesis_layer` | profiles.Pn1d | umol m-2 s-1 | 每层净光合 |
| `net_photosynthesis_cab_layer` | profiles.Pn1d_Cab | umol m-2 s-1 | Pn × Cab 比例 |
| `sensible_heat_layer` | profiles.Hc1d | W m-2 | 每层显热 |
| `latent_heat_layer` | profiles.lEc1d | W m-2 | 每层潜热 |
| `assimilation_layer` | profiles.A1d | umol m-2 s-1 | 每层同化 |

### 19.5 标量摘要 (1 字段)

| 字段 | 旧 SCOPE 名 | 单位 | 维度 | 物理含义 |
|---|---|---|---|---|
| `canopy_avg_temperature_3d` | profiles.Tchave | degC | time | 冠层平均温度（Tchave）|

---

## Group 20: solver_internal (17 字段) — 2026-05-15 新增

**输出**: `solver_internal.nc`
**注册文件**: [src/+diagnostics/+groups/registerSolverInternal.m](../src/+diagnostics/+groups/registerSolverInternal.m)

Picard 内部诊断 — 土壤水/汽/气通量分量 + 梯度 + 元素平均导度。从 `flux_vars.*` 17 个字段全部迁过来。dims 全部 `time × layer`（NL 元素）。

### 20.1 构成核液水通量重构及分量 (4 字段, 单位 cm s-1)

这些量来自接受状态下的热/干空气系数核，可包含冻结吸力头、温度和气压项；正向 = 向下。它们不是 Richards 守恒水方程接受后的层间质量通量，也不能替代 `surface_water.nc` 中的实际顶通量/入渗量。

| 字段 | 旧名 | 物理含义 |
|---|---|---|
| `constitutive_liquid_flux` | flux_vars.QLiquid | 构成核液水通量重构 |
| `constitutive_liquid_flux_h` | flux_vars.QLiquid_h | matric head / 冻结吸力头分量 |
| `constitutive_liquid_flux_t` | flux_vars.QLiquid_T | 温度梯度分量 |
| `constitutive_liquid_flux_a` | flux_vars.QLiquid_A | 空气压梯度分量 |

### 20.2 水汽通量及分量 (4 字段, 单位 cm s-1)

| 字段 | 旧名 | 物理含义 |
|---|---|---|
| `q_vapor` | flux_vars.QVapor | 水汽总通量 (h+T+A) |
| `q_vapor_h` | flux_vars.QVapor_h | h 梯度分量 |
| `q_vapor_t` | flux_vars.QVapor_T | T 梯度分量 |
| `q_vapor_a` | flux_vars.QVapor_A | A 梯度分量 |

### 20.3 干空气通量 (1 字段, lSoilairefc=1)

| 字段 | 旧名 | 单位 | 物理含义 |
|---|---|---|---|
| `q_air` | flux_vars.QAir | cm s-1 | 干空气通量（仅在 `lSoilairefc=1` 时有效）|

### 20.4 梯度 (3 字段)

| 字段 | 旧名 | 单位 | 物理含义 |
|---|---|---|---|
| `gradient_matric_head` | flux_vars.Grad_h | cm cm-1 | 基质势梯度 |
| `gradient_temperature` | flux_vars.Grad_T | K cm-1 | 温度梯度 |
| `gradient_air_pressure` | flux_vars.GradPg / 10 | Pa cm-1 | 内部 10 Pa 压力尺度换算后的空气压梯度 |

### 20.5 元素平均导度 (3 字段)

| 字段 | 旧名 | 单位 | 物理含义 |
|---|---|---|---|
| `cond_hydraulic_avg` | flux_vars.CondHydraulicAvg | cm s-1 | 元素平均水力导度 |
| `cond_thermal_avg` | flux_vars.CondThermalAvg | cm K-1 s-1 | 元素平均水的热致导度 |
| `cond_vapor_avg` | flux_vars.CondVaporAvg | cm s-1 | 元素平均水汽导度 |

### 20.6 饱和诊断 (2 字段)

| 字段 | 旧名 | 单位 | 物理含义 |
|---|---|---|---|
| `saturation_water_layer` | flux_vars.SatWater | 1 | 元素水有效饱和度 |
| `saturation_air_layer` | flux_vars.SatAir | 1 | 元素气有效饱和度 |

---

## 全局 NetCDF Attributes (每个 .nc 文件都带)

```text
:Conventions    = "CF-1.10"
:title          = "COMPASS output: <group>"
:source         = "COMPASS v33 (Zhu Jizixing, Northwest A&F University)"
:run_tag        = "<config 提取>"
:start_time     = "2021-01-01T00:00"
:end_time       = "2022-01-01T00:00"
:delt_t_seconds = 3600
:cfg_lFreezing  = 0|1
:cfg_lNefc      = 0|1
:cfg_lSoilairefc= 0|1
:institution    = "Northwest A&F University"
:contact        = "Zhu Jizixing <zjzx3100362274@gmail.com>, <3100362274@nwafu.edu.cn>"
:history        = "Created 2026-05-15 16:23:01 by diagnostics.flushToNetCDF"
:references     = "Zhu et al. 2026 GMD (in prep)"
```

---

## 加新变量的标准流程

1. **想清楚归哪个 group**（参考"五大变量分类"，state/aux/forcing/diag）
2. 打开对应 `register*.m`，加一段 `reg.register(struct(...))`
3. 在主循环里加**一行** `reg.write('your_name', value, KT);`
4. 更新本 md（在对应 group 表格加一行）
5. 跑 `tests/test_diagnostics_smoke.m` 验证管道不崩
6. 跑短窗口 COMPASS 验证 .nc 文件能 xarray 读

完成后**无需碰 writeCompassOutputs / flushToCsv**，注册中心自动处理。

---

## 兼容期 (1 个月过渡期, 至 ~2026-06-15)

`flushToCsv` 同时写出旧 CSV 命名，旧 plot 脚本（PlotResults/*.m）继续工作。

`Sim_output_<LegacyName>.csv` 的旧名映射定义在 [src/+diagnostics/flushToCsv.m](../src/+diagnostics/flushToCsv.m) 的 `legacyCaseFor()` 函数里（~80 行查找表）。

迁移完成后下个 minor version 砍掉 `flushToCsv`，所有 plot 改读 .nc。

---

## 参考

- `docs/netcdf_output_design_decisions.md` — 8 个设计决定
- `src/+diagnostics/README.md` — 包设计文档
- ClimaLand `define_diagnostics.jl` — 设计灵感来源
- CF Convention 1.10 — http://cfconventions.org/
- UDUNITS-2 — https://docs.unidata.ucar.edu/udunits/

---

## 2026-05-16 增量更新（N MM 改造 + O2 Bug A 修 + 真实 depth_node + 学校名修正）

| 改动 | 文件 | 内容 |
|---|---|---|
| **soil_nitrogen.nc**: +2 字段 | `registerSoilNitrogen.m` + COMPASS.m | `nit_mm_factor`, `nh4_substrate_total` (DAISY MM 诊断), Group 5 → 30 写入 |
| **soil_oxygen.nc**: Bug A 修 | COMPASS.m 主循环 + solveNitrogen.m | 5 字段从 100% NaN → 100% 有数据 (yearfull_v3 之前漏 reg.write) |
| **真实 depth_node** | COMPASS.m runMeta | NC 不再用 `linspace(0,1000,NN)` fallback, 改为 `cumsum(flipud(DeltZ))` 真实 Dtrmn_Z |
| **N 化学物理改造** | `nitrification.m` + `denitrification.m` + `assembleSoluteParameters.m` + `solveTransportStep.m` + `SoilProperties.m` + `setInitialNitrogen.m` + soilproperty.xlsx | 线性 K_nit → DAISY 显式源汇 MM (Snit 进 G1, 反硝化也改显式), 100KT smoke 守恒 0% |
| **学校名修正** | flushToNetCDF.m, output_catalog.md, netcdf_output_design_decisions.md | Tsinghua → Northwest A&F University, contact 加 `3100362274@nwafu.edu.cn` |

### 全年验证状态（进行中）
- `data/output_v33_nMM_yearfull/` (8760 KT) MATLAB 后台跑中, 23% 推进
- 验收: NH4 浓度全年中位 5-10 mg/kg, 硝化年通量 100-200 kg N/ha, N 守恒 < 1%

### 全局属性更新

```
:Conventions    = "CF-1.10"
:title          = "COMPASS output: <group>"
:source         = "COMPASS v33 (Zhu Jizixing, Northwest A&F University)"
:institution    = "Northwest A&F University"
:contact        = "Zhu Jizixing <zjzx3100362274@gmail.com>, <3100362274@nwafu.edu.cn>"
:run_tag        = "<config 提取>"
:start_time     = "2021-01-01T00:00"
:end_time       = "2022-01-01T00:00"
:delt_t_seconds = 3600
:cfg_lFreezing  = 0|1
:cfg_lNefc      = 0|1
:cfg_lSoilairefc= 0|1
:history        = "Created <datetime> by diagnostics.flushToNetCDF"
```

**注意**: 当前正在跑的 nMM_yearfull 任务在学校名修复**之前**启动, 仍会带"Tsinghua". 跑完后下次新跑就是 NWAFU.
