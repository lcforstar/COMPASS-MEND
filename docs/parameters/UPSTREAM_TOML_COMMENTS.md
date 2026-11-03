# 原作者通用TOML注释的保留范围

Editor: Wenjie Liu; Date: 2026-10-31.

原始配置及通用注释作者：Jizixing Zhu。长武冬小麦/MEND配置整理：Wenjie Liu。

用户指定参考文件：`COMPASS_MATLAB/data/sites/changwu_apple/site.toml`。
读取时原文件SHA256：`a45eaf056a11fcf590c8c1ba6a1c95e542db85a802f70ba85da71c9bc1b18042`。原文件保持原样，未被写入。

当前site.toml按适用范围保留原作者的单位、模式枚举、边界符号、肥料入口与分组说明。
原作者文件中的果树情景、PFT先验、冠层高度、深根与物候校准背景不迁入当前配置。
保留当前冬小麦/MEND、FAST10、TAN和初始池来源说明；所有399个配置值和section/键顺序不变。

## 恢复的注释块

| 当前section/key | 原文件行号 | 是否适配当前实现 |
|---|---|---|
| `simulation.spatial_mode` | 32–32 | 原说明保留 |
| `simulation.point_context` | 34–35 | 原说明保留 |
| `environment.wind_convert` | 79–79 | 原说明保留 |
| `environment.rho_thermal_soil_dry` | 84–84 | 原说明保留 |
| `environment.Ca` | 87–87 | 原说明保留 |
| `switches.solve_soil_heat` | 94–94 | 原说明保留 |
| `switches.water_retention_model` | 102–102 | 原说明保留 |
| `switches.wetting_heat_method` | 107–107 | 原说明保留 |
| `switches.temperature_driven_liquid_flow` | 114–114 | 原说明保留 |
| `switches.track_total_nitrogen` | 119–119 | 原说明保留 |
| `boundary.water_boundary_file` | 153–172 | 原说明保留 |
| `boundary.surface_water_capacity_cm` | 176–179 | 原说明保留 |
| `boundary.heat_top_kind` | 182–203 | 是 |
| `boundary.dry_air_top_kind` | 209–218 | 原说明保留 |
| `fertilization.pool_entry` | 248–253 | 原说明保留 |
| `vegetation.method` | 258–258 | 原说明保留 |
| `vegetation.LIDFa` | 269–269 | 原说明保留 |
| `vegetation.clumping_index_source` | 272–272 | 原说明保留 |
| `vegetation.canopy_layer_profiles` | 279–279 | 原说明保留 |
| `vegetation.canopy_air_profile_eb` | 281–281 | 原说明保留 |
| `vegetation.vegetation_data_file` | 286–286 | 原说明保留 |
| `vegetation.LAI_profile_mode` | 297–297 | 原说明保留 |
| `vegetation.canopy_storage_per_LAI` | 303–303 | 原说明保留 |
| `vegetation.hydraulic_scheme` | 307–312 | 是 |
| `vegetation.retention_slope_root` | 322–326 | 是 |
| `vegetation.root_dist_mode` | 340–340 | 原说明保留 |
| `vegetation.root_hydraulic_krmax` | 347–347 | 原说明保留 |
| `vegetation.photo_type` | 375–375 | 原说明保留 |
| `vegetation.gs_g0` | 384–385 | 原说明保留 |
| `vegetation.heat_gs_enable` | 389–389 | 原说明保留 |
| `vegetation.heat_gs_max` | 392–392 | 原说明保留 |
| `vegetation.heat_gs_T_on_C` | 395–395 | 原说明保留 |
| `model_library.directional_angles_file` | 437–437 | 原说明保留 |
| `leaf_optical.fluor_beta` | 452–452 | 原说明保留 |
| `soil_optical.BSM_brightness` | 463–463 | 原说明保留 |
| `numerics.NIT` | 473–473 | 原说明保留 |
| `numerics.tol_theta` | 476–476 | 原说明保留 |
| `numerics.soil_interface_closure_enable` | 487–487 | 原说明保留 |
| `numerics.hydraulic_kernel` | 497–497 | 原说明保留 |
| `numerics.extrapolate` | 500–500 | 原说明保留 |
| `numerics.checkpoint_interval` | 505–505 | 原说明保留 |

## 按当前实现修正的说明

- 日尺度辐射输入：true为日均W/m²；false为日累计MJ/m²/day。原模板的MJ/m²/s不适用于当前日资料。
- 原始forcing_start给第一行原始资料定时；模拟开始时间不应被用来重新贴原始行日期。驱动为区间起点，接受态为区间末。
- grid.NS保持“soil.csv参数采样行数”，不采用旧注释的“soil phases”。
- atmospheric热顶为上部20 cm响应匹配后的固定G*；surface_temperature是遗留Dirichlet，季节底温不是当前底边界。
- SCOPE force-restore G是经验/诊断量，不代替耦合土柱的匹配热通量。
- 氮边界当前驱动：顶部fallback为沉降通量mg N/m²/h，底部fallback为水相浓度mg N/L；逐小时强迫优先，不统一误标mol/m³。
- storage_dynamic当前仍共用诊断储水路径，不宣称已实现反馈式动态水容；组织线性持水式中的psi采用米水头。
- 根水/氮开关独立于PlantCN；当前PlantCN关闭时氮汇仍为被动质流，未启用植物需求驱动的主动吸氮。
- 不从原模板引入tol_air_mass_relative新键或0.01闭合目标；当前数值门槛保持原值。

校核依据：当前io/init、surface/radiation、soil水热气氮边界与vegetation/hydraulics活动代码；
热边界另核物理书§2.8.1.7和§5.7。Python非注释代码逐行及完整TOML对象一致；
MATLAB原生解析对象与序列化字节一致。未修改soil.csv、有机初始池、模型源码或后台任务。
