# COMPASS NetCDF 输出设计决定 (2026-05-15)

8 个关键选择已与用户确认。下一波重构按此实施。

## 决定记录

| 问题 | 选择 | 简述 |
|---|---|---|
| Q1 分文件粒度 | **B** | 按物理 group 分 ~13 个 .nc（ClimaLand 风格） |
| Q2 剖面对齐维度 | **B** | 统一插值到 `node` 维（NL×2 → NN 端点平均） |
| Q3 时间轴 | **A** | 外步 KT 末 8760 个时间点，CF datetime |
| Q4 深度方向 | **A** | 顶 → 底，depth_node=[0..1000] cm，CF `positive="down"` |
| Q5 变量命名 | **C** | 短名 (`soil_moisture`) + `standard_name` 双标 |
| Q6 单位 | **CF UDUNITS-2** | `m3 m-3`、`cm s-1`、`ug cm-2 s-1` 等 |
| Q7 浮点精度 | **C** | `single` (float32) + NetCDF4 deflate=5 压缩 |
| Q8 MATLAB API | **高层** | `nccreate / ncwrite / ncwriteatt` (不用 group) |

## 13 个 group

```
1. soil_water.nc          土壤水状态 (always)
2. soil_thermal.nc        土壤温度 (cfg.lThmrlefc=1)
3. soil_freezing.nc       冻土 (cfg.lFreezing=1)
4. soil_air.nc            气相 (cfg.lSoilairefc=1)
5. soil_nitrogen.nc       N 模块 (cfg.lNefc=1)
6. plant_canopy.nc        冠层 (LAI/leaf potential/fwet)
7. plant_hydraulics.nc    植物水力 (Kplant/Ksoilroot/RWU)
8. surface_water.nc       表层水通量 (precip/throughfall/infiltration/runoff)
9. evapotranspiration.nc  ET 拆分
10. conservation.nc       守恒残差 (water + N)
11. radiation.nc          辐射 (从 Scope_output)
12. numerics.nc           数值诊断 (StepFailed/PicardIter)
13. solver_internal.nc    Picard 内部 (从 flux_vars)
```

## Global attributes (每文件必带)

```
:Conventions    = "CF-1.10"
:title          = "<group description>"
:source         = "COMPASS v33"
:run_tag        = "<from config>"
:start_time     = "<datetime>"
:end_time       = "<datetime>"
:cfg_lFreezing  = 0|1
:cfg_lNefc      = 0|1
:cfg_lSoilairefc= 0|1
:institution    = "Northwest A&F University"
:contact        = "Zhu Jizixing <zjzx3100362274@gmail.com>, <3100362274@nwafu.edu.cn>"
:history        = "Created <datetime> by COMPASS.m"
:references     = "Zhu et al. 2026 GMD (in prep)"
```

## 实施工作量

- `src/+diagnostics/Registry.m` — 注册中心（变量元数据）4 h
- `src/+diagnostics/writeNetCDF.m` — 主写函数 3 h
- 99 个旧字段注册 4 h
- 旧 CSV 兼容写出（1 个月过渡期） 1 h
- plot 脚本改读 xarray 3 h
- 验证测试 2 h
**总计**：~17 小时 = 2 天

