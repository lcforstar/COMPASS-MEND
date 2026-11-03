# COMPASS Final Model Version Checklist

日期：2026-06-20

用途：在把当前 MATLAB 版本固定为可复现参考版本前，逐项检查主循环状态流、输入输出口径、数值稳定性和核心科学闭合。本文档接在 `docs/main_loop_state_flow_audit.md` 后使用；不替代状态流审计，也不讨论继续提速。

## 0. 冻结检查范围

最终检查先冻结“要验证的版本”，否则后续 NetCDF、图件和统计口径无法追溯。

必须记录：

- `git status --short`，明确哪些改动属于本次候选版本，哪些是无关工作区改动；
- `git rev-parse HEAD`；
- 主要配置文件快照：`data/config_main.txt`、`data/config_environment.txt`、`data/config_grid.txt`、`data/config_vegetation.txt`；
- 运行时覆盖项：`STEMMUS_START_TIME`、`STEMMUS_END_TIME`、`STEMMUS_NMSRMN`、`STEMMUS_OUTPUT`、`STEMMUS_DIAG_INTERVAL`；
- 当前数值 kernel：`HydraulicKernel`、`MatricKernel`。

通过标准：

- 能明确回答“最终检查针对哪一个 checkout、哪一套 config、哪一个输出目录”；
- 运行输出目录不复用旧目录，避免混入历史 NetCDF。

## 1. 主循环状态流 Gate

依据 `docs/main_loop_state_flow_audit.md`，先完成主循环状态清单，再检查公式。

交付物：

- `docs/main_loop_state_inventory.md`
- 一页流程图或表：`forcing -> canopy -> plant C/N -> root sink -> soil substeps -> nitrogen -> diagnostics`

必须覆盖的状态块：

- 土壤水：`SoilVariables.hh`、`h`、`Theta_LL`、`Theta_UU`、`Theta_II`、`KfL_h`
- 土壤热：`SoilVariables.TT`、`T`、`EnergyVars` 中热通量和热参数
- 土壤气：`AirVars`、`P_gg`
- 土壤氮：`NitrogenVariables`、`NitrogenFlux`
- 冠层能量和光合：`fluxes`、`thermal`、`profiles`、`rad`
- 植物水力：`plantHydro` 和可选动态植物水状态
- 植物 C/N：器官 C/N 池、NSC、allocation、actual/potential respiration
- 地表水和边界：`BoundaryCondition`、`DSTOR`、`RS`、入渗/径流项
- 纯诊断：`diagnostics.Registry` 和局部输出 history arrays

通过标准：

- 每个状态字段都有 owner、更新位置、rollback 策略和输出 group；
- 能区分 prognostic state、trial-only state、diagnostic-only state；
- 明确哪些字段是一小时滞后或日尺度更新，而不是误判为 bug。

## 2. Rollback 和 Adaptive Substep Gate

重点文件：

- `src/COMPASS.m`
- `src/+init/defineSoilVariables.m`
- `src/UpdateSoilWaterContent.m`
- `src/+soil/+water/solveSoilMoistureBalance.m`
- `src/+soil/+thermal/solveThermalEquations.m`
- `src/+soil/+gas/solveDryAirEquations.m`
- `src/+soil/+water/+dualporosity/dualPorExchange.m`

检查项：

- trial substep 失败后，水、热、气、边界和地表水 trial 项恢复到同一时间点；
- outer-step fatal rollback 时，跳过当前小时输出装配和写入，不产生半步诊断；
- Picard 内水、热、气使用同一个当前 guess，不混用接受前后的状态；
- 根吸水、双孔隙 exchange 等外部源项是否在 substep 内保持常数，要有明确物理解释；
- rollback 计数、Picard 迭代分布和 fatal failure 数量写入或可从日志追踪。

通过标准：

- 3 天冬季/冻融窗口没有 fatal rollback；
- 10 天困难窗口没有 fatal rollback；
- rollback 次数如有变化，必须能由步长、收敛阈值或物理开关解释；
- 被拒绝 substep 不改变最终 NetCDF 中的 prognostic state。

## 3. 输入 Schema Gate

重点文件：

- `run_stemmus_initialize.m`
- `src/+core/ModelSettings.m`
- `src/+core/EnvironmentSettings.m`
- `src/+core/VegetationData.m`
- `data/config_*.txt`

交付物：

- `docs/input_schema.md`

检查项：

- config 字段按 site 常数、forcing 时间序列、PFT/植被参数、管理措施、诊断开关分组；
- 每个字段列单位、默认值、允许范围和读取位置；
- `HydraulicKernel`、`MatricKernel` 等运行时选择来自 config，而不是临时改代码；
- 长武站点模板可以复制为新站点，而不需要改主循环。

通过标准：

- 从一个新输出目录的 metadata 可以追溯到完整配置；
- 任何关键参数变化都能在 config diff 或 manifest 中看到。

## 4. 输出 Schema Gate

重点文件：

- `src/+diagnostics/Registry.m`
- `src/+diagnostics/registerAllDiagnostics.m`
- `src/+diagnostics/+groups/*.m`
- `src/+diagnostics/flushToNetCDF.m`
- `docs/output_catalog.md`

建议先跑：

```powershell
& 'F:\Matlab\bin\matlab.exe' -batch "cd('F:\COMPASSV1.0'); addpath('src'); addpath('tests'); test_diagnostics_full_registry"
& 'F:\Matlab\bin\matlab.exe' -batch "cd('F:\COMPASSV1.0'); addpath('src'); addpath('tests'); test_diagnostics_smoke"
```

交付物：

- `docs/output_schema_audit.md`

检查项：

- registry 中每个变量都有单位、维度、source code 和物理含义；
- `docs/output_catalog.md` 与 `src/+diagnostics/+groups/register*.m` 一致；
- group 清单必须以 `src/+diagnostics/registerAllDiagnostics.m` 为当前事实来源，而不是以 `docs/output_catalog.md` 的旧统计数字为准；
- `deep_learning` 是训练/后处理诊断，时间维是 `epoch`，不能混入主循环 `KT` 状态流；
- `online_gs_validation` 只在 `STEMMUS_ONLINE_GS_DIAG=1` 时注册，是临时 write-only 验证输出，不能作为模型物理反馈；
- conservation group 保持小而权威；
- actual/potential respiration、GPP、NPP、Reco、NEE 的口径一致；
- legacy CSV 明确是兼容输出还是正式归档输出。

通过标准：

- full registry test 通过；
- diagnostics smoke test 通过；
- 没有注册但从未写入的正式主循环变量；允许存在明确标注为后处理或临时诊断的例外；
- `output_catalog.md` 的 group 数量、变量数和触发条件已经与当前 registry 对齐，或在 `docs/output_schema_audit.md` 中列为待修复差异；
- NetCDF attributes 没有丢失或类型漂移。

## 5. 短窗 A/B Gate

用途：在长运行前确认 array kernel、诊断输出和状态恢复没有破坏参考结果。

建议窗口：

- 冬季/冻融 3 天；
- 困难水热耦合窗口 10 天；
- 一个代表性生长季片段。

建议输出目录：

- `data/output_finalcheck_3d`
- `data/output_finalcheck_10d`
- `data/output_finalcheck_growth`

推荐命令模板：

```powershell
$env:STEMMUS_START_TIME='2021-01-01T00:00'
$env:STEMMUS_END_TIME='2021-01-04T00:00'
$env:STEMMUS_NMSRMN='301'
$env:STEMMUS_OUTPUT='F:\COMPASSV1.0\data\output_finalcheck_3d'
$env:STEMMUS_DIAG_INTERVAL='100'
& 'F:\Matlab\bin\matlab.exe' -batch "cd('F:\COMPASSV1.0'); addpath('src'); run('src/COMPASS.m');"
```

通过标准：

- 无 NaN/Inf；
- 无 fatal rollback；
- `conservation.nc`、`soil_water.nc`、`soil_thermal.nc`、`carbon_fluxes.nc`、`plant_hydraulics.nc` 写出完整时间轴；
- 水量残差、能量残差、关键剖面与 baseline 差异在预期范围内；
- 任何差异都有对应代码或配置解释。

## 6. 年尺度科学闭合 Gate

用途：确认模型作为最终参考版本不仅能跑完，还能在主要过程上闭合。

最小检查：

- 水量闭合：`conservation.nc/water_balance_residual`
- 能量闭合：`energy_balance.nc`
- 土壤温度剖面：`soil_thermal.nc`
- 土壤含水量剖面：`soil_water.nc`
- 土壤氮闭合：`soil_nitrogen.nc`、`soil_oxygen.nc`
- 植物 C/N 闭合：`plant_nitrogen.nc`、`carbon_fluxes.nc`
- 植物水力：`plant_hydraulics.nc`、`root_uptake.nc`
- 冠层光合和交换：`carbon_fluxes.nc`、`evapotranspiration.nc`、`canopy_profiles.nc`

通过标准：

- 年运行完成且时间轴完整；
- 关键守恒残差有定义清楚的单位和统计口径；
- 植物 C/N 池年内轨迹没有非物理跳变；
- GPP、ET、土壤水、土壤温度、叶水势若有观测，统计窗口和过滤规则明确。

## 7. 多站点验证准备 Gate

只有前六个 gate 通过后，再开始多站点验证。

每个站点必须准备：

- site config package；
- forcing 时间段、缺测处理和筛选规则；
- 观测变量与硬/软验证等级；
- 输出目录和 metadata；
- 指标定义：RMSE、KGE、bias、相关系数或守恒残差；
- 已知 caveat，例如管理措施缺失、观测深度不一致、叶水势采样时刻不一致。

第一批指标：

- 水量闭合；
- 能量残差；
- 土壤温度剖面；
- 土壤含水量剖面；
- ET 或 latent heat；
- GPP；
- 有数据时看 Reco/NEE；
- 有数据时看叶水势或冠层水分胁迫；
- 果园站点看植物 C/N 池是否合理。

## 8. Final Release Packet

最终版本归档时至少包含：

- candidate commit 或 patch diff；
- config snapshot；
- run manifest；
- NetCDF 输出目录；
- `docs/main_loop_state_inventory.md`；
- `docs/input_schema.md`；
- `docs/output_schema_audit.md`；
- 年尺度检查摘要；
- 多站点验证 checklist；
- 已知限制和不在本轮修复的风险。

最终通过标准：

- 代码、配置、输出、文档四者能互相追溯；
- 主循环状态流和 rollback 边界已经说明；
- NetCDF 输出是权威结果源；
- 后续 Julia 转写或多站点验证可以从冻结版本开始，而不是重新猜字段含义。
