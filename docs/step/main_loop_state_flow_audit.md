# COMPASS 主循环状态流审计计划

日期：2026-06-20

目的：把当前已经优化过的 MATLAB 模型整理成一个稳定、可解释、可验证的参考版本。这个文档不讨论继续提速，重点是：主循环怎么走、状态在哪里更新、输入输出如何固定、以后多站点验证从哪里开始。

## 范围

主入口：

- `src/COMPASS.m`

核心目标：

- 看清每个小时 `KT` 内部发生了什么；
- 区分“真正推进模型的状态量”和“只是输出用的诊断量”；
- 找到 adaptive substep 和 rollback 到底需要保存/恢复哪些状态；
- 在多站点验证或未来转 Julia 前，把输入和输出口径固定下来。

## 推荐顺序

不要一上来改某个公式。先把状态流画清楚，再整理输入输出，最后再逐个模块查公式和科学过程。

1. 主循环骨架
2. 土壤水-热-气内部求解
3. 输入 schema
4. 输出 schema
5. 冠层能量、辐射和光合
6. 植物水力和根吸水
7. 植物 C/N、物候、周转、采收和修剪
8. 多站点验证

## 1. 主循环骨架

目标：把模型看成一串“状态转换”，而不是一堆散函数。

`src/COMPASS.m` 当前大体结构：

1. 初始化
   - `run_stemmus_initialize.m`
   - `init.setVanGenuchtenParameters`
   - `init.defineSoilVariables`
   - `StartInit`
   - 滞回和双孔隙初始化
   - 诊断输出 registry 初始化
   - 边界条件初始化
   - 植物水力初始化

2. 小时外循环
   - 为当前 `KT` 选择 forcing 和 SCOPE/canopy 输入；
   - 更新 LAI、冠层高度、LUNA、叶片生化参数；
   - 计算辐射、冠层湿润、截留等输入；
   - 调 `ebal_clm45` 或 `ebal_clm45_profile` 计算冠层能量和光合；
   - 调 `runPlantCNStep` 推进植物 C/N；
   - 在土壤 Picard 前计算根吸水源项；
   - 进入 adaptive internal substeps，求解土壤水、干空气、热；
   - 在接受后的土壤状态上求解土壤氮；
   - 把当前小时的诊断量写进 `diagnostics.Registry`；
   - 运行结束后统一写 NetCDF。

3. adaptive internal substeps
   - 每个 trial substep 前保存状态；
   - 选择 trial `Delt_t`；
   - 如果开启双孔隙，在 Picard 前做一次 exchange；
   - Picard 内迭代求解土壤水、干空气、热；
   - 收敛则接受 substep，并累加到当前小时输出；
   - 不收敛则缩小步长重试；
   - 如果最小步长仍失败，则整个外步 rollback。

交付物：

- 一张一页图或表：`forcing -> canopy -> plant C/N -> root sink -> soil substeps -> nitrogen -> diagnostics`

要问的问题：

- 每个 `KT` 开始和结束时，哪些字段是真正的 prognostic state？
- 哪些字段只在 substep trial 里临时有效？
- rollback 必须恢复哪些字段？
- 哪些字段只是重新计算出来的诊断量，不应该手动恢复？
- 哪些字段本来就设计成一小时或一步滞后，比如 LAI 或植物 C/N 输入？

## 2. Prognostic State 清单

重构前先建状态表。建议列如下：

| 状态块 | 所属 struct | 典型字段 | 由谁更新 | rollback 需要恢复？ | 输出 group |
|---|---|---|---|---|---|
| 土壤水 | `SoilVariables` | `hh`, `h`, `Theta_LL`, `Theta_UU`, `Theta_II`, `KfL_h` | `UpdateSoilWaterContent`, `solveSoilMoistureBalance` | 是 | `soil_water`, `solver_internal` |
| 土壤热 | `SoilVariables`, `EnergyVars` | `TT`, `T`, `QET`, `QEB`, 热参数 | `solveThermalEquations` | 是 | `soil_thermal`, `energy_balance` |
| 土壤气 | `SoilVariables`, `AirVars` | `P_gg`, 气相变量 | `solveDryAirEquations` | 是 | `soil_air` |
| 土壤氮 | `SoilVariables`, `NitrogenVariables`, `NitrogenFlux` | `C_NH4_new`, `C_NO3_new`, 反应源项, 累积通量 | `solveNitrogen` | 部分需要，且在土壤状态接受后更新 | `soil_nitrogen`, `soil_oxygen` |
| 冠层能量 | `fluxes`, `thermal`, `profiles`, `rad` | `Tch`, `Tcu`, `GPP`, `LE`, `H`, 辐射诊断 | `ebal_clm45` | 不参与土壤 substep rollback，但必须和当前小时状态一致 | `energy_balance`, `carbon_fluxes`, `canopy_profiles` |
| 植物水力 | `plantHydro`, 可选 `W_plant_state` | 根导度、叶水势、水分胁迫 | ebal 水力求解器和可选动态 wrapper | 取决于 scheme | `plant_hydraulics` |
| 植物 C/N | `NitrogenVariables` | C/N 池、NSC、分配、actual/potential 呼吸 | `runPlantCNStep` | 外步失败时需要 | `plant_nitrogen`, `carbon_fluxes` |
| 地表水 | `BoundaryCondition`, 局部历史数组 | `DSTOR`, `RS`, 入渗/径流项 | 土壤水边界逻辑 | 是 | `surface_water`, `conservation` |
| 纯诊断 | `reg`, 局部 history arrays | 输出用字段 | `reg.write` 和诊断 helper | 否 | 所有 NetCDF groups |

交付物：

- `docs/main_loop_state_inventory.md`，或者把完整状态表作为本文档附录。

## 3. 土壤水-热-气内部求解

这是内部计算里最先该查的一块。原因：土壤水热气是底座，冠层、根吸水、氮迁移、C/N 都依赖它。

关键文件：

- `src/UpdateSoilWaterContent.m`
- `src/+soil/+water/calculateMatricCoefficients.m`
- `src/+soil/+water/calculateMatricCoefficientsArray.m`
- `src/+soil/+water/solveSoilMoistureBalance.m`
- `src/+soil/+thermal/solveThermalEquations.m`
- `src/+soil/+gas/solveDryAirEquations.m`
- `src/+soil/+water/+dualporosity/dualPorExchange.m`

一个 trial substep 内的过程：

1. 保存 trial 开始时的状态，供失败 rollback。
2. 准备 substep 初值。
3. 如果开启双孔隙，在 Picard 前做一次 `dualPorExchange`。
4. 每次 Picard 迭代：
   - 可选 Anderson 初值修正；
   - 用 `UpdateSoilWaterContent` 刷新当前 Picard guess 的水力/冻融性质；
   - 组装水分方程矩阵系数；
   - 求解土壤水分方程；
   - 更新地表储水/径流 trial 项；
   - 如果开启，求解干空气方程；
   - 求解热方程；
   - 计算水分和能量误差；
   - 如果本轮有效，提交 Anderson map。
5. 水和能量都收敛，则接受 substep。
6. 把接受的 substep 通量累加回当前小时输出。

要问的问题：

- `OLD`、`SAVE` 和 trial 字段在拒绝 substep 后是否恢复一致？
- 哪些源项在 Picard 内保持常数是物理设计，比如根吸水、双孔隙 exchange？
- 水、气、热是否都用同一个当前 Picard 状态？
- 边界条件 switch 只是诊断，还是会改变求解路径？
- array kernel 和 scalar reference 是否仍然逐值等价？

验证：

- 3 天冬季/冻融窗口 A/B；
- 10 天困难窗口 A/B；
- 本地 A/B 干净后，再做一年或多站点验证；
- 对比水量残差、能量残差、rollback 次数、Picard 迭代分布和关键剖面。

## 4. 输入 Schema 审计

主循环状态流清楚后再做输入。目标是让每个输入字段都能追到一个模型决策。

关键文件和来源：

- `run_stemmus_initialize.m`
- `src/+core/ModelSettings.m`
- `src/+core/EnvironmentSettings.m`
- `src/+core/VegetationData.m`
- `data/config_*.txt`
- `data/` 下 forcing 和 site 文件

任务：

- 列出所有 config 文件，以及它们写入哪些 struct 字段；
- 标清单位和允许范围；
- 分开 site 常数、forcing 时间序列、PFT 参数、管理措施和诊断开关；
- 清理或记录重复 alias；
- 确保 `HydraulicKernel`、`MatricKernel` 这类运行时选择是正常 config，而不是改代码；
- 做一个长武站点模板，之后可以复制成其他站点。

交付物：

- `docs/input_schema.md`
- 字段名稳定后，再考虑机器可读 schema。

## 5. 输出 Schema 审计

目标：NetCDF 是权威输出；每个诊断字段的定义、单位和来源都明确。

关键文件：

- `src/+diagnostics/Registry.m`
- `src/+diagnostics/registerAllDiagnostics.m`
- `src/+diagnostics/+groups/*.m`
- `src/+diagnostics/flushToNetCDF.m`
- `docs/output_catalog.md`

任务：

- 每个注册变量都要有单位、source code、维度和物理含义；
- 以 `src/+diagnostics/registerAllDiagnostics.m` 和 `src/+diagnostics/+groups/register*.m` 作为当前事实来源，检查 `docs/output_catalog.md` 是否滞后；
- 单独标注非主循环输出，例如 `deep_learning` 的 `epoch` 维度和 `online_gs_validation` 的临时 write-only 诊断属性；
- 区分 prognostic state、瞬时通量、累积通量、诊断 helper；
- actual/potential 呼吸口径必须明确；
- conservation 组要小而权威；
- 决定哪些 legacy CSV 还需要，哪些只是可选输出。

验证：

- diagnostics smoke test；
- `output_catalog.md` 和 registry 一致；
- 没有注册但未写的变量，除非明确标成 optional；
- conservation 字段必须和状态更新使用同一套通量定义闭合。

## 6. 冠层能量、辐射和光合

这块在土壤 inner solve 和 I/O 稳定后再查。原因：它直接控制 GPP、ET、叶温、冠层温度和植物水力胁迫，科学敏感度高。

关键文件：

- `src/ebal_clm45.m`
- `src/ebal_clm45_profile.m`
- `src/heatfluxes_clm45.m`
- `src/+plant/biochemical.m`
- `src/+Radiation/` 下辐射传输函数

要问的问题：

- 当前冠层假设是什么：连续冠层、sunlit/shaded、垂直 profile，还是 profile mode？
- 哪些温度是在 ebal 里求解的，哪些只是诊断？
- `Tch`、`Tcu`、叶温、`GPP`、`LE`、`H`、辐射输出如何平均？
- 哪些字段传给植物水力和气孔导度？
- 哪些量可以缓存，哪些必须每次迭代重算？

验证：

- GPP、ET、`Tch/Tcu`、叶温、气孔导度、水分胁迫和能量残差；
- 同时对比 baseline 和观测。

## 7. 植物水力和根吸水

关键文件：

- `src/+plant/+PlantHydraulic/Plantrootconductivity.m`
- `src/+plant/+PlantHydraulic/solveEbalHydraulics.m`
- `src/+plant/+PlantHydraulic/solveEbalHydraulicsLayered.m`
- `src/COMPASS.m` 中根吸水和 sink 逻辑

要问的问题：

- 根吸水是在外步算一次，还是进入 substep 内部？
- 土壤 node/layer 值怎么映射到根吸水层？
- 根导度、root-soil contact 诊断是否和水分方程 sink 一致？
- 可选 dynamic plant water state 是否反馈主求解，还是只是后处理/诊断？

验证：

- 根吸水剖面积分是否等于预期的蒸腾需求或实际蒸腾；
- 叶水势和水分胁迫不能依赖隐藏的输出字段；
- 对比叶水势观测。

## 8. 植物 C/N、物候、周转和管理

关键文件：

- `src/+plant/+integrated/runPlantCNStep.m`
- `src/+plant/+integrated/applyNSCCoupling.m`
- `src/+plant/+integrated/turnoverAndHarvest.m`
- `src/+plant/` 下 allocation、demand、LUNA、phenology helper

要问的问题：

- 植物池子扣的是 actual respiration 还是 potential respiration？
- NPP、CUE、Reco、NEE 和池子变化是否使用同一套 actual/potential 口径？
- harvest、prune、root turnover、senescence、litter 是否在植物 C/N 预算里闭合？
- 哪些过程是 daily，哪些是 hourly，哪些由物候触发？

验证：

- plant C residual；
- plant N residual；
- NSC/storage 年内轨迹；
- 各器官池年变化；
- GPP/Ra/CUE 同时对比硬生理靶和软通量靶。

## 9. 多站点验证 Gate

只有主循环、输入 schema、输出 schema 稳定后，再做多站点验证。

每个站点最小验证包：

- config/site package；
- forcing 时间段和筛选规则；
- 观测靶和硬/软等级；
- 模型 NetCDF 输出路径；
- 指标定义；
- 已知过程 caveat。

第一批建议指标：

- 水量闭合；
- 能量残差；
- 土壤温度剖面；
- 土壤含水量剖面；
- ET 或 latent heat；
- GPP；
- 有数据时看 Reco/NEE；
- 有数据时看叶水势或冠层水分胁迫；
- 果园站点看植物 C/N 池是否合理。

## 10. 重构规则

- 梳理状态流时不要顺手改方程。
- 移动代码前，先列清楚它读什么状态、写什么状态。
- 替换数值 kernel 时保留 scalar/reference 路径。
- 诊断口径要同步更新 registry、NetCDF metadata 和文档。
- 每个结构性 cleanup 都必须有一个验证目标：状态相等、守恒闭合，或观测指标不退化。

## 近期交付物

1. `main_loop_state_inventory.md`：主要状态字段和 owner 表。
2. `input_schema.md`：config 和 forcing 字段地图。
3. `output_schema_audit.md`：registry 和 output catalog 一致性报告。
4. 土壤 inner-solve 流程图，标清 rollback 边界。
5. 第一个多站点验证 checklist，使用当前稳定 array-kernel 版本作为 baseline。
