# 长武无机氮初始化与观测回算：零步核查

Editor: Wenjie Liu; Date: 2026-09-29.

同日后续：用户已授权修改，并确认水分、DIN、SOC 初态均为实测。
当前实施与验收见 [实测初态映射修复](CHANGWU_MEASURED_SOIL_INITIALIZATION_FIX_20260929.md)。
下文保留修复前诊断和授权状态，不代表当前生产仍使用旧公式。

## 范围与证据

本轮只新增独立诊断脚本 `tools/diagnostics/verifyChangwuDinObservationMapping.m`；未改生产初始化、无机氮方程/边界、评分器、站点、观测或历史 campaign。模型积分步数为 0。最终证据为 `outputs/din_mapping_audit_20260929/20260929_201930/verification.json`，7 项代数/不变量检查 PASS。PASS 指核查及候选代数通过，不是当前生产映射已经修复。

较早的 `20260929_201842` 证据采用 Y.Theta_L 进行库存比较，仅用于追踪诊断水量与求解器水量的差异；其库存结论由上述最终证据取代，保留原文件不覆盖。

## 三个不同的映射

1. `src/+init/buildInitialState.m` 用初始节点 θ、节点容重 Y.bd 和 `equations.interpolateElementToNode` 的节点 Kd，将实测总 N [mg/kg 干土] 换算为液相浓度 [ug/cm3 水]。该 Kd 插值内部权重为对侧层厚（单元中心的逆距离权重），不是控制体储量的本侧层厚权重。
2. 初始水热辅助量刷新后，`init.initializeSoilAux` 经 `soil.nitrogen.materializeAcceptedState` 将 `derived.Theta_LL` 保存到 `Y.soil_nitrogen.theta_liquid_prev`。`advanceAcceptedCandidate.seedInitialStorage` 用这份快照及逐层 RHOKG 计算首次库存；不是 Y.theta，也不是评分器采用的 Y.Theta_L 投影。本例快照与 Y.Theta_L 最大差为 0.01820990734 m3/m3。
3. `scoreChangwuWaterGate.extractAccepted` 在初始诊断总量为零时，用 Y.Theta_L 的节点投影及追加末层值的层容重/RHOKG 回算；演化后使用诊断总量。后者的生产生成函数 `advanceAcceptedCandidate.updateTotalNitrogen` 取各层第一端点加最顶端点，仍不是共享节点两侧控制体的守恒投影。`evaluateChangwuBaseline` 和 `diagnostics.writeStep` 也各有保留口径。

因此，不能仅把评分器 t0 改成初始化公式来宣称问题解决；那会掩盖运行时初始库存不一致。

## 数字

| 指标 | NH4 | NO3 |
|---|---:|---:|
| 原初始化按自身公式回算最大差 [mg/kg] | 8.88e-16 | 7.11e-15 |
| 旧门评分 t0 回算最大差 [mg/kg] | 0.5558866282 | 0.7267249641 |
| 求解器容量守恒投影相对目标最大差 [mg/kg] | 0.3963261409 | 4.0697770181 |
| 当前求解器初始库存 [ug N/cm2] | 4078.537376 | 6494.296295 |
| 目标剖面对应离散库存 [ug N/cm2] | 4084.489562 | 6774.905846 |
| 当前库存相对目标差 | -0.1457266% | -4.1418960% |

目标来自现有实测剖面的原插值/外推规则及模型干土质量积分，不是新增的全柱 N 实测。该差异发生在模拟开始之前；积分后闭合检查通过不能反证初态准确。

## 候选修复与所需范围

在每个共享节点，把相邻层各半厚度的容量求和：

- 干土质量容量 D = sum(dz/2 × rho)，单位 g 干土/cm2。
- NH4 容量 B4 = sum(dz/2 × (theta_effective + rho × Kd))，单位 cm。
- NO3 容量 B3 = sum(dz/2 × theta_effective)，单位 cm。
- 初始液相浓度 C4 = measured_NH4 × D/B4，C3 = measured_NO3 × D/B3。
- 观测回算使用相同容量：N4 = C4 × B4/D，N3 = C3 × B3/D。

离散库存严格等于 sum(measured_N × D)，与 `calculateStorage` 的逐层梯形积分一致。本次候选副本回算误差约 1e-15，库存误差不超过 1.82e-12 ug/cm2；非均匀层厚、两侧含水量跳变和均匀极限的手算检查也通过。没有把候选状态保存为正式初态。

实际修复必须在水热初始容量完成后再进行实测 DIN 相态换算，并使用同一接受态观测算子；不能只修 scorer，也不能在第一步积分中无台账重置浓度。动态冰相、表面积水及未溶解肥料是否属于观测可提取总 N，需要按现有观测定义明确，不由本次无冰/无肥 t0 代数测试自动放行。

此前用户要求保留无机氮接口，本轮止于证据与候选验证。若允许修改上述初始化换算及诊断/评分口径，下一阶段做零步实测约束、守恒、极限/失败检查和短窗直接积分/恢复验收；不修改输运方程、硝化反硝化动力学或边界类型，不启动长期计算。

历史敏感性样本继续保留并标记版本，不宣称白跑。评分口径变更可先从已有 checkpoint 重评分；初始化实态变更后，旧轨迹不能仅靠重评分变成新初态结果。是否补算取决于基准及少量关键点的版本对照，不能现在断言需要全部重跑。
