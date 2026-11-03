# DOC/DON 底部出流：实施与验收

Editor: Wenjie Liu; Date: 2026-09-29.

同日后续状态：微生物颗粒出口现已另行接通，见
[BIOPARTICLE 实施与验收](BIOPARTICLE_BOTTOM_BOUNDARY_IMPLEMENTATION_20260929.md)。
下文是 DOM 单项实施时的记录，其“颗粒仍闭合”和 tracer schema=2 描述只适用于
当时快照。同日后续编号边界与当前 tracer schema=4 见
[有机物边界类型](ORGANIC_BOUNDARY_TYPES_20260929.md)。下文原始数值与证据保持历史归属。

## 1. 交付范围

长武 `data/sites/changwu/site.toml` 已显式启用：

```toml
[boundary]
dom_bottom_kind = "advective_outflow"
dom_bottom_water_phase = "liquid_only"
```

该段位于原有 `[boundary]` 内，不要重复建表。缺省或显式 `closed` 保留闭柱模式。
Kd 和弥散长度仍归 `soil.csv`；本次不改任何动力学参数、观测、初态或管理输入。
输出跨过的是当前 500 cm 计算域底面，不宣称模拟了更深土层或到达地下水的过程。

适用范围：物理底面 node 1、自由排水、仅向外流。顶面 DOM 外部通量和底面扩散/
弥散通量为零。水的主方程未改；NH4/NO3 边界未改；BIOPARTICLE 共用输运核的
旧调用仍为闭边界，POC/MOC、吸附或附着颗粒没有被直接排走。

## 2. 显式水相合同（不是任意 QMB=液流）

`calculateCoupledBoundaryFlux` 给出接受的**总水**反力 QMB，不能一般地当成液态水。
本次在站点声明新的适用假设：自由排水底部无水汽交换（liquid_only）。在此条件下，
水 owner 的 `soil.water.bottomLiquidBoundaryForSolutes` 定义接受的液流为 QMB、
底部水汽流为 0。该分相是明确的边界闭合假设，不是从单元中心 QL 推算的水汽测量。

记录包含 physical_column_bottom、into_domain 正方向、cm/s、起止秒数、dt、
phase_contract。只有接受水子步进入 chronological schedule；DOM 按相同 dt/θ/QL
序列消费每个外部底面液流。正 QMB（向上入流）、非自由排水或升高的地下水接口
明确报错；没有暗设外部 DOC/DON=0，也没有实现地下水双向交换。

底部质量通量为 `(-QMB) * C_aqueous`，不是 `(-QMB) * total_mg/kg`。
线性吸附：`C_aqueous=S_total/(theta_liquid+rho_bulk*Kd)`；每个接受输运子步重分配。
动态吸附：只运移 total−sorbed，显式 sorbed 库保持不变。出流与内部面交换共同
进入非负性缩步；不事后扣库、不以 clamp 消除负值。

## 3. 端到端所有权与收支

| 项目 | 活动实现/合同 |
|---|---|
| 配置 | `soil.dom.normalizeBoundary`；`io.checkSite` 与 `init.buildSoilDOM` 均消费；未知 dom_* 名字报错 |
| 接受水驱动 | `solver.advanceSoilEquationOuterStep` → `recordAcceptedTransportTile` → `soil.dom.acceptedBottomDrivers` |
| 输运核 | `soil.dom.transportConservative` 默认闭；可选具名底部液流；DOC/DON 分别按自身浓度出流 |
| 总质量账 | `M_after−M_before+bottom_outflow=0`；新增按浮点尺度/子步数的闭合断言 |
| 累计记忆 | `Y.soil_dom.bottom_export_cumulative.{doc,don}_g_cm2`，仅开放模式拥有；候选上累计，原子 commit 才接受 |
| DON 标签 | 冻结总交换图包含 bottom_outflow，隐式混合给标签份额；step/cumulative.don_bottom_loss_g_cm2 是终端损失 |
| 标签闭合 | 标签留存量+标签底部出口闭合，0≤tag≤total；完整化肥累计账加入 DON 出口，tracer schema=2 |
| 诊断 | DOC/DON step、flux、cumulative；化肥 DON step/cumulative 和底面有向交换；这些新增出口账用 double NetCDF |
| 重启 | `checkpointCompatibilitySignature.dom_boundary` 绑定模式/相假设/账本版本；不同模式或旧无签名检查点拒绝 |

质量单位 g C 或 N/cm²；flux 为 g/cm²/s；乘 1e5 换为 kg/ha。出口不是呼吸、矿化、
数值损失或内部有机/无机氮交换。外部输入本版恒为零，不提供一个未消费的输入浓度键。
来源标签仍采用既有“完全混合、无源偏好”归因离散；未引入标签对总量物理的反馈。

## 4. 验证结果与限制

证据目录：`outputs/dom_bottom_boundary_20260929/`。
验证入口：`tools/diagnostics/verifyDOMBottomBoundary.m`。

- 单层手算、零流、零浓度、Kd=0/高 Kd、反向流拒绝、干底层拒绝、配置错误拒绝。
- 小库存自适应缩步和冻结子步重放；多层 C/N 质量与标签界限；动态吸附库保持。
- 两个非等长水子步与逐段调用一致；反复构造候选不改变输入、不重复累计。
- 单层解析解时间细化误差依次 1.212e-8、6.048e-9、3.021e-9 g/cm³，约一阶收敛。
- 48 h 新代码 closed 与冻结旧代码**接受态逐位一致**，只排除 tracer schema 升级及
  新增的零累计标签字段；开放/闭合配对的 h、T、Theta_L、theta 逐位相同。
- 开放 48 h 累计 DOC 出口 `1.1355265590690563e-6 g C/cm²`（0.113553 kg C/ha），
  DON `2.4995112464177282e-7 g N/cm²`（0.0249951 kg N/ha）。这是短窗模型预测，
  不是实测通量，不能线性外推全年。
- 独立用逐层输运改变量+出口重建质量账，DOC/DON 最大残差分别
  `1.216e-18 / 8.150e-20 g/cm²`。
- 再独立求和不重叠有机状态库，加呼吸、DIN 交换、管理/凋落输入及新出口，完整
  有机 C/N 库存账逐小时最大残差约 `1.295e-15 / 1.194e-16 g/cm²`。
- 保留原生产 C/N/无机 N/氧气断言；化肥源累计闭合最大残差 `8.314e-19 g N/cm²`。
- 专门用“初始 DON 的 30% 为标签”的**合成测试夹具**验证非零出口：标签损失
  `1.7889101126744534e-9 g N/cm²`，柱闭合残差 `2.972e-20`；未写入任何正式初态。
- KT24 检查点续跑 KT25 与连续运行逐位相同；闭/open 不匹配和旧无边界签名拒绝。
- NetCDF 出口、累计出口与内存接受账逐位一致。
- 合入后的正式站点不加测试覆盖直接初始化并积分 1 h，通过；既有接受子步
  DOM/BIOPARTICLE 重放、80 个随机有界来源输运及有向面诊断回归均通过。

`production_verification.json`、`label_probe_verification.json` 保存实际数字和状态；
`production_evidence.mat` 保留初末态、逐小时出口及独立库存残差。
模型级来源规范影响了实现：分相合同放在水 owner，出口进入候选/commit 与重启账本，
而非仅在输出端补一列通量。

验收限于单元/构造测试和 48 h 耦合窗口；不冒称全年或六年验证。未实现底部入流、
地下水浓度驱动、顶部灌溉 DOM、地表溶解态径流或微生物颗粒底部出口。

## 5. 五组机制试验和下一科学步骤

边界开发先在 staging 副本进行；五组旧边界 1464 h 对照全部 PASS、独立审计 PASS、
旧身份复核 PASS 后才更新主源码。旧 `frozen/`、全部结果、原始签名均保留。
主源码现已是新身份，旧 runner 的根目录哈希检查会拒绝追加，这是预期的保护。

五组结果见 `CHANGWU_CN_CONTRASTS_20260929.md` §6：初始 DON 分配控制了较多深层
NH4 增量；摄取相关速率组控制 DOC 耗尽速度；二者不能用同一调整替代。
新边界是补齐物理出口，不承诺提高观测拟合；DOC 原已偏低，增加出口可能更低。

后续应在新身份下将 open/closed 基线延伸到首个观测日，量化底部负荷与剖面差异，
再决定哪些已有敏感样本需要补算；没有重启 44 组，也没有启动年度 BO/六年碳率定。
正式 DOC/DON 率定还需要确认提取测定方法/土液比和 2000 初态来源。不能把 2006
DON 同时用作独立的 2000 初态和验证真值，也不能挑 total/free 中较好拟合的口径。
