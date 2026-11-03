# 微生物颗粒底部出流：实施与验收

Editor: Wenjie Liu; Date: 2026-09-29.

同日后续更新：当前站点已改为有机物编号边界 top=2 / bottom=3，支持的其他类型、
限制和最新验收见 [有机物边界类型](ORGANIC_BOUNDARY_TYPES_20260929.md)。下文保留底部出口
首次实现时的配置、tracer schema=3 和原始证据，不代表当前配置仍使用旧 kind 键。

## 1. 现在会排出什么

当前长武站点已在原有 `[boundary]` 下启用：

```toml
bioparticle_bottom_kind = "advective_outflow"
bioparticle_bottom_water_phase = "liquid_only"
```

本项独立于 DOC/DON 出口，且只有 `switches.microbial_transport=true` 时发生。
缺省 `closed` 保留旧闭边界；关闭运移时没有颗粒出口。边界是站点全局设置，
附着/脱附/截留和弥散等土层参数仍在 `soil.csv`，本次没有改率定参数或初始库存。

| 状态 | 是否能直接跨越 500 cm 底面 |
|---|---|
| MB_mobile：移动态活跃微生物 | 是，随底部液态排水 |
| BD_mobile：移动态休眠微生物 | 是，随底部液态排水 |
| ENZN_mobile：移动态胞内酶货载 | 是，与细胞使用相同输运算子和时间步 |
| attached：附着态 | 否，先脱附进入 mobile，才能在后续输运中排出 |
| retained：截留态 | 否；本次不改变既有不可逆截留设定 |
| POC、MOC、胞外土壤酶、吸附 DOM | 不通过这个微生物颗粒出口直接排出 |

这表示物质离开当前计算域，不表示已经模拟了 5 m 以下的运移或到达地下水。
DOC/DON 溶解态出口与微生物颗粒出口是两条独立路径，不能把同一微生物 C/N
再算一遍 DOC/DON 淋失。MB、BD、ENZN 本身也是三个不重叠的状态库。

## 2. 物理、数值及所有权合同

沿用 [DOM 底部出流合同](DOM_BOTTOM_BOUNDARY_IMPLEMENTATION_20260929.md) 的
自由排水、零底部水汽交换 `liquid_only` 假设。水 owner 记录每个接受水子步的
物理底面液流、起止时间、位置、单位及正方向。不是用最后一个单元中心 QL
替代底面流，也不把任意总 QMB 无条件等同于液态通量。

以进入计算域为正，接受底面 `q_liquid_cm_s <= 0`。对于 C 或 N：

```text
底部向外通量 = (-q_liquid) × (mobile_bulk / theta_liquid)
单位：cm/s × g/cm³水 = g/cm²/s
完整颗粒库存变化 + 颗粒出口 = 0
```

底部外部扩散/弥散通量为零，顶部外部颗粒通量为零。反向底部入流会报错，
因为没有外部微生物浓度；不暗设入流浓度为零。适用范围为物理底面 node 1
和自由排水，不支持升高地下水接口或双向地下水交换。

碳状态仍是 top-down mg C/g 干土，N 状态仍是 bottom-up g N/cm³土。
仅在既有适配层转换。开放模式先由 θ、厚度、内部液流和底面液流计算一个
与浓度无关的非负性时间步上限，使显式算子的最大单元流出系数乘 dt 不超过 0.9。
MB/BD/ENZN 的 C/N 因此共用同一接受时钟，不会因不同浓度分别缩步而拆散配对。
0.9 是数值安全系数，不是新增生态参数。闭边界沿用旧时钟以保证历史回归。
每个接受水段仍先输运，再按既有方法进行附着/脱附/截留；未重排反应流程。

活动路径：

- 初始化：`init.buildSoilMicrobes` → `soil.microbe.normalizeTransportBoundary`。
- 水驱动：`solver.advanceSoilEquationOuterStep` → 接受子步记录；DOM 或颗粒任一
  开放时记录底面，二者分别通过 `soil.dom.acceptedBottomDrivers` 校验只读驱动。
- C/N 候选：`soil.nitrogen.advanceAcceptedCandidate` →
  `soil.microbe.bioparticleTransport` → 共享保守输运核。
- 累计量：`Y.soil_microbes.bottom_export_cumulative.{MB,BD,ENZN}.{carbon,nitrogen}_g_cm2`，
  只在候选内更新，完整步接受后提交；重复构造候选不会重复累计。
- 肥料来源标签：冻结总 N 交换账包括底部出口，标签按既有有界混合方法分配；
  新增三个 `bioparticle_{mb,bd,enzn}_bottom_loss_g_cm2` 终端损失，进入柱闭合。
  标签 state schema 从 2 升为 3，不反馈到总量物理。
- 重启：`checkpointCompatibilitySignature.bioparticle_boundary` 绑定模式、相假设、
  数值时钟和出口账版本；闭/open 不同或缺少新合同的旧检查点拒绝直接续跑。

生产守恒容差未放宽，未靠事后 clamp 去除负库存。项目编码规范促使出口同时进入
物理扣库、候选/提交、来源归因、累计记忆、重启和诊断，而不是只增加输出通量。

## 3. 输出如何看

`nc/soil_organic.nc` 新增 24 个 float64 变量：

```text
bioparticle_{mb,bd,enzn,total}_{carbon,nitrogen}_bottom_outflow_{step,flux,cumulative}
```

step、cumulative 为 g C 或 N/cm²，flux 为 g C 或 N/cm²/s。
`total` 是三个组分的求和，不能再把 total 与三个组分相加。g/cm² 乘 `1e5`
变为 kg/ha。关闭出口或关闭运移时，出口为物理零。

`nc/soil_nitrogen.nc` 新增六个 float64 来源损失变量：

```text
fertilizer_tracer_bioparticle_{mb,bd,enzn}_bottom_loss_{step,cumulative}
```

原来源标签的有向底面交换也包含对应出口，顶面仍为零。
这些是“带化肥来源标签的 N”，并非全部颗粒 N；不能与总 N 出口重复相加。
做新边界下的化肥去向收支时，应同时纳入 DON 和这三个颗粒出口，不只统计无机氮淋失。

## 4. 已完成的验证

入口：`tools/diagnostics/verifyBioparticleBottomBoundary.m`。
证据：`outputs/bioparticle_bottom_boundary_20260929/`；`frozen/` 为修改前源/站点快照。

- 单层手算 C/N 出流；附着/截留库无直接扣除；零液流、零 mobile 库极限。
- 未知配置、缺少底面驱动、反向入流拒绝；重复候选不改输入、输出逐位一致。
- 非均匀 C:N、多层强排水、扩散/弥散及相间交换：各组共同时钟、非负、C/N 闭合，
  标签 0/30%/100% 及随机份额保持有界、留存量加出口闭合。
- 非等长/变 θ 液流序列与顺序调用一致到单位适配浮点误差，实测最大 1 ULP；
  不能把多一次 mg/g ↔ g/cm³ 转换误称为逐位相同。未改变生产容差。
- 单层解析排空解细化误差：`1.1971e-4 → 5.9464e-5 → 2.9636e-5 mg/g`，约一阶收敛。
- DOM 关闭/颗粒开放，以及颗粒运移关闭的独立开关测试通过。
- 48 h 颗粒 closed 与冻结修改前完整接受态逐位一致，仅排除新 tracer schema
  元数据和三个已断言为零的新来源累计字段；原有 DOM 开边界保持。
- 48 h 颗粒 open 耦合试验：独立求和不重叠的完整有机 C/N 库，加入 CO2、DIN
  交换、管理/凋落输入及 DOM/颗粒出口，逐小时最大残差为
  `8.4331e-16 g C/cm²`、`1.2403e-16 g N/cm²`。
- 开/闭配对的 h、T、Theta_L、theta 逐位相同；化肥源柱累计最大残差
  `7.6656e-19 g N/cm²`。
- KT24 检查点续跑 KT25 与连续运行逐位一致；不兼容/缺签名重启被拒绝。
- 24 个总量出口 step/flux/cumulative NetCDF 均与接受账逐位一致。
- 人工标记现存微生物 N 的 30%（仅测试标签，不改总库）后，MB、BD、ENZN 均有
  非零标签出流；柱闭合残差 `5.5013e-20 g N/cm²`，六个标签出口输出逐位一致。
- 原接受子步重放、80 个随机有界来源输运、有向面诊断、DOM 边界单元回归通过。

48 h 颗粒出口（模型短窗预测，非观测）：

| 组分 | C，g/cm² | N，g/cm² |
|---|---:|---:|
| MB | 1.5570581598e-8 | 1.9314675206e-9 |
| BD | 9.3260172231e-7 | 1.1657583342e-7 |
| ENZN | 8.8154485070e-15 | 2.9384828357e-15 |

合计约 **0.094817 kg C/ha、0.011851 kg N/ha**。当前初态/短窗下休眠态占主要
部分，不能据此宣称全年通量或完成率定，也不能直接线性外推。

五个 `*_verification.json` 分别保存 unit、production、label_probe、switches、
output_audit 结果；`production_evidence.mat` 保存初末态、逐小时出口和独立库存残差。

## 5. 与率定的关系及边界

本次完成的是开放土柱物理出口及短窗工程验收，不是率定最优参数。
未重跑 44 组、未启动年度优化或六年基线，旧 campaign、数据和签名保留。
新源码/边界必须使用新身份，不能接着旧 campaign 混跑。

下一科学步骤仍建议先做“新边界基线—旧边界对照”至首个观测日，查看剖面偏差
与出口负荷的实际影响，再决定需要补算哪些敏感性样本；不是先假定全部重跑。
DOC/DON 提取方法、土液比、2000 初态来源仍未解决，不能用新增出口替代这些确认。
