# 长武有机库周期交接与实测约束初态：实施记录

Editor: Wenjie Liu; Date: 2026-09-29.

同日后续：经用户授权，DIN 映射差异连同实测水分的生产端点回算已修复，
见 [实测水分/DIN/SOC 修复与验收](CHANGWU_MEASURED_SOIL_INITIALIZATION_FIX_20260929.md)。
本文下面保留先前接口实施状态与发现记录，不能将其旧映射限制理解为修复后的现状。

本记录补充 [初始化方案](CHANGWU_CONSTRAINED_INITIALIZATION_PLAN_20260929.md)，不覆盖历史试验或初值来源。仅完成接口与短窗验收；未启动年度预热、长期预热或正式率定，未改变无机氮方程、边界、站点输入或冻结 campaign。

## 1. 两个接口

`init.captureBiogeochemistryCycleCarry` / `init.applyBiogeochemistryCycleCarry` 用于代理驱动预热周期交接。它们交接完整 accepted Y、根系水状态和已接受气温历史尾段，重建派生缓存；新周期重置运行时钟并重新启用该周期管理事件。首步使用上一周期已接受地表温度，生长温度窗口使用连续气温历史。普通冷启动不走这些特殊分支。

身份检查覆盖参数、网格、边界、驱动、原始初态输入、实际加载的静态模型库数据和 src MATLAB 源码；输出路径和运行时间不参与身份。checkpoint 保存周期上下文并在恢复时复核。当前驱动仅支持显式恢复第一周期，不宣称任意后续周期都可通过同一入口恢复。

`init.materializeMeasuredBiogeochemistryInitialState` 用于正式模拟的零步初态构建，目标必须是尚未开始运行的新模型。保持该参数组重新初始化所得 NH4/NO3、水热和植被状态；逐层以实测 SOC 对应的初始化剖面约束有机 C 总量。该剖面的插值/外推层并不因此成为独立逐层实测。

具体做法是固定新初态 IOM，按一个逐层非负比例缩放供体可响应 C 分库及对应有机 N，保留相态比例与 C:N 配对。DOC 总量只计算一次；AOM 不进入当前 SOC 算子，单独保留并计入土壤加地表 AOM 的初始化库存。POC/EP 等合计字段重建，不重复相加。没有逐层 TN 实测约束，因此这些 N 分库仍属条件先验。

预热累计呼吸、出流及施肥来源标签不导入正式时期。初态条件化造成的 C/N 库存差额独立记录，不作为积分中的物理通量。`periodic_state` 必须是完整配置窗口末态且与正式起点同日历相位；`initial_guess` 和 `diagnostic_state` 仅供代数/技术测试，不证明预热收敛。初态接口尚未形成冻结 campaign 的完整来源清单，审计字段明确标记未完成来源/驱动哈希准入和科学放行。

## 2. 验收证据

- 实测约束初态：29 组检查 PASS，共 48 个已接受测试小时，最长连续窗口 24 h。包括逐层非均匀分库、单位/方向、相态恒等式、保护字段逐位不变、确定性、幂等性、非法供体拒绝，以及直接积分和断点恢复的完整 Y/管理状态逐位一致。证据：`outputs/measured_initialization_20260929/20260929_195327/verification.json`。
- 周期交接：PASS，共 10 个已接受测试小时。普通运行与周期重放 checkpoint 恢复均逐位一致，历史两小时状态回归一致，7 类非法合同被拒绝，恢复环境变量隔离并还原。以两小时技术周期检验缓存、根系水状态、地表温度及气温历史交接、管理重启和身份检查；这不是年度收敛试验。证据：`outputs/bgc_cycle_carry_check_20260929_195259/verification.json`，身份 schema 2。
- 测试继续使用原生产无机氮守恒容差 `1e-6 ug N/cm2`；未放宽容差。该单位不同于早期筛查报告所述 `mg N/cm2`，不能混用。

合成供体用于检验非均匀分库代数，另有 12 h 实际积分末态测试，但二者都不作为已预热的科学初态。

本次结束时，`site.toml`、`soil.csv`、`soil_initial_profile.csv`、`organic_initial_pools.csv` 与已完成桥接 `campaign.json` 的 SHA256 均与本轮开始记录一致。

## 3. 尚未解决的问题与下一步

同日零步后续核查：差异还包含初始化 θ 与氮求解器有效液态快照的不同。实际库存口径下
NO3 的初始库存比目标离散库存低约 4.142%；候选守恒换算已在独立副本通过代数检查，
未部署。详见 [DIN 映射核查及候选范围](CHANGWU_DIN_MAPPING_AUDIT_20260929.md)。

发现原无机氮初始化和水氮评分器 t0 的深度映射不同：初始化使用节点容重、厚度加权节点 Kd 和初始含水量；旧评分回算使用另一套层到节点映射/accepted 含水量。基准零步检查中，NH4/NO3 最大差分别为 0.5558866282 / 0.7267249641 mg N/kg。初始化按自身算子回算误差约为浮点舍入量级，不能据此声称评分算子也一致。

本轮按授权保留无机氮代码，仅将差异写入审计。正式率定前应先统一观测回算口径并重验，不应让优化器补偿这个映射差异。

随后仍需年度日历/作物管理相位验收、独立全柱 C/N 收支核查、带单位的分库收敛判据以及不同起始猜测的收敛分支比较。短周期交接 PASS 不等于年度周期稳定，更不等于未观测分库被唯一确定。正式率定前还需冻结完整初态来源与驱动身份。

## 4. 2026-09-30 预热准入核查

Editor: Wenjie Liu; Date: 2026-09-30.

本节取代前文同项工作的待办状态，不覆盖历史证据。按 `compass-coding-standard`
新增只读诊断、明确单位的比较器和初态文件接口，没有更改反应、输运、边界、植被方程、
站点输入或冻结年度 campaign。没有新增年度样本，尚未启动自然三周期预热或正式优化。

### 4.1 年度独立 C/N 预算：工程检查通过，物理碳准入未过

新增 `diagnostics.biogeochemistryInventory` / `diagnostics.biogeochemistryBudget`，
控制体为土柱、地表 AOM 和植物实体库。独立计算库存与实际外部累积输入/输出之差，
不用各组件已有 residual 倒推缺失通量。覆盖液态/吸附/冰相无机氮、尿素与颗粒肥料、
独立有机库、所有已启用的 DOM/颗粒边界、GPP、呼吸、收获、施肥及氮气体损失。
凋落、吸氮、矿化等内部交换不作为外源重复相加；来源标签不是新库存。

对冻结的新实测初态基准 41，重建同版本零步初态、读取 KT8760 末态，不再积分：

| 检查项 | 实际结果 |
|---|---:|
| 实体 C 库存起点 → 终点 | 30947.8503 → 30758.9976 g C/m² |
| 实体库存变化 − 物理外部净输入 | +29.6058623324 g C/m² |
| 显式纳入债务核销后的记账残差 | −1.5370×10⁻¹⁰ g C/m² |
| 总 N 库存变化 − 外部净输入 | +8.3397×10⁻¹² g N/m² |
| 年初 / 年末负 `xsmr` 债务 | 0 / 0 g C/m² |
| 全年记录的债务核销 | 29.6058623325 g C/m² |

代码事实：`vegetation.integrated.applyNSCCoupling` 的欠账分支允许维持呼吸先形成
负 `C_xsmr`；`turnoverCarbonMirror>winterWheatTurnoverC` 在收获清除未偿还欠账，
`runPlantCNStep` 将其计入 `cum_C_xsmr_writeoff`，且明确不把它当作负呼吸。
这项记账调整解释了上述物理差额，但不能当作真实光合或其他碳输入。
实体总 C 实际在下降，不能把差额直接表述成“凭空新增了 29.6 g C/m² 的实体库存”。
参考模型注释尚未独立核实，本节不据此宣称当前策略已获文献验证。

验证器 `tools/diagnostics/verifySpinupAnnualBudget.m` 的 **11 项检查 PASS**，包括
未记账 C/N 注入识别、缺失账本拒绝、DIN 库存与生产控制体对照、别名与内部交换检查。
冻结源码/输入/工具/设计匹配 campaign，历史文件前后 SHA256 相同。
证据：`outputs/spinup_annual_budget_20260930/20260930_095359/verification.json`。
N 使用原有数值门 `1e-6 ug N/cm² = 1e-12 g N/cm²`；C 的 `1e-10 g C/cm²`
仅为报告明确声明的数值诊断阈值，不是新批准的物理调整容限。
报告分别标记 `numerical_budget_passed=false`、`scientific_release=false` 和
`annual_launch_allowed=false`。不得通过删除核销记录、重命名为物理输入或放宽容差放行。

### 4.2 年度终点技术交接及管理重置

`tools/diagnostics/verifyAnnualBiogeochemistrySeam.m` 读取已有基准年度末态，
验证完整 accepted Y、根水历史、地表温度与气温尾段交接，时钟重置到 KT0；
下一周期前两小时积分成功，第二小时重启完整 Y 与管理状态逐位一致。
另用生产事件匹配器遍历两次完整年度日历，每次矿质施肥 1 次、有机管理事件 1 次，
窗内事件各触发一次，且与原年度已接受事件标志一致。

这只新增 **3 个诊断积分步**（第二小时重启计一次），不是第二年全年积分。
匹配器日历遍历也不是物理模拟。历史 checkpoint 没有预热来源身份，未将其追认成
首个预热周期；科学初态导出仍拒绝缺少积分前身份的供体。
证据：`outputs/annual_bgc_seam_20260930/verification.json`（PASS）。

### 4.3 分库收敛比较器

`soil.organic.spinupStateVector` 将独立 C/N 库按层、节点、端点或相态统一到
`g element/cm²`；聚合别名只核验、不重复计数，累积呼吸/出流不当作库存。
`compareSpinupCycles` 逐位置检查
`abs(end-start) <= abs_element_g_cm2 + rtol*max(abs(start),abs(end))`。
必须显式提供 C/N 绝对阈值、相对阈值与来源，没有混用原生单位的默认绝对下限。
通过此库检查不自动声明完整周期平衡。

29 组合成测试 PASS；真实当前初态和基准 KT8760 均解析为 **3233 行独立状态**，
有限值、单位、相态及别名检查通过。证据：
`outputs/spinup_convergence_check_20260930_095104/verification.json` 及同目录
`production_schema.json`，新增积分为零。
尚未制定/批准实际科学绝对阈值，也未替换正式预热 driver 的旧式比较分支。
不能将合成测试使用的浮点舍入阈值直接用于科学预热放行。

### 4.4 实测约束初态文件

新增 `init.exportMeasuredBiogeochemistryInitialArtifact` 和
`init.applyMeasuredBiogeochemistryInitialArtifact`。文件绑定供体 checkpoint 的 SHA256、
积分前源码/参数/输入/驱动/实际模型库身份、实测目标、条件化差额及验收证据。
参数发生变化不能近似复用；来源与收敛证据必须同时指向同一供体及身份。
缺少证据时明确 NOT_READY，默认拒绝用于正式模拟；显式诊断使用不改变未放行状态。
保留实测水分、NH4、NO3、SOC 约束，正式时期累积账本重新从零开始。
初态差额只记为初始化调整，不计入真实通量。

独立供体 1 h + 文件应用后 1 h 的验证由
`tools/diagnostics/verifyMeasuredInitialArtifact.m` 执行，**19 项检查 PASS，合计 2 个接受小时**。
证据为 `outputs/measured_initial_artifact/20260930_095136/verification.json`；
源码/输入身份首尾一致，未放行初态的默认正式应用确实被拒绝。
不伪造已通过的年度收敛证据来生成科学初态；本轮只产出诊断初态文件。

### 4.5 后续顺序与停止条件

1. 先单独核对冬小麦维持呼吸供碳、欠账偿还和收获清算策略，再决定范围明确的过程修正。
   当前未获本轮新增植物方程改变的具体科学选择；不直接切换呼吸方案或删除债务来凑账。
2. 修正若获授权，先做局部供碳/收获守恒与短窗回归，再补一条受影响的年度基准验证。
   现有五点结果保留为历史对照；此发现不自动要求全量重跑 44 点。
3. 全系统物理收支门通过后，冻结科学收敛政策并接入正式预热 campaign、初态证据出口。
   首批最多三个自然年度周期只是计算预算，不保证三年收敛；不得用任意调整加速凑出稳定。
4. 之后检验不同未观测初值分支、相同日历相位的环境与通量稳定性，形成实测约束正式 t0，
   才接水氮校准/独立验证、低维 SOC 率定及更长时期耦合验证。

## 5. 2026-09-30 冬小麦实体供碳修正

Editor: Wenjie Liu; Date: 2026-09-30.

用户授权修正后，本节取代 §4.5 中“尚未获过程修正授权”的状态；历史年度和零步审计
保留原状，不以新源码重新解释旧结果。本次仅修正冬小麦供碳及相应氮需求接口，
不修改实测初值、无机氮反应/输运方程、有机物边界或既定数值容差。

### 5.1 实际过程与适用范围

仅对 `winter_wheat_cn == 1` 启用实体供碳。保留当前 `resp_scheme = 0` 及潜在
维持呼吸的温度/组织 N 公式，不切换到苹果的经验底物响应、不增加率定参数。
单步可动用碳量为当步非负 GPP 输入加已有 labile、root storage、stem storage。
在同一 `g C/m²/step` 单位下：

```
available = GPP_step + labile_previous + root_storage_previous + stem_storage_previous
Rm_actual = min(Rm_potential, available)
Rg_actual = min(Rg_potential, max(available - Rm_actual, 0))
growth_basis = max(max(GPP_step - Rm_actual, 0) - Rg_actual, 0)
```

维持呼吸优先，生长呼吸共享剩余余额；碳分配端实际扣减对应库。已有可动用储备
可支付夜间呼吸，但本次不引入结构碳消耗、死亡过程或储备驱动的夜间结构生长。
微量 root/stem 储备转移按各自真实余额分配，避免旧分母下限在空 stem 中产生负数。
浮点资金检查仅允许有限次加减的舍入误差，不放宽物理守恒门。

`runPlantNDemand` 将同一冬小麦标志传入 `resp_only` 计算；氮需求和实际 C 分配
采用相同供碳基础。无机氮公式未改，但植物实际吸氮可以随有碳支持的生长而变化，
这是耦合响应，不应宣称整条 N 轨迹逐位不变。

冬小麦不再产生 `xsmr` 欠账；收获不再清除或核销非零欠账。旧非零债务显式报错，
而非静默置零。已有正 labile 收获清算仍作为真实 C 输出保留，收获器官、根系凋落
和保种划分不变。非冬小麦分支保持既有计算。

### 5.2 与参考模型的关系

本次不是声称恢复 CLM/CoLM 原算法：其原生实现确有带符号 `xsmr` 欠账和收获
清算。已核对固定提交的
[CoLM 供碳分配](https://github.com/CoLM-SYSU/CoLM/blob/626347a9156ab97db4ca846f47cd1b7f2728f705/main/BGC/MOD_BGC_Veg_NutrientCompetition.F90#L307-L348)
及 [CoLM 收获清算](https://github.com/CoLM-SYSU/CoLM/blob/626347a9156ab97db4ca846f47cd1b7f2728f705/main/BGC/MOD_BGC_CNCStateUpdate1.F90#L475-L481)。
COMPASS 此处选择独立实体库存必须由实际输入/输出闭合的合同，因此明确改变冬小麦
供碳规则，不能把旧欠账核销重新命名为物理碳输入。
[WOFOST72](https://github.com/ajwdewit/pcse/blob/67a28e56b0e34655f8d60b0b4a254a7c81efbb2f/pcse/crop/wofost72.py#L209-L214)
也区分潜在与可支付维持呼吸，但其日尺度同化限制不等同于本模型小时尺度的 NSC
储备实现。上述参照支持设计比较，不是长武生理表现已获验证的证据。

### 5.3 版本隔离与当前证据

新 checkpoint 对启用冬麦 C/N 的运行增加
`winter_wheat_carbon_funding = physical_nsc_v1`。缺少/不同标志的旧 schema-2
和旧 schema-1 checkpoint 均拒绝进入新过程；即使旧末态 `xsmr == 0` 也不例外。
历史年度资料继续作为旧过程证据，不能拼接成新版本轨迹。

`tools/diagnostics/verifyWinterWheatCarbonFunding.m` 已通过 **106 项纯函数/合同检查**，
新增模型积分为零。包括冻结旧源码的真实失败复现、昼夜/零碳/储备/高 N/微量库存
物理 C/N 守恒、收获路径、非冬麦各呼吸方案及实际 Apple/NET 的 exact parity、
空可选标志兼容性，以及 checkpoint 标志/旧版本拒绝。
证据：`outputs/winter_wheat_carbon_funding_final_20260930/verification.json`。

耦合 **48 h PASS**，MATLAB 正常退出：物理 C 残差 `-5.39499001028787e-16 g/cm²`，
N 残差 `-1.25333771139324e-16 g/cm²`；每步欠账及累计核销均严格为零。
48 步 DIN 最大残差 NH4/NO3 分别为 `3.865352482534945e-12` /
`5.10347319959692e-12 ug N/cm²`。实测水分/DIN/SOC 的初始化回算和源码/输入/工具
身份首尾检查通过。证据：`outputs/winter_wheat_carbon_funding_short_20260930/result.json`。
当前版本重启/周期交接回归 **PASS**，共 10 个已接受测试小时。普通运行和两小时技术
周期的 checkpoint 恢复均完整 Y 逐位一致，7 类非法合同拒绝；短周期不报告年度收敛。
证据：`outputs/winter_wheat_carbon_funding_restart_20260930/verification.json`。
旧两小时 reference 因初始化协议不一致未作新旧 Y 对比，不将此项写成历史逐位通过。
上述短窗检查不替代全年物理预算。
专用驱动 `tools/diagnostics/runWinterWheatCarbonFundingBaseline.m` 限定只验证基准
sample 41，保留七参数逐项等于 baseline 的断言；年度必须从新 `frozen` 副本运行，
使用 fresh 实测约束初态。每步检查零欠账/零核销和原 DIN 门，每 10 KT 独立检查
全系统物理 C/N 预算，保存初态、观测时刻和终点 checkpoint，并首尾复验身份。
年度验证不等于预热收敛，不自动放行 44 点重算、正式率定或六年碳基线。

2026-09-30 10:16（Asia/Shanghai）已启动一条新年度基准：
`outputs/changwu_crop_carbon_funding_20260930/annual`，固定 sample 41、8760 步，
源码/输入/设计独立冻结于同目录 `frozen`。899 个冻结文件与通过短窗的当前文件逐项
SHA256 完全一致。`campaign.json` 记录启动身份与前置证据，`annual/progress.json`
及 `annual.log` 为实时进度，`annual/result.json` 为最终验收；本记录写入时全年尚未完成。
受保护四个站点文件及旧 bridge campaign 哈希保持原值。

`compass-coding-standard` 在本次修正中落实为先复现、局部 C/N 物理闭合、保留非目标
PFT 行为、阻止跨版本恢复和独立年度复验，未通过放宽容差或删除历史证据获得通过。

### 5.4 新年度基准完成

2026-09-30 后续：新版本 sample 41 已完成 **8760/8760 h**，最终 `result.json`
为 `passed`，耗时 2765.6604 s。完整年度物理 C 残差
`-1.7180701306074297e-14 g C/cm²`（`-1.7180701306074297e-10 g C/m²`），
N 残差 `4.0202216555762504e-16 g N/cm²`（`4.0202216555762504e-12 g N/m²`）。
每步 `xsmr` 与累计核销严格为零；8760 步 DIN 数值门完整覆盖，NH4/NO3 最大残差
`3.637978807091713e-11` / `2.3987922759260982e-11 ug N/cm²`。
证据：`outputs/changwu_crop_carbon_funding_20260930/annual/result.json`。

另对最终 checkpoint 做零积分审查：`crop_live=false`、`crop_phase=4`，收获 C/N
分别为 `378.5936478` / `10.81681624 g/m²`，根残体 C/N 分别为 `44.81825324` /
`1.067101268 g/m²`。按冻结配置逐 8760 个小时间隔调用生产管理匹配器，矿质肥和
有机肥均各发生一次，与最终接受态的事件标志完全对应。未把匹配器遍历称为物理积分。

本年度的校准/验证 RMSE：水分 `0.0346615 / 0.0355243 m³/m³`；NH4
`7.94351 / 5.82082 mg N/kg`；NO3 `8.11763 / 6.98860 mg N/kg`。
校准各 45、验证各 30 个深度—日期配对。数值守恒通过与观测拟合达标是两回事，
NH4 拟合仍差，不借此宣布率定成功，也不据验证结果选择预热分支或参数。

## 6. 有界自然预热诊断：下一阶段合同

Editor: Wenjie Liu; Date: 2026-09-30.

用户要求年度完成后继续下一步。新年度物理预算门已过，下一项只做 **baseline 41、
最多三个自然年度周期** 的代理驱动诊断，不扩大到 44 点、低/高初态分支或正式优化。
`tools/diagnostics/runChangwuBiogeochemistryPilot.m` 复用现有公共运行、周期交接、
独立预算、分库向量和实测约束初态接口，不修改 `src` 或实测输入。

第一周期从同一 materializer 的 fresh 基准出发，在第一步之前绑定预热身份；之后
完整交接接受态与必要记忆，重放年度驱动及管理日历。原累计账本保留，各周期收支
取端点增量；仅第一周期的预算启用 `start_is_fresh`。不逐年重施实测值、不加酶、
不重分配微生物、不做慢库加速。此前年度是准入证据，因未在 runtime 捕获预热身份，
不追认成这次 pilot 的首个周期或初态供体。

每步检查 DIN 数值门、零欠账/核销；每 10 KT 与终点独立核对实体 C/N 预算。
每周期保存独立分库及各相的有单位原始漂移、真实外部通量、水热原生末态、作物与
管理事件以及 checkpoint；失败保留最后接受态并停止，不放宽门限继续。

科学绝对收敛阈值尚未批准。本批使用显式 `abs_C=abs_N=rtol=0` **仅作精确不变
诊断**，提取变化量，不将这一零门、内部 `state_gate_pass` 或无穷比值用作科学准入。
MAT 保留原始 Inf；JSON 显式写字符串 `"Inf"`，不伪装成零或无说明的 null。
顶层与逐周期均标记 `scientific_convergence_evaluated=false`、
`convergence_verified=false`、`scientific_release=false`、`convergence_status=not_assessed`。
不会因预算跑满或某一次状态接近就自动宣称平衡。

最终只导出 **NOT_READY 的实测约束候选初态**，绑定真实供体及来源。默认正式应用
必须拒绝它，显式诊断应用也不升级放行状态。正式率定仍需科学阈值/收紧检验、
同相位水热与通量稳定、初值分支依赖和 SOC 条件化后响应等证据。

### 6.1 技术验收与启动

`verifyChangwuBiogeochemistryPilot` 的 **22 项检查 PASS，5 个真实积分小时**：
两个各 2 h 的技术周期验证完整 Y 逐位交接、管理重启和周期/累计实体预算；
未放行初态默认应用拒绝、显式诊断应用保留实测及 fresh 账本，并再通过 1 h 物理预算。
已有输出、残留 resume 环境、超过三周期或零步窗口均被拒绝。
证据：`outputs/bgc_pilot_check_20260930/verification.json`。

短窗之后仅修正原生环境便捷摘录中一行字段名：`Theta_II → Theta_I`、`P_g → p_gas`。
完整 checkpoint/carry 原本已保存这两个真实状态，没有原始证据损失。反向替换该行的
SHA256 与通过 22 项测试的 driver 完全一致；零积分读取真实 checkpoint 验证
`Theta_I` 为有限 `58×2`、`p_gas` 为有限 `59×1`，`checkcode` 零提示。
未重复物理积分来测试这个纯输出字段修正，未将旧 driver 的 SHA 冒称新 SHA。

2026-09-30 11:10（Asia/Shanghai）已后台启动
`outputs/changwu_bgc_pilot_20260930/pilot`，**3×8760 = 26280 h** 为硬预算，
900 个冻结文件与当前验证文件逐项 byte hash 一致；生产源码、站点、模型库和率定
工具身份与刚完成的年度基准完全相同，只有专用 pilot driver 不同。
`campaign.json` 保存年度准入证据、短窗及纯输出变更边界、冻结身份和运行合同。
实时进度为 `pilot/progress.json`，完整日志为根目录 `matlab.log`；最终状态以
`pilot/result.json` 为准。本条记录为启动交付，不是三个周期已经完成的报告。

## 7. 2026-09-30 冬麦作物生命周期修正

Editor: Wenjie Liu; Date: 2026-09-30.

用户在正式推进前要求核清植被是否仍为苹果，随后明确授权修正。运行站点确实选择
独立 Winter Wheat PFT，冬麦物候、器官分配和收获已启用；但此前共用路径仍包含
树木的下一年储氮需求和木质 NSC 目标，且缺少灌浆器官 N 再转运与收获 mobile N
去向。旧年度物理守恒通过不证明这些作物生理假设正确，也不追认旧结果为新机制结果。

旧三周期 pilot 于 **11:52:11 +08:00** 精确停止，持久进度为 cycle 1 / KT8490，
真实最后可读取 checkpoint 为 KT7920（不是进度文件中过时的 KT0 指针）。没有删除
历史输入、日志、状态或冻结源码。详情见
`outputs/changwu_bgc_pilot_20260930/STOPPED_FOR_CROP_CN_REPAIR_20260930.json`。

### 7.1 本次局部物理合同

1. 冬麦根吸氮需求为 `max(0, D_growth - N_labile/dt)`，非活动期为零，不再加入
   树木 sapwood/next-year-regrowth 储备项。`N_labile` 仍是有质量的季内移动 N 库。
   非冬麦保留原式；本次不改变无机氮转化、输运或边界方程。
2. 冬麦不加载多年生根/茎木质目标储备；剩余同化碳留在 `C_labile`，继续支付呼吸。
   不是“作物不存在 NSC”。新的 fresh 轨迹四个多年生 C/N 储备应始终为零；六个
   perennial onset-transfer 库若被注入非零值则拒绝，不静默延续。
3. 在 accepted phase < 3 → pending phase 3 且 crop_live 的唯一边沿，以本步起点
   实际器官库计算 `dN_i=max(0,N_i-C_i/CN_fill_i)`，本步生长后器官扣除与 mobile
   N 贷入同时发生，下步先使用此 N 再向土壤请求。C 不随此 N 转运而凭空消失。
   叶/茎/细根参考 C:N 为 **65/100/40**，明确为 CoLM 参数默认值而非长武率定值。
4. 仅冬麦灌浆期的既有整器官 C:N 上界使用 `max(base_CN+offset,CN_fill)`，防止
   正常再转运后既有组织被错误视为超限。新组织 N 需求仍使用原基础 C:N。
   物候 `grnfill=0.40`、播期、LAI/株高/SAI 的 data 模式均不改变。
5. 收获执行一次，沿用已有 accepted `turnover_state.harvested_this_year`，下一次
   live pre-grain 阶段重新武装。种子配对保留：
   `seedC=min(initial_seed_c,C_grain,N_grain*CN_seed)`、`seedN=seedC/CN_seed`，
   `CN_seed` 保持既有叶 C:N 20，不凭空补齐不足的种子 C 或 N。
6. 继续本站地上秸秆及非保留籽粒离田、根残体入土的管理合同；收获剩余 mobile N
   归入离田植物 N，茎储备归离田、根储备归根残体，分别计入已有实际收获/凋落账本。
   `C_labile` 仍沿用此前明确的 crop cleanup 呼吸去向。**mobile N 离田是 COMPASS
   的年度管理适配，不是 CoLM 原版 retransn 清零公式，也不是新增实测管理证据。**

### 7.2 来源及适配边界

固定一手源码：CoLM-SYSU/CoLM commit
`626347a9156ab97db4ca846f47cd1b7f2728f705`。

- [NutrientCompetition:418–443、506–532](https://github.com/CoLM-SYSU/CoLM/blob/626347a9156ab97db4ca846f47cd1b7f2728f705/main/BGC/MOD_BGC_Veg_NutrientCompetition.F90#L418)：一次性器官 N 再转运、已有 retrans N 先满足当步需求。
- [MOD_Const_PFT:1482–1525](https://github.com/CoLM-SYSU/CoLM/blob/626347a9156ab97db4ca846f47cd1b7f2728f705/main/MOD_Const_PFT.F90#L1482)：冬麦为 0-based index 21，叶/茎/根灌浆 C:N 65/100/40；同文件 fcur2 冬麦为 1。

CoLM 原式用 `C*(1/CN_base-1/CN_fill)`；COMPASS FlexibleCN 的实际 N 不保证等于
`C/CN_base`，故改用实际 N 库，并限制为单向转出。CoLM 根 42→40 会使原式为负，
本实现不将该负数误作可供给 N，也不无来源地增加根 N。`N_labile` 代替独立
`retransn` 承载季内移动 N，属明确适配，不宣称完整 CoLM 逐行等价。

### 7.3 验证与版本隔离

冬麦 checkpoint 增加 `winter_wheat_crop_lifecycle=annual_reserve_grainN_v1`，
保留 `physical_nsc_v1`；旧生命周期 checkpoint 即使末态储备为零也不能接续。
先冻结修前 Apple/NET 的 10 组完整 demand/dispatcher/turnover 固定输入输出，
再检验非冬麦逐字段不变；冬麦单测覆盖内部转运、极限库存、阶段边沿、收获锁及
种子双元素限额。验证计数、正常退出状态和耦合结果在实际完成后追加，不能以前一版
106 项或旧年度 PASS 替代本次证据。

增强年度 driver 逐 accepted step 记录植物储备、mobile 和籽粒轨迹与最大值，不再
仅靠年末零库判定季内行为；真实收获以接受态累计收获账本增量识别，并检查非种子
库存去向。通过短窗与事件回归后，仅验证新基准 41 的完整年度，不自动重跑 44 点、
启动三周期预热或进行参数优化。实测水分、DIN、SOC 和站点输入保持不变。

### 7.4 局部验收结果

`verifyWinterWheatCropLifecycle` **73 项 PASS**，含 10 组 Apple/NET 修前真实函数
输出逐字段一致、冬麦不再申请树木储氮/不加载木质储备、灌浆边沿一次性与重试确定性、
下一步再转运氮用于籽粒、器官 N 双边转移守恒、收获所有非种子库存归账、种子 C/N
共同限额，以及六个非法多年生 transfer 库逐项拒绝。
证据：`outputs/winter_wheat_lifecycle_verify_20260930_r2/verification.json`。
原供碳回归 **106 项 PASS**，证据：
`outputs/winter_wheat_funding_postlifecycle_20260930_r2/verification.json`。
两组顺序执行 MATLAB 退出码 0；本轮局部验收不含土壤模型积分。

首次旧供碳回归发现空可选冬麦开关 `[]` 被新增逻辑误作标量布尔，已改用
`isequal(flag,1)`，随后上述两组完整通过；失败记录保留，未放宽物理或数值门。
单独少数 MATLAB 静态检查曾在输出后发生 ddux 退出清理线程异常，最终通过记录
以正常退出的完整回归为准，不将异常退出冒称成功。

### 7.5 耦合、重启验收与新年度启动

本次修正版 **10 个真实积分小时**的重启/周期交接验收 PASS，普通重启与 replay
重启均保持接受态逐位一致，7 类非法合同被拒绝；MATLAB 正常退出码 0。
证据：`outputs/winter_wheat_lifecycle_restart_20260930/verification.json`。
此处短周期只是技术交接验证，不代表年度预热收敛；不同初值协议的历史 2 h 结果
明确为不可比较，未据此宣称历史轨迹一致。

随后 **48/48 h 耦合短窗 PASS**，MATLAB 正常退出码 0；独立实体 C/N 残差分别为
`-5.39499001028787e-16 / -1.2533377113932431e-16 g/cm²`。NH4/NO3 逐步闭合
最大残差为 `3.865352482534945e-12 / 5.1034731995969196e-12 ug N/cm²`。
49 个初态及接受态观察点中四个多年生储备、呼吸欠账及核销均为零，初值实测核对和
入口/出口身份一致性通过；driver 静态检查零提示。证据：
`outputs/winter_wheat_lifecycle_short_20260930/result.json`。
短窗尚在出苗前，灌浆与收获事件均为零，不能据此宣称完整季节过程已通过。

**2026-09-30 13:15:57 +08:00** 启动新 sample 41 单条 8760 h 年度基准，launcher
PID 28660。独立冻结根为 `outputs/changwu_crop_lifecycle_20260930/frozen`，900 个
文件逐项 byte hash 与已验证版本一致；初始化后的源码、站点、模型库、率定工具、
专用 driver、design 六类身份均与预期相符。首次交付核查已达到 KT310，无失败。
实测水分、DIN、SOC、未观测有机库猜测及历史 campaign 原文件哈希未变。

运行清单与前置验收证据见 `outputs/changwu_crop_lifecycle_20260930/campaign.json`；
实时状态为 `annual/progress.json`，日志为 `annual.log`，最终以 `annual/result.json`
为准。本记录只确认已启动，不宣称年度通过；年度驱动还须核验一次真实灌浆边沿、
器官实际扣 N 与 mobile N 贷入匹配、真实收获后的非种子库存及全系统守恒。
本次不重跑 44 点，不自动接续旧 pilot，不启动正式预热或参数优化。
