# NH4 生物可及性保护阈值

Editor: Wenjie Liu; Date: 2026-11-01.

该机制是在现有线性平衡吸附架构中检验的经验生物可及性假设，不新增铵库、不创造铵，
不改变 NH4 初值分相、干土基观测算子或其液相加平衡吸附相的物理含义。
2006、2019 年 9 月的实测剖面用 2 mol/L KCl 提取铵。这个测量不能单独确定被保护的铵有多少。
[长武原始方法核查](../../data/verification/changwu_observation_provenance_20261103_01/长武处理、观测深度与参数依据.md)

## 配置与参数含义

```toml
[nitrogen]
nh4_protected_per_clay_mg_kg = 0.0
nh4_protection_scope = "nitrification_only"
```

缺省整个节、缺省 `nh4_protected_per_clay_mg_kg` 或显式设置该系数为 `0.0` 均关闭保护，保留原计算分支。
该系数必须为实数、有限、非负数值标量。布尔、字符串、数组和未知 `[nitrogen]` 键报错。
单位为 mg N/kg 干土，数值代表单位黏粒质量分数对应的经验保护量。
`13.95` 是历史 30 号使用过的候选系数。它并非长武直接测得的生理阈值。
以长武 2019 年收获后表土作 28 天培养回放，在较快有机周转且以文中表 1 矿质氮均值初始化时，表层 4.20 mg N/kg 的硝化专用阈值使净硝化低于实测均值减一个小区标准差。两种有机初态下，所试 0–1.5 mg N/kg 的离散点均满足该瓶内氮速率约束。瓶内实际起始氮和含水量没有逐瓶报告，不能据此选出田间阈值。
[同站培养中的保护量、净矿化与净硝化](../../data/verification/cand30_wang_june_assay_20261103_01/measured_DIN/六月长武培养对铵阈值的条件检验.md)
正式长武示例的系数缺省为零。实验需要正阈值时，在该运行独立的站点副本中设置系数与范围。
选 `nitrification_only` 时，阈值只在硝化处生效。选 `nitrification_and_mend` 时，
硝化和 MEND 均使用阈值。旧的正阈值配置缺省范围时沿用后者。
20 cm 原始样点黏粒为 `30.10%`，候选样点阈值约为 `13.95 * 0.301 = 4.20 mg N/kg`；
计算节点的实际阈值须使用下述映射和干质量投影，不直接套用样点数值。

## 映射、单位与节点预算

[`buildSoilParameters`](../../src/+init/buildSoilParameters.m) 已将 `clay_pct / 100` 转为
`s.clay` 质量分数，按站点深度映射到计算层；不可二次除以 100。
[`buildNitrogenSoilParams`](../../src/+init/buildNitrogenSoilParams.m) 使用最终活动容重构建
bottom-up `NL x 1` 只读参数：

\[
S_l = p\,f_{\mathrm{clay},l}\quad[\mathrm{mg\ N/kg}=\mu\mathrm{g\ N/g}],\qquad
M_l = \rho_{b,l}S_l\quad[\mu\mathrm{g\ N/cm^3\ soil}].
\]

节点保护面积质量为各相邻半单元的 `M_l * dz_l / 2` 之和；节点干土质量用同样的
`rho_b,l * dz_l / 2` 投影。两者相除得到节点干土基等效阈值，与
[`observationProfile`](../../src/+soil/+nitrogen/observationProfile.m) 的干质量投影同构。
两种范围都按共享节点上的剩余物理质量限制硝化。采用仅硝化范围时，
MEND 动力学使用根吸预留后剩余的液相底物，实际固持仍受节点供给量和原有守恒检查限制。
采用历史共用范围时，保护量也从 MEND 的动力学可用底物和固持额度中扣除。

## 活动范围与对照口径

支持经验矿化模式，以及显式微生物开启时的 `source_hourly` MEND。
非零阈值与活动显式微生物的其他反应模式组合会报错；零阈值保留原兼容分支。
根系吸氮、氨挥发和输运均可能使局部铵低于阈值，仅硝化范围下微生物也可以继续固持可用铵。

E1a 必须分别与各自模式的 E0 比较。既有 MEND 基线含 120 kg N/ha 尿素，以及
1326 kg C/ha、147.75 kg N/ha 有机肥（其中 TAN 7.3875 kg N/ha）；经验基线关闭了整个
有机肥入口，仅尿素生效。原有 75 Mg/ha 是材料施用量标签，不另行认定干重或湿重。
跨模式差异包含有机肥输入和有机氮周转差异，不能归因于矿化公式的单独作用。

## 输出与检查点

正阈值额外输出接受态的 `nitrification_protection_reduction_profile`、
`nitrification_protection_reduction_flux` 和 `nh4_protection_budget_profile`；
零阈值不注册这三个字段。`nh4_substrate_total` 继续表示原物理总底物，
正阈值的 MM 因子使用扣除保护量后的可及底物。
削减剖面及通量仅表示**阈值 MM 速率之后，共享节点预算造成的额外削减**，
不包含点态可及底物减阈值本身引起的 MM 减速。E1a 与 E0 总硝化量之差还包含
底物供源变化及 MEND 反馈，不等于上述额外削减通量的年度积分。

零阈值不增加检查点兼容签名字段，其他契约相同时允许原 schema-2 检查点续跑。
非零阈值签名绑定站点标量、实际层阈值和适用范围。历史共用范围保留原有签名，
仅硝化范围另有明确的签名字段。两个活动范围的检查点不能互相续算。

## 验证入口

- [`tests/testNh4ProtectionThreshold.m`](../../tests/testNh4ProtectionThreshold.m)：
  验证保护量单位、共享节点预算、零值兼容性和检查点适用范围。
- [`scripts/verify_nh4_protection_windows.m`](../../scripts/verify_nh4_protection_windows.m)：
  `verify_nh4_protection_windows(codeRoot,outputRoot,mode,variant,nSteps,baselineFile,resumeFile)`；
  `mode` 为 `mend`/`empirical`，`variant` 为 `default`/`zero`/`protected`（13.95）。
  每个窗口使用独立 MATLAB 进程和新的输出目录，只修改输出目录内的站点副本。
  `nSteps` 为目标绝对 KT。执行范围见下列日期记录。

以上入口及参数说明不构成测试已通过的声明。验收须保留各模块原生 C/N/DIN/O2 与
示踪容差，包含 O2 子步/重放门和接受态每步、累计严格门；短窗通过也不证明全年响应
或该经验阈值的科学有效性。

## 2026-11-01 执行记录

聚焦单元合同、缺省/显式零值各模式24小时原快照逐位数值对比、13.95阈值各模式24小时原生守恒门，以及24→28小时的旧检查点零值续跑比对/非零阈值续跑，均已通过。验收原件见 `data/output/changwu/nh4_protection_20261101_001826/validation_gate.json`。

两组E1a全年已按冻结代码副本启动，结果与自动评价位于 `data/output/changwu/nh4_protection_20261101_001826/`。全年完成前不据短窗结果宣称科学假设已验证或NH4失配已修复。原站点输入未改，当前试验只增加独立配置的保护标量13.95。

## 2026-11-03 仅硝化范围验收

正常站点入口在同一施肥窗口分别运行零阈值、正阈值仅限硝化、正阈值显式共用和旧版省略范围四种配置，各接受 6 小时。首小时正阈值两种范围的硝化均为 0.049843 kg N/ha，零值为 0.109624。仅硝化配置与零阈值的首小时全层 MEND 毛铵请求逐值相等。旧配置与显式共用的 1,433 个 NetCDF 数据集逐值相等。另在存在根吸铵的 2001 年 5 月真实断点上重演根后供给。65–70 cm 的物理余量为 5.8994 mg N/kg，仅硝化模式供给 MEND 的额度为 5.8994，共用模式为 1.5772。[四组公开入口与历史根吸小时](../../data/verification/nh4_scope_release_20261103_01/public/receipt.json)

对部署范围另核对 862 个运行文件，只有此次阈值接线涉及的四个源码文件不同，21 个正式站点文件与归档相同。取紧邻本次变更的源码快照，和当前正式入口分别以同一 MATLAB 单线程运行零阈值示例 6 小时，23 个 NetCDF 文件的 1,417 个变量数组及掩码逐值相等。两种正阈值之间的检查点互换已在聚焦用例中验证拒绝。以上检查覆盖程序行为，不能确定长武的保护系数。[正式默认站点的改前改后对照](../../data/verification/nh4_scope_release_20261103_01/default_parity.json) · [聚焦测试入口](../../tests/testNh4ProtectionThreshold.m)
