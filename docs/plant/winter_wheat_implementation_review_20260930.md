# 长武冬小麦地上植被修订与验证

> **版本状态：2026-09-30用户缩小范围后的历史实现记录。** 本文保留回退前的交付范围、来源及当时验证数值。按用户最新范围，APSIM茎WSC、站立枯叶、实际叶氮LUNA反馈及新增冠层fail-fast属于撤回项；地上保留冬小麦物候、CoLM植被C/N分配和原CoLM凋落物，其它沿用COMPASS。正文的完成表述和历史测试不能解释为当前生产状态；当前回退状态和后续检查以[本轮范围回退报告](winter_wheat_scope_rollback_20260930.md)为准。

Editor: Wenjie Liu; Date: 2026-09-30.

## 结论与交付范围

本次审阅后按用户授权实现了观测 LAI 收获、实际绿叶氮对 LUNA 的反馈、实际
PlantCN/AOM/MEND 呼吸、年度冠层日历修正、年度茎 WSC 与站立衰老叶，以及参数
来源和消费者核查。冬小麦不再沿用树木的来年储氮、木质碳储备和树龄冠层规则。
年度 WSC 是当季可移动碳库，不能与多年生储备混用。

当前交付验证的是过程实现、约束及短窗口耦合闭合；长武全年产量、GPP/NEE、叶级
容量的独立数量验证仍需新基准和相应观测。运行中的历史冻结目录
`outputs/changwu_crop_lifecycle_20260930/frozen` 未参与本次编辑。

## 六项修改的实际行为

| 授权项 | 实现与证据 | 科学边界 |
|---|---|---|
| LAI 到零才收获 | `crop_harvest_rule="observed_lai_zero"`；播种、出苗、生理成熟和收获分开；正 LAI→零的下降边沿在能量交换前形成候选，接受后一次性提交 | 当前输入给出的收获日期为 2001-06-09。不会把“6 月底”的文献日期写成现有 LAI 文件的日期。 |
| 实际叶 N→光合 | `luna_leaf_n_source="plant_pool"`；绿叶 N/当前绿叶 LAI；去氮立即约束功能容量，日优化使用接受时钟；休闲和新作季清空旧容量 | LUNA 的全局机制不等于长武品种率定；独立冬小麦数据只能支持目前已核的数量级和反馈方向。 |
| 实际呼吸接入 NEE | 最终接受端点的 `cum_R_a`、AOM 累积 CO2、MEND `IC` 增量合成 `NEE=Ra+Rh_aom+Rh_mend-GPP`；输出各分量 | 使用各库所有者的实际增量，包含实际 N 反馈后的微生物溢出呼吸，内部凋落物转移不另算 CO2。 |
| 年度冠层 | EPIC/热量与年度日历按冬小麦作季处理；无叶/收获步的根吸水候选为零，保留活作物有符号根通量；拒绝尚未闭合的模型株高/SAI 配置 | 当前 LAI、株高、SAI 仍使用原有数据。未把缺少闭合方程的预测配置标记为可用。 |
| WSC 与站立衰老叶 | 新茎生长按结构碳/WSC 分开，实际库存和 N 容量限制 WSC→籽粒，成长呼吸计一次；衰老叶从绿叶转到站立死叶，收获单独移除 | WSC 比例/动员率来自明确的作物模型参考；目前没有粒数与独立潜在粒重模块，不能声称已验证潜在产量。 |
| 参数来源/消费者 | PFT 与 site 覆盖编译为单一只读参数视图；需求、分配、初始库及周转消费一致参数；核查 CoLM、文献和独立叶片数据 | 有出处的参考默认与长武实测分开；缺少出处的光学/水力参数明确保留为待核先验。 |

## LAI 日期与文献的冲突

使用原始气象/LAI 驱动重放得到：2000-09-14 播种、2000-09-15 热量判据出苗、
2001-05-04 进入灌浆、2001-05-25 热量生理成熟（LAI 约 3.901）、2001-06-09
观测 LAI 到零收获。热量出苗不是实测出苗，热量成熟也不是管理收获。
旧默认 `grnfill=0.40` 对应 2001-04-05 灌浆；本次按本地 CoLM 冬小麦
`grnfill=0.67` 更正后为 2001-05-04，旧记录保存在 `colm_prior_parameter_calendar`。

长武长期田间研究说明冬小麦通常在 6 月底收获：
[Contribution of fertilisation, precipitation, and variety to grain yield…](https://www.tandfonline.com/doi/pdf/10.1080/09064710.2016.1149215)
和 [Responses of Wheat Yield…: A 35-Year Experiment](https://www.mdpi.com/2077-0472/12/9/1498)。
用户所述郭胜利论文原文尚未取得，本报告不冒充完成其原文核验。现有 LAI 的 DOY160
对应 2001-06-09，需查原始观测的年、DOY、插值和收获元数据后决定是否更正 forcing。
本次遵照用户指示让收获跟随 LAI 边沿，没有擅自改动观测。

## AOM、凋落物与根呼吸的去向

`AOM` 表示加入的有机物，不是根呼吸这个过程的别名。当前实现中，表面 AOM
分解有直接 CO2；入土 AOM 碳进入 POM/DOM，再由 MEND 的呼吸库记录。凋落物
可进入 DOM 或 AOM，死亡根残体进入土壤有机质。活根的自养呼吸已包含在 PlantCN
的 `Ra_plant`，不能再加入 `Rh_aom`。施入有机肥同样可能提供 AOM/MEND 底物。

最终 NEE 在土壤 N 与根实际吸氮校正之后计算。原经验温度/降水土壤呼吸不再作为
当前完整 C/N 耦合链的 NEE 来源。碳单位为 kg C m^-2 s^-1；`Resp` 为
µmol CO2 m^-2 s^-1，正 NEE 表示向大气排放。

## 参考参数与未解决的配置问题

用户提供的本地 CoLM 文件是冬小麦热量、春化和参考 C:N 参数的核查入口。
参数须核到初始化表及活动消费者，不能仅因同名就复制多年生转移逻辑。
已按同方程/同单位的 CoLM 冬小麦默认更正：

| 参数 | 旧 COMPASS 值 | 已核 CoLM/当前值 | 消费者 |
|---|---:|---:|---|
| `arooti` | 0.20 | 0.10 | 初始根分配及随 HUI 变化的年度分配 |
| `fleafi` | 0.80 | 0.90 | 营养生长初始叶分配比例 |
| `bfact` | -10.0 | 0.10 | 营养生长叶分配曲线的形状 |
| `grnfill` | 0.40 | 0.67 | 灌浆阶段和春化阈值，按 HUI 比例 |
| `mxmat` | 330 | 270 | 播种后的最长生理生育天数 |

参考来源为用户提供的本地 `main/MOD_Const_PFT.F90` index21，以及
`MOD_BGC_Veg_CNPhenology.F90`/`MOD_BGC_Veg_NutrientCompetition.F90` 的活动消费者。
采用这些参考默认不代表长武标定。CoLM `grperc=.11` 是成长呼吸/净新生物量比值，
现有 `R_g_fraction=.25` 是 LPJ/Sitch 的底物份额，两者定义不同，未直接互换。
详细映射、独立叶片数据、来源缺口见
[叶氮与参数审计](winter_wheat_leaf_n_and_parameter_audit_20260930.md)。
年度 WSC 的实验机制、APSIM 参考值、单位和库存守恒见
[年度叶片与移动碳合同](winter_wheat_annual_physiology.md)。

当前 `SLA0=0.035 m²/g C` 在已核的 119 片独立冬小麦叶片范围内；
`Cdm=0.012 g DW/cm²` 对应 LMA=120 g DW/m²，与同一 SLA 隐含 C/DW=0.238。
独立数据 C/DW 为 0.420–0.474，LMA 最大约 94.87。若两参数表示同批叶片，则存在
实质不一致。用独立均值换算的 Cdm 约 0.006356 仅是跨站候选，没有作为长武实测替换。
CoLM BGC 参数不能直接补齐 fluspect 的 Cw、N_mesophyll、LIDF、叶肉 CO2 导度的
长武观测来源；这些仍须相应光谱、叶角或气体交换数据。

## 检查记录

过程报告均保留原容差，不以降低标准换取通过。

| 验证 | 当前证据 |
|---|---|
| 实际绿叶 N/LUNA、生化消费者及参数覆盖 | 新 CoLM 默认后重跑 39 项通过；`outputs/review_aboveground_implementation_20260930/luna_feedback_after_colm_defaults/report.json` |
| 年度冠层/两作季/无叶水力/候选提交 | 新 CoLM 默认 38 项通过；`colm_parameter_canopy_final/verification.json`；旧参数 32 项报告另存 `colm_prior_parameter_calendar` |
| 年度 WSC、站立衰老叶与 C/N 守恒 | 47 项通过；`outputs/review_aboveground_implementation_20260930/annual_physiology_final/verification.json` |
| 原生命周期及碳供给回归 | 73/106 项通过；`original_lifecycle_final`、`original_funding_final`；原冻结输入/参考不变 |
| 实际碳交换及公开输出消费者 | 12 项通过；`outputs/review_aboveground_implementation_20260930/actual_carbon_exchange_final/verification.json` |
| 实际吸氮回写 | 29 项通过；`outputs/review_aboveground_implementation_20260930/actual_root_N_final/verification.json`；冬小麦以实际土壤吸氮从同一接受态重算整候选 |
| 联合灌浆/收获和真实 checkpoint | 120 个接受步通过；`outputs/wheat_coupled_aboveground_20260930_final/verification.json`；两窗口完整 Y 及 8 项碳通量逐位一致 |
| 真实初态短运行 | 48 小时通过，MATLAB exit0；`outputs/wheat_aboveground_fresh48_20260930_final/result.json`；入/退出源码及输入身份一致，非全年验收 |

短联合窗口使用人工正作物库存和冷土壤初态，以直接覆盖灌浆、衰老及收获事件。
它没有包含前期真实土壤/作物轨迹，不能当作该季田间模拟。真实 `beginRun` 重启检查
比较完整接受态，保留原先失败证据，不改为更宽容差。最终联合窗口的最大残差为
C=3.073225e-15 g/cm²、N=2.370066e-16 g/cm²、DIN=5.137668e-12 µg/cm²，
NEE 分量闭合误差为零。灌浆窗口有正的实际根系吸氮，并用独立的层端点/层厚积分
检查实际吸氮率与植物累计账本；收获窗口检查一次性收获、当步无叶交换及新库清理。

联合检验发现并修复了两个接受合同问题。根请求元数据曾在候选物化时丢失，实际
吸氮汇总又先加层再乘层厚，产生了两物种向量而非正确的地面积标量。现在先按层
加物种/平均端点，再按该层厚积分；冬小麦以实际吸氮从同一接受态重算器官分配、
C:N 限制和 WSC 供应，逐项确认已被土壤消费的凋落物/收获通量未变化。非作物候选
物化路径维持历史行为，未把本次年度重算扩大到苹果。

另一差异来自上一接受步的 `Throughfall`/`EVAP` 未写入 checkpoint。它们实际被下一步
突发湿润边界预测消费，不能当作可丢弃诊断。真实断点的 Y、Lambda/EHCAP/KfL/Theta
都相同；仅恢复这两个量就使下一步完整 Y 逐位相同，因果证据在
`outputs/wheat_restart_cache_probe_20260930_r2/harvest/history_causality.json`。
现在仅补保存/恢复这项小型接受记忆，没有保存整个 workspace 或更改水热方程。
4 个损坏历史拒绝检查和 1 个旧文件缺项标记检查通过。

测试使用 `driver` 诊断 profile，因此这些短窗不生成 `carbon_fluxes.nc`。实际 NEE
及各分量在接受的 scope fluxes 上逐步核对；公开 registry 写入另有 12 项消费者检查。
新基准若需要导出冠层/NEE 的 NetCDF，应使用现有 `full` profile。

## 新基准的冬小麦验收方案

1. 在独立输出目录用本次源码从一致初态启动；老 checkpoint 的参数/过程签名不兼容
   时重启拒绝，不能把中途变更拼成一条基准轨迹。
2. 检查作季事件：播种/出苗一次、成熟与收获分开、LAI 到零仅一次收获、休闲 GPP 与
   根吸水为零；核对收获前后各器官 C/N、WSC、死叶、输出/根残体/种子去向。
3. 每步和整季检查 PlantCN+DOM+AOM+MEND+DIN 的碳氮总账及
   `NEE=Ra+Rh-GPP`，再比较连续运行与事件前后重启。
4. 在独立数据上检验叶 Narea/叶绿素→Vcmax/Jmax/A–Ci，再检验 LAI、器官生物量及
   C/N、籽粒产量/氮或蛋白。区分数据驱动 LAI 与预测 LAI，不能把 prescribed LAI 的
   一致性当作物候预测准确度。
5. 用独立站点/年份或留出施氮处理检验响应。WSC 比例和动员率、叶光学参数、水力
   参数做有界敏感性分析后，才决定是否需要长武品种率定。没有观测的数据列保持
   “未验证”，不以参考模型数值替代实验结果。
