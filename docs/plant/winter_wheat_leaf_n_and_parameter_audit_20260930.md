# 冬小麦叶氮反馈、适用性验证与参数消费者核查

> **版本状态：2026-09-30用户缩小范围后的历史实现记录。** 本文保留回退前的实现、文献来源和当时验证数值。按用户最新范围，APSIM茎WSC、站立枯叶、实际叶氮LUNA反馈及新增冠层fail-fast属于撤回项；地上保留冬小麦物候、CoLM植被C/N分配和原CoLM凋落物，其它沿用COMPASS。本文旧39/38项检查不能解释为当前生产验证；当前回退状态和后续检查以[本轮范围回退报告](winter_wheat_scope_rollback_20260930.md)为准。

Editor: Wenjie Liu; Date: 2026-09-30.

## 本次实现与结论

长武的 `luna_leaf_n_source="plant_pool"` 现在使用接受态绿叶氮除以当前交换冠层的绿叶
LAI，即 `N.leaf / LAI_green`，而不再用固定 `1/(SLA0*CN_L)` 代替实际叶氮。
`N.leaf` 与 `C.leaf` 是绿叶库；站立衰老叶、籽粒种子和根不贡献光合酶或叶绿素氮。
本次未把树木的来年储氮、木质储备重新接回冬小麦。

灌浆或衰老移出的叶氮会立即限制同一步的生化容量，优先于 LUNA 的缓慢日驯化。
容量满足

    N_lc + N_resp + Vcmax25/Fc25 + Jmax25/Fj25 <= FNCa <= N.leaf / LAI_green.

结构氮仍沿用 LUNA 的 SNC 闭合，既有 FNCa=15 g N m^-2 leaf 上限保留。
投影只调整辅助酶容量/光学输入，不修改任何物理 C/N 库或守恒账本。日优化只在接受
外层步后提交，使用实际接受时钟的天文日长。无叶、收获休闲及新的作季清空旧作物的
酶容量；绿叶重新出现后用新叶实际氮进行优化。候选构造与序列化恢复检验覆盖了该记忆。

另修复了一处确定的光谱适配问题：`RTMo.rad.PAR` 已是入射 PAR，`computeLuna` 原本
将吸收 PAR 除以吸收比例来还原入射光。因此实际叶氮模式传入 `rabsorb=1`，表示输入
无需还原，绝不表示真实叶片吸收率为 1。热红外反射率/透射率不再作为 PAR 的吸收率。
显式 `pft_constant` 模式保留历史适配，用于其他 PFT 的兼容对照。

LUNA 的总叶氮与气候驱动机制来自
[Ali et al. (2016), GMD 9:587–606, DOI 10.5194/gmd-9-587-2016](https://gmd.copernicus.org/articles/9/587/2016/)。
该论文不是长武冬小麦标定；本次容量预算投影是 COMPASS 为实际 C/N 耦合添加的约束。

## 已执行的验证及其边界

`tools/diagnostics/verifyWinterWheatLunaFeedback.m` 的最终证据在
`outputs/review_aboveground_implementation_20260930/luna_feedback_final/`。
39 条纯过程检查通过，包括实际氮浓度输入、N 移除后的容量预算、无叶极限、新作季
重置、接受态不变、完整气候缓存恢复、实际日长、热红外参数独立性、活动 C3 生化消费，
以及初始化/根氮需求对显式 C:N 参数的实际消费。没有用该测试声称全年田间准确。

固定合成气候和水力条件时，叶氮响应如下。每项容量日优化 60 次；生化检验使用
Tleaf=Tgrowth=25°C、吸收光 1500 µmol photons m^-2 s^-1、Ca=400 ppm、诊断 Ci=300 ppm、
水力因子 1、当前小麦气孔/叶肉参数。它隔离 N 反馈的方向，不是田间拟合。

| 叶 N (g N m^-2 leaf) | Vcmax25 (µmol m^-2 s^-1) | Jmax25 (µmol m^-2 s^-1) | 净 A (µmol m^-2 s^-1) |
|---:|---:|---:|---:|
| 0.05 | 1.0523 | 1.7889 | 0.2248 |
| 0.30 | 11.0931 | 20.3045 | 2.3622 |
| 0.75 | 21.7924 | 40.0476 | 4.5327 |
| 1.50 | 37.0223 | 69.6658 | 7.4720 |
| 3.00 | 63.5863 | 128.7291 | 12.2147 |

独立冬小麦田间研究也表明叶总 N 与 Vcmax25 的关系受到氮分配和物候影响。
[Li et al. (2023), DOI 10.3389/fpls.2023.1225295](https://www.frontiersin.org/journals/plant-science/articles/10.3389/fpls.2023.1225295/full)
在禹城冬小麦观测中报告 Vcmax25 与叶 Narea 的 R²=0.36、与叶 Chlarea 的 R²=0.69。
这支持将实际 N 与功能氮分配/叶绿素区分，但不能据此宣称本模型已达到这些拟合度。

## 独立原始数据核查

已获取并保留
[Fang et al. (2022), 4TU 数据集 DOI 10.4121/21304917.v1](https://doi.org/10.4121/21304917.v1)，
相关论文为
[Fang et al. (2023), Global Change Biology 29:505–521, DOI 10.1111/gcb.16488](https://pubmed.ncbi.nlm.nih.gov/36300859/)。
该 CC BY 4.0 数据集来自 Wageningen 的 2019–2020 年温室/气候箱冬小麦
Thésée、Récital，含充分供水/干旱与生长温度处理。
论文 Methods 的气体交换/荧光小节明确叶材料在 70°C 烘至恒重，随后测质量和 N；
因此下表使用的是干质量，而非把未区分的鲜质量擅自当作干质量。

原始文件在 `outputs/review_aboveground_implementation_20260930/leaf_evidence/`，未改写。
`publisher_metadata.json` 保留发布方元数据；`source_file_manifest.json` 同时记录下载 URL、
发布方 MD5、实际验证 MD5、SHA-256 与字节数。三个原始数据文件全部通过发布方 MD5 检验。
`fetchWinterWheatLeafEvidence.py` 可复取固定 v1；`auditWinterWheatLeafTraits.py` 用原始
叶质量、面积和 C/N 百分数独立计算下表，保留 Excel 行号。完整 119 条记录和 10 个处理
组汇总作为证据 CSV，未导入长武的 forcing、参数或观测。

| 独立实测/单位换算量 | 最小值 | 均值 | 最大值 |
|---|---:|---:|---:|
| 叶 C 质量分数 (g C/g DW) | 0.41982 | 0.44949 | 0.47383 |
| 叶干质量面积比 LMA (g DW/m² leaf) | 40.2375 | 55.6559 | 94.8695 |
| SLA 按 C 质量换算 (m² leaf/g C) | 0.023202 | 0.041554 | 0.057962 |
| 叶 N 面积密度 (g N/m² leaf) | 1.07643 | 2.11035 | 3.03262 |
| Cdm 按干质量换算 (g DW/cm² leaf) | 0.004024 | 0.005566 | 0.009487 |

此独立范围支持检查当前 SLA 的单位和数量级。气体交换工作簿没有直接给出 Vcmax/Jmax
列，LeafCN 工作簿也没有重复编号，故不能伪造逐叶配对或把重新拟合的 apparent Vcmax
当作论文报告的、考虑有限叶肉导度的真实容量。没有报告缺少输入条件支撑的 LUNA RMSE。
论文的温度/干旱驯化结果还提示，即使 N 反馈正确，也需验证热响应、叶肉导度和水力限制。
Fang 等发现这些处理下的生化驯化不能仅由叶总 N 变化解释；本次修复叶 N 的实际消费，
未补造独立的冬小麦干旱生化驯化方程。

## 当前参数的来源与活动消费者

`compilePlantParameters` 把数值 PFT 表与显式 site 覆盖编译成一个只读 view；源标签区分
`PFT_table`、`site.toml:vegetation` 和原有的具名兼容常数。需求、碳分配、周转与 LUNA
消费同一组已编译数值，避免初始化使用覆盖值、后续又回到 PFT 原值。无效覆盖会报错。

| 项 | 当前值/来源等级 | 活动消费者和解释 |
|---|---|---|
| SLA0、SLAm | 0.035 m²/g C、0；已核 pinned CoLM 的冬小麦 `slatop`/`dsladlai` 数值，属于参考默认 | PlantCN 叶 C→LAI 与叶生长上限；固定 N 模式使用其计算 Narea。实际 N 模式在有绿叶时使用当前交换 LAI/C.leaf 推导面积质量比。当前实测 LAI 模式不以它替代 LAI 观测。 |
| CN_L、CN_target_root | 20、42；已核 CoLM `leafcn`/`frootcn` | 初始 N、器官生长 N 需求、C:N 上限；实际光合 N 浓度来自 N.leaf，不被目标 CN_L 固定。 |
| CN_target_stem | 50；参考默认 | 是否有 stem 分配的门和兼容 fallback。当前有独立 CN_lw/CN_dw 时，它不是 stem 数量方程中的唯一目标；仅改此项的正值不会改变活/死 stem 的目标比值。 |
| CN_lw、CN_dw、flivewd | 50、500、1；参考/既有 PFT 设置 | `buildInitialPlantPools`、`computeOrganDemand`、`computeCarbonAllocation` 分别消费活/死 stem 比值。冬小麦 flivewd=1，无新 heartwood；调年度结构茎的 C:N 应覆盖 CN_lw。 |
| CN_target_fruit | 50；已核 CoLM `graincn` | 冬小麦的 fruit 槽位表示 grain；消费于初始籽粒 N、生长 N 需求及 grain C:N 约束。名称不代表苹果果实逻辑被启用。 |
| CN_leaf/stem/root_grainfill | 65、100、40；已核 pinned CoLM `fleafcn`/`fstemcn`/`ffrootcn` | 供体可移出 N 的保留界限；不是凭比值生成新的 N。站立衰老叶沿用 leaf=65 是明确标记的 COMPASS surrogate，非长武实测。 |
| Jmaxb1 | 0.0936；LUNA 全局参考参数 | 日 N 优化；冬小麦不消费苹果的 0.244 反推值。尚未做长武品种率定。 |
| Vcmax25/manual | 100.7；Spring Wheat proxy | 非 LUNA 的生化基底。当前实际 N LUNA 路径在交换前覆盖成氮约束容量，不能用此数字证明当前模拟已校准。 |
| g_m | 0.295 mol m^-2 s^-1；C3-crop proxy，冬小麦来源未独立定量核实 | 活动 `biochemical` 的有限叶肉 CO₂ 导度。需冬小麦/品种数据检验，暂不伪造替代值。 |
| Medlyn 参数与 LUNA stomata | 交换使用站点 gs_g1；优化内部仍用 LUNA 的 Ball–Berry mp=9 | 两种模型的气孔参数不互为同义值；当前验收说明算法一致性边界，不能把 mp 当作 Medlyn g1 率定结果。 |
| Cdm | 0.012 g DW/cm²；既有站点配置，原始数值来源未核实 | `buildScope`→`precomputeScopeOptics`/动态 Cab 重算→fluspect 的 Kall；影响光学与吸收辐射，不决定 PlantCN 干质量。当前值相当于 LMA=120，超出本独立样本的最大 94.87。 |
| Cw | 0.009 cm；既有配置，未找到长武原始叶鲜/干质量或光谱来源 | 同上，水吸收光谱参数；不是 steady-SPAC 的器官储水库。当前仍作为待核先验。 |
| N_mesophyll | 1.4，无量纲；既有配置，数值来源未核实 | `buildScope` 映射为 leafbio.N，fluspect 紧凑叶层结构参数；它不是 g N/m²，也不是叶肉 CO₂ 导度 g_m。 |
| LIDFa、LIDFb | -1、0；配置明确为 legacy erectophile 值 | `buildVeg`→`precomputeScopeOptics`→`equations.leafangles`→RTMo/RTMt 的方向分布。未配置实测叶角文件，不能标为冬小麦实测，也不能只因 erectophile 就判为苹果专属。 |
| 热红外 rho/tau | 0.01、0.01；既有配置，数值来源未核实 | 热发射率与热 RT；本次实际 N LUNA 不再把它们当作 PAR 吸收率。 |

CoLM 数值证据固定为
[commit 626347a, main/MOD_Const_PFT.F90](https://github.com/CoLM-SYSU/CoLM/blob/626347a9156ab97db4ca846f47cd1b7f2728f705/main/MOD_Const_PFT.F90)，
数组为 0-based，冬小麦 index=21；本地原样副本 SHA-256 记录于 manifest。
当前 COMPASS 的 SLA 单位由 `LAI=SLA*C.leaf` 和实际 N/C 输入合同确认。
历史 `pft_clm5_audit.md` 有互相矛盾的 SLA dry-mass/carbon 单位结论；此处用当前生产合同
与独立叶 C 百分数重新检验，不据旧审计文字对整个 PFT 表再次乘 2。

## 仍需要的田间检验与具体方案

当前 SLA0=0.035 落在独立冬小麦范围，但若 Cdm=0.012 与该 SLA 描述同一批叶，它们隐含
的 C/DW=1/(0.035×0.012×10000)=0.2381，与上述独立实测 0.4198–0.4738 不相容。
用该数据均值 0.44949 换算的候选 Cdm≈0.006356，属于跨站参考候选而非长武实测。
本次保留原配置，明确建议后续先量/核对长武同批叶片的干质量、C 质量和面积，再统一
SLA 与 Cdm；有叶光谱时用 PROSPECT/fluspect 联合约束 Cw、N_mesophyll、Cab。

冬小麦数量验证需要独立的 Narea、叶绿素、A–Ci/荧光、叶肉导度或已估计的 Vcmax25/Jmax25，
并与同期气候/物候配对。固定光合/水力参数后比较施氮处理的叶 N→容量响应及季节衰老，
再检查 LAI、冠层 GPP/NEE、器官 C/N、籽粒产量和蛋白/氮含量。C/N 守恒与重启一致只能
证明过程和账本实现，不等同于这些独立观测的数量验证。当前缺少对应长武叶级数据，因此
不宣称冬小麦产量或通量已通过田间适用性验证。

## 本地 CoLM 冬小麦默认值复核与采用

按用户“可采用本地 CoLM 冬小麦不确定参数”的指示，只替换与当前活动方程具有相同
单位和含义的 PFT25 默认值。参考目录为
`C:/Users/Administrator/Desktop/CoLM/CoLM-master`；CROP 表中 index 21 是 winter_wheat，
index 22 是 irrigated_winter_wheat，index 15 是 generic C3 crop。本次使用 index 21。
`auditCoLMLocalWinterWheat.py` 只读解析字面量，并对每个数组确认 79 项，保留零起始
索引、声明行、具体数值行、单位、COMPASS 参数和活动消费者。

证据文件为
`outputs/review_aboveground_implementation_20260930/colm_local_audit/local_parameter_crosswalk.json`。
其中 30 项 crosswalk PASS 是来源/数值检查，不是田间模型验证；已记录以下原文件 SHA-256：

| 本地源文件 | SHA-256 |
|---|---|
| main/MOD_Const_PFT.F90 | 0dedc917e1de007a96e484196e82f4c4b647190f9a50af5bce23633e9e9b4f41 |
| main/BGC/MOD_BGC_Veg_CNPhenology.F90 | 56efffcaa788ae7c08c2d9b1ee45b37fc32b4ba8dbe2be2558892419bf1e9778 |
| main/BGC/MOD_BGC_Veg_NutrientCompetition.F90 | c4e006bed29471cb99dcb0bdf74cbd981c042d36ccfbb50f131b4201b0a7c936 |

本地 Const_PFT 与前述固定 commit 的原始下载文件 SHA-256 一致。以下 `L` 是该本地文件
的数值行；原始声明行和每个参数的永久链接另见 JSON。只改冬小麦行，未改其他 PFT。

| 参数 | 旧 COMPASS 默认 | 本地 CoLM / 本次默认 | 数值行、单位和活动消费者 |
|---|---:|---:|---|
| arooti | 0.20 | 0.10 | L1380；无量纲；vegetative/grain-fill 根分配式 |
| fleafi | 0.80 | 0.90 | L1410；无量纲；非根分配中的初始叶比例 |
| bfact | -10.0 | 0.10 | L1425；无量纲；HUI 叶分配指数形状 |
| grnfill | 0.40 | 0.67 | L1334；HUI 比例；灌浆事件、春化截止和分配分母 |
| mxmat | 330 | 270 | L1319；播种后日数；热量成熟判据 |
| baset / lfemerg | 0 / 0.03 | 保持 0 / 0.03 | L1350 / L1304；°C / HUI；积温与出苗 |
| arootf / astemf / allconss | 0 / 0.05 / 1 | 保持 0 / 0.05 / 1 | L1395 / L1365 / L1455；无量纲；根分配终值和灌浆茎衰减 |
| SLA0 / SLAm | 0.035 / 0 | 保持 0.035 / 0 | slatop L1273 / dsladlai L1258；面积/g C；叶 C→LAI |
| CN_L / CN_target_root | 20 / 42 | 保持 20 / 42 | leafcn L911 / frootcn L933；g C/g N；新叶/根 N 需求 |
| CN_lw / CN_dw / flivewd | 50 / 500 / 1 | 保持 50 / 500 / 1 | livewdcn L948 / deadwdcn L963 / flivewd L1228；结构茎活/死分配和 N 需求 |
| CN_target_fruit | 50 | 保持 50 | graincn L978；g C/g N；grain 槽位的 N 需求 |
| CN_leaf/stem/root_grainfill | 65 / 100 / 40 | 保持 65 / 100 / 40 | fleafcn L1486 / fstemcn L1501 / ffrootcn L1516；g C/g N；实际供体 N 保留界限 |
| croot_stem_alloc | 0 | 保持 0 | croot_stem L1198；根/茎分配比；年度小麦不分配木质粗根 |

这些参数的单位和消费方向逐式对照
[CoLM NutrientCompetition L398–454](https://github.com/CoLM-SYSU/CoLM/blob/626347a9156ab97db4ca846f47cd1b7f2728f705/main/BGC/MOD_BGC_Veg_NutrientCompetition.F90#L398)
与当前 `winterWheatPhenology`。当前 `computeOrganDemand` 和 `computeCarbonAllocation`
消费同一 pending HUI 分配系数；CN_target_stem 的正值只是 stem demand 的门和兼容
fallback，真正的年度结构茎比值由 CN_lw=50 决定。deadwdcn=500 因 flivewd=1 不产生
新的 dead stem。根供体原始 42→40 公式可以是负数，本模型仍维持实际 N 库存下的单向
非负转运，未把此来源数值误当作可负向转运的授权。

没有直接替换以下有语义差异或来源缺口的项：

- CoLM grperc=0.11 是增长呼吸/净新生物量的 factor；allocation 的 `1+g1` 在 L487。
  COMPASS `R_g_fraction=0.25` 是 LPJ/Sitch 对 maintenance 后可用底物的 fraction，
  当前 computeNPP 和 WSC→grain 都用这一底物 fraction，后者 yield=1−fraction。
  若要改成 CoLM 增长方案，等价底物比例为 0.11/1.11≈0.099099，而非直接赋 0.11；
  这需要单独验证和科学方案选择。本次保留 0.25，并纠正原错误的 CoLM g1 来源注释。
- CoLM fcur2=1（L1243）令当季新生物量不进入下年 display storage。
  COMPASS 年度路径已显式跳过树木 reserve loading；当前兼容字段 fcur=0.95 不控制
  该路径，不能据它宣称小麦仍储存 5% 木质碳。未导入多年生 storage/transfer pools。
- 检查本地 crop phenology、nutrient allocation、PFT stocks/fluxes 后，没有发现独立的
  站立枯叶库或无结构 N 的 stem WSC→grain owner。grainc_storage 是籽粒器官 storage，
  不能等同于 stem WSC。WSC 分配 0.30/日转运 0.50 保持原 APSIM 机制启发/明确参考默认
  的归属，不能改称 CoLM 原始系数。
- fleafcn=65 的源注释是“grain fill; leaf”。lflitcn=25（L993）消费于
  CNPhenology 的 `leafn_to_litter=leafc_to_litter/lflitcn`（L1356、L1410），它对应实际
  凋落通量，并非站立枯叶 C/N。当前 standing-dead CN65 保持 COMPASS surrogate 标签；
  未以 CoLM 的另一池 C/N 偷换这项。
- allconsl=3、declfact=1.05 虽然在表中声明，但所审 NutrientCompetition 方程未消费它们。
  未给 COMPASS 添加无消费者的参数。CoLM laimx=7 是通用参考；当前显式 site laimx=6.2
  和 LAI 观测保持原始归属。Cdm/Cw/N_mesophyll/LIDF/g_m 不由这份 crop C/N 表证明。

`verifyWinterWheatCanopyLifecycle` 最终通过 38 项检查（0 coupled soil steps），增加了
22 个编译默认值对应同一源表、4 个独立 HUI 分配/实际 N-demand 消费者数量检查和当前
forcing 的灌浆日期门。原有事件一次性、候选重试、checkpoint 和两季 cohort 检查保留。
旧实际日历的 32 项 PASS 证据在 `colm_prior_parameter_calendar/verification.json`；
新证据在 `colm_parameter_canopy_final/verification.json`，没有覆盖旧参数历史。

| 当前 forcing 事件 | 旧参数日历 | 本地 CoLM 默认后的日历 |
|---|---|---|
| 播种 | 2000-09-14 | 2000-09-14 |
| 出苗 | 2000-09-15 | 2000-09-15 |
| 灌浆开始 | 2001-04-05；LAI=3.504 | 2001-05-04；LAI=6.19 |
| 热量成熟 | 2001-05-25 | 2001-05-25 |
| 实测 LAI 正值→0 管理事件 | 2001-06-09 | 2001-06-09 |

这些日期是用当前 forcing 重放得到的模型事件，不能当作长武真实播种、开花或收获
观测。用户指定的 observed_lai_zero 仍优先决定物理收获；未复制 CoLM 把热量成熟
直接作为收获的管理分支，也未改写观测来满足文献的六月底收获日期。
