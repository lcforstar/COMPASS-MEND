# 长武站土壤碳氮参数来源与率定边界

## 当前本地配置：2026-10-12 分层参数与 CSV 自由 DOC 初态

Editor: Wenjie Liu; Date: 2026-10-12.

- 本轮按用户明确要求，将 MEND 8+2 候选参数全部放入 `soil.csv`，并从
  `site.toml` 删除对应十个键。18 个深度行沿用迁移前相同数值，未率定或人为制造深度差异；
  `mend_vertical_parameter_mode="sample_linear"`，水力/无机氮映射保持原样。
- 水力参数采用 2026-10-08 用户接受的 `n_restore` 组；无机氮参数采用
  `nitrogen_priority_rows_20261009_01/trial_00024`。迁移前后原 29 列全部不变。
- 本轮普通初始化读取 `soil_initial_profile.csv`，不导入历史预热 MAT。
  水分、NH4、NO3 和总 SOC 为用户确认的实测支持剖面；网格插值/外推不等于独立逐层实测。
  温度、气压和容重按同一输入文件读取，不扩大其观测来源声明。
- 用户本轮明确选择 CSV 中的 DOC 数值作为自由溶解态 D 初值。因此仅将
  `doc_mgkg` 列名改为 `doc_free_mgkg`，保留所有数值，单位为 mg C/kg 干土。
  这是一项明确的相态含义改变：旧列表示 D+Q，新列表示 D；不能称为旧初态逐位等价。
  在最终初始液态含水量下由吸附关系计算 QOC，从 MOC 中预留 D+Q，并保持各层实测总 SOC。
  DOC 仍为用户指定试算初值，不升级为同期实测或预热收敛值。
- DON 保持原 `don_mgkg` 总量输入含义；未观测有机分库仍按 `organic_initial_pools.csv`
  先验构建。本轮未执行预热。诊断采用 DOC=D、MAOC=MOC+QOC，SOC 账本不重复加入 MAOC。
- 当前服务器冻结的 15 组试验不因本轮本地变更而修改、停止或重启。
  迁移前快照、检查及本地年度运行记录位于 `outputs/local_mend_layer_annual_20261012_01/`。

以下日期记录保留历史配置与来源，不覆盖以上本轮选择。

## 0. 2026-09-29 用户澄清：实测无机氮/SOC 与预热起始猜测值

Editor: Wenjie Liu; Date: 2026-09-29.

- 用户明确并随后补充：初始无机氮（NH4、NO3）和 SOC 均以实测值为依据；
  当前 DOC/DON 为用户任意指定的试算初值，后续需要通过预热建立正式初态，
  不作为初始日期的实测 DOC/DON。
- “随机给的”在此记录为 **initial guesses（起始猜测值）**，不声称它们来自已定义
  的随机分布、抽样种子或独立田间测量。本次不改写或重新抽取任何初值。
- SOC 总量有实测依据不等于 POC/MOC、微生物、酶及各有机氮分库均已实测。
  未观测分库的状态/分配需要预热和必要的先验约束；预热结果属于模型推定，
  不升级为观测值。本澄清不补造 SOC 的采样日期、深度覆盖或外推部分的证据。
- 既有核查发现初始 DON 与 2006 年剖面数值对应，该数值核查记录保留；但不再
  以查找一份不存在的 2000 年 DOC/DON 实测初值作为预热设计的前置要求。
  预热仍需检验对不同起始猜测值的依赖，不能仅改来源标签就宣称消除了初值记忆。
- 正式初始化采用“实测 NH4/NO3/SOC 约束 + 未观测有机库预热”的待设计流程。
  具体驱动/管理循环、收敛准则和模拟起点的联合状态构建尚未冻结；不能将预热
  得到的 NH4/NO3/SOC 直接当作替代实测的最终初值，也不按步强制重置这些库。
  初始化条件化及预热到正式模拟的状态/账本衔接须显式记录，不能靠未记账的
  补池或库间转移隐藏 C/N 差额；预热不能保证唯一确定真实分库。
- DON 初态会影响实测起点之后的无机氮演变，因此最终无机氮率定需要与未观测
  有机库的预热初态一致；这不表示初始 NH4/NO3 缺乏观测，亦不表示硝化、
  反硝化等动力学参数已由初值观测确定。当前五点
  `changwu_watern_bridge_20260929` 保持冻结输入，仅用于相同起始猜测下的新旧
  版本对照，不因本澄清中断、重启或改写；它不是已预热的最终率定实验。
- 本次只更新来源/工作顺序说明；未修改模型、初态 CSV、参数、观测或冻结副本，
  未启动预热。后续 DOC/DON 观测的测定相态和实验口径仍须独立确认。

## 1. 配置归属

2026-10-10 更新（Wenjie Liu）：按用户批准，长武生产配置选择 MEND 原生
`mend_ihr = 1`，使生长与活跃微生物维持呼吸按实际 DOC 摄取分配。
`iHR=0` 与 `iHR=1` 均为原 MEND 选项；这次是配置选择，不是另写本构或
将原 `iHR=0` 定性为 COMPASS 移植错误。分配关系为
`Rg + Rm = (1 - Yg) * U`，休眠和氮溢流呼吸另计。

本次将 `soil.csv` 各行 `mend_ks_dom_mg_g` 从 `0.00010034` 改为
`0.01 mg C/g dry soil`（即 `10 mg C/kg`）。新值是 Jian et al. (2024),
DOI: 10.1029/2023JG007825，补充材料 Table S3 的先验下限，作为首轮
待验证站点先验，尚非长武率定最优值。保留其原温度响应及在摄取、休眠中的
共同使用；它不是 `kd_doc_cm3_g` 平衡吸附系数。原 pH 曲线与其他 MEND
动力学参数不变。下文 2026-09-29 的数值记载为历史迁移记录。

同初态对照试验保留先前冻结的服务器基准初始化模型及其水氮参数，显式替换
上述两个生效参数，另存新模型身份与参数快照；不得将旧参数检查点当作新配置
的普通续跑。当前站点表中其他水氮参数不因本次试验被覆盖。

2026-10-11 更新（Wenjie Liu）：依据 MEND 的 case/soil-series 参数语义及当前
长武缺少分层微生物动力学观测的事实，`mend_vd_pom_h`、
`mend_ks_pom1_mg_g`、`mend_vg_h`、`mend_ks_dom_mg_g` 恢复为
`site.toml [microbe]` 的全站标量。数值保持 50、60、0.00425194、0.01，
本次只改变输入所有者，不改变方程或参数值。逐层 pH、DOM 平衡吸附、弥散及
微生物介质输运参数继续由 `soil.csv` 持有。若未来获得分层动力学证据，适配器
仍支持完整四列 soil profile；site/soil 双定义和不完整定义继续拒绝。

2026-09-29 更新（Wenjie Liu）：固定的站点统一参数在 `site.toml`，需要按层
表示/率定的过程参数在 `soil.csv`，初始库在专用初态文件。CSV 按列名读取，
不依赖列的位置；`#` 行保留中文含义和单位。以下为历史迁移记录。

本次将 `mend_vd_pom_h`、`mend_ks_pom1_mg_g`、`mend_vg_h`、
`mend_ks_dom_mg_g` 从 site 迁入 soil。18 行分别沿用 50、60、0.00425194、
0.00010034 的原先验，不是新的率定值。现已允许各行不同，按站点选定的
`sample_linear`/`horizon_bounds` 合同映射到计算层。

其余 MEND 生理参数和 POM2/MOM 共享关系本次不改。旧站点只有完整 site 标量组
且没有 soil 动力学组时可兼容；双定义或不完整组会报错。完整合同见
`docs/soil/biogeochemistry/MEND_LAYER_KINETICS_20260929.md`。

下文保留历史来源记录；其中旧 MEMC 条目和历史数值不能代替现行配置快照。

## 2. 长武文献约束（固定观测或验收范围）

| 观测量 | 长武文献值 | 在本工程中的用途 |
|---|---:|---|
| 0–20 cm 微生物量碳 MBC | 254.65–745.26 mg C kg⁻¹ | 约束首次初始化与年度 MBC 范围，不直接替代 MEMC 速率 |
| 0–20 cm 微生物量氮 MBN | 31.70–120.54 mg N kg⁻¹ | 约束 MB 的 C:N 与有机氮闭合 |
| 初始 SOC / MBC / DOC | 6.5 g C kg⁻¹ / 113 mg C kg⁻¹ / 22.43 mg C kg⁻¹ | 独立的低肥力初值参照；DOC 质量浓度不能无换算地当作水相 mg L⁻¹ |
| 长武冬小麦表土 SOC / TN / MBC / DOC / pH | 8.9 g kg⁻¹ / 1.2 g kg⁻¹ / 163 mg kg⁻¹ / 32.4 mg kg⁻¹ / 7.9 | 全年结果数量级与碱性土壤背景验收 |
| 长武站气候与土壤背景 | 年均温约 9.2 °C，年降水约 572 mm，容重约 1.36 g cm⁻³，pH 约 8.4 | 固定场地背景交叉检查；不覆盖项目已有逐层实测表 |

当前 `organic_initial_pools.csv` 顶层 MB/SOC 比例约 0.0362；与
`soil_initial_profile.csv` 的 SOC 联用得到的 MBC 数量级落在长武长期试验报道范围内。
深层比例仍是初始化先验，必须用剖面观测或多年度预热稳定性检验，不应宣称为实测。

## 3. 源码固定常数

| 配置项 | 当前值 | 依据与处理 |
|---|---:|---|
| `respiration_fraction` | 0.55 | 本地 CoLM `MOD_Initialize.F90` 的 `rf_decomp`；固定结构常数 |
| `q10`, `frozen_q10` | 1.5, 1.5 | 本地 CoLM 初始化源码；固定环境响应形式 |
| `smpmax_mm`, `smpmin_mm` | −200, −200000 mm | 本地 CoLM 初始化源码；固定水分标量边界 |
| `depth_efolding_m` | 10 m | 本地 CoLM BGC 分解级联源码；固定深度标量 |
| `p_ep`, `p_em`, `r_ep`, `r_em`, `Q_max`, `dd_beta` | 0.001, 0.001, 0.001, 0.001, 1.2, 1 | 本地 MEMC `data-raw/memc_params.csv`；保持源方程参数合同 |

这些值不是“长武率定结果”。固定的理由是保持所耦合源码方程及其结构参数；若以后
开展结构敏感性试验，应另立方案，不能与场地率定混在一起。

## 4. 待率定参数及当前先验

AOM 和当前长武的 MEND 生理/动力学参数保存在 `site.toml`。DOM 吸附/弥散和
BIOPARTICLE 介质参数可随土层变化，保存在 `soil.csv`。
现行 source-hourly 不使用旧 MEMC 动力学字段。以下来源表不是已率定参数表。

| 组 | 字段 | 当前先验 | 主要约束数据 |
|---|---|---:|---|
| AOM | `k_litter_1d`, `k_straw_1d`, `k_manure_1d` | 0.050676, 0.013418, 0.070000 d⁻¹ | 凋落物/秸秆/有机肥质量损失与 CO₂ 通量 |
| AOM | `doc_leach_1d` | 0.0001 d⁻¹ | DOC 脉冲与 AOM→DOC 通量 |
| DOM（逐层） | `kd_doc_cm3_g`, `kd_don_cm3_g` | DOC: 0.5（20 cm）/3.98（≥40 cm）；DON: 3.98 cm³ g⁻¹ | 吸附试验或 DOC/DON 穿透曲线；2026-10-12 核查：表层 0.5 与 Camino-Serrano 2018 经验式（Kaiser 1996 数据）对长武黏粒/pH 的预测 0.512 数值巧合一致，深层 3.98 位于该数据库包络 0–4 cm³/g 上缘；DON 无公开线性 Kd 数值锚点。仍为待约束先验，不限制后续逐层率定；核查见 `docs/soil/biogeochemistry/KD_UNITS_LITERATURE_PARTITION_AUDIT_20261012.md` |
| DOM（逐层） | `dom_disp_cm` | 4 cm | 长武原状土柱氯离子平均弥散度的区域代理；DOC/DON 仍需敏感性分析 |
| MEMC | `vp_1d`–`mr_1d` | 当前 MEMC 源表/工程先验；`mr_1d=0.0052202` d⁻¹ | 时钟固定为 day/d⁻¹；POC/MOC/DOC/MBC、酶与 CO₂ 多目标率定 |
| MEND | `dormant_maintenance_ratio` | 0.001 | MEND 固定 β；无来源的独立休眠死亡支路已从代码和参数接口移除 |
| MEND 微生物 DIN 摄取 | `mend_vnup_mb_h`、`mend_ks_nh4_mb_mg_g`、`mend_ks_no3_mb_mg_g` | 0.001、0.001、0.001 | 上游 `MEND_namelist.nml` 的初始 `Pinitial(28:30)`；长武待率定先验，**不是**同文件最后覆盖的优化值 0.1、0.00018、0.00041。当前两个 Ks 在 COMPASS 中按 mg N/g 干土使用；论文与 Fortran 参数注释写 mg N/cm³，而 Fortran 计算直接与标为 mg N/g 的库存相除，因此不能声称数值相同就完成了物理单位对齐。 |
| BIOPARTICLE（逐层） | `attachment_1d`, `detachment_1d`, `straining_1d`, `microbe_disp_cm` | 0.0024, 0.0024, 0.0001 d⁻¹, 4 cm | PFLOTRAN 方程结构下的非零工程先验；当前各层同值不限制后续逐层率定 |

独立 `DON→NH₄⁺` Michaelis–Menten 支路未进入当前模型，因此原
`don_vmax_ugN_cm3_d` 和 `don_km_mgN_L` 占位字段已从配置合同删除。

## 5. 文献

1. 臧逸飞等. 26年长期施肥对土壤微生物量碳、氮及土壤呼吸的影响.
   *生态学报* 35, 1445–1451 (2015). DOI: 10.5846/stxb201305070967.
2. 长武长期氮肥梯度试验的 SOC、DOC 和 MBC 数据：
   DOI: 10.11674/zwyf.2023388.
3. Long-term fertilization and straw incorporation experiment at Changwu.
   *Frontiers in Microbiology* (2022). DOI: 10.3389/fmicb.2022.926506.
4. Gao, Y. et al. Effects of fertilization on wheat yield and SOC accumulation
   in rainfed loessial tablelands. DOI: 10.11674/zwyf.2009.0613.
5. 长武冬小麦根系研究场地说明（年均温、降水、质地、容重、pH、OM、TN）：
   University of Western Australia research repository accepted manuscript,
   `AAM_Wheat_cultivars_with_small_root_length.pdf`.
6. Kaiser, K., Guggenberger, G. & Zech, W. Sorption of DOM and DOM fractions
   to forest soils. *Geoderma* 74, 281–303 (1996).
   DOI: 10.1016/S0016-7061(96)00071-7（线性初始质量等温线 m 的源数据库）.
7. Camino-Serrano, M. et al. ORCHIDEE-SOM: modeling soil organic carbon (SOC)
   and dissolved organic carbon (DOC) dynamics along vertical soil profiles
   in Europe. *Geosci. Model Dev.* 11, 937–957 (2018).
   DOI: 10.5194/gmd-11-937-2018（Eq. 15 的 KD–黏粒–pH 经验式；KD 包络 0–4 cm³/g）.
8. Vandenbruwane, J., De Neve, S., Qualls, R. G., Sleutel, S. & Hofman, G.
   Comparison of different isotherm models for dissolved organic carbon (DOC)
   and nitrogen (DON) sorption to mineral soil. *Geoderma* 139, 144–153 (2007).
   DOI: 10.1016/j.geoderma.2007.01.012（线性 IM 等温线在深层不饱和层位失效）.
9. Kothawala, D. N., Roehm, C., Blodau, C. & Moore, T. R. Selective adsorption
   of dissolved organic matter to mineral soils. *Geoderma* 189–190, 334–342
   (2012). DOI: 10.1016/j.geoderma.2012.07.001（Langmuir Qmax/b 证据）.
10. Lilienfein, J., Qualls, R. G., Uselman, S. M. & Bridgham, S. D. Adsorption
    of dissolved organic carbon and nitrogen in soils of a weathering
    chronosequence. *Soil Sci. Soc. Am. J.* 68, 292–305 (2004).
    DOI: 10.2136/sssaj2004.2920.
11. Li, W. D. et al. 四川盆地农田土壤颗粒态和矿物结合态有机碳特征及影响因素.
    *环境科学* 46, 388–398 (2025). PMID: 39721644（表层 POC/SOC 32.5%、
    MAOC/SOC 67.5%，n=871；仅作分库先验量级参照，不是长武数据）.

文献值用于固定场地背景或给率定设置观测边界；缺少长武直接测量的过程参数仍明确
标为“待率定”，不因相邻地区或其他土类的数值而升级为固定参数。

## 6. 率定准入

当前观测覆盖、参数可辨识性、分阶段率定顺序和正式预热准入门见
`docs/soil/biogeochemistry/CHANGWU_BGC_CALIBRATION_READINESS.md`。
DOC/DON、MEMC/MEND 和 BIOPARTICLE 的逐参数来源核查、郭胜利论文可用边界及
首次正式预热前的来源归一化清单见
`docs/soil/biogeochemistry/CHANGWU_BGC_PARAMETER_SOURCE_AUDIT.md`。
在该准入门满足前，当前数值均保持“固定结构常数”或“待率定先验”身份。

## 7. 本次补充的NM处理微生物约束（0–20 cm）

根据用户提供的长期施肥试验图表，NM（氮肥+有机肥）处理的均值±标准差为：

| 观测量 | NM值 | 模型率定用途 |
|---|---:|---|
| 微生物量碳 MBC | 574.58±24.52 mg C/kg | 0–20 cm活跃+休眠微生物碳总量约束 |
| 微生物量氮 MBN | 107.52±2.08 mg N/kg | 微生物C:N与有机氮闭合约束 |
| 微生物熵（代谢商） | 5.16±0.62 % | 微生物呼吸相对于MBC的辅助约束 |
| 微生物量C:N | 5.34±0.13 | MB-C/MB-N化学计量约束 |
| 土壤基础呼吸 | 约26 mg CO₂-C kg⁻¹ h⁻¹（图1读图值） | 0–20 cm呼吸通量数量级辅助约束；保留读图不确定度 |
| 微生物代谢商 qCO₂ | 约0.95（图2读图值） | 与呼吸和MBC联合的辅助约束；保留读图不确定度 |

表中MBC/MBN为直接表格值；基础呼吸和qCO₂为图1/图2读图近似值，不作为高精度观测输入。模型输出须先转换到0–20 cm干土质量基准后比较。MBC不能直接等同于活跃库BA，应比较BA+BD；BA/BD分配仍由MEND状态方程决定。土壤呼吸比较时必须统一时间单位及C/CO₂质量基准。

## 8. 2026-09-14：源小时 C/N 生产模式

Changwu 已显式选择 `microbe.reaction_method = "source_hourly"`。
动态吸附参数使用对齐参考 namelist 的最后一次 Pinitial 赋值：Qmax=1.5 mg/g、
Kba=6 g/mg、Kdes_ref=0.006 mg/g/h、fpENZN=1；iHR=0。
MOD_MEND.F90 的 2008–2027、2329–2332 行提供吸附/解吸温度响应：
Kads=Kba*Kdes_ref，活化能分别为 5 和 20 kJ/mol。
这些数值是源码对比参数，不构成新增观测率定证据。

该模式中土层 Kd 仅用于初始化 DOC/DON 分相，后续由源小时更新推进动态吸附。
已有 reaction_substep_s 仅作用于旧 RK4 路径。新增胞内氮酶初始 C/N 均为零，
并独立保存、输运及计入有机 C/N 总量。旧模式检查点不自动迁移。
根优先、游离 DIN 动力学与剩余总量上限、未解析 NO2、供氧限制及空间运移仍为
明确的 COMPASS 耦合差异，详见
[生产接口契约](../../../docs/soil/biogeochemistry/MEND_SOURCE_PRODUCTION_CONTRACT.md)。

## 9. 2026-09-30：冬小麦地上过程修订

Editor: Wenjie Liu; Date: 2026-09-30。

- `crop_harvest_rule="observed_lai_zero"` 按当前输入 LAI 的正值→精确零值边沿
  触发一次收获；当前输入对应 2001-06-09。这是用户选定的管理代理，并非从
  LAI 独立验证了田间实际收获日期。长武长期试验原始论文记载翌年6月底收获
  （DOI [10.1080/09064710.2016.1149215](https://doi.org/10.1080/09064710.2016.1149215)）。
  观测与论文日期差异保留待核；生理成熟的热量条件单独保留。
- `luna_leaf_n_source="plant_pool"` 使用实际绿色叶片 N / 当前绿色 LAI，
  随氮再转运和衰老更新功能氮上限；LUNA 来源为 Ali et al. (2016)，
  DOI [10.5194/gmd-9-587-2016](https://doi.org/10.5194/gmd-9-587-2016)。
  全球氮分配机制来源不等于长武品种的光合参数率定。
- `crop_wsc_alloc_fraction=0.30` 与 `crop_wsc_retrans_daily_fraction=0.50`
  是有明确原始源码位置的 APSIM 参考默认值，按 COMPASS 小时步长与实际
  C/N 库存适配；不是冬小麦实测转运比例。年度 WSC 机制由 Takahashi et al.
  (2001) 冬小麦实验支持（DOI [10.1626/pps.4.160](https://doi.org/10.1626/pps.4.160)）。
  默认衰老组织 C:N 复用已编译 `CN_leaf_grainfill`，明确属于待验证的代理。
  源码固定版本、方程、限制和库存闭合见
  [年度叶片/WSC说明](../../../docs/plant/winter_wheat_annual_physiology.md)。
- NEE 使用最终接受的 PlantCN `cum_GPP/cum_R_a`、AOM 累计 CO₂ 和
  MEND `IC` 差量；无新增经验呼吸率参数。活根属于植物自养呼吸；表面
  AOM 直接呼吸与 MEND 土内呼吸分别消费，凋落物内转移不重复计为呼吸。
- 叶片光学 `Cdm/Cw/Nmesophyll`、叶倾角、SLA 和编译参数消费者的证据等级
  与后续校准边界，见 `docs/plant/winter_wheat_implementation_review_20260930.md`。
- 按用户提供的本地 CoLM 冬小麦参数表及同义消费者核实后，更正 PFT25 默认
  `arooti=.10`、`fleafi=.90`、`bfact=.10`、`grnfill=.67`、`mxmat=270 d`。
  当前实际 forcing 的灌浆事件由旧 2001-04-05 移到 2001-05-04；生理成熟和
  LAI-zero 收获日期保持 2001-05-25/06-09。原版 CoLM 的热量成熟即收获规则
  未覆盖用户选定的 LAI 收获策略。源文件、精确值/行号和消费者见
  [叶氮与参数审计](../../../docs/plant/winter_wheat_leaf_n_and_parameter_audit_20260930.md)。

原始观测、当前冻结年度基准及历史来源记录不因上述修订而被重写。


## 10. 2026-10-06：用户采用本地最终水力参数

Editor: Wenjie Liu; Date: 2026-10-06。

用户明确决定：“这组的参数可以作为最终参数了”。采用刚展示的本地
`foundation_alpha_rebalance` 完整参数表；全年生产结果已完成 8760 个接受步，
全部 75 对水分观测及预测来源已核验。用户采用后停止后续搜索，两项未完成的
200 cm 候选保留日志与检查点，不参与评分；六层逐行搜索没有全部完成。

采用档案：`outputs/local_water_rebalance_20261006/accepted_final_20261006_192026/`。
原始选定表 `soil.csv` 的 SHA256 为
`79f74b4965f0aa8057daae9fecbade8531aa801e105cb034a023e3fcb0396b36`。
本地活动表仅写入其 `ks_cmd`、`vg_alpha_1cm`、`vg_n`、`theta_r`、`theta_s`
五列；非水力参数、原始注释、编码及换行保留。修改前完整表及逐项应用凭证见
`outputs/local_water_rebalance_20261006/acceptance_control_20261006_192026/`。

小于 200 cm 的水力参数采用 2026-10-06 10:16 的服务器最优候选 20005
快照（其中 20 cm 的 `vg_n=1.15`）；这只是该时刻的服务器最优组。
200/220/240/260/280/300 cm 采用本地选定组，500/700/1000 cm 固定。
Ks 的表单位是 cm/d，α 是 1/cm，n、θr、θs 为无量纲；m 由 n 自动重建。

| 深度 (cm) | Ks (cm/d) | α (1/cm) | n | θr | θs |
|---|---:|---:|---:|---:|---:|
| 200 | 16.49114072 | 0.0065 | 1.597984717 | 0.1 | 0.407 |
| 220 | 12.62124411 | 0.0063 | 1.52300289 | 0.095 | 0.4074 |
| 240 | 8 | 0.0082 | 1.453357967 | 0.09 | 0.4078 |
| 260 | 8 | 0.01114285714285714 | 1.458112163 | 0.085 | 0.4082 |
| 280 | 10 | 0.0145 | 1.415218849 | 0.08 | 0.4086 |
| 300 | 9 | 0.0145 | 1.466234283 | 0.075 | 0.409 |

率定仅使用 2000 年 8/10/12 月三期；200–300 cm 合并 RMSE 为
0.0101174622862 m³/m³，最大绝对误差为 0.0199561709881 m³/m³。
各目标深度的 RMSE 和每期绝对误差均低于用户给定的 0.03 m³/m³。
2001 年 2/4 月只作保留期诊断，不参与选择：目标深度合并 RMSE 为
0.0144724502940 m³/m³，最大绝对误差为 0.0308399750614 m³/m³。
因此保留期不能表述为全部通过 0.03 阈值，也不升级为独立正式验证。

既有冻结年度结果作为本次采用证据；本次应用参数没有新增年度模拟。
原始观测、历史参数与评分、冻结源代码和服务器文件均保留。

## 2026-10-10：植物与腐熟厩肥直接输入 MEND 的本地试跑

Editor: Wenjie Liu. 用户明确选择本次路由与厩肥碳分配先验；不代表完成率定。

- 植物由 CoLM 产生实际凋落 C/N 与时序，入土后一次分到 POM1/POM2/DOM，
  碳比例沿用本地 MEND `Fraction_type1=[0.07,0.37,0.56]`。新增植物凋落不经
  CoLM 组分 AOM 前置分解；地上输入按最上层即时接收，根沿现有深度权重。
- 腐熟厩肥有机碳采用 DOC=0.10，其余按既有 POM1:POM2=7:37 分配，故为
  `[0.1431818181818182,0.7568181818181818,0.10]`。10%参考 Morvan & Nicolardot
  (2009), doi:10.1007/s00374-009-0355-1 的腐熟牛厩肥水溶碳近似值。
  水提条件、原料成熟度与本地材料不完全对应；水溶性也不证明快速生物可利用性。
  余下碳的7:37是用户选择的模型先验，不是该论文实测值。
- 先保留原有5% TAN/总氮先验并扣除 TAN，再分配有机氮。使用源 MEND 的
  输入纤维素 C:N=500 与 `rCN_LIG2LAB=2`（本地 `SOIL_INI.dat`），以实际输入
  总C、总有机N闭合；不会把碳的10%或56%复用为DON比例。此处外源组分C:N
  关系同样属于移用假设，不是腐熟厩肥DON实测证据。
- 混施0–20 cm全量按交叠厚度均匀分层；表施直接进入最上层。取消旧80%混施、
  20%地表AOM的隐式拆分。各形态C/N使用同一空间质量权重，存储单位分别转换。
- 独立秸秆还田保留straw AOM；长武当前还田事件为空。历史AOM库存不清零。
- 新增路由、分配与C:N参数进入重启兼容性签名，旧检查点不能静默当作新过程续跑。
  新入POM的碳较旧AOM更早进入现有SOC观测汇总，须同时检查完整C/N库存与通量。


## 2026-10-12：用户采用旧四参数开展一年 DOC 初值对照

Editor: Wenjie Liu.

按用户明确指令，18个 soil.csv 输入层统一采用 Vg=0.00212597 h^-1、KsDOM=0.30 mg C/g、CUE20=0.25、alpha=0.09545998321928836；Vm仍为0.000224362722223168 h^-1。
trial24水氮及其他参数不变，继续线性平衡吸附。关闭被动氮肥示踪；施肥输入与总氮过程保留。
本地活动初始剖面保留用户CSV自由DOC及原DON，零DOC/DON仅在独立试验副本。服务器另用有明确来源的历史预热自由DOC；不替换完整碳氮状态。
本次一年试验目录：outputs/old4_trial24_annual_20261012_01/。改前文件与参数来源保存在其provenance目录。61天证据不代表全年通过。

## 2026-10-30：用户确认收获期指定根残体输入

Editor: Wenjie Liu.

用户明确选择沿用工作簿给定数量、每年6月25日输入。新增 `[root_residue]`
独立土壤边界驱动；保持 `plant_cn_cycle=false`、`litter_input=false` 及既有
Jackson根剖面不变。冻结的ID41和历史模拟不重写。

- 数量来自 `NM冬小麦实测.xlsx` 的“生物量”表：`F3=1.75 t C ha⁻¹ yr⁻¹`、
  `F6=18.52 kg N ha⁻¹ yr⁻¹`。它们是工作簿给定的静态年输入，不是同一行年份的
  年度实测通量。原工作簿SHA256为
  `5b18a4da33785ab4a61f3a4e062098296f70e14cf58610e5ad31e80e7f626b50`。
  新事件直接使用1750 kg C/ha和18.52 kg N/ha，不由舍入C:N反算N；对应C:N为
  94.4924406048。根际沉积与死根的贡献在原表中未分开，本次采用用户指定的
  收获期单次输入简化。
- 臧逸飞等2015《26年长期施肥对土壤微生物量碳、氮及土壤呼吸的影响》
  （DOI [10.5846/stxb201305070967](https://doi.org/10.5846/stxb201305070967)，
  p1447 §1.2，同含NM处理）记载次年6月下旬收获。6月25日是用户确认的代表日，
  不冒充2000/2001年逐年核实的收获日期。原有LAI-zero物候代理日期不变。
- 郭胜利、吴金水、党廷辉2008《中国农业科学》41(3):744–751，
  DOI [10.3864/j.issn.0578-1752.2008.03.015](https://doi.org/10.3864/j.issn.0578-1752.2008.03.015)，
  记载秸秆不还田，地下根系为返田残体主体；其0–20 cm根茬量按籽实产量35%
  估算。该口径不用于替换本次用户选定的工作簿数量，也不声称工作簿数量由该式
  得到。不额外加入没有数量依据的秸秆或麦茬。
- 根残体沿当前 `vegetation.root.RAI` 的归一化剖面分配。当前为Jackson
  `root_type=6`、最大深度200 cm；不更改活根/水力分布，也不读取未启用的
  `litter.root_distribution_beta` 作为第二套根分布。根有机C/N先进入既有
  根残体AOM组分，随后沿既有分解流程供给MEND；没有直接NH₄或TAN补给。
- 指定根残体是土壤控制体的外部输入，独立于本次未跟踪植物库存的被动吸氮；
  不计为PlantCN内部凋落转移，不与原肥料输入重复。未来启用PlantCN/自动凋落时
  需先取消同一来源的指定根事件，防止重复归还。
- 当前年窗为2000-06-15至2001-06-15，实际命中的根事件是2000-06-25，
  不是窗外的2001-06-25。前一季根输入不自动解释为本年383 kg N/ha吸氮的
  固定归还份额。是否改善NH₄拟合必须另经对照模拟验证，不能由接口通过代替。

原文、工作簿单元格、事件配置与检查记录见
`data/verification/changwu_harvest_residue_input_20261030_01/`。


## 2026-10-31：采用已确认的K1无机肥0–20 cm全年参数基线

Editor: Wenjie Liu.

用户在查看完整水碳氮剖面后确认当前参数，并要求归纳对齐到本地。采用的实际算例为
`mend_k1_mineral_0_20cm_annual_20261031_16 / K1_mineral_0_20cm`，不是旧release附带的ID41_ref。
该算例8760步及原生水/热/C/N/O2和矿肥、有机肥有机N、有机肥TAN三来源守恒检查通过；
归档绑定release manifest SHA为`f3c2ade97fc2eb3515d0cdcbbc506850627df2ae9c5a880bf76ecd633ac4b1ae`。

- soil.csv采用实际已材料化的有效参数，18行×37列，保持8个MEND逐层列；
  不再次乘历史MOSPO倍率。与改前本地相比8列130格更新，水力、碳参考参数等其他29列不变。
- organic_initial_pools.csv采用已验证的B_ref分配；62.5–287.5 cm的22个计算层
  仅MON/DON/QON初始N三列共66格更新。全部初始C分配不变，各层总N守恒至双精度舍入。
- FAST10入口、矿质肥来源标签及有机肥有机N/TAN来源标签启用。
  有机肥10%有机C与纯有机N快入DOC/DON；TAN按原5%单列、不重复计算。
- 无机肥与有机肥均为0–20 cm；原9月15日10:00及120 kg N/ha化肥量保留。
- 本轮采用无指定根残体的完整全年方案，取消活动[root_residue]事件。
  上一节2026-10-30根残体数量、日期、来源及原授权仍为历史记录，不在本轮生效。
  原完整事件与站点24个文件保存在本轮before/site；未来恢复应作为另一套输入单独验证。
- Qmax=1.5、Kba=6、Kdes=0.006仍在TOML。前期讨论的三参数迁到soil接口未实施，
  不把文件位置调整当成本次数值采用的一部分。本地新版schema保留其可选接口能力。
- 本地src817个文件（含806个.m）与本算例r4逐字一致，本次不覆盖或改动生产源代码。
  气象、观测、初始水热/DIN、网格、求解器和本地相对路径保留；不引入服务器绝对路径。

原全年性能：水RMSE0.0241069241 m3/m3、NO3 RMSE6.8926827122 mgN/kg、
NH4 RMSE1.7668810635 mgN/kg；SOC2019参照RMSE0.9976968386 gC/kg，
DOC2006参照RMSE24.1168624295 mgC/kg。五期水氮此前参与筛选，碳为跨年参照，
不将这些指标改称独立验证。数值被用户采用，不等于全部参数已直接实测。

差异、原始SHA、改前副本、参数目录与本地验证记录集中在
`data/verification/mend_k1_accepted_local_alignment_20261031_21/`。
本轮不恢复服务器率定，也不启动新的全年模型计算。


## 2026-10-31：Qmax、Kba、Kdes迁为soil逐层单一所有权

Editor: Wenjie Liu.

在已采用的 `K1_mineral_0_20cm` 数值基线上，按用户确认把
`mend_qmax_mg_g`、`mend_kba_g_mg`、`mend_kdes_mg_g_h` 从 `[microbe]` 移到
`soil.csv`。18行分别保留1.5 mg C/g干土、6 g干土/mg C、0.006 mg C/g干土/h。
三键不再在TOML重复定义；旧版TOML标量站点仍由兼容接口支持。

当前soil18×40，含11个MEND列；microbe35键、TOML合计399叶键。
这次只变更所有权与向量接口，不改变初始有机C/N、管理事件或参数数值；
上文8/38、Qmax仍在TOML及r4源码逐字一致的描述均为迁移前历史记录。
新接口源码不再以迁移前r4源码哈希作为当前源码一致性的声明。

Kdes_ref既参与 `Kads=Kba*Kdes_ref*Arr(T,5)`，也构造
`Kdes=Kdes_ref*Arr(T,20)`，其单位为mg C/g干土/h，不能简写为h^-1。
三参数最初数值沿用上游示例先验；文件迁移不构成长武测量或新的场地率定。

迁移后参数总表以重新读取当前输入的目录为准。原8760步指标保留为历史已采用结果；
本轮接口短步回归结果单列，以实际验证回执为准，不冒称重新全年验证。
