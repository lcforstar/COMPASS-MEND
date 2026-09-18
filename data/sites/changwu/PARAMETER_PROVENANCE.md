# 长武站土壤碳氮参数来源与率定边界

## 1. 配置归属

`soil.csv` 按以下固定顺序组织：土壤基本属性 → 水力参数 → 无机氮参数 →
有机碳/AOM 参数 → DOC/DON 参数 → MEMC 参数 → MEND 参数 →
BIOPARTICLE 参数。第二行是以 `#` 开头的中文含义和单位说明，读取器会将其
作为注释忽略。

`site.toml` 只保留管理事件、数值步长、初始化选项和有明确源码依据的结构常数。
所有动力学强度、吸附、休眠和运移参数均在 `soil.csv` 中率定。当前实现要求这些
率定列在各深度重复同一个初值；这是“场地全局参数”约束，不是假定已有分层证据。

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

AOM 与 MEMC/MEND 参数在当前方程中不设置独立土层自由度，保存在
`site.toml`。DOM 吸附/弥散和 BIOPARTICLE 介质参数可随土层变化，保存在
`soil.csv`。

| 组 | 字段 | 当前先验 | 主要约束数据 |
|---|---|---:|---|
| AOM | `k_litter_1d`, `k_straw_1d`, `k_manure_1d` | 0.050676, 0.013418, 0.070000 d⁻¹ | 凋落物/秸秆/有机肥质量损失与 CO₂ 通量 |
| AOM | `doc_leach_1d` | 0.0001 d⁻¹ | DOC 脉冲与 AOM→DOC 通量 |
| DOM（逐层） | `kd_doc_cm3_g`, `kd_don_cm3_g` | 3.98, 3.98 cm³ g⁻¹ | 吸附试验或 DOC/DON 穿透曲线；当前各层同值仅为初始先验，不限制后续逐层率定 |
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
