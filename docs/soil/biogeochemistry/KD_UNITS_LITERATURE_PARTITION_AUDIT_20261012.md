# Kd 单位、文献取值范围与初始分库方式核查：2026-10-12

Editor: Wenjie Liu  
Date: 2026-10-12

## 范围

只读核查与文档补充；模型推进步数为 0，未修改任何方程、配置、初态、观测或历史结果。
本核查承接 `outputs/doc_partition_audit_20261012/AUDIT.md` 的下一步决定与
`outputs/doc_free_initial_profiles_20261012/QA_AND_FINDINGS.md` 提出的
"Kd、水/土浓度单位与供体库政策是否共同代表目标实测库"问题，按用户指定优先级
核查三件事：**kd 的单位、文献取值范围、初始分库方式**。

## 1. kd 单位核查

结论：DOC/DON 平衡吸附 Kd 的单位在全链路为 **cm³/g（数值等于 L/kg）**，方程量纲
闭合，未发现单位混用、1000 倍转换或剖面颠倒。

| 位置 | 检查内容 | 结论 |
|---|---|---|
| `soil.csv` 表头 | `kd_doc_cm3_g`、`kd_don_cm3_g`（"DOC/DON平衡吸附Kd(cm3/g)"） | cm³/g ✓ |
| [tableSpec.m:68](src/+io/+adapters/tableSpec.m#L68) | 单位元数据 `'kd_doc_cm3_g','cm3/g','kd_don_cm3_g','cm3/g'` | ✓ |
| [buildSoilParameters.m:85](src/+init/buildSoilParameters.m#L85) | soil.csv → `dom_kd_doc_cm3_g` 直读 | ✓ |
| [checkSite.m:318](src/+io/checkSite.m#L318) | 要求逐层 >0 | ✓ |
| [equilibrate.m:28](src/+soil/+dom/equilibrate.m#L28) | S=(θ+ρ·Kd)·Caq；q_sorbed=Kd·Caq；ρ·Kd 无量纲 | ✓ 量纲闭合 |
| [phaseState.m:31](src/+soil/+dom/phaseState.m#L31) | free=θ·Caq；sorbed_bulk=ρ·q | ✓ |
| [buildSoilDOM.m:81](src/+init/buildSoilDOM.m#L81) | mg/kg→g/cm³ 乘 ρ×1e-6；mg/L 旧口径乘 (θ+ρKd)×1e-6 | ✓ |
| [totalFromFree.m:21](src/+soil/+dom/totalFromFree.m#L21) | T=D·(1+ρKd/θ)，即 Q=D·ρKd/θ，不凭空造质量 | ✓ |
| [writeStep.m:129](src/+diagnostics/writeStep.m#L129) | D=θ·Caq、Q=ρ·Kd·Caq，D+Q=T 闭合断言 ≤1e-12 | ✓ |
| 审计回算公式 | D_mgkg=Caq_mgL·θ/ρ；Q_mgkg=Caq_mgL·Kd | ✓ 量纲验证 |

换算陷阱记录（防止未来搬运文献数值时出错）：

- L/kg = cm³/g，数值完全相同；
- m³/kg 的数值是 cm³/g 数值的 1/1000（如 0.001226 m³/kg = 1.226 cm³/g）；
- 文献 mg/L 是水相浓度，mg/kg 是干土质量基准，二者按 (θ+ρKd)/ρ 换算，
  代码已在 `buildSoilDOM` 拒绝同填两个口径；
- NH4 吸附系数 `kg_cm3g`（当前 0.668/1.159 cm³/g）是独立参数，不在本次 DOM Kd 范围。

## 2. kd 文献取值范围

现有 0.5/3.98 在来源表中标为"待约束先验"，此前未记录文献范围。本次核查结果：

| 来源 | 内容 | 数值 |
|---|---|---|
| Kaiser et al. 1996, Geoderma 74:281–303, DOI 10.1016/S0016-7061(96)00071-7 | 125 个森林土壤样品；线性初始质量等温线 RE=m·Xi−b；m 由 Alo/Fed/OC 控制（无碳酸盐土壤解释 75–77% 方差，含碳酸盐土壤 Fed+OC 解释 58–67%） | m 数值范围未在公开摘要给出 |
| Camino-Serrano et al. 2018, GMD 11:937–957, DOI 10.5194/gmd-11-937-2018 | 用 Nodvin et al. 1986 关系 KD=m/(1−m)·(V/M) 把 Kaiser 1996 的 m 换算为 KD；Eq. 15: **KD = 0.001226 − 0.000212·pH + 0.00374·clay**（clay 为 0–1 小数）[m³/kg]，adjusted R²=0.25；34 个欧洲剖面 | 数据库 KD 包络 **0–0.004 m³/kg（0–4 cm³/g）**（补充材料 Fig S1 轴范围，已核实原文与补充材料） |
| Vandenbruwane et al. 2007, Geoderma 139:144–153, DOI 10.1016/j.geoderma.2007.01.012 | 同一批土壤上 DOC 与 DON 等温线比较；线性 IM 等温线只适合低容量或近饱和层位，深层不饱和层位为曲线、修正 Langmuir 最优 | 线性常数 Kd 用于深层是最弱假设（定性证据；数值表未公开） |
| Kothawala et al. 2012, Geoderma 189–190:334–342, DOI 10.1016/j.geoderma.2012.07.001 | 43 个加拿大矿物层位，Langmuir Qmax/b | 吸附容量证据，不是线性 Kd |
| Lilienfein et al. 2004, SSSAJ 68:292–305, DOI 10.2136/sssaj2004.2920 | 年代序列上 DOC+DON 吸附，Langmuir + null-point，吸附随 allophane 增加 | 同上；DOC/DON 同批实验证据 |

长武代入（Eq. 15）：

- 表土黏粒 30.1%、pH 8.68 → KD=0.000512 m³/kg = **0.512 cm³/g**，与当前表层
  kd_doc=0.5 几乎重合（差 2.4%）；
- 40 cm 以下（pH 8.77、黏粒 30.9%）→ 预测约 0.52 cm³/g；当前深层 3.98 cm³/g
  位于 Kaiser 数据库包络的上边缘（0.004 m³/kg ≈ 4 cm³/g），没有 pedotransfer 支持；
- Eq. 15 的拟合域是酸性/中性欧洲森林土壤；长武是 pH 8.4–8.7 的碳酸盐黄土，
  超出拟合域（低黏粒时该式在 pH>5.8 已给负 KD）。表层 0.5 与预测重合只能记录为
  "数值巧合一致"，不能作为验证；
- DON：未取得公开的线性 Kd 数值锚点（Vandenbruwane 2007 数值表在付费墙内）；
  文献普遍认为 DON 吸附弱于或不同于 DOC。当前 kd_don=3.98（与深层 DOC 同值）
  保持"待约束先验"身份。

结论：**0–4 cm³/g 是当前可核实的文献包络**（Kaiser 1996 数据库经 Camino-Serrano 2018
换算）；0.5 在包络内且与 pedotransfer 预测一致；3.98 在包络上缘但无独立依据。
按 `tools/remote/doc_proxy_calibration/BOUNDS_AND_DOC_CONTRACT.md` 的原则，该包络
不得与动力学吸附 Kba/Qmax/Kdes 的范围混用，也不自动成为搜索边界。

## 3. 初始分库方式核查

原生冷启动链路：

1. `soil_initial_profile.csv` 实测 SOC（soc_total_mgkg）→ g/cm³
   （[buildSoilCarbon.m:181](src/+init/buildSoilCarbon.m#L181)）；
2. `organic_initial_pools.csv` 逐层先验比例（frac_POC/MOC/MB/IOM，每行和为 1）
   → pocC/mocC/mbC/iomC；qoc/ep/em 初始为 0
   （[buildSoilCarbon.m:33](src/+init/buildSoilCarbon.m#L33)）；
3. POM 按 `mend_lci_initial` 分 POM1/POM2；有机氮 = 碳/目标 C:N，MB 按
   `active_fraction_initial` 分活跃/休眠；
4. [buildSoilDOM.m](src/+init/buildSoilDOM.m) 读 DOC/DON 输入；新口径
   `doc_free_mgkg`=D 经 [totalFromFree.m](src/+soil/+dom/totalFromFree.m) 重建 D+Q；
5. 通用路径 [reserveInitialDOM.m](src/+init/reserveInitialDOM.m) 从 MOC 扣 D+Q、
   从 MON 扣 DON，SOC 守恒；
6. 长武代理驱动 `applyChangwuProxyDomInitialState`（doc_free_d_v2）不从 MOC 单扣，
   而是按原比例缩放全部响应库（POM1/2、MOC、MB、BD、ENZP1/2、EM、ENZN、NEC），
   配对有机氮同步缩放，SOC 逐层守恒；实测 DON 总量直接代入（口径 D+Q，
   已在 audit struct 显式声明 `don_input_basis='legacy_total_DON_plus_QON...'`）。

核查结果：

- 比例文件自标"先验，不是逐层实测值"，表头与代码口径一致 ✓；
- 表层有效 MB/SOC≈0.0362，与长武 MBC 文献范围 254.65–745.26 mg/kg 兼容
  （PARAMETER_PROVENANCE §2）✓；
- 表层 MOC 先验 0.72 与农田表层 MAOC/SOC 文献量级一致：四川盆地 871 个表层样
  POC/SOC 32.5%、MAOC/SOC 67.5%（Li et al. 2025, 环境科学 46:388–398,
  PMID 39721644）。但模型 MOC 是动力学库（新定义 MAOC=MOC+Q），不等于 <53 µm
  实测 MAOC，口径差异必须声明，不能直接对标；
- 深层比例（40–500 cm 的 POC 0.05–0.10、MB 0.005–0.02）无逐层实测或文献锚点；
- 全剖面 IOM=0 意味着把全部 SOC 当作响应库，是无惰性库数据的先验选择，
  未做过敏感性检验；
- 两种供体政策并存（通用 MOC-only 与代理驱动全响应库比例缩放）已记录，
  不是代码缺陷；新 doc_free_d_v2 筛选使用后者（最低响应库保留比例约 0.4183）。

## 4. 结论与遗留

| 项目 | 状态 |
|---|---|
| kd 单位 | 全链路 cm³/g 一致，量纲闭合，无换算错误；m³/kg↔cm³/g 陷阱已记录 |
| 文献范围 | DOC 线性 Kd 可核实包络 0–4 cm³/g（Camino-Serrano 2018/Kaiser 1996）；表层 0.5 与 pedotransfer 预测一致；深层 3.98 在包络上缘、无独立依据；DON Kd 无公开数值锚点 |
| 初始分库 | 方法正确（实测 SOC × 先验比例，D/Q 从响应库预留，SOC 守恒，闭合断言存在）；除表层 MB（长武 MBC 约束）与 MOC 量级（农田 MAOC 文献）外，比例均为先验；IOM=0 与深层比例无锚点 |

遗留（按证据需求排序）：吸附试验或 DOC/DON 穿透曲线仍是确定长武 Kd 的正途；
深层 Kd=3.98 需要解释或替换依据；DON Kd 缺文献；IOM=0 与深层分库比例敏感性未做。
辨识性扰动（1/0.3/0.1/0.03 倍 Kd 短窗试验）仅作为分配机制辨识工具，不作为
无先验自由率定，见 `outputs/doc_partition_audit_20261012/AUDIT.md` 下一步决定。

本核查不启动任何运行、不改变任何参数值；文档来源均可核实，未核实的部分已显式标注。
