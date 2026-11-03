# 六年 DOC/DON 与表层 SOC：参数候选与初始化审计

Editor: Wenjie Liu; Date: 2026-09-22.

## 1. 范围与证据归属

本文件记录对既有结果的只读审计，不代表候选参数已经完成六年率定。此次审计没有修改原站点、观测、历史输出或服务器文件。下面的参数方向以当前生产源码为机制依据，历史数值以对应冻结运行的接受态为依据，二者不混写成同版本参数敏感性结论。

- 已完成的六年运行：`outputs/mend_multiyear_20260917/vg050_2000_2006`。时间为 2000-06-15 00:00 至 2006-09-15 00:00，共 **54,792 小时**，不是整六个日历年。`result.json` 为 `passed=true`；全期最大逐小时有机 C/N 残差分别为 `2.18353108285955e-15`、`2.25836200209977e-16 g/cm²`。
- 运行定义：`outputs/mend_multiyear_20260917/runMultiYear.m`、`protocol.md`、`start.ps1`。历史源码位置由该启动脚本明确指定为 `docs/model_comparison/annual_runs/corrected-manure-tan05-annual-20260916/src`；不能用当前 `src` 的哈希替代它。
- 独立观测配对：`C:/Users/Administrator/Desktop/COMPASS_Report_Code_20260918/matched_observations.csv`、`metrics_by_period.csv`；生成方法见同目录 `prepare_figures.py` 及项目内 `outputs/mend_multiyear_20260917/prepare_figures.py`。
- 当前观测工作簿：`data/sites/changwu/NM冬小麦实测.xlsx`。本次重算 SHA256 为 `2c9392b50057de25bdb02ba2ded518cdc021b1e53a7f6038bf3ab084e0e6f9bd`，与上述观测配对的 `source_manifest.json` 一致。
- 库存/通量辅助核查：`outputs/mend_multiyear_budget_20260918/阶段核查.md`、`mass_balance_check.json`、`annual_layer_concentrations.csv`。其中部分重放只覆盖有限窗口；六年缺项差额不能假称独立测得的输运通量。

## 2. 观测算子及实测偏差

沿用已确认的观测处理，不为了接近实测而修改输出定义：

1. DOC/DON 比较的是**模型水相库存按干土质量折算**的 mg C 或 N/kg，即 `free_g_cm3_bulk / bulk_density_gcm3 * 1e6`；不是水溶液 mg/L，也不是包含吸附相的总 DOM。没有引入额外浸提比例。
2. 既有 DOC/DON、SOC 剖面沿模型单元中心深度线性插值到工作簿深度标签；0–20 cm SOC 时间序列另按各单元与 0–20 cm 的重叠厚度乘容重加权。不要把 20 cm 点值与 0–20 cm 平均值混作同一观测算子。
3. SOC 使用 `soc_observation_profile`，排除未分解 AOM。含 AOM 总有机 C 应另列，不能用它替换原 SOC 评分。
4. 匹配的是观测日期 00:00 的接受态；从开始时刻计算小时数后使用相应结束态记录。SOC 初始 2000-06-15 值属于初始化约束，不是独立验证点。

### 2.1 2006-09-15 DOC/DON

| 变量与范围 | 模拟，mg/kg | 实测，mg/kg | 说明 |
|---|---:|---:|---|
| DOC，20 cm 点 | 14.690616 | 243.800000 | 低约 94.0% |
| DON，20 cm 点 | 0.885439 | 12.389510 | 低约 92.9% |
| DOC，15 个观测深度的算术均值 | 3.874732 | 59.631333 | 不是土柱库存平均 |
| DON，15 个观测深度的算术均值 | 0.175150 | 4.893217 | 不是土柱库存平均 |

DOC 和 DON 剖面 RMSE 分别为 `78.844442` 和 `5.101291 mg/kg`。其剖面相关性虽高，量级仍明显偏低；相关性不能替代偏差和质量库存核查。

### 2.2 六个年度 0–20 cm SOC 端点

| 日期 | 模拟，g C/kg | 实测，g C/kg | 实测减模拟，g C/kg |
|---|---:|---:|---:|
| 2001-06-15 | 9.325180 | 9.830000 | 0.504820 |
| 2002-06-15 | 9.238305 | 10.050896 | 0.812591 |
| 2003-06-15 | 9.429696 | 10.072300 | 0.642604 |
| 2004-06-15 | 9.635961 | 10.317513 | 0.681552 |
| 2005-06-15 | 9.897894 | 10.360000 | 0.462106 |
| 2006-06-15 | 10.184058 | 10.650000 | 0.465942 |

六点 RMSE 为 `0.608720 g/kg`，平均模拟偏差为 `-0.594936 g/kg`。保存容重对应 0–20 cm 干土质量 `2.6e6 kg/ha`；最后一个差值相当于 `1211.449 kg C/ha` 的库存差，**不是允许额外施入的碳量**。2006-09-15 模拟 SOC 为 `10.479179 g/kg`，没有该日同算子的 SOC 实测配对，不能拿它与 6 月 15 日观测比较来缩小误差。

### 2.3 分相与总量的上限核查

本次直接读取 `continuous/nc/soil_organic.nc` 的末时刻各库，并用 `initial_state.mat` 容重与 `soil_water.nc:depth_node` 重算：

| 2006-09-15，0–20 cm 干土加权 | 游离相，mg/kg | 吸附相，mg/kg | 总量，mg/kg |
|---|---:|---:|---:|
| DOC-C | 12.564582 | 33.020995 | 45.585577 |
| DON-N | 0.956277 | 20.007004 | 20.963281 |

同一时刻 20 cm 插值点的总 DOM-C 为 `58.971148 mg/kg`，仍小于 DOC 实测 `243.8 mg/kg`；不能仅靠对当前端点库存重新分相达到该 DOC 观测。40 cm 总 DOM-C 为 `231.206658 mg/kg`，提示深度间情况不同。DON 在 20 cm 总量为 `19.411947 mg/kg`，但 60 cm 总量 `4.741564` 已低于实测 `5.375520 mg/kg`；更深处亦有总量不足。因此分相、供给、摄取与层间迁移必须分开解释。

### 2.4 数据局限

- 本次六年窗口只有一个 DOC/DON 观测日期；15 个深度不是 15 组时间序列，不能唯一辨识多个分解、摄取、吸附和运移参数。
- 工作簿现有数值按用户先前确认的观测身份保留，不重写来源；其深度标签是否严格代表点样或厚度复合样，没有在本次重新获得更细的采样支持信息，因此继续使用原算子，不自行更换。
- 初始 DON 数值与后期部分观测存在已有对应关系，初始 DOM 分配是先验，不可将它作为同日独立验证证据。
- 2010 年播前 MBC/MBN、微生物 C:N 与呼吸资料只能作跨年辅助约束；不能以其样本标准差设本六年每小时的硬阈值。工作簿 SMB 的历史解释是微生物残体 C，不等于活跃加休眠 MBC。
- 地上残体移除和模型根碳归还沿用原管理。工作簿外部估计根碳输入 `1750 kg C/ha/year` 不得在现有 PlantCN 返回之外再当成一笔新增碳源。

## 3. 旧六年有效参数不等于当前 TOML 文本

历史 `runMultiYear.m` 直接装入旧 `initial_state.mat` 的全部 Y/parameters。以当前站点重新初始化时，需要显式重现以下有效设置，而不是复制旧 Y 跨版本继续：

| 项目 | 当前 `data/sites/changwu` 文本/常规初始化 | 旧六年有效值或操作 |
|---|---|---|
| 时间末端 | 2001-06-15 | 2006-09-15 |
| `mend_vg_h` | 0.00425194 h⁻¹ | 0.00212597 h⁻¹ |
| `mend_ks_dom_mg_g` | 0.00010034 mg C/g | 0.30 mg C/g |
| `mend_ihr` | 0 | 1 |
| `mend_alpha`/派生 Vm | alpha=0.05012233；普通初始化按 Vg 派生 | 原试验只改 Vg，Vm 保持 0.00022436272222316797 h⁻¹ |
| 源小时酶初值 | ENZP1、ENZP2、EM 均为零 | 每类为所在层总 SOC 的 5e-5，全剖面播种；C 从 MOC、N 从 MON 守恒转移 |
| 0–20 cm MBC/SOC | 0.08 | 0.06，只改变中心深度 <20 cm 的单元；释放 C/N 回 MOC/MON |
| 初始总 DOM-C | 当前 CSV 原值 | 原值不变，0 cm=245、20 cm=195 mg/kg，0–20 cm平均=220 mg/kg |
| 示踪开关 | true | 历史运行未包含新来源标签；新对照和候选都保持 true |

`init.buildSoilMicrobes` 定义 `Vm=Vg*alpha/(1-alpha)`，并从 Vm 派生死亡和休眠维持率。若通过 fresh 配置保持旧 Vg 和旧 Vm，必须使用 `alpha=Vm/(Vg+Vm)=0.09545998321928836`；不是同时将 Vm 再减半。alpha 的不同仅用于表达旧有效参数，不能冒称新的独立率定结果。新构建后应逐项核对 Vm、rMORT、VmD、酶分配率和其余参数。

### 3.1 需要原样保留的输入

本次比较 `data/sites/changwu_vg050_2000_2006_20260917` 与当前站点：

- `soil_initial_profile.csv`、`organic_initial_pools.csv`、`plant_initial_pools.csv`、`forcing.csv`、`canopy_measurements.csv`、`root.csv`、`water_boundary.csv` **逐字节相同**。
- `soil.csv` 唯一文本差为首层 DOC Kd 的 `0.50` 与 `0.5`，数值全部相同。已率定水力/无机 N、容重、pH、SOC 背景、Kd、弥散、附着/脱附/截留参数均不需为重现旧初值而改动。
- 两个 `site.toml` 解析后的配置差异仅为 `time.end` 和当前新增的 `switches.fertilizer_source_tracer=true`；上述 MEND 有效差异来自旧驱动脚本运行时设置，不存在于旧 site 文本中。
- `soil_initial_profile.csv` SHA256：`2848d4981a1809708c70c9af22fa87a13807bcf1d3918e2faade5303cfd3db25`。
- `organic_initial_pools.csv` SHA256：`72daeabb5edfb81f50b30188637cd6fbf8509cf3ebf7afabb2dc1e328d7dbe92`。

### 3.2 可审计的 fresh 初始化重现规则

权威历史变更算法见 `outputs/mend_doc_annual_20260916/driver/runDocAnnualTrial.m:19-78`；绝对质量记录见该运行 `ks030_year/initialization.json`。当前原生初始化先由 `buildSoilCarbon` 分配 SOC，再由 `reserveInitialDOM` 从 MOC/MON 预留总 DOM，因此 **DOM220 不需要额外添加、替换或从 MOC 再扣一次**。

在当前版本原生初始化之后、运行与 tracer 统计起点建立之前，进行一次显式初始化适配：

1. 取 top-down 各层 `SOC=POM1+POM2+MOC+MB+BD+DOM+(IOM+NEC)`，所有项统一为 mg C/g 干土。要求 ENZP1、ENZP2、EM 原初值为零。
2. 全剖面 `seed=5e-5*SOC`；`seedVolume=seed*bulkDensity/1000`。MOC 减 `3*seed`，ENZP1、ENZP2、EM 各置 seed，EP=2*seed。MON 减 bottom-up `seedVolume`，各酶 N=对应 `seedVolume/3`。初始碳分配参数缓存中的 MOC、EP、EM 同步更新。
3. 仅对单元中心 <20 cm 且网格恰覆盖 0–20 cm 的单元，将 MB、BD 与 mobile/attached/retained 的 C/N 各相乘 0.75；保留 active:dormant、mobile:attached 比例。减少的 C 加 MOC，减少的 N 加 MON；深层 MB/BD 不变。同步初始 MB/MOC 参数缓存和 Tot/EP 等派生聚合值。
4. 验证各层总 SOC 与总有机 N、全柱总 C/N、非有机接受态、水热/无机 N、管理输入保持不变；并将有意义的逐池、逐相初值与历史 `initial_state.mat` 对照。新增 tracer 初始来源标签应为零；不能把历史无标签 Y 整体覆盖新 schema。

单改 CSV 首行 `frac_MB=.06, frac_MOC=.74` **不能**精确复现：旧表是点剖面，改首行会通过插值改变 20–40 cm；而按默认 CN 重建 MOC-N 与原 MB-N→MON 守恒归還不同。当前 CSV 也没有独立酶初值列。所以应保存原 CSV＋显式一次性初始化规则，而不是暗改输入后声称初始 Y 相同。

旧初始化验收数值：

- MBC→MOC：`491.9200000000002 kg C/ha`；对应 microbial N→MON：`61.48999999999998 kg N/ha`。
- MOC→三酶：`46.41727546875 kg C/ha`；MON→三酶：`15.47242515625 kg N/ha`。
- 每类酶初始 0–20 cm 平均：`0.473 mg C/kg`。
- 初始 0–20 cm SOC=`9460 mg/kg`，MBC=`567.6 mg/kg`，总 DOM=`220 mg/kg`。
- 操作前后全柱初始 C/N：`3.09448503125`、`0.26379449599875715 g/cm²`。

## 4. 候选选择：本轮只改 CUE，不同时改 Kd

建议同一修复版源码上运行两组 fresh 六年：**重现旧有效参数的当前版基线**，以及仅 `mend_cue_20=0.25` 的温和候选。其余输入、初始化适配、示踪设置和输出相同。两组在同版本运行能区分新输运实现的影响与参数变化影响；当前代码不能使用旧六年曲线充当未改参数的精确对照。

### 4.1 CUE 0.2010932 → 0.25 的机制与限度

活动源码 `src/+soil/+microbe/mendSourceHour.m:83-91`：

```text
U = min(((Vg + Vm)/Yg) * D * MBA / (KD + D), D)
iHR=1: growth+maintenance respiration = U*(1-Yg)
```

在同一状态、环境和非 donor-cap 条件下，较高 Yg 降低 DOC 总摄取及生长相关呼吸；`U*Yg` 不直接下降。20 °C 同态摄取比为 `0.2010932/0.25=0.8043728`。因此它有同时缓解 DOC 消耗与 SOC 呼吸损失的明确机制。实际 Yg 按 `mendEnvironmentalParameters.m` 的温度函数变化，N 限制/溢流、酶、休眠和运输会反馈，不能把同态比率说成六年误差改善百分比。

`0.25` 是本轮可检验参数候选，不是新实测值，也不是已辨识的最优值。需同时验收 DON、SOC、MBC/MBN、CO2、无机 N 和守恒，不能只择最高 DOC。

### 4.2 Kd 与其他参数：保留为下一轮单独候选

`src/+soil/+dom/phaseState.m` 和 `equilibrate.m` 的总量约束为 `S=(theta+rho_b*Kd)*Cwater`。DOC/DON Kd 各乘 0.5 会提高同态自由相比例，但也增加摄取和迁移，且第 2.3 节显示表层 DOC 有总库存不足。因此本轮**不与 CUE 一起改 Kd**；将来需要以当前版基线分别做 DOC Kd、DON Kd、配对 Kd 的可辨识比较，保留深度形状和总量守恒。

BIOPARTICLE 附着增强或脱附减弱可减少表层细胞迁出，但可能增强表层 DOC 摄取。旧总微生物柱库存由 `5496.80` 增至 `7506.29 kg C/ha`，同时表层下降，不可简单解释成全柱生长崩溃。闭合柱层间迁移不等于离开 500 cm 底部的淋失。

### 4.3 历史反证和冻结项

- `outputs/mend_vg_annual_20260917/全年Vg试验结果.md`：仅 built Vg×0.5 曾使 0–20 cm 年末 DOC 从 `5.927` 提高至 `12.108 mg/kg`，SOC 从 `9.2962` 提高至 `9.3252 g/kg`，但 MBC 从 `528.31` 降至 `512.48 mg/kg`；不能无限下调 Vg 当作无代价改善。
- `outputs/mend_supply_annual_20260917/全年供碳试验结论.md`：manure 速率 `0.07→0.035→0.0175 d⁻¹` 时，年末 DOC `5.93→5.44→4.90 mg/kg`，SOC `9.296→9.244→9.167 g/kg`，全年有机 CO2 增加。较多 DOC 供给未保证终点 DOC 增加，故本轮保持 manure 速率 0.07 d⁻¹。
- 冻结水力和无机 N 参数、pH、容重、初始实测 SOC、施肥量/日期、作物输入和有效 C:N 合同；不增加浓度下限、不加虚构碳源、不关闭微生物运移、不重设每年土壤库。
- `linear_equilibrium` 下源 MEND 的 Kads/Kdes 已在 `advanceSourceCandidate.m` 置零。不能把 `mend_qmax_mg_g`、`mend_kba_g_mg`、`mend_kdes_mg_g_h` 当成本运行有效 DOC 吸附率定旋钮。

## 5. 版本隔离与交付验收

历史冻结 `+soil/+nitrogen/advanceAcceptedCandidate.m:197` 将最后接受水子步 `solver_diagnostics.QL` 用于外小时 DOM/BIOPARTICLE。当前文件的 2026-09-21 改动改为按 `accepted_transport_schedule` 顺序重放每个 QL/theta 子步；无机 N 使用接受时长加权液相通量，而不是最后子步通量。该变化会改变物理运输；旧六年的迁移数值与候选之间不是纯参数效应。

2026-09-22 标签修复的合同则是不改变总量动力学，证据和有限窗口见 `TRACER_LOCAL_REPAIR_20260922.md`。本次新基线与新候选都必须使用同一份包含该修复的源码；旧标签结果不混入新结果。

新六年启动前应检查 fresh 初始化审计、所有启用来源的逐库/逐相标签上下界与总量闭合、短窗示踪开关物理一致性、输出字段存在性与单位/维度/时间轴、恢复一致性。正式六年仍需完成之后才能评价终点 DOC/DON 与 SOC，不以成功启动替代率定完成。

## 6. 新试验初始化实施后的独立复核

对 `tools/local/initializeChangwuSixYearTrial.m`、`Prepare-ChangwuSixYearTrials.ps1`、`verifyChangwuSixYearInitialization.m` 作有界只读复核。酶播种、表层微生物各相 C/N 缩减、MOC/MON 归还及初始缓存更新与上述历史算法一致。实际生成的两份站点 TOML 经解析只有 `microbe.mend_cue_20` 一项差异；所有非 TOML 输入文件与当前原站点逐字节相同。重构参数计算的 Vm 与历史保存值在 double 精度下精确相同。标签 `emptyState` 初始化为零、不携带物理池容量缓存，因此守恒初始化适配不会需要额外标签重分配。

`outputs/sixyear_preflight_20260922/initialization-verification.json` 记录：新基线与 CUE 候选初始 Y 完全相同，参数仅 Yg_ref 不同；对历史检查共有 29 个顶级 Y 字段完全相同，其中 soil_organic、soil_microbes、soil_dom 整个 struct 完全相同，逐池最大绝对差为零；两个窗口内均有 7 次矿质肥事件和 7 次有机肥事件。

对其余 `vegetation`、`plant_acclimation` 两个非整结构相同字段，独立启动一次 MATLAB fresh 初始化并递归比较，**只发现以下四个叶差异**：

| 字段 | 历史初始形状 | 新初始形状 | 共有前缀比较 |
|---|---|---|---|
| `Y.vegetation.lai_state.LAI` | 9636×1 | 60272×1 | 完全相同，最大差 0 |
| `Y.vegetation.lai_state.SAI` | 9636×1 | 60272×1 | 完全相同，最大差 0 |
| `Y.vegetation.lai_state.hc` | 9636×1 | 60272×1 | 完全相同，最大差 0 |
| `Y.plant_acclimation.leaf_resp.Ta_history` | 8760×1 | 54792×1 | 新旧全为 NaN，共有前缀完全相同 |

解释：历史六年加载的是一年试验的初始 Y，因此留有一年预分配数组。新试验按 54,792 小时 fresh 构建；`src/+init/buildTime.m:78` 设置 `Nmsrmn=ceil(total_hours*1.1)`，`buildVegState.m:48` 据此构建 LAI/hc/SAI 序列；`buildInitialState.m:205` 按 tS 分配尚无接受时步的温度历史 NaN 缓冲。差异是符合新时间窗口的数组扩展，不是土壤、作物 C/N、物候初态或既有冠层观测值的改变。所有其他递归叶均完全相同，未发现未解释的初始物理状态差异。

独立诊断脚本及结果均保存在 `outputs/sixyear_preflight_20260922/`：

- `auditHistoricalInitialYDifferences.m`；
- `historical-initial-Y-differences.json`、`.mat`、`.log`。

诊断 MATLAB 正常退出（exit code 0），仅进行初始化和旧 MAT 对照，没有推进模拟，也没有编辑模型或初始化 helper。该结论支持进入两组同版本短窗及六年验证，不替代长窗结果本身。
