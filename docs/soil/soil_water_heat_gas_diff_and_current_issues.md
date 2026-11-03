# 土壤水-热-气模块与 STEMMUS_SCOPE-main 原版差异清单

对比基准：`F:\STEMMUS_SCOPE-main\src`

当前代码：`F:\COMPASSV1.0\src`

记录日期：2026-05-05

本文只记录当前需要逐一验证的差异，不代表所有差异都必须回退。判断顺序建议是：先确认物理量定义，再确认数值稳定项，最后确认输出字段。

## 1. 当前优先验证目标

1. 水分三相关系是否守恒：
   `Theta_UU = Theta_LL + RHOI/RHOL * Theta_II`，当前常数约为 `RHOI/RHOL = 0.92`。
2. 气相孔隙定义是否合理：
   `Theta_g = POR - Theta_LL - Theta_II`。
3. 干空气压力场是否符合边界条件：
   底部定压时，深层气压应由 `BtmPg` 控制；如果启用静气压底边界，则底部应略高于表层。
4. 热方程是否仍被当前额外稳定项改变：
   特别是 `KTT/VTT` 的人工弥散和限幅。
5. 非收敛是否来自热方程相变热容、干空气方程压力矩阵，还是水分导数。

## 2. 水分模块差异

### 2.1 `UpdateSoilWaterContent.m`

当前文件：`src/UpdateSoilWaterContent.m`

原版文件：`F:\STEMMUS_SCOPE-main\src\UpdateSoilWaterContent.m`

主要差异：

1. 表面张力温度修正的实现方式不同。

原版临时覆盖主状态：

```matlab
SoilVariables.hh = COR .* SoilVariables.hh;
calculateHydraulicConductivity(...)
SoilVariables.hh = SoilVariables.hh ./ COR;
```

当前代码不直接覆盖 `hh`，而是写入：

```matlab
SoilVariables.hhU = COR .* hh;
SoilVariables.hU  = COR .* SoilVariables.h;
```

这能避免污染主状态，但不是逐句原版。

2. 无滞后时气相孔隙定义不同。

原版：

```matlab
Theta_V = POR - Theta_UU;
Theta_g = Theta_V;
```

当前：

```matlab
Theta_V = POR - Theta_LL - Theta_II;
Theta_g = Theta_V;
```

当前写法更符合冻土三相体积关系：

```text
POR = Theta_LL + Theta_II + Theta_g
```

但这不是原版。

3. 当前滞后修正后会同步更新冰相和导水率。

当前在滞后扫描线改变 `Theta_LL` 后，会保持 `Theta_UU` 不变并反算：

```matlab
Theta_II = max(0, (Theta_UU - Theta_LL) * RHOL / RHOI);
```

随后同步更新 `Se`、`KL_h`、`KfL_h`、`KfL_T`。原版没有这一整段同步逻辑。

4. 当前只在非饱和区执行滞后修正。

当前使用：

```matlab
if hh(MN) < 0
    % hysteresis
else
    % saturated, skip hysteresis
end
```

原版没有这个 `hh < 0` 限制，只按 `IH` 分支判断。

### 2.2 `calculateMatricCoefficients.m`

当前文件：`src/+soil/+water/calculateMatricCoefficients.m`

原版文件：`F:\STEMMUS_SCOPE-main\src\+soilmoisture\calculateMatricCoefficients.m`

当前已经按原版恢复导数交叉割线写法：

```matlab
DTheta_UUh = (Theta_LL - Theta_L) * COR / DhUU;
DTheta_LLh = (Theta_UU - Theta_U) * COR / DhU;
```

并且已经删除之前的小分母保护。后续验证重点：

1. `DTheta_LLh` 是否还会在 reject 或冻融段爆炸。
2. 冻融关闭时，`DTheta_LLh` 和 `DTheta_UUh` 是否退化一致。
3. 水分矩阵容量项和热方程潜热项是否使用了与原版一致的导数。

### 2.3 冰阻抗

当前 `src/+init/defineSoilVariables.m` 中：

```matlab
SoilVariables.Imped = repmat(3, NL+1, 1);
```

已从之前的 `7` 改为 `3`。这更接近原版/常用设置。后续需要验证冻融期 `KfL_h = KL_h * 10^(-Imped * Ratio_ice)` 是否导致水分变化过强或过弱。

## 3. 干空气模块差异

### 3.1 `calculateDryAirParameters.m`

当前文件：`src/+soil/+gas/calculateDryAirParameters.m`

原版文件：`F:\STEMMUS_SCOPE-main\src\+dryair\calculateDryAirParameters.m`

主要差异：

1. 储气项不同。

原版等价写法：

```matlab
stor = POR + (Hc - 1) * Theta_LL;
```

当前：

```matlab
theta_g_node = max(POR - Theta_LL - Theta_II, 1e-14);
stor = theta_g_node + Hc * Theta_LL;
```

这是当前建议保留的物理修正，因为冰占据孔隙后不能再当作气相空间。

2. 干空气方程使用了节点插值。

当前对以下单元量插值到节点：

```matlab
D_Vg_n = interpolateElementToNode(D_Vg, DeltZ);
V_A_n  = interpolateElementToNode(V_A, DeltZ);
QL_n   = interpolateElementToNode(QL, DeltZ);
```

原版直接使用层量 `D_Vg(i)`、`V_A(i)`、`QL(i)`。这会改变 `ND=1/2` 两侧的干空气矩阵系数。需要单独验证。

3. 当前对气压梯度限幅。

当前：

```matlab
DPgDZ = min(max(DPgDZ, -1e5), 1e5);
```

原版没有这个限幅。该项是数值保护，不是原版物理公式。

4. 当前对 `Kg` 和扩散系数设置正下限。

当前在非积水条件下使用：

```matlab
safe_Kg = max(Kg, 1e-15);
safe_diff = max(D_V + D_Vg_n, 1e-18);
```

原版没有这些下限。它们会增强矩阵可解性，但可能改变低气相孔隙条件下的压力响应。

5. 当前对干空气对流项加入人工弥散。

当前：

```matlab
Kah = max(Kah, abs(Vah) * DeltZ / 2);
KaT = max(KaT, abs(VaT) * DeltZ / 2);
Kaa = max(Kaa, abs(Vaa) * DeltZ / 2);
```

原版没有。需要验证它是否影响深层气压长期偏低或锯齿。

6. 当前输出兼容字段。

当前同时输出：

```matlab
AirVariables.QL_A
AirVariables.QL_a
```

其中 `QL_a` 是为了和最新版能量方程字段名一致，`QL_A` 是当前其它模块已有调用名。

### 3.2 `calculateBoundaryConditions.m`

当前文件：`src/+soil/+gas/calculateBoundaryConditions.m`

原版文件：`F:\STEMMUS_SCOPE-main\src\+dryair\calculateBoundaryConditions.m`

主要差异：

1. 当前加入了随当前表层气压变化的静气压底边界。

当前：

```matlab
BtmPg = top_pg_now + rho_air * g_cgs * depth_cm;
```

该值只有当底部气压边界为定压时才生效：

```matlab
if NBCPB == 1
    RHS(1) = BtmPg;
end
```

如果 `NBCPB == 2`，则是底部通量边界：

```matlab
RHS(1) = RHS(1) + BCPB;
```

此时前面计算的 `BtmPg` 不会施加到底部。

2. 当前表层边界 fallback 与原版略不同。

原版直接使用：

```matlab
RHS(n) = TopPg(KT);
```

当前先检查 `TopPg(KT)` 是否存在和有限，然后才用，否则回退到 `BtmPg`。

3. 当前对 RHS 做极值截断。

当前：

```matlab
RHS = min(max(RHS, -1e30), 1e30);
```

原版没有。该项主要防止数值溢出。

## 4. 热方程模块差异

### 4.1 `calculateEnergyParameters.m`

当前文件：`src/+soil/+thermal/calculateEnergyParameters.m`

原版文件：`F:\STEMMUS_SCOPE-main\src\+energy\calculateEnergyParameters.m`

主要差异：

1. 当前修复了原版函数化代码中的 `CTh` 断行问题。

原版函数化代码中 `else` 分支实际写成：

```matlab
EnergyVariables.CTh = thermal_term * DTheta_LLh;
+ vapor_term + air_term;
```

第二行不会加进 `CTh`。当前代码把三项合并为完整公式：

```matlab
CTh = thermal_term * DTheta_LLh + vapor_term + air_term;
```

这更接近老版 `EnrgyPARM` 的物理公式，但不是逐句最新版原版。

2. 严格冻融热容分支条件不同。

原版：

```matlab
if abs(SAVEDTheta_LLh - SAVEDTheta_UUh) ~= 0
```

当前：

```matlab
if lSFCC && lFreezing && abs(SAVEDTheta_LLh - SAVEDTheta_UUh) ~= 0
```

当前多了 `lSFCC` 和 `lFreezing` 限制。

3. `SFCC == 0` 分支条件不同。

原版：

```matlab
if ModelSettings.SFCC == 0
```

当前：

```matlab
if ~lSFCC && lFreezing
```

当前多了 `lFreezing` 限制。

4. `CTT <= 0` 回退保护已删除。

之前当前代码有：

```matlab
if CTT <= 0
    CTT = generic storage;
end
```

现已删除，行为更接近原版。

5. 当前对热方程通量/系数使用节点插值。

当前：

```matlab
D_Vgg = interpolateElementToNode(D_Vg, DeltZ);
V_AA  = interpolateElementToNode(V_A, DeltZ);
QLL   = interpolateElementToNode(QL, DeltZ);
QVV   = interpolateElementToNode(QV, DeltZ);
Qaa   = interpolateElementToNode(Qa, DeltZ);
```

原版直接使用 `D_Vg(i)`、`V_A(i)`、`QL(i)`、`QV`、`Qa`。特别是 `QL/QV/Qa` 本身是层间通量，插值到节点后会改变同一单元两端的 `VTT`，这是需要验证的重点。

6. 当前对 `KTT/VTT` 加了额外稳定项。

当前：

```matlab
KTT = max(KTT, Lambda_eff);
KTT = max(KTT, abs(VTT) * DeltZ / 2);
VTT = min(max(VTT, -vtt_limit), vtt_limit);
```

原版没有。该项可能改善收敛，但也会改变热对流强度。

7. 当前 `CTg` 与最新版一致，但与最老版 `EnrgyPARM` 不同。

当前/最新版：

```matlab
CTg = (c_L*RHOL + c_a*Hc*RHODA) * KfL_h * TT - c_L * Srt * TT;
```

最老版：

```matlab
CTg = (c_L*RHOL + c_a*Hc*RHODA) * KL_h * TT;
```

差异是当前使用含冰阻抗的 `KfL_h`，并扣除了根系吸水热项 `Srt`。

8. 当前输出字段已按最新版补齐。

当前输出：

```matlab
EnergyVariables.Qtot
EnergyVariables.QL_h
EnergyVariables.QL_T
EnergyVariables.QL_a
EnergyVariables.QVH
EnergyVariables.QVT
EnergyVariables.QVa
```

旧字段 `EnergyVariables.QV`、`QV_h`、`QV_T`、`QVA` 已不再作为输出字段保留。局部变量仍用于内部计算。

### 4.2 `solveThermalEquations.m`

需要后续检查：

1. 温度异常回退条件是否与原版一致。
2. `TT/T/TOLD`、`TT_CRIT/T_CRIT`、`hh_frez/h_frez` 的时步更新是否完整。
3. 子步回退时是否只回退独立状态量，还是连派生系数一起回退。

## 5. 输出字段状态

上一次已有结果目录 `data/output_airfix_hydrostatic_1month` 是旧命名输出，验证结果为：

```text
SoilUnfrozenWater = SoilMoisture + 0.92 * SoilIce
max residual = 1.39e-15
```

说明该旧结果中：

```text
Sim_output_SoilMoisture       = Theta_LL
Sim_output_SoilUnfrozenWater  = Theta_UU
Sim_output_SoilIce            = Theta_II
```

当前代码已经改为：

```text
Sim_output_SoilMoisture       = Theta_UU
Sim_output_SoilUnfrozenWater  = Theta_LL
Sim_output_SoilIce            = Theta_II
```

下一次运行后应验证：

```text
SoilMoisture = SoilUnfrozenWater + 0.92 * SoilIce
```

## 6. 后续逐项验证顺序

建议按下面顺序验证，不要一次性改多处：

1. 用当前代码跑 1 天，确认输出命名和三相守恒。
2. 跑 1 个月，检查气压剖面：
   - `P_bottom - P_top` 是否接近静气压边界；
   - 是否仍出现深层长期低压。
3. 关闭热方程额外稳定项，仅验证 `KTT/VTT` 对收敛和温度锯齿的影响。
4. 把热方程中的 `QL/QV/Qa` 节点插值改回层量，单独验证温度和水汽通量是否更平滑。
5. 对干空气模块做三组边界测试：
   - `NBCa_bottom=1` 定压；
   - `NBCa_bottom=2` 通量；
   - 固定表层/底部压力的静态测试。
6. 冻融期重点输出：
   - `Theta_UU`、`Theta_LL`、`Theta_II`；
   - `DTheta_LLh`、`DTheta_UUh`；
   - `KfL_h`、`Ratio_ice`；
   - `CTT`、`CTT_PH`、`CTT_LT`；
   - `P_gg`、`Caa`、`RHS`。

## 7. 当前判断

目前第一优先级仍是土壤水-热-气内部一致性，SCOPE 暂时不是第一嫌疑。SCOPE 输出可以继续观察，但如果三相含水量、气相孔隙和干空气压力场还没有稳定，不宜先校正冠层或辐射部分。

## 8. 2026-05-06 Picard 系数刷新试验与实测对比

本节记录一次专门针对 `COMPASS.m` Picard 内部水-气-热系数刷新顺序的试验。试验目的不是重新标定参数，而是检查更严格的系数同步是否值得作为默认主循环。

### 8.1 对比组合

观测数据：

```text
tests/changwu_obs_cleaned.xlsx
```

模拟窗口：

```text
2021-01-01 00:00 到 2021-02-01 00:00
```

四组模拟输出：

| 组别 | 输出目录 | 说明 |
|---|---|---|
| Original | `data/原版代码数据` | 原版函数输出 |
| Baseline | `data/output_compare_baseline_1month` | 当前推荐基线顺序 |
| StrictRefresh | `data/output_compare_strict_refresh_1month` | 水后、气后、热后均完整刷新 |
| MidRefresh | `data/output_compare_mid_refresh_1month` | 水后、气后刷新，去掉热后刷新 |

生成脚本：

```text
tests/compare_four_soil_moisture_obs.m
```

输出文件：

```text
docs/figures/four_run_obs_compare/four_run_obs_aligned_values.csv
docs/figures/four_run_obs_compare/four_run_obs_stats_by_depth.csv
docs/figures/four_run_obs_compare/four_run_obs_weighted_summary.csv
```

### 8.2 加权总体指标

| Run | N | Weighted RMSE | Weighted MAE | Weighted Bias | Mean KGE | Mean R |
|---|---:|---:|---:|---:|---:|---:|
| Original | 6135 | 0.07970 | 0.05470 | +0.01420 | -1.7347 | -0.1759 |
| Baseline | 6135 | 0.02896 | 0.01831 | +0.00736 | -0.5159 | -0.1269 |
| StrictRefresh | 6135 | 0.03028 | 0.01963 | +0.00716 | -1.8261 | -0.2497 |
| MidRefresh | 6135 | 0.02894 | 0.01828 | +0.00684 | -1.6302 | -0.1523 |

按单月土壤含水量 RMSE 看，`MidRefresh` 与 `Baseline` 几乎持平，`MidRefresh` 仅小约 `1.8e-5`。这个差异远小于数值稳定性代价，不足以支持切换默认主循环。

### 8.3 图形结果

四组时间序列与实测对比：

![four run obs timeseries](figures/four_run_obs_compare/four_run_obs_timeseries.png)

分深度 RMSE：

![four run obs rmse by depth](figures/four_run_obs_compare/four_run_obs_rmse_by_depth.png)

均值剖面与总体 RMSE：

![four run obs profile summary](figures/four_run_obs_compare/four_run_obs_profile_summary.png)

### 8.4 数值稳定性对比

1 个月运行稳定性：

| 组别 | 计算时间 | StepFailed | Reject trials | 涉及 reject 的小时数 | 备注 |
|---|---:|---:|---:|---:|---|
| Baseline | 约 320 s | 0 | 0 | 0 | 最稳定 |
| StrictRefresh | 约 1051 s | 0 | 52 | 18 | 更硬，局部温度/气压差异放大 |
| MidRefresh | 约 1402 s | 0 | 103 | 55 | 比 StrictRefresh 更难收敛 |

极值变化也说明刷新版本会显著改变冻结期内部状态：

| 变量 | Baseline 范围 | StrictRefresh 范围 | MidRefresh 范围 |
|---|---:|---:|---:|
| SoilMoisture | 0.166-0.284 | 0.179-0.368 | 0.180-0.397 |
| SoilTemperature | -25.72-15.31 C | -25.72-30.79 C | -25.72-75.41 C |
| AirPressure | 84.01-88.46 | 69.45-94.02 | 71.94-89.20 |

其中 `MidRefresh` 在 1 月末出现 60 C 以上的内部土温异常，不能作为正式科学结果。

### 8.5 结论

本次试验支持以下判断：

1. `Original` 与实测差距明显更大，不能作为当前推荐结果。
2. `Baseline` 相比 `Original` 显著改善土壤含水量 RMSE，并且一个月内没有 reject、没有 StepFailed。
3. `StrictRefresh` 和 `MidRefresh` 虽然尝试提高 Picard 内部水-气-热系数同步，但会明显增加非线性刚性，导致更多子步拒绝和局部物理异常。
4. `MidRefresh` 的总体 RMSE 略低于 `Baseline`，但优势只有约 `1.8e-5`，同时计算时间约为 `Baseline` 的 4.4 倍，并出现更强的温度异常，因此不采用。

最终决定：

```text
当前正式主循环继续使用 Baseline 顺序。
不采用 StrictRefresh。
不采用 MidRefresh。
保留 Tol_Water = 100。
保留 min-dt 失败回退到上一个完整外步状态的逻辑。
```

后续如果继续研究系数同步，建议只作为独立实验分支，并先关闭 Anderson 或增加专门的气压/冻结诊断输出，不应直接并入正式主循环。

## 9. 全年 Imped=7 运行中的难收敛段记录

后台全年运行目录：

```text
data/output_fullyear_imped7_20260506_194252
```

截至当前日志，`KT=1850-1851` 是明显难点，约对应 `2021-03-19` 凌晨。`KT=1850` 最终用 `641` 个内部子步完成一小时；`KT=1851` 继续出现大量 reject 和小步推进。

主要现象：

1. 非收敛原因主要来自水分方程，不是热方程。日志中常见 `Energy=1e-3~1e-2`，但 `Water` 在失败时可达 `1e3~1e4 cm`。
2. 成功接受的小步也常接近当前水分保护阈值，例如 `Water=60-99 cm`，说明该段只是被 `Tol_Water=100` 放行，并不是严格水分收敛。
3. 前序 `nonconv_diag_KT01843_dt1s.mat` 到 `nonconv_diag_KT01849_dt1s.mat` 显示土壤处在冻结临界和近饱和组合状态：`TT min` 约 `-0.08 C`，`iceMax` 约 `0.236`，`hh max` 接近 `0`，`theta max` 接近饱和含水量。
4. 气象上该段处于 3 月中旬湿冷过程之后，有高湿、低辐射、小雨或降雨日背景；但 `KT=1850-1851` 并不是简单的瞬时降雨峰值问题，更像是前几个小时冻结-融化-湿润状态累积后的水力边界/冻结水力参数硬点。

下一步最值得改的数值策略：

```text
在连续 reject 的小时内，不要两次成功就立即放大内部步长。
应临时锁住最近稳定的小步长，等通过冻结/饱和临界段后再恢复增长。
```

建议实现方向：

1. 在 `COMPASS.m` 内部子步逻辑中记录当前外步内的连续 reject 次数和最近 reject 次数。
2. 如果同一 `KT` 内频繁出现 `3.5 s` 成功、`7-28 s` 失败的模式，则进入 hard-substep 模式。
3. hard-substep 模式下将 `sub_dt_next` 固定为最近稳定步长，或要求更多连续成功次数后才允许增长，例如从当前 `2` 次提高到 `8-10` 次。
4. 当外步完成或连续较大步长稳定通过后，再退出 hard-substep 模式，恢复原增长规则。

该修改只改变内部时间步调度，不直接改变水热气物理公式。它的目标是减少冻结/近饱和临界段中的无效放大-失败-回退循环，提高全年运行效率，并避免长时间卡在少数小时。
