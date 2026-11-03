# N 模块重构 + 物理修正进度

> 维护：朱寄子星 · 启动 2026-05-13 · 最近更新 **2026-05-18**
> 新会话第一时间读本文件即可定位进度，需要细节时回看 §"已完成修复细节"。

---

## 当前状态（一句话）

阶段 1 重构完成 ✅ · 阶段 2 物理修正 + 守恒闭合完成 ✅ · **P3 全年 2021 跑验证完成** ✅ (NH4/NO3 守恒 0.0000%, 水 2.605%) · **阶段 3 (DAISY MM + A2 量纲修正 + NH4 dissolution pool) 完成** ✅ (Svol +60%)。下一阶段：**Phase C 根 N 主动吸收** (DSSAT NUPTAK 风格 + RTNH4/RTNO3 苹果偏好)。

### 2026-06-03 root_uptake demand-cap 收敛 bug 修复 ✅
[solveRootSurface.m:104](../src/+soil/+nitrogen/+root_uptake/solveRootSurface.m) 的 demand 模式通量收敛判据 `abs(I_mid − demand_per_root) < tol` 误用 `tol`（C_root **浓度**容差 [kg/m³]）去比**通量**差 [kg/(m·s)]——单位不匹配，通量差 (~1e-12) << 浓度容差 (~1e-6) → **第 1 次迭代即假收敛**，返回 C_root 中点的 I_DAISY（≈半土壤供应），与 demand 脱钩 → 根 N 吸收达 **10–50× demand** → 植物 N_labile 虚高。
- **修**（1 行）：通量判据改**相对容差** `abs(I_mid − demand_per_root) <= 1e-3 * demand_per_root`；C_root 区间判据 `(hi-lo)*0.5 < tol` 不动（该处 tol 单位正确，作兜底）。
- **配套**（runPlantNDemand.m）：demand 的 storage_deficit 改看**合并储备** `N_labile + N_storage`（合 LPJ-GUESS guess.cpp:1679 `nstore()`），reserve 满时 demand→0→uptake 真正归零。
- **全年验证** output_v25_F3fix（vs v24）：N_labile 末值 **47→5.1**；年根吸 **52.5→10.6 g N/m²**；supply_limit_factor（吸/需）mean **16.9→0.87**；NH4/NO3 守恒 **0.0000%**；不崩。
- **下游正确后果**：根吸砍 5× → 多约 42 g N/m² 留土壤（淋溶/反硝化，守恒账已平）。详见 [cn_coupling §12.30](cn_coupling_path_and_literature_audit.md) + P0 进度记忆。

### 2026-05-14 关键修复 — cv_init 锚定 (test28)

**问题**：1 月 jan2021 测试 NO3 累计残差 ~158 μg/cm² (6.99%)，KT-by-KT 单步残差却 ≈ 0。

**定位**：[no3-ledger] 输出显示 KT=24 时残差已经是 +157.67，到 KT=744 仍 +158.11 — 720 个 KT 中间严格闭合，**全部漂移都发生在 KT=1**。

**根因**：[solveNitrogen.m:91-101](../src/+soil/+nitrogen/solveNitrogen.m#L91-L101) 计算 `cv_init = cvno3OLD(KT=1)` 时用 `Theta_Solute_LL`（= KT=1 水量解之后的 θ），但配对的 C 是 IC C。两者隔着一次水量解，时步错位导致 KT=1 内的 θ 跳变（春初冻土初始化）被记成 cv_init 偏差，污染整个累计账本。

**修复**：`Theta_Solute_LL` → `Theta_Solute_L`（v11 用 Theta_LL_for_N_prev 给出的 t=0 真实 θ）。NH4 同步修。

**结果**：

| 物种 | 修前 (test27) | 修后 (test28) |
|---|---|---|
| NH4 | 0.21% | **0.0039%** |
| NO3 | 6.99% | **0.0197%** |

**通用教训**：任何"上一步储量"（X_OLD / cv_init / 残差锚）的 θ 和 C **必须取同一时步快照**，不能跨水量解混用。这条比反应速率公式还容易踩。

---

## 阶段 0: 包重构（已完成）

`src/+soil/+nitrogen/` 结构：

```
solveNitrogen.m              dispatcher (~250 行)
applyFertilizer.m            施肥即时事件
+chemistry/                  反应
  calculateRetardationFactors.m
  ammoniumVolatilization.m
  ureaHydrolysis.m
  nitrification.m            (O2 调制)
  denitrification.m          (O2 调制)
  freezeThawPartitioning.m   (冰相 N 包合/释放, K_segr=0.3)
+organic/                    有机 N
  mineralization.m
+transport/                  ADE
  effectiveLiquidWaterForSolute.m
  solveTransportStep.m       7-step 主驱动
  calculateBoundaryConcentrations.m
  calculateBoundaryFlux.m
  applyBoundaryConditions.m
  assembleSoluteParameters.m
  assembleMatrixCoefficients.m
  assembleEquation.m
  solveSoluteLinearSystem.m
  soluteDispersion.m
  pecletCourantNumber.m
  calculateStorage.m
+root_uptake/                植物 N 吸收
  nitrogenUptake.m           passive + Michaelis-Menten
+fertilizer/                 施肥事件
  initFertilizerController.m
  fertilizerFlux.m
```

清理：
- ~~`STEMMUS_N.m` (1387 行老 main loop)~~ 已删
- ~~`InNitrogenConcentrationDistribute.m` / `N_sub.m` / `N_Flux.m`~~ 已删
- ~~`sync_nitrogen_flux_aliases` / `get_nitrogen_flux_field`~~ 已删
- ~~`convertElementToNode.m`~~ 换用 `equations.interpolateElementToNode`

---

## 阶段 1: 物理修正（按完成顺序）

### ✅ P0a — Effective liquid water for solute
- 新文件 [effectiveLiquidWaterForSolute.m](../src/+soil/+nitrogen/+transport/effectiveLiquidWaterForSolute.m)
- 初版"扣冰"被撤回：Theta_LL 已经是液态水（含冰量在 Theta_II），无需再扣
- 当前实现：仅做 `min(θ_LL, POR)` 顶部 cap 和 ponding-aware 顶部节点处理
- 调用点 [solveNitrogen.m:31](../src/+soil/+nitrogen/solveNitrogen.m)

### ✅ P0a.5 — 冰相 N 包合 / 释放
- 新文件 [freezeThawPartitioning.m](../src/+soil/+nitrogen/+chemistry/freezeThawPartitioning.m)
- 冰生长 (Δθ_ice > 0): `entrap = K_segr=0.3 × C_liquid × Δθ_ice` 进 ice pool
- 冰融化 (Δθ_ice < 0): 按 melt 比例释放回液相
- 新 prognostic state: `SoilVariables.C_NH4_ice / C_NO3_ice` (NL×2)
- 每步守恒检查 `partition_diag.mass_conserved_check_*`

### ✅ P0b — O2 transport 冻土期 fallback
- 修改点 [solveNitrogen.m:102-117](../src/+soil/+nitrogen/solveNitrogen.m)
- 触发条件：`nnz(Theta_g < 0.02) / (2*NL) > 0.7`
- 行为：跳过 solveO2Transport，复用上一步 f_aerobic / f_anaerobic
- 首步 fallback：默认 `f_aerobic = f_anaerobic = 0.5`

### ✅ P1 — N 沉降时间序列输入
- Excel 加 4 列：
  - col 22 `NH4_dep` (mg/m²/day)
  - col 23 `NO3_dep` (mg/m²/day)
  - col 24 `NH4_bot` (mg/L)
  - col 25 `NO3_bot` (mg/L)
- 数据文件 `data/Meterology data2021-2022.xlsx` (368 行 × 25 列)
- [MeteoData.m](../src/+core/MeteoData.m) 读取 + lTimeScale 1=hourly/0=daily 切换
- [COMPASS.m:1665-1705](../src/COMPASS.m) NitrogenBoundary 用 MeteoData 时间序列，缺失列 fallback BoundaryCondition 标量
- 单位转换在 [calculateBoundaryConcentrations.m](../src/+soil/+nitrogen/+transport/calculateBoundaryConcentrations.m): FLUX → CONC via Precip × 3.6e4 L/(m²·h)
- ⚠️ 限制：当 Precip < 1e-7 cm/s（无雨）→ 干沉降丢弃，cvTop = 0

### ✅ P1.5 — 施肥事件 + 守恒账
- 3 个测试事件：`2021-04-05 NH4 50 kg/ha`、`04-15 NO3 50`、`04-25 urea 50`
- [applyFertilizer.m](../src/+soil/+nitrogen/applyFertilizer.m) 采用 **pre/post storage 差** 计算 `fert_diag.*_added`（单位与 `cv_nh4` / `cv_no3` 完全一致）
- COMPASS.m 累加到 `cum_nh4_in_mg_cm2` / `cum_no3_in_mg_cm2`

### ✅ P2 — 1 月 + 4 月跑验证
- test4 (修前): NH4 残差 131% → test6 (含 fert 修): 10.76% → test7 (含 Shys 修): **0.0026%**
- NO3 全程 ~4%

### ✅ P3 — 全年 2021 跑验证（完成 2026-05-18）
- StartTime=2021-01-01, EndTime=2022-01-01, Nmsrmn=8760
- 验收: 全年 N 守恒 < 5%, 冻融期不卡, N 轨迹符合物理（春融硝化峰 / 雨季淋溶 / 夏季根吸）
- **实际结果** (v33.2-dissolution-pool tag, output_v33_dissolve_fix):
  - **NH4 守恒残差 0.0000%** ✓ (远好于 5% 验收线)
  - **NO3 守恒残差 0.0000%** ✓
  - 水量残差 2.605% < 5% ✓
  - 冻融期不卡: 偶尔 thermal NOT-CONV 都被 SOR substep refinement 自动恢复
  - N 轨迹物理合理: Snit 5-9 月主峰 (温度+暖湿), Sden 7-10 月 (NO3 + WFPS), Smin 8 月单峰, Shys 4-5 月施肥脉冲, Svol 4-6 月跟着施肥
  - 非施肥期根区 NH4 中位数 7.8 mg/kg ✓ 落入 IC 5-8 目标
  - 苹果园 N 通量量级合理 (黄土高原高 N 输入 → Snit 1540, Sden 369, N2O 67 kg/ha)

### P4 — 详细规范化重写（计划）
对所有 +chemistry/+organic/+transport 按 Skill §2 doc 三件套 + ClimaLand 风格扫尾。

---

## 已完成修复细节（2026-05-13 晚 / 5-14 早）

### A. 速率常数批量重命名
caller 早已用可读名，但 callee 内部参数名仍是 K11/K21/K31 等老符号。这次把 callee 也对齐：

| 旧 | 新 | 文件 |
|---|---|---|
| K11 | K_nitrification | nitrification.m |
| K21 | K_denitrification | denitrification.m |
| K31 | K_mineralization | mineralization.m |
| Km1 | KM_denitrification | denitrification.m |
| KC1 | KC_denitrification | denitrification.m |
| Kurea1 | K_urea_hydrolysis | ureaHydrolysis.m |
| RHOKG1 | RHO_b_Kd | nitrification / denitrification / mineralization / calculateRetardationFactors / soluteDispersion / assembleSoluteParameters |

### B. ice partition 单位 bug（freezeThawPartitioning.m）
- 问题：`delta_NH4_to_ice += entrap_NH4` 累加的是 `[μg/cm³]`（浓度），ledger 期望 `[μg/cm²]`（深度积分质量）。隐性单位错配
- 修：签名加 `DeltZ` 参数，按 `entrap × DeltZ(ML) × 0.5` 求深度积分（每个 ND 贡献 0.5 × DeltZ，对应 cv_nh4 梯形积分定义）
- caller [solveNitrogen.m:45](../src/+soil/+nitrogen/solveNitrogen.m) 同步加 DeltZ

### C. cvShys 加入 NH4 input ledger（决定性守恒修复）
- 问题：`assembleSoluteParameters.m` `G1 = Shys + Smin - Svol - SinkS(:,:,1)` ——
  Shys（尿素水解 NH4 源）作为 ADE 源项进入传输方程，但 cumulative ledger 只加 `cvTop + cvBot + cvSmin`，Shys 隐性输入未记账
- 量化验证：test6 NH4 缺口 +649 μg/cm² ≈ test7 cvShys +648.91 μg/cm²，完美匹配
- 修：
  ```matlab
  % solveNitrogen.m per-step
  Shys_avg_now = (NitrogenVariables.Shys(:,1) + NitrogenVariables.Shys(:,2)) * 0.5;
  cvShys_now   = DeltZ' * Shys_avg_now;
  NitrogenFlux.cum_cvShys      += Delt_t * cvShys_now;
  NitrogenFlux.cum_nh4_in_mg_cm2 += Delt_t * cvShys_now;
  ```
- 在 [COMPASS.m](../src/COMPASS.m) breakdown print 加 `Shys (urea→NH4)` 行

### D. 守恒报告 print 精度
- `%+.2f` → `%+.4e`：原格式让 <0.005 μg/cm² 的值全打成 +0.00，误导排查（让人误以为 cvTop / ice partition 等过程 0）

### E. SinkS → NitrogenVariables 传递
- [COMPASS.m:935-946](../src/COMPASS.m)：`applyRootUptake` 返回的 `SinkS` 之前只用作局部变量，没赋给 `NitrogenVariables.SinkS` → N 输运永远看不到根吸
- 修后 test7 cvSinkSnh4 = 0.62、cvSinkSno3 = 26.96 ✓

### F. 根吸不能让 N 反流（HR 保护）
- [nitrogenUptake.m:77-86](../src/+soil/+nitrogen/+root_uptake/nitrogenUptake.m)：水量可在 HR 反向流回土（Srt < 0），但溶质不行
- 修：`Srt_uptake_only = max(Srt(ML, ND), 0)`，HR 反向流 N=0

### G. fertilizer 质量记账对齐 storage 单位
- [applyFertilizer.m](../src/+soil/+nitrogen/applyFertilizer.m)：施肥前后用相同 `Σ DeltZ × (θ + RHO_b_Kd) × C` 公式 snapshot，差值即 `fert_diag.nh4_added`，与 `cv_nh4` 完全同单位（μg/cm²）

### H. cTop/cBot 改名（与化学物质对齐）
- `cTop / cBot / cTop1 / cBot1` → `cTop_NH4 / cBot_NH4 / cTop_NO3 / cBot_NO3`（cTop1 不是"iter 1"，是 NO3）
- 6 个文件统一

---

## 验证结果对照

### test4（重构验证版，无任何物理修）
```
NH4 残差: 131%   NO3 残差: 97%   ← 完全不闭合
```

### test6（含 fert pre/post + 重命名）
```
NH4: ΔStorage=+911.5  Net=+262.8  Residual=+649 (10.76%) ✗
   Smin=0  Top=0  Bot=0  iceN=0  SinkS=0
   Snit=−2824.89  Cvnit_N2o=−57.65  Svol=−0.84
NO3: ΔStorage=+2841.3  Net=+2962.3  Residual=−121 (4.00%) ✓
   Snit(in)=+2824.89  Sden=−30.34
```
缺口 +649 μg/cm² 来自尿素水解 Shys 没记账。

### test7（含所有修，决定性）
```
NH4: ΔStorage=+911.23  Net=+911.40  Residual=−0.17 (0.0026%) ✓✓✓
   Top=0  Bot=0  Smin=0  Shys=+648.91  iceN=0
   Snit=−2824.6  SinkS=−0.62  Svol=−0.84  Cvnit_N2o=−57.65
NO3: ΔStorage=+2814.61  Net=+2935.17  Residual=−120.56 (3.95%) ✓
   Snit(in)=+2824.6  SinkS=−26.96  Sden=−30.25
```
✅ cvShys=648.91 完美匹配 test6 缺口（−1.5% 偏差仅由 KT 步内 Picard 收敛轻微差异）
✅ SinkS 起作用（NH4 0.62、NO3 26.96 反映冬季 → 4 月施肥后根吸增加）

---

## 仍有的小问题（不影响守恒，但用户问过）

1. **cvSmin = 0**：理论 4 个月应 ~50 μg/cm²，实测 +0.0000e+00 即便 %+.4e 精度也是 0
   - 排除：K_mineralization 数组非零（Excel col 24 = 0.05，÷86400 后 5.79e-7）
   - 排除：Theta_f 非零（van Genuchten @ h=-101.97 cm 计算得 ~0.3）
   - 待查：mineralization.m 内 eTheta / eT / ePH 中某个为 0
   - 验证方法：solveNitrogen 调完 mineralization 后加临时 fprintf 抓 Smin(NL,1) / Theta_Solute_LL / Theta_II
   - **不影响守恒**（即便 Smin 为 0，cvSmin 进 ledger 和 ADE G1 都是 0，自洽）

2. **Top deposition = 0**：calculateBoundaryConcentrations 通过 Precip 换算 FLUX→CONC，无雨期 cvTop=0
   - 物理设计（只算湿沉降），不是 bug
   - 干沉降需要 surface flux BC term，未来工作

3. **ice partition delta = 0**：理论上 1-2 月冻融期应有，实测 0
   - 待查：`Theta_II(new) - Theta_I(old)` 是否在 Picard 之间被刷新成相同值
   - 也可能 cfg.lFreezing 未真正打开（已确认 config 是 1）

4. **NO3 残差 3.95%**：在 5% 阈值内但偏大，可能是 ice partition NO3 漏算（一旦 #3 修通可能自然降）

---

## 测试历史

| 测试 | 时间 | 配置 | NH4 残差 | NO3 残差 | 关键变化 |
|---|---|---|---|---|---|
| test4 | 5-13 早 | 4 月跑 | 131% | 97% | 重构验证 |
| test6 | 5-13 晚 | 4 月跑 | 10.76% | 4.00% | + fert pre/post + 重命名 |
| test7 | 5-13 深夜 | 4 月跑 | **0.0026%** | 3.95% | + cvShys + ice DeltZ + SinkS 传递 + HR 保护 |
| year2021_test1 | 5-14 上午 | 全年 | TBD | TBD | 全年验证中 |

---

## 文件改动清单（按修改重要性）

| 文件 | 关键修改 |
|---|---|
| [src/+soil/+nitrogen/solveNitrogen.m](../src/+soil/+nitrogen/solveNitrogen.m) | + cum_cvShys 累计 + ice partition 传 DeltZ + per-process ledger |
| [src/+soil/+nitrogen/+chemistry/freezeThawPartitioning.m](../src/+soil/+nitrogen/+chemistry/freezeThawPartitioning.m) | 单位修：delta_*_to_ice 加 DeltZ 加权 |
| [src/+soil/+nitrogen/applyFertilizer.m](../src/+soil/+nitrogen/applyFertilizer.m) | pre/post storage snapshot 计算 fert_diag |
| [src/+soil/+nitrogen/+root_uptake/nitrogenUptake.m](../src/+soil/+nitrogen/+root_uptake/nitrogenUptake.m) | HR 保护 `Srt_uptake_only = max(Srt, 0)` |
| [src/+soil/+nitrogen/+chemistry/nitrification.m](../src/+soil/+nitrogen/+chemistry/nitrification.m) | K11/RHOKG1 → K_nitrification/RHO_b_Kd |
| [src/+soil/+nitrogen/+chemistry/denitrification.m](../src/+soil/+nitrogen/+chemistry/denitrification.m) | K21/Km1/KC1 → K_/KM_/KC_denitrification |
| [src/+soil/+nitrogen/+organic/mineralization.m](../src/+soil/+nitrogen/+organic/mineralization.m) | K31 → K_mineralization |
| [src/+soil/+nitrogen/+chemistry/ureaHydrolysis.m](../src/+soil/+nitrogen/+chemistry/ureaHydrolysis.m) | Kurea1 → K_urea_hydrolysis |
| [src/+soil/+nitrogen/+chemistry/calculateRetardationFactors.m](../src/+soil/+nitrogen/+chemistry/calculateRetardationFactors.m) | RHOKG1 → RHO_b_Kd |
| [src/+soil/+nitrogen/+transport/soluteDispersion.m](../src/+soil/+nitrogen/+transport/soluteDispersion.m) | RHOKG1 → RHO_b_Kd |
| [src/+soil/+nitrogen/+transport/assembleSoluteParameters.m](../src/+soil/+nitrogen/+transport/assembleSoluteParameters.m) | RHOKG1 → RHO_b_Kd |
| [src/+soil/+nitrogen/+transport/calculateBoundaryConcentrations.m](../src/+soil/+nitrogen/+transport/calculateBoundaryConcentrations.m) | FLUX→CONC 单位转换 + 字段改名 |
| [src/+soil/+nitrogen/+transport/effectiveLiquidWaterForSolute.m](../src/+soil/+nitrogen/+transport/effectiveLiquidWaterForSolute.m) | 新文件（最终等价 inline，保留作组织目的） |
| [src/COMPASS.m](../src/COMPASS.m) | breakdown print %+.2f→%+.4e + Shys 行 + SinkS→NitrogenVariables 传递 + fert 调用签名 |
| [src/+init/setInitialNitrogen.m](../src/+init/setInitialNitrogen.m) | 全 caller 用可读名 |
| [src/+core/MeteoData.m](../src/+core/MeteoData.m) | + NH4_dep/NO3_dep/NH4_bot/NO3_bot (Excel col 22-25) |
| [src/+core/FluxVariables.m](../src/+core/FluxVariables.m) | cTop/cBot/cTop1/cBot1 → cTop_NH4/cBot_NH4/cTop_NO3/cBot_NO3 |
| [data/Meterology data2021-2022.xlsx](../data/Meterology%20data2021-2022.xlsx) | 新文件 368 行 × 25 列 |

---

## 跨会话接力指令

新会话开始：
1. 读 `CLAUDE.md` 顶级警示 → invoke `compass-coding-standard` skill
2. 读本文件 §"当前状态" 和 §"测试历史" 表
3. 改 N 模块代码前先看 §"已完成修复细节" 避免重复踩坑（特别是 cvShys ledger、ice DeltZ、HR 保护）
4. 守恒挂了优先看 cum_*_in_mg_cm2 各分量是否有新源项漏记账（参考 Shys 修复模式）

---

## 关键命令

```bash
# 切 lNefc=1
sed -i 's/lNefc=0/lNefc=1/' data/config_stemmus_scope.txt

# 重置 IC
cp data/InitialConditions_jan01_backup.xlsx data/InitialConditions.xlsx

# 全年跑
rm -rf data/output_v33_N_apr2021 2>&1
F:/Matlab/bin/matlab.exe -batch "run_compass_window" > year2021_test1.log 2>&1

# 看守恒
grep -A 35 "N mass balance" year2021_test1.log | tail -45
```

---

## 阶段 3: 化学速率论修正 (2026-05-15 ~ 18, 完成)

| 子阶段 | 时间 | 内容 | tag |
|---|---|---|---|
| **3a** | 2026-05-15 | DAISY MM 改造 (硝化) — 线性 K_nit → MM | (无 tag) |
| **3b** | 2026-05-16 | Vmax_nit DAISY 校准 (0.5→5/3/1/0.5 μg/cm³/d) | (无 tag) |
| **3c** | 2026-05-17 | **A2 反硝化量纲修正** — 修 C² bug + 土基 MM 对称 | **v33.1-N-MM-A2** |
| **3d** | 2026-05-17 | NH3 挥发 L + pH 表层修正 | (3e 内一并) |
| **3e** | 2026-05-17/18 | **NH4 颗粒肥 dissolution pool** (APSIM-SoilN 风格) | **v33.2-dissolution-pool** |

### 3a — 硝化 DAISY MM 改造 (2026-05-15)

**问题**：旧线性 `Snit = K_nit·C` 在低 NH4 时让 Snit ≈ 0、高 NH4 时无饱和 → 文献偏离 + 数值刚性。

**修复**：`src/+soil/+nitrogen/+chemistry/nitrification.m` 重写为 DAISY Hansen 1991 风格：
- 显式源汇 (不走 SVnit 半隐式 — 因为 MM 在施肥脉冲后失配)
- $S_{nit} = V_{max}^{nit} \cdot e_T \cdot f_{aero} \cdot e_{pH} \cdot \dfrac{M}{K_m + M}$
- $M = (\theta + \rho_b K_d) \cdot C_{NH4}$ (土基, 含吸附)
- Excel 新列 `Vmax_nit` [μg/(cm³ 土·d)], `Km_nit` [μg/cm³ 土]
- ADE: $G_1 += -S_{nit}$ (NH4 汇), $G_2 += 0.98 S_{nit}$ (NO3 源 + 0.02 N2O 气逸出)

**校准** (2026-05-16, DAISY 源码 `nitrification_soil.C:111,114`):
- DAISY 默认: `k_10 = 5e-6 g N/cm³/d`, `k = 5e-5 g N/cm³`
- COMPASS Vmax_nit 从 0.5/0.3/0.1/0.05 (低 10-100×) 调到 5/3/1/0.5
- Km_nit = 50 = DAISY 默认

**结果**：非施肥期根区 NH4 中位数从 47 → 7.8 mg/kg ✓ 落入目标 5-8 区间。

### 3c — A2 反硝化量纲修正 (2026-05-17) 🌟

**bug 发现**：旧反硝化 `Sden = V × θ × mm × C` 实际 = `V × θ × C²/(C+Km)`, 不是 MM 而是 **C² 二阶反应**。

**修正**：[denitrification.m](../src/+soil/+nitrogen/+chemistry/denitrification.m) 重写为 DAISY 风格土体基 MM, 与硝化**对称**：
- $S_{den} = V_{max}^{den} \cdot e_T \cdot f_{anaer} \cdot \dfrac{M}{K_m + M}$
- $M = \theta \cdot C_{NO3}$ (NO3 不吸附 → 不含 ρKd 项, 与 NH4 唯一物理差异)
- Excel 重命名: `Km` → `Vmax_denit` [μg/(cm³ 土·d)], `KC` → `Km_denit` [μg/cm³ 土]
- 新数值: Vmax_denit = 1.5/1.0/0.3/0.1, Km_denit = 3 (= θ·10 mg/L 液)

**ADE 量纲一致性**（写入 `docs/cn_coupling_path_and_literature_audit.md §5.2`):

| 项 | 公式 | 单位 |
|---|---|---|
| 左 1 瞬态液相 | $\partial(\theta C)/\partial t$ | μg/(cm³ 土·s) ✓ |
| 左 2 瞬态吸附 | $\rho_b \partial s/\partial t$ | μg/(cm³ 土·s) ✓ |
| 弥散 | $\partial/\partial z(\theta D \partial C/\partial z)$ | θ 显式吸进 → 土基 ✓ |
| 对流 | $\partial(qC)/\partial z$ | q 本身土基, 不再 ×θ ✓ |
| **源汇** | $S_{nit}, S_{den}, S_{min}, ...$ | **必须土基** ✓ |

**全年验证 (2021)**：

| 物种 | 守恒残差 | Sden 年通量 |
|---|---|---|
| NH4 | **0.0000%** | — |
| NO3 | **0.0000%** | A1 → 497, **A2 → 374** kg/ha (-25%) |
| 总 N2O | — | A1 → 89.7, **A2 → 76.9** kg/ha (-14%) |

**贯通字段名**: `SoilProperties.m`, `setInitialNitrogen.m`, `solveNitrogen.m` 全改 `KM/KC` → `Vmax_denit/Km_denit`.

**论文意义**: 旧公式低浓度时高估反硝化（C² 退化为线性），高浓度时低估。A2 让低/高浓度区都符合 MM 物理。

### 3d — NH3 挥发参数修正 (2026-05-17)

[ammoniumVolatilization.m](../src/+soil/+nitrogen/+chemistry/ammoniumVolatilization.m) 2 处修：

| 项 | 改前 | 改后 | 理由 |
|---|---|---|---|
| **L** (扩散层厚) | 0.0025 m (裸地) | **0.05 m** (果园树冠下) | Sherlock & Goh 1985 原文针对植被覆盖 |
| **pH 节点** | `PH1(5)` (2cm depth) | **`PH1(1)`** (0-0.25cm 表层) | Svol 是表面通量, 应取表土 pH |

**数值影响**: L 改 → h_m × (0.0025/0.05)^0.2 = ×0.72; pH 几乎不变 (全剖面 pH 差 0.013).

**真正 Svol 偏低根因**: 施肥**深度 5 cm** 不是表面撒, 即使 dissolution pool 完美工作, 释放出 NH4 在 5 cm 处, Svol 公式只看表层节点 NN → 物理锁死. 修 L+pH 是物理正确但量级有限.

### 3e — NH4 颗粒肥 Dissolution Pool 🌟🌟 (2026-05-17/18)

**问题陈述**: 模型把 NH4 颗粒肥**瞬间均匀稀释**到 5 cm 施肥层 → 表层 NH4 液相浓度峰值被低估 (实际颗粒表面 ~1000 mg/L vs 模型 45 mg/L) → Henry 驱动 NH3 挥发通量低估 10-30 倍.

**实施 (APSIM-SoilN 风格, Probert 1998 Agric Systems 56(1):1-28)**:

| 项 | 含义 |
|---|---|
| **新模块** | [nh4FertDissolution.m](../src/+soil/+nitrogen/+fertilizer/nh4FertDissolution.m) — 60 行, 镜像 ureaHydrolysis 结构 |
| **新 prog state** | `SoilVariables.C_NH4_granule` [μg/cm³, NN×1] — 颗粒肥缓冲池 |
| **公式** | $\Delta G = C_{gran} \cdot (1 - \exp(-k_{eff} \cdot \Delta t))$, $k_{eff} = 5 \cdot \text{WFPS} \cdot K_{dis} \cdot Q_{10}^{(T-20)/10}$ |
| **路径** | 施肥 → granule pool → 一阶溶解 → dissolved NH4 → 走原有 ADE/反应 |
| **新参数** | Excel `K_nh4_dissolve` 列 [1/d], 7 行全填 0.3 (APSIM 区间 0.1-0.5 中位) |
| **ADE 耦合** | $G_1 += S_{dis}$ (NH4 源, 与 Shys/Smin 并列) |
| **守恒** | granule 纳入 cvnh4 储量; granule→dissolved 是**内部转移**, 自动闭合 |

**文献支撑**:
- **Probert M E et al. 1998** Agric Systems 56(1):1-28 — APSIM-SoilN urea-hydrolysis transient pool template
- **Sommer S G & Hutchings N J 2001** EJA 15(1):1-15 — Mean-field 假设是 NH3 挥发低估主因
- **Génermont S & Cellier P 1997** Agric Forest Meteorol 88:145-167 — Surface vs subsurface NH4 pool separation
- **Lloyd J & Taylor J A 1994** Func Ecol 8:315-323 — Q10=2.5 温度响应

**实施清单 (Step 1-7)**:

| Step | 文件 | 改动 |
|---|---|---|
| 1 | `+fertilizer/nh4FertDissolution.m` | 新文件 |
| 2 | `data/soilproperty.xlsx` col 27 | 加 K_nh4_dissolve, **7 行全填 0.3** |
| 3 | `+core/SoilProperties.m` | `K_nh4_dissolve_type/K_nh4_dissolve` 字段 + Excel 解析 + flip |
| 4 | `+core/StateVariables.m`, `+init/defineSoilVariables.m` | `C_NH4_granule` (NN×1) prog state |
| 5 | `+soil/+nitrogen/applyFertilizer.m` | NH4 路由到 granule; `fert_diag.nh4_added` 含 granule 增量 |
| 6 | `+soil/+nitrogen/solveNitrogen.m` | 调 `nh4FertDissolution` → Sdis 进 ADE G1; granule 扣减放 Picard iter **外** |
| 7 | `+transport/assembleSoluteParameters.m` | $G_1 += S_{dis}$; `nargin<43` backward-compat |
| 7b | `+transport/solveTransportStep.m` | 传 Sdis 给 assemble |
| 7c | `+transport/calculateStorage.m` | 4th varargin = `C_NH4_granule` 累加到 cvnh4 |

**踩过的两个大坑 (写入 cn_coupling_path md §12.6)**:

1. **`nargin < 45` bug**: assembleSoluteParameters 实际 43 个参数 (Sdis 是第 43). 误用 45 → Sdis 一直 reset 为 0 → granule 释放但 dissolved 不增 → mass balance 漏 651 μg/cm². **改 43 修复**.
2. **Excel 加列只填前 4 行**: 后 3 行 (300/500/700 cm) 留空 → Excel 解析 dropna 删整行 → 深层 VG 参数 (ks, n, alpha) 丢 → K 飙升 5× → 底部漏水 +110 mm/年 → 水量残差 6.59%. **7 行全填 0.3 修复, 水量回到 2.605%**.

**A2 (无 pool) vs Dissolution Fix 全年对比** (Changwu apple orchard, 2021):

| 指标 | A2 | Dissolution Fix | 变化 |
|---|---|---|---|
| 水量残差 | 2.605% | **2.605%** | ✓ 一致 (pool 不动水) |
| BottomFlux | 41.56 mm | **41.56 mm** | ✓ 一致 |
| NH4 守恒 | 0.0000% | **0.0000%** | ✓ |
| NO3 守恒 | 0.0000% | **0.0000%** | ✓ |
| Snit 年通量 | 1788 kg/ha | **1540 kg/ha** | ↓ 14% (NH4 被 Svol 拿走更多) |
| Sden 年通量 | 374 kg/ha | 369 kg/ha | ✓ 不变 |
| **Svol (NH3 挥发)** | **0.125 kg/ha** | **0.200 kg/ha** | **+60%** ✓ 物理收益 |
| NO3 末态储量 | 27133 μg/cm² | 25225 μg/cm² | ↓ (Snit 减少传导) |

**仍偏低**: 0.2 kg/ha vs 文献 15-94 区间. 根因不在 dissolution pool, 而在施肥**深度 5 cm**, NH4 在 5 cm 释放但 Svol 公式只看 0-0.25 cm 表层节点. 论文里写 limitation, 不再调.

---

## 阶段 4: Phase C 根系 N 主动吸收 (2026-05-18 启动, 进行中)

详细 plan 在 [`cn_coupling_path_and_literature_audit.md §12.7`](cn_coupling_path_and_literature_audit.md).

**M1 (1-2 周, 当前下一步)**: DSSAT-NUPTAK 风格根 N 主动吸收
- M1.1 新文件 `src/+plant/+root/calculateRootNUptake.m`: $U_N = \max(RWU, 0) \cdot C_l + RAI \cdot V_{max} C_l / (K_m + C_l) \cdot f_{water} \cdot f_{demand}$
- M1.2 物种偏好系数 RTNH4=0.7, RTNO3=0.3 (苹果偏氨)
- M1.3 集成到 SinkS (替代当前纯被动质量流)
- 期望: NH4:NO3 根吸收比例从 1:41 翻到 3:7 → NO3 累积下降

**M2 (2-3 周)**: leaf N pool + Vcmax 反馈 (CLM5 LUNA 风格)
**M3 (2-3 周)**: C-N 反分配 + 衰老 retranslocation

完成后是论文 Tier 2 (GMD / Agricultural Forest Meteorology), 第一个真正"SPAC + SCOPE + N 闭环"果园模型.

---

## Tag 历史

| Tag | 日期 | 内容 |
|---|---|---|
| v33.0-thomas-sor-fix | 2026-05-12 | Thomas TDMA log-SOR + Tol_Water=50 (水分) |
| **v33.1-N-MM-A2** | **2026-05-17** | **A2 反硝化量纲修正 + N 模块完整重构 + DAISY 校准** (acbfd72) |
| **v33.2-dissolution-pool** | **2026-05-18** | **NH4 颗粒肥 dissolution pool + Phase C 启动 plan** (10d2ce9) |
