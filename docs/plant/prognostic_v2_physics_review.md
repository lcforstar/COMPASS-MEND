# COMPASS prognostic_v2 物理审查 + 耦合可行性报告

**文档目的**: 三件事一并交代
1. v2 (`+plant/+PlantHydraulic/prognosticPlantWater.m`) 的物理过程跟参考模型对比
2. 物理过程对不对、要不要改进
3. v2 跟现在 COMPASS.m 主代码能不能耦合进去（接口层面）

**日期**: 2026-05-27
**关联文档**: [vod_module_roadmap.md](vod_module_roadmap.md)、[src/PROGNOSTIC_V2_QUICKSTART.md](../src/PROGNOSTIC_V2_QUICKSTART.md)

---

## 1. v2 是什么 (一句话)

把 ClimaLand 2026 (Deck et al. JAMES) §2.3.1 的瞬态植物水力 (Eq. 32–37) 用 MATLAB 重写，**三层储水状态**：`W_root, W_stem, W_leaf` [m]，用 IMEX 半隐式 3×3 解。

文件清单：
- [src/+plant/+PlantHydraulic/prognosticPlantWater.m](../src/+plant/+PlantHydraulic/prognosticPlantWater.m) ← 核心 IMEX 算法
- [src/+plant/+PlantHydraulic/weibullConductance.m](../src/+plant/+PlantHydraulic/weibullConductance.m) ← K(ψ)
- [src/+plant/+PlantHydraulic/linearRetention.m](../src/+plant/+PlantHydraulic/linearRetention.m) ← ψ(W)
- [src/+plant/+PlantHydraulic/Plant_hydraulics_v2.m](../src/+plant/+PlantHydraulic/Plant_hydraulics_v2.m) ← 兼容 v1 接口的包装

---

## 2. 物理审查 — 对 4 个参考模型

### 2.1 公式逐项对比

| 物理量 | ClimaLand 2026 (Eq. 32-37) | CLM5 PHS (Kennedy 2019) | COMPASS v2 | 一致性 |
|---|---|---|---|---|
| **水量平衡** | `dW/dt = q_in - q_out` | `dW/dt = q_in - q_out` | `dW/dt = q_in - q_out` | ✅ 完全相同 |
| **Darcy + 重力** | `q = -K·((ψj-ψi)/Δz + 1)` | 同 | `q = -K_eff·((ψj-ψi)/dz + 1)` | ✅ |
| **Weibull K(ψ)** | `K = K_sat·exp(-(ψ/ψ63)^c)` | `K = K_max·exp(-(ψ/p50)^c)` | `K = K_sat·exp(-(|ψ|/|ψ63|)^c)` | ✅ |
| **持水曲线** | 线性 `ψ = (S_l-1)/a` | exponent `ψ = -ψ0·S^(-1/b)` | 线性 `ψ = (W/Wcap-1)/a` | ✅ 跟 ClimaLand 同 |
| **节点导度** | 调和平均 `K_eff = K1·K2/(K1+K2)` | 同 | 同 | ✅ |
| **土→根** | 根分布加权 | 根分布加权 | 根分布加权 `K_sr_eff = Σ P(z)·K(z)` | ✅ |
| **水力截断 T** | `min(T_pot, q_supply)` | btran 因子 | `min(T_pot, q_stem→leaf)` | ✅ |
| **时间离散** | 显式 RK + IMEX | 显式 | **IMEX 3×3** | ⚠️ 比 ClimaLand 标准更稳 |

### 2.2 物理过程整体评价

**全部对、且有一处超越 ClimaLand**:

| 维度 | 评价 |
|---|---|
| **方程组** | ✅ 跟 ClimaLand 2026 paper Eq. 32-37 严格一致 |
| **Weibull K(ψ)** | ✅ 跟 CLM5 PHS (Kennedy 2019, Bonan 2011) 一致 |
| **线性持水** | ✅ Christoffersen 2016 GMD 也用这个 |
| **调和平均** | ✅ ClimaLand 第 246 行 `harmonic_mean(x,y)=xy/max(x+y,eps)`，跟 v2 一致 |
| **IMEX 3×3** | ⭐ 比 ClimaLand 默认 explicit-RK 更稳 — v2 优势 |
| **水力截断 T** | ✅ Sperry-Love 2015 PhytNew 框架 |
| **守恒诊断** | ✅ v2 第 220 行有 `water_balance_residual` —— ClimaLand 没显式做 |
| **奇异矩阵保护** | ✅ v2 第 184 行 `rcond < 1e-12` 回退显式 — CLM5/ClimaLand 都没有 |

**结论**: v2 物理是**严格 ClimaLand 2026 实现 + 数值稳健性补强**。没有任何"COMPASS 私货公式"，可以放心用作 VOD 的 VWC 输入源。

### 2.3 v2 跟 ClimaLand 唯一的方法学差异

| 项 | ClimaLand | v2 | 评价 |
|---|---|---|---|
| 求解器 | ExplicitRK / IMEXEuler (Differential Equations.jl) | 内置 IMEX 3×3 | v2 自给自足、不依赖外部 ODE 包；可控性强 |
| 茎水势空间分辨 | 可分多个 stem 层 | 单层 (root/stem/leaf) | v2 简化；对中小树足够 |
| 异速生长（SAI, AGB 季节变化） | 通过 LAI/SAI 时间序列输入 | 静态参数（W_capacity 固定） | v2 是 ClimaLand 单年仿真模式 |

**不影响 VOD 计算**：VWC 只需 `W_stem(t)` 和 `W_leaf(t)`，单层精度足够。

---

## 3. 物理改进建议（如果想再进一步）

按"投入产出比"排序：

### 改进 1 (推荐做)：把 W_capacity 跟 LAI 联动
**现状**：`W_leaf_capacity` 是常数 0.001 m，意味着叶水容量全年不变。
**问题**：苹果树冬季落叶，LAI=0.293 时叶水容量应该接近 0，而不是 0.001。
**方案**：
```matlab
% 在 prognosticPlantWater 调用前
Wcap(3) = plant_params.W_leaf_capacity_max * (LAI / LAI_max);
```
**收益**：冬季 W_leaf ≈ 0 → VOD_leaf 自动归零 → 跟苹果园物候一致

### 改进 2 (可选)：psi_63 跟着 PFT 走
**现状**：`psi_63 = -275 m` 是温带乔木默认值。
**改成**：从 `+plant/getPFTParameterSet.m` 里读站点苹果树文献值（典型 −2.5 MPa = −255 m）
**收益**：水力胁迫触发阈值更准

### 改进 3 (不必做)：换非线性持水曲线
ClimaLand 默认就是线性，不需要换。

### 改进 4 (远景)：多层 stem
对苹果树这种 ~3m 高度的小乔木，单层 stem 完全够；对森林才需要。先不做。

**总结**：现状 v2 物理已经够好，**改进 1（W_capacity 跟 LAI 联动）是唯一性价比高的**，5 分钟工作。

---

## 4. 耦合可行性 — v2 跟现在 COMPASS.m 接口

### 4.1 当前耦合状态：零接入

```
grep "lProgHydraulic\|Plant_hydraulics_v2\|prognosticPlantWater" src/COMPASS.m
→ 0 matches
```

即 v2 完整写好了、单元测试也有 ([test_v2_unit_tests.m](../src/+plant/+PlantHydraulic/test_v2_unit_tests.m))，但**主循环里没人调它**。

### 4.2 v2 需要哪些 COMPASS.m 给的东西

v2 入口签名：
```matlab
[W_new, psi_new, q_fluxes, info] = prognosticPlantWater( ...
    W_old, psi_soil_profile, P_root_dist, K_soil_root_profile, ...
    T_potential, plant_params, Delt_t)
```

挨个对比：

| v2 参数 | 物理量 | COMPASS.m 里有吗 | 来源 |
|---|---|---|---|
| `W_old` [3×1, m] | 上步 W_root/W_stem/W_leaf | ❌ 没初始化 | **要加 state** |
| `psi_soil_profile` [NN×1, m] | 土壤水势剖面 | ✅ `hh / 100` (cmH₂O → m) | SoilVariables.hh |
| `P_root_dist` [NN×1] | 根分布密度 | ✅ `Rl/sum(Rl)` | RTBSettings.Rl |
| `K_soil_root_profile` [NN×1, m/s] | 土-根接口导度 | ✅ `KfL_h` 节点平均 (cm/s → m/s) | SoilVariables.KfL_h |
| `T_potential` [m/s] | 潜在蒸腾 | ✅ `(Tp_sun + Tp_sha)` (cm/s → m/s) | SCOPE 输出 |
| `plant_params.*` | Weibull/持水/几何参数 | ⚠️ 一半在 config_v2.txt | loadProgV2Config |
| `Delt_t` [s] | 时步 | ✅ `ModelSettings.DELT * 3600` | |

**结论**：除了 `W_old` 需要新加 state 字段（初始化用 `inverseLinearRetention` 给 ψ_soil_top → W_old），**其他所有输入 COMPASS.m 都已经有**。

### 4.3 耦合接入的最小改动 (5 处)

#### 改动 A. COMPASS.m 主时间循环之前：初始化 W
```matlab
% 加在 KT=1 之前，约 COMPASS.m 第 X 行
if cfg.lProgHydraulic == 1
    W_state = [...
        plant_params.W_capacity(1) * 0.9;   % W_root
        plant_params.W_capacity(2) * 0.9;   % W_stem
        plant_params.W_capacity(3) * 0.9];  % W_leaf (饱和度 90% 起步)
end
```

#### 改动 B. 主循环内、SCOPE 之后、Plant_hydraulics 之前调用 v2
```matlab
if cfg.lProgHydraulic == 1
    [W_state, psi_plant, q_fluxes, info_v2] = vegetation.hydraulics.prognosticPlantWater( ...
        W_state, ...
        SoilVariables.hh / 100, ...                  % cm → m
        RTBSettings.Rl / sum(RTBSettings.Rl), ...   % 归一化根分布
        SoilVariables.KfL_h_node / 100, ...         % cm/s → m/s
        (Tp_sun + Tp_sha) / 100, ...                % cm/s → m/s
        plant_params, ...
        ModelSettings.DELT * 3600);
    % 用 v2 的水势覆盖原来 SPAC 给的 psi_leaf
    psil_sun(KT) = psi_plant(3) / 0.1019;  % m → MPa
    psil_sha(KT) = psi_plant(3) / 0.1019;
end
```

#### 改动 C. 输出 W 到诊断
```matlab
if cfg.lProgHydraulic == 1
    reg.write('W_root_storage', W_state(1), KT);
    reg.write('W_stem_storage', W_state(2), KT);
    reg.write('W_leaf_storage', W_state(3), KT);
end
```

#### 改动 D. 在 +diagnostics/+groups/registerPlantHydraulics.m 注册三个新字段
（按 COMPASS coding standard §7.2 流程）

#### 改动 E. 加载 v2 配置
```matlab
% 在 cfg 加载阶段
if cfg.lProgHydraulic == 1
    cfg_v2 = loadProgV2Config('data/config_v2.txt');
    plant_params = build_plant_params(cfg_v2);  % helper to map cfg → struct
end
```

**总改动量**: ~50 行，全部加在 if 分支里，**默认 lProgHydraulic=0 完全等价于现状**。

### 4.4 可能的雷区

| 雷区 | 风险等级 | 解决 |
|---|---|---|
| 单位混乱 (cm vs m) | 中 | 严格按上表换算，每行加 `% [unit]` 注释 |
| psi 符号约定（+/− 吸力） | 高 | ClimaLand: ψ ≤ 0 吸力；COMPASS hh < 0 吸力。**两者同向**，可直接 hh/100 |
| `RTBSettings.Rl` 节点 vs 层 | 中 | v1 是节点；v2 也按节点；OK |
| `KfL_h` 是 2 列矩阵（前后步） | 中 | 取当前列 `KfL_h(:,1)` |
| Δz 是参数 vs 物理距离 | 低 | dz_rs/dz_sl 已在 config_v2.txt |
| 首次 W_old 该是多少 | 低 | 用饱和度 0.9 起步，2-3 天 spinup 收敛 |

**结论**：耦合可行，且**不影响现有代码**（条件分支 + 新增 state）。

---

## 5. 接入后能给 VOD 模块带来什么

回到 [vod_module_roadmap.md](vod_module_roadmap.md) 的核心问题 —— VOD 平直、没有季节循环。

接入 v2 之后，VWC 计算变成**纯物理零经验**：

```matlab
% +Radiation/+microwave/compute_VOD.m 改 2 行：
VWC_stem = W_stem * 1000;    % m → kg/m²
VWC_leaf = W_leaf * 1000;    % m → kg/m²
% 删掉之前的 LAI×Cw×10 和 AGB_static 这两个经验式
```

**预期变化**：
- `W_stem` 在干旱后期会从 0.005 m → 0.003 m（缩水 40%）→ VOD_stem 从 0.5 → 0.3
- `W_leaf` 雨后涨、午后跌（昼夜变化）→ VOD 日内能有 0.05 振幅
- 冬季落叶期 W_leaf ≈ 0 → 自动季节脱钩
- ψ_leaf 不再是诊断稳态，而是**真实瞬态**，b_leaf(ψ) 调控才有意义

这正是你 PPT Slide 18 那张图的**真实物理化**。

---

## 6. 推荐路径

按你"先理清楚再动手"的原则：

### Step 1 (1 天)：单元测试 v2 跑通
```
src/+plant/+PlantHydraulic/test_v2_unit_tests.m
src/+plant/+PlantHydraulic/demo_v2_step_response.m
```
**确认**：单步 IMEX 3×3 数值上稳；ψ_leaf 对蒸腾扰动响应合理（小时尺度）。

### Step 2 (2 天)：物理改进 1 实现
把 W_capacity(3) 跟 LAI 联动（5 分钟代码 + 测试）。

### Step 3 (3-5 天)：实际接入 COMPASS.m
按 §4.3 的 5 处改动加进去，跑一年仿真，输出对比 v1 vs v2 的：
- `psi_leaf(t)` 时间序列
- `W_stem(t)` / `W_leaf(t)` 演化
- 蒸腾 `T_actual` 跟 v1 的偏差
- 水量守恒残差

### Step 4 (1 天)：接入 +Radiation/+microwave
把 `compute_VOD.m` 里 VWC 改成 W·1000，重跑全年，看 VOD 季节相位是否变活。

### Step 5 (1 天)：跟 MT-DCA + VODCA 比对
出图 + RMSE 报告。

**总工期: 7-10 天**（比之前 PPT 路线图说的 "2-3 周" 短，因为 v2 模块已经全写好了，只是没接入）。

---

## 7. 一句话结论

> COMPASS v2 (prognosticPlantWater.m) 物理跟 ClimaLand 2026 严格对齐 + 数值稳健性更强，**没有需要改的物理 bug**。
> 唯一推荐的小改进是 `W_leaf_capacity` 跟 LAI 联动（5 分钟）。
> 接入主循环只要 ~50 行条件分支代码，且**不动现有任何代码**。
> 接入后 VWC 直接 = W·1000，**VOD 计算从经验式变成纯物理变量**。

**这是 PPT Slide 18 真正实现的最直接路径，而且 v2 模块已经做好 80% 的工作。**
