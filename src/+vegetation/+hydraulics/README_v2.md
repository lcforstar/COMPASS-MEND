# Prognostic Canopy Energy & Plant Hydraulics (v2)

平行版（v2）实现，**不修改任何现有代码**。通过 `cfg.lProgEnergy` 和 `cfg.lProgHydraulic` 开关启用。

---

## 文件清单

| 文件 | 类型 | 职责 |
|---|---|---|
| `weibullConductance.m`     | 纯函数 | Weibull 木质部导度 K(ψ) |
| `linearRetention.m`        | 纯函数 | 线性持水曲线 ψ(W) |
| `inverseLinearRetention.m` | 纯函数 | 反函数 W(ψ)，初始化用 |
| `harmonicMean.m`           | 纯函数 | 调和平均（节点间 K_eff） |
| `prognosticCanopyEnergy.m` | 算法 | #6 冠层瞬态能量平衡（IMEX 标量） |
| `prognosticPlantWater.m`   | 算法 | #7 植物水量平衡（IMEX 3×3） |
| `Plant_hydraulics_v2.m`    | 接口 | 替换 `Plant_hydraulics.m` 的瞬态版（待写）|
| `tests_v2/`                | 测试 | 单元测试（待写）|

---

## 物理依据

| 模块 | 公式来源 |
|---|---|
| #6 冠层热容 | ClimaLand 2026 Eq. 39 |
| #7 植物水力 | ClimaLand 2026 Eq. 32–37 |
| Weibull K(ψ) | ClimaLand 2026 Eq. 34 |
| 线性持水曲线 | ClimaLand 2026 Eq. 35 |

---

## 数学规范

### #6 冠层瞬态能量

```
(ρch) · A_canopy · dT/dt = R_n - LE - H + S_re_e
```

IMEX 半隐式离散：
```
A_imp · ΔT = R_explicit
A_imp = (ρch·A)/Δt + ρa·cp/r_ah + (ρa·Lv/r_ae)·(Lv·qsat)/(Rv·T²) + 8εσT³
R_explicit = R_n - H - LE + S_re_e
```

### #7 植物水力

```
dW_root / dt = q_sr - q_rs
dW_stem / dt = q_rs - q_sl
dW_leaf / dt = q_sl - T_actual

q_ij = -K_eff · ((ψj - ψi)/Δz + 1)        Darcy + gravity
K(ψ) = K_sat · exp(-(|ψ|/|ψ63|)^c)        Weibull
ψ(W) = (W/W_cap - 1) / a                  Linear retention
```

IMEX 3×3 隐式系统，每步 O(1) 求逆。

---

## 调用示例

### #6 冠层能量

```matlab
[T_new, dT, info] = vegetation.hydraulics.prognosticCanopyEnergy( ...
    T_old, T_iter, Rn, H, LE, S_re_e, ...
    rho_c_h, A_canopy, ...
    rho_a, cp, Lv, r_ah, r_ae, ...
    qsat_curr, T_air, ...
    epsilon_canopy, sigmaSB, R_v, ...
    Delt_t);

% 守恒检查
disp(info.EB_residual)   % 应该接近 0
disp(info.tau_canopy)    % 热惯性时标 (s)
```

### #7 植物水力

```matlab
plant_params = struct( ...
    'K_sat',       2e-8, ...      % 木质部饱和导水率 m/s
    'psi_63',      -275, ...      % 63%损失水势 m
    'c_weibull',   4, ...
    'a_retention', [4.9e-4; 4.9e-4; 4.9e-4], ...   % [root;stem;leaf]
    'W_capacity',  [0.003; 0.005; 0.001], ...      % 各层容量 m
    'dz_rs',       0.5, ...       % 根→茎距离
    'dz_sl',       1.0);          % 茎→叶距离

[W_new, psi_new, q_fluxes, info] = vegetation.hydraulics.prognosticPlantWater( ...
    W_old, psi_soil_profile, P_root_dist, K_soil_root_profile, ...
    T_potential, plant_params, Delt_t);

% 守恒检查
disp(info.water_balance_residual)   % 应该 < 1e-6 m
```

---

## 默认参数（苹果园）

| 参数 | 值 | 来源 |
|---|---|---|
| ρch (无果) | 745 J/m²/K | ClimaLand US-Var |
| ρch (含果) | ~63000 J/m²/K | 苹果园估算 |
| K_sat 木质部 | 2e-8 m/s | ClimaLand 默认 |
| ψ63 | -275 m (-2.7 MPa) | ClimaLand 默认 |
| Weibull c | 4 | ClimaLand 默认 |
| a_retention | 4.9e-4 m⁻¹ | ClimaLand 默认 |

---

## 验证测试计划

1. **守恒测试**：跑 1 周，能量误差 < 0.1%，水量误差 < 1mm
2. **稳态退化**：Δt → ∞，结果应收敛到原 `Plant_hydraulics` 输出
3. **时标测试**：阶跃辐射 → τ_canopy；阶跃 VPD → τ_plant
4. **苹果园对比**：vs 通量塔实测（冠层温度日变化、中午水势）

---

## 论文 Claim

> **联合论文**：*"Prognostic canopy energy and plant hydraulics in STEMMUS-SCOPE: implications for sub-hourly to multi-day water-energy coupling"*
>
> 1. 冠层瞬态能量减小亚小时温度模拟 RMSE
> 2. 瞬态植物水力捕捉昼夜水势节律和干旱滞后
> 3. 保留多层 SCOPE + 自适应步长 + 三相耦合（vs ClimaLand 简化）
>
> 目标：JAMES / GMD

---

## 状态：全部完成 ✅

### 核心算法 (4 文件)
- [x] `weibullConductance.m`              Weibull K(ψ)
- [x] `linearRetention.m`                 线性持水曲线 ψ(W)
- [x] `inverseLinearRetention.m`          反函数 W(ψ)
- [x] `harmonicMean.m`                    节点 K_eff 调和平均

### IMEX 求解器 (2 文件)
- [x] `prognosticCanopyEnergy.m`          #6 IMEX 标量
- [x] `prognosticPlantWater.m`            #7 IMEX 3×3

### 高层接口 (1 文件)
- [x] `Plant_hydraulics_v2.m`             兼容 v1 接口的瞬态版

### 测试 & 演示 (2 文件)
- [x] `test_v2_unit_tests.m`              6 个单元测试
- [x] `demo_v2_step_response.m`           4 个物理 demo

### 配套（在 src/ 和 data/）
- [x] `src/ebal_clm45_v2.m`               冠层后处理包装
- [x] `src/AMainloop_v2_demo.m`           v2 演示驱动器
- [x] `src/compare_v1_v2.m`               v1 vs v2 对比脚本
- [x] `src/loadProgV2Config.m`            配置加载
- [x] `data/config_v2.txt`                v2 专用配置
- [x] `src/PROGNOSTIC_V2_QUICKSTART.md`   集成指南

## 快速验证

```matlab
cd 'f:/STEMMUS_Zhujizixing -苹果17 - 变步长/src'

% 1. 跑 6 个单元测试
test_v2_unit_tests

% 2. 跑 4 个物理 demo
demo_v2_step_response

% 3. 跑 v1 vs v2 对比表
compare_v1_v2

% 4. 演示完整时步推进
results = AMainloop_v2_demo();
```

最后更新: 2026-05-02
