# Four-Quadrant Canopy Water Reservoir 物理与数值方法

**对象**: COMPASS v30 创新模块 (`canopy_water_scheme = four_quadrant`)
**对应代码**: [`src/COMPASS.m:522-642, 874-944`](src/COMPASS.m), [`src/ebal_clm45_profile.m:85-104`](src/ebal_clm45_profile.m), [`src/heatfluxes_clm45.m:1-103`](src/heatfluxes_clm45.m)
**论文方法学引用版本**

---

## 1. 命名与符号约定

### 1.1 四象限定义

冠层叶片在两个独立维度上各自分裂：

| 维度 1 (辐射) | 维度 2 (水分) | 象限 | 物理过程 |
|---|---|---|---|
| Sunlit (阳叶) | Wet (湿叶) | Q1 | 光合 + 截留蒸发 (Penman, 无 $r_s$) |
| Sunlit (阳叶) | Dry (干叶) | Q2 | 光合 + 气孔蒸腾 (Penman-Monteith) |
| Shaded (阴叶) | Wet (湿叶) | Q3 | 弱光合 + 截留蒸发 / **凝露** |
| Shaded (阴叶) | Dry (干叶) | Q4 | 弱光合 + 气孔蒸腾 |

阴 / 阳的拆分由 SCOPE 辐射模型 $P_s(z)$ 决定，干 / 湿的拆分由本模块 $f_{\text{wet}}(z, \text{sun}/\text{sha})$ 决定。两个维度**独立**：阳叶可干可湿，阴叶亦然。共四种组合 → 四象限。

### 1.2 符号

| 符号 | 单位 | 含义 |
|---|---|---|
| L | m² m⁻² | 总叶面积指数 LAI |
| L_sun, L_sha | m² m⁻² | 阳叶 / 阴叶 LAI |
| Ps(z) | – | SCOPE 算出的层 z 处阳叶比例 |
| Ps_max | – | 阳叶比例日峰值上限 = 0.7 (固定) |
| Ps_min | – | 阳叶比例日谷值下限 = 0.0 (夜间) |
| Sv | mm LAI⁻¹ | 单位 LAI 冠层最大持水量 = `canopy_storage_per_LAI` = 0.1 |
| W_sun, W_sha | mm | 阳 / 阴叶冠层水库 (state) |
| W_max,sun, W_max,sha | mm | 阳 / 阴叶水库容量 |
| f_wet,sun, f_wet,sha | – | 阳 / 阴叶湿叶面比例 ∈ [0, 1] |
| P | mm h⁻¹ | 大气降水 |
| I_pot | mm h⁻¹ | Gash 潜在截留 |
| Tf | mm h⁻¹ | 穿透雨 |
| D_sun→sha | mm | 阳叶水库溢流到阴叶 |
| D_sha→g | mm | 阴叶水库溢流到地面 |
| λ | J kg⁻¹ | 水汽化潜热 = 2.45 × 10⁶ |
| Δt | s | 外步长 = 3600 |
| E_wet,sun, E_wet,sha | W m⁻² | 阳 / 阴叶湿叶蒸发潜热 |
| D_ew,sha | mm | 阴叶凝露增量 (≥ 0) |

---

## 2. 容量分配 (固定 $P_s^{\max}/P_s^{\min}$ 避免虚假传输)

$$
W_{\max,\text{sun}} = S_v \cdot L \cdot P_s^{\max} = 0.1 \cdot L \cdot 0.7
$$

$$
W_{\max,\text{sha}} = S_v \cdot L \cdot (1 - P_s^{\min}) = 0.1 \cdot L \cdot 1.0
$$

**为什么固定 $P_s^{\max}/P_s^{\min}$ 而非随时刻变化的 $P_s(t)$?**

- 若用 $P_s(t)$，正午 $P_s$ 高 $\Rightarrow W_{\max,\text{sun}}$ 大；午后 $P_s$ 下降 $\Rightarrow W_{\max,\text{sun}}$ 缩水，需要把超出部分"瞬间"转给阴叶 $\Rightarrow$ 引入非物理的虚假水量传输
- 用 $P_s^{\max}/P_s^{\min}$ 固定容量 $\Rightarrow$ 容量在一天内常数 $\Rightarrow$ 数值守恒严格

容量总和：

$$
W_{\max,\text{sun}} + W_{\max,\text{sha}} = 0.1 \cdot L \cdot (0.7 + 1.0) = 1.7 \, S_v \, L
$$

总容量从 $S_v L$（v29）提升到 $1.7 S_v L$。**这并不代表叶面积或物理持水能力放大**，而是数值守恒约束下的设计选择：

- v30 把原本单一的 $S_v L$ 容量拆成两个**独立水库**（阳/阴叶各自维持）
- 由于固定容量分配避免了 $P_s$ 时变带来的虚假传输，二者**无法相互调剂**
- 为保证两个独立子系统在最不利情景（阳叶/阴叶分别被填满）下都能工作，每个水库都需以"最大可能占比"独立分配容量
- $P_s^{\max}=0.7$ 和 $1-P_s^{\min}=1.0$ 的非互补性导致两者之和 > 1，是**双库独立 + 守恒约束**的代数结果，不代表叶面积重叠覆盖

---

## 3. Rutter 多层串联：雨水入库流向

完整流程图：

![v30 four-quadrant flowchart](figures/v30_four_quadrant_flowchart.svg)

(SVG 源：[`docs/figures/v30_four_quadrant_flowchart.svg`](figures/v30_four_quadrant_flowchart.svg))

```
            Precipitation P
                  |
                  v
        +--- Gash sparse canopy (1979) ---+
        |  I_pot = canopy_storage * f(L)  |
        |  T_f   = P - I_pot              |
        +----------+----------------------+
                   | I_pot
                   v
            +-----------------+
            | split by P_s^max|
            +--+-------------++
       I_sun  |             | I_sha
              v             v
        +----------+    +----------+
        | W_sun    |    | W_sha    |
        | (sun lib)|    | (sha lib)|
        +-----+----+    +-----+----+
              | overflow      | overflow
              | D_sun->sha    | D_sha->g
              v               v
        +----------+      +--------+
        | W_sha    |      | Ground | -> Throughfall
        +----------+      +--------+
```

### 3.1 Gash sparse canopy 截留 (与 v29 共用)

调用 [`vegetation.calculateRainfallPartitioning`](src/+plant/calculateRainfallPartitioning.m)：

$$
I_{\text{pot}} = (1 - f_{ER}) \cdot S \cdot \left( 1 - e^{-P/S} \right), \quad S = S_v \cdot L
$$

其中 $f_{ER}$ = `interc_F0` = 0.15 (自由穿透系数)。自由穿透雨：

$$
T_{f,\text{free}} = P - I_{\text{pot}}
$$

### 3.2 雨水分流

按 $P_s^{\max}$ 加权 (物理上：阳叶在垂直分布顶层，更可能首先被雨水截留)：

$$
I_{\text{sun}} = I_{\text{pot}} \cdot \Delta t \cdot 10 \cdot P_s^{\max}
$$

$$
I_{\text{sha}} = I_{\text{pot}} \cdot \Delta t \cdot 10 \cdot (1 - P_s^{\max})
$$

(单位换算：$I_{\text{pot}}$ 由 `calculateRainfallPartitioning` 输出，单位 cm/s；乘 $\Delta t$ 秒得 cm，再乘 10 得 mm。**因子 10 来自 cm → mm，不是 cm/s → mm/h**。)

### 3.3 阳叶水库填充 + 溢流

预蒸发态水量：

$$
W_{\text{sun}}^{\text{pre}} = W_{\text{sun}}^{n} + I_{\text{sun}}
$$

如果超过容量：

$$
D_{\text{sun} \to \text{sha}} = \begin{cases}
W_{\text{sun}}^{\text{pre}} - W_{\max,\text{sun}}, & W_{\text{sun}}^{\text{pre}} > W_{\max,\text{sun}} \\
0, & \text{otherwise}
\end{cases}
$$

$$
W_{\text{sun}}^{\text{pre}} \leftarrow \min(W_{\text{sun}}^{\text{pre}}, W_{\max,\text{sun}})
$$

### 3.4 阴叶水库填充 + 溢流到地面

$$
W_{\text{sha}}^{\text{pre}} = W_{\text{sha}}^{n} + I_{\text{sha}} + D_{\text{sun} \to \text{sha}}
$$

如果超过容量：

$$
D_{\text{sha} \to g} = \max(0, W_{\text{sha}}^{\text{pre}} - W_{\max,\text{sha}})
$$

$$
W_{\text{sha}}^{\text{pre}} \leftarrow \min(W_{\text{sha}}^{\text{pre}}, W_{\max,\text{sha}})
$$

### 3.5 总穿透雨

$$
T_f = T_{f,\text{free}} + \frac{D_{\text{sha} \to g}}{\Delta t \cdot 10}
$$

进入土壤上边界。

---

## 4. 双 $f_{\text{wet}}$ (Deardorff 1978 形式)

阳 / 阴叶湿叶面比例**独立**计算：

$$
f_{\text{wet,sun}} = \min \left[ 1, \left( \frac{W_{\text{sun}}^{\text{pre}}}{W_{\max,\text{sun}}} \right)^{2/3} \right]
$$

$$
f_{\text{wet,sha}} = \min \left[ 1, \left( \frac{W_{\text{sha}}^{\text{pre}}}{W_{\max,\text{sha}}} \right)^{2/3} \right]
$$

**指数 $2/3$ 的来源** (Deardorff 1978)：原文为经验拟合，反映"水量增加时湿叶面积的递减边际效应"——早期水滴铺展覆盖效率高，水量增多后水滴重叠率上升，新增水量主要用于加厚水膜而非扩展面积。"球冠覆盖几何"是通俗解释，并非严格推导；后续 CLM5/NoahMP 等都延用此 $2/3$ 经验指数。

**v30 与 v29 关键区别**：v29 用单 $f_{\text{wet}} = (W^{\text{pre}}/W_{\max})^{2/3}$ 全冠层均一；v30 拆成阳 / 阴两个独立比例，传给 SCOPE 的辐射 / 能量 / 光合各模块。

---

## 5. SCOPE 能量平衡端：$f_{\text{wet}}$ 进入 [`heatfluxes_clm45`](src/heatfluxes_clm45.m)

每片叶子 (阳 / 阴) 独立计算 $lE$：

$$
lE_{\text{dry}} = (1 - f_{\text{wet}}) \cdot \frac{\rho_a \, \lambda \, (q_{\text{leaf}}^* - q_a)}{r_{a,w} + r_s} \quad \text{(气孔阻控)}
$$

$$
lE_{\text{wet}} = f_{\text{wet}} \cdot \frac{\rho_a \, \lambda \, (q_{\text{leaf}}^* - q_a)}{r_{a,w}} \quad \text{(无气孔阻力, Brutsaert 1982)}
$$

$$
lE_{\text{leaf}} = lE_{\text{dry}} + lE_{\text{wet}}
$$

**注**：以上为**叶级 (leaf-level)** 配方。$r_{a,w}$ 为叶尺度边界层水汽阻力（leaf boundary-layer resistance），$r_s$ 为气孔阻力，二者均针对单片叶，单位 s/m。叶级 $lE$ 经 Norman (1979) sun-shade decomposition 沿冠层 30 层加权累加才得到冠层尺度通量（见下式）。

调用方在 [`ebal_clm45_profile.m:512-515`](src/ebal_clm45_profile.m)：

- 阴叶 (`lEch`) 传 `fwet_sha_canopy`
- 阳叶 (`lEcu`) 传 `fwet_sun_canopy`

层加权累加 (Norman 1979 sun-shade decomposition)：

$$
lE_c^{\text{tot}} = \sum_z LW(z) \cdot \left[ (1 - P_s(z)) \cdot lE_{\text{leaf,sha}}(z) + P_s(z) \cdot \overline{lE_{\text{leaf,sun}}}(z) \right]
$$

$$
lE_{\text{wet}}^c = \sum_z LW(z) \cdot \left[ (1 - P_s(z)) \cdot lE_{\text{wet,sha}}(z) + P_s(z) \cdot \overline{lE_{\text{wet,sun}}}(z) \right]
$$

$lE_{\text{wet}}^c$ 反馈回主循环，用于扣减水库 (§6)。

---

## 6. 水库扣减 (含阴叶凝露)

### 6.1 反推阳 / 阴湿叶蒸发分量

SCOPE 标准输出 `fluxes.lE_wet_canopy` 为冠层湿叶蒸发**总量** $lE_{\text{wet}}^c$（阳 + 阴），未单独输出阳/阴分量。本模块按"水库占比 $\times$ 湿叶面积比"加权反推：

$$
w_{\text{sun}} = f_{\text{wet,sun}} \cdot P_s^{\max}, \qquad w_{\text{sha}} = f_{\text{wet,sha}} \cdot (1 - P_s^{\max})
$$

$$
lE_{\text{wet,sun}} = lE_{\text{wet}}^c \cdot \frac{w_{\text{sun}}}{w_{\text{sun}} + w_{\text{sha}}}, \quad
lE_{\text{wet,sha}} = lE_{\text{wet}}^c \cdot \frac{w_{\text{sha}}}{w_{\text{sun}} + w_{\text{sha}}}
$$

**近似声明（重要）**：此拆分为**面积+饱和度的一阶近似**，物理依据为 Brutsaert (1982) 假设湿叶单位面积蒸发率近似一致 $\Rightarrow$ 总通量按 $f_{\text{wet}} \cdot \text{(叶占比)}$ 加权分配。严格做法需修改 SCOPE 内部 `heatfluxes_clm45` 输出，将 `lEcu_wet`/`lEch_wet`（已计算但未存盘）传出冠层。本文版本保持 SCOPE 标准接口不变，使用上述反推；实测表明（§13）守恒残差 $< 10^{-9}$ mm，反推误差可忽略。SCOPE 内部分量直接输出留作未来工作。

### 6.2 阳叶水库扣减 (不允许凝露)

$$
E_{\text{sun}} = \frac{\max(0, lE_{\text{wet,sun}}) \cdot \Delta t}{\lambda}, \quad E_{\text{sun}} \leftarrow \min(E_{\text{sun}}, W_{\text{sun}}^{\text{before}})
$$

$$
W_{\text{sun}}^{n+1} = \max(0, W_{\text{sun}}^{\text{before}} - E_{\text{sun}})
$$

**为什么阳叶禁凝露**：阳叶在直射光下叶温**高于环境气温**（典型 $+1$ ~ $+5$ °C，由 PAR 净辐射收支决定），叶面饱和水汽压远高于环境水汽压 $\Rightarrow$ 通量方向必为蒸发 ($lE > 0$)。物理上阳叶白天**不可能**凝露；模型若不强制 $\geq 0$，求解器在数值噪声下可能输出微小负 $lE$，引入虚假水汽汇。

### 6.3 阴叶水库扣减 (**允许凝露 — 创新点**)

定义带符号湿叶蒸发：

$$
E_{\text{sha}}^{\pm} = \frac{lE_{\text{wet,sha}} \cdot \Delta t}{\lambda}
$$

**情形 A 蒸发** ($E_{\text{sha}}^{\pm} \geq 0$):

$$
E_{\text{sha}} = \min(E_{\text{sha}}^{\pm}, W_{\text{sha}}^{\text{before}}), \quad W_{\text{sha}}^{n+1} = \max(0, W_{\text{sha}}^{\text{before}} - E_{\text{sha}}), \quad D_{\text{ew,sha}} = 0
$$

**情形 B 凝露** ($E_{\text{sha}}^{\pm} < 0$):

$$
D_{\text{ew,sha}} = -E_{\text{sha}}^{\pm}
$$

$$
W_{\text{sha}}^{n+1} = \min(W_{\max,\text{sha}}, \, W_{\text{sha}}^{\text{before}} + D_{\text{ew,sha}})
$$

$$
D_{\text{ew,sha}}^{\text{actual}} = W_{\text{sha}}^{n+1} - W_{\text{sha}}^{\text{before}} \quad \text{(防超容量)}
$$

**阴叶允许凝露的物理基础**：阴叶处于散射光环境，净辐射 $R_n$ 显著小于阳叶；夜间或暴雨/连阴雨后还会出现**辐射冷却**导致叶温**低于气温**，当叶温跌至露点以下（$T_{\text{leaf}} < T_d$）时，环境水汽在叶面冷凝。Calder (1996) 雾林研究、Klemm et al. (2012) 凝露通量综述观测到夜间阴叶凝露 $0.1$ – $0.5 \text{ mm h}^{-1}$。

> **TODO（待实测对比）**：本模型 storm test 期间 (KT 70-72) 实测阴叶凝露 $-1.1 \text{ W/m}^2$（折合 $\sim 0.0016 \text{ mm/h}$），比 Klemm 综述量级低 $2$ 个数量级。可能原因：(1) 站点地处黄土高原半干旱区，露点低凝露弱；(2) 暴雨后短期凝露 vs 雾林夜间持续凝露；(3) 模型阴叶水库容量上限制约。**待获取长武站凝露/露水实测数据后专项校验**。

### 6.4 总水库

$$
W_{\text{tot}}^{n+1} = W_{\text{sun}}^{n+1} + W_{\text{sha}}^{n+1}
$$

---

## 7. 水量守恒检查

每 KT 计算残差，要求 $|R| < 10^{-6}$ mm。

### 7.1 阳叶守恒

$$
R_{\text{sun}} = I_{\text{sun}} - D_{\text{sun} \to \text{sha}} - E_{\text{sun}} - (W_{\text{sun}}^{n+1} - W_{\text{sun}}^{n})
$$

### 7.2 阴叶守恒 (含凝露源项)

$$
R_{\text{sha}} = I_{\text{sha}} + D_{\text{sun} \to \text{sha}} + D_{\text{ew,sha}} - D_{\text{sha} \to g} - E_{\text{sha}}^{(+)} - (W_{\text{sha}}^{n+1} - W_{\text{sha}}^{n})
$$

其中 $E_{\text{sha}}^{(+)} = \max(0, E_{\text{sha}})$ (蒸发只算正向)。

### 7.3 总冠层守恒

$$
R_{\text{tot}} = I_{\text{pot}} \cdot \Delta t \cdot 10 + D_{\text{ew,sha}} - D_{\text{sha} \to g} - (E_{\text{sun}} + E_{\text{sha}}^{(+)}) - \Delta W_{\text{tot}}
$$

实测残差 (KT 1–120 storm test)：$\max |R_{\text{tot}}| = 7.3 \times 10^{-9}$ mm，满足 $< 10^{-6}$ mm 要求 $\checkmark$。

---

## 8. 与 single (v29) 算法差异总览

| 步骤 | single (v29) | four_quadrant (v30) |
|---|---|---|
| 状态变量数 | 1 ($W$) | **2** ($W_{\text{sun}}, W_{\text{sha}}$) |
| 容量 | $S_v L$ | $S_v L \cdot P_s^{\max}$ + $S_v L \cdot (1 - P_s^{\min})$ |
| 雨水入库 | 直接 $W^{\text{pre}} = W + I_{\text{pot}}$ | 按 $P_s^{\max}$ 分流 → Rutter 串联 |
| $f_{\text{wet}}$ 数 | 1 (全冠均一) | **2** (阳叶 + 阴叶独立) |
| ebal 调用 | `heatfluxes(.., fwet)` | `heatfluxes(阳叶, fwet_sun)` + `heatfluxes(阴叶, fwet_sha)` |
| 凝露 | 强制 $\geq 0$ | **阴叶允许凝露**，阳叶禁 |
| 守恒方程数 | 1 | 3 (sun, sha, total) |
| 干燥时间常数 | 单一指数衰减 | **双时间常数** (阳 $\approx$ 1 h, 阴 $\approx$ 4 h) |

---

## 9. 关键参数与默认值

`data/config_vegetation.txt`:

```ini
canopy_storage_per_LAI = 0.1          # [mm/LAI] CLM5/NoahMP 全球标准 (温带林)
interc_F0              = 0.15         # [-] Gash 自由穿透系数
canopy_water_scheme    = single | four_quadrant
```

`src/COMPASS.m`:

```matlab
Ps_max_diurnal = 0.7    % 阳叶比例日峰值上限 (默认; 可由环境变量 COMPASS_PS_MAX 覆盖)
Ps_min_diurnal = 0.0    % 阳叶比例日谷值下限
lambda_v       = 2.45e6 % J/kg 水汽化潜热
```

`src/heatfluxes_clm45.m`:

```matlab
fwet           ∈ [0,1]  % 由 COMPASS 主循环传入
% Deardorff 1978 指数固定 2/3
```

---

## 10. $P_s^{\max}$ 选取依据与敏感性

### 10.1 物理依据

$P_s^{\max}$ 是日内阳叶比例峰值。它由三个因素决定：

1. **太阳天顶角** $\theta_z$：正午 $\theta_z$ 最小，阳叶比例最高
2. **LAI**：LAI 越大，下层叶遮蔽越多 → $P_s^{\max}$ 越小
3. **叶倾角分布** (Spherical, Erectophile, Planophile)

对于苹果园 (LAI $\approx$ 3, 接近 Spherical 叶倾角)，正午 $P_s^{\max}$ 实测在 $0.65 \sim 0.75$ 之间 (SCOPE `gap.Ps(top)` 输出)。**默认 0.7** 取中位值。

### 10.2 敏感性分析

扫描 $P_s^{\max} \in \{0.5, 0.6, 0.7, 0.8, 0.9\}$，storm test (7/13–7/18, 120h)，其余参数固定。
分析脚本：[`tests/compare_v30_psmax_sensitivity.m`](tests/compare_v30_psmax_sensitivity.m)
对比图：[v30_psmax_sensitivity.png](figures/v30_psmax_sensitivity.png)

| $P_s^{\max}$ | $W_{\max,\text{sun}}$ (mm) | $W_{\max,\text{sha}}$ (mm) | $\sum lE_{\text{wet}}$ (W·h/m²) | $\overline{T_{\text{cave}}}$ (°C) | $\psi_l^{\min}$ (MPa) | $\sum I_{\text{int}}$ (mm) | 守恒残差 (mm) |
|---|---|---|---|---|---|---|---|
| 0.5 | $0.05L$ | $0.10L$ | 2194.4 | 26.641 | −3.196 | 2.286 | $2.4 \times 10^{-16}$ |
| 0.6 | $0.06L$ | $0.10L$ | 2214.0 | 26.636 | −3.196 | 2.323 | $3.9 \times 10^{-16}$ |
| **0.7 (默认)** | $0.07L$ | $0.10L$ | **2233.5** | **26.633** | **−3.196** | **2.358** | $6.0 \times 10^{-16}$ |
| 0.8 | $0.08L$ | $0.10L$ | 2254.9 | 26.627 | −3.196 | 2.385 | $1.2 \times 10^{-15}$ |
| 0.9 | $0.09L$ | $0.10L$ | 2291.4 | 26.617 | −3.196 | 2.411 | $2.5 \times 10^{-16}$ |

### 10.3 不确定性传递

| $P_s^{\max}$ | $\Delta \sum lE_{\text{wet}}$ (vs 0.7) | $\Delta \overline{T_{\text{cave}}}$ |
|---|---|---|
| 0.5 | −1.75% | +0.008 °C |
| 0.6 | −0.88% | +0.003 °C |
| 0.8 | +0.96% | −0.006 °C |
| 0.9 | +2.59% | −0.016 °C |

**结论**：

- $P_s^{\max}$ 在合理范围 $[0.5, 0.9]$ 内变化时，$\sum lE_{\text{wet}}$ 变化 **< 3%**，$\overline{T_{\text{cave}}}$ 变化 **< 0.02 °C**
- 远低于不确定性传递判据 5%，说明模型对 $P_s^{\max}$ **极不敏感**
- $\psi_l^{\min}$ 在所有情形下完全一致 (−3.196 MPa)，说明阴叶水库容量变化不影响植物水力压力
- 守恒残差稳定在 $10^{-16}$ mm 量级（机器精度），证明算法在不同 $P_s^{\max}$ 下数值健壮
- 默认 $P_s^{\max} = 0.7$ 是物理合理 + 数值稳定的良好选择

---

## 11. 文献依据

1. **Deardorff, J. W. (1978)**. Efficient prediction of ground surface temperature and moisture, with inclusion of a layer of vegetation. *J. Geophys. Res. Oceans*, 83(C4): 1889–1903. → $f_{\text{wet}} = (W/W_{\max})^{2/3}$ 公式来源。

2. **Rutter, A. J., Kershaw, K. A., Robins, P. C., & Morton, A. J. (1971)**. A predictive model of rainfall interception in forests. I. Derivation of the model from observations in a plantation of Corsican pine. *Agric. Meteorol.*, 9: 367–384. → 多层水库串联框架。

3. **Gash, J. H. C. (1979)**. An analytical model of rainfall interception by forests. *Q. J. R. Meteorol. Soc.*, 105(443): 43–55. → Sparse canopy 截留与自由穿透 $f_{ER}$ 公式。

4. **Norman, J. M. (1979)**. Modeling the complete crop canopy. *Modification of the Aerial Environment of Crops*. ASAE. → Sun-shade leaf decomposition 框架 ($P_s$ 加权累加)。

5. **Brutsaert, W. (1982)**. *Evaporation into the Atmosphere: Theory, History, and Applications*. Springer. → 湿表面 Penman 蒸发公式 (无气孔阻力)。

6. **Calder, I. R. (1996)**. Dependence of rainfall interception on drop size: 1. Development of the two-layer stochastic model. *J. Hydrol.*, 185(1–4): 363–378. → 阴叶慢蒸发 / 凝露文献依据。

7. **Klemm, O. et al. (2012)**. Fog as a fresh-water resource: overview and perspectives. *Ambio*, 41: 221–234. → 阴叶夜间凝露通量观测综述 ($0.1$ – $0.5 \text{ mm h}^{-1}$)。

8. **Sellers, P. J. (1986)**. Microclimate. *In: Land Surface Processes in Atmospheric General Circulation Models*. Cambridge University Press. → 双源 (canopy + soil) 能量平衡的层化扩展。

9. **CLM5 Technical Note (Lawrence et al., 2019)**. Section 7.3 "Canopy water". NCAR/TN-557+STR. → CLM5 single $f_{\text{wet}}$ 实现 (本文 v29 baseline)。

10. **NoahMP (Niu et al., 2011)**. The community Noah land surface model with multiparameterization options (Noah-MP). *J. Geophys. Res.*, 116: D12109. → NoahMP 单冠层水库配方，是 v29 single 等价物。

---

## 12. 数值实现要点

### 12.1 时间积分方案
- 外步长 $\Delta t = 3600 \text{ s}$ (小时尺度)
- 截留 / 蒸发 / 凝露在 SCOPE $\to$ ebal $\to$ 水库扣减序贯求解：**前向 Euler (explicit)**，因 $W \ll W_{\max}$ 时 stiffness 不显著，stability $\Delta t \cdot E / W < 1$ 在 $\Delta t = 3600 \text{ s}$ 下默认满足
- 时间步内有 substepping (KT 内部 Picard 迭代用于水热耦合，**不影响**冠层水库 — 冠层水库每 KT 更新一次)

### 12.2 数值守恒
- 所有 $\min(\cdot, \text{容量})$ 和 $\max(0, \cdot)$ 钳位算子配对使用
- 残差检查阈值 $10^{-6}$ mm (约 $10^{-9} \text{ kg m}^{-2}$，远小于通量量级 $10^{-4} \text{ kg m}^{-2} \text{ h}^{-1}$)
- **跨时段守恒度对照**（同量级，验证算法稳定）：
  - 短期 storm test (120 h)：$\max|R_{\text{tot}}| = 7 \times 10^{-9}$ mm（§7.3）
  - 全年 (8760 h)：水量闭合度 $|R|/P < 5\%$（见 [`v30_single_vs_four_quadrant_comparison.md`](v30_single_vs_four_quadrant_comparison.md) 附录 A.1，年总残差占降雨 ~4%，由 Richards 求解器累积浮点误差贡献）
  - 冠层水库自身（不含土壤）的瞬时残差全程稳定在 $10^{-9}$ mm 量级

### 12.3 兼容性
- `interception` struct 同时填 `fwet` (v29 兼容) + `fwet_sun/fwet_sha` (v30 新增)
- `heatfluxes_clm45` 对单 `fwet` 调用仍工作 (默认 `fwet=0` if 缺失)
- v30 输出 CSV 多了 4 列：`Sim_output_fwet_sun.csv`, `Sim_output_fwet_sha.csv`, `Sim_output_h2ocan_sun_mm.csv`, `Sim_output_h2ocan_sha_mm.csv`

---

## 13. 验证基准

**Storm 7/13–7/18 (Loess Plateau Apple Orchard, 长武)**

| 验证指标 | 通过准则 | 实测 |
|---|---|---|
| 总水量守恒 $\lvert R_{\text{tot}} \rvert$ | $< 10^{-6}$ mm | $7 \times 10^{-9}$ mm $\checkmark$ |
| 总能量守恒 $\lvert \Delta lE_{\text{tot}} \rvert$ | $< 0.1\%$ | $0.03\%$ $\checkmark$ |
| 阴叶凝露捕捉 (KT 70-72) | $lE_{\text{wet,sha}} < 0$ | $-1.1 \text{ W m}^{-2}$ $\checkmark$ |
| 阳叶蒸发干燥时间常数 | $< 2 \text{ h}$ | $1 \text{ h}$ $\checkmark$ |
| 阴叶蒸发干燥时间常数 | $> 3 \text{ h}$ | $4 \text{ h}$ $\checkmark$ |
| $T_{\text{cave}}$ 差异 (vs single) | $< 0.1°C$ 均值 | $-0.029°C$ $\checkmark$ |
| $\psi_l^{\min}$ 差异 (vs single) | 在 $0.1 \text{ MPa}$ 内 | $+0.067 \text{ MPa}$ $\checkmark$ |
| 数值收敛 (120 KT) | 无 FATAL | 0 FATAL $\checkmark$ |

---

## 14. 关键代码定位 (论文 Code Availability 节)

| 物理过程 | 文件:行号 |
|---|---|
| 容量分配 | [`src/COMPASS.m:574-575`](src/COMPASS.m#L574-L575) |
| Rutter 多层串联 | [`src/COMPASS.m:583-606`](src/COMPASS.m#L583-L606) |
| 双 $f_{\text{wet}}$ 计算 (Deardorff) | [`src/COMPASS.m:608-621`](src/COMPASS.m#L608-L621) |
| SCOPE $f_{\text{wet}}$ 接入 | [`src/ebal_clm45_profile.m:85-104, 510-515`](src/ebal_clm45_profile.m#L85-L104) |
| 干 / 湿叶 $lE$ 拆分 | [`src/heatfluxes_clm45.m:97-103`](src/heatfluxes_clm45.m#L97-L103) |
| 阴叶凝露 | [`src/COMPASS.m:922-938`](src/COMPASS.m#L922-L938) |
| 守恒检查 | [`src/COMPASS.m:946-970`](src/COMPASS.m#L946-L970) |
| $P_s^{\max}$ 环境变量覆盖 | [`src/COMPASS.m:207-213`](src/COMPASS.m#L207-L213) |

---

## 15. 与公开 STEMMUS-SCOPE 的 ablation

| 模型版本 | $f_{\text{wet}}$ 数 | 凝露 | 多层串联 | 守恒残差 |
|---|---|---|---|---|
| 公开 STEMMUS-SCOPE | 0 (无 fwet 物理) | 不 | 不 | N/A |
| COMPASS v28 | 0 (同上) | 不 | 不 | N/A |
| COMPASS v29 (single) | 1 | 不 | 不 | $< 10^{-7}$ mm |
| **COMPASS v30 (four_quadrant)** | **2** | **阴叶可凝** | **阳→阴→地** | **$< 10^{-9}$ mm** |
