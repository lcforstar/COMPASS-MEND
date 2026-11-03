# 热驱动气孔开度（Heat-induced \(g_s\)）设计与离线试验

**状态：** 设计 + 离线验证 + **已接入 ebal/biochemical**；长武 **`heat_gs_enable=false` 默认关**  
**日期：** 2026-07-19 / 2026-07-20  
**相关代码（离线）：** [`tests/plant/offline_heat_gs_demo.m`](../../tests/plant/offline_heat_gs_demo.m)  
**相关生产路径：**  
- `wt` → `sfactor` → `Vcmax/Jmax/TPU`：[`src/+plant/biochemical.m`](../../src/+plant/biochemical.m) §水力胁迫  
- 可选 \(g_{\mathrm{heat}}\)：同文件 Medlyn 后加法项（`heat_gs_*`）  
- ebal 传入 `sfactor` 与 leafbio：[`src/+plant/+energy/ebal.m`](../../src/+plant/+energy/ebal.m)  
- 配置：[`data/sites/changwu/site.toml`](../../data/sites/changwu/site.toml) `[vegetation]`  
**交接总览：** [`docs/surface/2026-07-20_atmospheric_matched_flux_guide.md`](../surface/2026-07-20_atmospheric_matched_flux_guide.md) §6  
**相关长文：** [`stomata_thermal_decoupling_model.md`](stomata_thermal_decoupling_model.md)（热–气孔正式计划 v3；本文件是其可落地的 **可选 C 路径**，不自动覆盖 v3 的 Stage 顺序）

---

## 0. 一句话

> **干旱路径（A/B）已有：`wt` 压 Vcmax → An↓ → Medlyn \(g_{s0}\)↓。**  
> **可选热路径（C）：在 \(g_{s0}\) 上加法附加 \(g_{\mathrm{heat}}\)（不放大 Medlyn \(g_0\)），且仅在有水、够热、VPD/An 允许时打开；默认关。**

---

## 1. 问题背景

### 1.1 已有 A/B（关孔路径）

| 步骤 | 行为 |
|------|------|
| SPAC / 脆弱性 | 得 `wt ∈ [0,1]` |
| `biochemical` | `Vcmax,Jmax,TPU × sfactor`（`sfactor=wt`）；**Rd 不乘** |
| Medlyn | \(g_{s0}=g_0+1.6\,g_1 A_n/(C_a\sqrt{D})\)，\(g_0\sim 10^{-4}\) 极小 |
| 能量平衡 | \(r_{cw}\) 进 H/LE，解叶温 \(T_c\) |

人话：**渴了关** → 蒸发少 → 叶温往往升高。  
这不是“热了主动开孔散热”。

### 1.2 为何不能靠放大 \(g_0\) 做热散热

Medlyn 残余导度 \(g_0\) 本就极小。干旱时 An→0 ⇒ \(g_s\to g_0\)，残余蒸发散热可忽略。

- 把 \(g_0\) 调大：污染夜间/极旱残余，**假机制**  
- \(g_s=g_{s0}(1+\phi)\) 在 \(g_{s0}\approx g_0\) 时仍≈0，**热开孔加不上**

因此 C 必须用 **加法** 有界项 \(g_{\mathrm{heat}}\)，且 **`wt` 一票否决**。

### 1.3 与 `stomata_thermal_decoupling_model.md` 的关系

| 文档 | 角色 |
|------|------|
| `stomata_thermal_decoupling_model.md` | 正式计划 v3：先分清热土≠热叶、诊断、再谈液压气孔等 |
| **本文** | 一条 **可开关的热附加 \(g_s\)** 规格 + 离线数字；**未授权默认上线** |

实施前仍须满足 v3 的“先诊断、热土≠热叶”约束；本文不替代 v3 的 Stage 0。

---

## 2. 规格公式

### 2.1 总式

\[
g_{s0}
=
g_0
+
1.6\,g_1
\frac{\max(A_n,0)}{C_a\sqrt{\max(D,D_{\min})}}
\]

\[
g_{\mathrm{heat}}
=
g_{\mathrm{heat,max}}
\cdot
\phi_T(T_c)
\cdot
\phi_w(wt)
\cdot
\phi_D(D)
\cdot
\phi_A(A_n)
\]

\[
g_s
=
\min\bigl(g_{s0}+g_{\mathrm{heat}},\, g_{s,\mathrm{cap}}\bigr)
\]

单位：\(g_s,g_{s0},g_{\mathrm{heat}}\) 为 \(\mathrm{mol\,m^{-2}\,s^{-1}}\)（与 Medlyn 一致，再换 \(r_{cw}\)）。

### 2.2 门函数（线性 ramp）

\[
\mathrm{ramp}(x;x_0,x_1)
=
\begin{cases}
0 & x\le x_0\\
(x-x_0)/(x_1-x_0) & x_0<x<x_1\\
1 & x\ge x_1
\end{cases}
\]

| 因子 | 定义 | 建议初值 |
|------|------|----------|
| \(\phi_T\) | \(\mathrm{ramp}(T_c;T_{\mathrm{on}},T_{\mathrm{full}})\) | \(T_{\mathrm{on}}=32^\circ\mathrm{C}\), \(T_{\mathrm{full}}=38^\circ\mathrm{C}\) |
| \(\phi_w\) | \(\mathrm{ramp}(wt;wt_{\mathrm{off}},wt_{\mathrm{on}})\) | \(wt_{\mathrm{off}}=0.4\), \(wt_{\mathrm{on}}=0.8\) |
| \(\phi_D\) | \(1-\mathrm{ramp}(D;D_{\mathrm{lo}},D_{\mathrm{hi}})\) | \(D_{\mathrm{lo}}=1.5\,\mathrm{kPa}\), \(D_{\mathrm{hi}}=3.5\,\mathrm{kPa}\) |
| \(\phi_A\) | \(\mathrm{ramp}(A_n;A_{\mathrm{lo}},A_{\mathrm{hi}})\) | \(A_{\mathrm{lo}}=1\), \(A_{\mathrm{hi}}=5\) µmol m⁻² s⁻¹ |
| \(g_{\mathrm{heat,max}}\) | 热项满额 | **0.08** mol m⁻² s⁻¹（先敏感性） |
| \(g_{s,\mathrm{cap}}\) | 总上限 | 0.6 mol m⁻² s⁻¹ |
| \(g_0\) | Medlyn 残余 | **保持 ~1e-4，禁止抬高** |

### 2.3 优先级（必须写死）

```text
1. enable=false → g_heat ≡ 0（生产默认）
2. wt ≤ wt_off → g_heat ≡ 0（干旱 = 现 A/B）
3. Tc ≤ T_on 或 D ≥ D_hi 或 An ≤ A_lo → 对应 φ=0
4. 永不通过放大 g0 实现热开孔
5. 不与土壤 Picard / 皮温–土温接口混 PR
```

### 2.4 与现有 `wt→Vcmax` 的共存顺序

```text
猜 Tc
  → biochemical：Vcmax×wt → An → Medlyn gs0
  → [可选 C] gs = gs0 + g_heat(Tc,wt,D,An)
  → rcw → H, LE
  → 能量残差更新 Tc
  → 水力内环更新 wt 后重算（含 g_heat）
```

阳/阴叶（及层）各自用自己的 \(T_c,wt,A_n,D\)。

---

## 3. 配置（已写入 site vegetation，默认关）

扁平键挂在 `[vegetation]`（与 `gs_g0`/`gs_g1` 同级），经 `init.buildVeg` →
`veg.photosynthesis.heat_gs_*` → `prepareScopeStep`/`leafbio` →
`ebal`/`ebalProfile` `apply_phase3_biochem_fields` → `vegetation.biochemical`。

```toml
# data/sites/changwu/site.toml [vegetation]  — production default OFF
heat_gs_enable   = false
heat_gs_T_on_C   = 32.0
heat_gs_T_full_C = 38.0
heat_gs_max      = 0.08    # mol m-2 s-1
heat_gs_wt_off   = 0.4
heat_gs_wt_on    = 0.8
heat_gs_D_lo_kPa = 1.5
heat_gs_D_hi_kPa = 3.5
heat_gs_A_lo     = 1.0     # umol m-2 s-1
heat_gs_A_hi     = 5.0
heat_gs_cap      = 0.6     # only applied when enable=true
```

**有实测再开**：`heat_gs_enable = true` 并按站点校准阈值；开前先做 E1–E5。

---

## 4. 离线玩具试验

### 4.1 如何运行

```matlab
cd F:\COMPASSV1.0
addpath('tests/plant')
offline_heat_gs_demo
```

脚本：[`tests/plant/offline_heat_gs_demo.m`](../../tests/plant/offline_heat_gs_demo.m)  
**不**调用 ebal / SPAC / 主循环。

### 4.2 工况表结果（2026-07-19 跑通）

| case | Tc | An | D | wt | gs0 | g_heat | **gs** | 解读 |
|------|---:|---:|---:|---:|---:|---:|---:|------|
| wet_mild | 28 | 15 | 1.0 | 1.0 | 0.281 | **0** | 0.281 | 不热，无热项 |
| wet_hot | 36 | 12 | 1.8 | 1.0 | 0.168 | **0.045** | **0.213** | 有水+热，gs↑~27% |
| wet_vhot | 40 | 10 | 2.2 | 1.0 | 0.126 | 0.052 | 0.178 | 热项接近满 |
| dry_hot | 36 | 1 | 2.5 | 0.25 | 0.012 | **0** | 0.012 | **wt 锁死热项** |
| dry_vhot | 40 | 0.2 | 3.0 | 0.2 | 0.002 | **0** | ~g0 量级 | 旱+An 塌 |
| mid_wt_hot | 36 | 6 | 2.0 | 0.55 | 0.080 | 0.015 | 0.095 | 半旱打折 |
| wet_hot_hiVPD | 36 | 12 | 4.0 | 1.0 | 0.113 | **0** | 0.113 | φ_D=0 |

### 4.3 敏感性摘要

- **扫 wt**（Tc=36, An 固定）：`wt≤0.4` → `g_heat=0`；`0.4–0.8` 线性开；`≥0.8` 热项满（在 φ_T、φ_D 允许下）。  
- **扫 Tc**（wt=1）：`Tc≤32` 无热项；`32–38` 上升；`≥38` 封顶。  
- **玩具能量平衡**（方向性）：`wt`↓ → LE↓ → 平衡叶温↑，与 A/B 故事一致。

### 4.4 离线结论

| 检查项 | 结果 |
|--------|------|
| 不靠抬 g0 | 通过 |
| 湿热 gs > gs0 | 通过 |
| 旱热 g_heat=0 | 通过 |
| 高 VPD 可关热项 | 通过 |
| 可与 A/B 共存（门控） | 通过（玩具层） |

**离线只证明门控逻辑；接入 ebal 后须重做 E1–E5 验收。**

---

## 5. 生产接入前验收（若将来实现）

| ID | 设置 | 期望 |
|----|------|------|
| E1 | `enable=false`，原干旱个例 | 与现网一致 |
| E2 | `enable=true`，湿润+高温 | gs↑、LE↑、\(T_c-T_a\) 略↓ |
| E3 | `enable=true`，干旱高温 | 与 false 几乎相同 |
| E4 | `enable=true`，温和 \(T_c<T_{\mathrm{on}}\) | 与 false 几乎相同 |
| E5 | An≈0 极端 | \(g_{\mathrm{heat}}\) 被 φ_A / φ_w 压住 |

建议诊断字段：`gs0_medlyn`, `g_heat`, `gs_final`, `phi_T/W/D/A`, `Tc`, `wt`, `LE_dry`。

---

## 6. 明确不做

| 不做 | 原因 |
|------|------|
| 默认 `enable=true` | 保护已测通 A/B |
| 放大 Medlyn \(g_0\) | 假热机制 |
| 强制 \(T_c\) 到设定温度 | 破坏能量平衡 |
| 在土壤 Picard 内算气孔 | 架构错误 |
| 与皮温–土温 \(K\) 消元同一 PR | 混线 |
| 用本离线表当标定 | 玩具 An/D，非长武时序 |

---

## 7. 与土接口议题的边界

| 议题 | 文档/状态 |
|------|-----------|
| 皮温 Dirichlet / KT8106 | 审计与接口设计；**另一条线** |
| resume `Y` 污染修复 | 已实现，保留 |
| 本热开孔 C | **可选叶生理扩展**；默认关 |

---

## 8. 建议决策

| 主科学问题 | 建议 |
|------------|------|
| 干旱、水力、产量 | **可不实现 C**；方法学写明“无热驱动开孔” |
| 有水热浪 + 叶温/gs 观测 | 按本文规格实现，默认关，先 E1–E4 |
| 仅想“叶温好看” | **不做** |

**当前仓库状态（2026-07-20）：** 生产链路已接线，**默认 `heat_gs_enable=false`**（`g_heat≡0`，行为与纯 Medlyn 一致）。开开关前须有实测再做 E1–E5 验证。

---

## 9. 修订记录

| 日期 | 内容 |
|------|------|
| 2026-07-19 | 讨论定稿：加法 \(g_{\mathrm{heat}}\)、门控、不抬 g0 |
| 2026-07-19 | `offline_heat_gs_demo.m` 跑通工况表与敏感性 |
| 2026-07-20 | 本文写入 `docs/plant/` |
| 2026-07-20 | 接入 ebal/biochemical + site vegetation 参数；开关默认关 |
