# 大气匹配通量顶热边界：原理、动机与代码

**配置名：** `heat_top_kind = "atmospheric"`（长武生产）  
**符号：** 通量进入土柱为正（W m⁻²）；内核 J cm⁻² s⁻¹ 为 `×1e-4`。

---

## 1. 为什么这么做

### 旧路：钉皮温 Dirichlet

ebal 解出皮肤温度 \(T_{ss}\)，土柱顶节点强制 \(T_N=T_{ss}\)。  
土柱反应通量 \(G_{\mathrm{soil}}\) 由矩阵解出，**不必等于** 地表可用能量 \(G_{\mathrm{sfc}}\)。

| 问题 | 后果 |
|------|------|
| \(G_{\mathrm{soil}}\neq G_{\mathrm{sfc}}\) | 界面能量账开裂 |
| \(T_{ss}\) 难步跳变 | Picard 难收敛、\(\Delta t\) 暴缩 |
| force-restore 经验 \(G\) 与土反应混用 | 诊断与权威边界不清 |

### 不做的路

| 方案 | 原因 |
|------|------|
| Picard 内反复重跑 ebal | 过贵、外层无界 |
| 裸地+雪聚合 Robin 直接进 \(h/U/P_g\) 顶行 | 部分雪盖数值不稳 |
| 抬高 Medlyn \(g_0\) 做“热散热” | 污染夜/旱残余（旁线气孔另论） |

### 本路图像

1. **表面**：ebal + 雪仍算一次，冻结仿射 \(G_s(T)\)（含裸地与雪，进土为正，斜率 \(a\le 0\)）。  
2. **土侧预测**：上层 **20 cm** 隐式热矩阵 → admittance \(K\) 与零通量顶温 \(T_0\)，得 \(G_g(T)=K(T-T_0)\)。  
3. **匹配**：\(G_s(T^*)=G_g(T^*)\) → 共同通量 \(G^*\)。  
4. **正式土柱**：整小时（含子步）用 **固定 Neumann \(G^*\)** 推 \(h/U/P_g\)；**不钉** \(T_N=T_{ss}\)。  
   皮温只作 ebal 诊断皮肤。

人话：**大气–冠层定“这小时该进多少热”；土柱自己演化顶温。**

---

## 2. 原理与公式

### 2.1 地表仿射 \(G_s(T)\)

雪盖率 \(f\in[0,1]\)，统一线性化温度 \(T_{\mathrm{ref}}\)：

\[
\begin{aligned}
G_b^{\mathrm{ref}} &= (1-f)\big[q_b + a_b(T_{\mathrm{ref}}-T_{b,\mathrm{lin}})\big],\\
G_{si}^{\mathrm{ref}} &= f\,q_{si}^{\mathrm{ref}} + q_{\mathrm{redistrib}},\\
a &= (1-f)a_b + f\,a_{si},\qquad a\le 0,\\
G_{\mathrm{ref}} &= G_b^{\mathrm{ref}} + G_{si}^{\mathrm{ref}},\\
G_s(T) &= G_{\mathrm{ref}} + a\,(T - T_{\mathrm{ref}}).
\end{aligned}
\]

- \(q_b,a_b\)：裸地外部能量平衡及其对土表温导数（进土为正时 \(a_b\le 0\)）。  
- \(q_{si},a_{si}\)：雪–土导热 \(C_{si}(T_{\mathrm{snow}}-T)\) 的参考值与斜率（\(a_{si}=-C_{si}\)）。  
- \(q_{\mathrm{redistrib}}\)：融尽等 **与 \(T\) 无关** 的焓重分配，只进截距。  
- **不含** SCOPE force-restore 经验 \(G\)。

### 2.2 上层 20 cm admittance \(K,T_0\)

接受态水/冰/气/热系数冻结。在 \(z=0\) 到 **第一个深度 ≥ 0.20 m 的节点** 上截断热系统；截断底节点 **Dirichlet 固定为接受态温度**。

| 探针 | 含义 |
|------|------|
| 顶通量 \(=0\) | 顶温 \(T_0\) |
| 顶通量 \(+1\,\mathrm{W\,m^{-2}}\) | \(\partial T/\partial G\) |

\[
K = \Big(\frac{\partial T}{\partial G}\Big)^{-1},\qquad
\frac{\partial T}{\partial G}>0,\qquad
G_g(T)=K\,(T-T_0).
\]

全柱浅于 20 cm 时不截断，保留配置底边界。  
**注意：** 20 cm 只服务 **匹配预测**；正式 \(h/U/P_g\) 仍在 **整柱** 上解。

### 2.3 解析匹配

\[
G_s(T)=G_{\mathrm{ref}}+a(T-T_{\mathrm{ref}}),\quad
G_g(T)=K(T-T_0).
\]

\[
T^*=\frac{K T_0 + G_{\mathrm{ref}} - a\,T_{\mathrm{ref}}}{K-a},\qquad
G^*=K\,(T^*-T_0).
\]

因 \(K>0\)、\(a\le 0\)，分母 \(K-a>0\)，交点稳定。

### 2.4 正式边界与残差

正式土柱顶：

\[
q_{\mathrm{top}}^{\mathrm{in}} \equiv G^*\quad\text{（固定 Neumann，进土为正）}.
\]

内核：`BCT = -q_in`，装配 `RHS(N) -= BCT` ⇔ `RHS += q_in`。

求解后评价温 \(T_{\mathrm{eval}}\) 上：

\[
r = G_{\mathrm{soil}}(T_{\mathrm{eval}}) - G_{\mathrm{sfc}}(T_{\mathrm{eval}}),
\]

\(G_{\mathrm{sfc}}\) 由冻结仿射再算（可含实现蒸发对潜热的修正）。  
\(|r|\le\varepsilon_G\)（量级 5 W m⁻²）接受；否则割线修正试验通量再解，次数有界。

---

## 3. 逻辑线路（一个大气小时）

```text
ebal / ebalProfile
  → 叶温、H/LE、皮温 Tss（诊断）
  → 不把 Tss 当土顶 Dirichlet

post-ebal + 雪
  → buildSoilTopHeatResponse
  → top_heat_response  {G_ref, a, T_ref, f, …}  冻结

advanceSoilWithInterfaceTemperatureClosure  (kind == atmospheric)
  → estimateDynamicSurfaceAdmittance   (上层 20 cm → K, T_0)
  → matchDynamicSoilHeatResponse       (解析 G*_pred, T*_pred)
  → soil_top_downward_energy_flux = G*_trial
  → materialize → atmospheric = 固定 Neumann
  → advanceSoilEquationOuterStep       (h/U/P_g Picard，可子步)
  → 评价 r = G_soil − G_sfc(T_eval)
  → |r| 过大则割线改 G*_trial 再试

接受后
  → 雪焓 reconcile（时间加权界面温）
  → commit / 诊断
```

**不变量：** 正常小时 **Picard 内不重跑 ebal**；仿射在 post-ebal 冻结。

---

## 4. 对应代码

| 环节 | 文件 |
|------|------|
| 小时编排 | `src/runOneStep.m` |
| 匹配 + Neumann 闭合 | `src/+solver/advanceSoilWithInterfaceTemperatureClosure.m`（`advanceDynamicFluxClosure`） |
| 20 cm admittance | `src/+soil/+thermal/estimateDynamicSurfaceAdmittance.m` |
| 解析交点 | `src/+surface/+energy/matchDynamicSoilHeatResponse.m` |
| 地表仿射 | `src/+surface/+energy/buildSoilTopHeatResponse.m` |
| 冻结图评价 \(G_{\mathrm{sfc}}\) | `src/+surface/+energy/evaluateFrozenSurfaceSoilHeat.m` |
| kind 规范化 | `src/+soil/+thermal/+boundary/normalizeSpec.m`（`atmospheric`） |
| 冻结 BC 请求 | `src/+soil/+thermal/+boundary/materializeRequest.m` |
| 试验 drivers | `src/+solver/materializeSoilBoundaryTrial.m` |
| 热矩阵顶行 | `src/+soil/+thermal/calculateBoundaryConditions.m` / `solveThermalTrial.m` |
| 站点配置 | `data/sites/changwu/site.toml` → `heat_top_kind = "atmospheric"` |

旁线（默认关，不进本路径主线）：可选 \(g_{\mathrm{heat}}\) 见 `docs/plant/stomata_heat_induced_gs_design.md` 与 `vegetation.biochemical`。
