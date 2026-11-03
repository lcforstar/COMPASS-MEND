# LAI 垂直分层与冠层耦合改进计划

作者：朱寄子星

日期：2026-04-19

## 目的

本文档记录当前 `ebal_clm45` 冠层结构改进的完成状态和后续计划。当前阶段的目标是先引入可配置的 LAI 垂直分布，并且不破坏已有 SCOPE 辐射传输和能量平衡流程。后续再逐步把 LAI 分层接入总通量积分、辐射传输、植物水力和根系吸水。

冠层层序遵循 SCOPE 约定：

```text
layer 1  = 冠层顶部
layer nl = 冠层底部 / 靠近土壤一侧
```

注意：这和 STEMMUS 一些土壤数组的层序不同。土壤变量里有些地方使用 `NL` 表示表层或靠近边界的层，需要单独判断，不能和冠层 `profiles` 的层序混用。

## 已完成改进

### 1. CLM4.5 能量平衡分支注释整理

已对 `src/ebal_clm45.m`、`src/resistances_clm45.m` 及相关辅助函数增加说明。当前 `ebal_clm45` 应理解为：

```text
SCOPE 辐射/冠层结构 + Wallace-Verhoef 冠层阻力网络 + CLM4.5/CLM5 风格空气动力和边界层改进
```

它不是纯 CLM big-leaf 模型，而是一个 SCOPE 与 CLM4.5 空气动力学改进的混合分支。

已完成内容：

- 更新 `ebal_clm45` 文件头注释。
- 增加动态土壤阻力、粗糙度、热发射率、空气动力阻力相关注释。
- 明确作者、日期和方法文档引用。

### 2. 使用 SCOPE `gap.k` 作为直射光消光系数

当前 `ebal_clm45` 主链路已经移除 Campbell-Norman fallback 消光系数。动态叶宽使用的直射光消光系数来自 SCOPE RTMo：

```matlab
canopy.Kb_direct = gap.k;
resist_in.Tao = canopy.Kb_direct;
```

物理含义：

- `gap.k` 是太阳直射方向的消光系数，已经按叶倾角分布 `lidf` 积分。
- `gap.K` 是观测方向 / view direction 的消光系数，不应用于直射光动态叶宽计算。

已完成内容：

- 删除 `src/+Radiation/LightExtinction.m`。
- 删除 `Campbell_Norman_fallback` 和 `Tao_dyn` 使用。
- `ebal_clm45` 不再从 `cfg.Veg.Tao` 读取静态 Tao。
- 增加 `thermal.Kb_direct`、`thermal.Kb_source`、`thermal.Tao` 诊断输出。

### 3. 叶宽开关统一

当前叶宽开关为：

```text
CropCode == 1  使用用户输入的 w_leaf_max
CropCode ~= 1  模型自动计算实际叶宽
```

自动计算公式：

```matlab
w_leaf_actual = w_leaf_max * (1 - exp(-Tao * LAI));
```

在 `ebal_clm45` 主链路中，`Tao` 通常就是 SCOPE 的 `gap.k`。

已完成内容：

- 更新 `src/resistances_clm45.m`。
- 更新 `data/config_vegetation.txt` 中的开关注释。
- 同步更新旧的优化阻力函数，使开关语义一致。
- 已用 MATLAB 最小算例验证开关行为。

### 4. 动态土壤和冠层边界参数

`ebal_clm45` 在能量平衡迭代前会更新以下参数：

- `soil.rss`：由表层土壤含水量动态计算土壤表面水汽阻力。
- `soil.rs_thermal` / `soil.epsilon_s`：由表层饱和度动态计算土壤热反射率 / 发射率。
- `soil.rbs`：可选的 LAI 相关土壤边界层阻力。
- `canopy.zo` 和 `canopy.d`：由 LAI/SAI/冠层高度动态计算粗糙度长度和零平面位移。

注意：

`soil.rss` 当前是在每个模型时间步进入能量平衡迭代前，根据当前表层含水量计算一次。它没有在 `while CONT` 温度迭代内部重复更新，因为当前能量平衡内循环并不更新表层土壤水分。

### 5. 热量和水汽空气动力路径一致性修正

CLM4.5 分支现在区分：

```matlab
rah  热量空气动力阻力
raw  水汽空气动力阻力
```

`ebal_clm45` 中温度更新分母已改为使用对应路径：

```matlab
Tch 潜热导数项使用 rac_w + rcwh
Tcu 潜热导数项使用 rac_w + rcwu
Ts  潜热导数项使用 ras_w + rss
```

同时，当 `SoilHeatMethod < 2` 时，土壤温度更新分母加入土壤热通量项对当前土壤表面温度的导数：

```matlab
dG_dTs = 2 * GAM / sqrt(pi * Deltat);
```

这使土壤能量平衡残差：

```matlab
EBers = Rns - lEs - Hs - G;
```

的温度更新更完整。

### 6. 植物水力初始值移到配置文件

原本 `Init_Plant.m` 里硬编码的植物水力初始值已移入：

```text
data/config_vegetation.txt
```

当前配置字段：

```text
hydraulic_init_psir
hydraulic_init_pstem
hydraulic_init_psil_sun
hydraulic_init_psil_sha
```

这些值是第一步暖启动或未求解时的兜底值，不是水力求解成功后的最终叶水势 / 根水势。

### 7. LAI 垂直分布第一阶段

已新增：

```text
src/+plant/calc_lai_profile.m
tests/check_lai_profile_shapes.m
```

配置位置：

```text
data/config_vegetation.txt
```

支持的 LAI 垂直分布模式（共 13 种，覆盖全球主要 PFT）：

```text
uniform              通用对照基线
pine_cone            温带针叶（松树、冷杉，圆锥形）
boreal_needleleaf    北方针叶（云杉、落叶松，圆柱形）  ← 新增
tropical_rainforest  热带雨林（双峰：主冠层+林下层）   ← 新增
broadleaf_tree       温带/北方落叶阔叶（Beta 分布）
fruit_ellipsoid      椭球形果树冠层
willow_drooping      垂枝柳冠层
orchard_open_center  开心形修剪果园
maize                玉米/高秆作物（峰值 0.40）
crop_upper           C3 作物/草地（指数衰减，k 可调）
c4_grass             C4 草地/稀树草原（陡指数 k=3.5）  ← 新增
shrub_low            温带灌木
tundra               苔原/极地矮灌木（近均匀）          ← 新增
```

`run_stemmus_initialize.m` 中新增的冠层字段：

```matlab
canopy.LAI_layer
canopy.layer_weight
canopy.LAI_profile_shape
canopy.LAI_profile_z
canopy.LAI_profile_mode
canopy.LAI_profile_reference_LAI
```

当前阶段这些字段主要用于诊断和为下一步耦合做准备。默认模式为 `uniform`，因此不会改变原模型行为。

## 当前 LAI 分层公式

所有 LAI 分布都会归一化，使其满足：

```matlab
sum(LAI_layer) = LAI
sum(layer_weight) = 1
```

层中心深度定义为：

```matlab
z = ((1:nl)' - 0.5) / nl;
```

其中：

```text
z = 0  冠层顶部
z = 1  冠层底部
```

### `uniform`

每层 LAI 相同。

适合：

- 缺少冠层垂直结构资料；
- 与原始模型保持一致；
- 作为对照组。

### `pine_cone`

松树 / 针叶树圆锥冠层。冠层顶部叶面积少，向下单调递增。公式为 `z^k * exp(-0.5*z)`，当默认 `k=1.5` 时峰值在 `z=3`（超出 `[0,1]` 范围），因此在整个冠层内无底部削弱，形状单调上升。`k` 越大，底部 LAI 集中程度越高。

适合：

- 松树；
- 杉树；
- 圆锥形针叶冠层。

### `fruit_ellipsoid`

果树 / 椭球形冠层。叶面积在冠层中部最大，顶部和底部较少。

适合：

- 苹果；
- 梨；
- 桃；
- 近似椭球冠层的果树。

### `willow_drooping`

柳树 / 垂枝冠层。叶面积在中上部较大，同时向下有拖尾。

适合：

- 柳树；
- 垂枝型树冠；
- 下垂枝叶明显的冠层。

### `broadleaf_tree`

普通阔叶树冠层。使用 Beta 分布形状，自然地在冠层顶部和底部趋于零。叶面积峰值深度由 `LAI_profile_peak` 控制，集中程度由 `LAI_profile_k` 控制（`k` 越大分布越窄）。默认 `peak=0.50, k=1.5` 给出对称分布；设 `peak=0.35` 可还原偏上峰的阔叶林特征。

适合：

- 普通阔叶林；
- 自然林冠；
- 未明确修剪形态的树冠。

### `crop_upper`

一般作物 / 草地冠层。叶面积上部较多，向下递减。

适合：

- 草地；
- 小麦 / 水稻等上层截光明显的作物；
- 简化作物冠层。

### `maize`

玉米 / 高秆作物。高斯分布，**硬编码峰值 peak=0.40、sigma=0.22**，表示玉米冠层叶面积集中在中上部的生态学默认值。与 `fruit_ellipsoid`（peak=0.50）有意区分。如需自定义峰位，改用 `fruit_ellipsoid` 并设置 `LAI_profile_peak`。

适合：

- 玉米；
- 高秆作物；
- 中上部叶面积集中的作物。

### `shrub_low`

灌木或低矮冠层。叶面积向下增加，中下部较多。

适合：

- 灌木；
- 低矮枝叶密集冠层；
- 下层叶面积占比较大的植被。

### `orchard_open_center`

果园修剪 / 开心形冠层。类似椭球分布，但略偏向中下部。

适合：

- 修剪后的苹果园；
- 开心形或纺锤形果树；
- 中部和中下部叶面积较明显的果园冠层。

## 验证方式

运行：

```matlab
run('tests/check_lai_profile_shapes.m')
```

测试内容：

- 每种模式返回 `nl` 层；
- `sum(LAI_layer) = LAI`；
- `sum(layer_weight) = 1`；
- 所有值有限且非负。

输出文件：

```text
tests/figures/lai_profile_shapes_overlay.png
tests/figures/lai_profile_shapes_panels.png
tests/figures/lai_profile_shape_summary.csv
```

## 尚未完成的改进

### 0. ~~时序模式下 `canopy.LAI_layer` 需按步更新~~（已完成）

`AMainloop.m` 中 `io.select_input` 调用后已新增：

```matlab
if isfield(canopy, 'layer_weight') && ~isempty(canopy.layer_weight)
    canopy.LAI_layer = canopy.layer_weight * canopy.LAI;
end
```

每个时间步 `canopy.LAI` 更新后，`canopy.LAI_layer` 随即按当前 LAI 重新缩放。`layer_weight` 仍保持不变（归一化，与 LAI 无关）。

### 1. ~~`ebal_clm45` 总通量接入 LAI_layer~~（已完成）

`ebal_clm45.m` Section 5 已替换为 LAI_layer 加权积分：

```matlab
Ps_nl = Ps(1:nl);
LW    = canopy.LAI_layer;    % [nl×1]  各层 LAI，sums to LAI
WW    = canopy.layer_weight; % [nl×1]  归一化权重，sums to 1

Xu1d   = equations.meanleaf(canopy, Xcu, 'angles');   % 先对角度平均 → [nl×1]
Xtot   = sum(LW .* ((1 - Ps_nl) .* Xsh + Ps_nl .* Xu1d));
```

覆盖的变量：`Rnctot`、`lEctot`、`Hctot`、`Actot`、`Tcave`、`Pntot`、`Pntot_Cab`、`Rntot_PAR`、`aPAR_Cab_eta`。

within-loop 的植物水力蒸腾需求（`Tp_sha_m_s`、`Tp_sun_m_s`）也同步改为 LAI_layer 加权，确保水力求解的蒸腾需求与最终输出通量一致。

当 `canopy.layer_weight` 不存在时自动回退到原始等层公式，向后兼容。

### 2. ~~RTMo 非均匀层厚~~（已完成）

`RTMo`、`RTMt_sb`、`RTMf`、`RTMz` 已全部完成非均匀层厚改造。

改造方式统一：

```matlab
% 有 canopy.layer_weight 时用实际层厚，否则退回 1/nl
if isfield(canopy, 'layer_weight') && ~isempty(canopy.layer_weight)
    dx_layer  = canopy.layer_weight(:);        % [nl×1] 非均匀，sums to 1
else
    dx_layer  = ones(nl, 1) / nl;              % [nl×1] 均匀 fallback
end
LAI_layer = LAI * dx_layer;                   % [nl×1]
```

各文件改动摘要：

| 文件 | 改动 |
|------|------|
| `RTMo` | `dx/iLAI` → `dx_layer/LAI_layer`；`xl` 由 `cumsum(dx_layer)` 得出；Ps/Po 有限尺寸修正逐层化；Pso 积分逐层化；`rho_sos`、`piLoc_` 按 `LAI_layer` 加权 |
| `RTMt_sb` | `dx/iLAI` → `dx_layer/LAI_layer`；`fHc` 改为 `[nl×1]` 向量；两流迭代步长逐层；方向热辐射贡献 piLo1/2/3 按 `LAI_layer` 加权 |
| `RTMf` | 删除 `iLAI`；`piLo1/2`、`fac_layer`、`dF1/dF2`、G1/G2 循环、`piLo3`、`Fem` 全部改用 `LAI_layer` |
| `RTMz` | 与 RTMf 相同处理（无 `Fem` 行） |

验证结果（`tests/check_lai_layer_integration.m`，2026-04-19）：
- RTMo 能量守恒通过（LAI 加权 Rnhc 差值 = 0.058 W/m²）✅
- ebal 总通量差值 < 1 W/m²（lEctot = −1.11 W/m²）✅
- EB 残差 < 0.006 W/m² ✅

旧的 54 W/m² 差异根因：RTMt_sb 热发射仍用均匀 1/nl 权重，与 RTMo 的 LAI_layer 不一致，叠加 Clausius-Clapeyron 放大。RTMt_sb 修复后消失。剩余 ~1 W/m² 为真实物理非均匀效应。

### 3. ~~分层植物水力~~（已完成 v2，LAI 比例导度）

#### 3.0 第一版（2026-04-19，需求比例分配）

2026-04-19 已接入第一版逐层植物水力胁迫。当前实现仍沿用现有 SPAC / `Plant_hydraulics` 核心求解器，但在 `ebal_clm45` 内把逐层潜在蒸腾需求传入新的 wrapper：

```matlab
vegetation.hydraulics.solveEbalHydraulicsLayered
```

新增输出：

```matlab
Tp_sun_layer_demand
Tp_sha_layer_demand
Tp_sun_layer
Tp_sha_layer
wt_sun_layer
wt_sha_layer
psil_sun_layer
psil_sha_layer
```

`ebal_clm45` 中的 `profiles` 和 `fluxes` 已输出逐层 `wt`；非均匀 `LAI_layer` 情况下，生化模块使用逐层 `wt_sun_layer / wt_sha_layer` 反馈到 `Vcmax`，同时保留需求加权的整体 `wt_sun / wt_sha` 作为汇总诊断。显式 uniform profile 会退回到原 two-leaf 结果并复制成逐层数组，因此不改变旧模型行为。

验证结果：
- `tests/check_lai_uniform_regression.m` 通过：显式 uniform 分层与旧 fallback 差异为机器精度。
- `tests/check_lai_profile_physical_impacts.m` 通过：输出逐层蒸腾需求、实际蒸腾、`wt` 和 `psil`。
- `tests/check_lai_layer_hydraulic_drought.m` 通过：在 `crop_upper / fruit_ellipsoid / pine_cone`
  三种冠型下构造干旱情景。复合干旱
  (`ea=4 hPa`, `KL_h × 1e-8`, `psi_soil=-1.4 MPa`) 使供给比降到约 0.039-0.044。
  `crop_upper` 顶部 20% 阳叶 `wt_sun≈0.204`、底部 20% `≈0.678`；
  `fruit_ellipsoid` 顶部 `≈0.336`、底部 `≈0.679`；
  `pine_cone` 顶部 `≈0.578`、底部 `≈0.246`。结果证明逐层水力胁迫不仅能分辨冠层上下差异，也能体现冠型差异。
- 极端供水受限时，`ebal_clm45` 会自动放宽能量平衡迭代上限：
  `supply_ratio < 0.2` 时最多 300 次，`supply_ratio < 0.05` 时最多 600 次。复合干旱最大能量残差从原先约 4.63 W/m² 降到
  `crop_upper≈0.40 W/m²`、`fruit_ellipsoid≈0.10 W/m²`、`pine_cone≈0.09 W/m²`。

#### 3.1 第二版（2026-04-20，LAI 比例导度分配）

**升级位置**：`src/+plant/+PlantHydraulic/solveEbalHydraulicsLayered.m`

**物理动机**：第一版假设每层 wt 等于聚合值（Phase 2），即所有层共享同一个 wt_sun / wt_sha。实际上，冠层上部辐射强、蒸腾需求高，每单位 LAI 的流量更大，叶片实际要经历更低的叶水势和更强的气孔关闭。LAI 比例导度方案将共享茎导度按各层 LAI 比例分配，使上部高需求层获得更高流量，从而解出逐层 ψ_leaf 和 wt。

**物理公式**：

茎端水势由聚合求解给出（与 v1 相同）：

$$\psi_{stem} = \psi_r - \rho_w g H / 10^6$$

每层叶水势满足：

$$K(\psi_{stem}) - K(\psi_{leaf,il}) = \frac{F_{il}}{f_{LAI,il}}$$

其中：
- $F_{il}$：第 $il$ 层实际流量（m s⁻¹），按各层需求比例从聚合实际流量中分配
- $f_{LAI,il} = \text{LAI}_{il} / \text{LAI}_{total}$：LAI 归一化权重（-）

当每层 LAI 均匀时，$f_{LAI,il} = 1/n_l$，退化为与聚合结果完全一致。

**实现**：`solveEbalHydraulicsLayered` 新增第 9 个参数 `mode`：

```matlab
% mode = 'LAI_proportional'（默认）→ Phase 3 v2
% mode = 'aggregate'           → 回退到 Phase 2 聚合标量重复
```

`ebal_clm45` 默认使用 `mode = 'LAI_proportional'`，可由 `options.use_layer_hydraulic_stress = false` 切换回聚合模式。

**验证结果**（`tests/check_layered_hydraulic_stress.m`，2026-04-20，k=108，Rin=759 W/m²）：

| 冠型 | 模式 | wt_sun_agg | wt_sun_min | wt_sun_max | wt_sha_agg | wt_sha_min | wt_sha_max |
|------|------|-----------|-----------|-----------|-----------|-----------|-----------|
| fruit_ellipsoid | LAI-proportional-v1 | 0.90 | 0.58 | 0.90 | 0.98 | 0.97 | 0.98 |
| pine_cone | LAI-proportional-v1 | 0.90 | 0.55 | 0.91 | 0.98 | 0.97 | 0.98 |
| crop_upper | LAI-proportional-v1 | 0.91 | 0.69 | 0.98 | 0.98 | 0.97 | 0.98 |

图像解读（`tests/figures/layered_hydraulic_stress_profiles.png`）：

- **wt_sun 逐层差异显著**：上部（depth≈0）蒸腾需求最强，wt_sun 最低（约 0.55–0.69），向下逐渐恢复至聚合值附近。下部层受到的水力胁迫明显小于上部层，聚合值（wt_sun_agg≈0.90，红色虚线）始终处于分层分布的中间位置。
- **wt_sha 几乎无变化**：阴叶蒸腾需求极低（接近 0），各层 wt_sha 均贴近 1.0，与聚合值（蓝色虚线）重合。这证明水力胁迫对阴叶几乎没有影响，主要从阳叶开始。
- **冠型影响**：`pine_cone` 上部 LAI 较多、需求集中，上层 wt_sun 最低（≈0.55）；`crop_upper` 上部需求同样高但 LAI 向顶部指数递减，上层虽然需求强但 LAI 较小，因此 f_LAI 较小，流量 / f_LAI 实际更大，wt_sun 最低仍在 0.69 左右；`fruit_ellipsoid` 居中。
- **对聚合通量的影响极小**（< 0.001 W/m²）：这是预期结果——Phase 3 改变了每层的 wt 分布，但 LAI 加权平均后聚合供水量不变，因此总 lE / H 与 Phase 2 基本一致。Phase 3 的价值在于**准确描述冠层哪个位置发生了水力胁迫**，为后续分层光合/蒸腾对干旱响应的精细诊断提供基础。

### 3.1 后续可进一步增强

```matlab
RWU_layer
```

`RWU_layer` 可以把当前“逐层需求驱动 stress”扩展成真正的逐层冠层-根系水力耦合。

### 4. 冠层空气剖面（Phase 4 热量/水汽/CO2 全部接入，2026-04-20 完成）

`ebal_clm45` 已新增逐层冠层空气输出：

```matlab
Ta_layer
ea_layer
CO2_layer
Kh_layer
z_layer_m
depth_lai
```

实现位置：

```text
src/+equations/calc_canopy_air_profile.m
src/ebal_clm45.m
src/resistances_clm45.m
```

当前做法：
- `resistances_clm45` 输出冠层顶热扩散系数 `Kh`、冠层顶风速 `uh` 和 Wallace-Verhoef 冠层衰减系数 `n_canopy`。
- `calc_canopy_air_profile` 用 `Kh` 和 `n_canopy` 生成逐层 `Kh_layer`。
- `Ta_layer / ea_layer` 由逐层热量和水汽源项通过 `Kh_layer` 的 steady flux-gradient 关系求得。
- `heatfluxes_clm45` 计算叶片通量时使用本层 `Ta_layer / ea_layer`；土壤通量使用冠层最底层空气。
- **`CO2_layer` 已修正并接入生化模块（2026-04-20）**，见 §4.1。
- 层序仍为 `1 = canopy top`，`nl = canopy bottom`。

叶片到空气的阻力目前仍保留原 SCOPE/CLM 有效阻力尺度，只替换空气状态为 `Ta_layer / ea_layer`。曾测试直接改成局地边界层阻力，但能量平衡温度更新过于刚性，在强辐射白天会发散。后续若要进一步物理化，应配套隐式冠层空气求解器，而不是只替换阻力。

#### 4.1 CO2_layer 公式修正与生化接入（2026-04-20）

**修正前问题**：旧公式使用叶片表面 `Cs`（近气孔 CO2 浓度）作为层内 CO2 "空气源"，物理上错误——叶片从空气吸收 CO2，空气不应以 Cs 为驱动。

**修正后公式**（与 `Ta_layer / ea_layer` 离散化方式一致）：

$$
\text{CO2\_layer}(il) = C_a - \int_{top}^{il} \frac{A_{mid}(z)}{K_h(z) \cdot \rho_{mol}} \, dz
$$

其中：
- $C_a$：冠层顶参考 CO2 浓度（μmol mol⁻¹）
- $A_{mid}$：中点处叶片净光合 + 土壤呼吸的累积通量（μmol m⁻² s⁻¹）
- $K_h$：湍流扩散系数（m² s⁻¹）
- $\rho_{mol}$：空气摩尔密度（mol m⁻³）

修正后梯度方向正确：光合作用消耗 CO2，冠层底部 CO2 浓度低于顶部。

**实现细节**：
- `air_source.A_layer_ground` 由 `ebal_clm45` 中当前迭代的 `Ahh / Auu`（GPP）逐层加权得出。
- `calc_canopy_air_profile` 中 CO2 积分与 Ta、ea 方向一致，从顶部向底部累积。
- `options.use_CO2_layer_biochem`（默认 `true`）控制 `CO2_layer` 是否作为每层光合的 $C_s$ 边界。

**生化反馈机制**：

```
CO2_layer(il) → Ca_leaf_sha(il) / Ca_leaf_sun(il)
              → heatfluxes_clm45 中的 Cs
              → vegetation.biochemical 的气孔阻力和 Ci
              → A / rcw / LE（影响整个 SCOPE 光合-气孔-蒸腾链路）
```

**修正前后对比**（k=108，Rin=759 W/m²，fruit_ellipsoid）：

| 指标 | 修正前（诊断松弛，错误 Cs 源） | 修正后（flux-gradient，正确光合源） |
|------|------|------|
| CO2 顶/底（μmol/mol） | 380 / 206（-174 梯度，物理错误） | 380 / 376（-4 梯度，合理） |
| dlEctot（clm45 - Path A） | ≈ +0 | +0.3–0.7 W/m² |
| EB 残差 | 不变 | 不变 |

**验证结果**（`tests/check_co2_layer_coupling_impact.m`，2026-04-20）：

- 修正后三种冠型（`fruit_ellipsoid / pine_cone / crop_upper`）的 CO2 垂直梯度约为 −4 至 −5 μmol/mol（从 Ca≈380 降至 375–376），物理合理。
- 开启 CO2 耦合使 lEctot 增加约 0.3–0.7 W/m²，Actot 增加约 0.003–0.01 μmol/m²/s，对 EB 残差无显著影响。
- 图 `tests/figures/co2_layer_coupling_profiles.png` 显示：三种冠型的 CO2 廓线均为单调递减（顶部 ≈380 → 底部 ≈375–376 μmol/mol），形状由冠型 LAI 分布驱动（`pine_cone` 底部 CO2 略低，因下部 LAI 集中；`crop_upper` 中部梯度较大）。

验证结果（旧测试仍通过）：
- `tests/check_lai_profile_physical_impacts.m` 已检查 `profiles.Ta_layer / ea_layer / CO2_layer / Kh_layer / z_layer_m` 存在、长度等于 `nl`、数值有限，且 `Kh_layer > 0`。
- `tests/check_lai_uniform_regression.m` 通过：显式 uniform 分层仍复现旧 fallback 路径，关键通量差异为机器精度。
- `tests/check_lai_layer_hydraulic_drought.m` 通过：复合干旱下 `crop_upper / fruit_ellipsoid / pine_cone` 的供给比仍约为 `0.039-0.044`，冠型差异保持；Phase 4 接入后严重供水受限情景的迭代上限提高到 900 次，最大能量残差约为 `crop_upper≈0.61 W/m²`、`fruit_ellipsoid≈0.10 W/m²`、`pine_cone≈0.20 W/m²`。
- 新增图：

```text
tests/figures/phase4_canopy_air_profiles.png
tests/figures/co2_layer_coupling_profiles.png
tests/figures/layered_hydraulic_stress_profiles.png
```

### 5. 根系吸水剖面输出

当前 `ebal_clm45` 里根系吸水输出仍是占位：

```matlab
RWU = zeros(size(PSIs));
```

新的植物水力求解器已经参与蒸腾供水限制，但显式土层 RWU 诊断还没有重新接入。

## 实现进度与后续计划

### 已完成（截至 2026-04-20）

| # | 任务 | 完成日期 |
|---|------|------|
| 1 | LAI 垂直分布形状（13 种 PFT）验证 | 2026-04-19 |
| 2 | `ebal_clm45` 总通量积分接入 `LAI_layer` | 2026-04-19 |
| 3 | RTMo / RTMt_sb / RTMf / RTMz 非均匀层厚改造 | 2026-04-19 |
| 4 | 分层蒸腾需求和逐层水力胁迫 Phase 3 v1（需求比例） | 2026-04-19 |
| 5 | 冠层空气剖面 `Ta_layer / ea_layer` 接入热量/水汽 | 2026-04-19 |
| 6 | Phase 3 v2：LAI 比例导度分配，逐层 ψ_leaf 和 wt | 2026-04-20 |
| 7 | `CO2_layer` 公式修正（flux-gradient 积分，正确光合源） | 2026-04-20 |
| 8 | `CO2_layer` 接入生化模块（`use_CO2_layer_biochem=true` 为默认） | 2026-04-20 |

### 待实现

| # | 任务 | 说明 |
|---|------|------|
| 9 | 显式逐层根系吸水 `RWU_layer` | 当前 RWU 为占位零向量；需重新接入土层 RWU 诊断 |
| 10 | 完整多冠层空气层 SPAC | 每层独立气孔-光合-水力闭合，需隐式冠层空气求解器 |
| 11 | 叶片局地边界层阻力物理化 | 当前阻力保留 SCOPE/CLM 有效尺度，替换需配套隐式求解 |

## 当前建议默认值

如果目标是苹果园 / 果树冠层，建议先测试：

```text
LAI_profile_mode=orchard_open_center
LAI_profile_peak=0.50
LAI_profile_sigma=0.24
LAI_profile_skew=0.8
LAI_profile_k=1.5
```

然后与以下模式对比：

```text
fruit_ellipsoid
uniform
```

`LAI_profile_skew=0.8` 是当前默认值，目的是让 `willow_drooping` 与 `fruit_ellipsoid` 区分开；如果果园开心形过度偏向冠底，可先测试 0.3-0.8 的敏感性。

## 2026-04-20 补充：分层情景验证与当前默认物理解释

### 1. 情景矩阵验证

新增脚本：

```text
tests/check_layered_scenario_transitions.m
tests/figures/layered_scenario_transition_summary.csv
```

验证了以下 7 个情景：

```text
wet_day_default
dry_day_default
wet_night_default
dry_night_default
wet_to_dry
dry_to_wet
dry_day_layer_exp
```

核心结果如下：

| 情景 | 时段 | 水力模式 | supply_ratio | wt_sun | EBmax (W m-2) |
|---|---|---|---:|---:|---:|
| wet_day_default | day | shared_path_scalar_repeated | 1.0000 | 0.9024 | 0.207 |
| dry_day_default | day | shared_path_scalar_repeated | 0.0471 | 0.2063 | 0.090 |
| wet_night_default | night | shared_path_scalar_repeated | 1.0000 | 0.9916 | 0.093 |
| dry_night_default | night | shared_path_scalar_repeated | 1.0000 | 0.7980 | 0.098 |
| wet_to_dry | day | shared_path_scalar_repeated | 0.0468 | 0.2064 | 0.099 |
| dry_to_wet | day | shared_path_scalar_repeated | 1.0000 | 0.9008 | 0.208 |
| dry_day_layer_exp | day | LAI_proportional_v1 | 0.0471 | 0.2063 | 0.096 |

物理判断：

1. 白天干旱相对湿润显著降低 `supply_ratio` 和 `wt_sun`，方向正确。
2. 夜间蒸腾需求远小于白天，模型可稳定收敛。
3. `wet_to_dry` 和 `dry_to_wet` 两类跃迁都可跑通，且最大能量残差保持在约 `0.21 W m-2` 以内。
4. 干旱夜间 `supply_ratio = 1` 但 `wt_sun < 1` 是当前模型设定下允许的结果：夜间需求很小，所以供需比不受限，但土壤与叶水势仍偏低，脆弱曲线仍给出胁迫。

### 2. 当前默认物理解释

根据 `docs/plant_hydraulic_method.md`，当前默认植物水力拓扑为：

```text
soil -> root -> shared stem node -> sunlit leaf / shaded leaf
```

因此当前主模型的默认设置是：

1. `LAI_layer` 参与辐射分层、能量积分、蒸腾需求分层、冠层空气剖面诊断。
2. 植物水力主闭合仍保持“共享上游路径 + 两个叶端分支”。
3. 默认 `layer_hydraulic_mode = shared_path_scalar_repeated`，即聚合 `wt_sun / wt_sha` 复制到各层。
4. 只有显式设置 `options.use_layer_hydraulic_stress = true` 时，才启用实验性的 `LAI_proportional_v1` 逐层 `wt` 方案。

### 3. 对逐层 wt 的解释边界

必须区分两类量：

1. **默认模式下的逐层 wt**
   - 只是共享路径聚合解在层维度上的重复输出；
   - 不会改变原 two-leaf 水力理论。

2. **实验模式下的逐层 wt**
   - 反映局部叶层在“需求/LAI 比值”下的相对受限程度；
   - 可用于诊断“哪一层更容易先受旱”；
   - 不能直接当作新的整株主风险指标，也不能替代原 `wt_sun / wt_sha` 的共享路径解释。

例如在 `dry_day_layer_exp` 中：

```text
aggregate wt_sun = 0.2063
top20 wt_sun     = 0.2026
bottom20 wt_sun  = 0.2120
```

这只能解释为“顶部局部更受限、底部略轻”，不能解释为整株水力主风险已经从 `0.2063` 改写成别的值。

### 4. 当前阶段结论

截至 2026-04-20，分层方案的结论是：

1. **有用且已生效的部分**
   - LAI 垂直剖面
   - 辐射分层
   - 能量与光合总量积分
   - 分层蒸腾需求
   - 冠层空气 `Ta_layer / ea_layer / CO2_layer / Kh_layer`

2. **仍保持保守默认的部分**
   - 植物水力主路径
   - 整株供水上限与 `wt` 主控制

3. **实验诊断用途**
   - `LAI_proportional_v1` 的逐层 `wt`
   - 用于敏感性分析，不作为默认主物理
