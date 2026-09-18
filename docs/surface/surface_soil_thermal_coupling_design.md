# COMPASS 地表–土壤热耦合早期设计（已由现役统一响应方案取代）

_历史设计冻结日期：2026-07-18。本文记录最初的 A0/A1 论证，不再作为实现或输出契约。现役权威说明见 `surface_soil_thermal_coupling_progress.md`、`partial_snow_two_tile_implementation_plan.md`、`output_catalog.md` 和 `compass_soil_water_heat_physics.md`；生产路径已经统一为裸地、部分雪和全雪共用的 `surface_energy_response`，旧输出名也已删除。_

## 1. 设计结论

第一优先固定为两步，且必须串行完成：

1. **A0：先建立无歧义的地表–土壤总能量账本。** 保留土壤顶部完整总内能反作用通量，拆出导热、液水、水汽和干空气携能，并把 SCOPE 的净辐射、显热、实际潜热、降水焓与土壤储能放进同一个接受步控制体。
2. **A1：再做无雪 `Tss–G_soil` 有界校正。** SCOPE 继续求裸土地表温度，土壤继续接收 Dirichlet `Tss`；土壤解出的总能量反作用用于同一步的地表校正。默认按需触发且最多增加一次土壤求解，不在两个模块之间开启无限 Picard。

COMPASS 的最终方法不是 Noah-MP、CLM5 或 ClimaLand 的复制，而是适配当前三方程结构的：

> **总能量反作用约束下的事务式分区校正（transactional bounded partitioned correction）。**

## 2. 当前接口及必须先澄清的量

### 2.1 当前无雪调用关系

当前接受步先运行 SCOPE 裸土地表能量平衡，再以其 `Tss` 作为土壤顶部温度边界，完整求解 `h/U/P_g`。土壤返回顶部热方程反作用 `QET`，但本步不再反馈给 SCOPE。因此问题是同一步接口滞后，不是土壤三方程内部失耦。

### 2.2 `G` 与 `QET` 的权威含义

SCOPE 中：

\[
R_{n,g}-H_g-LE_g-G=0,
\]

其中 `G` 为进入地表下方控制体的**净能量通量**，正值向下，单位 W m\(^{-2}\)。现有 force-restore 计算的 `G` 及 `dG_dTs` 只能作为第一次预测。

土壤热方程的公开边界约定为进入土壤正。当前内部 `QET` 顶部向外为正，因此实际总能量反作用是：

\[
\boxed{G_{\mathrm{soil,total}}
=q_{U,\mathrm{top}}^{\mathrm{in}}
=-10^4 QET}
\qquad [\mathrm{W\,m^{-2}}].
\]

它不是纯 Fourier 导热，而是：

\[
q_{U,\mathrm{top}}^{\mathrm{in}}
=q_{U,\mathrm{cond}}+q_{U,l}+q_{U,v}+q_{U,a}.
\]

因此 A1 的反作用目标应是**完整的顶部总能量通量**。不能把 `QET` 标成纯导热；也不能只把 Fourier 导热送回 `EBers`，否则土壤内的水汽、液水和空气携能会与地表潜热或降水焓失去统一控制体。

纯导热分量仍需输出，用于解释为什么总反作用与 force-restore 不同，但不取代总反作用成为 A1 权威量。

## 3. 从参考模型吸收什么

| 参考 | 借鉴点 | COMPASS 不照搬的部分 |
|---|---|---|
| Noah-MP | 地表 Newton 更新把土壤热导纳 `CGH` 放进分母；地表能量显式包含降水热；地温迭代次数有限 | 不改成 Noah 的多层温度方程，不使用固定 5 次循环 |
| CLM5 | 顶层热方程使用地表能量通量及其温度导数，接口采取半隐式线性化 | 不把当前 U-primary 土壤改成 CLM 的雪–土统一温度矩阵 |
| ClimaLand | 一个权威雪/土接口通量；两侧反号；水和能量储量与边界通量统一验收 | A0/A1 无雪阶段不引入其雪层结构和时间积分器 |
| 当前 SCOPE | 保留 `EBers`、裸土地表温度和现有 `dG_dTs` force-restore 热导纳 | force-restore `G` 降级为预测，不再作为接受步事实 |

参考依据：[Noah-MP 技术说明](https://www2.mmm.ucar.edu/wrf/users/physics/phys_refs/LAND_SURFACE/noah_mp_tech_note.pdf)、[CLM5 地表通量](https://escomp.github.io/CTSM/release-clm5.0/tech_note/Fluxes/CLM50_Tech_Note_Fluxes.html)、[CLM5 土壤与雪温度](https://escomp.github.io/CTSM/release-clm5.0/tech_note/Soil_Snow_Temperatures/CLM50_Tech_Note_Soil_Snow_Temperatures.html)、[ClimaLand 雪边界耦合](https://clima.github.io/ClimaLand.jl/stable/generated/integrated/handling_snow_fluxes/)。本地源码对照分别位于 `D:\OneDrive\应用\Model\NoahMP_Dual-main`、`D:\OneDrive\应用\Model\CTSM-release-clm5.0.37` 和 `D:\OneDrive\应用\Model\ClimaLand`。

## 4. A0：真正的地表+土壤总能量账本

### 4.1 统一符号和单位

公开诊断全部使用：

- 能量通量：W m\(^{-2}\)；**进入所述控制体为正**；
- 步累计能量：J m\(^{-2}\)；
- 温度：°C，温差：K；
- 土壤旧单位 `J cm^-2 s^-1` 仅在适配器内部出现，写出前乘 `1e4`；
- 每个量同时记录 `requested/predictor` 与 `actual/reaction`，两者不得共用字段。

### 4.2 三本账，而不是一个模糊残差

#### 土壤柱账

沿用当前精确储量定义：

\[
R_{\mathrm{soil}}
=\Delta E_{\mathrm{soil}}
-\left(G_{\mathrm{soil,total}}+G_{\mathrm{bottom}}\right)\Delta t
+E_{\mathrm{root}}.
\]

其中顶部、底部均以进入土柱为正，根携能为离开土柱的正量。该账验证 U-primary 土柱本身。

#### 零储能裸土地表账

无雪且暂无独立皮肤储能时，地表控制体的净外部输入为：

\[
F_{\mathrm{surface,ext}}
=R_{n,g}-H_g-LE_{g,\mathrm{actual}}
+H_{\mathrm{precip}}+H_{\mathrm{mass,other}}.
\]

表面账按“储量变化减净输入”定义：

\[
R_{\mathrm{surface}}
=\left(G_{\mathrm{soil,total}}
-F_{\mathrm{surface,ext}}\right)\Delta t.
\]

这里必须使用水活动集接受后的 `LE_actual`，不能在干旱受限时继续用 SCOPE 潜在蒸发。`H_precip` 使用到达地面的液态降水率和明确的降水温度；初版没有独立雨温观测时，以气温作为带标志的回退，不允许静默省略。`H_mass,other` 用于明确登记干空气等跨顶边界携能，A0 先量化后决定是否可忽略。

#### 地表+土壤总账

两账相加后，内部接口 `G_soil,total` 必须严格抵消：

\[
\boxed{
R_{\mathrm{total}}
=\Delta E_{\mathrm{soil}}
-F_{\mathrm{surface,ext}}\Delta t
-G_{\mathrm{bottom}}\Delta t
+E_{\mathrm{root}}
}
\]

并要求数值上：

\[
R_{\mathrm{total}}
=R_{\mathrm{surface}}+R_{\mathrm{soil}}
\]

在舍入容差内成立。总账不再次加入 `QET`，否则接口能量会双计。

### 4.3 顶部总通量分解

A0 必须同时给出：

\[
R_{\mathrm{decomp}}
=G_{\mathrm{soil,total}}
-\left(G_{\mathrm{cond}}+G_l+G_v+G_a\right).
\]

总反作用仍由完整耦合边界行恢复；四个机制分量从同一个接受状态、同一个顶端面和同一个时间平均规则恢复。不能把内部单元中点诊断直接冒充边界反作用。

### 4.4 A0 输出契约

建议归入独立的 `surface_energy.nc` 与现有 `conservation.nc`，字段至少包括：

| 字段 | 含义 |
|---|---|
| `ground_energy_flux_predictor` | force-restore 第一次预测，向下为正 |
| `soil_top_energy_flux_total_actual` | `-QET*1e4`，接受后的完整总能量反作用 |
| `soil_top_energy_flux_conductive_actual` | 顶端 Fourier 导热分量 |
| `soil_top_energy_flux_liquid_actual` | 液水携能分量 |
| `soil_top_energy_flux_vapor_actual` | 水汽显热与汽化能分量 |
| `soil_top_energy_flux_dry_air_actual` | 干空气携能分量 |
| `soil_top_energy_flux_decomposition_residual` | 总反作用减四分量之和 |
| `ground_latent_heat_flux_demand` | SCOPE 潜在裸土潜热 |
| `ground_latent_heat_flux_actual` | 水活动集接受后的实际潜热 |
| `precipitation_energy_flux` | 液态穿透雨携带的显热/焓 |
| `surface_external_energy_flux` | `Rn-H-LE_actual+Hprecip+Hmass` |
| `surface_energy_power_mismatch` | `F_surface,ext-G_soil,total` |
| `surface_energy_residual` | 地表步累计残差 |
| `soil_energy_residual` | 土壤柱步累计残差 |
| `surface_soil_total_energy_residual` | 合并控制体步累计残差 |
| `surface_soil_interface_cancellation_error` | 两本账相加与直接总账之差 |

另需输出 `scope_ground_ebal_residual`，保留 SCOPE 自身使用预测量时的代数闭合误差。它与使用实际土壤反作用后的 `surface_energy_power_mismatch` 不是同一个量。

### 4.5 A0 完成门槛

A0 只加诊断，不改变任何状态或边界。完成必须同时满足：

1. 模式 0 与当前全年基线状态和既有输出在数值容差内一致；
2. 顶部总能量反作用与四分量之和闭合；
3. 地表账、土壤账、直接总账的符号和单位测试通过；
4. 干旱蒸发受限时，`LE_demand-LE_actual` 可见且不会被 `G` 静默吸收；
5. 降雨时 `H_precip` 可见，降雨温度回退状态可追踪；
6. 接口抵消误差只允许达到舍入量级；
7. 全年报告分别给出累计残差、最大单步残差和按过程分组的能量贡献。

## 5. A1：无雪 `Tss–G_soil` 有界校正

### 5.1 求解问题

对同一个外层时间步，把土壤视为给定地表温度后返回总能量反作用的 Dirichlet-to-Neumann 映射：

\[
G_{\mathrm{soil,total}}=\mathcal D(T_s).
\]

需要求解的标量接口方程为：

\[
\boxed{
F(T_s)
=F_{\mathrm{surface,ext}}(T_s)
-\mathcal D(T_s)=0
}.
\]

每次计算 `D(Ts)` 时，土壤内部仍完整求解 `h/U/P_g`，水边界活动集、冻融、蒸汽和干空气都不从热校正中拆开。

### 5.2 默认算法：按需触发的一次校正

每个无雪外层步执行：

1. **冻结步初快照。** 保存接受态的冠层、地表、土壤、时间层和账本；所有候选都从该快照推进，不能重复累计降水、截留、蒸发或植物状态。
2. **预测。** 运行现有 SCOPE/`EBers`，用 force-restore `G_FR` 和 `dG_FR/dTs` 得到 `T_s,0`。
3. **第一次完整土壤解。** 以 `T_s,0` 为顶部 Dirichlet 温度，求解完整 `h/U/P_g`，得到实际 `G_0=D(T_s,0)`、实际蒸发、水边界状态和 A0 三本账。
4. **按需判定。** 若 `F(T_s,0)`、总账、三方程守恒和边界状态均已满足接受门槛，直接提交，不增加第二次土壤求解。
5. **构造有锚点的土壤反作用。** 若需校正，使用

   \[
   \widehat G_0(T_s)
   =G_0+C_0(T_s-T_{s,0}),
   \]

   而不是把 `G_0` 冻成常数。第一次斜率 `C0` 由 A0 选定的热导纳预测器提供，并受物理上下界约束。
6. **重算地表能量。** 从相同步初 SCOPE 状态重跑能量平衡，只重算温度相关的地表/冠层能量项，不重跑不依赖温度的光学传输；在 `EBers` 分母中使用 `C0`。得到受限的 `T_s,1`。
7. **第二次完整土壤解。** 从同一个接受态时间层以 `T_s,1` 重解 `h/U/P_g`，得到 `G_1=D(T_s,1)` 和新的实际水/能量通量。
8. **联合接受。** 只有 A0 总账、三方程守恒、水活动集稳定、`F(T_s,1)`、温差和反作用线性化误差全部通过，才一次性提交冠层、地表和土壤候选。

默认模式最多执行步骤 6–7 一次。任何未收敛候选都可丢弃；只有完整收敛的旧候选可以作为新候选的数值暖启动，且不得改写固定的步初时间层。

### 5.3 热导纳斜率

第一次校正不新增有限差分土壤求解。候选斜率为：

1. 当前 SCOPE force-restore 的 `dG_FR/dTs`；
2. 顶层离散热导纳估计，由顶层导热与该时间步有效储热共同给出；
3. 二者经 A0 离线有限差分试验比较后选定或有界混合。

约束固定为：

- `C0 >= 0`；
- 上界来自当前顶层热导率、网格尺度、热容量和时间步，不使用无来源常数；
- 若得到第二组完整反作用点，可计算

  \[
  C_{\mathrm{sec}}=
  \frac{G_1-G_0}{T_{s,1}-T_{s,0}},
  \]

  经正值、幅值和小分母保护后，仅供模式 2 的最后一次校正使用。

A0 的专项斜率 oracle 在选定时刻从同一快照分别求解 `Ts±δT`，用于评价斜率，不进入年度生产路径。

### 5.4 干旱活动集与实际潜热

当水边界为供水/蒸发受限状态时：

- 若实际蒸发等于 SCOPE 需求，沿用 SCOPE 的 `dLE/dTs`；
- 若实际蒸发受水量限制，校正方程使用实际潜热，并把其局部温度导数设为当前活动支可支持的值，初版可在同一活动支内取零；
- 若第二次土壤解改变了干湿活动支，旧线性化失效，必须重新判定或进入模式 2，不能让 `G` 吸收潜热差额；
- 液态入渗/蒸发通道始终保留，由水活动集自动选择，不设置人工“完全关闭液水”的热校正分支。

### 5.5 有界、回滚和失败策略

每次温度更新必须同时受四类边界约束：

- 最大 `|ΔTss|` 信赖域；
- `Tss` 物理范围与 U→T 可逆范围；
- 非负且有上界的热导纳；
- 新候选必须降低归一化接口残差，否则采用阻尼或区间内回退。

外层不无限循环：

| 模式 | 行为 |
|---:|---|
| `0` | 只运行 A0；保持当前单向路径 |
| `1` | 按需最多一次校正；生产默认候选 |
| `2` | 最多两次校正，第二次使用有保护的 secant/阻尼；用于强变天气和诊断 |

达到次数上限后：

1. 若存在一个完整收敛候选同时通过所有硬守恒门槛，只是未达到更严格的接口迭代软阈值，可提交残差最小者并标记 `bounded_soft_accept`；
2. 若总能量、水、干空气或非线性硬门槛不通过，整个候选事务失败，返回明确原因；不得提交半校正状态；
3. 土壤非线性失败沿用自适应 substep 缩步；若只是接口残差失败，先使用模式 2 的有界校正，不能用无脑 Picard 重试掩盖。

硬/软阈值必须由 A0 的舍入误差、现有能量容差、时间步敏感性和全年分布确定，设计阶段不写死经验数值。

## 6. 运行速度约束

A1 不以“每小时必做两遍”为默认：

- 第一次 SCOPE 与土壤解就是当前必需工作；
- A0 在第一次土壤解后计算接口残差；满足门槛则零额外土壤求解；
- 生产模式每个触发步最多增加一次土壤求解；
- 模式 2 最多再增加一次，只在强变或诊断时启用；
- 缓存 RTMo、光学和所有与 `Tss` 无关的量，只重算温度相关 ebal；
- 接受态快照不可变，但允许使用**已收敛**候选作为非线性初猜；未收敛候选绝不继承；
- 逐步记录 `soil_solve_count`、`ebal_solve_count`、是否触发校正和耗时。

若校正触发比例为 `p`，土壤主成本的年度上界近似从 `1` 增至 `1+p`，而不是无条件翻倍。实现验收必须报告触发比例和年度墙钟时间；在同一机器、同一配置、同一输出级别下与模式 0 比较，不能仅凭代码结构宣称比其他陆面模型快。

## 7. 验证矩阵

### A0 行为不变试验

- 无雨无雪白天与夜间；
- 强辐射裸土；
- 降雨入渗；
- 干旱蒸发受限；
- 冻融但无雪；
- 热、气开关组合；
- 模式 0 全年与当前 `8760/8760` 基线对照。

### A1 判别试验

- `G_FR` 已接近土壤反作用：不触发第二次土壤解；
- `G_FR` 明显偏离：一次校正降低接口残差和总账残差；
- 斜率过大/过小：信赖域与斜率夹紧生效；
- 干旱支：实际潜热替代需求潜热，不出现假能量；
- 降雨支：降水焓显式进入总账；
- 水边界同 dt 换支：每个热候选内的水活动集独立稳定；
- 候选土壤失败：回到不可变步初状态，无状态重复推进；
- 强变天气：模式 2 有界结束，无无限外层迭代。

### 完成判据

1. A0 的总通量分解、土壤账、地表账和合并账同时通过；
2. A1 模式 0 保持基线；
3. 模式 1 不增加非有限状态、热边界误切换或水/气守恒失败；
4. 校正后接口功率残差和总账残差相对模式 0 有系统性下降；
5. 年度额外成本由实际触发比例解释，单步求解次数严格有界；
6. 无雪 A0/A1 稳定前，不进入 10 m 初始化修正、全雪强制测试或 Stage B 两 tile 实现。

## 8. 实现落点（设计后的下一步）

预计最小改动面：

- `+soil/+thermal`：恢复接受态顶部总能量分解，保留 `QET` 反作用为权威总量；
- `+solver`：新增地表+土壤能量账本、热候选事务和联合接受状态；
- `+plant/+energy`：允许 `EBers` 接收有锚点的 `G(Ts)` 与导数，以及实际潜热活动支；
- `runOneStep`：只负责候选编排、回滚和一次性提交；
- `+diagnostics`：注册 A0/A1 字段并更新输出目录；
- `tests`：先完成 A0 行为不变和符号单位测试，再启用 A1 行为测试。

实现前仍需先写一个只读 A0 探针，确认顶部四分量能按与 `QET` 相同的端点、时间平均和参考能量恢复；若该恒等式不闭合，不得直接开始 A1。
