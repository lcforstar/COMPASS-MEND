# COMPASS 变饱和水流设计与实现

日期：2026-07-12
状态：water-only 新链已实现；水热全耦合验证待下一增量
适用范围：自然点位、实验土柱以及未来区域单元共用的 Richards 内核

## 1. 最终决策

不再采用“把剖面手工拆成饱和求解器和非饱和求解器”的旧 Route B
草案，也不把内部饱和区当成用户边界。最终方案是在同一个 mixed Richards
方程中自动处理变饱和状态：

\[
\Theta_u(h)=
\begin{cases}
\Theta_{VG}(h), & h<0,\\
\theta_s+S_s h, & h\ge 0,
\end{cases}
\]

\[
C_h=\frac{\partial\Theta_u}{\partial h}=
\begin{cases}
C_{VG}(h), & h<0,\\
S_s, & h\ge 0,
\end{cases}
\qquad
K(h)=
\begin{cases}
K_{MVG}(h), & h<0,\\
K_s, & h\ge 0\ \text{且未冻结}.
\end{cases}
\]

用户仍然只配置顶部和底部外边界。剖面内部的饱和区、非饱和区及
`h=0` 界面全部由当前迭代解自动决定。

这与 ClimaLand 的增广含水量思想一致：其
`soil_hydrology_parameterizations.jl:121-125` 在饱和分支使用
`psi=(vartheta-nu_eff)/S_s`，反解即 `vartheta=nu_eff+S_s*psi`；
`lines 149-150` 给出饱和导数 `dpsi/dvartheta=1/S_s`；
`lines 167-172` 在饱和时返回 `K=K_sat`。

## 2. 每个非线性迭代的 active-set 过程

1. 从已接受状态构造当前 Picard/Anderson 试值 `h^k`。
2. 逐端点判断 `h^k<0` 或 `h^k>=0`。
3. 使用对应的 `Theta_UU`、`DThetaWaterEquivDh` 和 `KfL_h` 组装同一
   Richards 方程。
4. 解得 `h^(k+1)`，保留所有正水头，不做全局 `-1e-6 cm` 截断。
5. 若任一节点跨过 `h=0`，active set 发生变化，本迭代不得接受；用新
   active set 再迭代。
6. 非饱和端点检查最大 `Delta theta`，饱和节点检查最大 `Delta h`，并用
   实际边界通量和内部根汇检查质量残差。
7. 只有三道门同时通过才接受；拒绝步不得提交 `h`、水相、时钟或诊断。

## 3. 冻融和总水储量

`Theta_UU` 是唯一的凝结态液态水当量储量，包含：

- 液态水；
- `rho_i/rho_l` 加权的冰水当量；
- 正压力储水增量 `S_s h`。

必须区分“总质量储量”和“真实孔隙相态”。在 (h\ge0) 时：

\[
\Theta_{UU}=\theta_p+S_sh,\qquad
\theta_p=\Theta_{LL}+\frac{\rho_i}{\rho_l}\Theta_{II}\le\theta_s.
\]

因此 `Ss*h` 只进入 Richards 储量与水量平衡，不进入液水、冰和气相体积。
这避免热方程把弹性储水误当成液水热容或冰相潜热。

mixed 时间项和 WBC 都使用：

\[
\theta_{eq}=\Theta_{UU}+\theta_v\frac{\rho_v}{\rho_l}.
\]

因此不能再用 `Theta_LL + Theta_II*rho_i/rho_l` 在诊断端重建总水；那会
遗漏正压力储量。相变可以改变液态和冰态分配，但不能改变同一时刻的
凝结态总质量。

已接受相状态现在由 `Y.Theta_U` 和 `Y.Theta_I` 保存。每个非线性迭代的
相变松弛都锚定同一个步初 `Y.Theta_I`，避免一次时间步内因多次 Picard
刷新而重复执行相变。接受后才提交新水相，拒绝后保持原相状态。

## 4. 严格水量平衡

统一土柱水量平衡为：

\[
R_W=(-QMT)\Delta t-Q_{root}\Delta t-(-QMB)\Delta t-\Delta W.
\]

- `-QMT`：求解后实际顶部入流，适用于定水头、定通量和自然边界；
- `-QMB`：实际底部出流；
- `Q_root`：对本次送入 Richards 的 `Srt` 做相同 FE 端点积分；
- `Delta W`：前后两端都按 `Theta_UU + vapor` 和相同 FE 求积计算。

`Infiltration_act` 和 `EVAP` 继续作为过程分量输出，但不再代替定水头的
实际顶部通量，也不在已包含净蒸发的 `QMT` 之外重复扣除蒸发。

默认精度：

- 长期自然模拟：`tol_water_mass_mm=1e-4 mm/step`；
- HYDRUS/ClimaLand 实验对标：`tol_water_mass_mm=1e-6 mm/step`。

方程、通量和守恒尺子完全相同，只有用户指定的验收精度不同。

## 5. h=0 界面输出

`computeZeroHeadDiagnostics` 对接受后的分段线性压力水头求全部 `h=0`
交点。它允许多个交点，因此可表示滞水层，而不是强制只有一个地下水面。

核心输出包括：

- `soil_saturated_node_mask`；
- `zero_head_crossing_depth`，从浅到深并用缺测值补齐；
- `zero_head_crossing_count`；
- `shallowest_zero_head_depth`；
- `saturated_profile_fraction`。

SV62 benchmark 另外保存全部压力水头剖面和
`zero_head_crossings.csv`，可直接绘制移动的零水头线。

## 6. 当前代码位置

- 变饱和本构：`src/+soil/+water/computeEndpointHydraulics.m`
- 总水相导数：`src/+soil/+water/computePhaseStoragePartials.m`
- mixed 储量：`src/+soil/+water/calculateTimeDerivatives.m` 的局部储量函数
- 保守系数：`src/+soil/+water/calculateWaterEquationCoefficients.m`
- 正水头提交：`src/+soil/+water/solveTridiagonalMatrixEquations.m`
- 实际通量 WBC：`src/+solver/computeWaterBalanceResidual.m`
- 根汇积分：`src/+solver/integrateRootSinkCmS.m`
- 零水头诊断：`src/+soil/+water/computeZeroHeadDiagnostics.m`
- 站点精度设置：各站点 `site.toml` 的 `[numerics]`

## 7. 已验证与未完成

已验证：

- SV62 顶部 `+5 cm` 在接受后仍精确为 `+5 cm`；
- `theta=theta_s+S_s h=0.46705000`；
- 饱和容量为 `S_s`，暖饱和导水率为 `K_s`；
- 30 s 定水头实际通量 WBC 残差为约 `-1.08e-7 mm`；
- 冻融液态/冰态闭合到孔隙水 `theta_p`，弹性储水单独保留在 `Theta_UU`；
- 固定温度的真实冻结 Richards 步已通过，WBC 残差约 `-1.16e-9 mm`；
- 多个 `h=0` 交点、滞水层和拒绝回滚；
- site/schema、初始化、NetCDF 以及 `runOneStep 109/109` 回归。

未完成：

- 水热全耦合矩阵在这些新物理导数上的整链对标；
- 饱和区气相体积归零后的气压方程处理；
- 完整大气互补边界、积水状态机、渗出面和动态地下水边界；
- 全部 W2 长案例与 HYDRUS/ClimaLand 图形比较。

后续顺序保持：先完成水热方程，再改气体方程。
