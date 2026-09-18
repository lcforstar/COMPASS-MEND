# MODFLOW 外部地下水边界契约

## 所有权

MODFLOW/BMI 是地下水位、饱和带水储量、地下水温度和地下水溶质的唯一所有者。COMPASS 只拥有水位以上包气带的 `h/U/P_g` 与土壤 N 状态，并向外发布实际交换 ledger。

## 冻结的 trial 输入

每个 Richards trial 从外部驱动冻结：

```text
active_node = water-table node (bottom-up)
h_gw        = fixed matric head at active_node
T_gw        = external groundwater temperature, if supplied
```

水边界为 `h(active_node)=h_gw`。节点以下的静水延拓不是 Richards 自由度，也不进入 Richards 水储量、根吸水或水量账本。

## 水量

水力反应必须在 `active_node` 恢复，而不是物理列节点 1。发布的方向为：

```text
QMB > 0: 外部地下水进入包气带
soil_to_groundwater.water = -QMB
```

水量验收使用水位以上有限元：

```text
Delta W_unsat = q_top_in - S_root - q_to_groundwater
```

## 热量（待实现）

热方程必须使用同一个 `active_node`。若 `T_gw` 有限，则在该节点施加内部 Dirichlet：

```text
T(active_node) = T_gw
```

热矩阵、储能积分、热反应恢复和 Thomas 活动范围均为 `[active_node, NN]`。地下水侧携热只可由同一 trial 的实际 `QMB` 计算；不得从旧节点 1 的通量或经验 force-restore 通量替代。

若 `T_gw` 缺失，热接口必须显式报告 `water_only_external_boundary`，不能静默把水位以下列继续当作自由包气带热域。

## 气体与氮（待实现）

水位以下没有包气带气相自由度。气压、O2 与氮不能在该段继续按普通土壤方程推进。要么截断至水位界面并给定外部边界，要么保持水力-only MODFLOW 模式并明确禁用外部地下水段的气/溶质耦合。

没有 MODFLOW 的地下水溶质浓度或通量时，禁止发布 `soil_to_groundwater.nitrogen` 为物理闭合结果。

## 验收

1. 固定水位节点和下方静水关系在接受后精确成立；
2. 水位以上水量账本闭合到 `tol_water_mass_mm`；
3. 内部温度边界实现后，水位以上能量账本和 `QMB` 携热项闭合；
4. rejection 后，`Y`、外部驱动快照和所有已接受 ledger 均不变；
5. Julia 参考窗口逐子步比较活动节点、QMB、h/U/P_g、接受决定及 ledger。
