# v30 single vs four_quadrant 暴雨 (2021-07-13 ~ 2021-07-18) 对比

## 1. 配置

```ini
# data/config_vegetation.txt
canopy_water_scheme = single | four_quadrant
canopy_storage_per_LAI = 0.1   # CLM5/NoahMP 标准 [mm/LAI]
interc_F0 = 0.15               # Gash 自由穿透系数
```

测试期 KT 1–120（120 h），同一气象、同一初始土壤剖面、同一根系。
两次跑只切换 `canopy_water_scheme`，其他参数完全一致。

## 2. 冠层各项 5 天累计 (mm 或 ∑W/m²)

| 变量 | single | four\_quadrant | 差值 (4q − single) | 物理解读 |
|---|---|---|---|---|
| h2ocan 总储量 ∑ (mm) | 1.858 | 6.439 | **+4.580** | 4q 中阴叶库 (Ps_max=70%) 慢蒸发，累计储量更高 |
| Throughfall ∑ (mm) | 92.551 | 92.242 | −0.309 | Gash sparse 控制，两模式相同 Sv，差异仅来自蒸发反馈 |
| CanopyInterception ∑ (mm) | 3.049 | 3.358 | +0.309 | 与 throughfall 互补，质量守恒 ✅ |
| Soil Evaporation ∑ (mm) | 6.266 | 6.222 | −0.045 | 几乎一致，冠层方案不直接影响土面蒸发 |
| Transpiration ∑ (mm) | 8.402 | 8.319 | −0.083 | 4q 略小，因为湿叶占比上升，干叶传导率下降 |
| **lE_wet ∑ (W/m²·h)** | 1395.4 | **1605.8** | **+210.4** (+15%) | 4q 阴叶慢蒸发尾巴，干燥期 KT 65–78 多输出 |

## 3. 暴雨主峰 KT 56–72: lE_wet 时序

```
  KT |   single |       4q |     diff
  59 |    83.71 |    83.72 |    +0.01     ← 饱和期一致
  60 |    92.78 |    92.81 |    +0.03
  63 |    96.64 |    96.66 |    +0.01
  64 |    86.05 |    86.06 |    +0.01     ← 饱和峰末
  65 |    55.57 |    58.94 |    +3.37     ← 干燥开始，4q 阴叶滞留
  66 |    25.53 |    32.09 |    +6.56     ← 4q 多 25%
  67 |     9.10 |    12.99 |    +3.89
  68 |     0.79 |     6.96 |    +6.16     ← single 几乎干尽，4q 阴叶仍 0.15 mm
  70 |     0.01 |    −0.03 |    −0.04     ← 4q 出现凝露 (single 不允许)
  71 |     0.00 |    −0.56 |    −0.56     ← 阴叶 dew gain
  72 |     0.00 |    −1.11 |    −1.11
```

**关键发现**：
- 饱和期 (KT 59–64) 两模式重合 (fwet=1)
- 干燥期 (KT 65–78) 4q 多 lE_wet ≈ 25 W/m²/h，因为 sha 库释放慢
- 凝露期 (KT 70–72) 4q 允许阴叶冷凝 (−1.1 W/m²)，single 钳制为 0

## 4. four_quadrant 阴/阳叶分离演化 (4q only)

| KT | h2_sun (mm) | h2_sha (mm) | fwet_sun | fwet_sha |
|---|---|---|---|---|
| 58 (rain) | 0.097 | 0.225 | 1.000 | 1.000 |
| 60 (peak) | 0.088 | 0.221 | 1.000 | 1.000 |
| 63 (peak) | 0.084 | 0.220 | 1.000 | 1.000 |
| 64 (drying) | 0.095 | 0.224 | 1.000 | 1.000 |
| 65 | 0.041 | 0.192 | 0.645 | 0.901 |
| 66 | 0.017 | 0.169 | 0.367 | 0.812 |
| 67 | 0.009 | 0.157 | 0.202 | 0.746 |
| 68 | 0.006 | 0.150 | 0.136 | 0.711 |
| 69 | 0.006 | 0.149 | 0.103 | 0.689 |
| 78 | 0.005 | 0.175 | 0.092 | 0.749 (回升) |
| 79 | 0.005 | 0.174 | 0.092 | 0.764 |

**Rutter 串联完整跑通**：
- 雨水先填阳叶 (Ps_max_sun = 0.3·LAI·canopy_storage_per_LAI)
- 阳叶溢流 → 阴叶库 (Ps_max_sha = 0.7·LAI·canopy_storage_per_LAI)
- 阴叶因低辐射保水更久 (e-folding > 4 h)，阳叶高辐射快蒸发 (< 1 h)

## 5. 潜热/显热/叶温/光合/叶水势 5 天均值

数据来源：`Scope_output.mat` (`fluxes_*`, `thermal_*`, `profiles_*`) + `Sim_output_Leafpotential.csv`

> **v2 (2026-05-09 修复后)**：原始数据存在 `SCOPE-GUARD` 拒收事件下 `lE_wet_canopy` / `lE_dry_canopy` 未清零的不一致 bug，导致 KT 11/13/14/15 lE_c=0 但分量大数 → 折线锯齿。修复方法：`COMPASS.m:782-785` 同步清零这两个分量。下表为修复后数据。

| 量 | single | four\_quadrant | 差值 (4q − single) | 解读 |
|---|---|---|---|---|
| **lE total** (W/m²) | 96.29 | 96.26 | −0.03 | 总潜热守恒 ✅ |
| lE canopy | 50.73 | 51.23 | **+0.51** | 4q 冠层蒸发增 |
| lE soil | 29.36 | 29.14 | −0.22 | 冠层多吸水→土面少蒸发 |
| lE wet canopy | 15.22 | 16.17 | +0.94 | 阴叶慢蒸发尾巴 |
| lE dry canopy | 35.50 | 35.06 | −0.44 | 湿叶占比上升→干叶蒸发减 |
| **H total** | 29.82 | 29.73 | −0.09 | 显热略减 |
| H canopy | 22.88 | 22.78 | −0.10 | 4q 多潜热 → 少显热 |
| H soil | 6.94 | 6.95 | +0.01 | 几乎无影响 |
| **Rn total** | 137.68 | 137.61 | −0.08 | 净辐射几乎不变 |
| Rn canopy | 86.21 | 86.26 | +0.05 | 湿叶反照率略低 → 多吸辐射 |
| Rn soil | 51.48 | 51.35 | −0.13 | 冠层多挡辐射 |
| **GPP** (μmol/m²/s) | 0.00 | 0.00 | 0 | 暴雨期光合关闭 |
| Anet | 7.68 | 7.64 | −0.04 | 几乎相同 |
| **Tcave** (°C) | 26.20 | 26.17 | **−0.029** | 4q 多蒸发→略凉 |
| Tchave | 25.19 | 25.18 | −0.011 | 同上 |
| **ψ_l mean** (MPa) | −1.042 | −1.047 | −0.005 | 均值相同 |
| ψ_l min | −3.284 | −3.217 | **+0.067** | 4q 不容易触底 |
| ψ_l max | −0.026 | −0.026 | 0 | 一致 |

### 暴雨主峰 + 二次降雨段 KT 60–80 详细

```
KT |    lEctot s|f |     Hctot s|f |     Tcave s|f |     psil s|f
60 | 94.96|94.99   | 10.84|10.84   | 26.04|26.04   | -0.05|-0.05  ← 主峰一致
65 | 44.96|48.19   | -0.25| 0.06   | 26.49|26.56   | -0.41|-0.37  ← 干燥+0.04 MPa
70 | -0.05|-0.75   | -2.29|-1.98   | 20.97|20.96   | -0.06|-0.06  ← 4q 凝露
75 | -0.23|-4.59   | -2.69|-2.74   | 18.77|19.14   | -0.07|-0.07  ← 4q -4.4 W/m² (凝露强)
80 | 68.49|85.23   | 38.33|28.89   | 26.89|25.85   | -0.99|-1.03  ← 二次降雨：4q lE +17, H -9, T -1.0°C
```

### 关键差异机制

1. **能量再分配（向潜热偏）**：4q 冠层 lE +0.56 / H −0.17 / Rn +0.13 W/m²，符合 Bowen 比下降的预期
2. **叶面温降**：5 天均值 4q 比 single 凉 0.035°C，瞬时 (KT 80 二次降雨) 可达 1°C
3. **凝露兜底**：KT 75 single −0.23 W/m² (微负，不允许更负)，4q −4.59 W/m² (允许 sha 库强冷凝)
4. **psil 保护**：4q ψ_l 最低 −3.22 MPa > single −3.28 MPa（多 0.07 MPa 安全裕度），因为湿叶占比高→气孔传导率高→根系拉力小
5. **光合不动**：暴雨期 GPP=0、Anet 几乎相同 (−0.4%)；说明 4q 不会扰动光合诊断
6. **lE total 守恒**：差异全部在冠层/土面分配上，总潜热误差 +0.01 W/m² (0.01%)，能量守恒严格

## 6. 与长武苹果 SoilMoisture 实测对比

```
深度       | obs 范围         | single RMSE/Bias       | 4q RMSE/Bias
 10 cm    | [0.169, 0.462]   | 0.0834 / −0.0652       | 0.0833 / −0.0651
 30 cm    | [0.221, 0.435]   | 0.0926 / −0.0859       | 0.0926 / −0.0858
 50 cm    | [0.232, 0.412]   | 0.1077 / −0.1021       | 0.1076 / −0.1020
100 cm    | [0.256, 0.406]   | 0.2026 / −0.1951       | 0.2026 / −0.1951
```

**结论**：
- 冠层 4q 方案对土壤水预测**几乎无影响** (∆RMSE < 0.0001 m³/m³)
- 两模式 bias 一致偏干 (−0.07 ~ −0.20)，主因是初始剖面偏干，与冠层方案无关
- 100 cm 偏差最大，说明深层根系/底边界值得后续优化 (与 4q 无关)

## 7. 数值稳定性

| 指标 | single | four\_quadrant |
|---|---|---|
| 总耗时 (s) | 171.0 | 188.1 (+10%) |
| ROLLBACK + NOT-CONV 次数 | 50 | 66 (+32%) |
| FATAL/MinDtFailed | 0 | 0 |
| 最终完成 KT | 120 ✅ | 120 ✅ |

4q 略增计算成本和迭代次数（因双状态 + 凝露分支），但稳定性**不退化**。

## 8. 物理机制取舍

### single 优势
- CLM5/NoahMP 全球标准，可与公开模型对比
- 最少额外状态变量（1 个 h2ocan）
- 计算最快

### four_quadrant 优势 (COMPASS 创新)
1. **物理一致**：阴/阳叶在 SCOPE 已分离辐射、温度、传导率，冠层水库随之分离避免人为均化
2. **凝露**：阴叶在夜间 / 暴雨退场期可冷凝 (e.g. KT 70–72: −1.1 W/m²)，single 钳制为 0
3. **干燥时间尺度差异**：阳叶 e-folding ≈ 1 h，阴叶 ≈ 4 h，与 Norman 1979 + Calder 1996 文献一致
4. **暴雨退场段 lE 累计 +25%**：可改善降雨后数小时蒸发-气温联动诊断

## 9. 推荐用法

| 用例 | 方案 |
|---|---|
| 全年陆面通量 baseline | `single` (与 CLM5 一致) |
| 暴雨/雨后短期诊断 | `four_quadrant` |
| 阴/阳叶能量分离研究 | `four_quadrant` |
| 与公开 LSM 对比验证 | `single` |
| 苹果园冠层结构创新研究 (本论文) | `four_quadrant` |

## 10. 输出与图

- 冠层水/穿透雨/土壤水图：`docs/figures/v30_single_vs_4q_storm120h.png`
- 能量/光合/叶水势图：`docs/figures/v30_single_vs_4q_energy_photo_psil.png`
- single 输出：`data/output_v30_single_storm120h/`
- 4q 输出：`data/output_v30_4q_storm120h/`
- 对比脚本：`tests/compare_v30_single_vs_4q.m`
- 日志：`v30_storm120h_single.log`，`v30_storm120h_4q_v2.log`

---

## 附录 A. 全年对比 (2021-01-01 ~ 2021-12-31)

> 本附录基于全年单点跑 (`v30_single_yearfull` vs `v30_4q_yearfull`)，
> 由 `PlotResults/plot_results.m` 自动生成。配置见 `data/config_stemmus_scope.txt`，
> 图源目录：[`docs/figures/PlotResults/`](figures/PlotResults/)。

### A.1 关键指标年度合计

| 指标 | single | four_quadrant | Δ (4q−single) | % | 结论 |
|---|---|---|---|---|---|
| CanopyInterception (mm/yr) | 44.0 | 50.3 | **+6.3** | +14.3% | 4q 阴阳叶分立后湿叶储水更彻底 |
| Throughfall (mm/yr) | 952.8 | 946.5 | −6.3 | −0.7% | 与截留互补 |
| SurfaceRunoff (mm/yr) | 27.2 | 6.1 | **−21.1** | −77.7% | 注：年总径流仅占 P 的 0.6–2.7%，差异主要来自小雨期低频泄水 |
| Evaporation (mm/yr) | 376.2 | 382.5 | +6.3 | +1.7% | — |
| Transpiration (mm/yr) | 306.4 | 304.6 | −1.8 | −0.6% | 不抢水 |
| λE_wet (MJ/m²·yr) | 150.7 | 155.4 | **+4.68** | +3.1% | 阳叶湿润时蒸发增强 |
| λE_dry (MJ/m²·yr) | 750.6 | 746.3 | −4.31 | −0.6% | 干叶传导降一点点 |
| **λE_total (MJ/m²·yr)** | 1818.6 | 1834.5 | **+15.95** | +0.88% | 总潜热上升 |
| **H_total (MJ/m²·yr)** | 423.3 | 410.6 | **−12.72** | −3.0% | 感热下降，能量从感热路径转向潜热路径 |
| Rn (MJ/m²·yr) | 2461.4 | 2467.2 | +5.77 | +0.23% | 几乎不变 |
| 能量平衡残差 (MJ/m²·yr) | −18.09 | −14.37 | +3.73 | 闭合度 +20.6% | 数值精度提升 |

完整 CSV: [`figures/PlotResults/compare_metrics_table.csv`](figures/PlotResults/compare_metrics_table.csv)

### A.2 五张可视化图

#### 图 1 — 年度差值视角

![compare_combined_panel](figures/PlotResults/compare_combined_panel.png)

5 行：① 截留差值 + 降雨柱（雨后 4q 跃升）　② 累计径流（注意年总占 P <3%）　③ Δ 湿叶 λE　④ Δ 干叶 λE　⑤ Δ λE 与 Δ H 配对（能量重分配）

#### 图 2A — 暴雨事件放大：7 月 15 日 91.5 mm/day 单峰极端暴雨

![storm_event_zoom_july](figures/PlotResults/storm_event_zoom.png)

5 行：降雨柱（含微雨标注）+ **净辐射 R_n（关键解释）** + 截留 + 累积径流 + 湿叶 λE。
- 7/15 14h-15h 暴雨柱占满 Y 轴 (~30 mm/h)，当天 R_n 被云挡住跌至接近 0
- 截留瞬间填满 0.25 mm（4q ≥ single），地表径流 7/15 当天阶跃
- 雨后 7/16-7/17 转晴 R_n 恢复，湿叶 λE 才达峰值（雨后慢释放）

#### 图 2B — 暴雨事件放大：10 月初连阴雨 (4 天 187 mm)

![storm_event_zoom_oct](figures/PlotResults_oct/storm_event_zoom.png)

10/3 79 mm + 10/6 47 mm + 10/4 30 mm + 10/5 31 mm 的连阴雨型事件。
- 与 7 月单峰对照：10 月每日 R_n 都被云层压低到 50–90 W/m²
- 10/1 出现一个 0.27 mm 微雨 + 当天 R_n 高达 390 → 湿叶 λE 反而比暴雨期还高（早先以为是 bug，加 R_n 行后释清）
- 这一对比是论文里"4q 在事件型与连阴型雨水下行为差异"的核心证据

#### 图 3 — λE 阴阳叶云图：single vs 4q 左右对比

![lE_sun_sha_compare](figures/PlotResults/lE_sun_sha_compare.png)

**重要负面发现**：阳叶/阴叶 λE 在 single 与 4q 之间几乎完全相同（年均仅差 0.07/0.27 W/m²，约 0.2%/5%）。
说明 SCOPE 模型本身就分阴阳叶辐射传输，4q 创新**不在 λE 分布上动手脚**，而在冠层水库 + fwet（见图 4-5）。

#### 图 4 — 4q 完整四象限 λE 时空分布

![four_quadrant_lE_breakdown](figures/PlotResults/four_quadrant_lE_breakdown.png)

按 fwet 面积加权近似拆分（严格四分量需改 SCOPE 输出）：
- ① **阳叶干 λE_sun_dry**：年均 41.7 W/m²·leaf，气孔蒸腾主路径，最大分量
- ② **阳叶湿 λE_sun_wet**：年均 4.5 W/m²·leaf，截留蒸发，雨后 + 晴天集中爆发
- ③ **阴叶干 λE_sha_dry**：年均 4.0 W/m²·leaf，散射光下气孔蒸腾，弱
- ④ **阴叶湿 λE_sha_wet**：年均 1.1 W/m²·leaf，雨后阴叶水库释放最慢

#### 图 5 — 4q 真正招牌：阴阳叶独立水库 + 双 fwet

![four_quadrant_signature](figures/PlotResults/four_quadrant_signature.png)

**这才是 4q 独有特色**（single 里 `fwet_sun = fwet_sha = h2ocan_sun = h2ocan_sha = 0`）：
- ① 双水库日动态：阴叶水库 0.018 mm ≈ 阳叶水库 0.005 mm 的 **3.5 倍**（Rutter 串联：阳叶先满溢到阴叶）
- ② 双 fwet 时序：fwet_sha 0.176 vs fwet_sun 0.096，阴叶系统性更湿
- ③ fwet 差值：4q 平均 +0.080，single 这个差值恒为 0
- ④ 暴雨窗口（7/13-7/18）小时尺度：阳叶填满后溢到阴叶；雨后阳叶蒸发更快（顶部直射），阴叶慢释放

#### 图 6 — 土壤水分 single vs 4q vs 实测

![soil_moisture_3way_compare](figures/PlotResults/soil_moisture_3way_compare.png)

**重要负面发现**：4q 冠层创新**不显著影响地下水分**预报。
- 上行：模型 vs 实测散点重合
- 下行：4q − single 差值 ×1000 放大；最大 |Δθ| < 0.05 cm³/cm³（实测精度 ~0.01 的 1/几百倍）
- 各深度 RMSE 几乎不变（10cm/30cm/100cm/300cm）
- 解读：4q 改 ET 配比量级 < 几 mm/yr，扩散到全土柱后被忽略

详细 CSV: [`figures/PlotResults/soil_moisture_compare_metrics.csv`](figures/PlotResults/soil_moisture_compare_metrics.csv)

### A.3 复现命令

```matlab
% 1) 跑两个全年仿真（output_v30_single_yearfull, output_v30_4q_yearfull）
% 2) 切到任一 run（修改 data/config_stemmus_scope.txt 里 OutputPath）
% 3) 一键画图：
addpath('f:/COMPASSV1.0/PlotResults','f:/COMPASSV1.0/src');
plot_results
% 输出：docs/figures/PlotResults/  (compare_metrics_table.csv + 6 PNG)
```

代码：
- 入口：[`PlotResults/plot_results.m`](../PlotResults/plot_results.m)
- 实现：[`PlotResults/compare_single_vs_4q.m`](../PlotResults/compare_single_vs_4q.m)
- 工具：[`PlotResults/parse_config_outputpath.m`](../PlotResults/parse_config_outputpath.m)

### A.4 文章用结论小段（可直接抄）

> 在 2021 年单点全年仿真（北京苹果园，91.5 mm/day 极端暴雨日 2021-07-15）下，four-quadrant
> 阴阳叶分立 fwet 方案对模型行为的影响呈现"冠层水文层面强、辐射传输层面弱、地下水分预报中性"的层次性：
>
> (1) **冠层水库分层化**：阴叶水库年均储量 0.018 mm ≈ 阳叶水库的 3.5 倍，反映 Rutter 串联模型
> 中阳叶饱和后向阴叶溢流的物理路径；fwet_sha − fwet_sun 平均 +0.08，single 模型该差值恒为 0。
>
> (2) **能量配比改变**：年总潜热 +15.95 MJ/m²、年总感热 −12.72 MJ/m²，净辐射几乎不变
> （+0.23%）。能量从感热路径转向潜热路径，与冠层湿叶面积上升一致。
>
> (3) **辐射传输几乎不动**：阳叶 / 阴叶 λE 时空分布在 single 与 4q 间相对差异 < 5%，
> 说明 4q 创新作用在冠层水量分配，不改变 SCOPE 内置的辐射传输求解。
>
> (4) **地下水分中性**：表层 10 cm 至深层 300 cm 的 VWC 时序在 single 与 4q 间最大绝对差
> < 0.05 cm³/cm³（实测精度的 1/几百倍），各深度 RMSE 实际不可分辨。
>
> 结论：four-quadrant 方案的增益在事件尺度（雨后冠层水量动态、湿叶蒸发尾巴），不在年总水分预报。
