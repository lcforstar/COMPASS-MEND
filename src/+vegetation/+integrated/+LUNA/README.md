# +vegetation/+integrated/+LUNA — CLM5 LUNA v1.0 移植

## 物理问题

叶片氮 (N_leaf) 如何动态决定光合容量 (Vcmax25 / Jmax25)。

闭环最后一环：N_leaf 涨 → Vcmax25 涨 → GPP 涨 → NPP/D_N 涨 → root uptake 涨 → N_leaf 涨。

## 算法 = CLM5 LUNA (Ali 2016 GMD)

5 池机理分配 N_leaf:
- N_lc: 光捕获 (chl + thylakoid)
- N_et: 电子传递 → Jmax
- N_cb: 羧化 (Rubisco) → Vcmax
- N_resp: 维持呼吸酶 → Rd
- N_store: 储存 (调控缓冲)

每日爬山优化: max(PSN − C_cost)
每日 ≤ 10% 变化率 (Q10Enz=2, Suzuki 2001)

## 模块组织

```
+LUNA/
  ├─ lunaParams.m              全 LUNA 常数集中 (1:1 from LunaMod.F90:55-78)
  ├─ computeLuna.m             主入口 (S4 待写)
  ├─ nitrogenAllocation.m      爬山优化核心 (S3 待写)
  ├─ nitrogenInvestments.m     N→Vcmax/Jmax 投资 (S2 待写)
  ├─ photosynthesisLuna.m      内嵌 Farquhar+BB (S2 待写)
  ├─ +nue/                     NUE 计算
  │   ├─ nueRef.m              25°C, 380ppm 参考 NUE
  │   └─ nueCurrent.m          当前条件 NUE
  ├─ +tempresp/                温度响应
  │   ├─ vcmxTKattge.m         Kattge & Knorr 2007
  │   ├─ jmxTKattge.m          
  │   └─ respTBernacchi.m      Bernacchi 2001
  └─ +accum/                   10 天滑动累加 (S4 待写)
```

## 输入 / 输出

**Inputs**:
- `NitrogenVariables.N_leaf(KT-1)` [g N/m² ground] — 上 KT 末叶 N pool
- 10 天滑动均气候 (PAR_10d, T_leaf_10d, RH_10d, CO2_10d, r_b_10d, dayl)
- canopy 状态 (LAI, SLA, nrad layer count)

**Outputs**:
- `NitrogenVariables.Vcmax25_LUNA` [μmol/m²/s] — 每日更新 Vcmax25
- `NitrogenVariables.Jmax25_LUNA` [μmol/m²/s] — 每日更新 Jmax25
- 5 子池诊断: PNlc / PNet / PNcb / PNresp / PNstore

## 跟外部耦合

```
+integrated/+LUNA → biochem.m (line 125):
    Vcmax25_base = NitrogenVariables.Vcmax25_LUNA
    Jmax25_base  = NitrogenVariables.Jmax25_LUNA  (不再用 R_JV)

biochem.m 后续:
    Vcmax = f_T(Vcmax25) × sfactor   ← sfactor 来自 Plant_hydraulics.m (不动)
    Jmax  = f_T(Jmax25)  × sfactor
    Rd    = f_T(Rd25)               (不缩放)
```

三层串联: **N 控容量 (LUNA 每日) → 温度响应 (Kattge 每 KT) → 水力胁迫 (sfactor 每 KT)**.

## 参考文献

- Ali AA et al. (2016) GMD 9:587-606 — LUNA v1.0 model description
- Xu C et al. (2012) New Phytol. 193:30 — N allocation algorithm origin
- Lawrence DM et al. (2019) CLM5 Tech Note Ch.9 — formula numbering
- Rogers (2014) Photosynthesis Research — Fc25 = 294.2
- Evans JR (1989) Oecologia 78:9-19 — Cb (chl N for light capture)
- Bernacchi CJ et al. (2001) PCE 24:253 — Kc25/Ko25/Cp25 + RespT
- Kattge J & Knorr W (2007) PCE 30:1176 — Vcmax/Jmax T-acclimation
- Suzuki T et al. (2001) — Rubisco ~7d half-life

源代码 (本地): `D:\OneDrive\CTSM-release-clm5.0.37\src\biogeophys\LunaMod.F90` (1352 行)

## 实施进度

- [x] S1: lunaParams + tempresp + nue (6 文件)
- [x] S2: nitrogenInvestments + photosynthesisLuna
- [x] S3: nitrogenAllocation (爬山优化)
- [x] S4: computeLuna 主入口 + 10 天 buffer
- [x] S5: biochem.m hook + COMPASS.m 调度
- [x] S6: 诊断 (Vcmax25/Jmax25/FNCa/PNlc/PNcb 逐小时输出)
- [x] S7 (2026-06-02): N_leaf → Cab 叶绿素转换 (nLeafToCab.m, 两路并列输出, 见下)

## N_leaf → Cab (叶绿素含量, S7)

`nLeafToCab.m` 把 LUNA 的叶 N 转成叶绿素含量 Cab [μg/cm²] (PROSPECT/Fluspect 输入)，
两条独立文献路径**并列输出**做值对比 (非开关切换):

| 路径 | 输入 | 公式 | 出处 | 特点 |
|---|---|---|---|---|
| `cab_from_nlc` (唯一) | N_lc = PNlc×FNCa [g N/m²叶] (LUNA 动态光捕获 N) | Chl = N_lc/M_N × (1/50) × M_Chl | Evans 1989 Oecologia 78:9 p.12 (50 mol thylakoid N/mol Chl) | **带 LUNA 氮分配信号**, 物理自洽 |

> ~~Route B (`cab_from_totaln`, 总叶N×3.8e-3)~~ **已删 (2026-06-02)**: lnc=1/(SLA×CN_L) 在 lnc_opt=false 下是固定常数, Route B 算出近常数 Cab (v21 全年恒 98.4≈固定100), 无动态无意义。保留唯一动态路径 cab_from_nlc。

常数: M_N=14.007 g/mol, M_Chl=893.5 g/mol (Chl a, C55H72MgN4O5)。Cab[μg/cm²] = Chl[g/m²叶] × 100。

**单元测试 (2026-06-02)**: 典型苹果叶片 (lnc 0.5-3, PNlc 0.2-0.4, FNCa 1-2) 两路 Cab 落在
20-80 μg/cm² 合理区间; 路径 A 随 LUNA 分配动态变, 路径 B 仅随 lnc 变。高 PNlc×FNCa 时路径 A
可冲到 ~100, 全年跑后需看实际范围。**端到端全年验证待下次全年跑**(产出两列对比)。

**尚未做**: Cab 耦合进 SCOPE 辐射 (fluspect leafbio.Cab)。届时再加 cfg 开关选用哪一路。
