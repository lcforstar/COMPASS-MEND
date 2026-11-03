# 输入表原始注释沿革（归档）

Editor: Wenjie Liu; Date: 2026-10-31.

以下说明从输入表开头移出；当前输入表只保留一行说明。原始表头、数值及历史署名均保留。

## soil.csv

```text
soil.csv: 基本属性 | 水力参数 | 无机氮参数 | DOC/DON参数 | BIOPARTICLE介质参数 | MEND分层反应参数
# 分层土壤、介质及11个当前MEND参数列位于本表；其他全站设置见site.toml。
# DOC/DON平衡吸附、弥散及微生物输运参数可随层变化；速率单位见列名。
# Change: 2026-10-11 Wenjie Liu -- four uniform MEND kinetic traits moved to site.toml without changing values.
# 深度(cm),黏粒含量(%),砂粒含量(%),土壤有机碳(mg/kg),容重(g/cm3),pH(-),饱和导水率(cm/d),饱和含水量(m3/m3),残余含水量(m3/m3),VG-n(-),VG-alpha(1/cm),铵态氮扩散(cm2/d),硝态氮扩散(cm2/d),无机氮弥散(cm),铵吸附系数(cm3/g),反硝化Vmax(ugN/cm3/d),反硝化Km(ugN/cm3),尿素水解(1/d),硝化Vmax(ugN/cm3/d),硝化Km(ugN/cm3),经验矿化(1/d),铵肥溶解(1/d),DOC平衡吸附Kd(cm3/g),DON平衡吸附Kd(cm3/g),DOM弥散(cm),微生物附着率(1/d),微生物脱附率(1/d),微生物截留率(1/d),微生物运移弥散(cm)
# Change: 2026-10-12 Wenjie Liu -- all 8+2 MEND candidates moved from site.toml to soil.csv; all 18 rows retain the existing local values; no calibration or depth differences introduced.
# MEND新增列单位: Vd[h^-1], KsPOM1[mg C/g], Vg[h^-1], KsDOM[mg C/g], ratio[-], alpha[-], mortality_multiplier[-], CUE[-], SWP_A2D[MPa], tau[-].
# Change: 2026-10-17 Wenjie Liu -- user-specified KsDOM prior 0.02 mg C/g; other values unchanged.
# Editor: Wenjie Liu; Date: 2026-10-18.
# Change: 2026-10-18 Wenjie Liu -- user-specified pH-disabled Vd prior 0.10 h^-1 in all 18 layers; all other cells unchanged.
# Change: 2026-10-18 Wenjie Liu -- freeze sample-39 Ks/alpha/top VG-n; set user nitrification prior at 20-120 cm; Vd stays 0.10 h^-1.
# Change: 2026-10-18 Wenjie Liu -- user smooth NH4 Kd and nitrification profiles at all 18 parameter depths; every other column unchanged.
# Change: 2026-10-18 Wenjie Liu -- restore only 20-120 cm nitrification to frozen csv_cold_regression values times sample-39 factor 0.580343248604; other cells unchanged.
# Change: 2026-10-20 Wenjie Liu -- freeze selected layered_04 Sample 39 DIN kinetics, NH4 Kd and split DIN/DOM dispersion; retain hydraulics, MEND and empirical mineralization cells.
# Change: 2026-10-21 Wenjie Liu -- Reduced-DOM diagnostic: CUE .42, Vd .25 h^-1, KsDOM 2.0 mg C/g, paired Kd .10 cm3/g; DIN/hydraulics unchanged.
# Editor: Wenjie Liu; Date: 2026-10-21.
# Change: 2026-10-21 Wenjie Liu -- user absolute targets Yg20=.45, KsDOM=.06 mgC/g, shared Vd=.15 h^-1; raw40-300cm DOC Kd=.50; raw<=40cm Km_nit=35; all other cells retained.
# Editor: Wenjie Liu; Date: 2026-10-21.
# Change: 2026-10-21 Wenjie Liu -- annual test: Yg20=.42, KsDOM=2 mgC/g, Vd=.15 h^-1, DOC Kd=.10 cm3/g throughout; raw<=40cm Km_nit=20. Other cells retained.
# Change: 2026-10-24 Wenjie Liu -- isolated S1; user-authorized source/sink and transport candidate; all other columns retained.
# Editor: Wenjie Liu; Date: 2026-10-25.
# Change: 2026-10-25 Wenjie Liu -- promote S1 after PASS at KT=8760; all 18 layers mend_vd_pom_h=0.25 h^-1 and mend_ks_dom_mg_g=10.0 mg C/g; all other parameter cells retained.
# Editor: Wenjie Liu; Date: 2026-10-25.
# Change: 2026-10-25 Wenjie Liu -- N5 isolated candidate; only raw20/40cm Vmax_nit/Km_nit/kg cells may differ; sample_linear and all deeper raw samples/MEND unchanged.
# Change: 2026-10-26 Wenjie Liu -- server Sample 0: promote N5 CN10 and carbon baseline; user-specified nitrification Vmax at raw 60-300 cm. Other cells retained.
# Change: 2026-10-31 Wenjie Liu -- sole-owner Qmax[mg C/g dry soil], Kba[g dry soil/mg C], Kdes[mg C/g dry soil/h] moved from TOML; all 18 rows preserve 1.5, 6.0, 0.006.
```

## organic_initial_pools.csv

```text
# Sample 0: measured 2000 SOC with user-specified cold-start priors.
# Existing MOC/QOC/EM names correspond to MOM/QOM/ENZM. AOM starts at zero.
# Editor: Wenjie Liu; Date: 2026-10-17.
# Change: 2026-10-17 Wenjie Liu -- replace historical warm-state fractions; preserve column contract.
# Literature attribution supplied by user; SOURCE UNVERIFIED 2026-10-17.
# Change: 2026-10-18 Wenjie Liu -- depth-dependent low DOC/QOC priors; transfer all removed C to MOC; MON/QON/DON C:N=10; exact decimal 12-pool closure.
# Change: 2026-10-21 Wenjie Liu -- 1.2% SOC total DOM diagnostic (.009 DOC + .003 QOC), funded from MOC; total SOC/N preserved at MON/QON/DON C:N=10. Actual phases follow equilibrium.
# Editor: Wenjie Liu; Date: 2026-10-21.
# Change: 2026-10-21 Wenjie Liu -- user DOM priors DOC=QOC=.005 at centers<=40cm and .002 at40-300cm; MOC is exact carbon remainder; paired MON/QON/DON retain CN10; >300cm unchanged.
# Editor: Wenjie Liu; Date: 2026-10-21.
# Change: 2026-10-21 Wenjie Liu -- 0-20cm initial litter C=.0003 g/cm3, paired N using Wei et al.2015 Table1 straw C:N=669.15/9.83; components .25/.50/.25. Existing 12 SOC fractions unchanged; extra unreacted AOM belongs to initial full C/N inventory.
# Change: 2026-10-23 Wenjie Liu -- 按用户确认清零初始麦茬全部组分碳氮；保留十二碳库比例及施肥输入。
# Change: 2026-10-24 Wenjie Liu -- user-authorized 2006-09-15 measured DOC profile initializes D at the 2000-06-15 start; temporal proxy, not same-date observation.
# Source: NM冬小麦实测.xlsx, DOC和DON!B3:B17, mg C/kg dry soil; 0-300 cm layer-bottom convention. See doc_initial_2006_proxy.json.
# Mapping: overlap-weighted measured horizons onto the existing 58 layer centers; >300 cm retains prior D. Q is independent and unchanged.
# Carbon: transfer D difference from MOC, preserving each layer SOC. All N fractions are unchanged; D/DON and MOC/MON C:N may therefore change.
# Change: 2026-10-24 Wenjie Liu -- isolated S1; Q0=Q*(D0) at 20 C, Q increment funded only from same-layer MOC; every N input unchanged.
# Editor: Wenjie Liu; Date: 2026-10-25.
# Change: 2026-10-25 Wenjie Liu -- promote S1 after PASS at KT=8760; 20 C Langmuir Q0 equilibrium (Kba=6.0, Qmax=1.5); extra QOC funded by same-layer MOC; free DOC, total C and all N fractions retained.
# Editor: Wenjie Liu; Date: 2026-10-25.
# Change: 2026-10-25 Wenjie Liu -- N5-N8 common CN10 ONLY in native elements 0-300 cm; exact N1 DON/QON/MON fractions in that interval; below300 cm S1 unchanged; carbon cells unchanged.
# Change: 2026-10-27 Wenjie Liu -- authorized isolated ablation B_partition; exact changed cells recorded in case audits; no model-code change.
```
