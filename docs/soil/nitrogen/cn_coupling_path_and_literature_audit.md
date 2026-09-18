# COMPASS 植物 C / N / 光合 / 分配 — 方法学 + 路线 + 待办

> 📍 **本文档已重写 (2026-05-29)**。正文 (§A-§F) 是论文级方法学整理 (闭环图 + 方程 + 路线 + 参数 + 待办 + 文献)。
> 历史实施档案 (M1~M2-B.8 怎么一步步做 + 踩的坑) 压缩进**附录 H** (原 §0-§12.33 原样保留, 仅查阅)。
> 备份: `docs/.cn_coupling_backup_pre_rewrite.md`。当前版本: **v10** (Changwu LMWE 苹果园站, PFT 23)。

---

## 📍 路线图 (全文唯一权威 — 严格按本表自上而下、块内自上而下执行)

> 全文只此一份路线图 (阶段 + 细项 + 优先级已三合一); §C / 附录 §0 / §12.33.5 / 旧 §E / 旧 §12.36 均已删, 防误读。
> **执行规则**: 块顺序 = 优先级。**P0 全做完才进 P1, P1 完才 P2**。块内按行从上到下。✅=完成, ⏳=待做, ✗=弃。
>
> **当前位置 (2026-06-02 更新)**: P0 主体推进完: F1✅/F5✅/13c✅(提交65c75f3)/f_alloc标定✅(提交74e3d21,产量18306→46868对标49315差5%)。
> **现役焦点**: CUE 偏低 0.43(实测0.71)根因链已锁定→SLA0 0.0125→0.028 修 lnc(v22验证中 output_v22_sla028)。F3(N_labile虚高)[**2026-06-03 真凶定位**: 全年验证 output_v24_F3combined 证明 N_labile 累积**不是 demand 公式/retransl 问题, 是根吸超 demand ~3.7×不听需求**(年根吸52.5 vs demand~14 gN/m²); fix#1(demand 看合并储备 labile+storage, 合 LPJ nstore guess.cpp:1679)已做、守恒0.0000%但**单独无效**(reserve仍12→58); **已修✅ 2026-06-03**: 真凶=solveRootSurface.m:104 收敛单位 bug(通量差[kg/m/s]比 C_root 浓度容差[kg/m³]→第1次迭代假收敛→根吸 10-50×demand), 通量判据改相对容差 1 行修复; 全年 output_v25_F3fix: **N_labile 末值 47→5.1, 年根吸 52.5→10.6 gN/m², supply factor mean 17→0.87, 守恒 0.0000%, 不崩**; 非建模选择、是 bug; fix#1(合并储备)+fix#2(收敛)两修都留; 详见 P0 记忆]/F2(LAI在线)仍挂。
> **遥感正演线 (并行)**: Cab 动态✅接通(cab_from_luna) → 含水量(W_leaf/VOD/Cw 叶照茎 M=C×WDR×rwc)设计定调待实施(§12.36.3) → 9 参数自动化(排含水量后)。
> **详细进展散见 §12.30(13c/呼吸)/§12.32(F2)/§12.36(fluspect 9参), 下方路线图已据此重排。**

### ✅ 已完成 (阶段 1-12 + 12.5·✓ + F1, 基线)
| 阶段 | 内容 | 出处 |
|---|---|---|
| 1-6 | 土壤N(DAISY+Nacry根吸) · 根吸demand-driven · 植物N5器官池+LPJ NPP+NEE · LUNA N→Vcmax · 植物C6池+NSC · Turnover/retransl/采收 | M1~M2-B.8 |
| 7-12 | Phase D物候(Sykes+CLM5+DAISY DS)+LAI双引擎 · LUNA strict(Jmaxb1=0.244) · R_m_stem只活木 · Rd光抑制 · C:N对齐FruitTree+gs_g1=8.2 | v3~v10 |
| 12.5·✓ | FlexibleCN(C:N带内浮动) · C主导init · C_labile初值统一C主导 · 植物N守恒残差(C/N对称) · harvest改DS_fruit(解耦frPHU) · 方案A fcur越冬储备 · §B文档补全(B.0池/B.9物候/B.10 NEE) · lLUNA移veg | §12.34序/§B |
| **12.5·F1 ✅** | **R_g 0.25→0.11 全24 PFT** (CLM5 grperc; 出处核实 Sitch PDF p8 vs CLM5 Tech Note)。Jan-Aug smoke 验 (同截 KT5831): **GPP 不变 1178.9** (变量隔离✓), **R_a −12.6%** (854→747), **C_fruit +3.4%** (111→115)。守恒 C 5.4e-12/N 1.4e-13 机器精度。**是产量低主因之一非全部** (省碳按 f_alloc 仅 30% 进果; 抬产量主力在 F2 GPP) | §12.34.1 |
| **12.5·F5 ✅** | **细根周转物候门控** fineRootTurnover 加 phen_gate=ps.phen (休眠期降速, CLM5 bglfr CNPhenologyMod.F90:2623)。Jan-Aug smoke: cum_N_root_litter 106.9→83.1 (−22%, 冬季少掉根), N_root_fine 末 7.05→7.74 (池保住), 守恒 N1.2e-13/C6.0e-12 | §12.34.1 |
| **树体几何 ✅** | DBH 4.2→**21.3** cm / canopy_height 2.0→**3.1** m / tree_CAH 4.13→**10.2** m² (旧值假, 改长武实测 Wang et al. 2024 Agronomy 14(5):963 Tbl1, 18龄乔砧大树 A15-A20插值)。佐证 C_stem_wood=1121 量级合理(非虚高) | config_vegetation |
| **13c 茎双池 ✅** | 提交 **65c75f3**。单 N/C_stem_wood 拆 livestem(CN_lw60,呼吸)/deadstem(CN_dw500,不呼吸), flivewd分配+lwtop周转, 解茎周转(N3)/边材呼吸/茎N需求。Jan-Aug smoke 守恒机器精度, 周转生效。CN_lw/flivewd 后按 Dabbaghi 实测改(60/0.15) | §12.30.3a |
| **f_alloc 标定 ✅** | 提交 **74e3d21**。叶0.40→0.12/果0.15→0.50/木0.20→0.30/根0.25→0.08 (Zanotelli 2013+Wu 2012 实测; 旧值叶果反了)。全年验: **产量 18306→46868 kg/ha** (对标实测 49315 差5%, 原低3倍根治), C_leaf 403→185, 守恒机精 | §12.32.9 |
| **C:N 出处补全 ✅** | 6 器官 C:N 全标 Dabbaghi 2022 Table C1 出处 (leafCN19.7/frootCN32/livewdCN60/grainCN139, Zanotelli 实测); CN_L=19.7 是实测正确值(非虚高) | getPFTParameterSet |
| **Cab 动态 ✅** | cab_from_luna: LUNA 光捕获 N → 叶绿素 Cab (Evans 1989 50:1), 替 SCOPE 固定100。Route B(总N)删(伪动态)。fluspect 9 参第1个过程化 | §12.36 |
| ✗弃 | cushion储备回补(v8/v9抽水泵bug, v10删) · Jmaxb1率定(0.244→0.18 LUNA补偿无效, 回退) · Cab Route B(lnc固定→伪动态) | — |

### 🔴 P0 — 现在做: 苹果 Bar B 主线 (阻塞"本周定模型物理"; 严格按序)
| 序 | 项 | 做什么 | 风险 | 出处 |
|---|---|---|---|---|
| **1 ✅** | ~~F1 生长呼吸~~ | R_g 0.25→0.11 (R_a−12.6%/C_fruit+3.4%) | 低 | §12.34.1 |
| **3 ✅** | ~~F5 细根门控~~ | phen_gate (cum_root_litter−22%) | 低 | §12.34.1 |
| **13c ✅** | ~~茎双池~~ | 提交 65c75f3 (见已完成块) | 中 | §12.30.3a |
| **f_alloc ✅** | ~~分配标定~~ | 提交 74e3d21 (产量对标实测差5%, 见已完成块) | 中 | §12.32.9 |
| **CUE ⏳现役** | **SLA0→修 lnc→降 Vcmax/Rd** | CUE 0.43(实测0.71)偏低: 根因链=喂 LUNA 的 lnc=1/(SLA×CN_L) 固定公式, SLA0=0.0125→lnc4.06(实测2-2.3的2倍)→Vcmax/Rd高→R_a超2.3倍。**改 SLA0 0.0125→0.028** (Dabbaghi/Zanotelli实测, lnc→1.81)。v22验证中(output_v22_sla028): 看 CUE 是否升向0.71+产量/LAI不崩。Jmaxb1率定失败已弃(LUNA补偿) | 中 | §12.30遗留2/§12.36.2 |
| **F3 ⏸️挂** | **N_labile 虚高** | 13c后仍 24.9(sap_term 2.17被regrow 11.6压死, max取regrow; 根因=N进labile无出口非公式)。代码保留不回退(守恒✓)。后续: F3需求加"超target停D_N根吸"或学CLM5过剩外排 | 中 | §12.30.3a遗留1 |
| **F2 ⏸️挂高风险** | **C_leaf→LAI 在线** | 范式转变非切开关: EPIC离线预算查表→leafc_sla在线算(叶碳每步变)。主循环现绕过computeLAIfromPhen(hotfix)。⚠️CLM5警告leafC↔LAI不收敛需滞后+peaklai。**离线验证已做: f_alloc标定后 leafc_sla LAI峰2.62 vs EPIC2.79(0.94)可切不崩**, 但生长季低15%+晚30d。架构改造压轴, 不擅自切 | **高** | §12.32.9/§14 |
| └诊断 | (随手) D1 ebal_residual · D3 GPP csv全0 · 夜间叶呼吸漏算(R_m_leaf夜=0, 补了CUE更低暂不动) | 低 | — |

### 🛰️ P0b — 遥感正演线 (并行, 碳氮水→叶光学/微波, GMD 特色, 反"缝合怪")
| 序 | 项 | 做什么 | 状态 | 出处 |
|---|---|---|---|---|
| RS1 ✅ | **Cab 动态** | LUNA 光捕获N→叶绿素(nLeafToCab, cab_from_luna), 替固定100。Route B删 | 完成 | §12.36 |
| RS2 ⏳ | **含水量统一 M=C×WDR×rwc** | 叶照茎(13c茎已此式): 碳驱动饱和含水×capacitance相对饱和度rwc(乘非加)。**RWC_leaf已算(wrapper:109)只需接出**: 注册rwc_leaf_v2+load_state读+VOD叶用+fluspect Cw。修当前叶茎不对称(茎有rwc叶无)。⚠️WDR_leaf 与 MC 须统一(现2.5 vs fluspect隐含0.43打架) | 设计定调待实施 | §12.36.3 |
| RS3 ⏳ | **fluspect 9参自动化** (排 RS2 之后) | 9参(Cab/Cca/Cdm/Cw/Cs/Cant/N/V2Z/fqe)能过程化的自动从碳氮水状态算: Cab(✅)+Cw/Cdm(含水量改好后, RS2带出)+Cca(0.25×Cab半动态)。N叶肉结构/Cant暂固定(无过程)。**依赖 RS2 含水量先理顺** | 待 RS2 | §12.36 |
| RS-时机 | 实施: v22(SLA)验完 + 森林告一段落再动手 (避免同时太多变动) | — | — |

### 🟡 P1 — 苹果验证后: 接土壤 + 物候动态 (P0 全通才做)
| 序 | 项 | 做什么 | 依赖 | 出处 |
|---|---|---|---|---|
| **5** | **12.5·F4 litter C:N** | litter C:N=60 喂 SOM (FruitTree; 非活叶19.7) + 解耦 C/N 回收 (turnoverCarbonMirror.m:78 临时件转正) | 接师弟 SOM | §12.34.1 |
| **6** | **阶段13 动态分配** | f_alloc 静态→动态: 13a 苹果(FruitTree Eq4-7, 含果实多阶段TSum+展叶期优先=P4) / 13b LPJ grass(leaf:root) / 13c LPJ tree(pipe+sap/heart双池=N3/GL2, DBH; 含 N1/N2 三通道充值; **边材双池收口包见 §12.30.3a**: 一并解决茎周转/F3储备/边材容量/边材呼吸) / 13d crop(DS) | 12.5全通 | §F |
| 6·arch | └ **统一 growth_form 标签**: 物候 4-flag (§B.9) 与分配 3-套分流共享同一 PFT 生长型分类, 建分配前先对齐, 别造第二套打架标签 | — | §B.9 |
| 7 | 阶段15 lnc_opt=.true. (土壤 N 真反馈光合, lnc=N_leaf/LAI; =P2) | F2 后 | — |
| └✅已解 | ~~P3 slatop 0.0125 vs 0.028~~: 已改 0.028 (Dabbaghi/Zanotelli 实测, 修 lnc/CUE, 见 P0 CUE 行); v22 验证中 | — | §12.36.2 |

### 🔵 P2 — 通用 LSM 多 PFT 结构 (阶段16, **按真实验证点驱动**, 非苹果阻塞; 不为"通用"提前堆砌无法验证的结构)
| 触发站点 | 项 | 出处 |
|---|---|---|
| 森林/乔木站 | GL1 粗/细根分离 · GL2 live/dead 木质+CWD(=N3 sap→heart) · PH3 常绿 bglfr 连续落叶 | §12.35 |
| 草地/作物站 | PH2 stress_decid 实现(4累积器) · PH1 crop 物候(出苗/DS多阶段/春化/光周期) | §12.35 |
| 落叶通用 | GL3/PH5 storage/xfer 转运三联体 · GL4 gresp_storage 生长呼吸储备 | §12.35 |
| 热带 | PH4 phenology ANY 混合型 (dispatcher 现强制 4 选 1) | §12.35 |
| 多年 spin-up | D4 N_storage 跨年闭环 · 植被动态(establishment/mortality, 远期) | — |

### ⚪ 远期 / 他人 / 中低优先 (不排队, 按需穿插)
- 阶段17 SOM 3池(修 NEE 反符号, 师弟) · 阶段18 多 patch 陆面扩展。
- D2 水残差确认+bottom_drainage · 文献 Ghimire 2016/Poorter 2012 下载 · 修剪模块 prune_fr。
- 中/低 (§12.34.2): 干旱 β 下调 R_m · 茎分配随 NPP · 回流 N 夏季可用 · 奢侈吸收存 reserve N · NSC cap 记债 · 根 R_m 分层温度。

---

## §A 模型概述 + 三大闭环

COMPASS = 土壤 (STEMMUS 水热 + STEMMUS_N) + 植物 (光合 LUNA-Farquhar + C/N 6 池 + 物候 + 水力) 全耦合单点陆面模型。
区别于 CLM5/LPJ-GUESS: **同时**有 DAISY+Nacry 单根 ADE 根吸 + CLM5 LUNA + DAISY 6 池 C/N 循环 + 植物水力 prognostic, 站点严格守恒 (机器精度), 多 PFT (24)。

### §A.1 三大闭环逻辑图 (水 / 碳 / 氮, v12 态: C 主导 + C:N 弹性)

```
   ┌──────────────────────── 气象 forcing (T, PAR Q, CO2 Ca, VPD, 降水) ────────────────────────┐
   │                                                                                            │
   ▼                                                                                            │
┌─────────┐  lnc=1/(slatop·leafcn)  ┌──────────────┐                                            │
│LUNA 每日 │ ──(lnc_opt=.false.)───▶ │Vcmax25/Jmax25│                                            │
│Ali 2016 │ ◀╌╌╌ N 反馈 (待 lnc_opt) │  (v5 接通)    │                                            │
│Jmaxb1   │       ╎                  └──────┬───────┘                                            │
└─────────┘       ╎                         ▼                                                   │
                  ╎          ┌──────────────────────────────────┐                               │
                  ╎          │  biochemical.m (Farquhar-Collatz) │◀──── PAR/CO2/VPD ─────────────┘
   N_leaf 池 ─────╯          │  Ac/Aj/Ap 共限制 → A, Ag          │
   (光合容量)                │  Rd×0.67 白天光抑制 (v7)          │◀══ ψ_leaf 气孔限制 (水力胁迫)
                             │  Medlyn gs (g1=8.2 v10)           │            ▲
                             └──────┬───────────────┬───────────┘            ║ (水力上行)
                              GPP   │            Rd_leaf │  蒸腾 T ═══════════╝
                                    ▼                │   (gs 驱动失水)        ║
                  ┌────────────────────────────┐    │                        ║
                  │ computeNPP (Sitch 2003)     │    │            ┌───────────────────────┐
                  │  R_m,leaf = Rd_leaf ────────┼────┘            │ 根系水力 (calSPAC/V2) │
                  │  R_m,root = r·N_root·g(T)   │                 │  ψ_soil→ψ_root→ψ_stem  │
                  │  R_m,stem = r·N_stem·       │ f_livewood      │  →ψ_leaf (gsr 串联)    │
                  │            f_livewood·g(T)  │ =0.22(v6)       │  根吸水 → 蒸腾流 ═══════╝
                  │  R_g=0.25(GPP-R_m)          │                 └───────────▲───────────┘
                  │  NPP=GPP-R_m-R_g            │                             ║ 土壤水 ψ_soil
                  └──────────────┬─────────────┘                             ║ (RWU 根水吸取)
                             NPP │ (=bminc)                                   ║
                                 ▼                                  ┌─────────╨──────────┐
   ┌──────────────────────────────────────────────┐                │ 土壤 (STEMMUS 水热) │
   │ applyNSCCoupling (NSC cap + overdraft)        │                │  hh / TT / Theta_LL │
   │  R_a cap=min(R_a_pot,GPP+NSC) [curmr 等价]    │                └─────────┬──────────┘
   │   R_a_pot = R_m+R_g (来自 NPP)                │
   │  C_labile<0 → C_storage overdraft 借          │                          ║ 土壤无机 N
   │  (cushion 已删 v10)                           │                ┌─────────╨──────────┐
   └──────────────┬────────────────────────────────┘                │ 土壤 N (STEMMUS_N)  │
        NPP_actual│                                                  │ NH4/NO3+硝化/反硝   │
                  ▼                                                  │ +施肥/淋溶          │
   ┌──────────────────────────────────────────────┐                └─────────┬──────────┘
   │ C 主导分配(v12): A_C=A_N×actual_CN            │                          ║ 根吸 N (U_N)
   │  C先定NPP·f→D_N=C/CN→f_N→CN弹性[+10]→A_C跟随   │                          ║ DAISY demand+Nacry
   │ 【f_alloc 待§F动态化(物候/grass/pipe)】        │                ┌─────────╨──────────┐
   │  草→LPJ grass (leaf:root,ltor)                │                │ runPlantNDemand     │
   │  自然树→LPJ tree (pipe Eq2/4+ltor 功能平衡)   │  ltor 受 ψ 调  │ +applyRootUptake    │
   └──┬──────────┬──────────┬──────────┬───────────┘  (缺水多分根) └─────────┬──────────┘
      ▼          ▼          ▼          ▼                                      ║ U_N
   C_leaf   C_root_fine C_stem_wood C_fruit                                   ║
      │          │          │          │     ┌──────────────────────────┐    ║
      │          │          │          └────▶│ allocatePlantN: 根吸 U_N  │◀═══╝
      ▼(LAI)     │          │(拆sap/heart待) │ → 6 N 池                   │
┌──────────┐    │   ┌─────────────────────┐ │ (labile/leaf/root/stem/   │
│LAI 双引擎 │    │   │turnoverAndHarvest    │ │  fruit/storage)           │
│epic(校准) │    └──▶│ 落叶/根周转/采收/春萌 │ └────────────┬─────────────┘
│leafc_sla  │◀───────│turnoverCarbonMirror  │              │ N_leaf 池
│(C_leaf×SLA│        │ (N→C, CN_target)     │              ╰────────▶ 回 LUNA/光合容量 (闭环)
└──────────┘        └──────────┬───────────┘
                               ▼ litter → 土壤 SOM (师弟 M3: litter→矿化→无机N, 闭环)
                          回土壤 N 池 (释放无机 N 供再吸收)

══ 双线 = 水/N 上行闭环 (根→冠层); ── 单线 = 碳下行 (光合→分配→池)
闭环1 (水): 土壤水→根吸→ψ_leaf→气孔→蒸腾→失水→拉低ψ→再吸  (水力串联 gsr)
闭环2 (N): 土壤无机N→根吸U_N→N_leaf→光合容量(Vcmax)→GPP→分配→落叶litter→矿化→土壤N
闭环3 (C): GPP→NPP→各器官C池→turnover litter→SOM→R_h→大气, NEE=R_a+R_h−GPP
```

### §A.2 核心设计原则 (贯穿全链)

1. **物种属性借文献, 站点属性保标定**: C:N/Vcmax/异速系数 → CLM5/FruitTree/LPJ 通用值; GDD 物候 + LAI 校准 → 长武 25yr/Wang 247obs。
2. **守恒第一**: 每通量进 cum_* 累积, ΔΣpool ≡ Σin−Σout 机器精度。
3. **§5.5 反幻觉**: 每参数引用具体到 文件:行号 + 本地 PDF; 无出处标 [site-calibrated]。
4. **不破不立**: 新方案走 config 开关 (LAI_engine / alloc_scheme), 旧路径保留对比。
5. **苹果优先**: 苹果是当前唯一验证站点 (FruitTree 同物种最对口), 自然树通用方案后做。

---

## §B 物理方程 (按数据流逐步, 公式 + 单位 + 出处)

> **范围**: §B 只覆盖**植物 C/N/光合链**(LUNA→光合→呼吸→NSC→分配→物候/turnover→LAI→根吸N→NEE)。
> 土壤水热/能量/辐射/植物水力/土壤N转化(矿化/硝化/反硝化/ADE)等**非 C/N 链**过程不在此文档,
> 见各自模块文档。先读 §B.0 认识池子, 再看 §B.1+ 的逐步流程。

### §B.0 植物 C/N 池子清单 (6 C 池 + 6 N 池, prognostic state)

植物碳氮各 6 个**预报池 (prognostic state)**, 一一对应(同名器官 C 池与 N 池配对), 加 N_labile/C_labile
共 12 个 + N_labile 标量 = 13 个状态字段。**C 主导**: C 池是驱动量(C_leaf→LAI→GPP, 生物量好测),
N 池 = C/CN 跟随(§B.5)。储备池(storage)无单一 C:N, C/N 各自独立给值。

| 池 | 单位 | 苹果初值 | 功能 (这个池干什么) | 出处 |
|---|---|---|---|---|
| **C_leaf** | g C/m² | 89.56 | 叶碳; 经 SLA 决定 LAI→光合(§B.7); 落叶时进 litter | Wang 2009 |
| **C_root_fine** | g C/m² | 150.49 | 细根碳; 决定根吸水/N 能力; 周转 100% 进 litter | — |
| **C_stem_wood** | g C/m² | 1121.26 | 茎/木质碳(多年生主体, 量最大); 慢周转 | Tagliavini 2003 |
| **C_fruit** | g C/m² | 53.56 | 果实碳; DOY150-290 累积, 成熟(DS_fruit>0.95)一次性采收 90% | Mao 2018 |
| **C_labile** | g C/m² | 40.97 | **当日可支配碳缓冲池**; GPP−R_a 进, 分配出; 余碳/透支都在这里平账(§B.4/B.5); C 主导(config C_labile_init, N=C/CN_L 派生) | CLM5 |
| **C_storage** | g C/m² | 204.92 | **跨年越冬储备碳**(NSC); fcur 充值(§B.5)+ 落叶再吸收进, 春萌释放回 C_labile | FATES/CLM5 |
| **N_leaf** | g N/m² | =C_leaf/19.7 | 叶氮; 决定 Vcmax/Jmax(经 LUNA, §B.1)→光合容量; 衰老 retransl 70%回 storage | Wang 2009 |
| **N_root_fine** | g N/m² | =C/32 | 细根氮; 跟 C_root 同步周转 | — |
| **N_stem_wood** | g N/m² | =C/60 | 茎氮 | Tagliavini 2003 |
| **N_fruit** | g N/m² | =C/139 | 果实氮; 跟 C_fruit 同步采收 | Mao 2018 |
| **N_labile** | g N/m² | 2.08 | **当日可支配氮缓冲池**; 根吸 U_N 进(§B.6), 分配出; 决定 f_N_growth(§B.5) | M2-A |
| **N_storage** | g N/m² | 10.17 | **跨年越冬储备氮**; 落叶再吸收进, 春萌释放回 N_labile | CLM5 |

**池间物质流(一句话)**:
```
[进] 根吸 U_N → N_labile;  光合 GPP−R_a → C_labile
[分配] labile → 4 器官 (leaf/root/stem/fruit) + storage (fcur 越冬)      §B.5
[周转] 器官 → litter(出地, 给土壤SOM) / storage(叶N再吸收) / harvest(果实出地)  §B.8/B.9
[越冬] storage ←(秋,落叶再吸收) →(春,萌发释放) labile                    §B.8
[出口] litter→SOM(土壤N上游) · harvest(产量) · NEE=R_h+R_a−GPP(碳通量)    §B.10
```

- 出处: 池定义 config_vegetation.txt §[9b]+[init] (line 326-386); 状态推进 runPlantCNStep.m;
  初始化 setInitialNitrogen.m / setInitialPlantBGC. **C 主导初值**: 各器官给 C_*_init, N=C/CN 反算(§B.5)。
- **[一致性 2026-05-30]** C_labile 初值也已统一为 C 主导: config 只给 `C_labile_init`, `N_labile = C_labile / CN_L`
  反算(与器官池一致, setInitialNitrogen.m; 旧 `N_labile_init` 字段已删)。运行时两 labile 池本就各自独立演化
  (C 由 GPP−R_a, N 由根吸 U_N), 初值统一只为范式一致, 不影响守恒/物理。
- 守恒(机器精度), 两个残差诊断**均已输出**(conservation group, 对称):
  - `PlantC_residual` = cum_GPP − cum_R_a − litter_C − harvest_C − ΔC_pool  (v12 实测 abs_max 3.98e-13 g ✅)
  - **`PlantN_residual`** = **cum_U_N** − litter_N − harvest_N − ΔN_pool  (retransl 是池间转移不计; [N 守恒 2026-05-30] 新增)
  - 实现: cum_U_N/N_pool_init_const (setInitialNitrogen.m), 残差 (COMPASS.m), 注册 (registerPlantNitrogen.m)。

### §B.1 LUNA — 叶氮→光合容量 (每日)

```
lnc = 1 / (slatop × leafcn)                                          [g N / m² 叶]
JmaxCoef = Jmaxb1 × (hourpd/12)² × [1 − exp(−relhExp·(relh10−minrelh)/(1−minrelh))]
→ hill-climb 优化 → Vcmax25_z, Jmax25_z                              [μmol m⁻² s⁻¹]
```
**参数定义**:
| 符号 | 含义 | 单位 | 值/来源 |
|---|---|---|---|
| lnc | 叶氮浓度 (面积基) | g N m⁻² 叶 | =1/(slatop·leafcn) |
| slatop | 冠层顶比叶面积 | m² gC⁻¹ | 苹果 0.0125 (Wang 2009) |
| leafcn | 叶 C:N 比 | gC gN⁻¹ | 苹果 19.7 (v10 FruitTree) |
| Jmaxb1 | LUNA PAR-响应电子传递投资系数 | — | 苹果 0.244 [site-cal, Ali 2016 范围 0.05-0.30] |
| hourpd | 日长 | h | forcing 派生 |
| relh10 | 10 日均相对湿度 | — | forcing |
| relhExp, minrelh | LUNA 湿度响应形状/下限 | — | CLM5 LunaMod 常数 |
| Vcmax25_z, Jmax25_z | 25℃ 最大羧化/电子传递速率 (LUNA 优化输出) | μmol m⁻² s⁻¹ | 每日更新 |
- 出处: Climate10dBuffer.m:234 (lnc); nitrogenAllocation.m:58 (JmaxCoef); CLM5 LunaMod.F90; Ali et al. 2016 GMD 9:587
- ⚠️ lnc_opt=.false.: lnc 是 PFT 常数, **土壤 N 暂不实际反馈光合** (待办 §E P2, 阶段15)
- 接入: leafbio.Vcmo + use_luna/Jmax25_LUNA 注入 biochem (v5, ebal_clm45.m:1106)

### §B.2 光合 (Farquhar-Collatz + Medlyn) — biochemical.m

**动力学常数温度响应** (Bernacchi 2001, biochemical.m:340-342):
```
Kc = 404.9e-6·p·f_Arr(79430,T)   Ko = 278.4e-3·p·f_Arr(36380,T)   Γ* = 42.75e-6·p·f_Arr(37830,T)   [bar]
f_Arr(Ha,T) = exp[ Ha/(Tref·R)·(1 − Tref/T) ]
```

**Vcmax/Jmax/TPU/Rd 温度响应** (peaked Arrhenius, biochemical.m:355-358; CLM5 Table 9.1):
```
X(T) = X25·f_Arr(Ha)·f_Peak(Hd,ΔS,T)
f_Peak(Hd,ΔS,T) = [1+exp((Tref·ΔS−Hd)/(Tref·R))] / [1+exp((T·ΔS−Hd)/(R·T))]
```
| 量 | Ha | Hd | ΔS | 25℃ 比例 |
|---|---|---|---|---|
| Vcmax | 72000 | 200000 | 668.39−1.07·Tgr (动态) | — |
| Jmax | 50000 | 200000 | 659.70−0.75·Tgr (动态, Kattge 2007 热适应) | — |
| TPU | 72000 | 200000 | 490 | TPU25 = **0.167·Vcmax25** (v3 修, CLM5 PhotosynthesisMod.F90:1405+1517) |
| Rd | 46390 | 150650 | 490 | Rd25 = 0.015·Vcmax25 (C3) |

**叶暗呼吸白天光抑制** (Kok, v7, biochemical.m:367; CLM5 1562, Lloyd 2010 + Metcalfe 2012):
```
Rd = Rd × (0.67 + 0.33·𝟙[Q≤0])          Q>0白天×0.67, 夜间×1.0
```

**电子传递** (biochemical.m:397-402; θ_psii=0.7 CLM5 1255):
```
Jpot = 0.5·po0·Q·fPAR     po0 ≈ 0.80
θ_J·Jact² − (Jpot+Jmax)·Jact + Jpot·Jmax = 0     θ_J=0.7, 取小根
```

**三限制共限制** (Collatz 1991 / CLM5; biochemical.m:539-575):
```
Vc = Vcmax·(Ci−Γ*)/(Ci+Kc(1+O/Ko))                Rubisco
Ve = Jact·(Ci−Γ*)/(Ci+2Γ*)·effcon                 光, effcon=1/4
Vs = 3·TPU                                          TPU
V  = smooth_min(Vc,Ve; θ_cj=0.98)   Ag = smooth_min(V,Vs; θ_ip=0.95)   A = Ag − Rd
```

**Medlyn 气孔 + Ci 迭代** (biochemical.m:583-593,441-454; Medlyn 2011; g1 FruitTree):
```
gs = g0 + 1.6·(1+g1/√VPD)·A/Cs        g1 = 8.2 [kPa^0.5] (v10, Zanotelli 实测)
Ci = Cs − 1.6·A/gs                     阻尼定点迭代 (damp 0.5, ~10 步)
```
边界层 Ca→Cs 在 ebal heatfluxes_clm45.m:107。**输出**: GPP=Ag·LAI·12e-9, Rd_leaf=Rdtot·12e-9 [kgC m⁻²s⁻¹] (ebal_clm45.m:990)。

**§B.2 参数定义**:
| 符号 | 含义 | 单位 | 值/来源 |
|---|---|---|---|
| Kc, Ko | CO₂ / O₂ Michaelis-Menten 常数 | bar | 25℃: 404.9 / 278.4e3 μmol mol⁻¹ (Bernacchi 2001) |
| Γ* | CO₂ 补偿点 (无暗呼吸) | bar | 25℃ 42.75 μmol mol⁻¹ |
| Ha | Arrhenius 激活能 | J mol⁻¹ | Vcmax/Jmax/TPU/Rd: 72000/50000/72000/46390 |
| Hd | 高温去活化能 | J mol⁻¹ | 200000/200000/200000/150650 |
| ΔS | 熵项 (热适应) | J mol⁻¹ K⁻¹ | Vcmax/Jmax 动态 (Kattge 2007); TPU/Rd=490 |
| Tref | 参考温度 | K | 298.15 |
| R | 气体常数 | J mol⁻¹ K⁻¹ | 8.314 |
| Tgr | 生长温度 (热适应输入) | ℃ | 10 日均气温 |
| Q | 吸收 PAR (APAR) | μmol photon m⁻² s⁻¹ | forcing/RTM |
| po0 | 最大暗光化学效率 Fv/Fm | — | Kp/(Kf+Kd+Kp)≈0.80 |
| Jpot, Jact | 潜在/实际电子传递速率 | μmol e⁻ m⁻² s⁻¹ | — |
| θ_J | 电子传递曲率 | — | 0.7 (CLM5 θ_psii) |
| effcon | CO₂ : 电子 转换效率 | mol mol⁻¹ | 1/4 (C3) |
| θ_cj, θ_ip | Vc-Ve / V-Vs 共限制曲率 | — | 0.98 / 0.95 (Collatz/CLM5) |
| Ci, Cs, Ca | 胞间/叶面/大气 CO₂ | ppm (或 bar) | Ci 迭代解 |
| gs, g0, g1 | 气孔导度 / 残余 / Medlyn 斜率 | mol m⁻² s⁻¹ / mol m⁻² s⁻¹ / kPa^0.5 | g0=1e-4, g1=8.2 (v10) |
| VPD | 叶面饱和水汽压差 | kPa | forcing 派生 |
| Vcmax25/Jmax25/TPU25/Rd25 | 各 25℃ 基值 | μmol m⁻² s⁻¹ | LUNA 或 PFT 表 |
| Ag, A | 总 / 净 光合 | μmol CO₂ m⁻² s⁻¹ | A=Ag−Rd |

### §B.3 维持/生长呼吸 + NPP (Sitch 2003 LPJ) — computeNPP.m

```
g(T) = exp[ 308.56·(1/56.02 − 1/(T_℃+46.02)) ]                    (Lloyd&Taylor 1994, Sitch Eq23)
R_m,leaf = Rd_leaf                                                  (避免双计)
R_m,root = r_m10·N_root·g(T_soil)
R_m,stem = r_m10·(N_stem·f_livewood)·g(T_air)                     只活木呼吸 (v6)
R_g = 0.25·max(GPP−R_m,0)     NPP = max(GPP−R_m−R_g,0)     R_a = R_m+R_g
```
**参数定义**:
| 符号 | 含义 | 单位 | 值/来源 |
|---|---|---|---|
| g(T) | Lloyd-Taylor 温度响应 | — | Sitch Eq23 (T clamp [−40,50]℃) |
| r_m10 | 10℃ 维持呼吸基率 | gC (gN·d)⁻¹ | 0.066 (Sitch Table3) |
| N_root, N_stem | 细根/树干 N 池 | g N m⁻² | prognostic |
| f_livewood | 活木(边材)占树干比例 | — | 苹果 0.22 [site-cal, CLM5 CNMRespMod:265 livestem 只活木] |
| T_air, T_soil | 气温 / 根区土温 | ℃ | forcing / 土壤 |
| R_m,{leaf,root,stem} | 各器官维持呼吸 | kgC m⁻² s⁻¹ | — |
| R_g | 生长呼吸 | kgC m⁻² s⁻¹ | =0.25·max(GPP−R_m,0), Ryan 1991 |
| R_a | 总自养呼吸 | kgC m⁻² s⁻¹ | =R_m+R_g |
| NPP | 净初级生产力 | kgC m⁻² s⁻¹ | =max(GPP−R_m−R_g,0) |
- 出处: computeNPP.m:92-114; Sitch 2003 GCB 9:161 Eq22-24

### §B.4 NSC 联动 (applyNSCCoupling.m)

```
R_a_potential = R_m + R_g           ← 来自 §B.3 (总自养呼吸 potential, cap 前)
NSC_avail = C_labile + C_storage
R_a_actual = min(R_a_potential, GPP·dt + NSC_avail)              碳赤字下调呼吸 (=CLM5 curmr_ratio)
NPP_actual = max(GPP·dt − R_a_actual, 0)
alloc 后 C_labile<0 → C_storage overdraft 借 (池间转移, 守恒)
```
**参数定义**:
| 符号 | 含义 | 单位 | 备注 |
|---|---|---|---|
| C_labile | 流动 NSC (现金, 过路池) | g C m⁻² | prognostic |
| C_storage | 储备 NSC (淀粉) | g C m⁻² | prognostic |
| NSC_avail | 可用 NSC = C_labile+C_storage | g C m⁻² | — |
| R_a_potential | 理论**总自养**呼吸 (cap 前) = R_m+R_g | g C m⁻²/dt | 来自 §B.3 computeNPP (R_m_leaf+root+stem + R_g); 代码 applyNSCCoupling.m:81 (v10 改名, 旧名 R_m_potential 误导) |
| R_a_actual | NSC-capped 实际呼吸 | g C m⁻²/dt | =min(R_a_potential, GPP·dt+NSC) |
| dt | 时间步 | s | 3600 |
| storage_drawdown_event | 动用储备标记 (0/1) | — | 𝟙[R_a_actual<R_a_potential], **非致死** |
- 出处: applyNSCCoupling.m:85-116; CLM5 NutrientCompetitionCLM45defaultMod.F90:709-720 (curmr); Sala 2012/Hartmann 2016
- FATES cmort (真碳饥饿致死, 基于 store/target_leaf 比) 对比见附录 §12.31; COMPASS 不接致死 (单点大叶, 人工果园)

### §B.5 C/N 分配 — C 主导 + C:N 弹性 (v12 2026-05-29, computeOrganDemand + allocatePlantN + computeCarbonAllocation)

**范式: C 主导 (carbon-dominant), N 跟随, C:N 在 [CN_target, CN_max] 弹性浮动。** 对标 CLM5 FlexibleCN
(`plant_calloc = plant_nalloc × c/n_allometry`, NutrientCompetitionFlexibleCNMod.F90:477; nscalar
弹性 :1660-1667) + FruitTree (Dabbaghi 2022 p.4 "C and corresponding N")。**不是 N 主导** (历史
v12 初版误做 A_C=A_N×CN_target 已纠正)。

```
[C 先定]   A_C_target_organ = NPP_actual × f_alloc_organ      Σf_alloc=1; 果实门 DOY150-290
[N 需求]   D_N_organ        = A_C_target_organ / CN_target_organ
[N 受限]   f_N_growth       = clamp(N_labile_supply / Σ(D_N·dt), 0, 1)
           A_N_organ        = D_N_organ × f_N_growth × dt       实际 N 分配 (N 守恒)
[C:N 弹性] actual_CN_organ  = min(CN_target_organ / f_N_growth, CN_target_organ + offset)
[C 跟随]   A_C_organ        = A_N_organ × actual_CN_organ
           ├ N 足 (f_N=1)        → actual_CN=CN_target → A_C=A_C_target  (全额)
           ├ N 微缺 (未撞顶)      → C:N 上浮, A_C≈A_C_target            (C 照长, 每 C 少配 N)
           └ N 严重缺 (撞 CN_max) → A_C=A_N·CN_max < A_C_target          (C 被砍, 余碳留 C_labile)
守恒: ΔC_labile + ΣΔC_organ = GPP·dt − R_a·dt (机器精度); ΣA_N + N_labile_new = N_labile_supply (N 守恒)
```

**参数定义**:
| 符号 | 含义 | 单位 | 值/来源 |
|---|---|---|---|
| f_alloc_organ | NPP→各器官分配分数 (C 先定) | — | 苹果 leaf0.25/root0.25/stem0.20/fruit0.30 [site], Σ=1 |
| CN_target_organ | 各器官目标 C:N | gC gN⁻¹ | FruitTree Zanotelli 2010 实测: leaf19.7/root32/stem60/fruit139 (Dabbaghi 2022 Table C1) |
| cn_flex_offset | C:N 弹性半幅 (CN_max=target+offset) | gC gN⁻¹ | 10 (CLM5 FlexibleCN:1660-1661 leafcn±10) |
| f_N_growth | 全局 N 限制因子 (0-1) | — | N_labile_supply / ΣD_N (allocatePlantN) |
| actual_CN_organ | 弹性实际 C:N (N 不足时上浮) | gC gN⁻¹ | clamp(CN_target/f_N, [CN_target, CN_max]) |
| 果实门 | fruit 分配仅 DOY150-290 | — | 非果实期 fruit 重分配 leaf(60%)/root(40%) |
- 出处: computeOrganDemand.m:82-104 (C 先定→D_N); allocatePlantN.m:85-110 (f_N + actual_CN 弹性);
  computeCarbonAllocation.m:73-83 (A_C=A_N×actual_CN); config cn_coupling_mode=c_follows_n / cn_flex_offset=10
- 物理: N 胁迫下叶 C:N 升高 (单位叶含氮↓→光合 N 限制经 LUNA 反馈), C 照长不立即减叶量 (保护 LAI 校准),
  仅 C:N 撞生理上限 CN_max 才弃碳。CLM5 nscalar 同理。
> 单一真相: f_alloc 仍固定 (功能平衡/异速待 §F 动态分配); C:N 弹性已根治旧 v10 无界漂移 (82-92→[target,max])。
> **→ 分配后数据流分三支**: N 侧供给来自 §B.6 根吸; 叶碳 C_leaf 经 §B.7 转 LAI 回喂光合; 器官池经 §B.8 turnover 流出。

### §B.6 根吸 N (上行供给, runPlantNDemand + applyRootUptake)

> ← 承 §B.5: 分配的 N 侧供给来自 N_labile; 本节讲 N_labile 怎么从土壤吸收补充。

```
[需求 need-based] demand = D_N_growth + max(0, N_labile_target − N_labile) / dt
   D_N_growth      = Σ_organ NPP·f_alloc_organ / CN_target_organ        各器官组织建造 N (= §B.5 的 D_N)
   N_labile_target = (C_leaf/CN_L + C_root/CN_target_root) · nrelocfrac  可支配池储备目标 (够重建明年叶+细根)
[吸收] Nacry 根膜双底物: U = Vmax·C_root/(Km+C_root) − Cmin floor;  DAISY 单根 ADE C_root 二分法解
       实际吸收 ≤ demand 且 ≤ 土壤供给 (供给不足 demand_saturated 封顶)
U_N → N_labile → (供 §B.5 allocatePlantN 分 5 器官) → N_leaf 回 §B.1 LUNA 光合容量 (闭环2)
```
**参数定义**:
| 符号 | 含义 | 单位 | 值/来源 |
|---|---|---|---|
| D_N_growth | 各器官生长 N 需求 (NPP·f_alloc/CN_target) | g N m⁻²/dt | computeOrganDemand.m:89 |
| N_labile_target | 可支配池储备目标 (重建明年叶+细根 N) | g N m⁻² | (C_leaf/CN_L+C_root/CN_root)·nrelocfrac; LPJ-GUESS growth.cpp:1479 (Smith et al. 2014 GMD 7:2261) |
| nrelocfrac | 衰老再吸收分数 | — | 苹果 0.7 / 草作物 0.5 (PFT 表) |
| U / U_N | 单底物/总根吸 N 速率 | g N m⁻²/dt | Michaelis-Menten |
| Vmax | 最大吸收速率 (NH4/NO3 各一) | μmol g⁻¹根 h⁻¹ | Nacry 根膜 (PFT 表) |
| Km | 半饱和常数 | mg/L | NH4/NO3 各一 |
| Cmin | 根膜物理下限 floor | mg/L | Nacry |
| C_root | 根表无机 N 浓度 | mg/L | 二分法解 (DAISY 单根 ADE) |
| SRL | 比根长 | m根 gC⁻¹ | 苹果 ~40 (Wells 2002) |
- 出处: runPlantNDemand.m, applyRootUptake.m, +PlantNitrogen/; LPJ-GUESS need-based (Smith 2014) + Nacry 2013 (Plant Soil 370:1) + DAISY 单根 ADE + DSSAT NUPTAK。

### §B.7 LAI (叶碳→冠层, 双引擎 computeLAIfromPhen.m)

> ← 承 §B.5: 分配出的叶碳 C_leaf 经 SLA 转成 LAI; LAI → 回喂 §B.2 光合 (上行闭环1)。

```
epic_phen_mix (默认,校准): LAI = (1−phen)·ALAI_MIN + phen·EPIC_curve   [247obs, RMSE 0.348]
leafc_sla_phen (碳驱动):   LAI = (1−phen)·ALAI_MIN + phen·(SLA0/SLAm)·(exp(C_leaf·SLAm)−1)
SLA(LAI) = SLA0 + SLAm·LAI                                            (CLM5 Thornton 2007)
```
**参数定义**:
| 符号 | 含义 | 单位 | 值/来源 |
|---|---|---|---|
| phen | 物候活跃度 (0=休眠,1=满叶) | — | Phase D 物候 (§B.9) |
| ALAI_MIN | 落叶期最小 LAI | m²/m² | 苹果 0.293 |
| SLA0 | 冠层顶比叶面积 | m² gC⁻¹ | 0.0125 (Wang; vs FruitTree 0.028 待核) |
| SLAm | SLA 随 LAI 梯度 dSLA/dLAI | m² gC⁻¹ /LAI | 0.00133 (CLM5) |
| C_leaf | 叶碳池 (面积基) | g C m⁻² | prognostic / =N_leaf·CN_L |
| EPIC_curve | 247-obs 校准 LAI 曲线 | m²/m² | RMSE 0.348 |
- 出处: computeLAIfromPhen.m:50-76; Sitch 2003 Eq5; CLM5 Thornton 2007 Eq3 (SLA 梯度)。两引擎并存 (epic 保校准, leafc_sla 接动态分配 §F)。

### §B.8 Turnover / 再吸收 / 采收 (器官出口, Phase D 物候)

> ← 承 §B.5+§B.7: 长成的器官 C/N 池经此流出 — litter(给土壤)/storage(越冬)/harvest(产量); 时序由 §B.9 物候管。

```
叶衰老 (CLM5 strict): rate=1/ndays_off≈6.7%/d; offset_counter≤1 末日全清
                       litter = d_N_leaf·(1−nrelocfrac); retransl → N_storage
根周转 root_turnover_rate (受物候门控 phen_gate, F5); 苹果采收 N_fruit·harvest_fraction 带走;
春萌 N_storage→N_labile (1/spring_release_days)
```
**参数定义**:
| 符号 | 含义 | 单位 | 值/来源 |
|---|---|---|---|
| ndays_off | CLM5 落叶 ramp 时长 | d | 苹果 15 (CLM5) |
| nrelocfrac | 衰老叶 N 再吸收比例 | — | 苹果 0.7 (Aerts 1996) |
| root_turnover_rate | 细根周转率 (×phen_gate 物候门控, F5) | yr⁻¹ | PFT 表; 门控 CLM5 bglfr |
| harvest_fraction | 采收带走果实 N 比例 | — | 苹果 0.90 (商业采摘) |
| spring_release_days | 春萌储备释放时长 | d | 苹果 60 |
| **DS_fruit** | **果实独立发育期** (驱动 harvest) | — | bloom 后 GDD/TSum_fruit; >0.95 触发收获 |
- harvest 由果实独立时钟 DS_fruit>0.95 触发 (与叶物候 frPHU 解耦, seasonDecidPhenology.m:315-325)。
- 出处: leafSenescenceRetransl.m (CLM5 CNPhenologyMod.F90:2450), turnoverAndHarvest.m, fineRootTurnover.m (F5 门控), springStorageRelease.m; 详见 phaseD_phenology_audit.md。
- **litter→SOM 出口** (exportLitterToSOM.m): leaf/root litter 的 C/N 打包交土壤 SOM 模块, 是 C/N 链交给土壤 N 的**唯一下游接口** (→ 续 §B.10 NEE 的 R_h)。

### §B.9 物候相位机 (决定何时长叶/落叶/采收, +integrated/+phenology/)

物候是 C/N 链的**时间开关**: 不直接动 C/N 池, 而是维护一个 **per-PFT 状态机** (PhenologyState, 23 个状态变量),
逐日推进, 产出 `phen∈[0,1]` (叶活跃度) + 事件 flag (onset/offset/harvested), 再驱动 §B.7 LAI 和 §B.8 turnover。
**CLM5 filter+flag-gate 架构**: 每 PFT 4 个二元物候 flag 恰好 1 个=1, dispatcher (phenologyDispatch.m) 路由到
对应分支。**站点级**: 每 PFT 各持一个 state struct (CLM5 是 _patch 数组, COMPASS 单点每 PFT 一份)。

```
                    ┌─────────────── phenologyDispatch (按 4 flag 路由) ───────────────┐
                    ▼              ▼                    ▼                     ▼
              evergreen      season_decid          stress_decid           crop
              (phen≡1)    (冷-GDD 耦合, 苹果)     (4 累积器 冷+旱)        (DS, 占位)
                              │
   ┌──────────────────────────┼──────────────────────────────────────┐
   │ 休眠 dormant ──冷月reset──▶ 萌发 onset ──Sykes ramp──▶ 满叶 ──日长──▶ 落叶 offset ──▶ 休眠
   │ chilldays↑gdd5=0          gdd5>gdd0          phen→1     dayl<crit    phen→0,落叶retransl
   └──────────────────────────────────────────────────────────────────┘
        └─(苹果)果实子时钟 budburst→bloom→DS_fruit→harvest, 与叶物候解耦
```

#### PhenologyState 状态变量 (23 个, PhenologyState.m)
| 变量 | 含义 | 类型/单位 | 出处 |
|---|---|---|---|
| chilldays | 累积冷日 (T_avg<5°C) | days | LPJ driver.cpp:550 |
| gdd5 | 累积 GDD (base 5°C) | °C·d | LPJ driver.cpp:548 |
| mtemp_last | 上一日月均温 (冷月 reset 用) | °C | — |
| ifsensechill | 最暖日 latch (gate 冷月 reset, 防早 reset) | bool | LPJ driver.cpp:540-545 |
| prev_dayl / ws_flag | 昨日日长 / 冬至-夏至升降标志 | h / 0,1 | CLM5 :802-807 |
| onset_gddflag / onset_gdd | 春 GDD 累积器开关 / 累积值 (stress 用) | 0,1 / °C·d | CLM5 :718-719 |
| **phen** | **叶活跃度 (0 休眠→1 满叶)** | — | LPJ growth.cpp:96-99 |
| aphen | 本季累积叶龄 (Σphen), 封顶 210 | days | LPJ growth.cpp:150 |
| onset_flag / onset_counter | 萌发期标志 / ramp 倒计时 | bool / days | CLM5 CNVegState:81 |
| offset_flag / offset_counter | 落叶期标志 / ramp 倒计时 | bool / days | CLM5 :87 |
| dormant_flag | 冬眠标志 (初始 true) | bool | CLM5 :79 |
| budburst_doy / senescence_doy | 萌发/落叶起始 DOY (诊断) | DOY | — |
| gdd5_at_budburst / bloom_doy / gdd5_at_bloom | 果实子时钟锚点 (苹果) | °C·d/DOY | seasonDecid:302-311 |
| **DS_fruit** / harvested_flag / harvest_doy | 果实发育期 / 已采标志 / 采收 DOY | — / bool / DOY | seasonDecid:315-326 |

#### 分支 1: season_decid — 季节落叶 (苹果 23 / 落叶树 3,7,8,11 / 落叶灌 10) — **苹果主用**
Sykes 1996 chilling-GDD 耦合萌发 + CLM5 日长落叶 (seasonDecidPhenology.m, LPJ growth.cpp:82-101 逐字移植):
```
① GDD/冷日累积:   gdd5 += max(0, T_avg−5);  if T_avg<5: chilldays += 1
② 最暖日 latch:    if mtemp == year_mtemp_max → ifsensechill = true   (driver.cpp:540)
③ 冷月 reset:      if mtemp_last≥5 且 mtemp<5 且 ifsensechill →
                      gdd5=0; chilldays=0; dormant_flag=true; offset_flag=false   (driver.cpp:583)
④ Sykes 阈值+萌发:  gdd0 = k_chilla + k_chillb·exp(−k_chillk·chilldays)         ← 冷量越多, 萌发所需 GDD 越少
                    if gdd5>gdd0 且 aphen<210:  phen = min(1, (gdd5−gdd0)/phengdd5ramp)
                       若 dormant: onset_flag=true; budburst_doy=doy
                    else: phen = 0
⑤ 日长落叶:        if dayl<crit_dayl_h 且 非dormant 且 aphen>30 → offset_flag=true (落叶+N再吸收 §B.8)
   aphen += phen   (累积叶龄, 最冷日 reset)
```
| 参数 | 含义 | 苹果值/来源 |
|---|---|---|
| k_chilla/b/k | Sykes 1996 冷量-GDD 系数 (a+b·exp(−k·chill)) | LPJ europe.ins; 苹果高需冷 k_chillb=1150 |
| phengdd5ramp | 萌发后 phen 0→1 的 GDD ramp 宽度 | °C·d (LPJ) |
| crit_dayl_h | 落叶日长阈值 | 10.9 h (CLM5 crit_dayl=39300s 全局) |
| APHEN_MAX / APHEN_GROW_MIN | 叶龄上限 210d / 落叶最小叶龄 30d | LPJ growth.cpp:50 |
- **半球感知** (LPJ guess.h:178): 最冷日 N 半球 DOY14 / S 半球 DOY195, 按 `forcing.latitude_deg` 符号切; aphen/果实状态在最冷日 reset → 跨年循环。

#### 分支 2: stress_decid — 胁迫落叶 (C3/C4 草 13,14 / 热带雨绿 6 / 落叶灌 12)
CLM5 CNStressDecidPhenology 4 累积器 (stressDecidPhenology.m, CNPhenologyMod.F90:974-1417 逐字移植):
```
① onset_fdd:  冻日累积 (dormant 且 onset_gddflag=0 且 T_soil<0)         → >15d latch onset_gddflag=1
② onset_gdd:  春 GDD 累积 (onset_gddflag=1 且 T_soil>0)
③ onset_swi:  土壤湿日累积 (soil_psi > −0.6 MPa)                        → >15d
④ 萌发条件:   onset_swi>15 且 prec10>20mm 且 (gddflag=0 或 gdd够) 且 dayl>6h → onset, ramp ndays_on
   落叶:      dayl<crit_dayl_h  或  严重干旱 (psi < −2.0 MPa)
```
常数 (CLM5 scalar): crit_onset_fdd/swi=15d, soilpsi_on=−0.6/off=−2.0 MPa, rain_thresh=20mm, dayl_min=6h。
⚠️ 苹果站未启用 (PFT 23 是 season_decid); 待 C3 草 / 多 PFT 共生站点验证。

#### 分支 3: evergreen (针叶常绿 1,2 / 阔叶常绿 4,5 / 常绿灌 9)
phen≡1 恒定, 无萌发/落叶事件; LAI 靠 leaf_long 连续周转驱动 (evergreenPhenology.m; CLM5 :502-666)。

#### 分支 4: crop (作物 15-22) — **占位**
年播种→DS 累积→收获 (cropPhenology.m, DAISY/LPJ 风格); 苹果站未启用, 待作物 PFT 站点完整实现。

#### 苹果果实独立子时钟 (DS_fruit, P5 fix 2026-05-26) — season_decid 内嵌, 与叶物候**解耦**
```
萌发(budburst) ─GDD_bloom(≈250)─▶ 开花(bloom) ─▶ DS_fruit = GDD_since_bloom / TSum_fruit(≈2132)
harvest: DS_fruit > harvest_DS(0.95) → 一次性触发 harvested_flag        seasonDecidPhenology.m:315-327
reset:   最冷日清零 gdd5_at_budburst/bloom_doy/DS_fruit/harvested_flag → 跨年
```
| 符号 | 含义 | 苹果值/来源 |
|---|---|---|
| GDD_bloom_after_budburst | 萌发→开花热需求 | ≈250 °C·d (Wang 2009 DOY80→130) |
| TSum_fruit | 开花→成熟热需求 | ≈2132 °C·d (Wang 2009 DOY130→285) |
| harvest_DS | 成熟采收阈值 | 0.95 (DAISY DS=2 + EPIC harvest index 0.9-0.95 中值) |
- **[关键修复]** harvest 改 DS_fruit 驱动后与 frPHU 天花板**彻底解耦** (旧版看 frPHU>0.90, 长武 max=0.648 永不触发);
  v12 实跑 cum_N_harvest=0.80 / cum_C_harvest=111.3 g/m² (=果实峰值 90%=harvest_fraction) ✅。详见 §10b.4 + §12.27.6。
- 出处: phenologyDispatch.m / seasonDecidPhenology.m / stressDecidPhenology.m / evergreenPhenology.m /
  cropPhenology.m / PhenologyState.m; 文献 LPJ-GUESS growth.cpp+driver.cpp (Sykes 1996) + CLM5 CNPhenologyMod.F90;
  详见 phaseD_phenology_audit.md。⚠️ 仅 season_decid (苹果) 实战验证过, 其余 3 分支待对应 PFT 站点验证。

### §B.10 NEE 碳通量闭合 (computeNEE.m, 碳环出口)

碳环最终收口: 整个生态系统对大气的净 CO₂ 交换。
```
NEE = R_h + R_a_plant − GPP            [kg C/m²/s]   (Chapin 2002 §6.4)
    = (R_a_plant − GPP) − (−R_h)
    = −NPP_eco + R_h
符号: NEE>0 生态系统是 CO₂ 源(净排放); NEE<0 是汇(净吸收)
```
| 符号 | 含义 | 单位 | 来源 |
|---|---|---|---|
| GPP | 冠层总初级生产 | kg C/m²/s | SCOPE biochem(§B.2) |
| R_a_plant | 植物自养呼吸 (R_m+R_g) | kg C/m²/s | computeNPP(§B.3) |
| R_h | 异养(土壤微生物)呼吸 | kg C/m²/s | fluxes.Resp(土壤模块, 师弟) |
- 出处: computeNEE.m; Chapin 2002 Ch.6 §6.4; CLM5 Tech Note §15.1。Phase C M2-B.1b S3 从 ebal.m 独立
  (旧版 NEE 漏了 R_a_plant, 现 LPJ NPP 提供完整 R_a 才能闭合)。
- 注: R_h 来自土壤 SOM 分解(非本文档范围), litter→SOM 入口见 §B.8 末 exportLitterToSOM。

---

## §B.11 🍎 苹果专用公式汇总 (FruitTree port, PFT23)

> 苹果 (PFT23) 相对通用 PFT 额外的 FruitTree 物理过程, **完全照搬 CLM5-FruitTree + 适配 COMPASS 架构**。
> 源: Dombrowski & Hoffmann 2022 GMD 15:5167; `D:/OneDrive/CTSM-release-clm5.0_FruitTree/`。
> 全部 **apple-scoped** (门控 `use_fruittree_alloc` / `use_cesaraccio_chill` / `has_storage_pool`+leafc_sla), 其他 PFT 不变。

### §B.11.1 4 相位 GDD 分配 + J-Norman 叶锥 (Group A, fruitTreeAllocation.m)
设 `g = gdd_budburst` (自 budburst 的 GDD, 基温 5°C)。器官分配比 a_leaf/a_root/a_stem/a_repr (∑=1):

**(0) 萌发前** (g 非有限 或 g ≤ 0): a_leaf = a_root = a_stem = a_repr = 0 (休眠门控)。

**(1) 采后/成熟** (harvested = 1 或 g ≥ hybgdd) — 储备回填 (CLM5-FruitTree :1742-1747):
```
a_leaf = aleafstor (0.3)
a_stem = a_root = (1 − a_leaf) / 2 = 0.35
a_repr = 0
```

**(2) Phase 1** (0 < g < grnfill) — 建冠, 无果实, prog = g / grnfill:
```
a_root = arooti − (arooti − arootf)·prog                    [0.7 → 0.2 线性]
a_leaf = max(0.05, (1 − a_root)·fleafi·τ_leaf(g))           [CLM5-FruitTree :1666 地板 0.05]
a_stem = 1 − a_leaf − a_root                                [= 余]
a_repr = 0
```

**(3) Phase 2-4** (grnfill ≤ g < hybgdd) — 茎衰减, 果实取余, prog2 = (g − grnfill)/(hybgdd − grnfill):
```
a_root = arootf − (arootf − arootf2)·prog2                  [0.2 → 0.08 线性]
a_leaf = max(1e−5, (1 − a_root)·fleafi·τ_leaf(g))           [taper → ~0 过 lfmat]
aleaf_p1end = max(0.05, (1 − arootf)·fleafi·τ_leaf(grnfill))   [Phase 1 末叶分配, 解析]
astemi      = 1 − aleaf_p1end − arootf                      [Phase 1 末茎分配]
a_stem = max(astemf, astemi − (astemi − astemf)·prog2)      [→ astemf=0.22]
a_repr = max(0, 1 − a_root − a_stem − a_leaf)               [果实 = 余, 随叶/茎关闭而升]
```

**(4) 归一**: tot = a_leaf+a_root+a_stem+a_repr; 若 tot>0 则各 a_* /= tot (保 ∑=1)。

**J-Norman 叶锥 τ_leaf(g)** (leafGDDtaper.m; Dombrowski 2022 GMD 15:5167 Eq.4; CLM5-FruitTree NutrientCompetitionFlexibleCNMod.F90:1663-1666), 归一到 g=0 时 1、g=gdd_lfmat 时 0:
```
g_c = clamp(g, 0, gdd_lfmat) ;   b = ft_bfact (−0.5) ;   eb = exp(−b)
τ_leaf = (eb − exp(−b·g_c/gdd_lfmat)) / (eb − 1)           [b≈0 → 线性 1 − g_c/gdd_lfmat]
τ_leaf = clamp(τ_leaf, 0, 1)
```
- 参数 (FruitTree nc citrus idx35): arooti=0.7, arootf=0.2, arootf2=0.08, fleafi=0.85, astemf=0.22, aleafstor=0.3, ft_bfact=−0.5, ft_GDDlfmat(=gdd_lfmat)=1350, grnfill=400, hybgdd=2880 [°C·d]。
- 注: a_leaf 峰 0.255 = (1−arooti)·fleafi 是 4 相位耦合结果, 非单参数 [[project_compass_fruittree_alloc_params]]。

### §B.11.2 Cesaraccio 三角冷激 (Group B, cesaraccioChillDay.m)
逐日给 Cd (chill-day, ≤0) 与 Ca (anti-chill-day, ≥0) 增量, dT = Tmax − Tmin, baset=4。
**完整 6 分支** (CLM5 CNPhenologyMod.F90:1824-1862; Cesaraccio et al. 2004 IJB 48:67):
```
dT ≤ 0 (无日较差, 退化为均温):
   Tmean ≤ baset 且 ≥ 0:  Cd = −(baset − Tmean)
   Tmean > baset:         Ca = Tmean − baset

Cd (chill-day):
   0 ≤ baset ≤ Tmin ≤ Tmax:        Cd = 0
   0 ≤ Tmin ≤ baset < Tmax:        Cd = −[ (Tmean − Tmin) − (Tmax − baset)²/(2·dT) ]
   0 ≤ Tmin ≤ Tmax ≤ baset:        Cd = −(Tmean − Tmin)
   Tmin < 0 ≤ Tmax ≤ baset:        Cd = −Tmax²/(2·dT)
   Tmin < 0 < baset < Tmax:        Cd = −Tmax²/(2·dT) − (Tmax − baset)²/(2·dT)
   else:                           Cd = 0

Ca (anti-chill-day):
   0 ≤ baset ≤ Tmin ≤ Tmax:        Ca = Tmean − baset
   0 ≤ Tmin ≤ baset < Tmax:        Ca = (Tmax − baset)²/(2·dT)
   0 ≤ Tmin ≤ Tmax ≤ baset:        Ca = 0
   Tmin < 0 ≤ Tmax ≤ baset:        Ca = 0
   Tmin < 0 < baset < Tmax:        Ca = (Tmax − baset)²/(2·dT)
   else:                           Ca = 0
```
累积 (seasonDecidPhenology): `chill_day += Cd` 至 crequ(−68); 达标后 `anti_chill_day += Ca`。
**budburst ⟺ anti_chill_day > 0  AND  anti_chill_day ≥ −chill_day**。
- 通用落叶 PFT 保留 Sykes 1996 (apple-only 切换, use_cesaraccio_chill=1)。

### §B.11.3 果实维持呼吸 (Group C, computeNPP.m)
```
R_m,fruit = r_m10 · N_fruit · g(T_air) · (1e−3 / 86400)        [kg C/(m²·s)]
g(T) = exp{ 308.56 · [ 1/56.02 − 1/(T_C + 46.02) ] }          (Sitch 2003 Eq 23 = Lloyd&Taylor 1994)
   T_C clamp [−40, 50];  r_m10 = 0.066 gC/(gN·d) @10°C (Sitch Table 3)
非果实 PFT: N_fruit = 0 → R_m,fruit = 0
```
- 源 CLM5 FruitTree CNMRespMod.F90:267-268 (O.Dombrowski): `grain_mr = grainn·br·tc` (perennial)。
- g(T) 与 R_m,root/stem 一致 (§B.3); 果实地上 → 用 g(T_air)。影响小 (果实 N 池小)。

### §B.11.4 采收 (harvested_flag, appleFruitHarvest.m + turnoverAndHarvest.m)
```
采收日 ⟺ harvest_signal (= phen_state.harvested_flag) ≥ 0.5  AND  ~harvested_this_year  AND  N_fruit > 0:
   C_harvest = harvest_fraction(0.90) · C_fruit           (turnoverCarbonMirror)
   N_harvest = harvest_fraction(0.90) · N_fruit           (turnoverAndHarvest)
   C_fruit −= C_harvest ;  N_fruit −= N_harvest           (留种 10%)
   harvested_this_year ← 1 (锁存)
重置 (re-arm): ~env.phen_state.harvested_flag (落叶/COLDEST_DAY 每年归零) → harvested_this_year ← 0
```
- **关键修复 2026-06-04**: 重置原挂 frPHU<0.01 (跨年不归零 → year-2+ 永不采收), 改挂物候 harvested_flag (每年 COLDEST_DAY 归零) [[project_compass_frphu_epic_lai_role]]。
- harvested_flag 由物候在 DS_fruit ≥ harvest_DS 时置 1 (seasonDecidPhenology)。

### §B.11.5 春季 onset-xfer 储备点火 (Group E, springOnsetXfer.m) — 解锁在线 leafc_sla LAI
苹果多年生春季展叶: budburst 时从 NSC 储备点火新组织 (光合启动前), 破在线 LAI 死亡螺旋。
fstor2tran=0.5, ndays_stor=50, α_k = 4 相位 budburst 分配比 (k∈{leaf,root,stem}, 归一到 ∑=1):
```
(1) budburst 一次性 mobilize (Cesaraccio 触发, onset_mobilized_flag 控一次):
      ΔC = fstor2tran · C_storage ;   ΔN = fstor2tran · N_storage
      C_k_xfer += α_k·ΔC ;  N_k_xfer += α_k·ΔN     (k = leaf,root,livestem)
      C_storage −= ΔC ;  N_storage −= ΔN

(2) onset 窗 (DOY ∈ [budburst, budburst+ndays_stor]) 每步 xfer→display:
      counter = max(dt, (budburst + ndays_stor − DOY)·86400)   [s]   (末步 counter=dt → f=1 排空)
      f = dt / counter
      ΔC_k = f · C_k_xfer ;  C_k += ΔC_k ;  C_k_xfer −= ΔC_k    (k 器官)
      ΔN_k = f · N_k_xfer ;  N_k += ΔN_k ;  N_k_xfer −= ΔN_k
```
- 源 CLM5 FruitTree CNPhenologyMod.F90:1946-1996 (perennial onset, 全器官 storage→xfer→display, C+N), :1890-1903 (xfer 按 onset_counter 排放)。
- **架构适配**: CLM5 每器官独立 storage; COMPASS 单 C_storage/N_storage → 按 α_k 分到 3 个器官 xfer 池。
- **C:N 一致 + 守恒**: C 和 N **全程同比例** (fstor2tran + α_k + f) 移动 → 移出 C:N = 储备 C:N, 不凭空造 C:N → 机器精度守恒; FlexibleCN 在分配步把器官 C:N 夹到目标带 [target, target+offset]。
- **xfer 池入守恒账** (COMPASS.m C_pool_now/N_pool_now + setInitialNitrogen 6 池 init 0): storage→xfer→display 全在池内流动。
- ndays_stor=50 (citrus idx35) ≠ 基础落叶树 ndays_on=30。触发用 Cesaraccio budburst_doy (§B.11.2)。

### §B.11.6 苹果尚未实现的 FruitTree 过程 (待办)
| 过程 | FruitTree 源 | COMPASS 状态 |
|---|---|---|
| 采后 storage_flag 储备回填 (Group D) | CNPhenologyMod:2035 (harvest→storage_flag=1) | 部分 (aleafstor), storage_flag 切换未做 |
| 修剪 prune (Group G) | prune_flag:1868 (offset 时剪枝去生物量) | 未做 |
| 有机肥 fert (Group H) | fert/manunitro:1967-1969 (budburst 后施 nfertdays_on 天) | 未做 (COMPASS 有施肥事件, 时序不同) |
| 移栽建植 transplant | :1744 (首年 transplant=5 gC 给叶/根/10%死木) | 未做 |
| 果园轮替 mxmat (3→20 龄→重植) | idpp/mxmat:1937 (达最大寿命→rotation) | 未做 |

#### 四个果树生理缺陷待办 (2026-06-09, 用户提出, 已核实现状)

**① 动态根系生长 (接入前必读 §B.11.8)** — 现状: 根分布用 Jackson 1996 固定 RAI 查表, 不随碳变, 不找水不
加深。要做动态根系需:
- (a) **拆细根/粗根 2 池** (不是 CLM5 的 3 池): C_root_fine=吸水+快周转+动态响应; C_root_coarse=储备+
  锚定+慢+几乎不吸水。动态响应的是细根, 粗根不参与快速找水。
- (b) **RAI 接细根碳** (打通两套根): RAI(层j) = C_root_fine(层j) × SRL × (2π·r_root) / 层底面积。
  即 碳→根长(×SRL)→根面积(×周长)。细根碳多→RAI 大→吸水强。替代现在的 Jackson 固定查表。
- (c) **分层 + 找水分配**: C_root_fine 按层存; 每步新根碳按层权重分 ∝ 层湿度 × 层已有根量 → 湿层多投根
  (树"找水"= 把碳投到湿土层)。对标 CLM5 dynamic root (Drewniak et al. 2019 GMD)。
- (d) **根深加深**: Z_root 随树龄/根碳渐近增大 (主线, 如 Z_root=Z_max·(1-exp(-k·根碳))); 或表层持续旱
  时加速下扎 (前沿)。加深 = RAI 分层数组增新深层格子。苹果根可达 2-3 m 找深层水。
- 闭环: 碳→长根→各层RAI→吸水→改土壤湿度→反馈下次找水分配。建议两步走: 先(a)+(b)让碳驱动吸水, 再(c)+(d)
  完整找水+加深。

**② LAI 峰形: 平顶 → 圆拱** — 现状: laimx=2.8 硬封顶 (peaklai cap, config_vegetation.txt:92) + ft_GDDlfmat
=1350 叶 taper→0, 造成 LAI 夏季**平台钉死 2.80** (yr1 DOY175-293 共118天平顶)。但长武实测成龄树是**带起伏的圆
拱** (DOY150-250 在 2.3-2.85 波动, 拟合非平顶), 见 data/苹果树LAI.xlsx 247obs + data/changwu_LAI_obs.png。
- 根因: laimx 硬 cap 一撞顶就锁死。CLM5 自己也是硬 cap (CNVegStructUpdateMod.F90:177 peaklai=1), 所以平顶
  其实**忠实 CLM5** — 改成圆拱 = **超越 CLM5**, 需引入叶龄周转 (老叶凋、新叶生, 盛夏自然起伏) 或软减速 (接近
  laimx 渐进降叶分配而非硬切)。属研究创新点非 bug 修复。是否做看是否要对齐实测 vs 对齐 CLM5。

**③ 大小年 (alternate bearing / biennial bearing)** — 现状: **完全没有** (grep 零命中)。苹果实测大小年 6 倍
(课题组靶: 2021大年49315 vs 2020小年8229 kg/ha, 见 memory apple_yield_benchmark)。机理: 大年果多→消耗储备
+抑制花芽分化→次年花少=小年→储备恢复→再大年。模型要做需: 花芽分化子模型 (当年果负荷反馈次年花量) + 储备耦合
(大年透支储备→小年低产)。这是果树**独有**、森林/作物没有的过程, 是 COMPASS 苹果对标实测逐年产量的关键缺口。
当前模型每年产量趋同 (无大小年), 是逐年验证的硬伤。

**④ 产量精度** — 现状: f_alloc 标定后**已对标实测差 5%** (提交 74e3d21: 18306→46868 vs 实测 49315 kg/ha,
原低 3 倍已根治)。剩余偏差可能来自: (i) 无大小年 → 逐年产量不对 (见③); (ii) 那 5% 系统偏差 (分配/果实膨大/采
收时序微调); (iii) 产量换算 C×148 (memory apple_yield_benchmark) 的系数。产量"不太对"更可能是**逐年模式不对
(无大小年)** 而非年均值, 优先级上③ > ④。

---

### §B.11.7 完整果树生理过程时序 (一年闭环) — 2026-06-09

苹果 (PFT 23) 一年的 C/N/物候完整链条, 由**冷激→积温**两套时钟串联驱动:

```
[休眠] 落叶后 dormant_flag=1
  │  累积冷激 chill_day (T < baset=4°C 的冷天)            ← Cesaraccio 2004
  │  冷激达 crequ = -68 °C·d → 解除休眠 ("睡够了")
  │  累积回暖 anti_chill_day → 抗冷激 ≥ 冷激缺口
  ▼
[发芽 budburst] onset_flag=1, DOY~80  ← 冷激链终点 = 积温链零点
  │  春萌储备点火: 各器官储备抽 fstor2tran=0.5 → xfer → display (springOnsetXfer)
  │  展叶: phengdd5ramp=200 °C·d 内 phen 0→1, ndays_on=30 d ramp
  │  分配 (fruitTreeAllocation 4 相位, GDD since budburst):
  │    根 arooti=0.7 → arootf=0.2 → arootf2=0.08 (三段台阶递减)
  │    叶 = (1-aroot)·fleafi=0.85·leafGDDtaper(g) (taper 1→0 over ft_GDDlfmat=1350)
  ▼ +GDD_bloom=369 °C·d
[开花 bloom] DOY~130
  ▼ +grnfill=400 °C·d (from budburst)  [注: 与开花仅差 ~31, 花谢即坐果]
[灌浆 fruit-fill] 碳大量切向果实, DS_fruit = GDD_since_bloom / TSum_fruit=2208 累积
  ▼ GDD=1350 (from budburst)
[冠层成熟] 叶 taper→0, LAI 封顶平台 (laimx=2.8 物理上限叠加)
  ▼ DS_fruit ≥ harvest_DS=0.95
[采收] DOY~285, 摘走 harvest_fraction=0.90 果实; 采后叶碳 aleafstor=0.3 → 储备
  ▼ GDD=hybgdd=2880 (from budburst)
[生理成熟]
  ▼ 日照 < crit_dayl_h=10.9 h, DOY~294
[落叶 senescence] offset_flag=1, ndays_off=45 d 线性落叶, 叶 N 按 nrelocfrac=0.7 回收→储备
  │  落叶若撞冷月重置 → leaf_flush_pending 强制排光残叶 (v63 修复)
  ▼ 月均温首跌破 5°C → 冷月重置 → dormant_flag=1
[休眠] (回到顶部, 循环)
```

**关键认识**: 发芽**无积温参数** — 苹果发芽由冷激 (crequ/baset) 触发, 防暖冬假醒; 发芽日是后续所有积温
(开花 369 / 灌浆 400 / 封顶 1350 / 成熟 2880) 的零点。TSum_fruit=2208 例外, 起点是**开花**不是发芽。

### §B.11.8 根系表示的三个局限 (2026-06-09, 接入动态根系前必读)

COMPASS 用**两套脱钩的根表示**: 碳循环一套, 吸水一套。三个局限:

1. **无显式粗根池**。CLM5 分 frootc / livecrootc / deadcrootc; COMPASS 压成单 `C_root_fine` (当年根, 名为"细
   根实为全部根的当年碳) + `C_root_storage` (粗+细根越冬储备)。粗根碳被并入这两池, **名不副实**。
   决策见 memory project_compass_coarse_root_decision: 多 PFT 长期跑也不需显式分 (LPJ-GUESS/FATES 都不分),
   保持单根池; 但命名误导是事实。

2. **碳根与吸水根脱钩**。吸水 (Plantrootconductivity) 用的 RAI (根面积指数) 和根长密度 Rl 来自 Jackson 1996
   全球根分布**固定查表** (calculateRootAreaIndex RAIlist + β), **完全不读 C_root_fine 碳池**。
   - ✅ 好处: 碳循环 C_root_fine 因 95% 根分配 (fcur) 都进它而偏高, **不影响吸水** (解耦保护, 吸水不会被高估)。
   - ⚠️ 局限: 碳长出多少根, 吸水的根不跟着变 — 两套根不一致。

3. **接入动态根系生长前必须先打通这两套**。要做根系动态 (分层/深度/根尖/碳驱动长根) 时, 现在的 Jackson 固定
   RAI 是死结 — 必须改成 RAI/根长由 C_root_fine 动态算 (carbon→root area 映射), 这是架构级改动。当前简化
   做不了真正的根系生长, 只能做"根分布固定、碳在固定分布上流动"。

**对储备衰减的关联**: 根分配 arooti=0.7 很大, 但其中 95% (fcur=0.95) 进 C_root_fine 当年细根 (吸水吸氮,
当年周转死掉变 litter, 是流量), 仅 5% (1-fcur) 进 C_root_storage 越冬 (存量)。春萌取走储备一半 (fstor2tran
=0.5) 却只存 5% → 衰减。修法 croot_stem 让粗根储备从茎流补 (见 §B.5 + memory clm5_xsmrpool_gap), 不动当
年根 (当年根是树活命必需, fcur 不能降)。

---

## §C 阶段路线图 → 见**文档开头**「📍 阶段路线图」(全文唯一权威, 含 12.5 C/N 耦合修复 F1-F5)。此处旧表已删, 防误读。

---

## §D 关键参数 (苹果 PFT 23, 含出处)

| 类 | 参数 | 值 | 出处 |
|---|---|---|---|
| 光合 | Vcmax25 (目标) | 105 μmol/m²/s | Wang 2009 苹果实测 |
| | Jmaxb1 | 0.244 [site-cal] | 反推, Ali 2016 范围内 |
| | TPU25 | 0.167·Vcmax25 | CLM5 1405/1517 |
| | Rd 光抑制 | ×0.67 白天 | CLM5 1562, Lloyd 2010 |
| | gs_g1 (Medlyn) | 8.2 | FruitTree Zanotelli 实测 |
| 呼吸 | r_m10 | 0.066 gC/gN/d | Sitch Table3 |
| | f_livewood | 0.22 [site-cal] | CLM5 CNMRespMod:265 |
| | R_g_fraction | 0.25 | Sitch p.169 |
| C:N | CN_L / root / stem / fruit | 19.7/32/60/139 | FruitTree Zanotelli 实测 (v10) |
| 分配 | f_alloc leaf/root/stem/fruit | 0.40/0.25/0.20/0.15 [site] | (待动态分配替代) |
| 物候 | GDD_bloom / TSum_fruit | 369/2208 °C·d | 长武 25yr 标定 |
| | k_chillb | 1150 | Sykes 1996 (苹果高需冷) |
| LAI | SLA0 | 0.0125 m²/gC | Wang 2009 (vs FruitTree 0.028, 待核 §E) |
| 异速 (LPJ, 自然树用) | k_latosa/k_allom1/2/3/wooddens/ltor_max | 6000/250/60/0.67/200/1.0 | LPJ global.ins group |

全 57 字段详见附录 §H.9 (旧 §9 PFT 表)。

---

## §E 待办 → 已并入**文档开头「📍 阶段路线图」**(细项 N/P/D + 通用 LSM 缺口 GL/PH 全挂对应阶段)

> 原 §E「待办合并清单」已删, 全部待办归到开头唯一权威路线图 (阶段表 + 「细项待办」+ 三层主线), 防分散打架。
> 多模型审计修复报告见 §12.34; 物候对标见 §12.35; 池结构对标见开头阶段 16 子项。

---

## §F 动态分配设计 (LPJ + FruitTree 三套分流, 阶段 13)

### §F.1 选型: LPJ-GUESS (vs FATES)
LPJ = 异速 (pipe model Eq2/4) + 环境功能平衡 (ltor 受 wscal/nscal 调, growth.cpp:1199)。
选 LPJ 理由: 全 PFT (树/草/作物) 全覆盖 + 功能平衡+异速双满足 + COMPASS R_m 已 Sitch 同源。FATES 偏木本 cohort, 草/作物无解。

### §F.2 三套分配分流 (用现有 PFT 字段)
```
crop_phenology==1 (15-22 作物) ─┐
PFT 23 苹果 (矮化果园)          ─┴→ [物候动态] FruitTree Eq4-7 (GDD: aleaf 衰减/arepr 增, 面积尺度)
CN_target_stem<=0 (12-14 草)   ──→ [LPJ grass] leaf:root 按 ltor (无茎)
其余 (1-11/24 自然树灌)         ──→ [LPJ tree] pipe model Eq2/4 + ltor (个体尺度, 需 tree_CAH/DBH 转换)
```

### §F.3 LPJ tree 核心方程 (growth.cpp:396-470)
```
(1) bminc = leaf_inc+root_inc+sap_inc                守恒
(2) leaf_area = k_latosa × sapwood_xs_area           pipe model (茎随叶)
(3) cmass_leaf = ltor × cmass_root                    功能平衡 (ltor 受 ψ 调)
(4) height = k_allom2 × diam^k_allom3                 树高异速
→ Eqn(13) 二分法求 leaf_inc
grass (无茎): leaf_inc = (bminc−leaf/ltor+root)/(1+1/ltor); root_inc = bminc−leaf_inc
```

### §F.4 漏洞决策 (5 个)
| # | 漏洞 | 决策 |
|---|---|---|
| 1 | 单位尺度 (LPJ 个体 vs COMPASS 面积) | 自然树用 tree_CAH/DBH 转换; **苹果改物候分配避开** |
| 2 | LAI 用 N_leaf×CN_L 非 C_leaf | 改直接用 C_leaf |
| 3 | height/diam 数值奇异 | 移植 LPJ clamp/cmass_debt |
| 4 | 苹果矮化≠通用乔木异速 | **苹果走 FruitTree 物候 (类作物), 不走 LPJ 乔木** |
| 5 | crop 分支 | DS 驱动 |

### §F.5 LPJ 全 PFT 异速参数 (group 继承, global.ins)
| PFT 类 | k_latosa | k_allom1 | k_allom2 | k_allom3 | wooddens | ltor_max |
|---|---|---|---|---|---|---|
| 阔叶树 | 6000 | 250 | 60 | 0.67 | 200 | 1.0 |
| 针叶树 | 5000 | 150 | 60 | 0.67 | 200 | 1.0 |
| 灌木 | — | 100 | 5 | 0.67 | 250 | 1.0 |
| 草 | — | — | — | — | — | 0.5 |

### §F.6 LAI 双路径并存 (LAI_engine 开关)
- epic_phen_mix (默认): EPIC 曲线, 校准不动, 分配用固定 f_alloc (不引入 Eq2)
- leafc_sla_phen (重构): LAI=C_leaf×SLA, 走完整 LPJ 异速 (含 Eq2 pipe model)
即 LAI_engine 同时决定 LAI 来源 + 是否走完整 LPJ 异速。

### §F.7 实施 (每步单独验证守恒)
步0 config alloc_scheme 开关 → 步1 苹果物候动态 (FruitTree, 无新状态, 先做最高ROI) → 步2 LPJ grass → 步3 LPJ tree (拆 sap/heart 双池+DBH 动态, 连锁 14 文件, 大工程) → 步4 ltor 功能平衡 (ψ 驱动) → 步5 苹果 fruit → 步6 crop。

*源码: LPJ-GUESS @ D:/OneDrive/LPJ-GUESS-v4.1.1/guess_4.1/modules/growth.cpp; FATES @ D:/OneDrive/fates; CLM5 @ D:/OneDrive/CTSM-release-clm5.0.37; FruitTree(Dabbaghi 2022)/CLM5 TechNote PDF @ data/references/allocation_allometry。*

### §F.8 C/N 主导范式调研 + C:N 弹性 (2026-05-29, 关键建模决策)

**起因**: v10 发现 C_fruit/N_fruit 实际 C:N=82-92, 但参数 CN_target_fruit=139 → 采收 C=N×139 反推超 C_fruit 池 → clamp → PlantC 守恒残差 -29.6 gC (KT6433 采收事件)。

**三模型 C/N 主导设计 (看实际代码)**:
| 模型 | 主导 | C:N | N 不足时 | 代码 |
|---|---|---|---|---|
| DSSAT/DAISY (农业) | **N 主导** (demand) | 可漂移 | N 限制生长 | — |
| LPJ-GUESS | C 主导 (光合) | cton_leaf [min,max] 浮动 | 限叶面积 | growth.cpp |
| **CLM5 FlexibleCN** | **C/N 协商** | actual_leafcn ∈ [leafcn, leafcn_max] **有界弹性** | downreg(下调GPP) 或 flexcn(C:N升高叶稀释) | NutrientCompetitionFlexibleCNMod.F90:38/433/439 |
| FruitTree (苹果) | **C 主导** (物候分配 C, N 跟随 = C/CN_target) | 严格守 CN_target | — | Dabbaghi 2022 p.4 "corresponding N" |

**CLM5 协商逻辑 (FlexibleCNMod.F90)**:
```
1. C 算 N 需求: GPP → plant_ndemand (gN/m²/s)
2. N 供给比 fpg = 根吸供给 / demand
3. fpg<1 (N 限制) 两选项 (downreg_opt 开关):
   downreg=true:  downreg = excess_cflux/gpp; psn×(1−downreg) 下调光合 → 保 C:N
   downreg=false: actual_leafcn 升高 (叶稀释), C 不丢 → flexcn
4. actual_leafcn ∈ [leafcn(目标), leafcn_max(上限)]  ← 有界, 不无限漂移
```

**COMPASS 现状问题**: N 主导 (demand) + C 按 f_alloc 独立 → **C:N 无约束自由漂移** (82-92)。比 CLM5 少了"C:N 有界弹性"。

**苹果靠什么 (非 LPJ — LPJ 无苹果/作物)**: 苹果走 **FruitTree** = C 主导 + N 跟随 (N=C/CN_target) → C:N 自动守 CN_target。

**结论 + 决策**:
1. **没有单一'C 主导 vs N 主导'答案** — CLM5 协商式 (C/N 双向 + C:N 有界弹性) 最成熟。
2. 苹果果实 C:N=139 是充分供 N 目标值 (果实糖多氮少, 合理); 82-92 = 果实分 N 偏多, 但**若长武 N 充足则真实** — 关键缺约束界, 不是值错。
3. **短期 (v11)**: 采收改 C_fruit×harvest_fraction (守恒, 尊重 C:N 弹性, 不强求 139)。**权宜但物理可接受** (摘整果, C/N 各按比例扣)。
4. **根治 (阶段 13a 苹果物候分配)**: 学 FruitTree **C 主导 + N 跟随 (N=C/CN_target)**, 或学 CLM5 **C:N 有界弹性 [min,max]**。让 C/N 绑定, 消除无约束漂移。保留 COMPASS N 限制特色 (CLM5 downreg/flexcn 思路)。

*源码: NutrientCompetitionFlexibleCNMod.F90:38(actual_leafcn)/251(leafcn_max)/433-490(downreg/flexcn); pftconMod.F90:135(leafcn)/225(fleafcn 灌浆期CN)/254(fun_cn_flex); LPJ growth.cpp; FruitTree Dabbaghi 2022 p.4。*

### §F.9 C-follows-N 架构重构 (2026-05-29 决策, 根治 C:N 漂移)

**起因 (细查 v10)**: C_fruit/N_fruit 实际 C:N=82-92 ≠ CN_target=139。逐步追根:
- 单步分配 A_C_fruit/A_N_fruit = **139** (正确, 公式对)
- **但初值 C:N = 旧值** (fruit C 11.56/N 0.385 = 30; leaf 25 vs 19.7; root 35 vs 32; stem 100 vs 60)
- 初值 30 + 增量 139 混合 → 池 C:N 漂到 82-92
- **根因**: 改 CN_target (v10) 时忘同步初值; 更深是 **C_init 和 N_init 是两个独立 config 硬值, 物理上应绑定 (N=C/CN)**, 各填一个必然脱节。

**更深的架构冗余 (用户洞察)**: 当前 C 和 N **各分各的**:
- N 侧: D_N = NPP·f_alloc/CN_target; A_N = D_N·f_N_growth (N 限制打折)
- C 侧: A_C = NPP·f_alloc (独立, 不打折, 没用 N 结果)
- → C 算了 N 需求 (D_N=C/CN), 但 C 分配又自己来一遍, **不自洽 → C:N 漂移**。

**已实现架构 — C 主导 + C:N 弹性 (v12, 2026-05-29 完成 + 单元验证)**:

> ⚠️ 纠错记录: 本节初版 (同日早些) 写成 **N 主导** (A_C=A_N×CN_target, C:N 恒=target)。
> 用户两次纠正"不是 N 主导是 C 主导"。读 CLM5 源码确认 CLM5 是 **C 主导**
> (availc 先定 → nlc=calloc/c_allometry → N=C/leafcn → C:N 经 supply/demand+nscalar 弹性,
> FlexibleCN:508/520/551/1660)。下为纠正后的最终实现:

```
1. A_C_target = NPP · f_alloc_organ                C 先定 (computeOrganDemand)
2. D_N_organ  = A_C_target / CN_target_organ        N 需求 (C 提需求)
3. f_N_growth = clamp(N供给 / Σ(D_N·dt), 0, 1)       N 限制比 (allocatePlantN)
4. A_N_organ  = D_N · f_N_growth · dt               实际 N 分配 (N 守恒不变)
5. actual_CN  = min(CN_target / f_N_growth, CN_target+offset)   ★C:N 弹性, clamp 到 CN_max★
6. A_C_organ  = A_N_organ · actual_CN  ★关键★        C 跟随 (N 足→全额; 缺→C:N 上浮 C 照长; 撞顶→砍 C)
7. 砍掉的 C (A_C_target − A_C) 经 dC_labile=GPP−R_a−ΣA_C 自动留 C_labile
```

**好处**:
| | |
|---|---|
| C:N 有界弹性 | actual_CN ∈ [CN_target, CN_target+10], 根治旧 v10 无界漂移 (82-92) |
| C 主导护 LAI | N 缺时 C 照长 (C:N 上浮), 不减叶量 → 不与 EPIC LAI 校准争控制权 |
| N 限制经光合反馈 | 叶 C:N↑ → 单位叶含氮↓ → LUNA 下调 Vcmax (闭环2), 物理真实 |
| N 守恒不破 | allocatePlantN Step1/4 不动, 仅加 actual_CN 输出 (单元验证 err=0) |
| 对标 CLM5/FruitTree | CLM5 nscalar (FlexibleCN:1660-1667); FruitTree C 主导 N 跟随 |

**已改动 (全部完成)**:
1. **VegetationData.m + config_vegetation.txt**: 加 cn_coupling_mode (c_follows_n/parallel) + cn_flex_offset=10
2. **allocatePlantN.m**: 加第5参 pftParams; Step3b 算 actual_CN (N 不足时上浮 clamp CN_max), 输出 .actual_CN_*
3. **computeCarbonAllocation.m**: c_follows_n 分支 A_C = A_N × actual_CN (A_N_g 带 .CN_* 字段)
4. **runPlantCNStep.m + runPlantNDemand.m**: allocatePlantN 传 pftP_demand; A_N_g 带 actual_CN; mode 开关透传
5. **初值**: config 删 5 organ C_*_init 硬值 → setInitialNitrogen fallback C=N×CN_target (自洽)

**单元验证 (2026-05-29 通过)**: 三场景 — N 足 f_N=1 actual_CN=19.7 A_C 全额; N 微缺 f_N=0.85
actual_CN=23.2 C 照长; N 严重缺 f_N=0.30 actual_CN=29.7 撞顶 A_C 0.016<目标 0.036 砍 C; N 守恒 err=0。

**保留兼容**: config cn_coupling_mode = c_follows_n(默认, C 主导+弹性) / parallel(旧 A_C=NPP·f_alloc 对比)。

**待**: 全年验证 C:N 落在 [target, CN_max] + Plant C/N 守恒机器精度 (跑中)。后续阶段 13 动态分配 (f_alloc 动态化, C 主导框架不变)。

*依据: CLM5 FlexibleCN nscalar (NutrientCompetitionFlexibleCNMod.F90:508/520/551/1660-1667);
FruitTree C 主导 (Dabbaghi 2022 p.4 + Table C1); 用户架构决策 2026-05-29 (C 主导 + CLM5 ±10 弹性)。*

---

## §G 参考文献

**物理基础**: Šimůnek 2008 (HYDRUS); van Genuchten 1980 (SSSAJ 44:892)。
**根吸**: Nacry 2013 (Plant Soil 370:1); Wells & Eissenstat 2002; DAISY 7.1.0。
**C/N 分配+物候**: Sitch 2003 (GCB 9:161); Mao 2018 (苹果果实 N); Aerts 1996 (Oikos 67:469, retranslocation); Tagliavini 2003 (苹果茎 N); Wang 2009 (苹果叶 N/Vcmax 长武); Sykes 1996 (chilling-GDD); Dabbaghi/Dombrowski 2022 (GMD 15:5167, CLM5-FruitTree)。
**光合+LUNA**: Farquhar 1980 (Planta 149:78); Medlyn 2011 (GCB 17:2134); Kattge & Knorr 2007 (PCE 30:1176); Bernacchi 2001 (PCE); Ali 2016 (GMD 9:587, LUNA); Lloyd & Taylor 1994 + Lloyd 2010 + Metcalfe 2012 (Rd 光抑制); Collatz 1991。
**异速**: LPJ-GUESS (Smith 2001/Sitch 2003); FATES (Fisher 2015); Shinozaki 1964 (pipe model); Ghimire 2016 (JAMES, CLM5 分配); Poorter 2012 (New Phytol, 全球分配 meta)。
**NSC/呼吸**: Sala 2012 (Tree Physiol 32:764); Hartmann 2016 (New Phytol 211:386); Ryan 1991。
本地源码/PDF: D:/OneDrive/{CTSM-release-clm5.0.37, DAISY, LPJ-GUESS-v4.1.1, fates}; data/references/。

---

# 附录 H — 历史实施档案 (原 §0-§12.33, 仅查阅)

> 以下是重写前的完整历史文档 (M1~M2-B.8 逐步实施 + 踩坑 + v3-v10 验证细节)。旧 M 编号映射见 §C 路线表。
> 当前状态/方程/待办**以上方正文 §A-§G 为准**, 本附录仅供追溯实施过程。

## 0. 阶段总览 → **旧 M1~M3 路线图已删** (2026-05-20 版, 严重过时)。当前唯一权威路线图见**文档开头「📍 阶段路线图」**(阶段 1-18 + 12.5 C/N 耦合修复)。下方附录 §1+ 仅留各阶段实施细节/踩坑追溯。

---

## 1. SPAC 全耦合闭环路径

### 1.1 完整数据流

```text
                              ┌─────────────────┐
                              │ Meteo (forcing) │
                              │  Ta, P, Rn, RH  │
                              └────────┬────────┘
                                       │
              ┌────────────────────────┼────────────────────────┐
              ▼                        ▼                        ▼
    ┌──────────────┐         ┌──────────────┐         ┌──────────────┐
    │ Energy balance│         │ Richards eq │         │ N chemistry  │
    │ ebal_clm45    │         │ θ, ψ, K(h)  │         │ 硝化/反硝化  │
    │ → Rn,H,LE,G,T │         │ → soil θ    │         │ → SinkS, S_min│
    └───────┬───────┘         └──────┬──────┘         └──────┬───────┘
            │                        │                        │
            └─────────┬──────────────┴────────────┬───────────┘
                      ▼                           ▼
            ┌──────────────────┐      ┌─────────────────────┐
            │ Plant hydraulic  │      │ Root uptake          │
            │ (Plant_hydraulics)│     │ DAISY+Nacry+demand   │
            │ → ψ_leaf, sfactor │     │ (M2-B.2)             │
            └─────────┬────────┘      │ → U_N, SinkS         │
                      │               └─────────┬───────────┘
                      ▼                         │
            ┌──────────────────┐                │
            │ LUNA (每日)       │                │
            │ N→Vcmax (M2-B.3)  │                │
            │ → leafbio.Vcmo    │                │
            └─────────┬────────┘                │
                      ▼                         │
            ┌──────────────────────────────────┴─┐
            │ biochem.m (Farquhar + Ball-Berry)  │
            │  Vcmax(T,ψ) → A_canopy → GPP        │
            └─────────────┬──────────────────────┘
                          ▼
            ┌────────────────────────────────────┐
            │ Plant N module (M2-B.1b)            │
            │  computeNPP / computeOrganDemand /  │
            │  allocatePlantN                     │
            │  → 5 池: labile/leaf/root/stem/fruit│
            └─────────────┬──────────────────────┘
                          ▼                          ⏳ 待 M2-B.4
            ┌─────────────────────────────────────┐
            │ Turnover + Harvest + N_storage     │
            │ → Litter (M3 SOM 入口)              │
            └─────────────────────────────────────┘
```

### 1.2 三层时间尺度

```text
KT (1h)  forcing 读 / Richards / ebal / biochem / SoilN ADE / root_uptake / Plant N alloc
day      LUNA 优化 Vcmax25 (KT % 24 == 0)
10-day   Climate10dBuffer 滑动均 (Ta/Tleaf/PAR/RH/CO2/rb/GPP)
year     M2-B.4 harvest event (frPHU > 0.90)
```

---

## 2. 模块组织 (current state)

按 Skill §4.1 ClimaLand standalone/integrated 风格:

```text
src/
├── +core/                            数据结构 + 共享 state
│   ├── Climate10dBuffer.m            10 天滑动均 (M2-B.3 新加)
│   ├── ModelSettings.m
│   ├── SoilProperties.m
│   └── ...
├── +soil/+nitrogen/                  土壤 N (单体)
│   ├── +chemistry/                   硝化/反硝化/矿化
│   ├── +organic/                     SOM (M3 待)
│   ├── +transport/                   ADE
│   └── +root_uptake/                 DAISY+Nacry+demand
│       ├── solveRootSurface.m       二分法 (M2-B.2 加 demand mode)
│       ├── nitrogenUptake.m         demand 拆分 + Nacry cap
│       ├── daisySingleRootFlux.m
│       └── nacryMembraneFlux.m
├── +plant/
│   ├── biochemical.m                 Farquhar (M2-B.3 加 use_luna hook)
│   ├── +PlantNitrogen/               5 organ pool (M2-B.1b)
│   │   ├── computeNPP.m             LPJ Sitch 2003
│   │   ├── computeOrganDemand.m     DOY 150-290 fruit gate
│   │   └── allocatePlantN.m         N-limited alloc
│   ├── +PlantHydraulic/              ψ_leaf + sfactor
│   ├── +root/                        applyRootUptake
│   ├── +integrated/                  跨子模块耦合
│   │   ├── computeNEE.m              NEE = R_h + R_a − GPP (M2-B.1b)
│   │   └── +LUNA/                    CLM5 LUNA v1.0 (M2-B.3)
│   │       ├── computeLuna.m         主入口
│   │       ├── nitrogenAllocation.m  爬山优化
│   │       ├── nitrogenInvestments.m N→Vcmax/Jmax
│   │       ├── photosynthesisLuna.m  Farquhar+BB ci 迭代
│   │       ├── lunaParams.m          14 常数
│   │       ├── quadraticRoots.m
│   │       ├── +nue/                 NUEref + NUEcurrent
│   │       └── +tempresp/            Kattge + Bernacchi
│   └── getPFTParameterSet.m         24 PFT × 57 字段
├── +diagnostics/                    输出注册
└── COMPASS.m                        主循环调度
```

---

## 3. M1 实施完成档案 — DAISY 单根 + Nacry 根膜联立 (2026-05-18)

**目标**: 替代 v33.2 passive+active 二分支根吸, 用单一物理机制 (单根稳态柱 ADE + 根膜 M-M).

### 3.1 物理框架

每层根, 同时求解:
- **DAISY 单根柱 ADE** (Kløffel 2026, DAISY 7.1.0 Appendix 7.2): 单根周围土壤 N 浓度场稳态分布, 给出"土壤端供给" $I_{DAISY}(C_{bulk}, C_{root}, q_w, D, r_r, L_{rld})$
- **Nacry 根膜 M-M** (Nacry et al. 2013 Plant Soil 370:1, Couturier 2007): 根表 N → 细胞质, $I_{Nacry}(C_{root}, V_{max}, K_m, C_{min}, SRL)$
- **耦合方程**: $I_{DAISY}(C_{root}) = I_{Nacry}(C_{root})$, 用二分法求 $C_{root} \in [C_{min}, C_{bulk}]$

### 3.2 关键文件

```
src/+soil/+nitrogen/+root_uptake/
  ├── nitrogenUptake.m           主调度 (M2-B.2 加 demand 模式)
  ├── solveRootSurface.m         二分法 (M2-B.2 加 demand mode)
  ├── daisySingleRootFlux.m      DAISY 单根稳态解析
  └── nacryMembraneFlux.m        Nacry M-M
```

### 3.3 已知边界状态 (M1)

- `cbulk_below_cmin`: 土壤 N 太低 → 不吸
- `membrane_limited`: 土壤供应充足 → 根膜瓶颈
- `soil_limited`: 根膜想吸但土壤跟不上
- `converged` / `failed`: 二分收敛或不收敛

### 3.4 M1 验证 (2026-05-18 全年 smoke)

- ✓ NH4/NO3 守恒 0.0% (机器精度)
- ✓ 全年总吸合理量级 (但偏高 ~800 kg/ha, 因无植物 demand 反馈)
- ⚠ 触发 M2-A → M2-B 系列改进

---

## 4. M2-B.1b 实施完成档案 — PFT-unified + 5 organ N pool + LPJ NPP + NEE (2026-05-19)

### 4.1 目标

- 拆 `STEMMUS_N.m` 单体, 整体 PFT-化 (24 PFT)
- 加 LPJ Sitch 2003 NPP / R_a 完整分解
- 5 organ N pool 动态 (labile + leaf + root + stem + fruit)
- NEE 完整化 (移出 ebal, R_h + R_a − GPP)

### 4.2 关键文件

```
src/+plant/+PlantNitrogen/
  ├── computeNPP.m              R_m_leaf + R_m_root + R_m_stem + R_g = R_a
  ├── computeOrganDemand.m      f_alloc × NPP / CN_target → D_N,i
  └── allocatePlantN.m          U_N + D_N + N_labile_prev → 4 organ + labile
src/+plant/+integrated/computeNEE.m   NEE = R_h + R_a − GPP
src/+plant/getPFTParameterSet.m       24 PFT × 53 字段 → 后扩到 57 (M2-B.3)
```

### 4.3 5 organ pool

```text
prognostic (per m² ground):
  N_labile     [g N/m²]  临时仓库, 根吸进, 4 organ 出
  N_leaf       [g N/m²]  LUNA 反馈用 (经 lnc, M2-B.3 改 CLM5 默认)
  N_root_fine  [g N/m²]
  N_stem_wood  [g N/m²]  annual PFT = 0
  N_fruit      [g N/m²]  DOY 150-290 累积 (Apple)
```

### 4.4 NPP / R_a 分解 (LPJ Sitch 2003 Eq. 8-14)

$$
R_m^{leaf} = R_{d,leaf}^{Farquhar} \quad (\text{已在 biochem 计算, M2-B.1a 抽出})
$$

$$
R_m^{root} = r_{m25}^{root} \cdot Q_{10}^{(T_{soil}-25)/10} \cdot N_{root,fine} \cdot 10^{-3}
$$

$$
R_m^{stem} = r_{m25}^{stem} \cdot Q_{10}^{(T_{air}-25)/10} \cdot N_{stem,wood} \cdot 10^{-3}
$$

$$
R_g = R_{g,frac} \cdot \max(GPP - R_m, 0)
$$

$$
NPP = \max(GPP - R_m - R_g, 0)
$$

$$
R_a = R_m^{leaf} + R_m^{root} + R_m^{stem} + R_g
$$

### 4.5 NEE 完整化

$$
NEE = R_h + R_a - GPP
$$

其中 $R_h$ 暂用 `fluxes.Resp` (SCOPE 简化土壤呼吸 — 等 M3 SOM 替换).

### 4.6 全年验证 (output_v40_M2B1b_full)

- ✓ 4 organ pool ΔN ≡ Σ A_i (机器精度)
- ✓ NH4/NO3 守恒 0.0%
- ⚠ NEE 反符号 (R_h SCOPE 简化偏大 1164 g/yr)

---

## 5. M2-B.2 实施完成档案 — Demand-driven 根吸 (2026-05-20)

### 5.1 目标

替代 M2-A 的 sink_scale 临时补丁. 用 **DAISY fixed_sink 反演**: 根据整株 demand 反算 $C_{root}$, 让根吸自然受 N_labile 反馈.

### 5.2 关键改动

**solveRootSurface.m 加 mode**:
```matlab
% M1: 解 I_DAISY = I_Nacry (耦合平衡)
% M2-B.2: 解 I_DAISY = demand_per_root (反演)
```

**nitrogenUptake.m demand 拆分**:
```matlab
demand_NH4 = D_N_total × Vmax_NH4 / (Vmax_NH4 + Vmax_NO3)
demand_NO3 = D_N_total × Vmax_NO3 / (Vmax_NH4 + Vmax_NO3)
demand_per_root_layer = demand × RLD_weight
```

**demand 公式 (修正后)**:
```
labile_deficit = max(0, N_target − N_labile_prev)
demand_total   = labile_deficit / dt
```

不加 D_N_growth (避免 N_labile 涨到 40+). alloc 内部从 N_labile 取 growth.

### 5.3 调试历程

**第一次全年 (v1)**: demand 公式错误 = D_N + (target − labile)/dt → N_labile 失控涨到 **40.55** g/m²

**修正后 (v2)**: demand = max(0, (target − labile))/dt → N_labile 稳定 **1 → 1.5** g/m²

### 5.4 全年验证 (output_v40_M2B2_fix_full)

| 指标 | 值 |
|---|---|
| N_labile 起 → 终 | 1.0 → 1.51 ✓ |
| Σ A_i ≡ ΔN_organ | 21.19 / 21.19 (残差 0) |
| Annual SinkS NH4 | 13.65 kg N/ha |
| Annual SinkS NO3 | 8.05 kg N/ha |
| 总根吸 | **21.7 kg N/ha** (M1 时代 800+) |
| supply_limit_factor | mean=3.23 (土壤供给>demand 3×) |

**注**: 总根吸 21.7 偏低 (苹果实测 200-400) — 因为 M2-B.4 未上, N_leaf 不下降, demand 不补.

---

## 6. M2-B.3 实施完成档案 — LUNA N→Vcmax 反馈 (2026-05-20)

### 6.1 目标

兑现 §11 钦定 "Ali 2016 LUNA 为 Vcmax/Jmax 主要模型出处" 承诺. 完整移植 CLM5 LunaMod.F90 (1352 行 Fortran).

### 6.2 算法 1:1 移植 (Ali 2016 + CLM5 Tech Note Ch.9)

**5 池机理分配** (`F_{NCa} = N_{lc} + N_{et} + N_{cb} + N_{resp} + N_{store}`):

| 池 | 决定 | 默认占比 (Apple) |
|---|---|---|
| $N_{lc}$ | 光捕获 (chl) → θ | ~0.20 |
| $N_{et}$ | 电子传递 → Jmax | ~0.04 |
| $N_{cb}$ | Rubisco → Vcmax | ~0.09 |
| $N_{resp}$ | 呼吸酶 → Rd | ~0.01 |
| $N_{store}$ | 储备 | ~0.66 |

**每日爬山优化** (LunaMod.F90:757-942):
$$
\max_{N_{lc}} \big[\text{PSN}(N_{lc}) - C_{cost}(N_{lc})\big]
$$

**容量公式** (LunaMod.F90:372-373):
$$
V_{cmax,25}^{opt} = N_{cb} \cdot F_{c,25}, \quad F_{c,25} = 294.2 \text{ (Rogers 2014)}
$$
$$
J_{max,25}^{opt} = N_{et} \cdot F_{j,25}, \quad F_{j,25} = 1257.0 \text{ (Xu 2012)}
$$

**每日变化率限速** (Suzuki 2001):
$$
\Delta V_{cmax} \le V_{cmax,old} \cdot Q_{10}^{enz, 0.1(T_l-25)} \cdot 0.1
$$

### 6.3 关键文件 (~750 行 MATLAB)

```
src/+plant/+integrated/+LUNA/
  ├── lunaParams.m              14 LUNA 常数 + Jmaxb1 override
  ├── computeLuna.m             主入口 + daily-rate + decay branch
  ├── nitrogenAllocation.m      爬山优化 (LunaMod.F90:757-942 1:1)
  ├── nitrogenInvestments.m     N → Vcmax/Jmax (LunaMod.F90:948-1028)
  ├── photosynthesisLuna.m      内嵌 Farquhar+BB ci 迭代
  ├── quadraticRoots.m
  ├── +nue/{nueRef, nueCurrent}.m
  └── +tempresp/{vcmxTKattge, jmxTKattge, respTBernacchi}.m

src/+core/Climate10dBuffer.m    10 天滑动均 (8 ring buffer)
```

### 6.4 调度 (COMPASS.m)

```text
KT 循环:
  forcing read → meteo / rad / canopy
  luna_clim.push(KT, forcing+thermal)        # 每 KT push 10 天 buffer
  if KT % 24 == 0 && luna_clim.isWarm():
      Jmaxb1 = vegetation.getPFTValue(pftType, 'Jmaxb1')
      lunaIn = luna_clim.assembleLunaInput(canopy)
      lunaOut = computeLuna(lunaIn, luna_state, lunaParams(Jmaxb1))
      luna_state ← lunaOut  (carry over)
  leafbio.Vcmo = luna_state.Vcmax25_z(1)     # 覆盖 PFT 表 105
  ebal_clm45(...)                             # ebal 完全不改, 自动用上
  biochem.m: Vcmax(T) = Vcmo × f_T × sfactor
```

### 6.5 关键设计决策

**CLM5 lnc_opt=.false. 路径** (PhotosynthesisMod.F90:1321):
$$
lnc = 1 / (SLA \times CN_L)
$$

lnc 是 **PFT 常数**, **不依赖 prognostic N_leaf 池**. 解决 mid-year 启动问题 (cold-start N_leaf=0 / 2 偏小):
- Apple: lnc = 1/(0.0125 × 25) = **3.20 g N/m² leaf** (Wang 2009 实测范围)
- FNCa = 3.20 − 0.32 (SNCa) = 2.88

**N_leaf prognostic 池** 仍跟 root uptake + NEE 耦合, 只是不直接喂 LUNA.

### 6.6 调试历程 (5 bug)

| Bug | 症状 | 修复 |
|---|---|---|
| #1 N_leaf 时序 | LUNA 看 N_leaf=0 | `N_leaf(KT-1)` (KT 时序在 alloc 之前) |
| #2 lnc 用 prognostic | mid-year 启动偏低 | **CLM5 默认 `1/(SLA×CN)`** |
| #3 PAR 单位 | LUNA 看 PAR≈0 | `rad.PAR [mol/m²/s] × 1e6 → μmol → W/m²` |
| #4 Apple SLA0 | lnc=1.33 偏低 | `0.030 → 0.0125` (Wang 2009 / CLM5 BDT) |
| #5 leafbio.Vcmo 未覆盖 | GPP 跟 baseline 100% 一致 | 直接 `leafbio.Vcmo = LUNA.Vcmax25` (0 改 ebal) |

### 6.7 1:1 算法审计

跟 CLM5 LunaMod.F90 严格 1:1 比对, 14 项一致 + 1 个 cold-start bug (PNlc=0.3→0). 设计选择:
- `lunaIn.tleaf10 = tleafd10` (CLM5 line 300 自带 bug, 1:1 复制)
- `Jmaxb1` per-PFT (CLM5 namelist scalar, 我们扩展)

### 6.8 PFT 表关键参数 (Apple, M2-B.3 修正)

| 字段 | 原值 | 改为 | 出处 |
|---|---|---|---|
| `SLA0` | 0.030 | **0.0125** | Wang 2009 / CLM5 BDT_temperate (Lawrence 2019 Table 8.1) |
| `SLAm` | 0.004 | **0.00133** | CLM5 BDT_temperate dsladlai (CLM4.5 Tech Note) |
| `Jmaxb1` | 0.1745 | **0.1745 (保留)** | [site] Apple 站点反推 (Wang 2009 Vcmax25=105); 其他 23 PFT 改 0.0936 (CLM5 默认) |

### 6.9 全年验证 (output_v40_M2B3_v3_full, 2021)

| 指标 | 值 |
|---|---|
| **Vcmax25_LUNA** | 50-93 (mean 75.4) μmol/m²/s |
| **Jmax25_LUNA** | 85-185 (mean 149) |
| FNCa | 2.88 (常数, CLM5 默认 lnc) |
| PNcb (Rubisco N) | 0.07-0.11 (mean 0.091) |
| PNlc (光捕获 N) | 0.18-0.24 (mean 0.20) |
| **GPP 年累计** | **1059.7** g C/m²/yr |
| **vs M2-B.2 baseline (PFT=105)** | **−15.07%** |
| NEE | +691.7 (反符号, 待 M3 SOM) |

**季节动态** (Vcmax25 月均):
```
Jan 61.5  Feb 89.0  Mar 83.9  Apr 76.0  May 78.9
Jun 74.6  Jul 72.0  Aug 70.8  Sep 68.6  Oct 68.2  Nov 81.7  Dec 79.9
```

冬末 (Feb 89) / 秋落叶 (Nov 81) 高: LAI 低 + PAR 强 → 单位叶 N 高 → LUNA 投资 Rubisco 多.
夏 (Aug 71) 低: 满冠 + 光合饱和 → 优化保守.

跟 Ali 2016 Fig 5 一致.

### 6.10 守恒 (机器精度)

| | 残差 |
|---|---|
| 4 organ pool ΔN ≡ Σ A_i | 0 (完美) |
| NH4_cum_residual | -6.3e-6 |
| NO3_cum_residual | +1.7e-10 |
| water_balance_residual | +0.006 mm/yr |

### 6.11 已知限制

| 限制 | 后续 |
|---|---|
| LUNA nrad=1 (单层, fVh 接力分层) | M3 LUNA 分层 (nrad=SCOPE 60) |
| N_leaf 只进不出 (M2-B.4 待) | M2-B.4 turnover |
| NEE 反符号 | M3 SOM |
| `lnc_opt=.false.` 失去 N_leaf 直接反馈 (CLM5 默认行为) | 未来 cfg.lnc_opt=true |

### 6.12 物理图景

```text
    LUNA (每日)                      ebal (每 KT)                biochem
                                                                
    Climate 10d:                                                
      Ta, T_leaf_d/n,                                           
      PAR, RH, CO2, rb                                          
            │                                                   
            ▼                                                   
    lnc = 1/(SLA×CN) = 3.20                                     
    爬山优化 5 池                                                
    Vcmax25_LUNA ────────────▶ leafbio.Vcmo ──┐                
                                              │                
                                              ▼                
                                  fVh × leafbio.Vcmo            
                                              │                
                                              ▼ biochem_in.Vcmo
                                       ┌──────────────────────┐
                                       │ Vcmax(T,ψ) =          │
                                       │   Vcmo × f_T(T_leaf)  │
                                       │       × sfactor(ψ)    │
                                       │ → GPP                 │
                                       └──────────────────────┘
```

3 层串联: **N 控容量 (LUNA, 日)** → **冠层 N decay (fVh, KT)** → **温度+水力 (biochem, KT)**

---

## 7. M2-B.4 实施计划 — Turnover + Retransl + N_storage + Apple Harvest (待实施)

### 7.1 目标

解决 M2-B.3 两个物理不真实问题:
1. N_leaf 全年只涨不降 (2 → 12.5 g/m²) — 没掉叶/没 retransl
2. N_fruit 永远积累 (0 → 2.04) — Apple 有果实但**没收获**

### 7.2 4 个 N 转移事件 (走现有 frPHU, 不引入 DOY hardcode)

**(a) Apple 果实收获** (一次性, frPHU > 0.90):
```
N_removed = N_fruit × harvest_fraction
cum_N_harvest += N_removed
```

**(b) 落叶 retranslocation** (每日, frPHU > DLAI):
```
N_leaf_senesce_rate = 1 / (leaf_long_days × (1 - DLAI))
d_N_leaf = N_leaf × rate × dt
N_storage += d_N_leaf × retransl_fraction
cum_N_litter += d_N_leaf × (1 - retransl_fraction)
```

**(c) 根周转** (每日, 全年):
```
d_N_root_litter = N_root_fine × root_turnover_rate × dt
cum_N_root_litter += d_N_root_litter
```

**(d) 春萌 N_storage 释放** (每日, 0 < frPHU < spring_release_endphu):
```
d_N_release = N_storage / spring_release_days × dt
N_labile += d_N_release
```

### 7.3 PFT 表新加 5 字段 (group [M])

| 字段 | Apple | 出处 | 校准等级 |
|---|---|---|---|
| `harvest_threshold` | 0.90 frPHU | EPIC harvest index (Williams 1989) | [tune] |
| `harvest_fraction` | 0.90 | **[site] 校准** (商业 80-95%) | ❌ |
| `root_turnover_rate` | 1/365 = 0.00274/day | Wells & Eissenstat 2002 (~年率) | ⚠️ [site] |
| `spring_release_endphu` | 0.20 | **[site] 校准** | ❌ |
| `spring_release_days` | 60 | **[site] 校准** | ❌ |

### 7.4 实施步骤 (5 步, 总 3 天)

```
S1 (0.5天) PFT 表 group [M] × 24 PFT + 头部注释
S2 (1天)   +plant/+integrated/turnoverAndHarvest.m + +turnover/ 4 子函数
S3 (0.5天) COMPASS.m 调度 hook + cum_* 累积
S4 (0.5天) 诊断 9 字段 + reg.write + CSV
S5 (0.5天) 单元测试 + 7天 + 全年 + §12.25 存档
```

### 7.5 验证标准

| 指标 | 期望 |
|---|---|
| N_leaf 季节 | 春萌 2 → 夏峰 14 → 秋冬 retransl 后 ~2-3 |
| N_fruit | DOY 150-290 升, frPHU=0.9 时一次性扣 ~10% |
| N_storage | 起 0 → 秋峰 ~7 → 春萌释放后 ~0-2 |
| 整株 N 守恒 | 机器精度 (Σpool + cum_harvest + cum_litter ≡ ΣU_N × dt + 初值) |
| 总根吸 | M2-B.2 21.7 → **200-400 kg N/ha** (Wang 2009 苹果实测) |

---

## 7b. M2-B.8 实施计划 — Plant C Pool (6 池, 跟 N 池 1:1) (待 M2-B.4 验证后)

### 7b.1 战略定位

**目标**: 完成 **W (water) + E (energy) + C (carbon) + N (nitrogen) 完整耦合** SPAC 真实循环.

**架构选择**: **DAISY 简化路线** (5-6 池, 单态) + **CLM5 关键算法** (PFT 表 / LUNA / Sitch alloc) 杂交.

```text
设计哲学:
  - 苹果园站点应用 (不是全球陆面) → 不需要 CLM5 三态 + 几百年 spin-up
  - 跟现有 N 池 1:1 平行 (DAISY DM 池路线)
  - 留 cfg.l* 开关给陆面扩展用 (M3+ 时可启 CLM5 三态)
```

参考文献:
- **CLM5** `CNVegCarbonStateType.F90`: 三态池 (storage/xfer/主, 19+ 池) — 过度复杂, 不采纳
- **DAISY** `production.C`: 5 个 DM 池 (单态), `WLeaf/WStem/WSOrg/WRoot/WDead` — **核心参考**
- **NoahMP** `module_Noahlsm.F.niu`: 4+1 mass 池, `LFMASS/RTMASS/STMASS/WOOD/GRAIN` — **辅助参考**

### 7b.2 6 个 C 池 (跟你 6 个 N 池 1:1)

```text
# Prognostic state (per m² ground)
C_labile       g C/m²  ← 跟 N_labile 平行 (NSC, Non-Structural Carbohydrates 临时仓库, CLM5 cpool 风格)
C_leaf         g C/m²  ← 跟 N_leaf
C_root_fine    g C/m²  ← 跟 N_root_fine
C_stem_wood    g C/m²  ← 跟 N_stem_wood (multi-year, 不周转)
C_fruit        g C/m²  ← 跟 N_fruit (DOY 150-290 累积 + frPHU>0.90 收获)
C_storage      g C/m²  ← 跟 N_storage (NSC + 木质 N 互转, 跨年缓冲)
```

### 7b.3 初值 (init from N × CN_target)

```text
C_labile_init    = N_labile_init   × CN_L          = 1.0 × 25 = 25 g C/m²
C_leaf_init      = N_leaf_init     × CN_L          = 2.0 × 25 = 50
C_root_fine_init = N_root_fine_init × CN_target_root = 5.0 × 35 = 175
C_stem_wood_init = N_stem_wood_init × CN_target_stem = 15  × 100 = 1500
C_fruit_init     = N_fruit_init    × CN_target_fruit = 0   × 30 = 0
C_storage_init   = N_storage_init  × CN_L          = 0 × 25 = 0
```

### 7b.4 演化方程 (每 KT, 跟 M2-B.1b 的 N 演化 镜像)

**Carbon flow 总图**:

```text
GPP (SCOPE biochem)
  ↓
C_labile 池 (积累)
  ↓ alloc (f_alloc_X 跟 N 一样)
[C_leaf, C_root, C_stem, C_fruit]
  ↓ Q10 维持呼吸 (R_m ∝ C, NoahMP 风格)
R_a → 出回大气 (NEE)

落叶 (frPHU>DLAI):
  C_leaf → retransl × C_storage  (跟 N retransl 同步, 用 retransl_fraction=0.7)
  C_leaf → litter × (1-retransl)  → 给师弟 SOM

根周转 (每日):
  C_root → litter                  → 给师弟 SOM

果实收获 (frPHU>0.90, Apple):
  C_fruit → cum_C_harvest (出地, 不给 SOM)

春萌 (0<frPHU<0.20):
  C_storage → C_labile  (跟 N spring_release 同步)
```

**公式**:

```matlab
% Alloc (类比 M2-B.1b allocatePlantN)
% [F2 stage-13a 已实现 2026-06-03, 验证 output_v41_F2leaftaper_2021]
% 叶分配不再用常数 f_alloc_leaf, 改加 FruitTree GDD-phase taper (解决 C_leaf 单调长到 ~403 失控):
%   g        = max(0, gdd5 − gdd5_at_budburst)                 [°C·d] 萌发后 GDD (base 5°C)
%   g_leaf   = (exp(−b) − exp(−b·g/GDDlfmat)) / (exp(−b) − 1)   [-] J-Norman 归一化递减, clamp[0,1]; g=GDDlfmat 时 =0
%   f_leaf   = f_alloc_leaf × g_leaf                            冠层成熟后叶分配→0, C_leaf plateau (不再到 403)
%   f_shaved = f_alloc_leaf × (1 − g_leaf)                      削掉份额 → 果窗内(DOY150-290)入 fruit / 窗外入 stem (Σf 守恒)
%   出处: Dombrowski et al. 2022 GMD 15:5167 Eq.4; CLM5-FruitTree
%         NutrientCompetitionFlexibleCNMod.F90:1663-1665; clm5_params.c171117__FruitTree.nc idx35 (苹果 GDDlfmat=1350, b=−0.5)
%   C 主导: 同一 taper 在 computeOrganDemand 镜像 → N 需求与叶 C 分配一致; 由 use_fruittree_alloc=1 门控 (当前仅苹果 PFT23)
A_C_leaf  = f_leaf   × NPP × dt   [g C/m²/KT]   (f_leaf  = f_alloc_leaf × g_leaf)
A_C_root  = f_alloc_root × NPP × dt
A_C_stem  = f_stem   × NPP × dt   (f_stem  = f_alloc_stem  + 窗外 f_shaved)
A_C_fruit = f_fruit  × NPP × dt   (f_fruit = f_alloc_fruit + 窗内 f_shaved; fruit 期外重分配, 跟 N 一致)

% C 池更新 (每 KT)
dC_labile  = GPP × dt − Σ A_C_X − R_a_extra      % NSC 平衡 (R_a 来自 GPP-NPP)
dC_leaf    = A_C_leaf  − R_m_leaf × dt − dC_leaf_senesce
dC_root    = A_C_root  − R_m_root × dt − dC_root_turnover
dC_stem    = A_C_stem  − R_m_stem × dt           % 无 turnover (多年生木质)
dC_fruit   = A_C_fruit − dC_fruit_harvest        % R_m_fruit 简化为 0
dC_storage = +dC_leaf_retrans − dC_spring_release

% 跟 N 池同步 litter/harvest (按 CN_target 换算)
dC_leaf_senesce  = dN_leaf_senesce × CN_L
dC_leaf_retrans  = dN_leaf_retrans × CN_L  (= dC_leaf_senesce × retransl_fraction)
dC_root_turnover = dN_root_turnover × CN_target_root
dC_fruit_harvest = dN_fruit_harvest × CN_target_fruit
dC_spring_release = dN_spring_release × CN_L

% 累积 (给师弟 SOM 接口)
cum_C_leaf_litter  += dC_leaf_senesce × (1 - retransl_fraction)
cum_C_root_litter  += dC_root_turnover
cum_C_harvest      += dC_fruit_harvest
```

### 7b.5 R_m 选择 (重要决策)

**A. R_m ∝ N (现状, M2-B.1b 已实现, LPJ Sitch 2003)**:
```matlab
R_m_X = r_m25_X_perN × Q10^((T-25)/10) × N_X × 1e-3
```

**B. R_m ∝ C (NoahMP/CLM5 严格物理)**:
```matlab
R_m_X = r_m25_X_perC × Q10^((T-25)/10) × C_X × 1e-3
% 需要重定 r_m25_X_perC, 不可直接用现 PFT 参数
```

**M2-B.8 决策**: **保留 A** (不动现有 R_m). 理由:
- 你 PFT 表 r_m25_root_perN=4.5e-7 已校准
- CN_target 守恒时 A/B 数值近似一致 (CN 偏离 <5%)
- 改 B 需要重校 5 个 PFT × 24 = 120 参数
- M3+ 阶段如果需要再切到 B

### 7b.6 实施步骤 (S1-S7, ~1.5 天)

```
S1 (0.5h) setInitialNitrogen.m 扩展为 setInitialPlantBGC.m (改名)
          - 加 C_labile/C_leaf/C_root_fine/C_stem_wood/C_fruit/C_storage 6 池
          - 初值 = N × CN_target

S2 (2h)   src/+plant/+PlantCarbon/computeCarbonAllocation.m (新建)
          - 类比 allocatePlantN.m, 但走 C 流, 用 f_alloc_X (PFT 表已有)
          - 输出 C 池更新 + cum_C_*

S3 (3h)   turnoverAndHarvest.m 加 C 同步
          - 4 子函数 (appleFruitHarvest, leafSenescenceRetransl, fineRootTurnover, springStorageRelease)
            全部加 C 版本 (按 CN_target 换算)
          - 主调度 turnoverAndHarvest.m 返回 6 C 池更新 + 4 cum_C_*

S4 (2h)   COMPASS.m 调度 hook
          - M2-B.4 N 块后跟 M2-B.8 C 块 (跟 N 同 KT)
          - 6 C 池写回 PlantBGCVariables.C_X(KT)

S5 (2h)   诊断: 11 字段
          - 6 prognostic: c_labile, c_leaf, c_root_fine, c_stem_wood, c_fruit, c_storage
          - 4 cum_*: cum_c_leaf_litter, cum_c_root_litter, cum_c_harvest, cum_C_storage
          - 1 守恒: plant_C_residual = GPP − R_a − cum_C_litter − cum_C_harvest − ΔC_pool

S6 (1h)   litter pack 接口 (给师弟 SOM 用)
          src/+plant/+integrated/exportLitterToSOM.m:
            litter_pack = struct( ...
              'C_leaf_litter', cum_C_leaf_litter_KT, ...
              'N_leaf_litter', cum_N_leaf_litter_KT, ...
              'C_root_litter', ..., 'N_root_litter', ...,
              'CN_leaf_litter', CN_L,    % 标识 litter 品质给 SOM 用
              'CN_root_litter', CN_target_root)

S7 (0.5天) 验证
          - 7 天 smoke (不破)
          - 全年 smoke (~3.5h)
          - 守恒检查: plant_C_residual < 1e-3 g C/m²/yr (机器精度)
          - C 池跟 N 池 CN_actual 监控 (应该 ≈ CN_target ± 5%)
          - 写 §12.26 完成存档
```

### 7b.7 验证标准

| 指标 | 期望 |
|:--|:--|
| C 池守恒 | plant_C_residual < 1e-3 g C/m²/yr (机器精度) |
| C-N 同步 | CN_leaf_actual ∈ [22, 28] (target 25 ± 12%) |
| GPP/NPP/R_a | 跟 M2-B.4 baseline 完全一致 (C 池只是镜像, 不影响 GPP) |
| C_storage 季节 | 秋积累 → 春释放 (跟 N_storage 同步) |
| C_fruit 收获 | DOY 245 frPHU>0.90 时一次性扣 90% |
| litter pack | 全年 cum_C_leaf_litter ≈ 50-100 g C/m²/yr (= cum_N_leaf_litter × 25) |

### 7b.8 给师弟 SOM 的接口

```matlab
% 师弟 SOM 模块每 KT 调:
[S_min_NH4_today, R_h_today, dCSOM] = soil.organic.SOMupdate( ...
    litter_pack_from_plant, ...    % 你 M2-B.8 输出
    SoilState, ...                  % 现有 SoilProperties + Theta + TT
    SOMpool_prev, ...               % 师弟自己的 SOM 池 state
    Delt_t);

% 师弟产出, 接到你模块:
NitrogenFlux.cvSmin(KT) = S_min_NH4_today  % 进 NH4 ADE (已有接口)
fluxes.Resp = R_h_today                     % 替代 SCOPE 简化 R_h, 修 NEE
```

### 7b.9 跟陆面扩展的兼容性

**预留接口** (M2-B.8 时只设字段, 不实现):

```matlab
% Trick state 三态 (CLM5 风格, M3+ 启用):
PlantBGCVariables.C_leaf_storage = []   % 字段占位
PlantBGCVariables.C_leaf_xfer    = []
% cfg.lLeafC_triState = 0 默认关, =1 启 CLM5 三态

% Multi-patch (陆面扩展):
% 当前 N/C 池都是 tS×1 标量 → 未来扩成 tS×Npatch
% 不动当前 setInitialNitrogen.m, 但写代码时注意 (:) ones-pattern
```

---

## 7c. 完整 SPAC W+E+C+N 真实循环 — 战略路线

### 终极目标

```text
真实陆面循环:
  大气 (Forcing P, Rn, T, RH, CO2, wind)
       ↓ ↑
  植物 (LAI, fwet, 光合, 蒸腾, alloc)
       ↓ ↑                            
  根系 (吸水 + 吸氮 + 释放呼吸 + 凋落根)
       ↓ ↑                            
  土壤 (水势, 温度, NH4, NO3, SOM, O2)
       ↓
  外界 (淋失, 蒸发, 渗漏)

闭环要素:
  ✓ Water 闭 (Richards + 蒸腾 + ETa)
  ✓ Energy 闭 (Rn = H + LE + G, 含冻融)
  ⏳ Carbon 闭 (GPP - R_a - R_h - NEE)
  ⏳ Nitrogen 闭 (input - leach - vol - harvest - litter ≡ pool change)
```

### 你的位置 vs 顶级模型

| 模型 | W | E | C | N | C-N 耦合 | 站点深度 | 陆面广度 |
|:--|:-:|:-:|:-:|:-:|:-:|:-:|:-:|
| **STEMMUS-SCOPE** (你上游) | ✓ | ✓ | ✓ GPP | ✓ N transport | △ 弱 | ⭐⭐ | ⭐⭐ |
| **DAISY** | ✓ | ✓ | ✓ DM 池 | ✓ NNI | ✓ 软限 | ⭐⭐⭐⭐ | ⭐ |
| **CLM5** | ✓ | ✓ | ✓ 19 池 | ✓ 19 池 | ✓ 独立 | ⭐⭐⭐ | ⭐⭐⭐⭐⭐ |
| **NoahMP** | ✓ | ✓ | ✓ 5 池 | ✗ | ✗ | ⭐⭐ | ⭐⭐⭐⭐ |
| **ClimaLand** | ✓ | ✓ | △ Pmodel | ✗ | ✗ | ⭐⭐ | ⭐⭐⭐ |
| **COMPASS 现状** | ✓ | ✓ | ✓ 6 池 (已完, 阶段5) | ✓ 6 池+完整循环 | ✓ LUNA | ⭐⭐⭐⭐ (有上 CLM5/DAISY) | ⭐ (站点) |
| **COMPASS 终目标** | ✓ | ✓ | ✓ 6 池 + 动态分配 | ✓ 6 池 | ✓ LUNA + NNI | ⭐⭐⭐⭐⭐ | ⭐⭐⭐ (多 patch) |

**COMPASS 独特优势**:
- 唯一同时有 DAISY+Nacry **单根 ADE** + CLM5 **LUNA** + DAISY 6 池 N 完整循环的模型
- Apple 站点率定 + 严格守恒 (机器精度) + multi-PFT 支持

### 完整路线

> ⚠️ 本处旧 "7 阶段 M1~M5" 路线图 (2026-05-20) **已过时并合并**。Plant C 6 池 (旧 M2-B.8)、turnover (旧 M2-B.4) 等
> **早已完成** (现到阶段 12/v10)。**统一路线图 (阶段 1-18 重排) 见文末 §12.33.5**, 以那个为准。
> 当前位置: 阶段 12 (v10 跑中) → 下一步 阶段 13a (苹果物候动态分配)。

**论文卖点**: 苹果园站点 **W+E+C+N 完整真实循环** + 动态分配 + (师弟 SOM 接通后) NEE 正符号 — 论文 1 (Tier 1 IF 4-6) 核心。

---

## 8. 物理推导精华 (保留)

### 8.1 ADE 方程 (HYDRUS / Šimůnek 2008 §3.4)

$$
(\theta + \rho_b K_d)\frac{\partial C}{\partial t}
= \frac{\partial}{\partial z}\!\left(\theta D_{sh}\frac{\partial C}{\partial z}\right)
- \frac{\partial(qC)}{\partial z} + S_{soil}
$$

其中:
- $C$: 液相浓度 [μg N / cm³ 液]
- $\theta$: 体积含水率 [m³/m³]
- $\rho_b K_d$: 阻滞因子 (NH4 有, NO3=0)
- $D_{sh}$: 弥散系数 [cm²/s]
- $q$: Darcy 通量 [cm/s]
- $S_{soil}$: 源汇 [μg N/(cm³ 土·s)] = 矿化 + 硝化 − 反硝化 − 根吸 − 挥发

### 8.2 DAISY 单根稳态柱 ADE (Kløffel 2026)

参考: DAISY 7.1.0 Appendix 7.2

单根周围 N 浓度场 $C(r)$ 满足径向 ADE, 解析解给出:
$$
I_{DAISY}(C_{root}) = \frac{2\pi r_r b_x D_x}{r_r - R_b} (C_{bulk} - C_{root})
$$
其中 $R_b = 1/\sqrt{\pi L_{rld}}$ (根间距).

### 8.3 Nacry 根膜 M-M (Nacry 2013)

$$
I_{Nacry}(C_{root}) = \frac{V_{max} \cdot (C_{root} - C_{min})}{K_m + (C_{root} - C_{min})}
$$

转 SI: Vmax 从 [μmol/(g 根·h)] → [kg N/(m 根·s)] = $V_{max} \times 10^{-6} \times \frac{1}{SRL \times 3600}$.

### 8.4 LUNA 5 池约束 (Ali 2016 Eq.7)

$$
F_{NCa} = N_{lc} + N_{et} + N_{cb} + N_{resp} + N_{store}
$$
$$
0.05 \le N_{lc}/F_{NCa} \le 0.95
$$
$$
N_{lc} + N_{et} + N_{cb} + N_{resp} \le 0.95 \cdot F_{NCa}
$$

### 8.5 LPJ NPP 完整分解 (Sitch 2003 Eq.8-14)

$$
NPP = \max\big(GPP - R_a, 0\big)
$$
$$
R_a = R_m^{leaf} + R_m^{root} + R_m^{stem} + R_g
$$
$$
R_m^i = r_{m25}^i \cdot Q_{10}^{(T_i - 25)/10} \cdot N_i \cdot 10^{-3}
$$
$$
R_g = R_{g,frac} \cdot \max(GPP - R_m, 0)
$$

### 8.6 NEE 完整化 (Chapin 2002 §6.4)

$$
NEE = R_h + R_a - GPP
$$

符号: NEE > 0 = 生态系统对大气是 CO₂ 源, NEE < 0 = 汇.

Apple 苹果园物理预期: **−400 ~ −800 g C/m²/yr** (净汇).
当前 +692 反符号, 因 SCOPE simple `R = 2.5 + 0.054·Ts` 全年 1164 g C/m²/yr 偏大. **M3 SOM 修复.**

---

## 9. PFT 表 (Apple PFT 23, 57 字段汇总)

详见 `src/+plant/getPFTParameterSet.m` 头部完整 doc (M2-B.3 重写完成).

### 9.1 已校准 Apple 关键值

```text
[A] Identification
  PFT = 'Apple', FullName = 'Apple (Domestic)'

[B] Attributes
  photo_type='C3', is_perennial=1, is_evergreen=0, has_fruit=1, has_storage_pool=1

[C] Stomatal (CLM Table 8.1)
  g1=4.86, g0=0.0001, m=9, alpha=0.06

[D] Structure (Wang 2009 / CLM5 BDT_temperate, M2-B.3 修正)
  SLA0=0.0125, SLAm=0.00133, CN_L=25, FLNR=0.1007

[E] Hydraulic (CLM5 Table 8.1)
  psi_o=-35000, psi_c=-224000

[F] Photosynthesis (Wang 2009 站点率定, legacy biochem 用)
  Vcmax25=105, Vcmax25_fN=33

[G] LUNA (Ali 2016 + COMPASS Apple 反推)
  Jmaxb0=0.0311 [const], Jmaxb1=0.1745 [site, Wang 2009 反推],
  tc_j0=0.8054 [const], H_LUNA=6.0996 [const],
  N_struct_per_LMA=0.002, N_photo_ref=1.6, Rd_fraction=0.015

[H] Organ C:N + alloc (Yuan 2011, Mao 2018, Sitch 2003)
  CN_target_root=35, CN_target_stem=100, CN_target_fruit=30
  f_alloc_leaf=0.40, f_alloc_root=0.25, f_alloc_stem=0.20, f_alloc_fruit=0.15

[I] Nacry root membrane (Nacry 2013, Wells 2002)
  Vmax_NH4_umol_g_h=12, Vmax_NO3_umol_g_h=5
  Km_NH4=0.3, Km_NO3=0.3, Cmin_NH4=0.5, Cmin_NO3=0.25, SRL=40

[J] Phenology (Niinemets 2010, Aerts 1996)
  leaf_long_days=200, retransl_fraction=0.7

[K] LPJ respiration (Sitch 2003)
  r_m25_root=4.5e-7, r_m25_stem=1.5e-7, Q10_resp=2.0, R_g_fraction=0.25

[L] Pool init (Wang 2009, Wells 2002, Tagliavini 2003, Mao 2018)
  N_leaf_init=2.0, N_root_fine_init=5.0, N_stem_wood_init=15.0, N_fruit_init=0

[M] Turnover + harvest (M2-B.4 待加)
  harvest_threshold=0.90, harvest_fraction=0.90,
  root_turnover_rate=0.00274, spring_release_endphu=0.20, spring_release_days=60
```

### 9.2 参数 Tunability 标签 (PFT 表 head doc)

| 标签 | 含义 |
|---|---|
| `[const]` | CLM5 Fortran `parameter`, 不可调 |
| `[tune]` | Per-PFT 校准 |
| `[site]` | COMPASS Apple 站点率定 |
| `[struct]` | 形态学, 随物种变 |
| `[reserved]` | M3 留位 |

---

## 10. 完整 SPAC 闭环现状 (2026-05-20)

```text
              GPP (SCOPE biochem + LUNA Vcmax25)
               ↓ ✓
         ┌─────────────┐
         │ Plant C/N   │  ✓ M2-B.1b (5 organ pool, ΔN 机器精度守恒)
         │  4 organ池  │
         └─────────────┘
            ↑       ↓
        N uptake  Turnover/Litter
       ✓ M2-B.2   ❌ M2-B.4 待
       demand
            ↑       ↓
         ┌─────────────┐
         │ Soil N      │  ✓ NH4/NO3 守恒机器精度
         │ Soil C/SOM  │  ❌ M3 SOM 待 (NEE 反符号)
         └─────────────┘
            ↑       ↓
         Water ✓  O₂ ✓

✓ N_leaf → Vcmax 反馈 (LUNA, M2-B.3)        ← 闭环最大缺口 已补
❌ Litter → SOM → S_min (M3)                  ← 下一阶段
```

---

## 10b. M2-B.4 完成存档 (2026-05-21) — Turnover + Harvest + N_storage

### 10b.1 实施总览

| S | 内容 | 计划 | 实际 |
|---|---|---|---|
| S1 | PFT 表 [M] 5 字段 × 24 PFT | 0.5 天 | 30 min |
| S2 | 4 子函数 + 主调度 + 单元测试 | 1 天 | 2h (6 单元测试机器精度 PASS) |
| S3 | COMPASS 调度 hook | 0.5 天 | 1h (含 struct field bug 修复) |
| S4 | 9 诊断字段 | 0.5 天 | 30 min |
| S5 | 7天 + 全年 smoke + 文档 | 0.5 天 | 1h (+v1→v2 frPHU bug 1 次重跑) |
| **合计** | | **3 天** | **同日完成** |

### 10b.2 关键文件改动

**新建**:
```
src/+plant/+integrated/+turnover/
  ├── appleFruitHarvest.m         一次性收获事件 (frPHU>0.90)
  ├── leafSenescenceRetransl.m    每日落叶 + 70% retransl 到 storage
  ├── fineRootTurnover.m          每日根周转 (1/365/day)
  ├── springStorageRelease.m      春萌储备释放 (0<frPHU<0.20)
  └── README.md
src/+plant/+integrated/turnoverAndHarvest.m  主调度 (caller-facing)
```

**修改**:
- `+plant/getPFTParameterSet.m`: 24 PFT × group [M] 5 字段 (harvest_threshold/fraction/root_turnover_rate/spring_release_endphu/spring_release_days). Apple 用 [site] 校准值, 其他用 EPIC/CLM5 默认.
- `+init/setInitialNitrogen.m`: 加 N_storage prognostic + cum_N_harvest/leaf_litter/root_litter 累积 + turnover_state + turnover_diag 预分配 + frPHU_diag.
- `COMPASS.m` line 1271-1330: M2-B.4 调度块. 在 M2-B.1b alloc 后调 turnoverAndHarvest, 写回 6 池 + 累积 cum_*.
- `COMPASS.m` line 2563-2577: 9 reg.write (含 frPHU 单独存数组).
- `+diagnostics/+groups/registerPlantNitrogen.m` line 173-220: 9 字段注册 (plant_n_storage + 4 dN_today + 3 cum_* + frphu).
- `+diagnostics/flushToCsv.m`: 9 CSV 兼容映射.

### 10b.3 调试历程 (2 bug)

| Bug | 症状 | 修复 |
|---|---|---|
| **#1 struct field 不匹配** | KT=1 后崩, "不同结构体之间下标赋值失败" | turn_diag 加 frPHU 字段 vs 预分配 struct array 4 字段冲突 → 拆开存 (turnover_diag 4 字段 + frPHU_diag 单独 tS×1) |
| **#2 frPHU 多除一次 PHU** | v1 全年 frPHU max=0.000165, 3 个 event 永不触发 | cfg.Veg.phen_HUI 已经 = CHU/PHU (归一化), COMPASS 不该再除 PHU. 修后 v2 max=0.648 |

### 10b.4 全年验证 (output_v40_M2B4_v2_full)

> ⚠️ **[已被取代 2026-05-30] 本节是 M2-B.4 frPHU 时代旧跑, harvest/spring_release=0 的结论已失效。**
> Phase D 重构 (2026-05-26) 把 harvest 改由**独立果实时钟 DS_fruit** 驱动 (`DS_fruit > harvest_DS=0.95`,
> seasonDecidPhenology.m:315-325 + turnoverAndHarvest.m:92 `[P5 fix]`), 与 frPHU=0.648 天花板**彻底解耦**。
> 最新 v12 全年跑 (`output_v40_phaseD_year2021_v12_cdom`) 实测: **cum_N_harvest=0.80 g/m²,
> cum_C_harvest=111.3 g/m²** (= 果实峰值 N_fruit=0.89/C_fruit=123.7 的 90%, 即 harvest_fraction ✅),
> N_harvest_today 1 次性触发 ✅。下表 0 值仅作历史记录, 勿据此判断现状。

#### 4 个 event 触发情况

| Event | KT 触发数 | 总量 | 评估 |
|---|---|---|---|
| **root_turnover** | 8760/8760 全年 | 4.94 g/m²/yr = **49.4 kg/ha** | ✓ daily, 跟理论 1/365/day × N_root 完美对齐 |
| **leaf_senesce** | 142/8760 | 0.031 g/m²/yr | ✓ 启动 (frPHU>0.45 时段), 量小因季节短 |
| **harvest** | 0/8760 | 0 | ⏳ frPHU max=0.648 没到 0.90 阈值 (留 known issue) |
| **spring_release** | 0/8760 | 0 | ⏳ N_storage=0.07 太小, 阈值 0<frPHU<0.20 时 N_storage 还=0 |

#### Pool 全年动态

| 池 | start | end | range |
|---|---|---|---|
| N_labile | 1.0 | 2.07 | 0.58-4.33 |
| N_leaf | 2.0 | 12.5 | 1.94-12.5 |
| N_root_fine | 5.0 | **4.81** | 4.11-5.68 (turnover 启动后小幅净降 ✓) |
| N_stem_wood | 15.0 | 16.2 | 15-16.2 |
| N_fruit | 0 | 2.05 | 0-2.05 (DOY 150-290 累积) |
| **N_storage** | 0 | **0.07** | 0-0.07 (leaf retransl 进 ✓ 首次启动!) |

#### Cumulative (送 SOM 师弟)

```
cum_N_harvest      = 0     g/m² = 0    kg/ha   (frPHU 没到 0.90)
cum_N_leaf_litter  = 0.031 g/m² = 0.3  kg/ha   (短季节 senesce)
cum_N_root_litter  = 4.94  g/m² = 49.4 kg/ha   (主导, 给 SOM 入口)
```

#### 整株 N 守恒

```
pool_start  = 23.00
pool_end    = 37.73
cum_lost    = 4.97
Total U_N   = 19.69 g/m² = 196.9 kg/ha    ← Wang 2009 苹果实测 ~180 kg/ha ✓
```

**完美** — 总根吸跟苹果实测吻合 (M2-B.3 v3 只有 22 kg/ha, M2-B.4 turnover 加上后自然回到合理量级).

#### Soil N + Water 守恒 (不破)

```
nh4_cum_residual   = -6.30e-6   (机器精度)
no3_cum_residual   = +2.37e-10
water_balance      = +0.006 mm/yr
GPP                = 1059.7 g C/m²/yr (跟 M2-B.3 baseline 完全一致 ✓)
NEE                = +679.4 (反符号待 M3 SOM)
```

### 10b.5 已知问题与限制

| 问题 | 影响 | 后续 |
|---|---|---|
| ~~frPHU max=0.648 < 1.0 → harvest 没触发~~ | ✅ **已修 (Phase D P5, 2026-05-26)** | harvest 改 DS_fruit 时钟驱动, 不再看 frPHU; v12 实跑 cum_N_harvest=0.80/cum_C_harvest=111.3 ✅ |
| N_storage 峰只 0.07 | spring_release 可释量小 | leaf retransl 进 storage 量小 (季节短); 多年生跨年累积后改善, 待 spin-up |
| Annual PFT 没验证 | Wheat/Maize 跑通 | M2-B.6 待 |
| 多年跨年没验证 | N_storage 跨年闭环 | 需跑 2-3 年 spin-up |

### 10b.6 数据归档

- 输出: `data/output_v40_M2B4_v2_full/` (8760 KT × ~254 diag fields, 9 个 turnover 字段新加)
- 失败版本: `output_v40_M2B4_S3/`, `output_v40_M2B4_S4/`, `output_v40_M2B4_full/` (v1 frPHU bug)
- 分析脚本: 见 §10b.4 内嵌

### 10b.7 下一步

**M2-B.8** (Plant C 6 池, DAISY 简化路线) — 详见 §7b. 跟 N 池 1:1 平行, 完成 W+E+C+N 完整 SPAC 真实循环.

---

## 11. 参考文献

### 11.1 物理基础 (土壤/水/N transport)

- **Šimůnek J, Šejna M, Saito H, Sakai M, van Genuchten MT** (2008). The HYDRUS-1D Software Package for Simulating the One-Dimensional Movement of Water, Heat, and Multiple Solutes in Variably-Saturated Media. *Department of Environmental Sciences, University of California Riverside*.
- **van Genuchten MT** (1980). A closed-form equation for predicting the hydraulic conductivity of unsaturated soils. *SSSAJ* 44:892-898.

### 11.2 根吸 (M1 + M2-B.2)

- **Nacry P, Bouguyon E, Gojon A** (2013). Nitrogen acquisition by roots: physiological and developmental mechanisms ensuring plant adaptation to a fluctuating resource. *Plant and Soil* 370:1-29.
- **Couturier J et al.** (2007). The expanded family of ammonium transporters in the perennial poplar plant. *New Phytologist* 174:137-150.
- **Wells CE, Eissenstat DM** (2002). Beyond the roots of young seedlings: the influence of age and order on fine root physiology. *Journal of Plant Growth Regulation* 21:324-334.
- **Kløffel et al.** (2026) DAISY 7.1.0 — Single-root cylindrical ADE framework (Appendix 7.2).

### 11.3 植物 C/N 分配 + 物候

- **Sitch S et al.** (2003). Evaluation of ecosystem dynamics, plant geography and terrestrial carbon cycling in the LPJ dynamic global vegetation model. *Global Change Biology* 9:161-185.
- **Mao L et al.** (2018). Apple fruit N budget (Mao series, citations from Apple-orchard literature).
- **Aerts R** (1996). Nutrient resorption from senescing leaves of perennials: are there general patterns? *Oikos* 67:469.
- **Niinemets Ü** (2010). A review of light interception in plant stands from leaf to canopy in different plant functional types and in species with varying shade tolerance. *Ecological Research* 25:693-714.
- **Tagliavini M et al.** (2003). Apple stem N budget.
- **Wang Y** (2009). Apple leaf N + Vcmax site-specific data (Loess Plateau apple orchards).

### 11.4 光合 + LUNA (M2-B.3)

- **Farquhar GD, von Caemmerer S, Berry JA** (1980). A biochemical model of photosynthetic CO₂ assimilation in leaves of C3 species. *Planta* 149:78-90.
- **Ball JT, Woodrow IE, Berry JA** (1987). A model predicting stomatal conductance and its contribution to the control of photosynthesis under different environmental conditions. *Progress in Photosynthesis Research* 4:221-224.
- **Medlyn BE et al.** (2011). Reconciling the optimal and empirical approaches to modelling stomatal conductance. *Global Change Biology* 17:2134-2144.
- **Evans JR** (1989). Photosynthesis and nitrogen relationships in leaves of C3 plants. *Oecologia* 78:9-19.
- **Kattge J, Knorr W** (2007). Temperature acclimation in a biochemical model of photosynthesis: a reanalysis of data from 36 species. *Plant Cell Environment* 30:1176-1190.
- **Bernacchi CJ et al.** (2001). Improved temperature response functions for models of Rubisco-limited photosynthesis. *Plant Cell Environment* 24:253-259.
- **Ali AA et al.** (2016). A global scale mechanistic model of photosynthetic capacity (LUNA V1.0). *Geoscientific Model Development* 9:587-606. **核心引用**.
- **Xu C et al.** (2012). Toward a mechanistic modeling of nitrogen limitation on vegetation dynamics. *New Phytologist* 193:30.
- **Lawrence DM et al.** (2019). CLM5 Technical Note. *National Center for Atmospheric Research*. **核心引用 Ch.9 LUNA**.
- **Rogers A** (2014). The use and misuse of $V_{c,max}$ in Earth System Models. *Photosynthesis Research* 119:15-29.
- **Suzuki T et al.** (2001). Rubisco half-life ~ 7 days (cited in LUNA enzyme turnover).

### 11.5 EPIC / DSSAT / phenology

- **Williams JR** (1989). EPIC: The Erosion-Productivity Impact Calculator. *Trans. ASAE* 32:497-511.

### 11.6 CLM5 (源代码参考)

- **`D:\OneDrive\CTSM-release-clm5.0.37\src\biogeophys\LunaMod.F90`** (1352 行) — LUNA 完整 Fortran 实现
- **`PhotosynthesisMod.F90`** — biochem hook (`lnc_opt=.false.` 路径)
- **`SoilWaterPlantSinkMod.F90`** — 根吸耦合 (CLM5 HydStress branch)
- **`namelist_defaults_clm4_5.xml`** line 430 — `Jmaxb1 = 0.093563` (CLM5 全球标定)

### 11.7 NEE / Carbon cycle

- **Chapin FS III, Matson PA, Mooney HA** (2002). *Principles of Terrestrial Ecosystem Ecology*. Springer. Ch.6 §6.4 — NEE.

### 11.8 SOM (M3 阶段待用)

- **Hansen S et al.** (1991). DAISY: Simulation model for water and N dynamics. *Fertilizer Research* 27:245-259.
- **Davidson EA et al.** (2012). The Dual Arrhenius and Michaelis-Menten kinetics model for decomposition of soil organic matter. *Soil Biology & Biochemistry* 56:175-179. (DAMM model, ClimaLand 用)
- **Koven CD et al.** (2013). The effect of vertically-resolved soil biogeochemistry on the carbon cycle in CLM4. *Biogeosciences* 10:7109-7131.
- **Lloyd J, Taylor JA** (1994). On the temperature dependence of soil respiration. *Functional Ecology* 8:315-323.

---

## 12. 历史归档

完整设计讨论 + 失败方案 + M2-A 中间设计 + §12.1-§12.20 早期 plan: 
`docs/_archive_cn_coupling_pre_rewrite_20260520.md` (5919 行).

---

## §12.26 M2-B.8 + M2-B.8b + M2-B.8c 完成存档 (2026-05-21 ~ 2026-05-22)

### §12.26.1 三阶段总览

| 阶段 | 日期 | 工作内容 | 状态 |
|---|---|---|---|
| **M2-B.8 v1** | 2026-05-21 | Plant C 6 池 DAISY 简化路线 (alloc + mirror + litter pack) | ✓ |
| **M2-B.8b NSC 联动** | 2026-05-22 上午 | R_a NSC cap + C_storage overdraft + 7 诊断字段 | ✓ |
| **M2-B.8c 二期重构** | 2026-05-22 下午 | COMPASS.m 主循环 325 行 → ~30 行 (函数抽象) | ✓ |

### §12.26.2 M2-B.8 v1 实施 (6 C 池 DAISY 简化路线)

#### 设计选择

3 个候选路线对比 (cn_coupling §7b):

| 选项 | 池数 | 难度 | 选 |
|---|---|---|---|
| A: CLM5 三态 (storage/transfer/display) | 18 | 高 | ✗ |
| B: DAISY 单态 + n 池 | 5-6 | 中 | **✓** |
| C: 不做 (跳 M3 SOM) | 0 | 低 | ✗ |

**理由**: DAISY 简化路线既保证守恒、又跟 N 池 1:1 mirror, 给师弟 SOM 提供清晰接口, 论文图可画 6 池季节动态.

#### 关键文件改动

**新建**:
- `+plant/+PlantCarbon/computeCarbonAllocation.m` (114 行) — 6 C 池 alloc, fruit gate 同 N
- `+plant/+integrated/turnoverCarbonMirror.m` (93 行) — N→C via CN_target
- `+plant/+integrated/exportLitterToSOM.m` (46 行) — 师弟 SOM 接口

**修改**:
- `+init/setInitialNitrogen.m`: 6 C 池 + cum_C_* 累积器初始化
- `COMPASS.m`: C alloc/mirror 调度 + 11 reg.write
- `+diagnostics/+groups/registerPlantNitrogen.m`: 11 字段注册
- `+diagnostics/flushToCsv.m`: 11 CSV 映射

#### 初值 (Apple, init from N × CN_target)

| 池 | N_init [g N/m²] | CN_target | C_init [g C/m²] |
|---|---|---|---|
| labile | 1.0 | 25 (CN_L) | 25 |
| leaf | 2.0 | 25 | 50 |
| root_fine | 5.0 | 35 | 175 |
| stem_wood | 15.0 | 100 | 1500 |
| fruit | 0 | 30 | 0 |
| storage | 0 | 25 | 0 |

#### 全年验证 (output_v40_M2B8_full)

- ✅ N 闭合机器精度 (NH4 res=-6.30e-6, NO3 res=+2.36e-10)
- ✅ 6 C 池演化正常 (C_leaf 50→312, C_fruit DOY 150-290 累积到 61.6)
- ✅ CN 比例守恒 (C/N_leaf=25.00, C/N_root=35.00, C/N_stem=100.00)
- ❌ **Plant C 守恒残差 = -104.32 g C/m²/yr (9.84% of GPP)** ← 触发 M2-B.8b 修复

### §12.26.3 M2-B.8b NSC 联动 (Bug 修复)

#### 根因

`+plant/+PlantCarbon/computeCarbonAllocation.m:96` 的 `max(C_labile, 0)` clamp:
- 5084/8760 步 (58%) 触发 — 夜间 R_a > GPP 时 C_labile 应该负
- clamp 凭空创造碳, 年累计 +104 g C/m² 假账
- 物理含义: 模型"赖账" — 真实植物会调用 C_storage 储备或下调 R_a (NSC starvation)

#### 物理修复路径 (生态生理学文献)

| Level | 条件 | 真实植物反应 | COMPASS 实现 |
|---|---|---|---|
| 0 | GPP > R_m | 全 R_m 由 GPP, 剩余 NPP 入 alloc | computeCarbonAllocation |
| 1 | GPP < R_m, C_labile 够 | GPP + C_labile 抽 | applyNSCCoupling 自然处理 |
| 2 | C_labile=0, C_storage 够 | C_storage → C_labile overdraft | **新增**: applyNSCCoupling Step 5 |
| 3 | C_storage=0 | 下调 R_m_actual (acclimation) | **新增**: applyNSCCoupling Step 2 cap |
| 4 | 持续 NSC=0 | C_starvation 死亡 | M5+ 不做 |

#### 实施 (5 处改动)

**新建**:
- `+plant/+integrated/applyNSCCoupling.m` (109 行) — Step 1-6 完整实现:
  1. R_m_potential = R_m_total + R_g
  2. NSC_avail = C_labile_prev + C_storage_prev
  3. R_a_actual = min(R_m_potential, GPP×dt + NSC_avail) — cap
  4. 调 computeCarbonAllocation 传 actual R_a + R_g=0
  5. C_labile<0 时从 C_storage 借 (overdraft, 守恒严格)
  6. 6 诊断字段 (R_m_potential / R_a_actual / starvation_event / nsc_total / c_storage_drawdown)

**修改**:
- `computeCarbonAllocation.m` (line 96): 删 `max(C_labile, 0)` clamp + 注释解释
- `COMPASS.m` (~1294-1306): 直接调用 applyNSCCoupling, 累积 cum_R_a (actual) + cum_R_m_potential + cum_starvation_count
- `+init/setInitialNitrogen.m`: 加 3 个 NSC 累积器 + nsc_diag struct array + C_pool_init_const 守恒锚
- `+diagnostics/+groups/registerPlantNitrogen.m`: 7 个新字段 (r_m_potential, r_a_actual, starvation_event, nsc_total, c_storage_drawdown, cum_r_m_potential, cum_starvation_count)
- `+diagnostics/flushToCsv.m`: 7 个 CSV 映射

#### Bug 副产: 守恒残差 KT=1 偏 -0.024 g

诊断中发现 `C_pool_init` 用 `NitrogenVariables.C_labile(1)` 读取, 但 KT=1 后此位被覆写 → 残差永远差一个 KT=1 step (-0.024 g/m², 这是 M2-B.8 v1 全年 -104 g 的次要源, 主源仍是 clamp 创造).

**修复**: `setInitialNitrogen.m` 加 `NitrogenVariables.C_pool_init_const = Σ C_*_0` (一次性, 不被覆写). COMPASS.m 残差公式改读这个字段.

#### 全年验证 (output_v40_M2B8b_nsc_full, 2.7h 跑)

| 指标 | M2-B.8 v1 | M2-B.8b NSC | 改善 |
|---|---|---|---|
| Plant C 残差 abs_max | -104.32 g/m² | **9.15e-12 g/m²** | 11 数量级 |
| 相对 GPP | 9.84% | **8.6e-13 %** | 机器精度 |
| C_labile clamp | 5084/8760 步 | **0** (overdraft 严格转移) | ✓ |
| starvation 诊断 | 隐藏 | **2578 events/yr 可见** | 物理透明 |

#### 苹果 C_storage 全年峰值仅 1.79 g (Tagliavini 实测 30-150 g)

不是 NSC bug, 是 **M2-B.4 phenology frPHU max=0.648 < 0.90 导致 leaf_senesce 没全启动** → leaf retransl 入金通道堵塞. 这是早就标记的 M2-B.4 known issue, 待 phenology 参数校准 (M3 阶段处理).

### §12.26.4 M2-B.8c 二期重构 (函数抽象)

#### 动机

M2-B.1b/M2-B.2/M2-B.3/M2-B.4/M2-B.8/M2-B.8b 累积下来, COMPASS.m 主循环 1085-1410 区段有 325 行 C/N 代码块, 跟主循环数值层混在一起. 违反 compass-coding-standard §4.1 (按耦合度切层) 和 §10 "签名 ≤8 参数".

#### 实施

**新建** (2 个 integrated 函数):

| 函数 | 行数 | 职责 |
|---|---|---|
| `+plant/+integrated/runPlantNDemand.m` | 145 | M2-B.2 demand 计算 (PFT 参数加载 + LAI target + NPP/D_N + NH4/NO3 拆分 + Nacry params) |
| `+plant/+integrated/runPlantCNStep.m` | 158 | 主调度: 串调 NPP + N alloc + NSC C alloc + N turnover + C mirror, 7 step 全部写回 NitrogenVariables |

**改动**:
- `COMPASS.m` 1085-1410 段 (325 行) → ~30 行调度块:
  ```matlab
  [demand_ctx, root_uptake_ctx] = vegetation.integrated.runPlantNDemand(...);
  [..., SinkS, SPUptake, ...] = vegetation.root.applyRootUptake(..., root_uptake_ctx.daisyParams, root_uptake_ctx.nacryParams);
  [NitrogenVariables, R_a_plant_kg] = vegetation.integrated.runPlantCNStep(NitrogenVariables, demand_ctx, U_N_total_g_m2_s, env_turn, Delt_t, KT);
  if isfield(fluxes, 'GPP'), [fluxes.NEE, ~] = vegetation.integrated.computeNEE(...); end
  ```

#### Bit-identical 验证

| 测试 | KT 数 | max |diff| vs M2-B.8b | 结论 |
|---|---|---|---|
| 7 天 smoke (Jun 1-8) | 168 | **0.0000e+00** | ✓ 位级一致 |
| 全年 smoke (2021) | 8760 | **0.0000e+00** | ✓ 位级一致 |

5 个 NC group (plant_nitrogen, carbon_fluxes, conservation, soil_nitrogen, soil_water, energy_balance) 全部 max|diff|=0. 重构纯结构改动, 零物理影响.

### §12.26.5 +plant/+integrated/ 模块全景 (M2-B 完成态)

```
src/+plant/+integrated/
├── +LUNA/                         M2-B.3 LUNA N→Vcmax 日调度
├── +turnover/                     M2-B.4 4 子函数
│   ├── appleFruitHarvest.m
│   ├── leafSenescenceRetransl.m
│   ├── fineRootTurnover.m
│   └── springStorageRelease.m
├── runPlantNDemand.m              【M2-B.8c】 N demand 计算
├── runPlantCNStep.m               【M2-B.8c】 主调度 (alloc/NSC/turnover/Cmirror)
├── applyNSCCoupling.m             【M2-B.8b】 NSC cap + C_storage overdraft
├── computeNEE.m                   【M2-B.1b】 NEE = R_h + R_a − GPP
├── turnoverAndHarvest.m           【M2-B.4】 turnover 主调度
├── turnoverCarbonMirror.m         【M2-B.8】 N→C via CN_target
└── exportLitterToSOM.m            【M2-B.8】 师弟 SOM 接口
```

### §12.26.6 数据归档

- 输出: `data/output_v40_M2B8c_refactor_full/` (8760 KT, 全 105+ 诊断字段, 含 7 NSC 字段)
- 早期版本: `output_v40_M2B8_full/` (v1, 守恒 -104), `output_v40_M2B8b_nsc_full/` (b, 守恒机器精度), `output_v40_M2B8c_refactor_smoke7d/` (c 7d 验证)
- 分析脚本: `tools/analyze_M2B8_full.py`

### §12.26.7 下一步

按 cn_coupling §7c (完整 SPAC W+E+C+N 战略路线):

| 阶段 | 内容 | 阻塞依赖 |
|---|---|---|
| **M3 SOM** | 师弟 3 池 SOM 模块 (litter→active→passive), 接 exportLitterToSOM 接口 | 师弟开发中 |
| **M2-B.5** | Apple 多年生**跨年 reset** 验证 (harvested_this_year 标志 + N_storage 跨年闭环) | harvest 单年已通 (DS_fruit, P5); 待 spin-up 验证跨年清零/累积 |
| **M2-B.6** | Wheat (annual PFT) 全年验证 | 跑 M2-B.5 完 |
| **多年 spin-up** | 2-5 年 cyclic forcing 跑稳态 | M3 完 |
| **M4 LUNA lnc_opt=.true.** | 精细化 LUNA 5 池约束 | 多年 spin-up 后 |
| **论文图** | 6 C 池季节动态 + NSC starvation 时序 + 5 苹果园 cross-site | 上面全 |

---

---

## §12.27 M2-B.8d frPHU 索引 bug 修复 (2026-05-23)

### §12.27.1 Bug 发现

`§12.26 M2-B.8c` 全年验证后画 frPHU 时序诊断图发现:

| 现象 | 数据 |
|---|---|
| frPHU max | 0.6482 @ **DOY 14.4 (1月中!)** ← 物理不可能 |
| frPHU end | **0.0** (年末归零) ← 全年波动错误 |
| frPHU > DLAI=0.45 step 数 | 142/8760 (4.7 天) ← 应该是秋季 60+ 天 |

苹果 1 月中旬不可能 frPHU=0.65, 长武冬季 T_air<0°C 全部. 时序物理违背 → bug.

### §12.27.2 根因

`COMPASS.m:1119` (M2-B.4 引入):
```matlab
frPHU_now = cfg.Veg.phen_HUI(KT);  % BUG: KT 是小时索引 (1-8760)
```

但 `+core/VegetationData.m:169` 注释明确写:
```
phen_HUI    % 逐日 frPHU (-)
```
长度 = n_days (365), 不是 8760. `VegetationData.m:478` 注释:
```
% 主循环中访问: day_idx = ceil(KT / steps_per_day), 再索引 phen_*(day_idx)
```

但 COMPASS.m 没遵守这个规则:
- KT=1..24 (DOY 1) 实际读 phen_HUI(1..24) → 把 24 天 frPHU 增量在 1 天内放完
- KT=365 (~DOY 15) 已把 phen_HUI 365 个值全部读完
- KT=366..8760 越界 → 被 `KT <= numel(phen_HUI)` 守卫拦截 → frPHU=0

观察到的 max=0.648 @ DOY 14.4 完全匹配: 14 天内累积了原本属于 14×24=336 天的 frPHU 增量, 但由于年内 turnover 触发也错位, 真实 frPHU 末值 = 0.648 (PHU 校准本身就是 2021 实际 GDD 不够全 PHU=3921).

### §12.27.3 影响范围 (干净的局部 bug)

| 模块 | 受影响? | 原因 |
|---|---|---|
| **LAI 模拟** | ❌ 不受影响 | 走独立路径 `set_step_veg → obj.LAI` (小时数组), 跟 phen_HUI(KT) bug 无关 |
| **能量平衡/SCOPE** | ❌ 不受 | 用 canopy.LAI 不用 frPHU |
| **N/C alloc** | ❌ 不受 | fruit gate 用 DOY (150-290), 不用 frPHU |
| **NSC 联动** | ❌ 不受 | 跟 frPHU 无关 |
| **turnover 4 事件** | ✅ **全瘫** | 全部依赖 frPHU |
| **C_storage 充值** | ✅ 间接 | leaf_senesce 未触发 → retransl 不进 storage |

LAI 校准 (RMSE=0.348 vs n=247 obs) **完全不受影响**, 这是用户原有的正确架构.

### §12.27.4 修复

`COMPASS.m:1116-1125` (3 行 → 5 行):

```matlab
% [M2-B.4 fix 2026-05-22] phen_HUI 是日数组 (n_days), KT 是小时索引 →
%   必须用 day_idx = ceil(KT / steps_per_day). 见 VegetationData.m line 478 注释.
%   旧 bug: 直接 phen_HUI(KT) 让 frPHU 在 DOY 14 (KT=336) 就 "跑完" 365 天,
%   max 卡 0.648, 之后全年归零 (越界 0). 4 个 turnover 事件 99% 不触发.
frPHU_now = 0;
steps_per_day = round(86400 / Delt_t);
day_idx = ceil(KT / steps_per_day);
if isprop(cfg.Veg, 'phen_HUI') && ~isempty(cfg.Veg.phen_HUI) && day_idx <= numel(cfg.Veg.phen_HUI)
    frPHU_now = cfg.Veg.phen_HUI(day_idx);
end
```

### §12.27.5 全年验证 (output_v40_M2B8d_frphufix_full)

| 指标 | 修复前 (8c) | 修复后 (8d) | 改善 |
|---|---|---|---|
| frPHU max | 0.648 @ DOY 14 | **0.648 @ DOY 324** | 时序物理正确 |
| frPHU end | 0.0 | **0.648** | 持续, 不再跳变 |
| cum_C_leaf_litter | 0.77 g C | **65.4 g C** | **85x** |
| cum_N_leaf_litter | 0.031 g N | **2.62 g N** | **85x** ✓ Wang 2009 |
| C_storage 峰值 | 1.79 g | **135.54 g** | **76x** ✓ Tagliavini 30-150 |
| C_leaf 末值 | 312 g (虚高) | **96.6 g** | 秋季正常落叶 |
| C starvation events | 2578 / 8760 | **1190 / 8760** | -54% (storage 救援) |
| C 守恒残差 | 9.15e-12 g | **9.27e-12 g** | 仍机器精度 ✓ |
| N 总质量 (pool+lost) | 42.69 g | **42.69 g** | 完全相同 ✓ |

### §12.27.6 仍未完成项

- ✅ **harvest 已修 (Phase D P5, 2026-05-26)**: ~~frPHU max=0.648 < 0.9 → 没触发~~。Phase D 把 harvest
  改由**独立果实时钟 DS_fruit** (bloom 后单独累 GDD, `DS_fruit=GDD_since_bloom/TSum_fruit`) 驱动,
  `DS_fruit > harvest_DS(0.95)` 一次性触发 (seasonDecidPhenology.m:315-325; turnoverAndHarvest.m:92 `[P5 fix]`),
  与 frPHU=0.648 天花板**彻底解耦** (不再依赖 PHU=3921 那个不可达阈值)。v12 实跑 cum_N_harvest=0.80 g/m²
  / cum_C_harvest=111.3 g/m² (= 果实峰值的 90% = harvest_fraction ✅)。
- ⏳ **spring_release 仍 0**: 春萌窗口 0<frPHU<0.20 时, N_storage 还=0 (上一年没充进来). 需要多年 spin-up 跑稳定循环.

### §12.27.7 数据归档

- 输出: `data/output_v40_M2B8d_frphufix_full/` (8760 KT, frPHU 索引修复版, 当前最新权威结果)
- 早期: `output_v40_M2B8c_refactor_full/` (索引 bug 对照版)
- 诊断图: `docs/figures/M2B8d_frphufix_validation.png` (6 子图)
- 分析脚本: `tools/analyze_M2B8_full.py`

---

---

## §12.28 M2-B.8e + M2-B.8f turnover 阈值校准 + harvest 单次触发 (2026-05-23)

### §12.28.1 问题诊断

§12.27 M2-B.8d 修了 frPHU 索引 bug 后, harvest 仍未触发 (cum_N_harvest=0). 调查发现两个独立问题:

#### Issue A: frPHU 阈值不匹配长武实际范围

长武苹果 (站点-气候-PHU 三者结合):

```
config_vegetation.txt:
  PHU = 3921 °C·day (GA+fmincon 拟合 247 obs LAI, RMSE=0.348, 2026-05-08)
  T_BASE = 5.59

长武 2021 实际:
  年累积 CHU ≈ 2510 °C·day
  → frPHU_year_end = 2510 / 3921 = 0.648 (永远不到 1.0)
```

**这不是 bug 是设计**: PHU 在 EPIC 框架里同时控制 LAI dHUF 增长曲线 + frPHU 端值. 如果减小 PHU 让 frPHU→1.0, LAI 会 6 月就达峰 (实测 9 月达峰), RMSE 飙升. 所以 GA 校准必然给 LAI 拟合优先, frPHU 端值≈0.65.

EPIC 默认 turnover 阈值 (假设 frPHU∈[0,1]):
- harvest=0.90, spring_release=0.20

在长武 frPHU∈[0,0.65] 下永远不触发.

#### Issue B: harvest 跨年误触发 (M2-B.8e 测试发现)

把阈值改成 0.60 后, harvest 首次 DOY 268 触发成功. 但发现 **KT=8760 又触发一次**.

根因: `current_date(KT)` 在 KT 步进**后**计算 (`current_date(KT) = start_time + seconds(TIME)`). EndTime=2022-01-01T00:00 时, KT=8760 (last) 实际 current_date = 2022-01-01. `turnoverAndHarvest.m` 用 `env.year` 判断 year rollover → `harvested_this_year=0` 重置 → 又允许 harvest.

### §12.28.2 M2-B.8e Fix A: PFT 表阈值缩放 (站点-specific)

**重要警告**: 缩放后的值仅适用长武站点, 别的苹果产区不能直接套用!

`getPFTParameterSet.m` line 1879-1908 (Apple PFT 23):

| 参数 | 原 EPIC 值 | 长武缩放值 | 计算公式 |
|---|---|---|---|
| `harvest_threshold` | 0.90 | **0.60** | EPIC 0.90 × 长武 frPHU 上限 0.65 |
| `spring_release_endphu` | 0.20 | **0.13** | EPIC 0.20 × 0.65 |
| `DLAI` (config_vegetation.txt) | 0.45 (EPIC) | **0.592** | LAI 校准选的, 跟 turnover 解耦, 不动 |

迁移到新站点操作清单 (写在 PFT 表注释):
1. 跑一年 smoke 拿到该站点 frPHU 年末值 (= CHU_year / PHU)
2. spring_release_endphu = 0.20 × frPHU_yr_end
3. harvest_threshold = 0.90 × frPHU_yr_end
4. DLAI 跟 LAI 校准保持一致, 不动

长期方案: 等多年 spin-up + 跨站点验证完成后, 用 frPHU_norm = frPHU / max(frPHU_yr) 归一化, 阈值统一回 EPIC 标准 (站点无关).

### §12.28.3 M2-B.8f Fix B: harvest 用 frPHU 下降检测重置

`turnoverAndHarvest.m` line 57-72 改动: **不用日历跨年判断, 用 frPHU 物理下降判断**.

```matlab
% [M2-B.8e fix 2026-05-23] 用 frPHU 春萌前重置 (而不是日历跨年):
%   1. 避免 EndTime=YYYY+1-01-01T00:00 时 KT=Nmsrmn 跨年误触发
%   2. 多年 spin-up 时, 真正的"新生长季"语义是 frPHU 从 0 重新启动
%   逻辑: state.frPHU_prev 跟踪上一 KT frPHU; 当前 frPHU<0.01 且上次>0.1 → 新生长季
if ~isfield(state_in, 'frPHU_prev') || isempty(state_in.frPHU_prev)
    state_out.frPHU_prev = env.frPHU;
    if env.frPHU < 0.01
        state_out.harvested_this_year = 0;
    end
elseif env.frPHU < 0.01 && state_in.frPHU_prev > 0.1
    state_out.harvested_this_year = 0;
    state_out.frPHU_prev = env.frPHU;
else
    state_out.frPHU_prev = env.frPHU;
end
```

`setInitialNitrogen.m` 加 `frPHU_prev: []` 字段.

### §12.28.4 全年验证 (output_v40_M2B8f_singleharvest_full, 2026-05-23 13:22 完成)

| 指标 | M2-B.8d (旧阈值) | M2-B.8e (改阈值) | M2-B.8f (修跨年) |
|---|---|---|---|
| Harvest 触发次数 | 0 | 2 (DOY 268 + KT 8760 误触发) | **1 (DOY 268 = 9月25日)** ✓ |
| cum_N_harvest | 0 g N/m² | 2.02 g | **1.71 g** ✓ (Wang 2009 长武苹果) |
| C_fruit 末值 | 61.6 g (堆积) | 1.04 g (二次触发后) | **10.45 g** (10% 籽实残留, 物理对) |
| C 守恒残差 | 9.27e-12 g | 9.27e-12 g | **9.27e-12 g** ✓ (机器精度) |
| N 总质量 | 42.69 g | 42.69 g | **42.69 g** ✓ |
| C_storage 峰值 | 135.5 g | 135.5 g | **135.5 g** (Tagliavini 30-150) |
| C starvation | 1190/8760 | 1190/8760 | **1190/8760** |

物理路径完整:
```
Apple 2021 长武全年 SPAC W+E+C+N 循环:
  GPP                = 1059.7 g C/m²
  NPP                = 486.4 g C/m² (NSC cap 后)
  Harvest (DOY 268)  = 1.7 g N (商业采摘 90% N_fruit)
  Leaf litter        = 2.6 g N (秋季落叶, 70% retransl 进 C_storage)
  Root litter        = 4.9 g N (全年周转)
  C_storage 秋充冬定 = 0 → 135 g (跟 Tagliavini 实测 30-150 g 对得上)
  C starvation       = 1190 events (冬夜 NSC 不足触发 cap)
  
  C 守恒残差: 9.27e-12 g (机器精度)
  N 总质量守恒: 42.69 g (位级)
```

### §12.28.5 剩余 1 issue (留 M2-B.6 / spin-up 处理)

`spring_release` **仍 0 触发**: 春萌窗口 (0 < frPHU < 0.13) 时 N_storage 还是 0 (上一年没充进来). 需要多年 spin-up 让 N_storage 跨年闭环 (秋季 leaf retransl 充进, 春季 release 放出).

预期效果: 跑 2-3 年 cyclic forcing 后, year 2 春萌应该看到 N_storage → N_labile 释放, 真正启动跨年储备路径.

### §12.28.6 数据归档

- 输出: `data/output_v40_M2B8f_singleharvest_full/` (8760 KT, **当前最新权威结果**)
- 关键文件改动:
  - `+plant/getPFTParameterSet.m` Apple PFT 23 (M2-B.8e 阈值 + 30 行长武警告注释)
  - `+plant/+integrated/turnoverAndHarvest.m` line 57-72 (M2-B.8f frPHU 下降检测)
  - `+init/setInitialNitrogen.m` line 155 (turnover_state 加 frPHU_prev)
- 删除: M2-B.1b/B.2/B.3/B.4/B.8a-d/B.8e 所有中间产物 (保留 v32-v39 baseline)

### §12.28.7 M2-B 系列完整闭环达成 (6 阶段)

| 阶段 | 核心 | 验证 |
|---|---|---|
| M2-B.8 v1 | 6 C 池 DAISY 简化 | ✓ 6 池演化, CN 比例守恒 |
| M2-B.8b | NSC 联动 | ✓ C 残差 11 数量级 → 机器精度 |
| M2-B.8c | 二期重构 | ✓ COMPASS.m 325→30 行 bit-identical |
| M2-B.8d | frPHU 索引 fix | ✓ C_storage 1.79→135 g 进实测区间 |
| M2-B.8e | PFT 阈值缩放 | ✓ harvest 首次 DOY 268 触发 |
| **M2-B.8f** | **harvest 跨年 fix** | ✓ **单次触发, 守恒机器精度** |

---

*下一个完成存档: §12.29 (M3 SOM 阶段, 师弟开发) 或 §12.30 (多年 spin-up).*

---

## §12.30 自养呼吸校正 — R_m_stem sapwood + R_m_leaf 光抑制诊断 (2026-05-28, v6)

### §12.30.1 起因 — R_a/GPP = 109% 碳收支为负

v5 全年验证 (LUNA strict + Jmaxb1=0.244) 后查 C 池, 发现:

- **C_labile 全年钉 0** (DOY 63 起), 全年 42% 时段 = 0
- **1278 次 starvation event**, 全集中在 7-10 月生长季 (非冬季)
- 日碳收支: **R_a/GPP = 109%** (GPP 3.28 < R_a 3.58 g C/m²/day) — 呼吸比光合还多, 必然抽干 NSC

按器官分解 7-10 月维持呼吸 (g C/m²/day):
| 项 | v5 | 占 GPP |
|---|---|---|
| R_m_stem | 2.73 | **54%** ← 最大单项 |
| R_m_leaf | 1.98 | 39% |
| R_m_root | 0.44 | 9% |

### §12.30.2 根因 1 — R_m_stem 用全树干 N (含心材)

`computeNPP.m` 旧式:
```matlab
R_m_stem = r_m10 × N_stem_wood × g(T)    % N_stem_wood = 19.24 g N/m² (活+死全部)
```

物理错误: **只有边材 (sapwood / livestem) 活组织呼吸, 心材 (heartwood / deadstem) 是死组织不耗 O₂**。
大树 80-90% 木质是心材。C_stem_wood = 1924 g C/m² (38.5 t/ha, 生物量本身合理), 但全部拿去算呼吸 → R_m_stem 高 4.5 倍。

**文献验证 (SKILL §5.5 合规)**:
- CLM5 `CNMRespMod.F90:265`: `livestem_mr(p) = livestemn(p)*br*tc` — 只用 livestem N, deadstem 不进式 (已核实源码行)
- CLM5 `NutrientCompetitionCLM45defaultMod.F90:281`: `f4 = flivewd(ivt(p))` — 活木分配比例 (已核实)
- Sitch 2003 GCB 9:161 状态变量表 p.163: maintenance respiration ∝ "sapwood C", 区分 sapwood-heartwood (PDF data/references/.._Sitch_..pdf 第 3 页, 已核实)

### §12.30.3 修复 (B 方案: 单池 × fraction, 压缩态)

引入 `f_livewood` (活木占树干比例), 等价于 CLM5 livestem/total 比值:
```matlab
% computeNPP.m:99-106
f_livewood  = pftParams.f_livewood (缺字段 fallback 1.0)
R_m_stem_kg = r_m10 × (N_stem_wood × f_livewood) × g(T)
```

PFT 表 (`getPFTParameterSet.m`) 24 个 PFT 加 `f_livewood`:
- **苹果 PFT 23 = 0.22** [site-calibrated] — 反推 R_a/GPP≈0.70 (果园密植小树边材高占比); 无苹果直接文献; 介于 CLM5 flivewd 乔木~0.1 与果园 0.3 之间
- 其余乔木 = 0.10 (CLM5 flivewd 温带乔木默认)
- 草本/作物 (PFT 15/20/21/24) = 1.0 (无木质 stem 占位)

透传链: `runPlantNDemand.m:62` pftP_demand 加 f_livewood → computeNPP。

**升级路径 (A 方案, 留轮 3 = 阶段 13c)**: 把单 N_stem_wood 池显式拆成 livestem/deadstem (sap/heart) 双池,
f_livewood 变成初始分配比例, 加 Sitch Eq 27 活木→死木周转。

> **[F3 关联 2026-05-31] 边材 N 储备容量项 (待双池补全)**: F3 根吸储备目标用 LPJ max_n_storage =
> max(边材容量项, 来年萌发项) (von Bloh 2018; growth.cpp:1490)。**边材容量项 = C_sap × fnstorage / CN_L**
> 需边材碳 C_sap。当前**无 sap/heart 双池**, F3 用**压缩态代替**: C_sap ≈ C_stem_wood × f_livewood (0.22),
> 即 1121×0.22≈247 g C/m² → 容量项 ≈3.76 g N/m² < 萌发项 6.47, 故苹果 max 取萌发项。
> **拆双池 (13c) 后**: C_sap 成独立预报池, 边材容量项改用真 C_sap (随树龄/边材累积动态变, 大树边材多→容量项升,
> 可能超萌发项主导), 替掉 f_livewood 压缩态。fnstorage=0.3 (LPJ 落叶树 europe.ins)。
>
> **[F3 整体并入此处 2026-06-01]** F3 need-based 三版均未压住 N_labile (LAI桶5.5 / kstore21.8 / need-based24.6);
> 失败根因: 萌发项储备目标 `(C_leaf+C_root)/CN×nreloc` 用**当季膨胀** C_leaf (盛期364.9 vs 初值89.6) 致目标飙到~18,
> 且漏 LPJ `×turnover_leaf` 缩放。**根治需 13c 拆 sap/heart 双池**: max_n_storage 用真 C_sap (年尺度稳定量, 不随当季叶碳膨胀)
> + 补 turnover 缩放, 才能让储备目标稳在实测量级 (~几 g N/m²)。**F3 暂挂, 与边材双池 13c 一并解决**。
> **代码保留 need-based 版不回退** (用户决定 2026-06-01): N_labile 虚高 24.6 但守恒机器精度✓, 不阻塞 F1/F5/后续;
> 量级偏大记为**已知项**待 13c 根治。后续验证看到 N_labile 偏高时, 知道是此待解项, 非新 bug。

#### 📦 §12.30.3a 13c 边材双池收口包 (2026-06-01) — 同根问题一并解决
拆 sap/heart (livestem/deadstem) 双池是**以下 4 个缠在一起问题的共同根**, 13c 时一次性收口:
| 症状 | 现状 | 双池后根治 |
|---|---|---|
| **茎周转缺失** (N3) | N_stem_wood/C_stem_wood **只进不出**, 逐年虚涨 +3.5%/yr (turnoverAndHarvest/turnoverCarbonMirror 无 stem 分支) | Sitch 2003 Eq27 边材→心材转化 + 死木 litter |
| **F3 储备目标** | 萌发项用当季膨胀 C_leaf → N_labile 虚高 24.6 | max_n_storage 边材项用真 C_sap (年稳定量) |
| **边材容量** (GL2) | f_livewood=0.22 压缩态代替 (单池×分数) | C_sap 独立预报池, 随树龄/累积动态变 |
| **边材呼吸** (§12.30 已用压缩态) | R_m_stem = N_stem×f_livewood (已生效, B方案) | 改用真 livestem N (A方案), 心材不呼吸 |
- 共同前置: 把单 `N/C_stem_wood` 拆成 `livestem`(边材,活,呼吸+储备容量) + `deadstem`(心材,死,不呼吸) 双池, f_livewood(0.22) 变初始分配比 + Sitch Eq27 周转率。连锁 ~14 文件 (大工程)。
- 出处: Sitch 2003 GCB 9:161 Eq27 (sap→heart); CLM5 CNLivewoodTurnover / livestem_mr; LPJ growth.cpp cmass_sap/heart。
- **在此之前**: 茎虚涨/N_labile偏高/边材压缩态 均记为**已知项**, 不单独修 (修了也是局部补丁, 双池才根治)。

##### ✅ 13c 实施完成 + Jan-Aug smoke 验证 (2026-06-01, output_v18_13c_smoke)
拆双池全部实施 (Tasks 1-10, 11 文件 + VOD): getPFTParameterSet (CN_lw=50/CN_dw=500/flivewd/lwtop=0.7) /
setInitialNitrogen (7池 C主导初值) / computeCarbonAllocation+allocatePlantN (flivewd 拆) /
turnoverAndHarvest+turnoverCarbonMirror (lwtop livestem→deadstem 周转) / computeNPP (R_m_stem 用真 N_livestem) /
runPlantCNStep+runPlantNDemand (双池时序 + F3 sap_term) / COMPASS+register+flushToCsv (7池守恒+诊断) /
compute_VOD+load_compass_state (live/dead 加权 WDR_live=1.5/dead=0.4)。checkcode 全 CLEAN。

Jan-Aug smoke 四条验收:
| # | 标准 | 结果 | 判定 |
|---|---|---|---|
| 1 守恒 | PlantC/N 机器精度 | C_res 5.68e-12, N_res 1.63e-13 | ✅ 完美 |
| 2 茎不虚涨 | live+dead 不暴涨 | 1121→1250 (+11.5%/8mo); **lwtop 周转生效** (livefrac 0.22→0.14, deadstem 875→1071) | ✅ 机制对 |
| 3 N_labile | <~10 | init 2→end 24.9 (mean 18.5) | ❌ **未达标 (见下)** |
| 4 R_a/GPP | 0.7-0.85 | 0.513 | ⚠️ 偏低 (见下) |

**结构性目标全达成** (双池/周转/守恒/边材呼吸正确化), 提交 13c。两个量级问题**另开后续任务**:

###### 遗留-1: N_labile 仍虚高 24.9 (F3 老账, 13c sap_term 未压住)
诊断 (smoke 数据): sap_term 峰值仅 **2.17** (C_livestem 181 / CN_L 25 × fnstorage 0.3), regrow_term **11.6** (当季 C_leaf 367 膨胀),
`max(sap, regrow)` 全程 100% 取 regrow → **sap_term 形同虚设**, 设计 §6 期望落空。
更关键: N_labile_target 峰值才 11.6, 但实际 N_labile 囤到 25 (>target 2 倍), storage_deficit 早 =0,
→ **N_labile 虚高根因不是 F3 公式, 是结构性**: N 进 labile (根吸 U_N + retransl 回收) > 器官当季消耗, 多 N 囤 labile **无出口**
(不像 C 撞顶剩碳留 labile 是刻意; N 这边没有过剩外排机制)。f_N_growth 49% 步数=1 (N 充足)。
**修法 (后续任务)**: ① F3 需求加"N_labile 超 target 停 D_N_growth 根吸" (现只停 storage_deficit); 或 ② 学 CLM5 retransn 上限 + 过剩不吸。

###### 遗留-2: R_a/GPP 实测判据 = Zanotelli CUE 0.71 (⚠️ 纠错: 旧"0.7-0.85"无出处)
**判据纠错 (2026-06-02)**: 之前 §12.30.4/.7 标的 "健康果园 R_a/GPP 0.7-0.85" **是拿模型自身旧输出 (v6=78%/v7=71%) 反凑的循环论证, 无文献出处** (违 §5.5)。已作废。
**真实测判据**: Zanotelli et al. 2013 Biogeosciences 10:3089 p.3091 报该成龄苹果园 **CUE = NPP/GPP = 0.71 ± 0.12** → **R_a/GPP = 1−CUE = 0.29 ± 0.12 (≈0.17-0.41)**。
**对照**: v19_falloc 全年 R_a/GPP = **0.497** (CUE=0.503), 对实测 CUE 0.71 **偏低 → 即 R_a 偏高/呼吸吃碳过多** (我此前误判为"偏低", 方向反了)。
注: v6=78%/v7=71% 那时被夸成"健康", 实际比 Zanotelli 实测的 R_a/GPP≈0.29 高很多 (CUE 偏低), 一直是问题, 只是判据错把它当达标。
**修法 (后续任务)**: 模型 CUE 0.50 < 实测 0.71 = 呼吸偏高/GPP 偏低。查 R_m 量级 (livestem/root N×r_m10) 或 GPP (LUNA Vcmax)。目标 CUE→0.71 即 R_a/GPP→0.29。

###### ⚠️⚠️ R_a/GPP 口径错误 (2026-06-02, 影响所有历史判断)
**发现**: 历史 R_a/GPP (含 v6=78%/v7=71%, 及本人 diag 脚本) 都用 **GPP ≈ Σ(NPP+R_a)** 当分母, **错**。
因 computeNPP `NPP=max(GPP−R_m−R_g, 0)` 有 max 截断: 冬季/夜间 GPP≈0 但 R_a=R_m 仍累加 → NPP=0 而 R_a>0 → (NPP+R_a) **虚高真 GPP**。
v19 实测: Σ(NPP+R_a)=1884 vs 真 GPP (cum_gpp_plant 内部账)=1523, **虚高 361 (24%)**。
**正确口径** (用 cum_gpp_plant=ebal 真光合积分):
| 量 | 错口径(NPP+Ra) | **真口径(cum_GPP)** | Zanotelli 2013 实测 |
|---|---|---|---|
| GPP gC/m²/yr | 1884 | **1523** | 1346 (偏高13%) |
| R_a/GPP | 0.497 | **0.614** | 0.29 |
| CUE | 0.503 | **0.386** | 0.71±0.12 |
| R_m_leaf/GPP | 28.7% | **35.5%** | — |
**真相**: CUE 真实只 0.386 (非 0.50), 离实测 0.71 **更远**; R_a 是实测 2 倍; **R_m_leaf(叶 Rd 暗呼吸) 吃 GPP 35.5%** 是头号问题。
历史 v6/v7 "71-78% 健康果园" 判断**双重错** (①0.7-0.85无出处 ②NPP+Ra 错口径), 实际一直离实测远。
**Rd 链查实 (biochemical.m:317-367)**: Rd25=0.015×Vcmax25 (Farquhar C3 标准, 有据); 温响应对标 CLM5; v7 已加白天光抑制×0.67。**公式无错** → Rd 偏高 ⟺ **Vcmax25 偏高** (Rd∝Vcmax, GPP 也∝Vcmax, 同源)。LUNA 跑出 Vcmax25 均 95.6/峰 122。
**⚠️ 缺苹果叶 Vcmax25 实测靶** (搜到苹果 Asat~31 口径不同; 冬小麦 Vcmax25 50-133 非苹果) → 不能拍"Vcmax 偏高", 待找长武/文献实测靶再定 (同 f_alloc: 先找靶再标)。
**待办**: ① 修所有 R_a/GPP 诊断改用 cum_gpp_plant 分母; ② 找苹果 Vcmax25 实测靶; ③ 靶到后判 LUNA Vcmax 是否偏高 / Rd 光抑制是否够。

###### 🎯 CUE 偏低真根因 = N_leaf 偏高 (2026-06-02, Greer 2018 亲验)
**纠错教训**: 本人先凭印象认定"苹果 Vcmax25 应 60-90, 模型 95.6 偏高", 差点据此调 LUNA。查 **Greer 2018 Env Exp Bot (1-s2.0-S0098847218300625, 果园实测 cv. Cox's Orange, A/Ci)** 后发现: 苹果叶 **Vcmax 量级百级** (正文 p4 "carboxylation ranged 25-225 μmol m⁻² s⁻¹", 32°C 最优), **模型 95.6 不偏高甚至偏低**。印象错, §5.5 防的坑 (没出处别当判据)。
**R_a 拆项 (真口径 cum_gpp_plant=1523, R_a=cum_r_a_potential=936 两口径一致可靠)**:
| 项 | gC/m²/yr | %GPP | vs 实测 |
|---|---|---|---|
| R_m_leaf | 540 | **35.5%** | 头号超标 |
| R_m_root | 135 | 8.8% | |
| R_m_stem | 144 | 9.4% | |
| R_g | 117 | 7.7% | 正常 |
| **R_a 合计** | **936** | 61.4% | 实测~415, **超 2.3×** |
**真根因链 (闭合)**: 模型 **lnc(叶N)=3.64 gN/m²叶 vs Greer 实测 2.0-2.3 (1.7-2.2%DM, p9)** = **偏高 1.69×** → lnc 高 → LUNA Vcmax∝lnc 算高 + Rd∝Vcmax (Rd25=0.015×Vcmax) 高 → R_m_leaf 吃 35.5%GPP → R_a 超 2.3× → **CUE 0.386 偏低**。
**上游 = 叶 C:N (CN_L=19.7) 偏低**: 峰值 C_leaf 185/N_leaf 9.66 → 叶 C:N=19.1 (=设定 CN_L 19.7)。Greer 叶 N 1.7-2.2%DM, 叶 C 木本~45-50%(待出处) → 实测叶 C:N≈22-28, **模型 19.7 偏低** → 同叶碳配 N 多 1.5-1.7×。
⚠️ **CN_L=19.7 PFT 表无出处注释** (违 §5.4)。**但 2026-06-02 查实: COMPASS 6 个器官 C:N (叶19.7/根32/活木50/死木500/果139) 极可能本来就抄自 Dabbaghi 2022 FruitTree Table C1 (Appendix C, Zanotelli 2010 六次测量 + 2013 实测, 同行评审苹果专用)** — 只是 PFT 表漏写出处注释, 非凭空拍。
⚠️⚠️ **CUE 诊断待重审**: 若 CN_L=19.7 真是 Dabbaghi 苹果实测 (而非本人合成靶 24), 则"CN_L 偏低需改 24"**作废**, N_leaf 偏高 1.69× 的根得另找 (N 分配 A_leaf / lnc 口径 / LAI)。Table C1 PDF OCR 行列错乱 (leafCN 可能 19.7 也可能 32, 读不准), **待用户亲读 PDF 确认每个器官准确值再定**。
**叶 C:N 实测靶 (2026-06-02, 两端亲验, 合成靶)**:
- 叶 N = **1.95% DM** (Greer 2018 中值, 1.7-2.2%, 亲验 Table; 成熟叶多源 2.0-2.4%)
- 叶 C = **46.59% DM** (Ma et al. 2018 Biogeosciences 15:693 Table, **Deciduous broad-leaved 46.59±3.55%**, 亲验 line 129; PDF webfetch 已存 tool-results; 全球叶均 46.9%, 非传统 50%)
- → **叶 C:N = 46.59/1.95 = 23.9 ≈ 24** (范围 N 1.7-2.4% → C:N 19-27)
- 对照模型 **CN_L=19.7 偏低 ~18%** (落范围低端), 印证 N_leaf 偏高。
- ⚠️ 合成靶 (Greer 叶N × Ma 叶C, 不同研究相除, 非单篇苹果直测 C:N); 但两端均亲验实测, 远胜凭印象。
**待办**: ① 标 CN_L 19.7→~24 (有 Greer+Ma 出处); ② 验 N_leaf 降→lnc 降→Vcmax/Rd 降→CUE↑ (目标 0.71); ③ 注意 CN_L 同时影响 LUNA 输入 + 叶 N 需求, 改后看 GPP/产量是否仍合理 (CN_L↑→同叶碳配 N 少→可能 N_leaf 降但 LAI/GPP 受影响, 需全年验)。

###### ⚠️ CN_L 标定作废 + Jmaxb1 率定失败 (2026-06-02, 亲读 Dabbaghi Table C1)
**CN_L=19.7 是 Dabbaghi 苹果实测正确值** (Table C1 leafCN=19.7, Zanotelli 2010 六次测量), 上面"该改 24"的合成靶 (Greer×Ma) **作废**。同表亲读纠正: **CN_lw 50→60** (livewdCN=60), **flivewd 0.22→0.15** (=BDT CLM5)。6 个器官 C:N 已补 Dabbaghi 出处注释。
**Jmaxb1 率定 (0.244→0.18) 基本失败**: 目标 Vcmax25~78 (Dabbaghi 公式 5.2+34×lnc)。但全年跑 (output_v21_cab_jmaxb1, cab_from_luna=1 叠加) 实测 Vcmax25 仅 95.6→**89.8** (降 6%), 远没到 78 — **LUNA 非线性补偿**: Jmaxb1 降 26% 时 LUNA 把更多 N 分给羧化补回, Vcmax 几乎不动。**结论: 靠调 Jmaxb1 压 Vcmax 这条路走不通** (LUNA 自有 N 优化平衡点)。
**v21 四项 (cab_from_luna=1 + Jmaxb1=0.18, 真口径)**:
| 指标 | v19 | v21 | 目标 | |
|---|---|---|---|---|
| Cab 动态(Route A) | — | 66.8-95.5 均83.9 | ~40-55 | ⚠️ 比离线估高1倍但不崩/不触guard(<120), 接近原固定100故光谱冲击小 |
| Vcmax25 | 95.6 | 89.8 | ~78 | ⚠️ LUNA补偿, 没压下 |
| **CUE (真cum)** | 0.386 | **0.429** | 0.71 | ⚠️ 微升, 远不够; R_a/GPP 0.614→0.571 |
| R_m_leaf/GPP | 35.5% | 34.4% | — | 仍头号 |
| 产量(×148) | 45789 | **46868** | 49315 | ✅ 差5% |
| C_leaf峰 | 185 | 186.7 | ~185 | ✅ |
| 守恒 C/N | — | 1.2e-11/1.4e-13 | 机精 | ✅ |
**归因受限**: cab+Jmaxb1 叠跑, CUE 微升 0.386→0.429 分不清是哪个 (cab 改光谱 vs Jmaxb1)。
**新方向 (2026-06-02 决策)**: CUE 偏低根**不在 Jmaxb1** (压不动 Vcmax)。改查 **Rd/光抑制** (R_m_leaf=Rd∝Vcmax, 占 GPP 34%): Rd25/Vcmax 比 0.015 是否偏高 / 白天光抑制 ×0.67 够不够 / 或 R_m_root+stem。f_alloc(产量) + cab 动态(不崩) 已成功, 待提交; CUE/呼吸作独立线继续查。
Greer 2018 PDF 已存 data/references/fruit_tree/ (Table 2 OCR 行列错位, 精确 Vcmax25/Jmax25 不可靠配对, 仅用正文定性"百级"+叶N%DM 硬值)。

###### ✅ v22 验证 SLA0=0.028 修 lnc (2026-06-02, output_v22_sla028, cab_from_luna=1)
**SLA→lnc→Vcmax 链验证成功**: SLA0 0.0125→0.028, lnc=1/(SLA×CN_L)=1/(0.028×19.7)=1.81 (从4.06降, 落 Greer 实测 2.0-2.3)。
| 指标 | v19(SLA0.0125) | **v22(SLA0.028)** | 目标/实测 |
|---|---|---|---|
| Vcmax25 | 95.6 | **均61.7/峰84.2** | Dabbaghi公式~67-78 ✅落区 |
| Jmax25/Vcmax25 | — | 1.97 | Walker 1.6-1.7 (偏高但 Jmax 不缺) |
| CUE(真cum) | 0.386 | **0.464** | 0.71 (升但不够) |
| R_a/GPP | 0.614 | 0.536 | 0.29 |
| R_m_leaf/GPP | 35.5% | 27.4% | — |
| GPP | 1523 | **1178** | 1346 (⚠️偏低12%) |
| 产量(×148) | 45789 | **40898** | 49315 (差17%) |
| C_leaf峰 | 185 | 169.6 | — |
| Cab动态 | — | 均73.8(47-94) | 40-60 (v21 84, 降了点) |
| 守恒 | — | 1.0e-11/1.3e-13 | 机精 ✅ |
**关键: Rd 排除 (用苹果实测直接校, 不用比值)**: 模型 R_m_leaf(Rd) per 叶面积 = 峰 **1.34**/均 0.64 μmol/m²/s, vs **苹果实测 Braeburn 0.96/Fuji 1.31 (PMC2759240, Rd+Vcmax 实测)** — **模型 Rd 量级正常甚至偏低, 不偏高!** 排除"Rd 偏高致 CUE 低"。(Rd/Vcmax=0.015 是 Collatz 标准, 苹果实测比 0.006-0.017 跨幅大对 Vcmax 口径敏感, 不动。)
**⚠️ 判据重大纠错 (2026-06-02, 亲读 Zanotelli 原文 zan2013.txt): GPP 实测靶是 1263±189 不是 1346!** 之前一直挂的 "实测 1346" 来源不明 (疑记串), 违 §5.5。亲验 Zanotelli 2013 BG 10:3089 原文 (行 57/478/726, Table 6):
- **GPP = 1263 ± 189 g C/m²/yr** (EC 法; 另一独立法 1375±12, 取 EC 为主靶)
- **R_a = 303 ± 202 g C/m²/yr** (行 555, Eq 8); NPP = 380 ± 30 (行 57/478); **CUE = 0.71 ± 0.12** (行 42)
- max LAI = 2.8 (行 121, ≈模型 2.79 ✅)

**模型 v22 重新对标 (用真靶 1263±189)**:
| 量 | 模型 v22 | 实测 Zanotelli | 评判 |
|---|---|---|---|
| GPP | **1178** | **1263±189** (带宽 1074–1452) | **✅ 落误差带内!** (近中心) — **GPP 不偏低!** |
| R_a | **631** | **303±202** (上限 505) | ❌ 偏高 ~2× (即便实测上限 505, 模型 631 仍超) |
| CUE | 0.464 | 0.71±0.12 (下限 0.59) | ❌ 偏低 |
| LAI | 2.79 | 2.8 | ✅ |

**结论彻底翻转: GPP 根本没偏低 (1178 落 1263±189 误差带), 此前 "GPP 偏低 12%" 是错靶 1346 导致的误判, 作废。CUE 偏低的全部责任在 R_a 偏高 (631 vs 303±202)。** 与 Greer 2018 (fp18093.pdf, 苹果光合温度响应实测) 一致: 苹果 30°C 以下卡 Rubisco/Vcmax (Greer 行 273-284, Vcmax 实测 141.7 低N / 峰~225 @32°C); 模型 Vcmax 均 61.7 < 实测 141.7, 说明**模型靠较低 Vcmax + 配置凑出对的 GPP, 但呼吸侧 R_a 算高了**。
**待办 (聚焦 R_a, 不再碰 GPP)**: 拆 R_a = R_m_leaf 322.6(27.4%)/R_m_root 130.1(11%)/R_m_stem 81.3(6.9%)/R_g 96.8(8.2%) = 630.7。要降到实测 303 (砍 ~328), Rd(叶)已验正常 → 偏高在 **R_m_root/stem (合 211, r_m10=0.066 Sitch 通用值, 苹果是否该更低待查 Greer/Zanotelli Rd+暗呼吸) 或 R_g**。**注意 §5.5: r_m10/R_g_fraction 都有 Sitch/CLM5 出处, 改前必须先找苹果专属呼吸实测靶, 没找到就标 SOURCE UNVERIFIED 不擅改。**

###### ✅ R_a 终判 + r_m10 反推验算 (2026-06-02): 光调现有参数救不了 CUE, Rd 是真瓶颈但正常
**反推 r_m10 (纯计算, Rm_root+Rm_stem = k×r_m10, k=3203; Rm_leaf=Rd 来自 SCOPE biochem 不受 r_m10 控制):**
| 目标 | 需要 r_m10 | 可行性 |
|---|---|---|
| CUE=0.71 (中心) | **−0.0264** | ❌ 负数, 不可能 |
| Ra=372 (Dabbaghi 实测) | −0.0157 | ❌ 负数 |
| Ra=505 (Zanotelli 上限带, CUE~0.57) | 0.031 | ✅ 可行 (0.066→0.031 砍半) |

**关键: 即便 r_m10=0 (根/茎完全不呼吸), CUE 最多 0.646, 仍到不了 0.71** — 因为 **Rm_leaf(Rd)=323 占 GPP 27%, 是 SCOPE biochem Rd25=0.015×Vcmax 给的, 不归 r_m10 管**。
**Rd 单点已验正常**: Rd per 叶面积生长季 = 均 0.687/峰 1.186 μmol/m²leaf/s, 落苹果实测 Braeburn 0.96/Fuji 1.31 (PMC2759240) 区间内甚至偏低。Rm_leaf 年累积高 = Rd × LAI 2.26 × 全年时长, 非系数错。
**模型 R_a=631 落在 Zanotelli 实测 303 与 Dabbaghi CLM5-FruitTree 模型值 693 之间** (比专门的苹果 CLM5 模型还好一点)。
**终判 (§5.5)**: ① **Rdparam 验证过正常, 不能为凑 CUE 下调 (会破坏单点匹配)**; ② r_m10 是唯一合法可调旋钮且有 Dabbaghi "成熟果园维持成本高估" 佐证, 但 **Bepete & Lakso 1997 PDF 不在本地, 不据未读论文改**; ③ Rg_fraction=0.11 已 F1 改对。**所以现阶段无任何参数可在有出处前提下改降 R_a。CUE 0.46 是 CLM5/Sitch 系模型对苹果的系统性偏差 (Dabbaghi 自承 CUE 0.59 同偏), 论文如实写, 待 Bepete&Lakso 出处再修。**

###### ✅ g1/g0 死值清理 + 出处确认 (2026-06-02)
**git 溯源**: PFT 表 g1=4.86/g0=0.014 源自 ec84225 (2025-12-26), 原始注释 `gsP.g1 = 4.86; % 调参参数` — **手调凑值, 无文献出处** (§5.5 典型坑)。
**运行时真相**: biochemical.m 消费 `leafbio.medlyn_g1/g0`, 由 COMPASS.m:610-611 从 **config_vegetation.txt gs_g1/gs_g0** 读 (get_veg_property), **不读 PFT 表 g1/g0** → PFT 表那 48 个 (24×2) 是死值。
**g1 早已改对**: config gs_g1=**8.2** (你 2026-05-29 v8.1 改, 出处 Dabbaghi 2022 GMD 15:5167 Table C1 medlynslope=8.2 + Zanotelli 2015 苹果园 GPP/ET 实测调), 与本次亲读 Dabbaghi 原文挖到的 8.2 吻合。gs_g0=0.0001 与全 PFT 一致。
**已删** getPFTParameterSet.m 全部 48 行死值 g1/g0 + doc-string 标注 "g1/g0 read from config not PFT table" (checkcode 0 msgs)。免再被死值误导 (本轮就先被 4.86 误导过)。

###### ✅ v23 cab-off 对照: 拆 cab 动态对 GPP 的贡献 (2026-06-02, output_v23_caboff = cab_from_luna=0 + SLA0.028)
**目的**: v22(cab 开) 减 v23(cab 关, 退回 SCOPE 固定 Cab=100) 隔离 cab 动态 (Cab 100→均74) 对 GPP 偏低的责任。
| 指标 | v22 (cab ON) | v23 (cab OFF) | Δ (cab 动态贡献) |
|---|---|---|---|
| GPP_cum (cum_gpp_plant) | 1177.7 | 1163.8 | **+13.9 gC/m² (+1.20%)** |
| R_a_cum (cum_r_a_potential) | 630.7 | 628.4 | +2.3 |
| CUE = (GPP−R_a)/GPP 真口径 | 0.464 | 0.460 | +0.004 |
| 产量 C_fruit_max×148 | 40898 | 40147 | +751 kg/ha |
| Vcmax25 均/峰 | 61.74/84.24 | 61.76/84.24 | ≈0 (cab 不动 lnc/Vcmax, 符合预期) |
| Cab 均 (47-94 范围) | 73.8 | 73.7 | 残留 ~0 (caboff 仍写 cab_from_nlc 诊断但不喂光谱) |
| 守恒 C/N 残差 | 1.0e-11/1.3e-13 | 1.3e-11/2.6e-13 | 机精 ✅ |

**结论 (锁定)**: **cab 动态对 GPP 仅 +1.2% (+13.9 gC/m²), 不是 GPP 偏低 12% 的原因。** 之前预期 "Cab 100→74 减吸光大幅拉低 GPP" 被推翻 — Cab 73.8 仍在叶绿素吸收饱和段 (PROSPECT 吸收对 Cab>40 已近饱和, 边际吸收很小), 故降 Cab 对冠层吸光影响微弱。**GPP 偏低 12% (1178<实测1346) 的责任几乎全在 LUNA/气孔侧, 与遥感正演 cab 解耦** → 后续查 GPP 只盯 Medlyn g1 / 叶温 / Jmax 电子传递, 不必再纠结 cab。**且: cab 动态可安全保留** (产量/CUE 影响均 <2%, 不崩不触 guard), 遥感正演线 RS1 (Cab 动态) 与碳氮主线 (CUE/GPP) 互不干扰, 可独立推进。

**⚠️ NPP 字段口径坑 (再确认, §5.5 防伪)**: 分析中发现 `net_primary_productivity` 输出字段累积 = 782.8 gC/m² ≠ GPP−R_a = 547.0 (差 236)。根因 = `computeNPP.m:63` 的 `NPP = max(GPP−R_m−R_g, 0)` 在冬季/夜间 GPP<R_m 时把负 NPP 截成 0, 累积丢失净碳损失 → NPP 字段虚高 24% (与旧 "冬季 max(NPP,0) 虚高" 坑同源)。**真 CUE 必须用 (GPP−R_a)/GPP = 0.464, 不能用积分 NPP 字段 (0.665 虚高)。** R_a 各分量 (Rm_leaf 322.6/Rm_root 130.1/Rm_stem 81.3/Rg 96.8 = 630.7) 与 cum_r_a_potential 精确自洽。

###### F2-③ 架构改造决策: 暂挂记档 (2026-06-02)
leafc_sla 已验证**可切不崩** (v19 标定 f_alloc 后, leafc_sla LAI 峰 2.62 vs EPIC 2.79 = 0.94; 生长季比值 0.74-0.92; 峰晚 30d; r=0.777; 图 output_v19_falloc_cal/f2_lai_compare.png)。
**有意识排序决策 (非搁置)**: F2-③ 是把 LAI 从预计算搬到在线的架构大改 (解 hotfix spin-up + C_leaf 滞后喂), 高风险压轴。现在切 = 用大风险换暂时更差的 LAI (生长季低 15%+晚 30d), 且会与未解的呼吸偏差 (CUE 0.386) 耦合, 难定位。
**先清单变量**: ②呼吸标定 (CUE 0.386→0.71, 最大物理偏差, 有 Zanotelli 实测靶) → 修峰值晚 30d (senescence 相位) → **最后**上 F2-③ 合 C/N→GPP→Vcmax 闭环。
当前卡点: 呼吸标定缺苹果 Vcmax25 实测靶 (Rd∝Vcmax)。

### §12.30.4 v6 全年验证 (output_v40_phaseD_year2021_v6)

| 指标 | v5 | v6 | 评判 |
|---|---|---|---|
| R_m_stem (KT 1-719 对照) | — | **0.220×v5** | ✅ 精确生效 (=f_livewood) |
| **R_a/GPP 全年** | 109% | **78%** | ✅ 进健康果园区 (0.7-0.8) |
| GPP 日均 | 3.28 | 3.80 | ↑16% (呼吸降→NPP 升) |
| Starvation 全年 | 1278 | **776** | ↓39% |
| C_labile 砸 0 占比 | 42% | 33% | 改善 |
| N_leaf 年末 | 0.107 | 0.306 | ↑ |
| N_storage 年末 | 11.21 | 14.41 | ✅ |
| **C_storage 年末** | 151.57 | **244.85** | ✅ ↑62% (储备增厚) |
| 守恒 NH4/NO3/PlantC/水 | — | -6.3e-6/5.8e-11/-5e-4/0.007mm | ✅ 全机器精度 |

### §12.30.5 残留问题 — R_m_leaf 光抑制缺失 (待修, 下一轮)

R_m_stem 修好后, 最大呼吸项变成 **R_m_leaf = 54% GPP** (7-8 月果实膨大期 4.05 g C/m²/day)。
诊断发现 **白天 R_m_leaf (4.63) > 夜间 (3.28)** — 异常, 健康模型白天 Rd 应被光抑制 < 夜间。

根因: COMPASS `biochemical.m:358` `Rd = Rd25 × f_Arr × f_Peak` 只有温度响应, **缺白天光抑制**。

**CLM5 标准做法 (待引入, 已核实源码)**:
- `PhotosynthesisMod.F90:1562-1564`:
  ```fortran
  ! 0.67 from Lloyd et al. 2010 & Metcalfe et al. 2012
  if ( light_inhibit .and. par_z(p,1) > 0._r8 ) then  ! 白天有光
     lmr_z(p,iv) = lmr_z(p,iv) * 0.67_r8               ! 暗呼吸 ×0.67
  end if
  ```
- COMPASS biochem Rd 温度响应参数 (DeltaS_Rd=490, lmrha=46390, lmrhd=150650) 与 CLM5 line 1239/1248/1249 完全一致, 只差最后这步 ×0.67 光抑制。
- biochem 已有 `Q` (吸收 PAR) = CLM5 `par_z` 等价物, 可直接判白天。

**[已完成 v7, 见 §12.30.7]**: biochemical.m:360-367 已加 `Rd = Rd × (0.67 + 0.33×(Q<=0))` 光抑制。v7 验证 R_a/GPP 78%→71%, 白天 R_m_leaf≈夜间 (反常修正)。

### §12.30.6 数据归档

- 代码: `computeNPP.m` (f_livewood), `runPlantNDemand.m` (透传), `getPFTParameterSet.m` (24 PFT f_livewood), `ebal_clm45.m`+`ebal_clm45_profile.m` (v5 LUNA strict 3 字段透传)
- 输出: `data/output_v40_phaseD_year2021_v6/` (全年 8760 KT, 18:52 完成)
- v5 baseline 保留: `data/output_v40_phaseD_year2021_v5/`

### §12.30.7 R_m_leaf 光抑制实施 (v7, 2026-05-29 完成)

§12.30.5 诊断的"白天 R_m_leaf > 夜间"反常已修复。

**实施** (`biochemical.m:360-367`):
```matlab
Rd = Rd .* (0.67 + 0.33 .* double(Q <= 0));   % Q>0(白天)→×0.67, Q≤0(夜间)→×1.0
```
- 出处: CLM5 `PhotosynthesisMod.F90:1562-1564` light_inhibit 分支, 系数 0.67 来自 Lloyd et al. 2010 (Funct Ecol) + Metcalfe et al. 2012 (New Phytol) (已核实源码行)
- 用乘法 (0.67 + 0.33×夜间标志) 代 CLM5 的 if 分支, 避免 Q/Rd 维度不一致时逻辑索引错位; 数学等价
- Q = 吸收 PAR [μmol m⁻² s⁻¹] (biochem line 114), >0 即白天有光

**Ci 求解方式核实 (附带, 不改)**: COMPASS 解 Ci 用阻尼定点迭代 (biochemical.m:441-454, damp=0.5),
非纯解析。测试验证迭代精度 0.01% (误差 0.03 ppm, 10 次收敛), CLM5 也是 brent 数值求根 (非纯解析),
两者本质都是数值法, COMPASS 现状够用, 不改。边界层 gb 修正在 ebal `heatfluxes_clm45.m:107`
`Cc = Ca − (Ca−Ci)·ra_w/(ra_w+rs)` 实现 (Ca→叶面 Cs), biochem 内 Cs 已是叶面值, 故 Ci=Cs−1.6A/gs
只需气孔项, 不缺 gb (与 CLM5 合并式 Ci=Ca−An·(1.4/gb+1.6/gs) 拆分等价)。

**v7 全年验证 (output_v40_phaseD_year2021_v7, 02:02 完成)**:

| 指标 | v5 | v6 | v7 | 评判 |
|---|---|---|---|---|
| **R_a/GPP 全年** | 109% | 78% | **71%** | ✅ 进健康果园核心区 0.65-0.75 |
| R_m_leaf 白天 (7-8月) | — | 5.36e-8 | **3.65e-8** | ✅ ↓32% (=×0.67 精确) |
| R_m_leaf 夜间 | — | 3.80e-8 | 3.59e-8 | 基本不变 ✓ |
| 白天 vs 夜间 | — | 白天>夜间 ❌ | **白天≈夜间** ✅ | 反常修正 |
| Starvation 全年 | 1278 | 776 | **728** | ↓ |
| C_labile 砸0 | 42% | 33% | 32% | 小改善 |
| N_storage 年末 | 11.21 | 14.41 | **15.47** | ✅ |
| C_storage 年末 | 151.57 | 244.85 | **270.15** | ✅ ↑10% vs v6 |
| GPP 峰 | — | 29.68 | 29.53 | ✓ 几乎不变 (光抑制不动 GPP) |
| 守恒 NH4/NO3/PlantC | — | — | -6.3e-6/1.7e-11/-5e-4 | ✅ 全机器精度 |

**三轮呼吸校正完整链**: v5 R_a/GPP=109% → [R_m_stem sapwood] → v6 78% → [R_m_leaf 光抑制] → v7 71% ✅

### §12.30.8 残留 — 果实膨大期日内 C_labile 缓冲不足

v7 仍 32% 时段 C_labile=0 + 728 starvation, 集中 7-8 月果实膨大期 (碳需求峰)。
R_m_stem/R_m_leaf 已修, 残留来源候选:
1. **C_storage→C_labile 日内缓冲通道缺失** — 当前 C_storage 仅春萌 (spring_release) + overdraft (砸0后) 才流向 C_labile, 缺"低于阈值主动回补"。CLM5 cpool/storage 双向。
2. R_m_root 或果实分配 (A_C_fruit) 偏猛
3. 夜间 GPP=0 但 R_m 持续, C_labile 当天耗尽 (应靠 storage 缓冲)

**下一步候选**: (a) C_storage→C_labile 日内缓冲 (applyNSCCoupling 加目标库存回补, ~15 行); (b) 轮 3 livestem/deadstem 双池 + Sitch Eq 27; (c) v8 lnc_opt=.true. 土壤 N→光合。

*数据: v5/v6/v7 baseline 均保留 data/output_v40_phaseD_year2021_v{5,6,7}/ (各 ~1.1G, periodic cutoff 已删).*

---

## §12.31 vs FATES 碳代谢对比 + 碳饥饿致死决策 (2026-05-29)

FATES 源码已下载到 `D:/OneDrive/fates/` (NGEET/fates GitHub, shallow clone)，可合法引用 (§5.5)。

### §12.31.1 storage_drawdown_event 改名 (旧名 starvation_event)

旧名 starvation_event 易误解为"碳饥饿致死"。**实际只标记"本步 R_m>GPP, 动用 NSC 储备"** — 苹果树正常 NSC 周转 (白天存糖/夜间-果实膨大期取储备), 不触发任何死亡。

改名映射:
| 旧 | 新 |
|---|---|
| `starvation_event` | `storage_drawdown_event` |
| `cum_starvation_count` | `cum_storage_drawdown_count` |
| CSV `StarvationEvent` | `StorageDrawdownEvent` |
| CSV `cum_StarvationCount` | `cum_StorageDrawdownCount` |

改动文件: applyNSCCoupling.m / runPlantCNStep.m / setInitialNitrogen.m / flushToCsv.m / registerPlantNitrogen.m / COMPASS.m (诊断字段重命名, 不改物理)。

### §12.31.2 决策: COMPASS 不加 FATES 碳饥饿致死

**FATES 碳饥饿致死** (EDMortalityFunctionsMod.F90:201-228):
```fortran
frac = store_c / target_leaf_c                          ! 储备碳 / 满冠目标叶碳
cmort = mort_scalar_cstarvation * max(0, (upthresh-frac)/upthresh)   ! 线性
   或 = mort_scalar_cstarvation * exp(-frac/upthresh)               ! 指数
```
判断量是 **储备/目标叶碳比 frac 持续低于阈值** → 致死概率 cmort (减少 cohort 数量 n)。

**COMPASS 不加, 理由 3 条**:
1. **无 cohort/数量概念**: FATES 是个体木队列模型 (cohort 有 dbh/age/n), cmort 减少 n (林分疏开)。COMPASS 是单株大叶模型, 模拟 1 棵定植苹果树, 没有"死几棵"的作用对象。
2. **人工管理果园**: 苹果园灌溉/施肥/修剪, 树不会碳饥饿自然死亡, 跟 FATES 模拟的自然林分演替/竞争/更新是完全不同场景。
3. **CLM5 大叶模型自己也不加**: CLM5 主干 CNGapMortalityMod.F90 用固定 gap mortality 速率, 不做碳饥饿致死。碳饥饿致死只存在于 FATES 这类动态植被模型。COMPASS 对标 CLM5 大叶, 跟着不加。

**FATES 全部 7 个致死项** (EDMortalityFunctionsMod.F90:60): bmort(背景)/cmort(碳饥饿)/hmort(水力失败)/frmort(冻害)/smort(尺寸衰老)/asmort(年龄衰老)/dgmort(损伤) — 全部是 cohort 尺度, COMPASS 单株大叶都不适用。

### §12.31.3 vs FATES 碳代谢差异清单 (COMPASS 可改进项)

| 机制 | FATES | COMPASS 现状 | 评判 |
|---|---|---|---|
| 维持呼吸 cap | curmr_ratio=gpp/mr (CLM5 同) | R_a cap = min(R_m, GPP+NSC) | ✅ 等价 |
| 叶呼吸光抑制 | (CLM5 ×0.67) | v7 已加 ×0.67 | ✅ 已对齐 |
| sapwood 呼吸 | 只活木 | v6 f_livewood=0.22 | ✅ 已对齐 |
| **储备目标量 target_storage** | `bstore_allom(dbh)` 按异速生长定目标储备 (EDCohortDynamicsMod.F90:1442) | ❌ **无储备目标概念** | COMPASS 缺: C_storage 只被动充放, 无"目标库存", 故无法做"低于目标主动回补" |
| **春萌储备调用** | store_c → 按 deficit 比例 flush 各器官 (phenflush_fraction, EDPhysiologyMod.F90:1694) | spring_release 固定速率 1/60d 释放到 labile | ⚠️ COMPASS 简化 (固定速率 vs deficit 驱动), 可接受 |
| NSC 池数 | store (单池, 异速目标) | C_labile + C_storage (双池) | 架构不同, 各自自洽 |
| 致死 | 7 项 cohort 致死 | 无 (人工果园, 见 §12.31.2) | ✅ 故意不加 |

**唯一实质缺失**: COMPASS 无 **target_storage (储备目标量)** 概念。FATES 用异速生长 `bstore_allom(dbh)` 给每棵树定一个"该存多少碳"的目标, 用于:
- 碳饥饿致死判断 (frac=store/target) — COMPASS 不需要 (不致死)
- 分配时优先补足 storage 到目标 — **COMPASS 可借鉴**: 让 C_storage 在低于目标时优先回补, 解决 §12.30.8 的 C_labile 砸 0 (果实膨大期)

但 target_storage 的异速公式 `bstore_allom` 是 FATES 个体木异速 (dbh→储备), 苹果园单株大叶无 dbh 动态, 需 site-calibrate 一个固定目标 (如初始 C_storage 205 gC/m²)。属可选改进, 非 bug。

*数据/源码: FATES @ D:/OneDrive/fates/ (可验证引用).*

---

## §12.32 动态分配方案设计 — LPJ-GUESS 异速+功能平衡 (2026-05-29, 设计完成待实现)

> **导航**: 本节 §12.32.1-.7 是初版设计 (LPJ 一统), §12.32.8 细化, **§12.32.9-.11 是最终定案**
> (漏洞排查后改为**三套分配分流**: 苹果/作物走 FruitTree 物候动态, 草走 LPJ grass, 自然树走 LPJ pipe model)。
> **以 §12.32.11 最终架构为准**。初版 §12.32.1-.7 "LPJ 一统 4 步" 已被 §12.32.10 #4 决策修正
> (苹果不走 LPJ 乔木异速, 改物候分配)。实施优先级见 §12.32.11 末。

### §12.32.1 问题 — 固定 f_alloc 违反两条植物分配规律

COMPASS 现状: `C_organ = NPP × f_alloc_organ`, f_alloc 固定 (苹果 leaf 0.40/root 0.25/stem 0.20/fruit 0.15)。
**两条真实分配规律都违反**:
1. **功能平衡** (Brouwer 1963): 缺水→多分根, 缺光→多分叶。固定比例不响应环境。
2. **异速生长** (pipe model): 树越大→茎占比越高 (DBH 3→26 cm, stem 分配应升)。固定比例不随树龄变。

### §12.32.2 选型 — LPJ-GUESS vs FATES (选 LPJ-GUESS)

源码: LPJ-GUESS v4.1.1 `guess_4.1/modules/growth.cpp` allocation() (D:/OneDrive/LPJ-GUESS-v4.1.1/)。
FATES: `D:/OneDrive/fates/biogeochem/` (PARTEH)。

| 维度 | LPJ-GUESS | FATES | 选 |
|---|---|---|---|
| 功能平衡 (缺水多分根) | ✅ 内置 ltor (leaf:root) | ⚠️ 弱 | LPJ |
| 异速 (树龄→茎) | ✅ pipe model k_latosa | ✅✅ dbh 异速 | 平 |
| **草本** | ✅ grass 简化 (无茎, 只 leaf:root) | ⚠️ 偏木本 cohort | **LPJ** |
| **作物** | ✅ cropallocation.cpp (DS 驱动) | ❌ 无 | **LPJ** |
| 全 PFT 通用 | ✅ group 继承 | ⚠️ 偏木本 | LPJ |
| COMPASS 已有基础 | ✅ R_m 已用 Sitch/LPJ, sfactor 现成 | ⚠️ 要新建 dbh+PARTEH | LPJ |

**决策: 选 LPJ-GUESS**。三决定性理由: (1) 全 PFT (树/草/作物) 全覆盖 — COMPASS 是陆面模型有三类; (2) 同时满足功能平衡+异速两条规律; (3) COMPASS R_m 已 Sitch 2003 同源, 一脉相承。
FATES 唯一强项是纯 dbh 异速更直接, 但代价是加 dbh 动态+PARTEH 整套且草/作物无解。

### §12.32.3 LPJ 三类 lifeform 分配 (解决树/草/作物)

| lifeform | 机制 | 器官 | 公式 (growth.cpp) |
|---|---|---|---|
| **tree** (苹果/林木) | 异速+功能平衡耦合 | leaf+root+sap/heart(茎) | 解 Eqn(13) 数值求 cmass_leaf_inc |
| **grass** (小草) | 只 leaf:root 平衡, **无茎** | leaf+root | `leaf_inc=(bminc−leaf/ltor+root)/(1+1/ltor)`; `root_inc=bminc−leaf_inc` |
| **crop** (作物) | DS 物候驱动 (类 FruitTree) | leaf+root+收获器官 | cropallocation.cpp |

tree 核心异速关系 (growth.cpp:399-410):
```
(2) 叶面积 = k_latosa × 边材横截面积       ← pipe model (茎随叶/树大)
(3) cmass_leaf = ltor × cmass_root         ← 功能平衡 (leaf:root 比, ltor 受水分胁迫调)
(4) height = k_allom2 × diam^k_allom3       ← 树高异速
→ 解 Eqn(13) f(cmass_leaf_inc)=0 数值求根 (二分法 NSEG=20/JMAX=40/XACC=0.0001)
```

### §12.32.4 全 PFT 固定参数 (LPJ-GUESS group 继承, 有出处)

LPJ-GUESS `data/ins/global.ins` group 继承 (tree/shrub/grass + broadleaved/needleleaved):

| PFT 类 | k_latosa | k_allom1 | k_allom2 | k_allom3 | wooddens | ltor_max |
|---|---|---|---|---|---|---|
| 阔叶树 (苹果 PFT23) | 6000 | 250 | 60 | 0.67 | 200 | 1.0 |
| 针叶树 | 5000 | 150 | 60 | 0.67 | 200 | 1.0 |
| 灌木 | — | 100 | 5 | 0.67 | 250 | 1.0 |
| 草 C3/C4 | — (无茎) | — | — | — | — | 0.5 |

(单位: k_latosa [m²/m²] 叶面积:边材面积; wooddens [kgC/m³]; ltor_max [-] leaf:root 上限)
源: global.ins group "tree"(L230) / "broadleaved"(L329) / "grass"(L297) / "needleleaved"(L342)。

### §12.32.5 接入 COMPASS 设计

1. **lifeform 三分类已现成** (Phase D 加的字段, 不需新增):
   ```matlab
   if crop_phenology == 1        → crop   (PFT 15-22 Corn/Wheat/...)
   elseif CN_target_stem <= 0    → grass  (PFT 12-14 C3/C4 grass, 无茎)
   else                          → tree   (PFT 1-11/23/24 木本, 有茎)
   ```
   验证: crop_phenology PFT15-22=1 其余=0; CN_target_stem PFT12-14=0 其余>0; is_perennial 木本/草=1 作物=0。
   三个 flag 组合即可干净区分, computeCarbonAllocation.m:81 已用 CN_target_stem<=0 判草本。
2. **PFT 表加** k_latosa/k_allom1-3/wooddens/ltor_max (group 默认值, 上表)。
3. **computeCarbonAllocation.m 按 lifeform 三分支**:
   - tree → 解 Eqn(13) (新函数 allocTreeLPJ, 二分法)
   - grass → leaf:root 解析式
   - crop → 保留现有 f_alloc 或 DS 驱动 (后续)
4. **ltor 受水分胁迫调** (功能平衡): COMPASS 已有 sfactor (ψ_leaf 木质部脆弱性), 缺水时 ltor↓ (多分根)。LPJ 用 wscal (水分胁迫标量) 调 ltor, COMPASS 可用 sfactor 等价。
5. **config 开关** `alloc_scheme = fixed / lpj_dynamic` (默认 fixed 保兼容)。
6. **守恒**: bminc = leaf_inc + root_inc + sap_inc 严格 (LPJ Eqn1), 总 C 不变, 机器精度。
7. **苹果果实**: tree 异速底座 + 叠加 fruit 分配 (FruitTree GDD 思路, 物候门 DOY 150-290)。

### §12.32.6 物种 vs 站点参数区分 (重要, 沿用 §2.3 原则)

- **异速系数** (k_latosa/k_allom/wooddens/ltor) = **物种/生活型属性** → 用 LPJ group 通用值 (全 PFT)
- **物候 GDD 阈值** = **站点气候** → 保持 COMPASS 长武标定, 不动
- **苹果 C:N** = 物种属性 → 已用 FruitTree 实测 (§12.31)

### §12.32.7 实现风险 + 分步

- **风险**: 动 C 分配核心 (computeCarbonAllocation + 涉及 24 PFT 分类 + 数值求根)。大重构。
- **分步建议**:
  - 步 1: PFT 表加 lifeform + LPJ 异速参数 (不改逻辑, 加字段)
  - 步 2: grass 分支 (最简单, leaf:root 解析) 先验证
  - 步 3: tree 分支 (Eqn13 数值求根) — 苹果验证
  - 步 4: ltor 水分胁迫耦合 (sfactor 驱动)
  - 每步单独跑全年验证守恒 + 不破坏现有 PFT
- **保留兼容**: config alloc_scheme=fixed 默认, 旧 f_alloc 不删。

*源码: LPJ-GUESS v4.1.1 @ D:/OneDrive/LPJ-GUESS-v4.1.1/guess_4.1/modules/growth.cpp:352-470 (allocation), global.ins (PFT 参数, group 继承)。状态: 设计完成, 待实现 (大重构, 分 4 步)。*

---

## §12.32.8 完整 LPJ 动态分配实施设计 (2026-05-29, 待实现)

**决策**: 不做简化版, 直接上完整 LPJ (异速 pipe model + 环境功能平衡), 全 PFT (树/草/灌木/作物)。

### A. LPJ 三问核实 (看实际源码, 非印象)

| 问 | 答 | 代码证据 |
|---|---|---|
| 是异速+环境? | ✅ 是 | growth.cpp:1199 `ltor = min(wscal_mean, nscal) × ltor_max` (水分+N胁迫调 leaf:root); pipe model Eq2 k_latosa + Eq4 height=k_allom2·diam^k_allom3 |
| 全 PFT 参数齐? | ✅ group 继承 | TeBS = tree+broadleaved+shade_tolerant+summergreen; 所有 PFT 经 group 拿全参数 |
| 草/灌木/作物? | ✅ 三套 | grass: 无茎 leaf:root (growth.cpp:673); shrub: tree lifeform 小系数; crop: 独立 cropallocation.cpp (DS 物候) |

### B. LPJ 核心方程 (tree, growth.cpp:396-470)

异速约束 (个体尺度, kgC/棵):
```
(1) bminc = leaf_inc + root_inc + sap_inc           守恒
(2) leaf_area = k_latosa × sapwood_xs_area          pipe model (茎随叶)
(3) cmass_leaf = ltor × cmass_root                   功能平衡 (ltor 受 wscal/nscal 调)
(4) height = k_allom2 × diam^k_allom3                树高异速
→ 联立消元得 Eqn(13) f(leaf_inc)=0, 二分法数值求根 (NSEG=20, JMAX=40, XACC=1e-4)
```
grass (无茎, growth.cpp:673): `leaf_inc = (bminc − leaf/ltor + root)/(1+1/ltor)`, `root_inc = bminc − leaf_inc`

### C. 全 PFT 参数表 (LPJ-GUESS group 继承, 有出处 global.ins)

| group | lifeform | ltor_max | k_latosa | k_allom1 | k_allom2 | k_allom3 | wooddens | crownarea_max |
|---|---|---|---|---|---|---|---|---|
| tree | tree | 1.0 | (broadleaved 6000 / needleleaved 5000) | (bl 250/nl 150) | 60 | 0.67 | 200 | 50 |
| shrub (低灌木) | tree | 1.0 | — | 100 | 5 | 0.67 | 250 | 10 |
| grass | grass | 0.5 | — (无茎) | — | — | — | — | — |

COMPASS 24 PFT 映射:
- PFT 1-9 乔木 + 23 苹果 + 24 Robinia → tree (broadleaved/needleleaved 按种)
- PFT 10-11 灌木 (BDS) → tree (shrub 系数)
- PFT 12-14 草 (C3/C4) → grass (ltor_max=0.5, 无茎)
- PFT 15-22 作物 → crop (保留现有 f_alloc 或 DS, 不走 LPJ tree)

### D. 需新建的状态 (完整 tree 必须, = 轮3 "拆 stem 双池")

| 新状态 | 含义 | 现状 |
|---|---|---|
| **C_sapwood + N_sapwood** | 边材 (活, 呼吸, pipe model) | C_stem_wood 单池拆分 |
| **C_heartwood + N_heartwood** | 心材 (死, 不呼吸, 结构) | 同上拆分 |
| **height** | 树高 (从 cmass_sap 异速反算) | 仅 canopy_height 常数 → 改 prognostic |
| **diam** | 胸径 (从 sap+heart 反算) | tree_DBH=4.2 常数 → 改动态 |
| (cmass_debt 碳债) | LPJ 借贷池 | 可选, 初期可省 (设 0) |

sap/heart 拆分后: f_livewood (v6) 直接变成 sap/(sap+heart) 比例, R_m_stem 只算 sap (CLM5 livestem 一致), 物理更干净。

### E. 改动文件清单 (14 文件涉及 stem 池, 拆双池连锁)

| 文件 | 改动 |
|---|---|
| getPFTParameterSet.m | 24 PFT 加 ltor_max/k_latosa/k_allom1-3/wooddens; C_stem_wood→C_sapwood+C_heartwood init |
| VegetationData.m | C_sapwood/heartwood/height/diam 状态字段 + 加载列表 |
| setInitialNitrogen.m | 池初始化 (sap/heart 分) |
| **computeCarbonAllocation.m** | 核心: lifeform 三分支 (tree 调 allocTreeLPJ / grass 解析 / crop 保留) |
| **+plant/+PlantCarbon/allocTreeLPJ.m** (新) | LPJ Eqn13 二分法求根 + height/diam 异速 |
| allocatePlantN.m | N 镜像拆 sap/heart |
| computeNPP.m | R_m_stem 只算 sapwood (f_livewood→sap/heart) |
| turnoverAndHarvest.m | stem turnover 分 sap→heart 转化 (Sitch Eq27 活木→死木) |
| runPlantCNStep.m / runPlantNDemand.m | 透传新池/参数 |
| registerPlantNitrogen.m + flushToCsv.m + COMPASS.m | 诊断输出 sap/heart/height/diam |
| config_vegetation.txt | alloc_scheme 开关 (fixed/lpj_dynamic) |

### F. 守恒保证

- C: `bminc = leaf_inc + root_inc + sap_inc (+ fruit_inc)` 严格 (LPJ Eqn1), 总 plant C 不变, 机器精度
- N: 镜像 C 分配, 按各器官 CN_target 换算, 总 plant N 守恒
- sap↔heart 转化 (turnover): 池间转移, 不影响总量

### G. 分步实施 (每步单独跑全年验证守恒)

| 步 | 内容 | 验证 |
|---|---|---|
| 0 | config alloc_scheme 开关 (默认 fixed, 兼容) + PFT 表加 LPJ 参数 (闲置) | 加载不崩 |
| 1 | **拆 stem 双池** C_sapwood/C_heartwood + N 镜像 + height/diam 动态状态 (fixed 分配下先跑通双池) | 守恒 + sap/heart 比合理 + R_m_stem 用 sap |
| 2 | grass 分支 (leaf:root 解析, 最简单) | C3/C4 grass PFT 守恒 |
| 3 | tree 分支 allocTreeLPJ (Eqn13 二分法) + height/diam 异速更新 | 苹果守恒 + leaf:root + 树高增长合理 |
| 4 | ltor 功能平衡 (wscal/nscal 调, COMPASS 用 sfactor ψ_leaf + N 胁迫) | 旱季根分配↑ |
| 5 | 苹果 fruit 分配叠加 (tree 底座 + 物候门 DOY 150-290) | 果实期分配对 |
| 6 | crop 分支 (作物 DS 驱动, 或保留现有 f_alloc) | 作物 PFT 守恒 |

**关键风险**: 步 1 (拆 stem 双池) 是连锁最广的 (14 文件), 且改 prognostic 状态结构, 必须先单独验证守恒不破再进步 3。height/diam 异速反算要防数值奇异 (幼树 diam→0)。

**物种 vs 站点**: 异速系数 (k_latosa/k_allom/wooddens/ltor_max) 用 LPJ group 通用值 (物种属性, 全 PFT); 物候 GDD 阈值保 COMPASS 长武标定 (站点属性, 不动)。

*源码: LPJ-GUESS growth.cpp:352-470 (allocation), :1199 (ltor 环境调), :673 (grass), cropallocation.cpp (作物); global.ins group 参数。状态: 完整设计完成, 待实现 (7 步大重构, 步1 拆 stem 双池连锁 14 文件)。*

### §12.32.9 LAI 路径决策 — 双引擎并存, 各走各的 (2026-05-29 关键决策)

**核心冲突**: LPJ Eq2 `leaf_area = k_latosa × sapwood_xs_area` (pipe model) 让叶面积受边材异速锁定 →
C_leaf 受边材约束 → LAI 由分配决定。这会**顶掉现有 EPIC LAI 经验曲线** (247 obs 校准, RMSE 0.348)。

**决策**: 用 config `LAI_engine` 开关让**两条 LAI 路径并存, 互不干扰**:

| LAI_engine | LAI 来源 | 分配方案 | 校准 | 状态 |
|---|---|---|---|---|
| **`epic_phen_mix`** (默认) | EPIC 经验曲线 (自由拟合) | 固定 f_alloc (或不含 Eq2 的简单分配) | 247 obs RMSE 0.348 | **保持原样, 不动** |
| **`leafc_sla_phen`** | C_leaf × SLA (碳驱动) | **完整 LPJ 异速 (含 Eq2 pipe model)** | 重做 (LAI 由异速决定) | **重构目标** |

**链条**:
- epic_phen_mix: 物候 + EPIC 曲线 → LAI; 碳分配独立 (不反馈 LAI)。LAI/C_leaf 双轨 (现状妥协, 保校准)。
- leafc_sla_phen: 分配 (Eqn13 受 Eq2 约束) → C_leaf → LAI = C_leaf×SLA。LAI/C_leaf 物理一致 (pipe model)。

**好处**:
1. EPIC 路径 (校准过的) 完全不受 LPJ 影响, RMSE 0.348 安全。
2. LPJ 完整异速 (含 Eq2) 在 leafc_sla_phen 路径单独实现/验证, 不破坏现有结果。
3. 可直接对比两路径 (同站点 epic vs lpj-physical)。

**实施约束**: §12.32.8 的 LPJ tree 分支 (allocTreeLPJ + Eq2 + height/diam 异速) **仅在 LAI_engine=leafc_sla_phen 时激活**;
epic_phen_mix 路径保持现有固定 f_alloc (computeCarbonAllocation 现状), 不引入 Eq2。
即: `LAI_engine` 同时决定 LAI 来源 + 是否走完整 LPJ 异速分配。

### §12.32.10 设计漏洞排查 + 决策 (2026-05-29)

完整 LPJ 设计认真挑出 5 个漏洞, 逐个决策:

| # | 漏洞 | 严重 | 决策 |
|---|---|---|---|
| 1 | **单位尺度**: LPJ Eq2/Eq4 pipe model 是个体尺度 (m²/棵, kgC/棵), COMPASS 碳氮是面积 (gC/m²地面) | 🔴 | 自然树用 tree_CAH (每棵冠面积 4.13 m²) + tree_DBH 动态做个体↔面积转换; **苹果改走物候分配避开 pipe model (见#4)** |
| 2 | leafc_sla_phen 引擎 LAI 用 N_leaf×CN_L 反推 C_leaf, 非直接 C_leaf 池 | 🔴 | **改直接用 C_leaf 池** (CN 本就是 C↔N 转化, 直接用 C_leaf 更准, LPJ 分配产出 C_leaf 直接喂 LAI) |
| 3 | height/diam 异速反算数值奇异 (幼树 diam→0 除零/虚根) | 🟡 | **移植 LPJ clamp + cmass_debt 借贷保护** (growth.cpp NSEG/JMAX/XACC/CDEBT_MAXLOAN) |
| 4 | 苹果矮化密植 ≠ 通用乔木异速 (LPJ k_latosa/k_allom 针对自然林大树) | 🟡 | **苹果不走 LPJ 乔木异速, 改走 FruitTree 物候动态分配 (类作物)**; 自然树才用 LPJ pipe model |
| 5 | crop 分支未定 | 🟢 | 作物用 DS/物候驱动 (FruitTree Eq4-7 / cropallocation 思路) |

**关键转折 (#4 决策)**: 苹果归"物候动态分配"路线 (与作物同), 不走 LPJ 乔木 pipe model。
好处: (a) 避开苹果矮化不适用通用乔木异速 (用 FruitTree Zanotelli 苹果实测更准); (b) 避开漏洞 1 单位尺度 (FruitTree 是面积尺度物候分配, 不用个体 pipe model); (c) 苹果是当前唯一验证站点, FruitTree 同物种最对口。

### §12.32.11 最终架构 — 三套分配按 lifeform/物候分流

```
分流判据 (现有 PFT 表字段):
  crop_phenology==1 (PFT 15-22 作物) ──┐
  PFT 23 苹果 (has_fruit, 矮化果园)   ──┴→ [物候动态分配] FruitTree Eq4-7 (GDD/DS 驱动, 面积尺度)
                                            aleaf 季初高指数衰减, arepr 果实膨大增长 plateau
  CN_target_stem<=0 (PFT 12-14 草)   ────→ [LPJ grass] leaf:root 按 ltor (无茎, 面积尺度直接)
  其余 (PFT 1-11/24 自然树/灌木)      ────→ [LPJ tree] pipe model Eq2/4 + ltor 功能平衡
                                            (个体尺度, 需 tree_CAH/DBH 动态转换 + sap/heart 双池)
```

| 路线 | PFT | 机制 | 单位 | 需新状态 | 参数源 |
|---|---|---|---|---|---|
| 物候动态 | 23 苹果 + 15-22 作物 | FruitTree GDD 阶段变系数 | 面积 gC/m² | 无 (复用现有 C 池) | FruitTree Table C1 (苹果) + crop |
| LPJ grass | 12-14 草 | leaf:root (ltor) | 面积 gC/m² | 无 | LPJ grass group (ltor_max=0.5) |
| LPJ tree | 1-11/24 树灌 | pipe model + ltor | 个体↔面积 | C_sap/C_heart 双池 + height + DBH 动态 | LPJ tree group (k_latosa/k_allom/wooddens/ltor_max=1) |

**实施优先级** (苹果是验证站点, 先做苹果路线):
1. **物候动态分配 (苹果)** — 先做, FruitTree Eq4-7, 面积尺度无新状态, 验证苹果。这是当前最高 ROI。
2. LPJ grass — 简单, 验证草本 PFT。
3. LPJ tree (pipe model + 拆 sap/heart 双池 + DBH 动态) — 最大工程, 自然树用, 后做。
4. crop — 作物 DS 驱动。

**注意**: §12.32.8 的"拆 stem 双池 + DBH 动态"**只自然树 (LPJ tree) 需要**, 苹果走物候分配**不需要** (复用现有单 C_stem_wood + f_alloc 风格但物候动态系数)。这大幅降低苹果路线工作量。

**LAI 路径** (§12.32.9): 苹果物候动态分配产出 C_leaf → leafc_sla_phen 引擎 LAI=C_leaf×SLA (#2 改直接用 C_leaf)。epic_phen_mix 仍保留 (EPIC 校准不动)。

---

## §12.33 全景梳理 + 逻辑图 + 待办合并 (2026-05-29 整理)

> 本节是 §12.26-§12.32 的统一总览。前面各 § 是细节归档, 本节合并待办 + 物理过程时间线 + 逻辑图, **以本节为最新状态入口**。

### §12.33.1 植物 C/N + 光合 物理过程时间线 (v5→v10)

| 版本 | 改动 | 出处 | 验证结果 | § |
|---|---|---|---|---|
| v5 | LUNA strict path 接通 (Jmax 走 LUNA 独立, ebal 透传 use_luna/Jmax25_LUNA) | CLM5 PhotosynthesisMod.F90:1515-1517 | Vcmax 峰 107.93 达标 | §12.30 前序 |
| (标定) | Jmaxb1 0.1745→0.244 (LUNA Vcmax→105) | Ali 2016 范围内, 反推 | Vcmax 年均 82 | — |
| v6 | R_m_stem sapwood 修正 (f_livewood=0.22, 只活木呼吸) | CLM5 CNMRespMod.F90:265 + Sitch 2003 | R_a/GPP 109%→78% | §12.30.2-4 |
| v7 | R_m_leaf 光抑制 (白天 Rd×0.67) | CLM5 PhotosynthesisMod.F90:1562 (Lloyd 2010) | R_a/GPP 78%→71%, 白天≈夜间 | §12.30.7 |
| (改名) | starvation→storage_drawdown (非致死澄清) | FATES EDMortalityFunctionsMod.F90:201 对比 | — | §12.31.1 |
| v8/v9 | cushion 储备回补 (FATES 1.2) — **失败** | FATES — | 抽水泵 bug, C_labile 砸0 32%→88%, **已删** | §12.31.3 → v10 删 |
| v10 | C:N 对齐 FruitTree (leafCN 19.7/frootCN 32/livewdCN 60/grainCN 139) + gs_g1 8.2 + C 分配诊断 (cushion 删) | FruitTree Table C1 (Zanotelli 实测) | 跑中验证 | §12.31.3 + §12.32 序 |

**累积成果**: R_a/GPP 从不合理 109% (碳收支负) → 健康果园 71%; C/N 比全用苹果实测; 储备增厚; 守恒机器精度。

### §12.33.2 待办合并 (去重 + 真实状态, 2026-05-29)

旧 § 里的 M2-B.4/B.5/B.8 等标记多已完成 (现已到 v10), 下表是**当前真实未完成项**:

| 优先 | 待办 | 出处/方案 | 风险 | 依赖 |
|---|---|---|---|---|
| 🔴🔴 最高 | **C/N 土壤-植被耦合修复 F1-F5** (生长呼吸2.5×/C_leaf→LAI/D_N驱动吸收/litter C:N=60/细根门控) — **先于物候, 耦合无误才做动态分配** | **§12.34** | 低-高 (逐条) | — |
| ✅ 完成 | v12 FlexibleCN + C主导init + 方案A(fcur储备) | §12.34 序 | — | smoke验证中 |
| 🔴 高 | **动态分配实施** (苹果物候 FruitTree Eq4-7) — **C/N 耦合修复(§12.34)之后** | cn_coupling §12.32.11 | 中 (苹果路线无新状态) | §12.34 完 |
| 🟡 中 | slatop 0.0125 vs FruitTree 0.028 核查 | 查 Wang 2009 原文 | 中 (动 LAI 校准 RMSE 0.348) | 不盲改 |
| 🟡 中 | C_labile 果实膨大期砸 0 (v7 32%) | 动态分配里作分配优先级解决 (非事后 cushion) | — | 动态分配 |
| 🟢 低 | LPJ tree pipe model + 拆 sap/heart 双池 + DBH 动态 (自然树用) | §12.32.11 步3 | 高 (连锁 14 文件) | 苹果路线后 |
| 🟢 低 | lnc_opt=.true. (土壤 N 真影响光合, lnc 用 N_leaf/LAI) | Climate10dBuffer.m:234 | 中 (需 N 池稳) | N 池稳后 |
| 🟢 低 | energy_balance_residual 诊断字段 masked | ebal 主循环加 reg.write | 低 | — |
| 🟢 低 | 水残差定位 (v3 9.2mm; v5+ 已降 0.007mm) | 加 bottom_drainage 累积 | 低 | 多半已解 |
| 🟢 低 | Ghimire 2016 / Poorter 2012 下载 | escholarship 限流 / Wiley 付费 | — | 次要 |
| 🟢 低 | 修剪模块 (prune_fr=0.85) | FruitTree 果园管理 | 中 | 未来 |

### §12.33.3 逻辑图 — 植物 C/N/光合/分配 全链 (v10 态)

```
   ┌──────────────────────── 气象 forcing (T, PAR Q, CO2 Ca, VPD, 降水) ────────────────────────┐
   │                                                                                            │
   ▼                                                                                            │
┌─────────┐  lnc=1/(slatop·leafcn)  ┌──────────────┐                                            │
│LUNA 每日 │ ──(lnc_opt=.false.)───▶ │Vcmax25/Jmax25│                                            │
│Ali 2016 │ ◀╌╌╌ N 反馈 (待 lnc_opt) │  (v5 接通)    │                                            │
│Jmaxb1   │       ╎                  └──────┬───────┘                                            │
└─────────┘       ╎                         ▼                                                   │
                  ╎          ┌──────────────────────────────────┐                               │
                  ╎          │  biochemical.m (Farquhar-Collatz) │◀──── PAR/CO2/VPD ─────────────┘
   N_leaf 池 ─────╯          │  Ac/Aj/Ap 共限制 → A, Ag          │
   (光合容量)                │  Rd×0.67 白天光抑制 (v7)          │◀══ ψ_leaf 气孔限制 (水力胁迫)
                             │  Medlyn gs (g1=8.2 v10)           │            ▲
                             └──────┬───────────────┬───────────┘            ║ (水力上行)
                              GPP   │            Rd_leaf │  蒸腾 T ═══════════╝
                                    ▼                │   (gs 驱动失水)        ║
                  ┌────────────────────────────┐    │                        ║
                  │ computeNPP (Sitch 2003)     │    │            ┌───────────────────────┐
                  │  R_m,leaf = Rd_leaf ────────┼────┘            │ 根系水力 (Plant Hydraul│
                  │  R_m,root = r·N_root·g(T)   │                 │ + calSPAC/V2)          │
                  │  R_m,stem = r·N_stem·       │ f_livewood      │  ψ_soil→ψ_root→ψ_stem  │
                  │            f_livewood·g(T)  │ =0.22(v6)       │  →ψ_leaf (gsr 串联)    │
                  │  R_g=0.25(GPP-R_m)          │                 │  根吸水 → 蒸腾流 ═══════╝
                  │  NPP=GPP-R_m-R_g            │                 └───────────▲───────────┘
                  └──────────────┬─────────────┘                             ║ 土壤水 ψ_soil
                             NPP │ (=bminc)                                   ║ (RWU 根水吸取)
                                 ▼                                  ┌─────────╨──────────┐
   ┌──────────────────────────────────────────────┐                │ 土壤 (STEMMUS 水热) │
   │ applyNSCCoupling (NSC cap + overdraft)        │                │  hh/TT/Theta_LL     │
   │  R_a cap=min(R_m_pot,GPP+NSC) [curmr 等价]    │                └─────────┬──────────┘
   │  C_labile<0 → C_storage overdraft 借          │                          ║ 土壤无机 N
   │  (cushion 已删 v10)                           │                ┌─────────╨──────────┐
   └──────────────┬────────────────────────────────┘                │ 土壤 N (STEMMUS_N)  │
        NPP_actual│                                                  │ NH4/NO3 + 硝化/反硝 │
                  ▼                                                  │ + 施肥/淋溶         │
   ┌──────────────────────────────────────────────┐                └─────────┬──────────┘
   │ computeCarbonAllocation 分配 (现固定 f_alloc) │                          ║ 根吸 N (U_N)
   │ 【待重构 §12.32.11 三套分流】:                │                          ║ DAISY demand + Nacry
   │  苹果/作物→FruitTree 物候动态 (aleaf衰减/arepr)│                ┌─────────╨──────────┐
   │  草→LPJ grass (leaf:root,ltor)                │                │ runPlantNDemand     │
   │  自然树→LPJ tree (pipe Eq2/4+ltor 功能平衡)   │  ltor 受 ψ 调  │ +applyRootUptake    │
   └──┬──────────┬──────────┬──────────┬───────────┘  (缺水多分根) └─────────┬──────────┘
      ▼          ▼          ▼          ▼                                      ║ U_N
   C_leaf   C_root_fine C_stem_wood C_fruit                                   ║
      │          │          │          │     ┌──────────────────────────┐    ║
      │          │          │          └────▶│ allocatePlantN: 根吸 U_N  │◀═══╝
      ▼(LAI)     │          │(拆sap/heart待) │ → 6 N 池 (labile/leaf/    │
┌──────────┐    │   ┌─────────────────────┐ │   root/stem/fruit/storage)│
│LAI 双引擎 │    │   │turnoverAndHarvest    │ └────────────┬─────────────┘
│epic(校准) │    └──▶│ 落叶/根周转/采收/春萌 │              │ N_leaf 池
│leafc_sla  │◀───────│turnoverCarbonMirror  │              ╰──────────────▶ (回 LUNA/光合容量, 闭环)
│(C_leaf×SLA│        │ (N→C, CN_target)     │
└──────────┘        └──────────┬───────────┘
                               ▼ litter → 土壤 SOM (师弟 M3, 接 litter→矿化→无机N, 闭环)
                          回土壤 N 池 (释放无机 N 供再吸收)

══ 双线 = 水/N 上行闭环 (根→冠层); ── 单线 = 碳下行 (光合→分配→池)
闭环1 (水): 土壤水→根吸→ψ_leaf→气孔→蒸腾→失水→拉低ψ→再吸  (水力串联 gsr)
闭环2 (N): 土壤无机N→根吸U_N→N_leaf→光合容量(Vcmax)→GPP→分配→落叶litter→矿化→土壤N
闭环3 (C): GPP→NPP→各器官C池→turnover litter→SOM→(R_h)→大气, NEE=R_a+R_h−GPP
物种参数借文献(异速/C:N/光合 CLM5/FruitTree/LPJ); 站点参数保标定(GDD物候 长武; LAI Wang247obs)
```

### §12.33.4 关键设计原则 (贯穿全链)

1. **物种属性借文献, 站点属性保标定**: C:N/Vcmax/异速 → CLM5/FruitTree/LPJ; GDD 物候/LAI 校准 → 长武 25yr/Wang 247obs。
2. **守恒第一**: 每个通量进 cum_* 累积, ΔΣpool ≡ Σin−Σout 机器精度。
3. **§5.5 反幻觉**: 每个参数引用具体到文件:行号 + 本地 PDF, 无出处标 [site-calibrated]。
4. **不破不立**: 新方案走 config 开关 (LAI_engine/alloc_scheme), 旧路径保留对比。
5. **苹果是当前唯一验证站点**: 优先苹果路线 (FruitTree 同物种最对口), 自然树通用方案后做。

### §12.33.5 统一路线图 → **已删, 见文档开头「📍 阶段路线图」(全文唯一权威, 含 12.5 C/N 耦合修复)**。
**辅助 (非主线, 非路线图)**: energy_balance 诊断字段接通 / 水残差定位 / slatop 核查 / 文献补下 (见 §12.33.2)。

---

## §12.34 C/N 土壤-植被耦合 物理过程审计 + 修复报告 (2026-05-30)

> **多模型对比审计**: 4 并行 agent 对比 CLM5 / FATES / DAISY / LPJ-GUESS / JULES-ClimaLand 源码 (均在 OneDrive)。
> **原则 (用户决策 2026-05-30)**: **先确保 C/N 土壤-植被耦合物理无误, 再做物候**。下列 F1-F5 修复**优先级高于 stage 13 物候**——耦合物理不对, 物候动态分配是建在错地基上。
> **背景**: 长武苹果园**过量 N 肥 → N 不限制生长 (f_N 二值=NPP有无, 非N胁迫)**, 故 C:N 锁 target 正确、FlexibleCN 对苹果休眠 (留作 N 限制 PFT 地基)。v12 已实现 FlexibleCN(带内浮动) + 方案A(fcur 越冬储备分配)。

### §12.34.1 高优先修复 (已亲核代码)

**F1 — 生长呼吸 R_g 偏高 ~2.5× (产量低 3× 的候选主因)**  〔严重度 高 / 核对 ✓ / 风险 低〕
- 现状: `computeNPP.m:122` `R_g = R_g_fraction(0.25) × max(GPP−R_m,0)` (Sitch 2003 "25% of remainder")。
- 参照: CLM5 `CNGRespMod.F90` `grperc=0.11` 作用在**新生长碳** `cpool_to_Xc` (等效 0.099 of (GPP−R_m)); FATES `EDPhysiologyMod.F90:800` 同式 `cost×(grperc+1)`。→ COMPASS 0.25 vs CLM5 0.099, **生长呼吸多算 2.52×** → 抬 R_a、压 NPP → 直接拉低产量。
- 修复步骤:
  1. PFT 表 PFT23 `R_g_fraction` 0.25 → **0.11** (CLM5 grperc), 注释引 `CNGRespMod.F90`; 保留 parallel 对比。
  2. (可选更严) 改 R_g 作用对象为分配的新生长碳而非 (GPP−R_m): `R_g = grperc × C_growth`。先只改系数验证产量响应, 再决定是否重排。
  3. 全年跑: 产量 (期望向 49315 kg/ha 靠) + R_a/GPP (期望降但仍合理) + 守恒。

**F2 — C_leaf 接通 LAI (让 C/N 真反馈 GPP; N→Vcmax 闭环前提)**  〔严重度 高 / 核对 ✓ / 风险 高〕
- 现状: LAI 走 EPIC 热单位 S 曲线 (`computeLAIfromPhen.m` `LAI_engine=epic_phen_mix`), `C_leaf` 不被读; `lnc` 冻成 PFT 常数 (`Climate10dBuffer.m:234`)。→ **C 主导/FlexibleCN/C:N 对 GPP 零影响 (纯账本层)**, 这解释了 v12 改分配 GPP/LAI 不变。
- 参照: CLM5 `CNVegStructUpdateMod.F90:152` `tlai=f(leafc·slatop)`; LPJ `lai=cmass_leaf×SLA`。
- 修复步骤:
  1. 切 `LAI_engine=leafc_sla_phen` (已存在但 inactive, `computeLAIfromPhen.m:59-76`): `LAI=C_leaf×SLA×phen`。
  2. ⚠️ EPIC LAI 是 Wang 247obs 校准 (RMSE 0.348, 峰 2.79)。切 leafc_sla 会动 LAI 校准基线 → 先**双引擎并跑对比** C_leaf×SLA 曲线 vs EPIC, 量级/峰值对得上再切, 否则 LAI/GPP 崩。
  3. 接通后同步开 `lnc_opt`=真 (lnc=N_leaf/LAI 喂 LUNA, stage 15) → N→Vcmax→GPP 闭环。
- 注: 原属 stage 14/15, 因是"C/N 反馈 GPP"前提, 提到耦合修复。高风险, 排 F1 之后。

##### 🔬 F2 深挖诊断 + f_alloc 标定 (2026-06-01)

**架构真相 (核实)**: 主循环当前**没调** `computeLAIfromPhen`。`LAI_mode='model'` 实际走 `COMPASS.m:550 canopy.LAI=cfg.Veg.LAI(KT)` (VegetationData 预计算的 EPIC 曲线), 一个 hotfix (`COMPASS.m:538-549`) 为解决"中途起跑 phenology 累加器未 spin-up"把在线 LAI 路径整个绕过。
→ 改 config `LAI_engine=leafc_sla_phen` **当前无效**。**矛盾本质**: EPIC 只依赖气象 forcing 可预计算(offline); leafc_sla 的 C_leaf 是主循环逐步算的(online, 闭环 GPP→C_leaf→LAI→GPP), 没法预计算。F2 真正要做的是"把 LAI 从预计算搬到在线", 是架构改动, 这才是它高风险的真相。

**离线双引擎对比 (用 v12_cdom 全年 C_leaf, 零风险诊断, tools/diag_f2_lai_engine_compare.m)**:
| 指标 | EPIC(校准基线) | leafc_sla(碳驱动) | 判读 |
|---|---|---|---|
| 峰值 LAI | 2.79 | **6.67** | ❌ 高 2.39× |
| 峰值 C_leaf | — | 403.2 gC/m² (要 195 才匹配 2.79) | C_leaf 高 ~2× |
| 峰值时相 | KT6313(~9月初) | KT7025(~10月下, 晚 30d) | ⚠️ 落叶晚 |
| 相关 r | — | 0.852 | ✅ 形状对, 量级崩 |

**根因 (核实, tools/diag_f2_cleaf_budget.m)**: 苹果 f_alloc 错配 + 果窗外倒给叶, 三处推高叶碳:
1. 基础 f_alloc_leaf=0.40 本就偏高 (实跑年均 0.43)
2. f_alloc_fruit=0.15 偏低, 且只 DOY150-290 窗口 (年均果占比仅 ~0.058)
3. `computeOrganDemand.m:74-76` 果窗外把果份额 **60/40 倒给叶/根** → 果窗外叶分配飙到 0.49 (成龄树非果期叶子不该长)
→ C_leaf 偏高 2× (F2 LAI 崩) + C_fruit 偏低 (产量低 3×, 对账 yield_benchmark 16490 vs 49315) + 碳氮整体偏多 (连带 N_labile 虚高), **三症同根**。

**实测标定靶 (亲验原文, PDF 在 data/references/fruit_tree/)**:
- **Zanotelli et al. 2013 Biogeosciences 10:3089** (Fuji/M9 成龄, 总 NPP 960 gC/m²/yr): **fruit 52%** (p.3098, 495±35 gC) / **leaf 11%** (p.3097, 106±5) / **woody AG 18%+BG 1%** (p.3097, 173+13) / fine-root standing ~1% (大部分 124 shed 由周转处理)。abstract: "largely allocated to fruit: 49%"。
- **Wu et al. 2012 PLoS ONE 7(6):e38883** (北京 18龄, 对口长武 18龄乔砧, 交叉验证): fruit 46/woody 37/leaf 14/fine-root 3%。
- (ScienceDirect NECB apple 付费墙未下, 盗版不碰; Zanotelli 同组同园开放版已覆盖。)

**路 2 改动 (2026-06-01, checkcode CLEAN)**:
- `getPFTParameterSet.m` PFT23: leaf 0.40→**0.12** / fruit 0.15→**0.50** / stem 0.20→**0.30** / root 0.25→**0.08** (出处 Zanotelli 2013 逐项注释)
- `computeOrganDemand.m:73-79`: 果窗外不再倒给叶, 改 **倒给 stem** (成龄树非果期碳进木质/储备, Zanotelli p.3097)
- ⚠️ 量纲注: f_alloc_fruit=0.50 是果窗内瞬时, DOY 门控会摊薄年均 (用户定: 果不细抠凑年均, 先看 C_leaf 降/产量升)
- 验证: 全年跑 output_v19_falloc_cal, 看 C_leaf 是否降到 ~195 / 产量 / R_a-GPP / 守恒。**通过后 leafc_sla 才可能切 (F2-③ 架构改造)**。

**F3 — 根吸 N 需求改 need-based (生长 D_N + 储备补到目标即停)**  〔严重度 高 / 核对 ✓实测 / 风险 中 / **已实现 2026-05-31**〕
- 现状(旧): `runPlantNDemand.m` 吸收需求 = `N_target_per_LAI×LAI − N_labile` (补 labile 桶, N_target_per_LAI=2.0 拍值); C 定义的 `D_N` 算了没喂 uptake。
- 决策历程 (反复对标, 见下):
  1. v1 `D_N + LAI兜底`: 守恒✓但 N_labile 5.5→21.8 (兜底目标 N_target≈10 太大)。
  2. v2 `D_N×kstore` (LPJmL von Bloh 2018 树1.15): 与 v1 结果**完全相同** (N_labile 仍 21.8) — 因 kstore 恒定放大, 不看 labile 满没满, 一直多吸。
  3. **实测核对** (fruit_tree/ 下 PDF): Cheng 2021 Frontiers (矮砧整株 N 4-7.6 g/m²) + Wang 2024 Agronomy 14:963 (长武实测 DBH 20-22cm/树高3.1m 乔砧大树)。结论: 长武乔砧大树 C_stem=1121 量级合理 (不虚高); N_labile=21.8 是 demand 公式推高 (与初值无关)。
  4. **三模型机制核对**: CLM5 npool 同样虚高 (nscalar 反馈弱); LPJ need-based (器官饱和 demand=0) + max_n_storage; DAISY 无独立 labile 池 (组织浓度封顶)。→ 学 LPJ need-based。
- **最终方案 (need-based, 对标 LPJ canexch.cpp:1163 + growth.cpp:1479; Smith et al. 2014 GMD 7:2261)**:
  ```matlab
  % runPlantNDemand.m:119+ (已实现)
  N_labile_target = (C_leaf/CN_L + C_root/CN_target_root) * nrelocfrac;  % 储备目标 = 重建明年一茬叶+细根 N (LPJ retransn_nextyear)
  storage_deficit = max(0, N_labile_target − N_labile_prev);            % labile 满了 = 0 (自停)
  demand = D_N.total + storage_deficit / Delt_t;                         % 生长 + 储备补满即停
  ```
  - 苹果算: 储备目标 ≈ (89.56/19.7 + 150.5/32)×0.7 = **6.47 g N/m²** → N_labile 预期稳 ~6.5 (对得上 F1=5.5 + 实测量级), 不再涨到 21.8。
  - **零新参数**: C_leaf/C_root (NitrogenVariables 时序), CN_L/CN_target_root, nrelocfrac=0.7 (PFT 表) 全现有。
  - 单位 g N/(m²·s): D_N.total (computeOrganDemand.m:89) + storage_deficit[g N/m²]/dt; 下游 applyRootUptake 内部转 kg N/(m-root·s) (不碰)。
  - **舍弃边材容量项** `C_sap×fnstorage/CN`: 需 sap/heart 双池 + 新参 fnstorage; 苹果该项(用 f_livewood 压缩态算)3.76 < 6.47 萌发项, max 取后者, 故现不需要。待 13c 拆双池补 (见 §12.30.3 [F3 关联])。
  - **全 PFT 通用 (三类植被的 N_labile 都改了)**: demand 公式全 PFT 共用一条 (runPlantNDemand 所有 PFT 走), 储备目标用的量 (C_leaf/C_root 时序, CN_L/CN_target_root, nrelocfrac) **三类都有**: 草(13-14 nreloc=0.5/CN_root=42)、作物(15-22 nreloc=0.5)、树(23 nreloc=0.7/CN_root=32, 24)。
    ⚠️ **但仅苹果(PFT23)实测验证过**; 草/作物的 nrelocfrac=0.5 是默认值非站点标定, 储备目标量级**待对应 PFT 验证点核** (同 GL/PH 验证驱动策略, 见阶段16)。
- 风险 (中): **必跑 smoke 验**: PlantN_residual 机器精度 + N_labile 回落到 ~6.5 (不再 21.8) + 果实期 U_N。详见记忆 [[project_compass_f3_demand_dual]]。

**F4 — litter C:N 用凋落物富集值 (而非活叶 C:N) [接师弟 SOM 时]**  〔严重度 高 / 核对 ✓+实测源 / 风险 低〕
- 现状: `exportLitterToSOM.m:41-42` `CN_leaf_litter=CN_L(19.7)`、`CN_root_litter=CN_target_root`。→ 喂 SOM 偏富氮 litter, 高估矿化——偏的正是 COMPASS 要闭合的土壤 N 供给。
- 参照: CLM5 单独 `lflitcn` (N 回收后更高, FUN 限 1.5-3× leafcn, `CNFUNMod.F90:1730`)。**FruitTree Table C1 `lflitCN=60` (Zanotelli 2010 实测)** vs 活叶 `leafCN=19.7` → **差 3×, 印证 N 回收富集**。
- 修复步骤:
  1. PFT 表加 `lflitcn` (苹果=**60**, FruitTree Zanotelli 实测, 标出处)。
  2. `exportLitterToSOM`: `CN_leaf_litter=lflitcn(60)`。
  3. 同步 **decouple 叶 C/N 回收分数** (现共用 retransl_fraction=0.7 → litter C:N=活叶 19.7): N 回收分数调到使 litter C:N→60 (N 回收多于 C); C 回收→litter (见 §12.34 注: 也解 v12 叶C回收 stopgap)。

**F5 — 细根周转加季节/温度门控**  〔严重度 中-高 / 核对 ✓ / 风险 中〕
- 现状: `fineRootTurnover.m:27` 恒 1/yr、无门控, **冬天/休眠也掉根** (0.27%/天)。→ 地下 litter N 时序错 (本该秋季脉冲)。
- 参照: CLM5 细根脱落同步落叶/温度 (`CNPhenologyMod.F90:2453,2627` offset-ramp + bglfr·lgsf); FATES `fates_phen_fnrt_drop_frac` 冷/旱触发。
- 修复步骤: `fineRootTurnover` 加温度/物候门 (休眠期/低温降速, 秋季同步落叶脉冲)。部分沾物候, 排 F 末。

### §12.34.2 中/低优先 (agent 标记, 实施时再亲核)
| 项 | 现状 vs 参照 | 严重 |
|---|---|---|
| 茎分配与 NPP 无关 | CLM5 低产时 f_stem→0 保护叶根 (FlexibleCN:424); COMPASS 恒 0.2 | 中 |
| 干旱无 β 下调 R_m | JULES `Rpm×β` (autotrophic_respiration.jl:240); COMPASS 半干旱夏季照付全额 → 加剧 NSC 抽干 | 中 |
| 无活木→死木 周转 + N 回流 | CLM5 `CNLivewoodTurnover` livestem→deadstem + livestemn_to_retransn + CWD litter; COMPASS 茎 N 永久锁 | 中 |
| 回流 N 进 storage 春季才放 | CLM5 retransn→npool 优先花掉抑制吸收 (FlexibleCN:1742); COMPASS 夏秋不能用、过度依赖土壤吸收 | 中 |
| 无奢侈吸收/存 reserve N | DAISY luxury uptake (crpn.C:67); COMPASS labile 卡 LAI 目标, 过量施肥站点存不了 reserve N | 中 |
| cn_flex_offset 绝对±10 | 苹果低 C:N 19.7 → ±51% 带 (CLM5 ±10 on 高 base = ±17-40%); 应取相对带 | 中(现muted) |
| NSC cap 丢未付维持呼吸 | CLM5 curmr+xsmrpool 记债 (FlexibleCN:1419); COMPASS 直接 cap 丢 | 中 |
| 根 R_m 用全柱均温 | CLM5 root-fraction 加权 tcsoi (CNMRespMod:289); COMPASS `mean(TT)` 拉入深冷层 | 低 |
| 细根周转率全 PFT 1/yr | FATES 1-2yr 分 PFT; 苹果 OK 其余偏快 | 低 |
| 果实 gate 60/40 硬分 / leaf_long_days 死参 / 常绿零落叶 / M-M demand 模式旁路 | — | 低 |

### §12.34.3 优先级 (明确提到物候之上)
**阶段 12.5 (C/N 耦合修复 F1-F5) 先于 阶段 13 (物候动态分配)。** 顺序建议: **F1(生长呼吸, 低风险高收益验产量) → F3(D_N 驱动吸收) → F4(litter C:N, 接 SOM 时) → F5(细根门控) → F2(C_leaf→LAI, 高风险双引擎对比后切)**。F2 接通后 N→Vcmax 闭环, 才真正"C/N 耦合无误", 然后进物候。中/低项按需穿插。

---

## §12.35 物候完整性对标 (CLM5 + LPJ-GUESS + DAISY, 2026-05-30)

三模型逐生长过程对标 §B.9 物候相位机。**苹果 (season_decid): 核心等价** (Sykes chilling-GDD + 日长落叶 + aphen
+ 冷月 reset, 三方确认结构等价, 甚至比 LPJ 原版多 onset_gddflag 门控 + 固定冷月 reset 防误触)。**通用 LSM
有 5 个真缺口**, 全部"对应 PFT 验证点激活时才暴露", 按验证驱动补 (同 GL1-GL4 池缺口策略)。

| # | 物候缺口 | 触发 PFT | 现状 | 出处 (对标行号) |
|---|---|---|---|---|
| PH1 | **crop 物候占位** (缺出苗 emergence/DS 多阶段/春化/光周期, 实现度<30%) | 作物 15-22 | 占位 | DAISY phenology_TSum.C:94-103(出苗)/vernalization.C:57-136; CLM5 CropPhenology:1420-2031 |
| PH2 | **stress_decid 占位** (4 累积器 onset_fdd/gdd/swi 有 docstring 无实现) | 草13,14/雨绿6/灌12 | 占位 | CLM5 CNStressDecidPhenology:974-1417 |
| PH3 | **background litterfall / 叶寿命** (常绿连续落叶 bglfr=1/leaf_long) | 常绿 1,2,4,5,9 | 简化(落叶靠 turnover 模块, 物候解耦) | CLM5 CNEvergreenPhenology:603-605 |
| PH4 | **phenology 混合型 ANY** (同 PFT 既 GDD 萌发又旱落叶; dispatcher 强制 4 选 1 建不了) | 热带季雨混合林 | 缺失 | LPJ growth.cpp:77-116 |
| PH5 | **storage→xfer→display 三阶段展叶** (CLM5 fstor2tran 缓冲春萌储备碳释放) | 所有落叶(含苹果) | 简化(线性 phen ramp, C/N 转移在 turnover 独立做) | CLM5 :923-941 |

**苹果当前验证点**: PH5 影响苹果但 Bar B 不阻塞 (功能未塌, 仅耦合度低); PH1-PH4 苹果不沾。**苹果物候 Bar B 已够。**
**非物候 (剥离)**: 探子提的 establishment/mortality/树龄/DBH = 植被动态 (vegdynam, 另议); allocation solver/ltor/果实多阶段
TSum/展叶期优先分配 = 动态分配 (阶段 13, 见 §F), 不算物候缺口。

> **全局修复优先级 → 见文档开头「📍 阶段路线图」的「三层主线」** (统一 F/GL/PH/N/arch, 全文唯一权威)。
> 此处不再单列, 防与开头分散打架。物候缺口 PH1-PH5 已并入开头阶段 16 (按验证点驱动) + 本节明细。

> **一句话**: 现在只做**第 1 层 F1→F3→F5→F2** (苹果 Bar B); 做完进第 2 层 (土壤 litter + 物候动态);
> 第 3 层全部等森林/草/作物**真实验证点**到位再按需补 — 不为"通用"提前堆砌无法验证的结构。

---

## §12.36 Fluspect 9 叶片参数 — 过程模拟可行性评估 (2026-06-02)

**背景**: Fluspect (叶片辐射传输, SCOPE 子模块) 需 9 个叶光学参数算反射/透射/SIF。COMPASS 默认从 SCOPE Excel 喂**固定值** (Cab=100 等)。问题: 这 9 个能否改由 COMPASS **过程模型动态算** (而非固定输入)? 这是 COMPASS"碳氮水状态→叶光学→冠层遥感(反射/SIF/VOD)"端到端过程驱动正演的核心 — 别的陆面模型喂固定叶参或经验 LAI 代理, 这是反"缝合怪"的真先进性 (GMD 卖点, 见 docs/model_comparison/)。

Fluspect 9 参数 (`fluspect_B_CX_PSI_PSII_combined.m:60-68`): Cab/Cca/Cdm/Cw/Cs/Cant/N/V2Z/fqe。

| # | 参数 | 含义 | COMPASS 现状 | 能否过程模拟 | 路径 |
|---|---|---|---|---|---|
| 1 | **Cab** | 叶绿素 a+b [μg/cm²] | **已做** (`nLeafToCab.m`, cab_from_luna=1) | ✅ **已模拟** | LUNA N_lc → Cab (Evans 1989 50:1) |
| 2 | Cca | 类胡萝卜素 [μg/cm²] | select_input 0.25×Cab | 🟡 半模拟 | Cab 动态后按比例跟动; 可接 LUNA |
| 3 | **Cdm** | 干物质 [g/cm²] | 固定 0.012 | ✅ **可模拟 (最干净)** | = 叶干重/面积 = C_leaf/0.466/LAI/1e4。手算 (2026-06-02) 0.0144 vs 固定 0.012, 量级对+口径清晰+严格物理。本质=1/SLA 体现, 跟 SLA0=0.028 自洽 |
| 4 | Cw | 叶含水量 [cm] | 固定 0.009 | 🟡 **暂不能干净接** | W_leaf (植物水力 v2 W_capacity=[.003;.005;.001]m) 是**经验诊断滞后量**, 无 LAI 归一化 (initWaterStorage LAI unused), 非真叶含水。手算 W_leaf×100=0.0097≈固定0.009 是**量级巧合非物理对应**, 拿它驱动 Cw = 伪过程化。**要接先把 W_leaf 升级成真叶含水** (每叶面积=鲜−干重, 按LAI缩放) |
| 5 | Cs | 棕色/衰老色素 [-] | 固定 | 🟡 可模拟 | 物候 senescence/phen → 秋季叶黄 |
| 6 | Cant | 花青素 [μg/cm²] | 固定 | 🔴 难 | 胁迫/低温诱导, 模型暂无机制, 固定合理 |
| 7 | N | 叶肉结构参数 [-] | 固定 | 🔴 难 | **叶解剖结构 (非 N 含量!易混)**, 无过程可算, 固定合理 (CLM5/SCOPE 也固定) |
| 8 | V2Z | 紫黄质循环 [-] | 固定 0 (=关) | 🟡 可模拟 | 光保护, 与 PAR/NPQ 相关 (LUNA kNPQs); 优先级低 |
| 9 | fqe | 荧光量子效率 [-] | 固定 | 🟡 部分 | 与 LUNA 电子传递 N_et 相关 |

**结论 (2026-06-02 手算夯实后修正)**:
- **Cab(已做) + Cdm(可干净接, 严格物理) = 真过程化的两个**。
- **Cw 暂不接**: W_leaf 是经验诊断滞后量 (非真叶含水, 无 LAI 归一化), 量级凑巧≈固定值不等于物理对应, 接它=伪过程化。**要接先升级 W_leaf 成真叶含水** (鲜−干重/叶面积)。这是植物水力模块的活, 不是拿现成 W_leaf 凑。
- Cca 已半动态 (0.25×Cab, Cab 动它跟动)。Cs/V2Z/fqe 半经验可接, 优先级低。
- **N(叶肉结构)/Cant 暂固定合理** — 无物理过程, 强行模拟=伪过程化 (CLM5/SCOPE 同样固定)。
- **教训**: 起初说"Cw 最有价值(W_leaf 现成)", 手算夯实发现 W_leaf 物理太弱 — 凑量级 ≠ 过程化。Cdm 才是真严格可接的。
**下一步推荐**: Cdm ← C_leaf/0.466/LAI (严格物理, 跟 SLA0=0.028 自洽); Cw 待 W_leaf 升级。
**引用**: Fluspect Vilfan et al. 2016 RSE 186:596; Cab-N Evans 1989 Oecologia 78:9 (见 nLeafToCab.m)。

###### §12.36.1 W_leaf 升级 = 真叶含水 (VOD + Cw 通吃, 2026-06-02 设计)
**动机**: VOD (compute_VOD.m:20-37) 已把 W_leaf_m/W_stem_m 当含水量 (VWC=W×1000 kg/m²地面), 但 W_leaf 是植物水力 v2 的**经验诊断滞后量** (W_capacity=[.003;.005;.001]m 拍值, 无 LAI 归一化)。输出 W_leaf~9.7e-5 m → VOD 叶含水 0.097 kg/m², 而真叶含水应 ~0.65 kg/m² → **VOD 叶含水偏低 ~7×**。Cw (fluspect) 同样卡这。**真叶含水是 VOD 的命门, 必须修。**
**真叶含水物理**: 叶含水[kg/m²地面] = (C_leaf/0.466) × MC/(1−MC), MC=叶含水率(鲜重基)。
- 每叶面积 Cw[cm] = 叶含水/LAI/10 (喂 fluspect); VOD 用每地面总量。
**含水率 MC 起步 = 0.43**: fluspect 原默认 Cw/(Cw+Cdm)=0.009/0.021=0.43 自洽反推 (有现值出处, 非拍; 用它算的稳态 Cw≈原固定0.009 不突变)。**偏低 (阔叶通常0.55-0.65), 待长武实测换** (vod_field_sampling 采样, [[project_compass_vod_field_sampling]])。
**两路并行 (用户定: 动态水库是模型亮点)**:
- **路线A (稳态真值, 先做基线)**: 叶含水 = C_leaf/0.466 × MC/(1−MC), 碳驱动稳态。快、物理清晰。
- **路线B (动态水库=卖点)**: 升级 Plant_hydraulics_v2 — W_capacity 从拍值 [.001] 改成**A 的稳态真值** (碳驱动物理容量), W_leaf 在其附近按蒸腾抽充 (动态滞后)。**B 用 A 当容量基线 → 自洽**, 且给"叶含水日变化/干旱响应" = VOD 日内信号 (SMAP/VODCA 看的就是动态含水, 别家无 prognostic 叶水库)。
- A vs B 之差 = 量化蒸腾日变化对 VOD 贡献。茎同理 (C_stem×WDR×rwc, 已是 13c VOD 思路)。
**实施时机**: 待 v22(SLA) 验完 + 森林告一段落再动手 (避免同时太多变动)。Cdm 可先接 (不依赖 MC, =C_leaf/0.466/LAI)。

###### ⚠️⚠️ §12.36.2 核 ClimaLand 源码后重大反转 (2026-06-02)
**用户质疑"W_leaf 模仿 ClimaLand 不该有问题" → 核 ClimaLand 源码, 我之前诊断错了, 真相反转**:
ClimaLand PlantHydraulics.jl 的 prognostic 水 = **ϑ_l (体积含水率) × AI × dz** (PlantHydraulics.jl:461/494 ∂ϑ_l/∂t=1/(AI·dz)·Δfa), 容量上限 = **ν × AI × dz**。
**关键: ClimaLand ν (plant_nu) = 1.44e-4 m³/m³** (toml/default_parameters.toml:92, "Porosity") — **极小!** 不是叶含水率 0.6-0.7。
→ ClimaLand 的 W **不是"叶片总含水量", 是"水势驱动的弹性储水 (capacitance)"** (量级 ~1e-7 m, 只算可交换/释放的那点水, 服务水势↔气孔反馈)。**ClimaLand 自己不拿它算 VOD/总含水** (Vegetation/ 无 VOD 模块, grep 证实)。
**所以**: COMPASS 拿 W_leaf 喂 VOD 是**概念错配** — capacitance 储水(~1e-7m) ≠ 叶总含水(~0.65kg/m², 结构+自由水60%鲜重), 差 6 量级。**不是"移植丢 LAI/ν 偏低 7 倍"** (那是我之前的错判), 是"用错了量"。用户对的是: ClimaLand W 本身没问题; 错在拿它喂 VOD。
**正解 (VOD 含水必须碳驱动, 同 13c 茎双池 M_stem=C×WDR×rwc 思路)**:
```
VOD 叶含水[kg/m²] = (C_leaf/0.466)×MC/(1−MC)  ×  [1 + capacitance相对波动]
                   └─ 总含水(大, 碳驱动, 路线A) ─┘   └ 动态波动(小, W_leaf capacitance 给) ┘
```
- 总含水量(主) 从碳算 (路线A, MC=0.43 起步)
- capacitance W_leaf 的日变化(次) 作动态调制 → 给 VOD 日内/干旱信号 (这才是"动态水库亮点"的正确用法: 不是 W 当总含水, 是 W 的波动调制碳驱动总含水)
- 茎同理 (13c 已是此思路)。fluspect Cw 同样用碳驱动总含水/LAI。
**作废**: §12.36.1 "W_capacity 改 ν×LAI×dz 当容量喂 VOD" 的方向 (基于 ν~0.65 误判); ClimaLand ν=1.44e-4 是 capacitance 不是总含水, 不能当总含水容量。

###### ✅ §12.36.3 最终设计: 叶照茎 M=C×WDR×rwc 统一 (2026-06-02 定调)
**统一公式 (叶茎一致, = 13c 茎双池思路)**:
```
M_organ[kg/m²] = (C_organ_kgm2 × WDR_organ) × rwc_organ
                 └─ 碳驱动饱和含水(静态基线,大) ─┘  └相对饱和度(动态,0~1)┘
```
**为何乘不加** (用户问): capacitance 是**相对波动**, 不是绝对补水。
- 乘: 总含水在碳定基线上下按比例波动 (干旱 rwc 0.9→0.7 = 总含水降10-30%), 物理对。
- 加: capacitance 绝对储水 ~1e-7m=1e-4kg/m² vs 碳驱动总含水 ~0.65kg/m², **加 = 被淹没** (0.65+0.0001≈0.65), 动态信号丢失, 且口径不一 (capacitance 储水 vs 总含水)。
- 茎已是乘 (M_stem=C×WDR×rwc_stem), 叶照抄。
**rwc 来源 = capacitance 相对饱和度** (动态水库的正确用法, 不浪费 ClimaLand 那套):
- `RWC = W_new / W_capacity` (储水/容量, Plant_hydraulics_v2_wrapper.m:108-109)
- **wrapper 已算 RWC_stem **和** RWC_leaf** (line 108/109)! 茎 RWC_stem 已注册(rwc_stem_v2)喂 VOD; **叶 RWC_leaf 算了但丢 diag 没接** — 这就是当前叶茎不对称 (茎 C×WDR×rwc 有动态, 叶仅 C×WDR 缺 rwc)。
**改动清单 (3 处, 叶照茎补 rwc, 不用改算法)**:
1. COMPASS.m:2249 旁 + registerPlantHydraulics.m: 注册 `rwc_leaf_v2` = v2_diag.RWC_leaf (照 rwc_stem_v2)
2. load_compass_state.m: 读 rwc_leaf_v2 (照 rwc_stem)
3. compute_VOD.m: 叶 `M_leaf = C_leaf×WDR_leaf×rwc_leaf` (现 agb_leaf_proxy=C×WDR 缺 rwc, 补上); fluspect Cw 同 (=M_leaf/LAI/单位)
**待核参数 (实施时, 有据再改)**:
- **WDR_leaf 与 MC 须一致**: 现 WDR_leaf=2.5(注"70%MC"=0.7/0.3=2.33), 但 fluspect 隐含 MC=0.43→WDR=0.43/0.57=0.75 — **打架! 同一含水率两种表达不能各设各的**, 统一到一个 MC 源 (起步 0.43, 待长武实测)。
- rwc_leaf 的 W_capacity 口径 (capacitance 容量, 与茎一致)。
**特色齐全**: 碳(C_leaf)→干生物量, capacitance(RWC_leaf=W/W_cap)→动态饱和度, 相乘→总含水→VOD+fluspect Cw。capacitance 动态水库成 rwc 物理来源 (=ClimaLand prognostic water 的正确用途), 叶茎统一, 别家无此"碳氮水→微波/光学"端到端动态链。
**实施时机**: v22(SLA) 验完 + 森林告一段落。

