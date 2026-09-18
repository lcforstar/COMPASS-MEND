# COMPASS vs CLM5 / LPJ-GUESS / Noah-MP / ClimaLand — 全模块先进性对比

> **用途**: GMD 文章用，论证 COMPASS 的先进性，反"缝合怪"。
> **生成方式**: 2026-06-01 由 6 个并行 agent 读 5 个模型本地源码 + 1 个 DAISY 旁证 整理；关键论断已亲验。
> **下次怎么用**: 先读本文档拿结论 + 源码出处坐标，只在要"精细化某一处"时再去读对方那几行源码，不必重读整仓库。
>
> **核实状态标注**:
> - ✅ = 已读真源码逐行核实出处（文件:行号确认）
> - ⚠️ = agent 报告，本人未逐行复核（写论文前需再查那几行）
> - ❌ = agent 报告有误，已纠正（见纠错栏）

---

## 0. 对比对象 + 本地源码路径（下次直接按此坐标查）

| 模型 | 语言 | 本地源码根 | 备注 |
|---|---|---|---|
| **COMPASS** | MATLAB | `f:/COMPASSV1.0/src/` | 本模型 |
| **CLM5** (CTSM 5.0.37) | Fortran | `D:/OneDrive/CTSM-release-clm5.0.37/src/` | 全功能对照标杆 |
| **LPJ-GUESS** v4.1.1 | C++ | `D:/OneDrive/LPJ-GUESS-v4.1.1/guess_4.1/` | DGVM, need-based N |
| **Noah-MP** | Fortran | `D:/OneDrive/NoahMP_Dual-main` | 业务化 LSM, 无 C/N 循环 |
| **ClimaLand.jl** | Julia | `D:/OneDrive/应用/Model/ClimaLand/src/` | 现代 Julia LSM |
| **DAISY** 7.1.0 | C++ | `D:/OneDrive/DAISY/daisy-7.1.0/src/` | 仅 C/E N 化学旁证 |
| LPJmL | — | 无本地源码 | 仅论文 von Bloh 2018 GMD 11:2789 |

**13 模块**: A土壤水 / B土壤热+冻融 / C土壤N / D根吸水RWU / E根吸N / F光合 / G气孔+植物水力 / H自养呼吸 / I C/N池+分配 / J茎周转+边材 / K N需求+储备 / L物候 / M微波VOD。

---

## A. 土壤水 (Richards)

| 维度 | COMPASS | CLM5 | LPJ-GUESS | Noah-MP | ClimaLand |
|---|---|---|---|---|---|
| 方程形式 | 混合 θ-h ⚠️ | θ 形式(混合可选但禁用) | 两层 bucket(非Richards) | bucket+简化层 | 纯 h 形式 Richards |
| 时间积分 | Picard 外 + 自适应内 | 半隐式线性化单步 | 日 bucket | 隐式 bucket | IMEX(ARS111) Newton |
| 线性解 | Thomas(零主元守卫) | Thomas(TridiagonalMod) | 层间直接重分配 | Thomas(SMFLX) | Julia `\`/GMRes |
| 对数空间/自适应步 | **是**(头跳>200cm对数阻尼ω=0.5) ⚠️ | 无显式对数 | 无 | 无 | 用户指定dt |
| 边界条件 | 顶通量+底地下水静水 | bc_flux顶/bc_aquifer底 | 隐式顶入渗零通量底 | 顶入渗 | FluxBC/MoistureStateBC |
| **多相(水汽气)** | **三相显式(RHOV/D_V/Khh_vapor)** ⚠️ | 单相液+冰 | 单相液 | 单相 | 单相液+冰 |

**源码坐标**:
- COMPASS: `src/+soil/+water/calculateMatricCoefficients.m`(多相系数,~L19-122 ⚠️未逐行验); `solveTridiagonalMatrixEquations.m:23-68`(Thomas+对数阻尼); `solveSoilMoistureBalance.m:38-88`(Picard)
- CLM5: `SoilWaterMovementMod.F90:968`(moisture_form), `:449-466`(线性化+Thomas), `:50-55`(bc); `TridiagonalMod.F90`
- LPJ-GUESS: `modules/soilwater.cpp:160-400`(两层 Gerten)
- Noah-MP: `module_sf_*` SMFLX(行号 agent 未取到 ⚠️)
- ClimaLand: `docs/tutorials/standalone/Soil/richards_equation.jl:13-26`

**COMPASS 优势(真)**: 三相(水-汽-气)Richards + 显式汽相迁移，5 个里唯一(STEMMUS 遗产)；水热同步隐式 vs CLM5 顺序。

---

## B. 土壤热 + 冻融

| 维度 | COMPASS | CLM5 | LPJ-GUESS | Noah-MP | ClimaLand |
|---|---|---|---|---|---|
| 热方程耦合 | 与水**同步**隐式(交叉导数CTh/KhT) ⚠️ | Crank-Nicolson,水后能量**顺序** | 无能量PDE(仅诊断) | bucket | 隐式耦合(IMEX) |
| 冻融潜热 | 表观热容(DTheta_LLT)+Clausius-Clapeyron ⚠️ | 表观热容(PhaseChange_beta) | 无 | bucket冰隐式 | 表观热容(推测) |
| 过冷水 | Clausius-Clapeyron(lFreezing,hh_frez) ⚠️ | C-C 经基质势 h(T) | 无 | 无 | VG+基质势(推测) |
| 导热参数 | **de Vries** 孔隙网络 ⚠️ | Johansen(Farouki 1981) | 无 | 经验TDFCND | 待查 |

**源码坐标**:
- COMPASS: `src/+soil/+thermal/solveThermalEquations.m:1-135`; `Latent.m`(L_f表); `+conductivity/`(de Vries,目录 ⚠️未列具体文件)
- CLM5: `SoilTemperatureMod.F90:126-141`(C-N), `:429-476`(PhaseChange), `:305-311`(SoilThermProp Johansen)
- ClimaLand: 热属性 agent 未读到 ⚠️

**COMPASS 优势(真)**: 水-热同步隐式(vs CLM5 顺序单向)；de Vries 机理导热(vs Johansen 半经验)；冻融过冷水 + 多相。

---

## C. 土壤 N 循环

| 维度 | COMPASS | CLM5 | LPJ-GUESS | Noah-MP | ClimaLand | DAISY |
|---|---|---|---|---|---|---|
| N 池 | NH4/NO3(含**冰相**)+有机分层+尿素+颗粒肥 | NH4/NO3+分解级联 | 仅有机C,无独立N池 | **无** ✅ | **无** ✅ | NH4/NO3(无冰相) |
| 硝化 | MM+温湿+pH ⚠️ | 一阶 k×[NH4]+温湿 | 无 | 无 | 无 | MM+温湿 |
| 反硝化 | MM+缺氧+T ⚠️ | f(resp)×f(NO3)+缺氧 | 无 | 无 | 无 | 两阶段+温湿 |
| **N 随水迁移** | **对流-弥散(Theta_Solute驱动)** ⚠️ | 层内扩散无对流 | 无 | 无 | 无 | 仅单根ADE |
| **冻融 N 分配** | **有(液-冰相K_segr=0.3)** ⚠️ | 无 | 无 | 无 | 无 | 无 |

**源码坐标**:
- COMPASS: `src/+soil/+nitrogen/solveNitrogen.m`(总调度,冻融分配~L136-178,ADE~L273); `+chemistry/nitrification.m`,`denitrification.m`,`freezeThawPartitioning.m`; `+organic/mineralization.m`; `+transport/solveTransportStep.m`
- CLM5: `SoilBiogeochemNitrifDenitrifMod.F90:35-46`(硝化), `:132-177`(反硝化); `SoilBiogeochemNitrogenStateType.F90:29-40`(N池)
- LPJ-GUESS: `modules/somdynam.cpp:48-66`(有机C,无独立N)
- ClimaLand: `src/standalone/Soil/Biogeochemistry/Biogeochemistry.jl`(仅CO2,无N) ✅
- DAISY: `src/daisy/chemicals/nitrification_soil.C:56-80`, `reaction_denit.C:69-135`

**COMPASS 优势(真)**: N 对流-弥散直接挂 STEMMUS 水迁移 + 冻融液冰相 N 分配，5 个里唯一。✅ 已验 ClimaLand/Noah 确无土壤 N。

---

## D. 根系吸水 RWU

| 维度 | COMPASS | CLM5 | LPJ-GUESS | Noah-MP | ClimaLand |
|---|---|---|---|---|---|
| 范式 | Sperry 土-根-茎-叶串联 ⚠️ | Feddes α+HydStress(经验) | 供需经验 | Feddes BTRAN | **Sperry**(诊断型) ✅有PlantHydraulics |
| 根导度拆分 | 可选root_krmax(katul/clm5/radial) | k_soil_root整合 | 无显式 | 无(仅BTRAN标量) | 土Ks缩放 |
| 重力项 | 显式逐节点(soil_gravity_MPa) | 显式(grav2) | 无 | 隐含 | 显式q=-K(Δψ/dz+1) |
| 气孔反馈 | **双向**(psir⇄psil bisection) | 单向down-regulate | 供需单向 | 单向BTRAN | 双向(Sperry) |

**源码坐标**:
- COMPASS: `src/+plant/+root/assembleRWU.m:41-44`(重力); `+plant/+PlantHydraulic/Plantrootconductivity.m:46-92`(根导度); `calSPAC.m`(SPAC双向)
- CLM5: `SoilWaterPlantSinkMod.F90:281`(k_soil_root), `:310`(grav2), `:242-339`
- LPJ-GUESS: `modules/canexch.cpp:1895-2022`(aet_water_stress)
- ClimaLand: `src/standalone/Vegetation/PlantHydraulics.jl:199-622` ✅ 确认存在

**COMPASS 优势(部分真)**: Sperry **非独有**(ClimaLand 也有)；独特在 Sperry + prognostic W_stem 储水(见 G/M)。

---

## E. 根系吸 N

| 维度 | COMPASS | CLM5 | LPJ-GUESS | Noah-MP | ClimaLand | DAISY |
|---|---|---|---|---|---|---|
| 被动(随蒸腾) | 有(q_w×C_bulk) ⚠️ | 有(Npassive∝rwu) | 无 | 无 | 无 | 有 |
| 主动(MM) | **Nacry膜双组分(Vmax/Km分NH4/NO3+Cmin)** ⚠️ | Nactive/Nnonmyc MM | 无 | 无 | 无 | 单根膜 |
| 土-根耦合 | **DAISY圆柱ADE+Nacry膜迭代求C_root** ⚠️ | 简单供需,无土壤物理 | 无 | 无 | 无 | DAISY ADE-膜 |

**源码坐标**:
- COMPASS: `src/+soil/+nitrogen/+root_uptake/nitrogenUptake.m`(Nacry MM ~L55-76,121-126); `solveRootSurface.m:63-120`(bisection求C_root)
- CLM5: `NutrientCompetitionCLM45defaultMod.F90:247-250`(Npassive/Nactive/Nnonmyc)

**COMPASS 优势(真)**: DAISY 圆柱 ADE + Nacry 膜 MM 迭代求根表面浓度，与 STEMMUS 水热完全耦合(CLM5 是分离的被动+主动双路径,无 C_root 反馈)。

---

## F. 光合

| 维度 | COMPASS | CLM5 | LPJ-GUESS | Noah-MP | ClimaLand |
|---|---|---|---|---|---|
| 模型 | FvCB-Farquhar(CLM5 port) | FvCB-Farquhar | Haxeltine&Prentice协调 | 经验(源码未取⚠️) | FvCB-Farquhar |
| Vcmax | **LUNA日优化(5池N分配)** ⚠️ | PFT常数 OR LUNA(可选) | 静态 N→Vmax线性 | 无 | 参数/输入驱动 |
| sun/shade | big-leaf+LUNA每层_z | sun/shade显式(nlevcan) | big-leaf(pop)/分层(cohort) | 无 | 单层 |
| 气孔 | Medlyn 2011+光抑制 | Ball-Berry OR Medlyn | 协调假设内嵌 | 经验 | Medlyn 2011 |
| 温度响应 | Kattge&Knorr 2007 peaked | Kattge&Knorr 2007 | Q10查表 | 无 | Arrhenius |

**源码坐标**:
- COMPASS: `src/+plant/photosynthesisLuna.m:1-3`(CLM5 LunaMod port); `computeLuna.m:1-55`; `nitrogenAllocation.m:5-150`(Ali 2016); `biochemical.m:583-595`(Medlyn)
- CLM5: `PhotosynthesisMod.F90:1488-1517`(LUNA path), `:1404-1405`(legacy Kattge), `:1532-1548`(温度)
- LPJ-GUESS: `modules/canexch.cpp:507-541`(Haxeltine N→Vmax), `:50-65`(Q10查表)
- ClimaLand: `photosynthesis_farquhar.jl:1-145`

**COMPASS 优势(对齐先进)**: LUNA 动态 Vcmax **CLM5 也有(可选)**，属"对齐 CLM5 最先进项"，非独有。LPJ/Noah/ClimaLand 无 LUNA。

---

## H. 自养呼吸

| 维度 | COMPASS | CLM5 | LPJ-GUESS | Noah-MP | ClimaLand |
|---|---|---|---|---|---|
| 维持R_m | Lloyd&Taylor Arrhenius(Sitch Eq23) ✅ | **Q10=1.5**(Ryan 1991) ✅ | Sitch Lloyd&Taylor | 无 | JULES |
| 生长R_g | Sitch 0.25 ✅ | **无R_g(NPP=GPP-R_m)** | Sitch 0.25 | 无 | JULES Rel参数 |
| **茎呼吸用哪池** | **仅livestem N(死木不呼吸)** ✅ | **仅livestemn** ✅ `livestem_mr=livestemn*br*tc` | sapwood(无活/死区分) | 无 | ηsl倒推无分离 |
| 温度形式 | exp{308.56×[1/56.02-1/(T+46.02)]} ✅ | Q10^((t-20)/10) ✅ | 同Sitch | 无 | JULES隐含 |

**源码坐标**:
- COMPASS: `src/+plant/+PlantNitrogen/computeNPP.m:92-93`(g(T)), `:97`(R_m_root), `:99-103`(R_m_stem用N_livestem,标[13c]), `:112-113`(R_g=0.25)
- CLM5: `CNMRespMod.F90:265` ✅ `livestem_mr(p)=livestemn(p)*br*tc`(死木不在式中); `:208-242`(br=2.525e-6,Q10=1.5); CLM5 **无 R_g**
- LPJ-GUESS: `growth.cpp`(Sitch r_m10 PFT表)

**COMPASS 优势(综合)**: Lloyd&Taylor 比 CLM5 Q10=1.5 物理；有 R_g(CLM5 无)；边材只活木呼吸 ✅ 对标 CLM5:265 准确。

---

## I. C/N 池 + 分配 ⭐

| 维度 | COMPASS | CLM5 | LPJ-GUESS | Noah-MP |
|---|---|---|---|---|
| C 池 | 7池:labile/leaf/root/**livestem/deadstem**/fruit/storage | leafc/rootc/livestemc/deadstemc/cpool | cmass_leaf/root/**sap/heart**/nstore | 仅总drymass,无器官 |
| 茎分化 | **双池 livestem(CN50)/deadstem(CN500)** | livestemc/deadstemc(CN见J) | sapwood/heartwood(CN~50/~2000) | 无 |
| C 分配 | **C主导** A_C=NPP×f_alloc;新茎×flivewd→live | C主导 | C主导(allocation.cpp) | 无 |
| **C主导 vs N主导** | **C主导**(NPP先定,N=C/CN跟随) | C主导 | C主导 | 无 |
| C:N 弹性 | **FlexibleCN**[target,target+offset=10],N不足时主动clamp砍C | **FlexibleCN**(nscalar,被动降ndemand) | **刚性C:N** | 无 |

**源码坐标**:
- COMPASS: `src/+plant/+PlantCarbon/computeCarbonAllocation.m:80-128`(cn_mode), `:154-155`(flivewd拆), `:180-194`(clampPoolCN); `+PlantNitrogen/allocatePlantN.m:98-126`(actual_CN弹性+双池N拆)
- CLM5: `NutrientCompetitionFlexibleCNMod.F90:1658-1670`(nscalar+leafcn clamp)
- LPJ-GUESS: `growth.cpp`(刚性cton_*_avr)

**COMPASS 定位(综合)**: FlexibleCN 学 CLM5 思想(实现更激进:直接clamp砍C vs CLM5被动降需求);C主导学 CLM5/LPJ。

---

## J. 茎周转 + 边材 ⭐ (13c 核心)

| 维度 | COMPASS | CLM5 | LPJ-GUESS |
|---|---|---|---|
| 活/死木分化 | livestem(CN50,呼吸)/deadstem(CN500,不呼吸) | livestemc/deadstemc | sapwood/heartwood(CN50/2000) |
| **活→死周转** | **lwtop=0.7/yr** ✅对标 | **lwtop_ann=0.7/yr** ✅ | turnover_sap物种0.05-0.4/yr |
| 转换公式 | dC=livestem×lwtop×dt; dN_to_dead=dN×(CN_lw/CN_dw) | `ctovr=livestemc*lwtop`;`livestemn_to_deadstemn=ctovr/deadwdcn` ✅ | cmass_sap→heart |
| 死木N回收 | CN_lw/CN_dw比例留死木,余retransl→storage | N按deadwdcn留,余retransl ✅ | 转换前完全retransl |
| 活木呼吸 | livestem(死木不呼吸) | livestemc(CNGRespMod:64) | sapwood贡献R_m |

**源码坐标**:
- COMPASS: `getPFTParameterSet.m`(CN_lw=50/CN_dw=500/flivewd/lwtop=0.7, PFT23~L2143-2152); `turnoverAndHarvest.m:165-180`(N活→死,retransl); `turnoverCarbonMirror.m:117-127`(C镜像)
- CLM5: `CNPhenologyMod.F90:401`(lwtop/31536000), `:2754`(ctovr=livestemc*lwtop) ✅, `:2757`(livestemn_to_deadstemn=ctovr/**deadwdcn**) ✅; `CNMRespMod.F90:265`(活木呼吸) ✅
- LPJ-GUESS: `growth.cpp:219-278`(turnover_sap)

**COMPASS 定位(综合)**: 结构最像 LPJ sapwood/heartwood,速率用 CLM5 lwtop=0.7 ✅;CN_dw=500 对标 CLM5 deadwdcn(✅见纠错)。

---

## K. N 需求 + 储备 ⭐

| 维度 | COMPASS | CLM5 | LPJ-GUESS |
|---|---|---|---|
| N需求公式 | 生长D_N + (N_labile_target-N_labile)/dt 储备缺口 | plant_ndemand_opt(0-3) | need-based: leaf+root+sap+store ndemand |
| **储备目标** | **max(sap_term, regrow_term)** | 年均retransn(无sap term) | **max(cmass_sap×fnstorage/cton_leaf, retransn)** ✅照搬源 |
| sap容量项 | C_livestem×fnstorage/CN_L (fnstorage=0.3) | **无** | cmass_sap×fnstorage/cton_leaf_bg(fnstorage=0.3) |
| 越冬储备 | fcur=0.5;(1-fcur)→storage;春萌release | fstor2tran=0.5 | scale_n_storage隐含 |

**源码坐标**:
- COMPASS: `runPlantNDemand.m:130-144`(N_labile_target=max(sap,regrow),fnstorage=0.3); `turnoverAndHarvest.m:182-189`(春萌release)
- CLM5: `NutrientCompetitionFlexibleCNMod.F90:1694-1765`(plant_ndemand opt), `CNPhenologyMod.F90:926`(fstor2tran)
- LPJ-GUESS: `growth.cpp:1491` ✅ `max_n_storage=max(cmass_sap*fnstorage/cton_leaf_bg, retransn_nextyear)`; `data/ins/europe.ins:288`(fnstorage 0.3 tree)

**COMPASS 定位(综合,sap term 照搬LPJ)**: sap 容量项逐字照搬 LPJ growth.cpp:1491(✅已亲验源),只改变量名;CLM5 缺此项。

---

## L. 物候

| 维度 | COMPASS | CLM5 | LPJ-GUESS | Noah-MP | ClimaLand |
|---|---|---|---|---|---|
| 叶展GDD | GDD5+CLM5日长门控 ⚠️ | GDD5(ws_flag) | GDD5全年累积 | 无 | 无(处方MODIS LAI) |
| 落叶 | 日长offset(ws_flag) | 日长<crit_dayl | raingreen wscal | 无 | 无 |
| 冷胁迫 | Sykes 1996 gdd0 | 无 | Sykes 1996 | 无 | 无 |
| APHEN上限 | 210天 | 无 | APHEN_MAX=210 | 无 | 无 |
| **独立果实物候** | **DS_fruit(0→1)** ⚠️ | 无 | 无 | 无 | 无 |

**源码坐标**:
- COMPASS: `src/+plant/+integrated/+phenology/seasonDecidPhenology.m`(GDD门控~L189-196,offset~L230-238,Sykes~L208-209,APHEN~L93/139,**果实DS_fruit~L299-338**)
- CLM5: `CNPhenologyMod.F90:669-971`(SeasonDecid,onset_gdd:905-908,offset:957)
- LPJ-GUESS: `modules/growth.cpp:80-178`(leaf_phenology,Sykes:82-99,APHEN:50/152)

**COMPASS 优势(真差异化)**: 独立果实物候 DS_fruit 是纯自然植被模型(CLM5/LPJ/Noah/ClimaLand)都没有的;且融合 Sykes+CLM5+LPJ 三套。

---

## M. 微波 VOD 前向 ⭐⭐ (头号卖点)

| 维度 | COMPASS | CLM5 | LPJ-GUESS | Noah-MP | ClimaLand |
|---|---|---|---|---|---|
| **内置微波/VOD** | **完整L波段前向** | **无(源码零)** ✅ | **无(源码零,仅参考PDF)** ✅ | **无(源码零)** ✅ | 实验骨架(_vod_smoke,无实现) |
| VOD=β·M·G | 是(VOD_organ=b·VWC) ⚠️ | — | — | — | 实验 |
| M_stem驱动 | prognostic W_stem / rwc_stem / empirical 三选 ⚠️ | — | — | — | — |
| **live/dead加权** | **WDR_live=1.5/WDR_dead=0.4** ⚠️ | — | — | — | — |
| ψ_leaf调b_leaf | Momen 2017 g(ψ) ⚠️ | — | — | — | — |
| 碳氮水闭环 | **C池→干生物量,rwc→相位,ψ_leaf→散射** ⚠️ | 无 | 无 | 无 | 无 |
| 物理b(CMEM) | Kirdyashev/Wigneron/Momen多源 ⚠️ | — | — | — | — |

**源码坐标**:
- COMPASS: `src/+Radiation/+microwave/compute_VOD.m`(VOD~L53-70,M_stem~L146-185,**双池加权~L118-126**,Momen ψ~L194-196,freeze禁用~L203-207,Kirdyashev~L209-258,H/V偏振~L276-290); `load_compass_state.m:78-152`(读C池/双池)
- CLM5/LPJ/Noah: **源码扫 microwave/vod/tau_omega/brightness 全空** ✅(2026-06-01 已 grep 验证)
- ClimaLand: 仅 `_vod_smoke_noplot.jl` 实验,无源码

**COMPASS 优势(真原创/唯一)**: ✅✅ 已亲验——5 个模型里**唯一**有完整"陆面C/N/水状态→L波段VOD"端到端前向算子。CLM5/LPJ-GUESS/Noah-MP 源码零微波;ClimaLand 仅实验骨架。这是 GMD 头号卖点。
**诚实**: 微波物理本身借了 Kirdyashev(CMEM)/Wigneron(L-MEB)/Momen 多源理论,但"陆面prognostic状态→VOD"的端到端集成 + 双池加权 + rwc驱动相位是原创组合。

---

## ★ agent 报告纠错栏 (本人亲验后修正)

| # | agent 原报告 | 纠正 | 依据 |
|---|---|---|---|
| 1 | (Agent5) "CLM5 无明确 CN_dw 参数,可能固定/隐含" | **错**。CLM5 **有 deadwdcn(死木C:N)**;茎周转 N 用它 | ✅ `CNPhenologyMod.F90:2757` `livestemn_to_deadstemn=ctovr/deadwdcn` |
| 2 | (Agent3/DG) ClimaLand "有Sperry"但又"诊断型"自相矛盾 | ClimaLand **有 PlantHydraulics.jl**(确存在);prognostic 程度待细查 | ✅ `src/standalone/Vegetation/PlantHydraulics.jl` 存在;⚠️ 细节未读 |
| 3 | (多 agent) 暗示 COMPASS C/N 多处"先进/独特" | **诚实定位**: I/J/K 是"综合 CLM5 FlexibleCN + LPJ need-based",**非原创算法**;真原创/唯一在 **M 微波** | ✅ K 的 sap_term 逐字照搬 LPJ growth.cpp:1491 |
| 4 | (Agent4/FH) Noah-MP 光合"经验协调" | Noah-MP 光合源码 agent **未取到**,标 ❓待补 | ⚠️ 未验 |

---

## ★ GMD 卖点排序 (反"缝合怪"的论证骨架)

1. **头炮(真唯一/原创组合)**: 陆面 C/N/水 → 微波 VOD 端到端(M) + prognostic W_stem 储水链(D/G)。**别家源码零微波**,别人复制不了的整合价值。
2. **硬底子(真优势)**: 三相土壤水 + 水-N 对流耦合 + 冻融 N 分配(A/B/C/E)。STEMMUS_N 深度,5 个里唯一。
3. **果树适配(真差异化)**: FruitTree 果实物候 DS_fruit(L)。自然植被模型都没有。
4. **C/N 机制(诚实定位"综合最优")**: CLM5 FlexibleCN × LPJ need-based 的首次耦合(I/J/K)。零件借的,组合是新的——别硬说原创。

**一句话框架**: COMPASS 不靠单模块比别人新,而靠"多相土壤物理 + 完整水-碳-氮耦合 + 微波可观测量"这条端到端链是 5 个模型里**唯一打通的**。缝合怪是拼了不通;COMPASS 是拼了能直接产出别人产不出的可观测量。

---

## ★ 写论文前待补验清单 (⚠️项)

- [ ] A/B 多相系数 `calculateMatricCoefficients.m` 逐行 + de Vries `+conductivity/` 具体文件
- [ ] C/E COMPASS 硝化/反硝化/Nacry 具体行号(agent 给的是文件名,部分行号粗)
- [ ] ClimaLand PlantHydraulics 是否真 prognostic ψ(vs 诊断)
- [ ] Noah-MP 光合/呼吸源码(agent 未取到行号)
- [ ] LPJmL 只能引论文 von Bloh 2018(无本地源码),对比表里 LPJmL 暂缺
- [ ] M 微波 compute_VOD 各物理项(Momen/Kirdyashev/Wigneron)对标原文献编号核对
