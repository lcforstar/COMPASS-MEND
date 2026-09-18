# COMPASS 碳氮–植被–土壤–光合耦合图

> COMPASS v30-2026 当前实现的碳（C）、氮（N）、植被、土壤、光合作用相互影响关系。
> 节点标签后附关键变量名/公式，便于回溯源码。本图据 `src/` 实际代码绘制，非教科书版。

---

## 🔗 总耦合图 / Full coupling diagram

四层（光合 → 碳分配 → 氮分配 → 土壤氮）自上而下为主物质流；粗箭头 `==>` 标出**①叶氮闭环**，虚线 `-.->` 为驱动/需求信号。

```mermaid
flowchart TB
    accTitle: COMPASS C-N Vegetation-Soil-Photosynthesis Coupling
    accDescr: Carbon, nitrogen, vegetation and soil coupling in COMPASS v30. Leaf nitrogen drives Vcmax through LUNA and sets GPP; GPP minus respiration gives NPP which fills the NSC pool and is allocated to organs; carbon allocation sets the nitrogen demand that drives root uptake and depletes soil mineral nitrogen; root uptake refills the labile nitrogen pool and closes the loop back to leaf nitrogen; seasonal reserves buffer growth across years.

    subgraph photo ["光合作用层 · Photosynthesis (GPP)"]
        lai["LAI 叶面积指数"]
        lnc["叶氮浓度 lnc = N_leaf / LAI"]
        luna["⚙️ LUNA 日步优化 → Vcmax25 / Jmax25"]
        sfactor["水分胁迫 sfactor = K(ψ_leaf)"]
        medlyn["⚙️ Medlyn 气孔导度 gs (VPD)"]
        gpp["📊 GPP — Farquhar 协同限制"]
        lnc --> luna
        luna -->|"Vcmax ↑"| gpp
        sfactor -->|"× Vcmax/Jmax/TPU"| gpp
        medlyn -->|"gs, Ci"| gpp
        lai -->|"光截获 / 分母"| lnc
    end

    subgraph carbon ["碳分配层 · Carbon (C-dominant + 弹性C:N)"]
        npp["NPP = GPP − R_m − R_g"]
        resp["呼吸 R_m(器官N, T) + R_g(0.25)"]
        c_labile["速效碳库 NSC C_labile"]
        alloc["⚙️ 碳分配 NPP × f_alloc (GDD 物候)"]
        organs_c["器官碳库 C: 叶/活茎/死茎/细根/粗根/果"]
        c_storage["💾 储备碳 stem/root_storage (FATES taper)"]
        resp -->|"− R_a"| npp
        npp --> c_labile
        c_labile --> alloc
        alloc --> organs_c
        c_labile -->|"盈余 surplus"| c_storage
    end

    subgraph nitro ["氮分配层 · Nitrogen (demand-driven)"]
        d_n["N 需求 D_N = NPP × f_alloc / CN_target"]
        f_n["⚠️ N 限制 f_N = N_labile / D_N"]
        flex_cn["弹性 C:N → [CN_target, +offset]"]
        n_labile["速效氮库 N_labile"]
        organs_n["器官氮库 N: 叶/茎/根/果"]
        n_storage["💾 储备氮 + 回流 retransn"]
        d_n --> f_n
        n_labile --> f_n
        f_n -->|"× A_N"| organs_n
        f_n -->|"f_N < 1"| flex_cn
        organs_n --> n_storage
    end

    subgraph soil ["土壤氮循环层 · Soil N"]
        org_n["土壤有机氮 SOM-N"]
        miner["矿化 mineralization"]
        nh4["NH4 铵盐 (μg/cm³)"]
        nitri["硝化 nitrification"]
        no3["NO3 硝酸盐 (μg/cm³)"]
        denit["反硝化 → N2 / N2O 损失"]
        uptake["⚙️ 根系吸收 DAISY ADE + Nacry MM"]
        fert["📥 施肥输入 fertilization"]
        org_n --> miner --> nh4
        nh4 --> nitri --> no3
        no3 --> denit
        fert --> nh4
        fert --> no3
        nh4 -->|"NH4"| uptake
        no3 -->|"NO3"| uptake
    end

    %% ---- 跨层反馈环 / cross-layer feedback loops ----
    gpp --> npp
    organs_n -.->|"驱动 R_m"| resp
    npp -.->|"× f_alloc / CN_target"| d_n
    flex_cn -->|"actual_CN 钳制"| alloc
    organs_c -->|"② C_leaf → LAI"| lai
    uptake ==>|"① U_N 吸收"| n_labile
    n_labile ==>|"N_leaf 分配"| organs_n
    organs_n ==>|"① 叶氮闭环"| lnc
    d_n -.->|"③ 需求驱动吸收"| uptake
    n_storage -->|"④ 春季 fstor2tran"| n_labile

    classDef photo_s fill:#dbeafe,stroke:#2563eb,stroke-width:2px,color:#1e3a5f
    classDef carbon_s fill:#dcfce7,stroke:#16a34a,stroke-width:2px,color:#14532d
    classDef nitro_s fill:#ede9fe,stroke:#7c3aed,stroke-width:2px,color:#3b0764
    classDef soil_s fill:#ffedd5,stroke:#ea580c,stroke-width:2px,color:#7c2d12

    class lai,lnc,luna,sfactor,medlyn,gpp photo_s
    class npp,resp,c_labile,alloc,organs_c,c_storage carbon_s
    class d_n,f_n,flex_cn,n_labile,organs_n,n_storage nitro_s
    class org_n,miner,nh4,nitri,no3,denit,uptake,fert soil_s
```

---

## 🔄 四大反馈环 / Four feedback loops

| # | 反馈环 | 路径 | 物理含义 |
| --- | --- | --- | --- |
| **①** | 叶氮–光合闭环 | `lnc → LUNA → Vcmax → GPP → NPP → D_N → 吸收 U_N → N_labile → organs_n → lnc` | 缺氮 → Vcmax↓ → GPP↓ → 需求↓ → 吸收↓ → 叶氮持续低，系统锁定在**氮限制态**，直到矿化/回流补给[^1][^2] |
| **②** | 碳–LAI 正反馈 | `GPP → NPP → 分配 → C_leaf → LAI → 光截获 → GPP` | 叶碳长出更多 LAI，截获更多光增 GPP；受 `peaklai`（LAI≥laimx）和休眠物候门控封顶[^2] |
| **③** | 弹性 C:N 缓冲 | `f_N < 1 → actual_CN↑ → 器官碳被钳制 → 剩余留 C_labile → 储备碳` | 氮缺时器官 C:N 上浮到 `[target, target+offset]`，多余碳不强塞器官而入储备，C 不浪费[^3] |
| **④** | 储备越冬多年信号 | `秋季落叶回流 retransn → N_storage → 春季 fstor2tran → N_labile → 叶萌发` | 上一年氮状况经储备库跨年传递，决定次年春季萌发强度（苹果大小年）[^4] |

> **水分旁路**：`sfactor = K(ψ_leaf)` 来自植物水力，乘在 Vcmax/Jmax/TPU 上（不乘 Rd），干旱直接压低 GPP→C_labile，与①②环叠加。

---

## 📋 各节点对应源码 / Source map

| 层 | 关键文件 | 关键函数 / 变量 |
| --- | --- | --- |
| 光合 | `src/+plant/biochemical.m`；`src/+plant/+integrated/+LUNA/computeLuna.m` | Farquhar+Medlyn 迭代；`lnc`→`FNCa`→`Vcmax25/Jmax25`；`sfactor` |
| 碳分配 | `src/+plant/+PlantCarbon/computeCarbonAllocation.m`；`computeNPP.m` | `NPP=GPP−R_m−R_g`(Sitch 2003)；`C_labile`；`clampPoolCN`；FATES `refillTo` |
| 氮分配 | `computeOrganDemand.m`；`allocatePlantN.m`；`runPlantNDemand.m` | `D_N=NPP·f_alloc/CN_target`；`f_N_growth`；`actual_CN` |
| 土壤氮 | `src/+soil/+nitrogen/{+organic/mineralization, +chemistry/nitrification, +root_uptake/nitrogenUptake}.m` | `Smin`、`Snit`、`SinkS`(DAISY+Nacry)；`solveNitrogen.m` 调度 |
| 主循环 | `src/COMPASS.m` | `runPlantNDemand → applyRootUptake → runPlantCNStep → solveNitrogen` |

---

[^1]: Sitch, S. et al. (2003). LPJ 动态全球植被模型——NPP 与氮需求/限制公式 (Eq. 12, 24, 41). *Global Change Biology* 9(2): 161–185. https://doi.org/10.1046/j.1365-2486.2003.00569.x
[^2]: Ali, A.A. et al. (2016). LUNA 模型——叶氮投资优化 Vcmax/Jmax (CLM5 LunaMod). *Geoscientific Model Development* 9: 587–606. https://doi.org/10.5194/gmd-9-587-2016
[^3]: Medlyn, B.E. et al. (2011). 最优气孔导度模型 (g1/√VPD). *Global Change Biology* 17(6): 2134–2144. https://doi.org/10.1111/j.1365-2486.2010.02375.x
[^4]: FATES / Massoud, E.C. et al. (2019). 非结构性碳储备目标驱动平滑充注 (smooth taper). https://doi.org/10.5194/bg-16-3499-2019
