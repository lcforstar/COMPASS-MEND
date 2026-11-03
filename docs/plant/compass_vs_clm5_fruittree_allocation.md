# COMPASS vs CLM5-FruitTree 分配/参数对比

**目的**: COMPASS 苹果 PFT 23 的 C/N 分配、C:N、物候、储备参数与 CLM5-FruitTree (Dombrowski/Dabbaghi et al. 2022 GMD 15:5167) 对标。
FruitTree 是**专门给落叶果树 (苹果验证)** 的 CLM5 子模型，与 COMPASS 同物种、同 GDD 物候框架，参数全部有 Zanotelli 苹果园实测出处 —— 是 COMPASS 苹果参数的金标准参照。

**源文件**: `data/references/allocation_allometry/Dabbaghi_2022_CLM5_FruitTree_GMD.pdf` (27 页, 已下载)
**调研日期**: 2026-05-29
**状态**: [更新 2026-05-29] C:N 4 项 + gs_g1 已改 (v10 验证中); slatop 待核; 动态分配走 cn_coupling §12.32.11 三套架构 (苹果走 FruitTree 物候)。

---

## 1. 核心发现 — 固定分配 vs GDD 动态分配

**COMPASS 现状**: `C_organ = NPP × f_alloc_organ`, f_alloc 固定 (leaf 0.40/root 0.25/stem 0.20/fruit 0.15), 果实仅 DOY 150-290 门控。

**FruitTree**: 分配系数随 **GDD 物候阶段动态变** (Eq 4-7):
- 叶分配 `aleaf`: 季初 fleafi=**0.85** → 指数衰减 → 末 aleaff=**0.01** (Eq 4, 指数因子 bfact=−0.5)
- 根分配 `afroot`: 初 **0.22** → 渐降至终值
- 茎分配 `alivestem` = 1 − arepr − afroot − aleaf (Eq 5, 余数; 衰减曲线 declfact=4, allconss=1.5)
- 果实分配 `arepr`: 果实期增长到 plateau (Eq 6-7)
- 采后: 全部新碳 → 储备池

**FruitTree 物候阶段** (GDD 从 budburst 起算, baset=4°C):
| 阶段 | GDD 阈值 (°days) | COMPASS Phase D 对应 |
|---|---|---|
| (1) bud break | GDDleaf=0 | budburst (Sykes chilling-GDD) |
| (2) 果实开始 | GDDfruit=400 | GDD_bloom_after_budburst=369 ✓接近 |
| (3) 果实成熟 | GDDripe=1100 | — |
| (4) 冠层成熟 | GDDlfmat=1350 | — |
| (5) 采收 | GDDmat=2880 | TSum_fruit=2208 |
| (6) 落叶起始 | crit_temp=278.15K(5°C) | crit_dayl / crit_temp |

**结论**: FruitTree 证实"固定 f_alloc 会出错" — 真实苹果叶分配季初高(0.85)随物候衰减，果实膨大期增长，不是恒定比例。COMPASS 已有 GDD 物候 (Phase D Sykes)，技术上可对接 Eq 4-7。

---

## 2. 参数逐项对比 (FruitTree 苹果实测 vs COMPASS 现状)

### 2.1 C:N 比 (⚠️ 几项差异大, FruitTree 有 Zanotelli 实测出处)

| 参数 | FruitTree 苹果 | COMPASS PFT23 | 差异 | 出处 (FruitTree) |
|---|---|---|---|---|
| leafCN | **19.7** | CN_L=25 | COMPASS 偏高 27% | Zanotelli 2010 实测 6 样本均值 |
| frootCN | **32** | CN_target_root=35 | 接近 | Zanotelli 2010 |
| **grainCN (果实)** | **139** | CN_target_fruit=30 | ⚠️ **差 4.6×** | Zanotelli 2010 实测 |
| livewdCN | **60** | CN_target_stem=100 | COMPASS 偏高 | Zanotelli 2010 |
| lflitCN (凋落物) | 60 | — | — | Zanotelli 2010 |

**最大问题: grainCN=139 vs COMPASS 30**。苹果果实高糖低氮 (C:N 高达 139), COMPASS 用 30 严重高估果实含氮 → 高估果实 N 需求/采收 N 移除。**建议核查并改向 139** (有实测出处)。

### 2.2 分配系数 (FruitTree Eq 4-7, COMPASS 无动态对应)

| 参数 | FruitTree 值 | 含义 | COMPASS |
|---|---|---|---|
| fleafi | 0.85 | 叶初始分配系数 | 无 (固定 0.40) |
| aleaff | 0.01 | 叶最终分配系数 | 无 |
| afroot (init) | 0.22 | 根分配初值 | 固定 0.25 |
| bfact | −0.5 | 叶分配指数因子 (Eq 4) | 无 |
| declfact | 4 | 茎分配衰减形状 (Eq 6) | 无 |
| allconss | 1.5 | 茎分配曲线幂 | 无 |
| fcur | 0.95 | 进显示池比例 (→5% 进储备) | cushion 等价 (✓已加 v8) |
| aleafstor | 0.3 | 采后叶储备分配 | 无 |

### 2.3 物候 GDD 阈值 — ⚠️ 不可移植! 仅交叉验证, 不可照搬

**重要区分 (2026-05-29)**: C:N 是**物种生理属性** (苹果叶 C:N≈20 跨站点稳定) → 可借 FruitTree 实测。
但 **GDD 物候是站点气候属性** (取决于当地温度积累) → **绝不能借 FruitTree**!
- COMPASS GDD 是**长武 25 年实测标定** (GDD_bloom=369 来自 2000-2024 Changwu 日 Tmax/Tmin, σ=38.9)
- FruitTree GDD 是 **Zanotelli 意大利 South Tyrol 果园** (不同站点/品种/气候)
- 两者**只能交叉验证量级合理性, 不可互相覆盖**。COMPASS 物候保持长武标定值不动。

| 参数 | FruitTree (意大利) | COMPASS (长武 25yr) | 评判 |
|---|---|---|---|
| baset (GDD 基温) | 4°C | T_base=5.59 (Wang 2009) | 都苹果园, 站点差异正常 |
| GDDfruit | 400 | GDD_bloom=369 | ✓ 量级一致, 各自站点标定, **不互改** |
| GDDmat | 2880 | TSum_fruit=2208 | 站点气候差异 (意大利 vs 长武), **COMPASS 保持长武值** |
| crit_temp | 278.15K | crit 落叶阈值 | 量级一致 |
| ndays_stor | 50 d | spring_release_days=60 | 接近, 各自站点 |
| laimx | 3 (意大利) | Wang 247obs ~5 (长武) | 长武密植, **保持长武实测** |
| flivewd | 0.15 | f_livewood=0.22 (v6 site-tune) | FruitTree 用 CLM5 BDT 默认 0.15 |

### 2.4 光合/结构 (⚠️ 部分差异)

| 参数 | FruitTree 苹果 | COMPASS | 差异 |
|---|---|---|---|
| slatop (SLA) | 0.028 m²/gC | SLA0=0.0125 | ⚠️ 差 2.2× (影响 LUNA lnc!) |
| medlynslope (g1) | 8.2 | g1_medlyn=4.86 | ⚠️ FruitTree 偏高 |
| i_vcad/s_vcad (Vcmax-N) | 5.2 / 34 | — | LUNA 相关 |
| nstem (密度) | 0.33 株/m² (3×3m) | — | 长武株距待核 |
| prune_fr (修剪比例) | 0.85 死茎 | — | COMPASS 无修剪模块 |
| ztopmx (树高) | 3.6 m | canopy_height | 核对 |
| taper (高:胸径比) | 120 | — | 异速相关 |

---

## 3. 可借鉴优先级 (状态更新 2026-05-29)

| 优先 | 项 | 改动 | 风险 | 出处 | 状态 |
|---|---|---|---|---|---|
| 🔴 高 | grainCN 30→139 | PFT 表 1 值 | 低 | Zanotelli 实测 | ✅ 已改 (v10) |
| 🔴 高 | leafCN 25→19.7 | PFT 表 1 值 | 低 (但影响 LUNA/N需求) | Zanotelli 实测 | ✅ 已改 (v10) |
| 🟡 中 | livewdCN 100→60 (CN_target_stem) | PFT 表 1 值 | 低 | Zanotelli 实测 | ✅ 已改 (v10) |
| 🟡 中 | frootCN 35→32 (CN_target_root) | PFT 表 1 值 | 低 | Zanotelli 实测 | ✅ 已改 (v10) |
| 🟡 中 | gs_g1 (medlynslope) 4.86→8.2 | config | 中 (影响 GPP/ET) | Zanotelli 实测 | ✅ 已改 (v10) |
| 🟡 中 | slatop 0.0125→0.028 核查 | 影响 LUNA lnc + LAI 标定 | 中 (动 LAI 校准 RMSE) | Zanotelli vs Wang 2009 冲突, 需查 | ⏳ 待核 (不盲改) |
| 🟢 低 | GDD 动态分配 (Eq 4-7) | 苹果走物候动态分配 | 高 | Dabbaghi 2022 | ⏳ 设计完成 (cn_coupling §12.32.11), 待实现 |
| 🟢 低 | 修剪模块 (prune_fr) | 新增 | 中 | 果园管理 | ⏳ 未做 |

**注意冲突点**:
- **slatop**: FruitTree 0.028 vs COMPASS 0.0125 (Wang 2009)。COMPASS LAI 校准 (247 obs, RMSE 0.348) 基于 0.0125, 改 slatop 会动 LUNA lnc 和 LAI。需先查 Wang 2009 原文 SLA 到底多少, 不能盲改。
- **grainCN**: 改 139 影响果实 N 池 + 采收 N 移除 (harvest_fraction), 要同步检查 N 守恒。
- **GDD 动态分配**: 物理最优但工作量大 (重构 computeCarbonAllocation + allocatePlantN), 且 FruitTree Eq 4-7 系数 (fleafi/bfact/declfact) 都是 Zanotelli 单站标定, 移植需重标。

---

## 4. 其他参考 (本文件夹)

- `CLM50_Tech_Note_Lawrence_2018.pdf` — CLM5 全 PFT 分配 (froot_leaf/stem_leaf, Ghimire 2016 源), 大叶尺度
- (待补) Ghimire 2016 JAMES — CLM5 分配原始出处 (escholarship 限流, TODO)
- (待补) Poorter 2012 New Phytol — 全球 1200 种 leaf/stem/root 分配 meta (Wiley 付费, TODO 找开放版)

**全 PFT vs 苹果专用**: COMPASS 是全 PFT 模型。FruitTree 参数只适用苹果 PFT 23; 其他 23 个 PFT 的分配仍应走 CLM5 全 PFT 方案 (Tech Note froot_leaf/stem_leaf) 或保持现有 f_alloc。苹果可单独用 FruitTree 实测值。
