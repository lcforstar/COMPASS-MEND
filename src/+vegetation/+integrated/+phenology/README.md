# +vegetation/+integrated/+phenology

## 物理问题

为 COMPASS 提供统一的多 PFT 物候触发框架，覆盖：
- **常绿** (evergreen): 持续 leaf turnover, 无 budburst/senescence
- **季节性落叶** (season_decid): 冷激+GDD 触发萌发, 日长触发落叶 (e.g., Apple, Que_rob)
- **胁迫落叶** (stress_decid): 温度/水分胁迫驱动 (e.g., C3/C4 grass, 热带 raingreen)
- **作物** (crop_phenology): 一年生 sowing→harvest 单周期

## 架构来源 (源码验证)

学 **CLM5 `CNPhenologyMod.F90`** (`D:/OneDrive/CTSM-release-clm5.0.37/src/biogeochem/`)：
- **filter + flag gate** 模式 (CLM5 lines 295-307): 主驱动顺序调用 4 个独立子程序, 每个子程序内部通过 PFT flag 决定是否处理该 patch
- **每个 PFT 独立 state** (`_patch` 后缀, CLM5 CNVegStateType.F90:79-88): onset_flag, onset_counter, onset_gdd, offset_flag 等

学 **LPJ-GUESS v4.1.1 Sykes 1996 公式** (`D:/OneDrive/LPJ-GUESS-v4.1.1/guess_4.1/modules/growth.cpp:91-99`)：
- **季节性落叶萌发公式**: `gdd0(c) = k_chilla + k_chillb·exp(-k_chillk·c)`, `phen = min(1, (gdd5-gdd0)/phengdd5ramp)`
- 引用源: Sykes MT, Prentice IC, Cramer W (1996) J Biogeogr 23:203-233 (通过 LPJ-GUESS 源码二次引用, 原始 PDF 未下载)

学 **CLM5 senescence**: 显式日长触发 `dayl < crit_dayl=10.9h` (CNPhenologyMod.F90:957)

## 输入

- **parameters**: 每 PFT 4 个 binary flag (evergreen/season_decid/stress_decid/crop_phenology), 季节性落叶 PFT 4 个 Sykes 参数 (k_chilla/k_chillb/k_chillk/phengdd5ramp), 全 PFT crit_dayl_h + nrelocfrac
- **state** (per-PFT): chilldays, gdd5, aphen, onset_flag, offset_flag, dormant_flag, phen
- **forcing**: T_air daily, dayl (day length, hours)

## 输出

- **phen** [0, 1] : 当前叶量比例 (per-PFT). 物理代码读取:
  - LAI 模块: LAI = LAImax × phen × dHUF(shape)
  - turnover 模块: 触发春萌 N release (phen 跨过 0.05), 叶 senescence (aphen 接近 APHEN_MAX), 苹果采摘 (DS_fruit > 1.0)
- **diagnostics**: budburst_doy, senescence_doy, chilling_completed_doy 每年

## 耦合点

- 调用入口: `runPlantCNStep.m` 主循环每日开头
- LAI 写入: `+canopy/+phenology/` 或 `src/+vegetation/phenology.m` (现有 CalculateLAI)
- N 触发: `+turnover/appleFruitHarvest.m`, `leafSenescenceRetransl.m`, `springStorageRelease.m`

## 文件清单

```
phenologyDispatch.m       主入口 (filter + flag gate 模式)
seasonDecidPhenology.m    Sykes 1996 + CLM5 daylength (Apple etc.)
stressDecidPhenology.m    温度/水分胁迫 (C3/C4 grass)
evergreenPhenology.m      持续 turnover (常绿)
cropPhenology.m           一年生作物 (一期占位, 复用现有 sowing→DS)
PhenologyState.m          per-PFT state container (CLM5 _patch 风格)
```

## 状态机 (季节性落叶, COMPASS 苹果走这条)

```
冬季 dormant:
  if T_day < T_chill_base (5°C):
    chilldays++
  gdd5 += max(0, T - T_gdd_base)

春季 (Sykes 1996 检验):
  gdd0_threshold = k_chilla + k_chillb × exp(-k_chillk × chilldays)
  if gdd5 > gdd0_threshold and aphen < APHEN_MAX:
    phen = min(1, (gdd5 - gdd0_threshold) / phengdd5ramp)
    if onset_flag == false: onset_flag = true, 触发春萌 N release
  
生长期 (phen 推进):
  aphen += phen   (累计 leaf-on days)
  LAI(t) = LAImax × phen × dHUF(shape)
  
秋季 (CLM5 daylength 触发):
  if dayl_hours < crit_dayl_h (10.9h):
    offset_flag = true
    phen 渐降 → 0 (in ndays_off=15 days)
    触发 leaf senescence (N retransloc)

冷月 reset (LPJ-GUESS driver.cpp:583-586):
  if mtemp_yesterday >= 5 and mtemp_today < 5 and ifsensechill:
    gdd5 = 0; chilldays = 0
```

## 守恒考虑

物候模块本身**不创造或销毁 N/C 质量**，只决定何时触发：
- 春萌 release: N_storage → N_labile (`springStorageRelease.m`, 既有质量守恒)
- 叶 senescence: N_leaf → nrelocfrac × N_storage + (1-nrelocfrac) × N_litter
- 苹果采摘: N_fruit → harvest_export

phen [0, 1] 是 auxiliary (不影响 N 守恒)。
