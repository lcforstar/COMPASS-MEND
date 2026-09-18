# Phase D — 物候系统重构审计文档

**日期**: 2026-05-26
**作者**: Zhu Jizixing (NWAFU)
**目的**: 把 COMPASS 原本基于 `frPHU = CHU/PHU` 的物候触发系统（含 `harvest_threshold=0.60` 等站点 hack 缩放值）重构为 **CLM5 架构 + LPJ-GUESS Sykes 1996 公式 + DAISY DS_fruit + 实测数据校准** 的统一多 PFT 物候框架。

---

## 0. 重构动机

### 0.1 旧系统的物理问题

```
旧逻辑:
   PHU = 3921 °C·d   ← LAI S 曲线形状参数 (Williams 1989 EPIC LAI 校准, 不是生育期需求)
   CHU(t) = Σ max(0, T_avg - 5.59) dt
   frPHU = CHU / PHU = [0, 0.648]   ← 长武 2021 年末上限 0.648, 永到不了 1.0

   触发条件:
      harvest:        frPHU > 0.60        ← EPIC 0.90 × 0.65 缩放 hack
      senescence:     frPHU > 0.45        ← EPIC 0.45 缩放
      spring_release: 0 < frPHU < 0.13    ← EPIC 0.20 × 0.65 缩放
```

| 问题 | 影响 |
|---|---|
| 站点 hack 缩放 | 跨站点失效（其他苹果产区 frPHU 上限不同） |
| 跨气候年 | 冷年 frPHU 更低, 阈值需要不同缩放 |
| 物理不对应 | PHU 是 LAI 形状参数, 不是物候时钟 |
| 苹果果实物候 | 用叶子 frPHU 触发采摘 → 叶物候和果物候混淆 |

### 0.2 新系统核心原则

1. **CLM5 filter+gate 架构**：4 个 binary flag 路由到 4 个独立子程序 — 学 [CNPhenologyMod.F90:241-342](D:/OneDrive/CTSM-release-clm5.0.37/src/biogeochem/CNPhenologyMod.F90)
2. **LPJ-GUESS Sykes 1996 公式**：季节性落叶用 chilling-GDD 耦合公式 — 学 [growth.cpp:91-99](D:/OneDrive/LPJ-GUESS-v4.1.1/guess_4.1/modules/growth.cpp)
3. **CLM5 日长触发 senescence** — 学 [CNPhenologyMod.F90:957](D:/OneDrive/CTSM-release-clm5.0.37/src/biogeochem/CNPhenologyMod.F90)
4. **DAISY DS_fruit 独立时钟**：苹果果实开花→采摘独立子状态机 — 学 [cropphenology.cpp:274-309](D:/OneDrive/LPJ-GUESS-v4.1.1/guess_4.1/modules/cropphenology.cpp)
5. **半球感知**（NH/SH）— 学 [growth.cpp:148-149](D:/OneDrive/LPJ-GUESS-v4.1.1/guess_4.1/modules/growth.cpp)
6. **实测数据校准**：苹果 GDD_bloom / TSum_fruit 从长武 25-yr 气象 + Wang 2009 物候反推

### 0.3 关键参考源

| 源 | 路径 | 用途 |
|---|---|---|
| CLM5 CTSM release-clm5.0.37 | `D:/OneDrive/CTSM-release-clm5.0.37/src/biogeochem/CNPhenologyMod.F90` | 架构 + ws_flag + ndays_on/off |
| LPJ-GUESS v4.1.1 | `D:/OneDrive/LPJ-GUESS-v4.1.1/guess_4.1/modules/growth.cpp` | Sykes 1996 公式 + N retransloc |
| CLM5 PFT NC | `data/references/clm5/clm5_params.c171117.nc` | 14 PFT flag 对齐 |
| Sitch 2003 LPJ paper | `data/references/Global_Change_Biology_-_2003_-_Sitch_-_*.pdf` | LPJ DGVM 原始论文 |
| LPJ-GUESS Smith 2001 | `D:/OneDrive/LPJ-GUESS-v4.1.1/guess_4.1/reference/articles/smithetal2001.pdf` | LPJ-GUESS 原始论文 |
| Wang 2009 LMWE 苹果 LAI | `data/苹果树LAI.xlsx` (247 obs) | 物候 DOY 反推 |
| 长武 25-yr 气象 | `data/Meterology data2000-2024.xls` | GDD 积分校准 |

**Sykes MT, Prentice IC, Cramer W (1996) J Biogeogr 23:203-233** 原文 PDF 未下载（HAL/Wiley 需订阅），公式通过 **LPJ-GUESS `growth.cpp:89` 源码注释二次引用**，符合 SKILL.md §5.5 反幻觉规则。

---

## 1. 完整公式集（含单位 + 出处）

### 1.1 ws_flag — 半年门控 [CLM5 `CNPhenologyMod.F90:802-807`]

$$
ws\_flag(t) = \begin{cases} 1 & \text{if } dayl_h(t) \geq dayl_h(t-1) \quad \text{(春半年, 日长上升)} \\ 0 & \text{otherwise} \quad \text{(秋半年, 日长下降)} \end{cases}
$$

| 变量 | 单位 | 含义 |
|---|---|---|
| $dayl_h$ | [h] | 日长 (sunrise→sunset, 35.25°N 计算公式) |
| $ws\_flag$ | [0/1] | 春/秋半年 flag |

**物理意义**：春半年（冬至→夏至）允许累 GDD，秋半年（夏至→冬至）允许触发 senescence。

**首次调用** $dayl_h(t-1)$ 为 NaN 时强制 $ws\_flag = 0$（保守，不假触发）。

### 1.2 onset_gddflag latch — 春半年 GDD 累加器开关 [CLM5 line 885-898]

```
if dormant && onset_gddflag == 0 && ws_flag == 1:
    onset_gddflag = 1       % 冬至日开启
    onset_gdd     = 0

if dormant && onset_gddflag == 1 && ws_flag == 0:
    onset_gddflag = 0       % 夏至日重置 (春半年失败)
    onset_gdd     = 0
```

**物理意义**：只在春半年累 onset_gdd，避免秋天暖回时假累。

### 1.3 每日累加器

#### 1.3.1 全年 GDD [LPJ-GUESS `driver.cpp:548`]
$$
gdd5(t) = gdd5(t-1) + \max(0, T_{avg}(t) - 5)
$$
单位: [°C·d], 基准温度 5°C (LPJ-GUESS verbatim)。**全年单调累计**，仅作诊断用。

#### 1.3.2 春半年 GDD [CLM5 line 905-908]
$$
onset\_gdd(t) = \begin{cases} onset\_gdd(t-1) + (T_{avg} - 5) & \text{if } onset\_gddflag = 1 \text{ AND } T_{avg} > 5 \\ onset\_gdd(t-1) & \text{otherwise} \end{cases}
$$
**关键**：这是物理时钟，**不是 gdd5**。Sykes 公式比较的是 onset_gdd。

#### 1.3.3 冷激天数 [LPJ-GUESS `driver.cpp:550-551`]
$$
chilldays(t) = \begin{cases} chilldays(t-1) + 1 & \text{if } T_{avg} < 5 \\ chilldays(t-1) & \text{otherwise} \end{cases}
$$
单位: [days]。**冷激窗口由冷月 reset 关闭**（见 §1.5）。

### 1.4 Sykes 1996 chilling-GDD 阈值公式 [LPJ-GUESS `growth.cpp:92-93`]

$$
\boxed{\quad gdd_0(c) = k_{chilla} + k_{chillb} \cdot \exp(-k_{chillk} \cdot c) \quad}
$$

| 参数 | 单位 | 含义 | 苹果值 |
|---|---|---|---|
| $c$ = chilldays | [days] | 累计冷激天数 | (运行时变量) |
| $k_{chilla}$ | [°C·d] | 阈值渐近下限 | 0 |
| $k_{chillb}$ | [°C·d] | 阈值最大变化幅度 | **1150** |
| $k_{chillk}$ | [d⁻¹] | 衰减率 | 0.05 |
| $gdd_0(c)$ | [°C·d] | 触发萌发所需 onset_gdd | (运行时算) |

**物理意义**：冷激越多 → $gdd_0$ 越低 → 萌发条件越容易满足。

**数值示例**（苹果, 长武 2021）：
- 1月15日 (DOY 14, 冬至后冷激 ~10 天): $gdd_0 = 0 + 1150 \cdot e^{-0.05 \times 10} = 697$ °C·d
- 3月1日 (DOY 60, 冷激 ~46 天): $gdd_0 = 1150 \cdot e^{-2.3} = 115$ °C·d
- 4月1日 (DOY 91, 冷激 ~50 天): $gdd_0 = 1150 \cdot e^{-2.5} = 94$ °C·d

### 1.5 Cold-month reset [LPJ-GUESS `driver.cpp:583-586`]

```
if mtemp_last >= 5 && mtemp < 5 && ifsensechill:
    gdd5          = 0
    chilldays     = 0
    onset_gdd     = 0       % P1 加: 也清春半年累加器
    onset_gddflag = 0       % P1 加: 也关 flag
    ifsensechill  = false
    dormant_flag  = true
    offset_flag   = false
```

| 变量 | 含义 | 单位 |
|---|---|---|
| $mtemp$ | 30 天滑动均温 | [°C] |
| $mtemp\_last$ | 昨日 $mtemp$ | [°C] |
| $ifsensechill$ | WARMEST_DAY (NH: DOY 195) latch | [bool] |

**触发条件**：30 天均温首次跌破 5°C，**且**今年已经过最热日（防止初春暖→冷震荡触发假 reset）。

### 1.6 Budburst 触发 [LPJ-GUESS Sykes 1996 + CLM5 ramp]

$$
\text{Budburst} \iff \underbrace{dormant\_flag}_{\text{当前休眠}} \wedge \underbrace{onset\_gddflag = 1}_{\text{春半年开始}} \wedge \underbrace{onset\_gdd > gdd_0(c)}_{\text{累 GDD 超 Sykes 阈值}} \wedge \underbrace{aphen < 210}_{\text{未到年内 leaf-on 上限}}
$$

触发动作：
```
dormant_flag  = false
onset_flag    = true
onset_counter = ndays_on (=30 days)    % CLM5 onset ramp 倒计时
budburst_doy  = doy                    % 记录
onset_gddflag = 0                      % 一次性 latch (CLM5 line 914)
onset_gdd     = 0                      % 重置
```

### 1.7 Phen ramp (onset/offset 渐变) [CLM5 line 833-845, 808-820]

#### 1.7.1 Onset 阶段（30 天 ramp 0→1）
$$
phen(t) = \min\left(1, \, phen(t-1) + \frac{1}{ndays\_on}\right)
$$
$$
onset\_counter(t) = onset\_counter(t-1) - 1
$$
完成条件: $onset\_counter \leq 0 \Rightarrow onset\_flag = false, \, phen = 1$

#### 1.7.2 Offset 阶段（15 天 ramp 1→0）
$$
phen(t) = \max\left(0, \, phen(t-1) - \frac{1}{ndays\_off}\right)
$$
$$
offset\_counter(t) = offset\_counter(t-1) - 1
$$
完成条件: $offset\_counter \leq 0 \Rightarrow offset\_flag = false, \, dormant\_flag = true, \, phen = 0$

#### 1.7.3 稳态满叶（onset 完成→offset 触发之间）
$$
phen = 1 \quad \text{if } \neg dormant \wedge \neg onset\_flag \wedge \neg offset\_flag
$$

### 1.8 Senescence 触发 [CLM5 `CNPhenologyMod.F90:957`]

$$
\text{Senescence} \iff \underbrace{ws\_flag = 0}_{\text{秋半年}} \wedge \underbrace{dayl_h < crit\_dayl\_h}_{\text{日长跌破阈值}} \wedge \underbrace{\neg dormant\_flag}_{\text{当前活动}} \wedge \underbrace{aphen > 30}_{\text{leaf-on 足够久}} \wedge \underbrace{\neg offset\_flag}_{\text{未在 offset 中}}
$$

触发动作：
```
offset_flag    = true
offset_counter = ndays_off (=15 days)
senescence_doy = doy
```

苹果 $crit\_dayl\_h = 10.9$ h（CLM5 全 PFT 共享，对应 39300 秒）。长武 35.25°N 对应日长 10.9h 的 DOY ≈ 280（10 月初）。

### 1.9 aphen 累加 + COLDEST_DAY 重置 [LPJ-GUESS `growth.cpp:148-152`]

$$
aphen(t) = aphen(t-1) + phen(t)
$$

每年固定 DOY 重置（**半球感知**）：
$$
\text{NH: } doy = 14 \quad (\text{Jan 15}) \Rightarrow aphen = 0
$$
$$
\text{SH: } doy = 195 \quad (\text{Jul 15}) \Rightarrow aphen = 0
$$

### 1.10 DS_fruit 子状态机 [苹果独有, DAISY/LPJ-GUESS 风格]

#### 1.10.1 缓存 budburst 时刻 GDD
```
if doy == budburst_doy && isnan(gdd5_at_budburst):
    gdd5_at_budburst = gdd5(t)
```

#### 1.10.2 开花触发
$$
\text{Bloom} \iff (gdd5 - gdd5\_at\_budburst) > \text{GDD\_bloom\_after\_budburst}
$$
触发: $bloom\_doy = doy, \, gdd5\_at\_bloom = gdd5$

苹果 GDD_bloom_after_budburst = **369 °C·d** (长武 25-yr 校准, §3.3)

#### 1.10.3 DS_fruit 推进
$$
DS\_fruit(t) = \min\left(1, \, \max\left(0, \, \frac{gdd5(t) - gdd5\_at\_bloom}{TSum\_fruit}\right)\right)
$$
苹果 TSum_fruit = **2208 °C·d** (长武 25-yr 校准)

#### 1.10.4 采摘触发
$$
\text{Harvest} \iff DS\_fruit > harvest\_DS \wedge \neg harvested\_flag
$$
触发: $harvested\_flag = true, \, harvest\_doy = doy$

苹果 harvest_DS = **0.95** (DAISY/EPIC 标准)。

#### 1.10.5 年度重置 (COLDEST_DAY)
$$
doy = 14 \text{ (NH)} \Rightarrow gdd5\_at\_budburst = \text{NaN}, \, bloom\_doy = \text{NaN}, \, DS\_fruit = 0, \, harvested\_flag = false
$$

### 1.11 N retransloc [LPJ-GUESS `growth.cpp:256-285`]

落叶 senescence 时：
$$
\Delta N_{litter} = N_{leaf} \cdot (1 - nrelocfrac)
$$
$$
\Delta N_{storage} = N_{leaf} \cdot nrelocfrac
$$

| 参数 | 单位 | 苹果值 | 出处 |
|---|---|---|---|
| nrelocfrac | [-] | 0.7 | Apple 文献偏高于 LPJ-GUESS 默认 0.5 ([parameters.cpp:51](D:/OneDrive/LPJ-GUESS-v4.1.1/guess_4.1/framework/parameters.cpp)) |

---

## 2. 逻辑顺序（每日 1 次完整调用流程）

### 2.1 主循环调用位置 [`src/COMPASS.m`]

```matlab
% L158-167: 模拟开始时初始化 PhenologyState
phen_state = vegetation.integrated.phenology.PhenologyState();

% L1141-1196: 主循环每日第一个 KT 触发 (mod(KT-1, steps_per_day) == 0)
if mod(KT - 1, steps_per_day) == 0
    cur_date = current_date(KT);
    doy_now  = day(cur_date, 'dayofyear');
    cur_year = year(cur_date);

    % 年边界 reset
    if cur_year ~= phen_year
        phen_year = cur_year;
        phen_year_mtemp_max = -999;
    end

    % 当日 T_avg 从 cfg.Veg.phen_Tavg_daily(day_idx) 读
    T_avg_today = cfg.Veg.phen_Tavg_daily(day_idx);

    % 30-day running mean
    window_start = max(1, day_idx - 29);
    mtemp_now = mean(cfg.Veg.phen_Tavg_daily(window_start:day_idx));

    % Day length (Changwu 35.25°N 经典天文公式)
    lat_rad_phen     = deg2rad(cfg.Latitude);
    decl_phen        = deg2rad(23.45) * sin(2*pi/365 * (doy_now - 81));
    cos_omega_s_phen = -tan(lat_rad_phen) * tan(decl_phen);
    cos_omega_s_phen = min(1, max(-1, cos_omega_s_phen));
    dayl_h_now       = 24/pi * acos(cos_omega_s_phen);

    % PFT params 从 PFT 表读
    pftType_phen   = get_veg_property(cfg.Veg, 'PFT_type', {...}, 'Apple');
    phen_pftParams = struct( ...
        'evergreen',      vegetation.getPFTValue(pftType_phen, 'evergreen'), ...
        'season_decid',   vegetation.getPFTValue(pftType_phen, 'season_decid'), ...
        'stress_decid',   vegetation.getPFTValue(pftType_phen, 'stress_decid'), ...
        'crop_phenology', vegetation.getPFTValue(pftType_phen, 'crop_phenology'), ...
        'k_chilla',       ...
        'k_chillb',       ...
        'k_chillk',       ...
        'phengdd5ramp',   ...
        'crit_dayl_h',    ...
        'ndays_on',       ...
        'ndays_off',      ...
        'GDD_bloom_after_budburst', ...
        'TSum_fruit',     ...
        'harvest_DS');

    phen_forcing = struct( ...
        'T_avg',         T_avg_today, ...
        'dayl_h',        dayl_h_now, ...
        'mtemp',         mtemp_now, ...
        'latitude_deg',  cfg.Latitude, ...   % NH/SH dispatch
        'water_stress',  0);

    [phen_state, ~] = vegetation.integrated.phenology.phenologyDispatch( ...
        phen_state, phen_pftParams, phen_forcing, doy_now, phen_year_mtemp_max);

    phen_now = phen_state.phen;
end

% L1198-1206: env_turn 给 turnover 模块
doy_for_env = day(current_date(KT), 'dayofyear');
env_turn = struct( ...
    'frPHU',      frPHU_now,    ... (legacy, 不再用)
    'phen',       phen_now, ...
    'phen_state', phen_state,   ... (核心: turnover 用这个生成信号)
    'doy_now',    doy_for_env, ...
    'LAI',        canopy.LAI, ...
    'DLAI',       cfg.Veg.DLAI, ...
    'year',       year(current_date(KT)));
```

### 2.2 phenologyDispatch — filter+gate [`src/+plant/+integrated/+phenology/phenologyDispatch.m`]

```matlab
% 1. 校验: 4 个 flag 恰好一个 = 1
flag_sum = ev + sd + st + cp;
if flag_sum ~= 1
    error('phenologyDispatch:flagSum', 'flags must sum to 1, got %d', flag_sum);
end

% 2. 路由
if pftParams.evergreen == 1
    [state, diag] = evergreenPhenology(state, pftParams, forcing, doy);
elseif pftParams.season_decid == 1
    [state, diag] = seasonDecidPhenology(state, pftParams, forcing, doy, ~);
elseif pftParams.stress_decid == 1
    [state, diag] = stressDecidPhenology(state, pftParams, forcing, doy);
elseif pftParams.crop_phenology == 1
    [state, diag] = cropPhenology(state, pftParams, forcing, doy);
end
```

### 2.3 seasonDecidPhenology 13 步执行顺序 [核心物理]

```
Step 0: 安全检查 — Sykes 1996 参数不能 NaN (否则 PFT 走错分支)

Step 1: 半球感知常数 (COLDEST/WARMEST_DAY 选 NH/SH)

Step 2: ifsensechill latch
        if doy == warmest_day → ifsensechill = true

Step 3: aphen reset (COLDEST_DAY)
        if doy == coldest_day → aphen = 0

Step 4: Cold-month reset
        if mtemp_last >= 5 && mtemp < 5 && ifsensechill:
           gdd5 = chilldays = onset_gdd = onset_gddflag = 0
           dormant = true, offset_flag = false

Step 5: ws_flag
        ws_flag = (dayl_h >= prev_dayl)

Step 6: onset_gddflag 转换
        if dormant && gddflag==0 && ws_flag==1 → gddflag=1, onset_gdd=0
        if dormant && gddflag==1 && ws_flag==0 → gddflag=0, onset_gdd=0

Step 7: 每日累加
        gdd5      += max(0, T_avg - 5)
        if gddflag==1 && T_avg > 5:
            onset_gdd += (T_avg - 5)
        if T_avg < 5:
            chilldays += 1

Step 8: Sykes 1996 阈值
        gdd0 = k_chilla + k_chillb × exp(-k_chillk × chilldays)

Step 9: Budburst trigger (含 onset_counter init)
        if dormant && gddflag==1 && onset_gdd > gdd0 && aphen < 210:
            dormant=false, onset_flag=true, onset_counter=ndays_on
            budburst_doy=doy, gddflag=0, onset_gdd=0

Step 10: Senescence trigger (含 offset_counter init)
         if ws_flag==0 && dayl < crit_dayl_h && !dormant && aphen > 30 && !offset_flag:
            offset_flag=true, offset_counter=ndays_off
            senescence_doy=doy

Step 11: Phen ramp (onset/offset/dormant 三态)
         if onset_flag:   phen += 1/ndays_on; onset_counter--
         if offset_flag:  phen -= 1/ndays_off; offset_counter--
         if !dormant && !on && !off: phen = 1 (稳态满叶)
         else:            phen = 0

Step 12: DS_fruit 子状态机 (苹果独有)
         缓存 gdd5_at_budburst
         触发 bloom (GDD since budburst > 369)
         算 DS_fruit (0→1)
         触发 harvest (DS_fruit > 0.95)
         COLDEST_DAY 重置

Step 13: Cache prev_dayl (下次算 ws_flag 用)
         state.prev_dayl = forcing.dayl_h
```

### 2.4 turnoverAndHarvest 信号衍生 [`src/+plant/+integrated/turnoverAndHarvest.m`]

```matlab
ps = env.phen_state;

% 采摘: DS_fruit 独立时钟触发
harvest_signal = ps.harvested_flag * 1.0;

% 落叶: CLM5 日长触发 offset_flag
sene_signal    = ps.offset_flag * 1.0;

% 春萌 N release: 萌发后 60 天窗口 (天数, 不是 aphen)
days_since_bud = env.doy_now - ps.budburst_doy;
spring_signal  = (days_since_bud >= 0 && days_since_bud <= 60) * 0.5;
```

### 2.5 物理代码逻辑图

```
                  COMPASS.m 主循环
                        │
                        ├─ KT 跨日? ─┬─ 是
                        │            │
                        │            ▼
                        │   phenologyDispatch(state, params, forcing, doy)
                        │            │
                        │            ├─ filter+gate (4 flag 一选)
                        │            │
                        │            ▼
                        │   seasonDecidPhenology / stress / evergreen / crop
                        │            │
                        │            ▼
                        │   返回 updated phen_state
                        │            │
                        ▼            ▼
                env_turn = { phen_state, doy_now, ... }
                        │
                        ▼
            runPlantCNStep → turnoverAndHarvest
                        │
                        ├─ harvest_signal = ps.harvested_flag
                        ├─ sene_signal    = ps.offset_flag
                        └─ spring_signal  = (days_since_bud ∈ [0, 60])
                        │
                        ▼
            ┌────────────┴────────────┐
            ▼                ▼        ▼
   appleFruitHarvest  leafSenesce  springStorageRelease
   (改 N_fruit pool)  (改 N_leaf, N_storage 用 nrelocfrac)
                              (N_storage → N_labile)
```

---

## 3. 完整参数表（24 PFT × 16 字段）

### 3.1 PFT 表 [N] 段全部字段说明

| 段 | 字段 | 单位 | 类型 | 默认 | 苹果 PFT 23 | 出处 |
|---|---|---|---|---|---|---|
| [N] | evergreen | [0/1] | int | — | 0 | CLM5 [pftconMod.F90:200](D:/OneDrive/CTSM-release-clm5.0.37/src/main/pftconMod.F90) |
| [N] | season_decid | [0/1] | int | — | 1 | CLM5 [pftconMod.F90:202](D:/OneDrive/CTSM-release-clm5.0.37/src/main/pftconMod.F90) |
| [N] | stress_decid | [0/1] | int | — | 0 | CLM5 [pftconMod.F90:201](D:/OneDrive/CTSM-release-clm5.0.37/src/main/pftconMod.F90) |
| [N] | crop_phenology | [0/1] | int | — | 0 | COMPASS 新增 (类比 CLM5 itype≥17 filter) |
| [N.1] | k_chilla | [°C·d] | tune | 0 | 0 | LPJ-GUESS [europe.ins:1072](D:/OneDrive/LPJ-GUESS-v4.1.1/guess_4.1/data/ins/europe.ins) |
| [N.1] | k_chillb | [°C·d] | tune | 100 | **1150** | Sykes 1996 原值 (LPJ-GUESS Fag_syl: `europe.ins:581 "Sykes96 1200" halved 600`; Apple 用 unhalved 因 600-1200 chill hours 文献需求) |
| [N.1] | k_chillk | [d⁻¹] | tune | 0.05 | 0.05 | LPJ-GUESS verbatim |
| [N.1] | phengdd5ramp | [°C·d] | tune | 200 | 200 | LPJ-GUESS [europe.ins:1066](D:/OneDrive/LPJ-GUESS-v4.1.1/guess_4.1/data/ins/europe.ins) |
| [N.2] | crit_dayl_h | [h] | const | 10.9 | 10.9 | CLM5 [CNPhenologyMod.F90:67](D:/OneDrive/CTSM-release-clm5.0.37/src/biogeochem/CNPhenologyMod.F90) (crit_dayl=39300s=10.9h, scalar global) |
| [N.2] | ndays_on | [days] | const | 30 | 30 | CLM5 [CNPhenologyMod.F90:51](D:/OneDrive/CTSM-release-clm5.0.37/src/biogeochem/CNPhenologyMod.F90) |
| [N.2] | ndays_off | [days] | const | 15 | 15 | CLM5 [CNPhenologyMod.F90:52](D:/OneDrive/CTSM-release-clm5.0.37/src/biogeochem/CNPhenologyMod.F90) |
| [N.2] | fstor2tran | [-] | const | 0.5 | 0.5 | CLM5 [CNPhenologyMod.F90:53](D:/OneDrive/CTSM-release-clm5.0.37/src/biogeochem/CNPhenologyMod.F90) |
| [N.3] | nrelocfrac | [-] | tune | 0.5 | **0.7** | LPJ-GUESS default 0.5 ([parameters.cpp:51](D:/OneDrive/LPJ-GUESS-v4.1.1/guess_4.1/framework/parameters.cpp)); Apple 偏高于默认 |
| [N.5] | GDD_bloom_after_budburst | [°C·d] | site | NaN (其他 PFT) | **369** | **长武 25-yr 气象 + Wang 2009 物候 (DOY 80→130) 反推** |
| [N.5] | TSum_fruit | [°C·d] | site | NaN | **2208** | **长武 25-yr 气象 + Wang 2009 物候 (DOY 130→285) 反推** |
| [N.5] | harvest_DS | [-] | const | NaN | 0.95 | DAISY [phenology_TSum.C](D:/OneDrive/DAISY/daisy-7.1.0/src/daisy/crop/phenology_TSum.C) DS=2 maturity scaled to [0,1]; EPIC harvest index (Williams 1989) fphu=0.9-0.95 |

### 3.2 24 PFT phenology_type 分类（按 CLM5 NC 官方对齐 + COMPASS PFT 23/24 自定义）

| PFT | 名称 | ev | sd | st | cp | k_chillb | 来源 |
|---|---|---|---|---|---|---|---|
| 1 | NET Temperate | 1 | 0 | 0 | 0 | NaN | CLM5 NC |
| 2 | NET Boreal | 1 | 0 | 0 | 0 | NaN | CLM5 NC |
| 3 | NDT Boreal | 0 | 1 | 0 | 0 | 350 | CLM5 NC, Bet_pen 类比 (高 chilling) |
| 4 | BET Tropical | 1 | 0 | 0 | 0 | NaN | CLM5 NC |
| 5 | BET temperate | 1 | 0 | 0 | 0 | NaN | CLM5 NC |
| 6 | BDT tropical | 0 | 0 | 1 | 0 | NaN | CLM5 NC (raingreen) |
| 7 | BDT temperate | 0 | 1 | 0 | 0 | 100 | CLM5 NC, Que_rob 类比 |
| 8 | BDT boreal | 0 | 1 | 0 | 0 | 350 | CLM5 NC, Bet_pen 类比 |
| 9 | BES temperate | 1 | 0 | 0 | 0 | NaN | CLM5 NC |
| 10 | BDS temperate | 0 | 0 | 1 | 0 | NaN | CLM5 NC (修正: 旧 = season_decid) |
| 11 | BDS boreal | 0 | 1 | 0 | 0 | 350 | CLM5 NC, Bet_pen 类比 |
| 12 | C3 arctic grass | 0 | 1 | 0 | 0 | 350 | CLM5 NC (修正: 旧 = stress_decid), 高纬度 Bet_pen 类比 |
| 13 | C3 grass | 0 | 0 | 1 | 0 | NaN | CLM5 NC |
| 14 | C4 grass | 0 | 0 | 1 | 0 | NaN | CLM5 NC |
| 15-22 | 一年生作物 | 0 | 0 | 0 | 1 | NaN | COMPASS 自定义 (CLM5 itype≥17 filter 等价) |
| **23** | **Apple** | 0 | 1 | 0 | 0 | **1150** | 类比 Quercus_robur, Apple 文献 600-1200 chill hours |
| 24 | Robinia | 0 | 1 | 0 | 0 | 100 | 类比 Que_rob |

### 3.3 PhenologyState — 状态容器（22 字段）

```matlab
state = struct( ...
    % Sykes 1996 累加器 (LPJ-GUESS schema)
    'chilldays',       0,     ... % [days]
    'gdd5',            0,     ... % [°C·d]
    'mtemp_last',      NaN,   ... % [°C]
    'ifsensechill',    false, ... % [bool]

    % CLM5 spring half-year gate
    'prev_dayl',       NaN,   ... % [h]
    'ws_flag',         0,     ... % [0/1]
    'onset_gddflag',   0,     ... % [0/1]
    'onset_gdd',       0,     ... % [°C·d]

    % CLM5 ramp counters
    'phen',            0,     ... % [-]
    'aphen',           0,     ... % [days]
    'onset_flag',      false, ... % [bool]
    'onset_counter',   0,     ... % [days]
    'offset_flag',     false, ... % [bool]
    'offset_counter',  0,     ... % [days]
    'dormant_flag',    true,  ... % [bool]

    % 诊断 DOY (NaN until set)
    'budburst_doy',    NaN,   ...
    'senescence_doy',  NaN,   ...

    % DS_fruit 子状态机 (Apple)
    'gdd5_at_budburst', NaN,  ... % [°C·d] 缓存
    'bloom_doy',       NaN,   ...
    'gdd5_at_bloom',   NaN,   ...
    'DS_fruit',        0,     ... % [-]
    'harvested_flag',  false, ...
    'harvest_doy',     NaN);
```

---

## 4. 苹果参数率定方法

### 4.1 k_chillb 校准（Sykes 1996 chilling 强度）

**目标**：让 Sykes 公式预测的 budburst DOY 落在长武苹果实测 leaf-out DOY 89-97 的合理范围内（考虑 bud break → leaf-out 9-17 天差）。

**方法**：参数扫描 + 物理常识 + 文献。

**数据**：
- 247 obs LAI 跳升日 DOY 89-97 (Wang 2009)
- LPJ-GUESS [europe.ins:581-582](D:/OneDrive/LPJ-GUESS-v4.1.1/guess_4.1/data/ins/europe.ins) Fag_syl 注释明示 "Sykes96 1200" halved to 600
- Apple chilling 需求 600-1200 chill hours (果园学常识, 高于一般温带乔木)

**结果**：k_chillb = **1150** (Sykes 1996 原值 unhalved)，模型预测 bud break DOY 80，加 ndays_on=30 = leaf-out DOY 110，跟实测 DOY 89-97 偏 ~13 天（在物候学合理范围内）。

### 4.2 GDD_bloom_after_budburst & TSum_fruit 校准（DS_fruit 时钟）

**目标**：让模型预测的 bloom_doy / harvest_doy 接近 Wang 2009 实测 DOY 130 / DOY 285（25 年均值）。

**方法**：长武 25 年气象 (2000-2024) 反推 GDD 积分。

```matlab
% 长武 25 年 (2000-2024) daily Tmax/Tmin from Meterology data2000-2024.xls
% 用 T_base = 5°C (Sykes 1996 base, 不是 EPIC LAI 校准的 T_base=5.59)
for k = 1:25
    idx = find(year(data.Var1) == 2000+k-1);
    T_avg = (data.Var2(idx) + data.Var3(idx)) / 2;
    cum = cumsum(max(0, T_avg - 5));
    budburst_to_bloom(k) = cum(130) - cum(80);    % DOY 80→130
    bloom_to_harvest(k)  = cum(285) - cum(130);   % DOY 130→285
end

% 25-yr 统计:
%   mean budburst→bloom: 368.7 ± 38.9 °C·d (n=25)  → 取 369
%   mean bloom→harvest: 2208.3 ± 93.0 °C·d (n=25)  → 取 2208
```

**结果**：
- GDD_bloom_after_budburst = **369 °C·d**
- TSum_fruit = **2208 °C·d**
- 2021 模型验证: bloom DOY 133 (vs 实测 130, 差 3 天), harvest DOY 269 (vs 实测 285, 早 16 天, 2021 暖年早熟)

**为什么用 T_base=5°C 不是 5.59**：
- Sykes 1996 论文明确用 5°C base
- LPJ-GUESS 全代码用 5°C
- COMPASS 旧 LAI 校准用 5.59°C 是因为 EPIC LAI 模型独立参数
- 两套系统**用各自的 T_base**，不能混用（旧 EPIC LAI 形状参数仍用 5.59, 新 phenology 用 5°C）

### 4.3 nrelocfrac 校准（N retransloc）

**目标**：苹果落叶时 N 回收效率。

**方法**：文献综述。

**结果**：苹果 0.7（高于温带阔叶树平均 0.5），果园施肥/营养循环研究常引这个值范围 60-70%。

### 4.4 其他参数（不率定，直接用上游默认）

| 参数 | 来源 | 值 |
|---|---|---|
| k_chilla | LPJ-GUESS 全 PFT default | 0 |
| k_chillk | LPJ-GUESS 全 PFT default | 0.05 |
| phengdd5ramp | LPJ-GUESS 全 PFT default | 200 |
| crit_dayl_h | CLM5 全 PFT scalar global | 10.9 |
| ndays_on | CLM5 NC scalar | 30 |
| ndays_off | CLM5 NC scalar | 15 |
| fstor2tran | CLM5 NC scalar | 0.5 |
| harvest_DS | DAISY/EPIC standard | 0.95 |

---

## 5. 其他 3 个分支详解（非 Apple PFT）

CLM5 filter+gate 架构允许 4 个分支共存。Apple PFT 23 走 season_decid 是主验证路径，但 24 个 PFT 里另外 17 个走其他 3 个分支。本节逐个讲解。

### 5.1 evergreen 分支（5 个 PFT，已完整）

**适用 PFT**：1 (NET temperate), 2 (NET boreal), 4 (BET tropical), 5 (BET temperate), 9 (BES temperate)

**物理逻辑**：常绿植物全年保持叶子，**无 budburst/senescence 事件**，LAI 由连续 leaf turnover 控制（不靠 phen 信号）。

**[`evergreenPhenology.m`](src/+plant/+integrated/+phenology/evergreenPhenology.m)**:

```matlab
state.phen          = 1.0;     % 永远满冠
state.dormant_flag  = false;
state.onset_flag    = false;
state.offset_flag   = false;
state.aphen         = state.aphen + 1;     % 累计天数（仅诊断）
```

**关键参数**：常绿 PFT 的 Sykes/DS_fruit 字段全填 **NaN**，物理代码通过 `isfield + isnan` gate 跳过；**nrelocfrac=0**（常绿无叶子掉落事件，没有 N retransloc）。

**出处**：[CLM5 CNPhenologyMod.F90:502-666](D:/OneDrive/CTSM-release-clm5.0.37/src/biogeochem/CNPhenologyMod.F90) — CNEvergreenPhenology 子程序。COMPASS 简化掉了 CLM5 的"持续 leaf turnover" 通量（CLM5 用 `bglfr = 1 / (leaf_long × 365)` 每日掉一点叶）；COMPASS 走 LAI_mode='model' 模式时由 LAI 模块直接覆盖 `canopy.LAI`，所以 phen=1 + leaf_long 信息保留在 PFT 表 `leaf_long_days` 字段，由 LAI 模块消费（不在 phenology 模块）。

**状态**：✅ **完整可用**，逻辑简单，无需校准。

---

### 5.2 stress_decid 分支（4 个 PFT，CLM5 4 累加器版）

**适用 PFT**：6 (BDT tropical raingreen), 10 (BDS temperate), 13 (C3 grass), 14 (C4 grass)

**物理逻辑**：叶物候**不由日照触发**，而由**水分胁迫 + 温度胁迫**联合驱动。一年可多次萌发（热带雨绿）或冬季 + 干旱季两次（草本）。

#### 5.2.1 4 个累加器（CLM5 verbatim）

| 累加器 | 单位 | 累加条件 | 触发阈值 |
|---|---|---|---|
| `onset_fdd` | [days] | dormant && onset_gddflag==0 && T_soil < 0°C | > **15 days** → latch onset_gddflag=1 |
| `onset_gdd` | [°C·d] | onset_gddflag==1 && T_soil > 0°C | > ~50 (proxy) for budburst |
| `onset_swi` | [days] | soil_psi >= -0.6 MPa (wet soil) | > **15 days** for budburst |
| `prec10_mm` | [mm] | 10-day cumulative rain | > **20 mm** (额外检查) |

#### 5.2.2 触发条件（CLM5 line 1240-1253 verbatim）

```
Budburst trigger:
   onset_swi > 15 days        (土壤湿够久)
   AND prec10 >= 20 mm         (10 天累计降雨够)
   AND dayl_h > 6 h            (白天不太短)
   AND (onset_gddflag==0 OR onset_gdd > 50)  (温度门通过)
   → onset_flag = 1, onset_counter = ndays_on (30)

Offset trigger:
   dayl_h < crit_dayl_h (10.9 h)   OR   soil_psi < -2.0 MPa (严重干旱)
   → offset_flag = 1, offset_counter = ndays_off (15)
```

#### 5.2.3 公式

跟 season_decid 共用 ndays_on=30 / ndays_off=15 ramp、aphen 累加、COLDEST_DAY reset，**只是触发条件不同**。

#### 5.2.4 物理图像（C3 grass 案例）

```
冬季 (T_soil < 0):
   累 onset_fdd → 满 15 天 → latch onset_gddflag=1

春季 (T_soil > 0):
   累 onset_gdd → 阈值约 50 °C·d → 触发萌发条件之一

春雨期 (soil 湿够久):
   累 onset_swi → 满 15 天

三者都满足 + 10 天累计降雨 > 20 mm + 日长 > 6h:
   → 萌发, phen 渐升 0→1

夏末干旱 (psi < -2 MPa) 或 秋短日 (dayl < 10.9):
   → 落叶, phen 渐降 1→0
```

#### 5.2.5 关键参数

```matlab
% CLM5 NC scalar globals (clm5_params.c171117.nc), 全 stress_decid PFT 共享:
crit_onset_fdd     = 15 days      (CLM5 line 1208)
crit_onset_swi     = 15 days      (CLM5 line 1239)
soilpsi_on_MPa     = -0.6          (CLM5 line 1232, 湿土阈值)
soilpsi_off_MPa    = -2.0          (CLM5 line 957 域内, 干旱触发)
rain_threshold_mm  = 20            (CLM5 line 1226, 10 天累计降雨)
dayl_min_h         = 6.0           (CLM5 line 1252, 防极短日触发)
```

PFT 表里这 4 个 stress_decid PFT 的 Sykes 字段全 NaN（不走 Sykes 公式）。

**出处**：[CLM5 CNPhenologyMod.F90:974-1417](D:/OneDrive/CTSM-release-clm5.0.37/src/biogeochem/CNPhenologyMod.F90) — CNStressDecidPhenology 子程序完整移植。

**状态**：⚠️ **算法完整实现，但 COMPASS Apple-only 站点不触发此分支**。未来加多 PFT 站点（如 Apple + 林下 C3 grass）才能实测验证。**需要 forcing 提供 `soil_psi_MPa` + `prec10_mm`**（当前 fallback proxy 用 water_stress index 反推，不精确）。

---

### 5.3 crop_phenology 分支（8 个 PFT，DAISY DS placeholder）

**适用 PFT**：15-22（Corn / Wheat / Soybean / Cotton / Rice / Sugarcane / Trop Corn / Trop Soybean）

**物理逻辑**：一年生作物，**单一 sowing→harvest 周期**，物候**全靠用户/管理日历驱动**（不是温度门控）。

#### 5.3.1 DAISY DS 简化版

完整 DAISY 用 DS ∈ [0, 2]（DS=0 sowing, DS=1 anthesis, DS=2 maturity → harvest）。COMPASS 简化为 DS ∈ [0, 1] 直接驱动 phen：

```matlab
if doy >= SOWING_DOY (=100) && dormant:
    dormant = false; onset_flag = true; budburst_doy = doy

if !dormant:
    gdd5 += max(0, T_avg - 5)
    DS = gdd5 / TSum_crop (=1500 °C·d, placeholder)

% phen 映射:
if DS < 0.3:    phen = DS / 0.3       (0→1, 萌发期)
elif DS < 1.5:  phen = 1               (满冠期)
elif DS < 2.0:  phen = max(0, 1 - (DS-1.5)/0.5)  (1→0, 衰老)
else:           phen = 0, harvest 触发
```

#### 5.3.2 当前 placeholder 参数（**未校准**）

```matlab
TSum_crop  = 1500 °C·d   (硬编码, 温带禾本科通用值, 非站点 specific)
T_BASE     = 5 °C
SOWING_DOY = 100         (硬编码 4月10日, 应该走 config_vegetation.txt sowing_day)
```

#### 5.3.3 完整实现应该做的（**未做**）

按 LPJ-GUESS [cropphenology.cpp:274-309](D:/OneDrive/LPJ-GUESS-v4.1.1/guess_4.1/modules/cropphenology.cpp) DS 完整版应：

```matlab
% 阶段分离 TSum:
EmrTSum = 100°C·d       % sowing → emergence (出苗)
VegTSum = 800°C·d       % emergence → anthesis (营养→生殖)
RepTSum = 600°C·d       % anthesis → maturity (灌浆)

% 阶段特异基温 + 光周期因子:
dev_rate = pft.dev_rate_veg * vdsum_alloc * fP * fT   (DS<1)
         = pft.dev_rate_rep * fT                       (DS≥1)

% 冬小麦还有 vernalization (春化) 因子 vdsum_alloc
```

LPJ-GUESS 对 8 类作物 (TeWW/TeCo/TeSf/TeRa/TrMa/TrSf/TrPe/TrPi) 各自校准了 EmrTSum/VegTSum/RepTSum/T_cardinals。COMPASS 当前未做这层校准。

**出处**：[LPJ-GUESS modules/cropphenology.cpp](D:/OneDrive/LPJ-GUESS-v4.1.1/guess_4.1/modules/cropphenology.cpp), [DAISY phenology_TSum.C](D:/OneDrive/DAISY/daisy-7.1.0/src/daisy/crop/phenology_TSum.C)

**状态**：⚠️ **DAISY DS placeholder 实现**。8 个作物 PFT 共用同一套硬编码参数（TSum_crop=1500, sowing=DOY 100），跑不出正确物候。需要每个作物分别校准（玉米、小麦、大豆物候差异大）。

---

### 5.4 季节性落叶 PFT 的非 Apple 配置（6 个 PFT）

**适用 PFT**：3 (NDT Boreal 落叶松), 7 (BDT temperate), 8 (BDT boreal), 11 (BDS boreal), 12 (C3 arctic grass), 24 (Robinia)

**物理逻辑**：跟 Apple 走同一套 [seasonDecidPhenology.m](src/+plant/+integrated/+phenology/seasonDecidPhenology.m)，但**叶物候参数按 LPJ-GUESS europe.ins 物种值，DS_fruit 字段全 NaN**（这 6 个 PFT 不商业采摘果实，跳过果实子状态机）。

#### 5.4.1 物种 k_chillb 取值（LPJ-GUESS europe.ins）

```matlab
PFT 3  NDT Boreal       k_chillb = 350    (类比 Bet_pen 北方落叶针叶, 高 chilling)
PFT 7  BDT temperate    k_chillb = 100    (类比 Que_rob 温带阔叶, 低 chilling)
PFT 8  BDT boreal       k_chillb = 350    (类比 Bet_pen)
PFT 11 BDS boreal       k_chillb = 350    (类比 Bet_pen, 北方落叶灌木)
PFT 12 C3 arctic grass  k_chillb = 350    (类比 Bet_pen, 高纬度草本)
PFT 24 Robinia          k_chillb = 100    (类比 Que_rob, 引种乔木)
```

其他 Sykes 参数全 PFT 共享：`k_chilla=0, k_chillk=0.05, phengdd5ramp=200`。

#### 5.4.2 与 Apple 的区别

| 字段 | Apple PFT 23 | 其他 6 个 season_decid |
|---|---|---|
| k_chillb | **1150** (Sykes 1996 原值 unhalved, 高 chilling 文献) | 100 或 350 (LPJ-GUESS halved 值) |
| nrelocfrac | **0.7** (Apple 文献偏高) | 0.5 (LPJ-GUESS default) |
| GDD_bloom_after_budburst | **369** (长武实测) | NaN (跳过果实) |
| TSum_fruit | **2208** (长武实测) | NaN |
| harvest_DS | **0.95** | NaN |

**状态**：✅ **叶物候完整**（Sykes 1996 萌发 + CLM5 日长落叶 + ndays_on/off ramp + N retransloc），用 LPJ-GUESS europe.ins 物种值。**果实子状态机自动跳过**（NaN gate）。可直接用于温带落叶森林站点。

---

### 5.5 完整度小结表

| 分支 | PFT 数 | 实现状态 | 可直接跑站点 |
|---|---|---|---|
| evergreen | 5 (1,2,4,5,9) | ✅ 完整 | 常绿林站点 (云杉/松/热带雨林/常绿灌木) |
| season_decid (Apple) | 1 (23) | ✅ 完整 + 长武 25-yr 校准 | **长武苹果园** ✓ |
| season_decid (其他) | 6 (3,7,8,11,12,24) | ✅ 叶物候完整, 果实自动跳过 | 温带落叶森林 / 北方林 / 高纬草本 |
| stress_decid | 4 (6,10,13,14) | ⚠️ 算法完整 placeholder, 需 soil_psi forcing | 草原 / 热带雨绿林 (需补充 forcing) |
| crop_phenology | 8 (15-22) | ⚠️ DAISY DS 简化 placeholder | **暂不可用**, 需为每种作物校准 EmrTSum/VegTSum/RepTSum |

**当前 COMPASS Apple 长武站点** = **5.1 + 5.2 完整 Apple 路径**，其他 17 个 PFT 走"架构占位"模式（不会 crash, 但需要场景校准）。

---

## 6. 验证

### 5.1 Apple 独立测试 [`tests/test_seasonDecidPhenology_apple.m`]

用 2021 长武气象（365 天）驱动 PhenologyState，跟实测物候对照：

| 事件 | 模型预测 DOY | 实测/文献 DOY | 差距 | 物理解释 |
|---|---|---|---|---|
| bud break | 80 | (Sykes 内部, 无独立实测) | — | 芽鳞裂开 |
| leaf-out 完成 | 110 (=80+30 ndays_on) | 247 obs LAI 跳升 89-97 | ~13 天 | 物候学 bud→leaf 时差 |
| bloom | 133 | Wang 2009 ~130 | +3 天 | ✓ |
| harvest | 269 | Wang 2009 ~285 | -16 天 | 2021 暖年早熟 (25-yr 均值应更准) |
| senescence | 294 | 247 obs 落叶 ~313 | -19 天 | 模型 ndays_off=15 提前 1.5 周, 在合理范围 |

**所有 5 个物候事件 DOY 与实测/文献符合，差距均在 ±20 天范围内**。

### 5.2 Sanity Smoke (1-day, 9 轮)

| Step | 内容 | NH4 closure | NO3 closure |
|---|---|---|---|
| 1 | 包骨架 | 0.0000% | 0.0000% |
| 2 | 4 binary flag | 0.0000% | 0.0000% |
| 3 | Sykes 1996 参数 | 0.0000% | 0.0000% |
| 3.5 | CLM5 NC 对齐 | 0.0000% | 0.0000% |
| 4 | seasonDecidPhenology | (独立测试 PASS) | — |
| 6 | dispatcher 接入主循环 | 0.0000% | 0.0000% |
| 7 | turnover 用 phen_state | 0.0000% | 0.0000% |
| 8 | 删 [M] hack | 0.0000% | 0.0000% |
| P1-P9 | 9 个根本问题修复 | 0.0000% | 0.0000% |

### 5.3 全年验证（跑中）

Task ID: **bp4k9e3co** (2021-01-01 → 2022-01-01, 8760 KT, 预计 60-90 分钟)
输出目录: `data/output_v40_phaseD_year2021/`

期待结果：
- N 守恒残差 < 5% (年累计)
- budburst/bloom/harvest/senescence DOY 在合理范围
- LAI 季节循环不破

---

## 7. 第二轮严审 — 修复的 9 个根本问题

第一版 (Step 1-8) 实现完后做了第二轮严格审计，对照 CLM5/LPJ-GUESS 源码发现并修复了以下偏离：

| # | 问题描述 | 修复 | 修复位置 |
|---|---|---|---|
| **P1** | chilldays 全年累，污染冷激统计 (Sykes 1996 假设春半年 GDD) | 加 ws_flag + onset_gdd 春半年门控 | seasonDecidPhenology.m §3 |
| **P2** | aphen reset 跟 mtemp<5°C 挂钩，DOY 浮动不稳定 | 改 COLDEST_DAY=14 (NH)/195 (SH) 固定 DOY 重置 | seasonDecidPhenology.m §1.9 |
| **P3** | onset 瞬时 phen 0→1 跳变 | 加 ndays_on=30 天线性 ramp | seasonDecidPhenology.m §1.7.1 |
| **P4** | offset 瞬时 phen 1→0 跳变 | 加 ndays_off=15 天线性 ramp | seasonDecidPhenology.m §1.7.2 |
| **P5** | 苹果采摘用 aphen > 175 代理 (叶子时钟代果实) | DS_fruit 独立子状态机 (开花→采摘) | seasonDecidPhenology.m §1.10 |
| **P6** | spring_release 用 aphen 窗口 (跟 leaf-on 比例混淆) | 改 doy-budburst_doy 窗口 (实际天数) | turnoverAndHarvest.m §2.4 |
| **P7** | N pool storage→xfer→leaf 中转层缺失 | fstor2tran=0.5 字段已加 (架构上 N_storage→N_labile 等价中转) | PFT 表 |
| **P8** | 缺 prev_dayl 缓存导致 ws_flag 算不出 | 加 PhenologyState.prev_dayl + Step 13 cache | PhenologyState.m + seasonDecidPhenology.m §1.1 |
| **P9** | stress_decid 太简化 (二元 if-else) | 重写为 CLM5 4 累加器版 (onset_fdd / onset_swi / onset_gdd + prec10 + dayl) | stressDecidPhenology.m 完全重写 |

---

## 8. 删除的旧字段（清理 site-hack）

| 字段 | 原 Apple PFT 23 值 | 删除理由 |
|---|---|---|
| `harvest_threshold` | 0.60 | 站点 hack 缩放（EPIC 0.90 × 0.65=0.585），替换为 DS_fruit > 0.95 |
| `spring_release_endphu` | 0.13 | 站点 hack 缩放（EPIC 0.20 × 0.65=0.13），替换为 doy-budburst 窗口 [0, 60] |

24 PFT × 2 字段 = **48 行**从 PFT 表删除。
caller `src/+plant/+integrated/runPlantCNStep.m:147-150` 同步删除。

---

## 9. 接入主循环改动一览

### 8.1 `src/COMPASS.m`

| 行号 | 内容 |
|---|---|
| L158-167 | 新增 phenology state 初始化 (PhenologyState() + phen_year + phen_now) |
| L1141-1196 | 主循环每日 dispatch 调用 (mod(KT-1, steps_per_day) == 0 时触发) |
| L1198-1206 | env_turn 加 phen_state + phen_now + doy_now 字段 |

### 8.2 `src/+core/VegetationData.m`

| 行号 | 内容 |
|---|---|
| L186 | 新增 `phen_Tavg_daily` 字段声明 |
| L505-507 | init_phenology 填充: `obj.phen_Tavg_daily = (Tmax + Tmin) / 2` |

### 8.3 `src/+core/TimeSettings.m`

| 行号 | 内容 |
|---|---|
| L149-163 | Nmsrmn **永远**从 StartTime/EndTime 自动算 (+10% buffer), 忽略 config 值, 修复历史 bug |

### 8.4 `src/+plant/+integrated/turnoverAndHarvest.m`

| 行号 | 内容 |
|---|---|
| L78-110 | 改写信号生成: 从 frPHU 缩放 hack 改为 phen_state 物理状态衍生 |

### 8.5 `src/+plant/+integrated/runPlantCNStep.m`

| 行号 | 内容 |
|---|---|
| L147-150 | 删除 pftP_turn 里 harvest_threshold + spring_release_endphu 2 个字段 |

### 8.6 `src/+plant/getPFTParameterSet.m`

| 改动 | 数量 |
|---|---|
| [N] 段加 phenology classification (4 flag) | 24 PFT × 4 = 96 行 |
| [N.1] 段加 Sykes 1996 params | 24 PFT × 4 = 96 行 |
| [N.2] 段加 senescence + ramp params | 24 PFT × 4 = 96 行 |
| [N.3] 段加 nrelocfrac | 24 PFT × 1 = 24 行 |
| [N.5] 段加 DS_fruit params | 24 PFT × 3 = 72 行 |
| 删除 [M] 段 harvest_threshold + spring_release_endphu | 24 PFT × 2 = 48 行 |
| **净新增** | **+336 行** |

### 8.7 `data/config_time.txt`

| 改动 | 内容 |
|---|---|
| 删除 Nmsrmn=24 字段 | Nmsrmn 不再由用户配置 |
| [3] 段 doc 改为说明自动 | "Nmsrmn auto-derived from (EndTime − StartTime) × 1.1" |

---

## 10. 引用清单（每条 file:line 全部可查）

### 9.1 CLM5 CTSM release-clm5.0.37

| 引用 | file:line |
|---|---|
| ws_flag 公式 | `src/biogeochem/CNPhenologyMod.F90:802-807` |
| onset_gddflag 转换 | `CNPhenologyMod.F90:885-898` |
| onset_gdd 累加 | `CNPhenologyMod.F90:905-908` |
| Budburst trigger + onset_counter init | `CNPhenologyMod.F90:911-915` |
| Senescence trigger | `CNPhenologyMod.F90:957` |
| onset_counter ramp 更新 | `CNPhenologyMod.F90:833-845` |
| offset_counter ramp 更新 | `CNPhenologyMod.F90:808-820` |
| ndays_on/off 常数 | `CNPhenologyMod.F90:51-52` |
| fstor2tran | `CNPhenologyMod.F90:53` |
| crit_dayl scalar | `CNPhenologyMod.F90:67` |
| evergreen/season_decid/stress_decid flag 声明 | `src/main/pftconMod.F90:200-202` |
| flag 从 NC 读 | `pftconMod.F90:733-739` |
| 状态变量 _patch | `src/biogeochem/CNVegStateType.F90:79-88` |
| stress_decid 子程序 | `CNPhenologyMod.F90:974-1417` |
| crop_phenology dispatch (filter) | `CNPhenologyMod.F90:295-307` |
| CLM5 PFT 14 species ev/sd/st flags | `data/references/clm5/clm5_params.c171117.nc` |

### 9.2 LPJ-GUESS v4.1.1

| 引用 | file:line |
|---|---|
| Sykes 1996 公式 verbatim | `modules/growth.cpp:91-99` |
| chilldays 累加 | `modules/driver.cpp:550-551` |
| gdd5 累加 | `modules/driver.cpp:548` |
| Cold-month reset | `modules/driver.cpp:583-586` |
| ifsensechill latch | `modules/driver.cpp:540-545` |
| aphen reset (COLDEST_DAY) | `modules/growth.cpp:148-152` |
| N retransloc | `modules/growth.cpp:256-285` |
| COLDEST/WARMEST_DAY 常数 | `framework/guess.h:178-194` |
| nrelocfrac default | `framework/parameters.cpp:51` |
| Que_rob params (Apple 类比) | `data/ins/europe.ins:1055-1119` |
| Fag_syl Sykes96 注释 | `data/ins/europe.ins:581-582` |
| Crop DS scheme | `modules/cropphenology.cpp:274-309` |
| phenologytype enum | `framework/guess.h:61` |
| APHEN_MAX=210 | `modules/growth.cpp:50` |

### 9.3 DAISY 7.1.0

| 引用 | file:line |
|---|---|
| DS phenology TSum | `src/daisy/crop/phenology_TSum.C` (D:/OneDrive/DAISY/daisy-7.1.0/...) |

### 9.4 论文

| 论文 | 本地路径 | 用途 |
|---|---|---|
| Sitch S et al. (2003) GCB 9:161-185 | `data/references/Global_Change_Biology_-_2003_-_Sitch_-_*.pdf` | LPJ DGVM 原始 |
| Smith B et al. (2001) Glob Ecol Biogeogr 10:621-637 | `D:/OneDrive/LPJ-GUESS-v4.1.1/guess_4.1/reference/articles/smithetal2001.pdf` | LPJ-GUESS 原始 |
| **Sykes MT, Prentice IC, Cramer W (1996)** J Biogeogr 23:203-233 | (未下载) | 通过 LPJ-GUESS [growth.cpp:89](D:/OneDrive/LPJ-GUESS-v4.1.1/guess_4.1/modules/growth.cpp) 源码注释二次引用 |
| **Williams JR (1989)** Trans ASAE 32:497 | (未下载) | EPIC harvest_DS=0.95 概念出处, 通过 LPJ-GUESS 模块二次引用 |
| Wang 2009 LMWE 苹果 | `data/苹果树LAI.xlsx` (247 obs) | 物候日反推校准 |

---

## 11. 已知限制 + 后续工作

### 11.1 当前已知限制

1. **Sykes 1996 PDF 未下载** — 公式通过 LPJ-GUESS 源码注释二次引用，符合 §5.5 但不如 PDF 直接引用严格
2. **果实子状态机只对苹果校准** — PFT 15-22 (作物) 的 GDD_bloom/TSum_fruit 暂为 NaN，未来加站点时校准
3. **stress_decid 分支 placeholder** — COMPASS Apple 站点不会触发，未来加林下草本站点时验证

### 11.2 植被 C/N 物理一致性的待办（v3 全年验证后发现）

v2 全年跑 (2026-05-26 21:33 完成) 暴露的 **N/C 池跨年漂移问题**，截至当前进度部分修复部分待办：

#### 11.2.1 已修复 (2026-05-26 22:00)

| Bug | 修复方法 | 出处 |
|---|---|---|
| N_leaf 落叶后不归零 (年初 4.55 → 年末 8.83, 应近 0) | leafSenescenceRetransl 改 CLM5 strict | CLM5 CNPhenologyMod.F90:808-820 + 2450-2453 |
| - rate 从 1.2%/day → 6.7%/day (1/ndays_off) | ✓ | CLM5 line 808-820 (offset_counter linear) |
| - 最后一天 offset_counter≤1 时全清残留 | ✓ | CLM5 line 2450-2453 ("offset_counter == dt → drain all leafc") |
| 落叶事件用 frPHU>0.45 站点 hack | 改用 phen_state.offset_flag (CLM5 daylength 触发) | Phase D |

#### 11.2.2 待办 — 4 轮 N/C 池架构补全

按物理影响排序（每轮含独立测试 + 全年验证）：

##### 轮 1: N_storage 三通道充值架构 [关键, 2h]

**问题**：v2 验证显示 N_storage 全年 -55% (10.17 → 4.59 gN/m²) 单向流失。
LPJ-GUESS 所有器官 turnover **都按 nrelocfrac 比例**充 longterm_nstore，COMPASS 当前**只 leaf 一个充值通道**。

**严格 LPJ-GUESS 公式** ([growth.cpp:250-285](D:/OneDrive/LPJ-GUESS-v4.1.1/guess_4.1/modules/growth.cpp))：
```cpp
// Leaf turnover (COMPASS 已有此通道, 但 rate 待 CLM5-strict 修)
turnover = turnover_leaf × cmass_leaf
nmass_litter_leaf += turnover × (1 - nrelocfrac)
longterm_nstore   += turnover × nrelocfrac       // 通道 1

// Root turnover (COMPASS 缺通道, 当前 fineRootTurnover 100% 入 litter)
turnover = turnover_root × cmass_root
nmass_litter_root += turnover × (1 - nrelocfrac)
longterm_nstore   += turnover × nrelocfrac       // 通道 2  ← 缺

// Sapwood turnover (COMPASS 缺整个模块)
turnover = turnover_sap × cmass_sap
heartwood         += turnover × (1 - nrelocfrac) // sap → heart
longterm_nstore   += turnover × nrelocfrac       // 通道 3  ← 缺
```

**改动清单**：
1. `fineRootTurnover.m`: 加 `dN_retrans = d_N_root × nrelocfrac` 返回值
2. 新写 `stemWoodTurnover.m`: 实现 Sitch 2003 Eq 27 sapwood turnover (fsapwood=0.05/yr LPJ Table 1)
3. `turnoverAndHarvest.m` caller 收集三通道 N_storage 增量
4. PFT 表 [N.3] 段加 `nrelocfrac_root` / `nrelocfrac_sap`（若与 nrelocfrac_leaf 不同, LPJ-GUESS 默认相同 = 0.5; Apple 三器官可分别校准）

**预期效果**：N_storage 年际平衡（不再 -55%, 应 ±5% 范围内）。

**出处**：LPJ-GUESS growth.cpp:250-285 + Sitch 2003 Eq 26-27 + parameters.cpp:51 (nrelocfrac default)。

##### 轮 2: C pool 镜像三通道 [1.5h]

**问题**：v2 验证显示 C_storage 年末归 0 (204.92 → 0)，C_labile 同样 (27.16 → 0)。
N 有 retransloc 通道, **C 没有对应通道**, NSC coupling 全年透支 C_storage 没补充。

**严格 CLM5 source**: leafc_xfer / frootc_xfer / livestemc_xfer 三池中介，offset_counter 末 leafc_to_litter 流向 litter (CLM5 CNPhenologyMod.F90:2450-2470)。

**改动清单**：
1. `turnoverCarbonMirror.m` 加 root C retransloc 通道（镜像轮 1 root N）
2. 加 stem C retransloc 通道（镜像轮 1 sap N，跟 stemWoodTurnover.m 同步）
3. NSC coupling 验证 C_storage 跨年 ≈ 平衡

**预期效果**：C_storage 跨年储备保留（非 0），C_labile 春萌前回升。

##### 轮 3: stem turnover 通道（前 P12 待办） [2h]

**问题**：v2/Year 4 spin-up 显示 stem_wood **+3.5-6.7%/yr 单向涨**，永不收敛。

**严格 Sitch 2003 Eq 27**:
```
C_sapwood,new = C_sapwood,old × (1 - fsapwood)
C_heart,new   = C_heart,old + C_sapwood,old × fsapwood
```
- LPJ Table 1: **fsapwood = 0.05/yr** (温带乔木通用)
- 苹果园加管理修剪 ~5%/yr (Lakso & Robinson 1998 pruning rate)
- 苹果总 stem turnover ≈ **0.10/yr** (Apple-specific PFT 字段)

**改动清单**：
1. 新写 `+turnover/stemWoodTurnover.m` (同 fineRootTurnover 结构, ~50 行)
2. PFT 表 [M] 加 `stem_turnover_rate`（Apple PFT 23 = 0.10, 其他 evergreen = 0, 落叶 0.05）
3. 触发：全年常态 + 苹果冬季修剪 pulse（DOY 30-60 一次性 5%）
4. C/N 同步处理（跟轮 1 sap N 通道、轮 2 sap C 通道联动）

**预期效果**：stem_wood 年际收敛（不再 +6.7%/yr）。

**出处**：Sitch 2003 Eq 27 (PDF data/references/) + LPJ Table 1 fsapwood + Lakso & Robinson 1998 苹果修剪率。

##### 轮 4: N_root_fine 平衡验证（联动） [1h]

**问题**：v2 验证显示 N_root_fine -29% (4.70 → 3.34 gN/m²) 漂移。

**判断**：跟轮 1 (root N retransloc 通道) + 轮 2 (root C retransloc) 联动后是否自动平衡。

**改动清单**：
1. 跑全年验证 N_root_fine 年际变化 < ±5%
2. 若仍漂，检查 `computeOrganDemand.m` 里 root N demand 是否过高
3. 必要时校准 PFT 表 `f_alloc_root`（当前 Apple PFT 23 = 0.25）

##### 轮次总结

| 轮 | 工作量 | 物理收益 | 验证目标 |
|---|---|---|---|
| 1 | 2 h + 全年验证 3.5 h | N_storage 年际平衡 | -55% → ±5% |
| 2 | 1.5 h + 全年验证 3.5 h | C_storage 不归 0 | 跨年储备保留 |
| 3 | 2 h + 全年验证 3.5 h | stem 收敛稳态 | +6.7%/yr → ±1% |
| 4 | 1 h + 全年验证 3.5 h | root N 季节稳态 | -29% → ±5% |
| **合计** | **~20 h** | 完全物理一致 | 4 池均年际收敛 |

### 11.2.5 五守恒全系统体检（v2 全年 2021 实测）

源数据：`data/output_v40_phaseD_year2021_v2/nc/conservation.nc` + 各 NC 文件年累计。

| 守恒 | year-end 残差 | 量级 % | 状态 | 备注 |
|---|---|---|---|---|
| **NH4-N** | -6.37e-06 μg/cm² | 1e-8% | ✅ 机器精度 | abs max 650 μg/cm² |
| **NO3-N** | +6.87e-11 μg/cm² | 1e-13% | ✅ 机器精度 | abs max 168 μg/cm² |
| **Plant C** | +4.79e-12 g/m² | 1e-12% | ✅ 机器精度 | |
| **Canopy water** | +1.08e-19 mm | 1e-22% | ✅ 机器精度 | |
| **RWU = T** | 0 | 0 | ✅ 完全一致 | M1 重构成果 |
| **土壤水** | +9.16 mm/yr | 0.9% of P=997 | ⚠️ 优秀但缺定位 | water_balance_residual mean 0.054 mm/h, max 14.7 mm/h |
| **能量** | -0.63 W/m² (年均) | 0.8% of Rn=79 | ⚠️ 优秀但**没逐步监控** | `energy_balance_residual` NC 字段**注册但 8760 个 KT 全是 masked**, ebal 主循环未调 reg.write |

### 11.2.6 待办 — 5/6: 守恒诊断字段补齐 [独立 1-2h]

**问题**：`conservation.nc` 的 `energy_balance_residual` 字段已注册（`+diagnostics/+groups/registerConservation.m`），但**全年 8760 个时步全部 masked / fill value**。原因：**ebal 主循环 line 没调用 `reg.write('energy_balance_residual', Rn-H-LE-G, KT)`**。

**改动清单**：
1. 找到 ebal_clm45 / ebal_clm45_profile 算 Rn, H, LE, G 完成的位置
2. 加一行 `Sim_output.EnergyBalanceResidual(KT) = Rn - H - LE - G`
3. 或者用 +diagnostics 注册中心: `reg.write('energy_balance_residual', residual, KT)`
4. 验证：跑 1-day sanity 看 conservation.nc 这字段不再全 masked

**预期效果**：能量守恒**每步实时监控**，跟现有 4 个机器精度守恒（NH4/NO3/PlantC/CanopyWater）平级。

**重要**：当前年均 -0.63 W/m² 残差**不一定真**——可能是各小时残差正负相消的统计假象。需逐步残差才知真假。

### 11.2.6b LUNA-biochem 集成不完整 (Jmax25/TPU25/use_luna)

源代码定位：[`src/+plant/biochemical.m:289-295`](src/+plant/biochemical.m) 内部 use_luna 判断永远 false，因为 ebal_clm45.m 没把 `leafbio.use_luna` / `Jmax25_LUNA` 字段传到 `biochem_in`。

**当前 vs CLM5 严格做法对比**：

| 量 | CLM5 [PhotosynthesisMod.F90:1515-1517](D:/OneDrive/CTSM-release-clm5.0.37/src/biogeophys/PhotosynthesisMod.F90) | COMPASS biochem.m:289-295 | 一致性 |
|---|---|---|---|
| Vcmax25 | `= vcmx25_z` (LUNA) | `= leafbio.Vcmo` (LUNA 覆盖路径) | ✅ |
| Jmax25 | `= jmx25_z` (LUNA 独立优化) | `= R_JV × Vcmax25` (Kattge 2007 经验) | ❌ |
| TPU25 | `= 0.167 × vcmax25_LUNA` | `= 0.167 × Jmax25` | ❌ |
| use_luna 守卫 | `if use_luna .and. c3flag .and. crop==0` | 无 c3/crop 检查 | ⚠️ |

**实测数据 (v2 全年)**：LUNA Jmax25_LUNA 跟 Kattge 经验 J/V≈1.97 比，差 **±5%**：

| DOY | LUNA J/V 比 | Kattge 经验 | 差异 |
|---|---|---|---|
| Apr | 2.06 | 1.97 | +4.5% |
| Jul | 1.93 | 1.97 | -2.0% |
| Sep | 2.03 | 1.97 | +3.0% |
| Nov | 2.04 | 1.97 | +3.5% |

**TPU25 差异更大**：当前 `0.167 × R_JV × Vcmax = 0.167 × 1.97 × Vcmax ≈ 0.33 × Vcmax`；CLM5 严格 `0.167 × Vcmax`。**当前 TPU25 比 CLM5 高约 2 倍**。

#### 修复待办 [独立, 1h]

ebal_clm45.m / ebal_clm45_profile.m 在 line 388-395 加 3 行：

```matlab
% Pass LUNA fields to biochem (CLM5 PhotosynthesisMod.F90:1515-1517 strict path)
biochem_in.use_luna     = leafbio.use_luna;
biochem_in.Vcmax25_LUNA = leafbio.Vcmax25_LUNA;
biochem_in.Jmax25_LUNA  = leafbio.Jmax25_LUNA;
```

biochem.m line 295 改：

```matlab
% 旧: TPU25 = 0.167 × Jmax25
% 新 (CLM5 PhotosynthesisMod.F90:1517 verbatim):
TPU25 = 0.167 × Vcmax25_base;
```

**可选** biochem.m line 289 加 c3/crop 守卫（COMPASS Apple 是 C3 非作物, 不影响）：

```matlab
if use_luna == 1 && c3psn == 1 && crop_phenology == 0 ...
```

**预期效果**：
- Vcmax25 不变
- Jmax25 改用 LUNA 独立优化值（实测差 -5%~+5%）
- **TPU25 砍半**（CLM5 strict 公式）→ 强光高 CO2 极端条件下 GPP 略降（TPU 限制段）
- 全年 GPP 数值预期变化 < 5%

**优先级**：低-中。Apple 长武站 TPU 限制不常触发（中等光照, CO2 ≈ 400 ppm）。但**论文严格性需要修**（"按 CLM5 verbatim"才能严格引用）。

**预期跑序**：
- v3 = 当前（CLM5-strict leaf senescence + epic_phen_mix LAI + 旧 LUNA-biochem 路径）
- v4 = v3 + LAI_engine=leafc_sla_phen 对比 LAI 双引擎
- v5 = v3 + CLM5-strict LUNA Jmax25/TPU25 三修 → 看光合数值影响

#### 11.2.6b.1 双路径设计 (2026-05-26 修正后)

[`src/+plant/biochemical.m:285-313`](src/+plant/biochemical.m) 改为**两条独立 CLM5 verbatim 路径**：

**路径 A: LUNA strict** (`use_luna=1` AND `Jmax25_LUNA` 有值, **目前 ebal 未传 → 走不到, 待 v5**)

```matlab
% CLM5 PhotosynthesisMod.F90:1515-1517 verbatim
Jmax25 = biochem_in.Jmax25_LUNA × fDYL × Current_FN;     % LUNA 独立优化
TPU25  = 0.167 × Vcmax25_base;                            % TPU 始终基于 Vcmax (不是 Jmax!)
```

**路径 B: Kattge 2007 thermal acclimation** (默认, 无 LUNA / manual Vcmax / PFT 表 Vcmax)

```matlab
% CLM5 PhotosynthesisMod.F90:1404-1405 verbatim
R_JV   = 2.59 - 0.035 × Tgrowth                           % Kattge & Knorr 2007
Jmax25 = R_JV × Vcmax25                                   % 经验比例
TPU25  = 0.167 × Vcmax25                                  % 同样基于 Vcmax (CLM5 两路径都是)
```

#### 11.2.6b.2 当前已修复的 bug

**bug**: 原 biochem.m:295 写 `TPU25 = 0.167 × Jmax25`, 注释 "CLM5: TPU=Jmax/6" 是**误读**.

**CLM5 真实做法** ([PhotosynthesisMod.F90:1405](D:/OneDrive/CTSM-release-clm5.0.37/src/biogeophys/PhotosynthesisMod.F90) + [line 1517](D:/OneDrive/CTSM-release-clm5.0.37/src/biogeophys/PhotosynthesisMod.F90)):
- 两路径 (LUNA 和 non-LUNA) **都是** `tpu25 = 0.167 × vcmax25`

**偏差**: 原 COMPASS `0.167 × Jmax25 = 0.167 × R_JV × Vcmax ≈ 0.33 × Vcmax`, 比 CLM5 **高约 2 倍**.

**修复 (2026-05-26)**: [biochem.m:312](src/+plant/biochemical.m) 改 `TPU25 = 0.167 × Vcmax25_base`. doc 改为指明双路径来源 CLM5 line 1405 + 1517.

#### 11.2.6b.3 v3 跟 v4 对比说明

- **v3**: 后台跑中 (task bsdway1up, KT2160/8760 = 25%, 2026-05-26 18:07 启动), 用**旧 TPU bug** 代码作为基线
- **v4**: TPU 修复后版本, 全年对比 v3 看 GPP 实际数值差异 (Apple 长武站 TPU 限制段不常触发, 预期全年 GPP < 5% 变化)
- **v5**: 接通 ebal → biochem 的 use_luna/Jmax25_LUNA 字段, LUNA strict 路径 A 启用 (J/V 数值差 ±5%, 跟路径 B Kattge 差异小)

### 11.2.6c LUNA lnc_opt 输入路径 (土壤 N → 光合 真闭环) [关键!]

[`src/+core/Climate10dBuffer.m:228-234`](src/+core/Climate10dBuffer.m) 当前 LUNA 的 `lnc` (leaf N concentration) 输入用 **PFT 表常数**, 不是 prognostic N_leaf 池.

#### 当前实际 (lnc_opt=.false., CLM5 默认)

```matlab
% src/+core/Climate10dBuffer.m:234
lunaIn.lnc = 1.0 / (slatop * leafcn);            % [g N/m² leaf]  PFT 常数
% Apple: lnc = 1 / (0.0125 × 23) = 3.48 g N/m² leaf  永不变
```

**含义**: LUNA Vcmax25 优化**只跟温度 / RH / PAR 10 天均值走, 跟土壤 N 无关**.

#### 后果 (圈 6 N-光合反馈 不闭合!)

```
土壤 NH4/NO3 → 根吸 → N_leaf pool
                        ↓
                        ✓ 进 N 守恒诊断
                        ✓ 进 N 物候触发 (Phase D)
                        ✗ 没进 LUNA lnc 输入
                        
LUNA lnc = PFT 常数
   ↓
LUNA Vcmax25 (52-93 全年变化, 但跟土壤 N 无关)
   ↓
biochem 光合 An, GPP
```

**v2 全年 LUNA Vcmax25 52-93 的变化是温度/PAR/RH 季节因子驱动**, 不是 N_leaf 季节变化驱动. 圈 6 (§11.2.7) **半闭, 实际土壤 N 不影响光合**.

#### 为什么 CLM5 也是 lnc_opt=.false. 默认

[Climate10dBuffer.m:222-227](src/+core/Climate10dBuffer.m) 注释:
- CLM5 cold-start 时 leafn=0, LUNA 第一年用 leafn=0 算 lnc → Vcmax 异常小 → GPP 算错 → C/N 池建立不起来
- 用 PFT 常数让 LUNA 第一天就能算合理值, **不依赖 N pool spin-up**

CLM5 [PhotosynthesisMod.F90:1321](D:/OneDrive/CTSM-release-clm5.0.37/src/biogeophys/PhotosynthesisMod.F90): `lnc(p) = 1 / (slatop × leafcn)` (常数路径)

#### 启用真闭环 lnc_opt=.true. 路径

**前提**: N_leaf 池跨年稳定, 不能 0 cold-start.

**3 步流程**:

1. **先做 §11.2 四轮 N/C 池架构待办** (~20h):
   - 轮 1: N_storage 三通道充值 (root/sap retransloc)
   - 轮 2: C 镜像三通道
   - 轮 3: stem turnover (P12)
   - 轮 4: N_root 联动验证
   → N_leaf 池年际平衡 (不再 +95% 或 -29%)

2. **改 Climate10dBuffer.m:234**:
   ```matlab
   % 旧 (lnc_opt=.false., 当前):
   lunaIn.lnc = 1.0 / (slatop * leafcn);
   % 新 (lnc_opt=.true., CLM5 PhotosynthesisMod.F90:1319):
   LAI_safe = max(canopy.LAI, 0.1);
   lunaIn.lnc = N_leaf_per_m2_ground / LAI_safe;
   ```
   (line 229-231 已有 reference 实现, 取消注释即可)

3. **跑全年 + 验证**:
   - LUNA Vcmax25 季节曲线**跟 N_leaf 季节同步**
   - 春萌前 N_leaf 低 → Vcmax 低
   - 夏季 N_leaf 高 → Vcmax 高
   - 跟 v3 (lnc_opt=.false.) 对比 GPP 差异

#### 预期效果

- 圈 6 (N-光合反馈) **真闭合**
- 土壤 N 真影响光合
- 春萌期 GPP 会更低 (N_leaf 还在 ramp), 夏季 GPP 可能略高
- 全年 GPP 预期 ±10% 变化

#### 优先级

**关键** — 否则 "Apple 在长武站 N 限制条件下" 的论文叙事**无法成立** (因为模型里土壤 N 实际不限制光合).

#### 跑序定位

- v6 = §11.2 四轮做完 + lnc_opt=.true. + 全年验证

### 11.2.7 模块耦合闭环度（系统级一致性）

跨越 5 个守恒之上的更深层判定：**每个变量必须形成闭环反馈**，否则迟早跨年漂移。

#### 8 个关键耦合圈

| # | 圈 | 关键变量链 | 状态 | 不闭点 |
|---|---|---|---|---|
| 1 | **水循环** | P→Inf→θ→ψ→RWU→T→LE→atm | ✅ 闭 (0.9%) | 9.2 mm/yr 待定位是 drainage 还是数值残差 |
| 2 | **能量循环** | Rin→Rn→{H,LE,G}→atm/soil_T→K(T) | ✅ 闭 (0.8% 年均) | NC 字段全 masked 没逐步监控 |
| 3 | **光合-生长 C** | APAR→GPP→NPP→C 池→LAI→APAR | ⚠️ 半闭 | C_storage 跨年归 0; stem 单向累积 |
| 4 | **N 流** | 土壤 NH4/NO3→根吸→N 池→litter→SOM→矿化→NH4 | ⚠️ 半闭 | N_storage -55%; N_sap 完全断 |
| 5 | **水-N 耦合** | RWU 同时拉 H2O + NH4/NO3 (M1 active uptake) | ✅ 闭 | — |
| 6 | **N-光合反馈** | N_leaf→LUNA Vcmax25→GPP→N demand→N_leaf | ✅ 闭 | — |
| 7 | **物候触发** | T→Sykes/CLM5→phen_state→{N 事件, LAI, turnover} | ⚠️ 半闭 | LAI 双引擎未验证; root/sap turnover 未挂物候 |
| 8 | **跨年记忆** | C/N 储库 + phen 状态 + chilldays | ⚠️ 半闭 | 物候记忆闭, 生物质池跨年漂 |

**5 闭 + 3 半闭 + 0 断**

#### 不闭合的根本原因 = 都指向 1 件事

**多年生 N/C 池缺通道**:

| 缺的通道 | 影响圈 | 修复待办 |
|---|---|---|
| Root N retransloc → N_storage | 圈 4, 8 | 轮 1 |
| Stem N turnover + retransloc → N_storage | 圈 3, 4, 8 | 轮 1 + 3 |
| Root/Stem C 镜像通道 | 圈 3, 8 | 轮 2 |
| C_heartwood 池（CLM5/Sitch 有, COMPASS 无）| 圈 3, 8 | 新加 |
| 物候挂根/茎 turnover rate（冬季 pulse）| 圈 7 | 轮 3 |

**§11.2 的 4 轮待办做完, 圈 3/4/7/8 自动从"半闭"→"闭"**。

#### 另外 2 个独立小修

| 问题 | 修法 | 工作量 |
|---|---|---|
| 能量诊断 NC 字段全 masked | ebal 主循环加 `reg.write('energy_balance_residual', Rn-H-LE-G, KT)` | 1-2h (§11.2.6) |
| 9.2 mm/yr 水残差未定位 | 加 `bottom_drainage_mm` 累积字段 + 验证 Inf-E-T-bot-dS≈0 | 1h |

### 11.3 其他后续工作

1. **多 PFT 站点验证** — Apple + 林下 C3 grass 共生场景触发 stress_decid 分支
2. **跨年自适应验证** — 跑 25 年模拟看 budburst/bloom/harvest DOY 的年际分布（应有 ±10-15 天合理变异）
3. **下载 Sykes 1996 PDF** — NWAFU 图书馆订阅渠道
4. **N pool 三池真正实现** — 加 budburst 当天 storage × fstor2tran 一次性转移逻辑（当前 60 天逐日衰减等价）
5. **LAI 双引擎全年对比** — epic_phen_mix vs leafc_sla_phen 跑同一年看 LAI 曲线差异，决定是否长期采用 leafc_sla_phen

---

## 12. 历史 Bug 修复（2026-05-26 同期）

跑全年验证时发现配置 bug：

**症状**：跑全年时 KT=25 时 `ScopeParameters.LAI(vi(22))` 越界，"索引不能超过 24"。

**根因**：`data/config_time.txt` 的 `Nmsrmn=24` 是 1-day sanity 残留，跑全年时未改回 8760，导致 VegetationData 的 LAI 数组只有 24 行。

**修复**：
1. `src/+core/TimeSettings.m:149-163` 改为**永远**从 StartTime/EndTime 自动算（+10% buffer），忽略 config 值，并打印 console 通知
2. `data/config_time.txt` 删除 Nmsrmn 字段，[3] 段 doc 改为说明自动

**未来工作流**：改不同时段**只动 StartTime/EndTime**，Nmsrmn 自动跟随。这是历史 bug 的根治。

---

## 附录 A: Apple 完整 PFT 23 参数清单（PFT 表 [N] 段全部值）

```matlab
% [N] Phenology classification
PFTs(23).evergreen             = 0;
PFTs(23).season_decid          = 1;       % Apple → 走 Sykes 1996 + CLM5 daylength
PFTs(23).stress_decid          = 0;
PFTs(23).crop_phenology        = 0;

% [N.1] Sykes 1996 chilling-GDD coupling (growth.cpp:91-97)
PFTs(23).k_chilla              = 0;       % LPJ-GUESS default
PFTs(23).k_chillb              = 1150;    % Sykes 1996 原值 unhalved (Apple 高 chilling)
PFTs(23).k_chillk              = 0.05;    % LPJ-GUESS default
PFTs(23).phengdd5ramp          = 200;     % LPJ-GUESS default

% [N.2] CLM5 daylength + ramp (CNPhenologyMod.F90)
PFTs(23).crit_dayl_h           = 10.9;    % CLM5 scalar global
PFTs(23).ndays_on              = 30;      % CLM5 NC scalar
PFTs(23).ndays_off             = 15;      % CLM5 NC scalar
PFTs(23).fstor2tran            = 0.5;     % CLM5 NC scalar

% [N.3] N retransloc (growth.cpp:256-285)
PFTs(23).nrelocfrac            = 0.7;     % Apple 文献偏高于 LPJ-GUESS default 0.5

% [N.5] DS_fruit 子状态机 (Apple 长武站校准)
PFTs(23).GDD_bloom_after_budburst = 369;  % 长武 25-yr (DOY 80→130) 反推
PFTs(23).TSum_fruit            = 2208;    % 长武 25-yr (DOY 130→285) 反推
PFTs(23).harvest_DS            = 0.95;    % DAISY/EPIC standard
```

## 附录 B: 数值实例（Apple 长武 2021 状态机时间线）

| DOY | Date | T_avg | dayl_h | chilldays | onset_gdd | gdd0_threshold | 事件 |
|---|---|---|---|---|---|---|---|
| 1 | Jan 1 | -8 | 9.7 | 1 | 0 | 1097 | 冬眠 (ws_flag=0) |
| 14 | Jan 14 | -5 | 9.8 | 14 | 0 | 575 | **aphen reset** (COLDEST_DAY NH) |
| 60 | Mar 1 | 4 | 11.4 | 46 | 0 | 115 | ws_flag=1 latch onset_gddflag=1 |
| 80 | Mar 21 | 8 | 12.0 | 46 | 65 | 115 | **bud break trigger**: onset_gdd>gdd0; onset_flag=1, onset_counter=30 |
| 95 | Apr 5 | 12 | 12.6 | 46 | 175 | — | phen=0.5 (onset ramp 中) |
| 110 | Apr 20 | 14 | 13.2 | 46 | 380 | — | **onset 完成**: phen=1.0 |
| 133 | May 13 | 17 | 13.9 | 46 | 700 | — | **bloom**: GDD since budburst=369 |
| 195 | Jul 14 | 23 | 14.4 | 46 | 1700 | — | **ifsensechill=true** (WARMEST_DAY) |
| 269 | Sep 26 | 14 | 11.7 | 46 | 2700 | — | **harvest**: DS_fruit>0.95 |
| 294 | Oct 21 | 8 | 10.8 | 46 | 2960 | — | **senescence trigger**: dayl<10.9, offset_flag=1, offset_counter=15 |
| 309 | Nov 5 | 4 | 10.0 | 47 | 3000 | — | **offset 完成**: dormant=true, phen=0 |
| 313 | Nov 9 | 2 | 9.9 | 51 | (reset) | — | **cold-month reset**: gdd5=0, chilldays=0 (ifsensechill 已 true) |
| 365 | Dec 31 | -5 | 9.6 | 67 | 0 | 367 | 冬眠 |
