# CLM5/LPJ PFT 生命周期补全操作单

目的：把 COMPASS 已有植被 C/N 循环中“多 PFT 生命周期出口没抄全”的部分，按 CLM5/LPJ-GUESS/DSSAT-CSM 源码补齐到可审状态。本文是给 Claude 执行的施工单；不要求自由发明新机理。

本文只处理 vegetation C/N lifecycle，不处理 SOM。SOM/litter 分解侧由其他工作线处理。

## 0. 总原则

1. 只按本地 CLM5 / LPJ-GUESS / DSSAT-CSM / DAISY 源码补缺，不为了追某个 CUE、NEE 或产量指标硬调参数。
2. 优先修“出口机制”：leaf/froot 背景凋落、作物 harvest/cleanup 的完整移植边界、通用 gap mortality 参数覆盖。
3. 苹果 PFT23 主线不能被破坏。任何通用补丁必须保持 PFT23 现有结果可解释。
4. 不改 storage/xfer 池结构。COMPASS 现在的 dual reserve / NSC / apple-specific pool design 是项目内决策，不向 CLM5 的 leaf/froot/live/dead/gresp 全套 storage/xfer 回退。
5. 不碰 SOM 内部过程。vegetation 出口只需要把 C/N 质量从植被池正确转出，并进入现有 litter/mortality/diagnostic 汇总。
6. 所有新增公式必须在代码注释里给出本地源码依据：文件名、行号或子程序名。

## 1. 依据源代码

本地参考源码路径：

- CLM5: `D:/OneDrive/应用/Model/CTSM-release-clm5.0.37`
- LPJ-GUESS: `D:/OneDrive/应用/Model/LPJ-GUESS-v4.1.1`
- DSSAT-CSM: `D:/OneDrive/应用/Model/CSM_Source_22Dec2015`

关键 CLM5 源码证据：

- Evergreen background litterfall:
  - `src/biogeochem/CNPhenologyMod.F90:603-604`
  - `bglfr = 1 / (leaf_long * dayspyr * secspday)`
- Stress-deciduous / long growing season background litterfall:
  - `src/biogeochem/CNPhenologyMod.F90:1370-1373`
  - `bglfr = (1 / (leaf_long * dayspyr * secspday)) * lgsf`
- Background litterfall application:
  - `src/biogeochem/CNPhenologyMod.F90:2623-2627`
  - if `bglfr > 0`, both `leafc_to_litter = bglfr * leafc` and `frootc_to_litter = bglfr * frootc`
- Crop phenology / harvest:
  - `src/biogeochem/CNPhenologyMod.F90:1420-2031`
  - harvest date set around `CNPhenologyMod.F90:1962-1964`
  - grain fill background litter around `CNPhenologyMod.F90:1994-1997`
- Crop state cleanup:
  - `src/biogeochem/CNCStateUpdate1Mod.F90:467-483`
  - after harvest, CLM5 drains small residual crop pools / xsmrpool / cpool / frootc to avoid carryover artifacts.
- Gap mortality:
  - `src/biogeochem/CNGapMortalityMod.F90:183-199`
  - mortality fluxes apply to leaf, froot, livestem, livecroot, deadstem, deadcroot.
  - wood mortality goes to CWD aggregation at `CNGapMortalityMod.F90:409-413`.

关键 LPJ-GUESS 源码证据：

- Crop development stage:
  - `guess_4.1/modules/cropphenology.cpp:260-309`
  - development rate uses temperature response, vernalization and photoperiod gates.
- Crop phenology / harvest control:
  - `guess_4.1/modules/cropphenology.cpp:311-315` onward.

关键 DSSAT-CSM 源码证据：

- Main crop driver:
  - `CROPGRO.for:243-250` calls `PHENOL` and receives `MDATE`, stages and vegetative/reproductive clocks.
  - `CROPGRO.for:265-279` calls `DEMAND`.
  - `CROPGRO.for:290-293` calls `NUPTAK`.
  - `CROPGRO.for:329-341` calls `VEGGR`.
  - `CROPGRO.for:346-358` calls `SENES` and `ROOTS`.
  - `CROPGRO.for:362-389` calls `GROW` to update organ biomass, N, senescence and LAI.
  - `CROPGRO.for:402-410` and `638-646` call `OPHARV`.
- Harvest timing and harvest fractions:
  - `AUTHAR.for:18-21` defines automatic harvest inputs/outputs.
  - `AUTHAR.for:119-120` maps harvest product/byproduct fractions.
  - `AUTHAR.for:132-183` supports maturity-date, fixed-date, days-after-planting, growth-stage and automatic-window harvest.
- Seasonal crop outputs and actual harvested yield:
  - `OPHARV.for:24-32` harvest output interface.
  - `OPHARV.for:508-522` computes harvested seed and byproduct using `HARVFRAC`.
- Crop-specific modules:
  - `PHENOL.for:21-28` generic CROPGRO phenology interface.
  - `CSP_PHENOL.for:20-30` sugarcane-specific phenology interface.
  - `CSCERES_Interface.for:123-152` cereal interface returns residue/senescence to CSM residue objects.
  - `CSCRP_Interface.for:192-237` wheat/barley interface returns residue/senescence and harvest residue.

关键 DAISY 源码证据：

- DAISY path: `D:/OneDrive/应用/Model/DAISY/daisy-7.1.0`
- Crop interface:
  - `include/daisy/crop/crop.h:92-100` crop `tick` couples time, bioclimate, soil, soil heat, soil water, chemistry and organic matter.
  - `include/daisy/crop/crop.h:102-116` harvest returns yield and residual C/N to topsoil/soil layers.
  - `include/daisy/crop/crop.h:143-151` exposes DS, shoot DM, storage-organ DM, total N/C and effective root density.
- Standard crop composition:
  - `src/daisy/crop/crop_std.C:61-80` composes seed, root system, canopy, harvesting, production, phenology, partition, vernalization, photosynthesis and crop N.
  - `src/daisy/crop/crop_std.C:291-307` updates dead organic matter, canopy, root density and nitrogen content after initialization.
  - `src/daisy/crop/crop_std.C:720-733` calls harvest with root density, phenological dates, production pools and stress counters.
- DAISY simple phenology / harvest:
  - `src/daisy/crop/phenology_TSum.C:34-49` defines thermal-sum crop phenology with emergence/vegetative/reproductive phases.
  - `src/daisy/crop/phenology_TSum.C:70-89` advances DS and marks flowering/ripe.
  - `src/daisy/crop/harvesting.C:49-76` harvest interface.
  - `src/daisy/crop/harvesting.C:80-203` removes stem/leaf/storage-organ C/N/DM and separates yield from field losses.

## 2. 当前 COMPASS 落点

优先阅读这些文件，不要先全库乱改：

- `src/+plant/+integrated/runPlantCNStep.m`
- `src/+plant/+integrated/turnoverAndHarvest.m`
- `src/+plant/+integrated/turnoverCarbonMirror.m`
- `src/+plant/+integrated/+phenology/evergreenPhenology.m`
- `src/+plant/+integrated/+phenology/cropPhenology.m`
- `src/+plant/+integrated/+phenology/phenologyDispatch.m`
- `src/+plant/getPFTParameterSet.m`
- `src/+plant/getPFTValue.m`
- `src/COMPASS.m`
- relevant diagnostics under `src/+diagnostics/+groups/registerPlantNitrogen.m`

已知现状：

- PFT23 apple 已有 fruit tree allocation、harvest、spring reserve、senescence、pruning、background mortality，是当前主线。
- PFT1-14/24 有部分 CLM5 allometry 参数，但生命周期出口不完整。
- PFT15-22 标了 `crop_phenology=1` 和 `has_fruit=1`，但 `cropPhenology.m` 明确是 placeholder，`turnoverAndHarvest.m` 对 annual PFT 早退，crop harvest/cleanup 没真正接上。
- `evergreenPhenology.m` 注释已经提到 CLM5 `bglfr`，但实际只设置 `phen=1`。
- `runPlantCNStep.m` 已有 Step 10 background gap mortality，但默认 `r_mort=0`，目前主要对 PFT23 生效。

## 3. 不改项

### 3.1 不改 storage/xfer 池结构

不要补 CLM5 的全套：

- `leafc_storage/xfer`
- `frootc_storage/xfer`
- `livestemc_storage/xfer`
- `deadstemc_storage/xfer`
- `livecrootc_storage/xfer`
- `deadcrootc_storage/xfer`
- `gresp_storage/xfer`

理由：COMPASS 当前池子是面向苹果 PFT 和 NSC/reserve 机制重构过的，项目判断是当前池结构更适合本模型。本文只修出口闭合，不把池结构改回 CLM5。

### 3.2 不改 SOM

不要新建或重写 litter/SOM 分解模块。vegetation 侧只需要把 `leaf_litter_today`、`root_litter_today`、`cum_C_mortality`、`cum_N_mortality` 等出口正确记账。

### 3.3 不为结果硬调

不要通过随意改 `r_mort`、`leaf_long_days`、`root_turnover_rate` 来追 CUE/NEE/LAI。参数改动必须来自 CLM5/LPJ/FruitTree 参数表或项目已有 calibration 文档。

## 4. Patch A：补 CLM5 bglfr 背景 leaf/froot turnover

### 4.1 目标

让 evergreen / grass / stress-deciduous 等非苹果或非 storage PFT 不再只“停分配”，而是按 CLM5 的 `bglfr` 把 leaf 和 fine root 以 leaf longevity 控制的速率转出。

CLM5 逻辑：

```text
bglfr = 1 / (leaf_long * dayspyr * secspday)
leaf_to_litter  = bglfr * leaf
froot_to_litter = bglfr * froot
```

注意单位转换：

- CLM5 `leaf_long` 单位是 year。
- COMPASS `leaf_long_days` 单位是 day。
- COMPASS 每步转出分数可写成：

```matlab
turnover_frac = min(1, dt / (leaf_long_days * 86400));
```

如果采用 stress-deciduous long-growing-season factor：

```matlab
turnover_frac = turnover_frac * lgsf_or_phen_gate;
```

但第一版建议先不要发明复杂 `lgsf`。COMPASS 没有 CLM5 原生 `lgsf` 状态时，可以用已存在 `phen_state.phen` 作为保守 gate，并在注释中明确这是 COMPASS 映射，不是 CLM5 原变量。

### 4.2 建议实现位置

优先在 `turnoverAndHarvest.m` 里实现，不要把实际落叶逻辑放在 `evergreenPhenology.m`。原因：

- phenology 函数只产 `phen/onset/offset` 状态。
- 真实 pool 扣减和 litter 诊断都集中在 turnover 层。
- C side 通过 `turnoverCarbonMirror.m` 镜像 N 诊断，减少双写风险。

### 4.3 具体步骤

1. 在 `turnoverAndHarvest.m` 初始化 `diag` 时新增诊断字段：
   - `leaf_background_litter_today`
   - `root_background_litter_today`
   - `bglfr_turnover_frac`

2. 在 annual 早退之前判断作物。不要让 crop placeholder 走自然 PFT 背景 turnover。

3. 在 leaf senescence block 之后、fine-root turnover 之前，新增一个 background litterfall block：

```matlab
is_crop = isfield(pftParams, 'crop_phenology') && pftParams.crop_phenology == 1;
is_storage_deciduous = isfield(pftParams, 'has_storage_pool') && pftParams.has_storage_pool == 1;

use_bglfr = ~is_crop;

% 对 apple/storage-deciduous，已有 autumn offset flush；背景 bglfr 可以先只用于 evergreen/grass。
% 如果要对 storage-deciduous 也开，必须证明不会和 offset leafSenescenceRetransl 双算。
use_bglfr = use_bglfr && ~is_storage_deciduous;
```

4. 计算 leaf/froot 同步背景转出：

```matlab
if use_bglfr && isfield(pftParams, 'leaf_long_days') && pftParams.leaf_long_days > 0
    bglfr_frac = min(1, dt / (pftParams.leaf_long_days * 86400));
    if isfield(env, 'phen_state') && isfield(env.phen_state, 'phen')
        phen_gate = max(0, min(1, env.phen_state.phen));
    else
        phen_gate = 1.0;
    end
    bglfr_frac = bglfr_frac * phen_gate;

    dN_leaf_bg = min(pools_out.N_leaf, pools_out.N_leaf * bglfr_frac);
    pools_out.N_leaf = pools_out.N_leaf - dN_leaf_bg;
    diag.leaf_litter_today = diag.leaf_litter_today + dN_leaf_bg;
    diag.leaf_background_litter_today = dN_leaf_bg;

    dN_root_bg = min(pools_out.N_root_fine, pools_out.N_root_fine * bglfr_frac);
    pools_out.N_root_fine = pools_out.N_root_fine - dN_root_bg;
    diag.root_litter_today = diag.root_litter_today + dN_root_bg;
    diag.root_background_litter_today = dN_root_bg;
    diag.bglfr_turnover_frac = bglfr_frac;
end
```

5. 检查 `turnoverCarbonMirror.m` 是否已经用 `diag.leaf_litter_today` 和 `diag.root_litter_today` 镜像 C。如果是，则不要额外写 C 扣减，避免双扣。

6. 更新 `evergreenPhenology.m` 注释：说明它只设置 `phen=1`，实际 `bglfr` 出口在 `turnoverAndHarvest.m`。

### 4.4 禁止事项

- 不要把 `leaf_long_days` 改成 365 来“看着稳定”。参数表已有按 CLM5 audit 更新的 PFT 值。
- 不要让 PFT23 apple 默认走背景 bglfr，除非单独设计和秋季 offset flush 的防双算逻辑。
- 不要在 `evergreenPhenology.m` 直接扣 pool。

### 4.5 验收

静态验收：

- `leaf_long_days` 至少在一个真实 turnover 公式中使用。
- `evergreenPhenology.m` 不再只是注释声称有 bglfr，而实际无人消费。
- `turnoverAndHarvest.m` 中 leaf 和 froot 使用同一个 `bglfr_frac`。

行为验收：

- PFT1 evergreen 单步运行时，`phen=1` 且 `leaf_background_litter_today > 0`。
- PFT13 grass 在生长季 `phen>0` 时有 leaf/froot 背景 turnover。
- PFT23 apple 的 harvest、senescence、reserve 指标不因 Patch A 明显改变。

## 5. Patch B：作物主干选择与完整移植路线

### 5.1 结论

COMPASS 收尾必须有作物。推荐作物主干是：

```text
主干：DSSAT-CSM 作物系统
辅证/接口参考：DAISY
全球 PFT 外壳参考：CLM5 + LPJ-GUESS
```

不要把 CLM5 crop 或 LPJ-GUESS crop 当成 COMPASS 的最终作物主干。它们适合提供 land-model 级别的 sowing/harvest/date/heat-unit 外壳，但不适合承担“特定作物可用”的产量、生物量、器官分配和 N 胁迫主模型。

### 5.2 三个候选源的定位

#### DSSAT-CSM：主干

优点：

- 专门做特定作物。PFT15-22 这类 corn/wheat/soybean/cotton/rice/sugarcane，与 CSM 的作物模型族直接对应。
- 有完整 daily crop loop，不只是 phenology：
  - `PHENOL` / `CSP_PHENOL`：物候、成熟日期、阶段。
  - `DEMAND`：器官需求。
  - `NUPTAK` / `CSP_NUPTAK`：作物 N uptake。
  - `VEGGR` / `GROW` / `CSP_GROW`：器官生长、分配、生物量、LAI。
  - `SENES` / `CSP_SENES`：衰老和残体。
  - `ROOTS` / `CSP_ROOTS`：根系生长。
  - `AUTHAR` / `OPHARV` / `CSP_OPHARV`：成熟、收获、产量、残留。
- 输出变量能对接常见农业观测：yield、biomass、LAI、grain N、stover、root、stage、harvest date。

缺点：

- Fortran 代码复杂，变量体系是 DSSAT 风格，不适合“直接整包塞进 MATLAB”。
- 需要做一层 COMPASS crop state adapter，把 CSM 的 `WTLF/STMWT/RTWT/SDWT/XLAI/MDATE/SENESCE/HARVRES` 映射到 COMPASS 的 C/N pools 和 diagnostics。

判断：这是唯一适合做 COMPASS 作物生产主干的候选。

#### DAISY：接口与管理参考

优点：

- 土壤-作物-氮-水耦合接口很清楚。
- `Crop::tick` 明确把 bioclimate、soil、soil water、chemistry、organic matter、residual C/N 串起来。
- harvest 设计好：把 stem/leaf/storage organ 分成 harvested yield 和 field loss，并把 C/N/DM 残留送回 topsoil/soil layers。
- `CropStandard` 的组合结构清晰：seed、root、canopy、harvesting、production、phenology、partition、vernalization、photo、crop N。

缺点：

- 更像可配置 generic crop framework，不是 DSSAT 那种作物族专用模型库。
- 若直接拿 DAISY 主干，PFT15-22 的特定作物参数和产量过程仍要另补。

判断：DAISY 不做主干，但要借它的接口思想，特别是 harvest/residual C/N 返回、root density 与 soil-N/water 连接、管理动作结构。

#### CLM5 / LPJ-GUESS：全球陆面外壳参考

优点：

- 与 COMPASS 当前 PFT/land-model 框架更接近。
- CLM5 有 crop seed、grain、harvdate、cleanup；LPJ-GUESS 有 sowing/harvest date、development stage、harvest index、double-crop 管理。
- 适合给 COMPASS PFT 表和大尺度默认逻辑提供最低限度框架。

缺点：

- 作物机理相对粗，目标不是把每个具体作物做成 DSSAT 级别。
- 只抄 CLM5/LPJ 会得到“能收获的全球 crop PFT”，但不是“特定作物模型”。

判断：CLM5/LPJ 只做外壳和 land-model 一致性参考，不做最终作物机理主干。

### 5.3 最终方案

采用分层方案：

1. COMPASS 保留自己的主循环、土壤水热氮、辐射、SPAC、photosynthesis 框架。
2. 作物内部生长主干按 DSSAT-CSM 移植成 MATLAB crop adapter。
3. DAISY harvest/residual 接口用于约束 C/N 守恒和残体返回。
4. CLM5/LPJ 只负责 PFT 表、global crop PFT 默认状态、sowing/harvest 外壳、seed/grain/cleanup sanity check。

### 5.4 建议新模块边界

不要把作物塞进现有 `cropPhenology.m` 一个函数。建议新增：

```text
src/+plant/+crop/
  runCropStep.m                 主调度，对应 CSM daily crop loop
  initCropState.m               播种/出苗/季节初始化
  cropPhenologyCSM.m            PHENOL / CSP_PHENOL 映射
  cropDemandCSM.m               DEMAND 映射
  cropGrowthCSM.m               GROW / VEGGR / CSP_GROW 映射
  cropNitrogenUptakeCSM.m       NUPTAK / CSP_NUPTAK 映射到 COMPASS soil N uptake
  cropSenescenceCSM.m           SENES / CSP_SENES 映射
  cropRootCSM.m                 ROOTS / CSP_ROOTS 映射到 root density/root pools
  cropHarvestCSM.m              AUTHAR / OPHARV / CSP_OPHARV 映射
  mapCropPFTToCSM.m             PFT15-22 -> CSM crop family
```

`src/+plant/+integrated/+phenology/cropPhenology.m` 只能保留为 dispatcher，不再承载完整作物生长。

### 5.5 PFT15-22 到 CSM 主干的建议映射

第一版映射：

- PFT15 corn：CSM CERES-maize / CSCER 类主干，或 CSM maize-specific file if present。
- PFT16 wheat：CSCER/CSCRP wheat/barley family。
- PFT17 soybean：CROPGRO grain legume family。
- PFT18 cotton：CROPGRO cotton family。
- PFT19 rice：ORYZA or CSM rice family if the local source path provides the active interface。
- PFT20 sugarcane：CSP / CASUPRO family，`CSP_PHENOL`, `CSP_GROW_CANE`, `CSP_OPHARV`。
- PFT21 tropical corn：same corn model with tropical cultivar/ecotype parameters。
- PFT22 tropical soybean：same soybean model with tropical cultivar/ecotype parameters。

Claude 必须先确认本地源码中对应 crop family 的 active interface，不要只按文件名猜。

### 5.6 第一阶段最小可用作物

不要一口气移植 8 个作物。第一阶段只做 2 个代表：

1. Corn：代表 cereal/grain crop。
2. Soybean：代表 CROPGRO legume/N-fixation crop。

若必须覆盖水稻/甘蔗，再单独作为第二阶段，因为 rice/sugarcane 的物候和生物量结构与普通 grain crop 不同。

第一阶段最低闭环：

1. sowing/emergence/maturity/harvest state。
2. development stage and stage dates。
3. leaf/stem/root/grain or storage-organ C/N pools。
4. LAI from leaf biomass or crop-specific LAI state。
5. root depth / root density passed to COMPASS water and N uptake。
6. N demand and uptake from COMPASS soil N。
7. senescence to vegetation litter/residue diagnostics。
8. harvest yield C/N removed from system。
9. non-harvested residue C/N returned to existing residual/litter path。
10. post-harvest cleanup so annual pools do not carry unrealistic biomass into next season。

### 5.7 临时护栏

在上述最小闭环完成前，PFT15-22 不能静默运行 placeholder。可以 fail-fast，但错误信息必须写成临时护栏：

```matlab
error(['COMPASS:CropCoreNotSelected', ...
       ' crop_phenology=1 requires the DSSAT-CSM crop adapter. ', ...
       'The old placeholder cropPhenology.m is not a production crop model.']);
```

这不是长期禁用作物，而是防止当前 placeholder 伪装成可用结果。作物是 COMPASS 收尾必做模块。

### 5.8 验收

方案级验收：

- 文档和代码注释明确：CSM/DSSAT 是 crop process 主干，DAISY 是 harvest/residual/interface 参考，CLM5/LPJ 是 land-model 外壳参考。
- PFT15-22 不再被描述为“可先不用”，而是被列为收尾必需模块。

第一阶段代码验收：

- Corn 和 soybean 至少能完成完整 sowing-to-harvest season。
- 输出 yield、harvest C/N、residue C/N、LAI、root depth/root density、N uptake、stage date。
- 收获后 annual crop leaf/stem/root/grain pools 清零或进入明确残留/harvest sink。
- C/N 守恒：harvest removal + residue return + remaining pools 与生长输入闭合。
- PFT23 apple 不受 crop adapter 影响。

## 6. Patch C：通用 gap mortality 参数覆盖

### 6.1 当前状态

`runPlantCNStep.m` 已有 Step 10 background gap mortality：

- 读取 `r_mort = vegetation.getPFTValue(pftType, 'r_mort', 0)`
- 对 `C_livestem/C_deadstem/C_livecroot/C_deadcroot` 和对应 N pool 按同一 fraction 扣减
- 累加到 `cum_C_mortality/cum_N_mortality`

这条逻辑主要服务 PFT23，因为 PFT23 有 `r_mort=0.02`。其他 woody PFT 如果没有 `r_mort`，就默认 0。

### 6.2 是否要改公式

第一版不改公式，只补参数覆盖和诊断说明。

CLM5 原版 mortality 包括 leaf/froot/live/dead wood；COMPASS 已经有 leaf/froot phenology turnover，若再把 leaf/froot 放入 gap mortality，可能与 bglfr 双算。第一版只维持 wood mortality 是可接受的，但必须在注释里写清楚：

- CLM5 source includes leaf/froot and wood mortality.
- COMPASS maps leaf/froot turnover through phenology/bglfr.
- COMPASS Step 10 only handles perennial woody structural mortality.

### 6.3 具体步骤

1. 在 `getPFTParameterSet.m` 为 woody PFT 添加或确认 `r_mort`。

建议范围：

- PFT1-11 woody natural PFT
- PFT23 apple
- PFT24 Robinia

不要给 grass/crop 添加 wood mortality。

2. 参数来源优先级：

   1. CLM5/FruitTree 参数表中已有 `r_mort` 或 equivalent annual mortality。
   2. CLM5 default `am` / gap mortality setting。
   3. 如果没有可靠来源，先不填，保持 0，并在文档/注释里标 `unverified`。

3. 不要为了把 wood pool 压到目标值而随意填大。

4. 在 Step 10 注释里补一句：

```matlab
% COMPASS applies gap mortality to woody structural pools only.
% Leaf/froot exits are handled by phenology/bglfr turnover to avoid double counting.
```

5. 如果 diagnostics 注册了 `cum_C_mortality/cum_N_mortality`，确认 PFT1-11/24 开启后能输出。

### 6.4 验收

- PFT23 `r_mort=0.02` 仍存在。
- 非 woody PFT 不应有 wood mortality。
- Woody PFT 如果没有 verified 参数，不能偷偷用 apple 参数。
- 多年 woody run 中 `C_deadstem/C_deadcroot` 不应无限单调积累而无出口。

## 7. Patch D：诊断和防双算审查

### 7.1 必查双算点

每个 patch 完成后，检查以下流是否重复扣：

- `diag.leaf_litter_today`
- `diag.root_litter_today`
- `diag.stem_turnover_today`
- `cum_C_mortality`
- `cum_N_mortality`
- `turnoverCarbonMirror.m` 对 C pool 的镜像扣减
- `COMPASS.m` 里 plant residual / litter / mortality 汇总

特别注意：

- 如果 background bglfr 已经加到 `diag.leaf_litter_today`，C mirror 不应另有单独 bglfr 扣 C。
- 如果 crop adapter 尚未完成，annual early return 不得 silently skip；必须临时 fail-fast。adapter 完成后，annual PFT 必须进入 `+plant/+crop/runCropStep.m`。
- 如果 PFT23 不走 bglfr，apple autumn senescence 不应变。

### 7.2 推荐新增诊断

如果现有 diagnostics 架构允许，增加：

- `plant_leaf_background_litter_N`
- `plant_root_background_litter_N`
- `plant_bglfr_turnover_frac`

如果 C mirror 可以提供 C 侧：

- `plant_leaf_background_litter_C`
- `plant_root_background_litter_C`

第一版不强求新增 netCDF 变量；但至少 `turn_diag` 内部字段要存在，便于调试。

## 8. 最小测试计划

### 8.1 静态检查

运行：

```powershell
rg -n "leaf_long_days|bglfr|leaf_background_litter|root_background_litter|CropCoreNotSelected|runCropStep|cropHarvestCSM|r_mort" src
```

验收：

- `leaf_long_days` 出现在真实 turnover 计算中。
- crop placeholder 不再静默 production run；若 CSM adapter 未完成，应能搜到临时护栏；若已完成，应能搜到 `runCropStep` 和 CSM crop adapter。
- `r_mort` 注释说明 wood-only mapping。

### 8.2 MATLAB 语法检查

优先用项目已有 MATLAB batch 方式。示例：

```powershell
& 'F:\Matlab\bin\matlab.exe' -batch "cd('F:\COMPASSV1.0'); addpath('src'); which vegetation.integrated.turnoverAndHarvest; which vegetation.integrated.phenology.cropPhenology;"
```

如果有现成 test runner，优先跑现成 runner。

### 8.3 PFT23 回归

用已有 Changwu apple 配置跑短期或已有 regression：

验收：

- Apple run 不应进入 crop adapter 或触发 crop 临时护栏。
- harvest/senescence/reserve 不崩。
- C/N conservation 无明显新 residual。
- PFT23 不应因为 Patch A 开始全年持续背景落叶，除非用户另批。

### 8.4 Evergreen smoke test

构造或切换一个 evergreen PFT，例如 PFT1：

验收：

- `phen=1`
- `leaf_background_litter_today > 0`
- `root_background_litter_today > 0`
- 每步转出比例约等于 `dt/(leaf_long_days*86400)`。

### 8.5 Grass smoke test

构造或切换 PFT13：

验收：

- 生长季 `phen>0` 时 leaf/froot 有背景 turnover。
- 如果 `phen=0` 且使用 phen gate，background turnover 应接近 0。
- 不需要 storage pool 也能有 leaf exit。

### 8.6 Crop adapter test

分两种状态验收：

未完成 CSM adapter 时：

- 选择 PFT15 时必须明确报 `CropCoreNotSelected` 或等价临时护栏。
- 错误信息说明作物收尾必须走 DSSAT-CSM adapter，不能继续 placeholder。

完成第一阶段 CSM adapter 后：

- PFT15 corn 可完成 sowing-to-harvest season。
- PFT17 soybean 可完成 sowing-to-harvest season。
- 两者输出 yield、harvest C/N、residue C/N、LAI、root depth/root density、N uptake、stage date。
- 收获后 annual crop pools 不应跨年残留未解释生物量。

## 9. Claude 执行顺序

按以下顺序改，不要并行乱改：

1. 阅读本文和第 2 节列出的 COMPASS 文件。
2. 对照 CLM5 源码确认 `bglfr`、mortality 行号；对照 DSSAT-CSM/DAISY 源码确认 crop 主干和接口边界。
3. 实施 Patch A：`bglfr` background leaf/froot turnover。
4. 跑 PFT23 短测，确认 apple 不受影响。
5. 实施 Patch B：先写 CSM/DAISY/CLM5-LPJ 作物主干选择说明和临时护栏；若进入代码实现，则从 corn + soybean 第一阶段 adapter 做起。
6. 实施 Patch C：补 woody `r_mort` 参数覆盖或明确留空 unverified。
7. 补必要注释和 diagnostics。
8. 跑第 8 节测试。
9. 输出改动摘要，必须列出：
   - 改了哪些文件
   - 哪些 PFT 行为变了
   - 作物 adapter 是临时护栏状态还是已实现 corn/soybean 第一阶段
   - PFT23 回归是否通过
   - 哪些参数来源仍 unverified

## 10. Codex 审查清单

Claude 改完后，Codex 按以下清单审查：

1. `storage/xfer` 池结构是否被改动。若被改动，除非用户明确批准，否则退回。
2. PFT23 是否被意外纳入 background bglfr。
3. `leaf_long_days` 是否真实驱动 leaf/froot litter，而不是只改注释。
4. `turnoverCarbonMirror.m` 是否发生 C 双扣。
5. crop PFT 是否仍在 placeholder 假跑。未完成 CSM adapter 时必须 fail-fast；完成后必须能跑 corn/soybean season。
6. `r_mort` 是否只给 verified woody PFT，是否误给 grass/crop。
7. 所有新增公式是否有 CLM5/LPJ/DSSAT-CSM/DAISY 源码依据。
8. 是否触碰 SOM。若触碰，单独拎出来，不并入本补丁。
9. MATLAB 是否至少通过语法/短测。
10. 结果解释是否仍围绕“生命周期出口补全”，没有伪装成 CUE/NEE 调参补丁。

## 11. 一句话边界

本补丁的科学边界是：补齐 CLM5/LPJ/CSM/DAISY 已有的 vegetation 与 crop 生命周期出口，使非苹果 PFT 不再“只进不出”，并把作物明确推进到 DSSAT-CSM 主干移植路线。它不是 SOM 改造，不是 storage/xfer 池重构，也不是用参数硬追观测结果。
