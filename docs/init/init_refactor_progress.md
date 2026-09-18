# 初始化重构进度备忘

日期: 2026-06-22

## 当前目标

逐步替换 `run_stemmus_initialize`，把站点文件夹输入整理成清晰的四层流程：

```text
io.loadSite(siteDir) -> io.checkSite(raw) -> init.buildState(raw) -> init.initialize(siteDir)
```

责任边界固定为：

- `io.loadSite`: 只读 `site.toml` 和 CSV，保留原单位，不做物理加工。
- `io.checkSite`: 只做结构、必需列、枚举、行数校验。
- `init.buildState`: 把 `raw` 加工成模型可用的 `built`。
- `init.initialize`: 未来顶层编排入口，替代老初始化脚本的 base workspace 输出。

## 已完成并验证

`buildState` 当前阶段：

```text
CFG init.buildConfig(raw.config)                           规范化配置 -> built.switches/numerics/grid/env/time
    (built.switches/numerics/grid/env/time 直接进 built 顶层, 不套 built.config)
A  init.buildDomain(built.grid)                            网格(用规范化 grid)
B  init.buildSoilParameters(raw.soil, domain)              VG/质地/N运移速率/派生
BD init.buildBulkDensity(raw.initial, soil, domain, cfg)   定容重来源 -> bd.node/layer/porosity
   -> soil.bulk_density = bd.layer; soil.porosity = bd.porosity
B2 init.buildMoistureChar(soil)                            特征含水量(field cap/wilt/air-dry)
B3 init.buildSoilThermal(soil)                             静态热参(需 B2 wilting_point + bd)
C  init.buildInitialState(raw.initial, built, cfg, bd)     初值 Y(用同一份 bd.node)
```

> 注：阶段标签按执行序排（B2 先于 B3，因为热参 GB2 依赖 B2 的 wilting_point）。
> `BD` 是**前置**容重决策：一次定死来源（`bulk_density_source`），soil/thermal/Y/N
> 全吃同一份 bd —— 取代了早先"事后补救"的 `reconcileBulkDensity`（已删）。

验证结果：

- `domain` 与老 `Dtrmn_Z` 逐元素一致。
- `built.parameters.soil` 已等价老 `core.SoilProperties` 的全部 soil 参数：VG、质地、N 运移/速率、派生量、特征含水量、热参。
- `BD`（buildBulkDensity，前置）：默认 `bulk_density_source="initial"` 用 IC 容重（复刻老 [7b]：节点→层均值 clamp、porosity=1-bd/2.65 clamp）；`"soil"` 用 soil.csv 容重。soil/thermal/Y/N 共用同一份 bd，无事后 reconcile。
- `test_buildState`: 76/76 PASS（soil gate-A 用 `bd_source="soil"` 对 standalone SP，+ BD `initial` 路径单测）。
- 读取/校验测试通过：
  - `test_readToml`: 15/15
  - `test_readCsvTable`: 11/11
  - `test_readExcel`: 9/9
  - `test_loadSite`: 11/11
  - `test_checkSite`: 12/12
- `src/+init` 当前相关函数和 `src/+io/checkSite.m` 已 `checkcode` CLEAN。
- 根目录空残留 `buildState.m` 已删除，MATLAB 只解析 `init.buildState`。

## 当前重要设计点

`buildInitialState` 已把 N 初值改为运行时需要的液相浓度：

- `nitrogen_input_phase = "total"`: `soil_initial_profile.csv` 的 NH4/NO3 当作土壤总量 `[mg/kg]`，转换成液相浓度。
- `nitrogen_input_phase = "dissolved"`: 输入直接当液相浓度 `[mg/L]`。

这与当前 N 模块的 storage 公式一致：

```text
NH4 storage = (theta + RHOKG) * C_NH4
NO3 storage = theta * C_NO3
RHOKG = bulk_density * KG
```

新增配置开关（site.toml [switches]，checkSite 已加枚举校验）：

- `nitrogen_input_phase = "total" | "dissolved"`（N 初值相位）。
- `bulk_density_source = "initial" | "soil"`（容重来源，控制 R 步是否覆盖）。

## built 容器约定（2026-06-22 重定）

**核心原则**：`parameters` 只放**物理/生物过程参数**，不是"所有不变量"。运行控制、几何、
边界、管理源各有自己的顶层家（学 ClimaLand：parameters 不含 time/grid/control/boundary）。
**单写者**：每个容器**只由一个 builder 写**。

| 顶层容器 | 内容 | builder |
|---|---|---|
| `built.domain` | 网格几何(dz/NL/NN/z_node/sumdz) + `.grid`(网格配置) + `.site`(lat/lon/elev) | buildGrid+buildDomain |
| `built.control` | `.switches` / `.numerics` / `.time`（求解+运行控制，非物理） | buildSwitches/buildNumerics/buildTime |
| `built.parameters` | `.soil`（后续 `.veg` / `.scope`）—— 仅过程参数 | buildSoilParameters... |
| `built.boundary` | 顶/底 BC 类型码 + 值（NBCh_top...，flat legacy 名） | buildBoundary |
| `built.sources` | `.fertilization`（管理源汇；后续灌溉等） | buildSources |
| `built.env` | **过渡容器**：环境标量待按角色再分配（见下） | buildEnv |
| `built.Y` | prognostic state | buildInitialState |
| `built.p` | aux cache（待） | — |
| `built.forcing` | 驱动时间序列（待 buildForcing） | — |

- **过渡 alias**（老主循环短路径，迁移完删）：`built.switches/numerics/time == built.control.*`，
  `built.grid == built.domain.grid`。
- **`built.env` 是临时停放**：字段已查代码角色 —— `Tav/At`→热底边界
  ([calculateBottomTempBoundary.m:39](../src/+soil/+thermal/calculateBottomTempBoundary.m#L39))、
  `Rho_Thermal_Soil_*`→土壤热(表层热阻 COMPASS.m:870)、
  `Z_measure/Wind_*/Rin_unit_is_watt + Ca/Oa`→forcing(MeteoData/光合)。
  buildForcing/buildVeg 落地时各归各位，env 容器随之消失。site 已镜像到 `built.domain.site`。

  > **built.env is transitional.** It exists only until buildForcing /
  > thermal-boundary / soil-runtime thermal builders consume and relocate its
  > fields. Do NOT treat `built.env` as final architecture.

- `raw.config` 只保留原始 TOML，不进计算逻辑。

## 代码组织约定（2026-06-22 定）

- 命名：camelCase 动宾，无人工前缀（私有性靠 file-local / `private/` 表达）。
- 公共阶段函数：`init.buildXxx`。
- init 层共享小工具：`src/+init/private/*.m`（+init 下所有函数可调，外部不可见）：
  `getSwitch` / `getSwitchString` / `getTableColumn` / `clampWithWarn` /
  `mapLayerToNodeMean` / `findNearestNode`。
- 物理公式专用 helper：留各自 file-local（如 `mapLayerToHorizon` / `classifyTexture` /
  `recomputeHeadFromTheta` / `validatePhysical` / `moistureViaMatricHead`）。
- 两种层→节点映射**故意并存**：`mapLayerToNodeMean`（简单平均，bd legacy 复刻）≠
  `equations.interpolateElementToNode`（厚度加权，RHOKG/物理场）。
- 边界：`io` 不依赖 `+init/private`（`io.checkSite` 暂留自己的 `getSwitch`；跨层配置访问器
  以后再提 `core.config`，不急）。
- 注释：local/private helper 一行用途 + 单位/出处，不写大段 docstring。

## config 家族拆分 (2026-06-22 更新)

- **`buildConfig` 已拆为 5 个单写者 builder**（删 `buildConfig.m` + `private/deriveTime.m`）：
  `init.buildSwitches`(→control.switches，含级联+enumToInt) / `init.buildNumerics`(→control.numerics) /
  `init.buildGrid`(→domain.grid，几何仍归 buildDomain) / `init.buildEnv`(→built.env 过渡) /
  `init.buildTime`(→control.time，折进原 deriveTime)。
- **`buildBoundary` 拆出 `buildSources`**：BC 码留 buildBoundary(→built.boundary)，施肥事件移到
  `init.buildSources`(→built.sources.fertilization)。
- 测试随之拆：`test_buildSources` 8 项（原 buildBoundary 的施肥段），`test_buildBoundary` 24 项(纯 BC)，
  `test_buildState` +6 容器布局守卫(=82)。`test_buildConfig` 41 项改为分别调 5 个 builder 组装。

历史细节（仍适用）：
- 开关 §8 名字映射 + enum→整数(代码验证: tian→3, peters_lidard/johansen→2, farouki→4; prunty→1; van_genuchten→1)+ setup_conditional_switches 级联(lSoilairefc=0→{lKaT,lKaa,lDVa,lKLa}=0; lThmrlefc=0→{lKLT,lDVT}=0)。
- 时间派生 `init/private/deriveTime.m`：start/end 解析 + TEND/tS/Nmsrmn(永远自动算+10%)/JN/lTimeScale。
- gate-A `test_buildConfig` 41/41(对老 ModelSettings 直比，含 lThermCond/时间）。
- **配置已同步**(2026-06-22)：site.toml 与 config_*.txt 语义一致 —— 热导率 `de_vries`(lThermCond=1)、整年(EndTime=2022-01-01)。规则：当前阶段两边保持一致；长期 site.toml 为唯一输入，config_*.txt 退役或仅作 gate-A reference。新增枚举 `de_vries→1`(原 toml 枚举漏了 de Vries 默认)。
- **运行期接线提醒**：tol_water/tol_energy 已进 `built.control.numerics`；老 `COMPASS.m`
  现用局部 `Tol_Water/Tol_Energy`，不是 `xERR/TERR`。新 `runOneStep`/主循环接入时
  应直接读 `built.control.numerics`。

## buildBoundary 完成 (2026-06-22)

- `init.buildBoundary(raw.config)` → `built.parameters.boundary`(canonical, 无 alias)。
- BC 类型码 + Dirichlet/通量值 friendly→legacy 名(NBCh_top/NBCT_*/NBCP_*/NBC_{NH4,NO3,UREA}_* + 值),全 passthrough;urea BC 用 BoundarySettings 默认(2/3/0/0,site.toml 不设 urea BC)。
- 施肥事件:site.toml 数组 → `bc.fertilization`(events_time/amount/type/ratios + depth_start/end);单位 kgN/ha(与 config_boundary.txt 一致)。
- site.toml 补了 T_top/P_top(占位,type 1/3 运行时覆盖)+ schema。
- gate-A `test_buildBoundary` 32/32(24 BC 字段 + 8 施肥;施肥解析老分号串比对)。

## forcing 读半完成 (2026-06-22)

驱动格式定为 **forcing.csv 命名列**(发布版 nc 留后)。读/校验这半已通,物理处理留给 buildForcing。

- **forcing.csv 生成**:从 `Meterology data2021-2023.xls` 一次性转(脚本用完即删,溯源写在 csv 头注释)。
  730 行日尺度(2021-2022),21 命名列,**保留 Excel 原值原单位**(rh=%, pg=hPa, rn=MJ/m²),
  缺测列(pg/rn/ts/tb/h_sur/wind_u/v 长武未测)写 `NaN`→读成 numeric NaN。单位换算/缺测填补
  全部留到 buildForcing(搬 MeteoData.load),gate-A 才干净。
- **tableSpec('forcing')**:required={date,tmax_c,tmin_c,tmean_c,rh_mean_pct,wind_ms,precip_mm,
  rin_wm2,rli_wm2};optional 含 ts/tb/pg/rn/sunshine/ca + 未来 N 沉降/底BC 4 列。按列名校验。
- **loadSite**:`raw.forcing` 改读 forcing.csv(替原 .nc TODO 占位)。
- **checkSite**:加 forcing 为 mandatory wide-table 校验。
- 测试:`test_loadSite` 12(+1 forcing 730 行)、`test_checkSite` 15(+2 forcing present/缺列报错)。

## buildForcing(物理半)— 搬 MeteoData，分 F1-F4 每步 gate-A

输出结构(用户拍板):`built.forcing.{raw, hourly, derived, solar}`,flat alias 按需再加。
gate-A 参照源 = 保留的 `Meterology data2021-2023.xls`(站点文件夹内),`core.MeteoData.load(Excel)` 取 `*_msr`。

- **F1 raw ✓ (2026-06-22)**:`init.buildForcing(raw.forcing)` → `forcing.raw`(legacy `*_msr` 等价层)。
  命名列→base 名(Tmax/Ta/RH/Pg/.../NO3_bot)+ 单位换算(RH×0.01、Pg hPa×100→Pa)+ **缺失 N 列补 0**
  (非 NaN,legacy load 行为);present-but-empty 列(pg/ts/tb/wind_u 长武未测)保持 NaN。
  **不碰小时展开**。forcing.csv 改 `%.15g` 全精度,round-trip Excel double。
  gate-A `test_buildForcing` **25/25**(24 字段 NaN-aware + 行数,逐元素对 MeteoData.load)。
- **F2 hourly ✓ (2026-06-22)**:`forcing.hourly` = 搬 `to_hourly`。日→时正弦气温(峰15:00)/地表温(峰16:00)、
  RHmax/RHmin 正弦或 ea 守恒湿度、自适应风速日较差、FAO-56 短波(`surface.radiation.distribute_radiation_fao56`)、
  SWAT 降水降尺度(`equations.downscale_rain_swat(.,6,0.35,10)` 全参→确定性,非随机)、N 沉降/24+底BC平铺。
  风速高度换算复用 `core.MeteoData.convert_wind_height`(吃 built.env)+ 最小 1 m/s;含 start_time 偏移截断。
  两 lTimeScale 分支都搬。gate-A `test_buildForcing` **50/50**(F1 25 + F2 25,逐元素对 to_hourly)。
  `buildForcing(rawForcing, built)`:built 缺省只出 F1;给了才出 F2(读 control.time/domain.site.latitude/env)。
- **F3 derived ✓ (2026-06-22)**:`forcing.derived` = 搬 `compute_derived`。es/Δ/qs(`equations.Calc_Esat_Flatau`)、
  ea/P_Va、qa、RHOV_A、位温/虚位温(`core.PhysicalConstants`)。**Pg 缺测时按高程 barometric 估算**(长武 pg 全 NaN→触发,
  须先于 Esat),有效 Pg/Ca 存进 derived(hourly.Pg 保持 to_hourly 原值不动)。gate-A +11(61/61)。
- **F4 solar ✓ (2026-06-22)**:`forcing.solar` = 搬 `compute_solar`。逐小时 tts/sun_altitude/azi/sun_R + 日出/日落/
  日长(每日缓存),复用 `surface.radiation.calculatesunpos_ultimate_azi_atan22`(method 0)。gate-A +7(**68/68**);全年 17520 步 7s。

**buildForcing 完成**:F1-F4 全过,`test_buildForcing` 68/68 逐元素对 `MeteoData`(load/to_hourly/compute_derived/
compute_solar)。已接进 buildState:`built.forcing = init.buildForcing(raw.forcing, built)`(用 domain.site lat/lon/elev、
control.time、env)。test_buildState 82 仍绿(含 forcing 构建,+7s)。

## 时间基重构 S1-S5 (2026-06-25)

时间基重构 S1-S5 已完成。S4 温度与 S5 测试已并入 S3 实现。

落地结果:
内部时间统一为 UTC;forcing 的本地时区只在读取口转换;太阳几何只使用 UTC + 经纬度,不再使用 round(lon/15)。

S1 buildTime:
三时基 {start,end,forcing_start}_{clock,utc};
clock_basis / utc_offset 接入 site.toml、schema、checkSite;
clockToUtc 接通。

S2 buildTimeAxis:
统一生成逐步 utcTime / clockTime / utcHour / clockHour / clockDayIndex;
trueSolarHour = mod(hour_angle/15 + 12, 24)。

S3 buildSolar + buildHourly:
buildSolar 使用 timeAxis.utcTime + lat/lon,输出 hour_angle / solarHour;
sunrise_anchored 温度使用 clockDayIndex + clockHour + sunrise_clock / sunset_clock;
dispatcher 顺序为 F1 -> buildTimeAxis -> buildSolar -> buildHourly -> derived。

关键修正:
长武太阳几何不再被 round(lon/15) 错算成 UTC+7;
sunrise_anchored 的 day/hour 统一到 clock 基,凌晨 00:00 正确归入当天前夜段;
F4 不再对老 MeteoData 做 bug 级复刻,改为物理判据;F1-F3 仍保留 legacy_sine 对 MeteoData 的逐元素 gate-A。

验收:
test_buildForcing 74/74;
当前完整输入层回归 304/304;
checkcode CLEAN;
jsonschema VALID。

遗留但不阻塞:
built.env 的 Z_measure/Wind/Ca/Oa 仍作为过渡容器被 buildForcing 消费,待 forcing/env 最终归位后移除;
start > forcing_start 时 timeAxis/solar 的截断对齐还需补测;
全球日尺度 Tmax/Tmin 的"日"归属仍需单独设计,全球优先使用小时 forcing。

> 注:forcing 链对**长武当前 `start == forcing_start` 的单点日尺度 forcing** 完整;
> `start > forcing_start`(offset>0)与全球 forcing 情况未触发验证(见上"遗留")。

## buildVeg I/O 契约 (锁定 2026-06-25, 实现前定)

分工:**buildVeg** = 植物长什么样 / 叶片怎么光合 / 冠层·根·水力参数 / 物候状态初值;
**buildInitialPlantPools** = C/N 初池(含 hydraulic_init_psi*);**buildScope** = 光谱/辐射;
**C/N 开关参数** → 单独 `built.parameters.plant_cn`。

**输入(消费)**:`raw.config.vegetation` + PFT 表(`vegetation.getPFTParameterSet(PFT_type)`)合并;
`built.domain`(根系按层);`raw.root`(root.csv,root_dist_mode=1);`built.control.time`(日序/树龄);
`built.forcing`(逐日 Ta/Tmin/Tmax → model 物候积温)。**buildVeg 在 buildForcing 之后**。

**输出**:
```
built.parameters.veg.{structure, canopy, phenology, hydraulic, root, photosynthesis}
built.vegetation_state.{phen_driver(forcing-only 预算), lai_state(主循环逐日推进), hydraulic_daily}
built.parameters.plant_cn.{cn_coupling_mode, cn_flex_offset, allocation, respiration}   ← 后续 builder
```
- 光合/气孔(photo_type/vcmax_source/vcmax_manual/luna/photo_Rd_*/leaf_resp_method/gs_*)→ `veg.photosynthesis`(叶生理,不进 scope)。
- prognostic 暂放 `built.vegetation_state`,**不进 built.Y**(主循环未统一 solver;built.Y 暂为土壤态)。以后重构再迁 `built.Y.{soil,veg,plant_cn}`。
- **LAI 不预算**:model 模式(长武)LAI 依赖水力胁迫+主循环态,逐日在循环算;buildVeg 只出物候参数+初值。constant→标量;data→预设时序。
- `cn_coupling_mode/cn_flex_offset` 归 `plant_cn`,不进 veg。

**实现分步(照 forcing F1-F4 打法,每步 gate-A)**:
- **V1 ✓ (2026-06-25)**:静态参数 → `parameters.veg.{structure,canopy,phenology,hydraulic,root,photosynthesis}`。
  site.toml [vegetation] **friendly 名 → legacy VegetationData 属性**映射 + enum/bool 转换(leaf_width→w_leaf_max、
  psicor→Psicor、vcmax_source→Photo_Vcmax_flag(manual=2)、luna→lLUNA、leaf_resp_method→Photo_LeafResp_Method、
  retention_slope_*→dyn_a_retention_*、path_*→dyn_dz_*...)。**PFT 表不进 V1**(它是 plant_cn/池的料;长武 vcmax_source=manual)。
  `built.parameters.veg`(暂 standalone,接 buildState 待 V4 完)。gate-A `test_buildVeg` **82/82** 对 VegetationData
  (load_from_txt config_vegetation.txt,已同步)。`root_measured_file` 不比(路径机制 new=raw.root vs legacy=路径串;数据等价由 V2 验)。
- **V2 ✓ (2026-06-25)**:`buildVeg(raw, built)` 加 `veg.root.RAI`[NL×1, bottom-up]。
  Jackson 路(root_dist_mode=0)**复用** `vegetation.hydraulics.calculateRootAreaIndex1`,gate-A 逐元素对它(验接线)。
  measured 路(root_dist_mode=1,长武)按 **HEAD-correct** 逻辑(RL=domain.tot_depth、Lmax=root_max_depth、
  raw.root RLD 重叠加权→层 → `RAI=2π·r_root·RLD·dz`)。**注:measured 暂不逐元素 gate-A 工作区 `calculateRootDistribution`**
  ——该函数工作区版 measured 分支裸用未定义的 RL/Lmax(HEAD 有定义,被 git-M 编辑删了,疑 Jizixing Zhu 重构 root_dist_mode);
  改用物理判据(RAI 有限/非负、sum 合理、随 raw.root RLD 线性=确实消费 root.csv);待该函数稳定再补逐元素。
  gate-A `test_buildVeg` **88/88**(V1 82 + V2 6)。不碰 git-M 的 calculateRootDistribution。
- **V3 ✓ (2026-06-25)**:`init.buildVegState(built, veg)` → `vegetation_state.phen_driver`(model 模式)。
  日 GDD 积温曲线复用 `vegetation.phenology.CalculateHeatUnit`(== init_phenology),吃 `forcing.raw.Tmax/Tmin`(日,== legacy *_msr);
  出 phen_CHU/HU/HUI/dHUF + phen_Tavg/Tmin/Tmax_daily。**LAI 不预算**(model 逐日靠水力胁迫,留主循环)。
  gate-A `test_buildVegState` **8/8** 对 init_phenology。
- **V4 ✓ (2026-06-25)**:`vegetation_state.lai_state` —— LAI/hc/SAI 步序列(长度 Nmsrmn,搬 to_hourly:
  constant 常数 / model=ALAI_MIN 占位[主循环逐日覆写] / data 暂不支持[InPlant 未入新站])+ 逐日 LAI 状态初值
  (lai_day_prev/dormant/LAI_D/HUI_D=DLAI/yr + lai_model_current,搬 init_day_state)。**model LAI 日更留主循环**(不碰)。
  gate-A `test_buildVegState` **17/17**(V3 8 对 init_phenology + V4 9 对 to_hourly/init_day_state)。

**buildVeg 完成(V1-V4)+ 已接入 buildState**:`built.parameters.veg`(buildVeg(raw,built))+ `built.vegetation_state`
(buildVegState(built, veg))。test_buildState 84(+2 veg 布局守卫)。`compute_model_lai` 整段 LAI **不搬**(逐日靠水力胁迫,主循环)。

## D-plant 契约 (锁定 2026-06-25, 实现前定)

`setInitialNitrogen` 是**整个 N 模块**初始化(4 块:①土壤 N 静态参→buildNitrogenSoilParams ②N 求解器态+③tS 诊断数组→
allocateWorkspace ④**植物 C/N 池**=D-plant)。`Init_Plant` 是植物**水力**(参数已在 veg;AUX=VC/KS;STATE=psi*)。
拆**三个 builder**,顺序 **P1→P3→P2**(P2 需 P3 的 CN/alloc/storage 参数):

```
P1 init.buildPlantHydraulicState(raw, built)   ✓ (2026-06-25, test 10/10)
     -> built.Y.plant_hydraulic  {psi_root, psi_stem, psi_leaf_sun, psi_leaf_sha, W_root/stem/leaf(dynamic)}
     -> built.p.plant_hydraulic  {PSIL_RANGE, VC, KS, gpmax}        ← 可重算缓存(不进 parameters)
     复用 VulnerabilityCurve/initWaterStorage/calculate_gpmax;组件级 gate-A(不跑 Init_Plant 脚本,measured 坏)
P3 init.buildPlantCNParams(raw, built)   ✓ (2026-06-25, test 22/22)
     -> built.parameters.plant_cn.{coupling, cn_ratio, allocation, respiration, turnover, storage}
     优先级 site.toml > PFT 表(getPFTValue,缺字段返默认)> 兜底;字段名核对 getPFTParameterSet。
     (stem_leaf 等 Apple 未定义→默认;只整理参数,不初始化池)
P2 init.buildInitialPlantPools(raw, built)   ✓ (2026-06-25, test 26/26)
     -> built.Y.plant_cn.{C, N, transfer, baseline}  C-primary 初值。读 raw.plant(plant_initial_pools.csv,精度同步到
        config 全精度)+ P3 的 plant_cn 参数;C→N=C/CN(叶 CN_L/根 CN_target_root/果 CN_target_fruit/活木·活粗根 CN_lw/
        死木·死粗根 CN_dw);stem live/dead 按 flivewd;粗根=croot_stem×stem;storage 按 f_storage_root 拆;
        N_storage 独立给(非 C/CN);xsmr 仅 C 侧(碳债,0,计入 C_pool_init_const);transfer(leaf/root/livestem)=0。
        gate-A 取径 A:真调 setInitialNitrogen(构 SP/VG/MeteoData/最小 state)只比第④块池 (1) 值(26 项)。

**D-plant 完成(P1/P3/P2)+ 已接入 buildState**(顺序 Y(soil)→P1→P3→P2):
`built.Y.plant_hydraulic`/`built.p.plant_hydraulic`、`built.parameters.plant_cn`、`built.Y.plant_cn`。
built.Y = 扁平土壤态 + 嵌套 plant_hydraulic/plant_cn(过渡,主循环重构再统一 Y.soil)。test_buildState 88。
```
- VC/KS 是参数派生的查表缓存 → `built.p.plant_hydraulic`(不放 parameters.veg.hydraulic)。
- prognostic 暂放 `built.Y.plant_hydraulic`/`built.Y.plant_cn`(土壤态仍扁平,不现在重构 Y.soil)。
- gate-A:P1 **不能跑 Init_Plant 脚本**(其 measured calculateRootDistribution 工作区坏),走组件级
  (复用 VulnerabilityCurve/initWaterStorage 复刻公式);P2 逐元素对 setInitialNitrogen 第④块池 (1) 值。

## D-soil 完成 (2026-06-25)

`init.buildNitrogenSoilParams(built)` → `built.parameters.soil_n`(单独容器,不塞进 soil 数组堆)。重排
setInitialNitrogen 第①块的**静态土壤 N 参数**:RHO_bulk、RHOKG(=KG·bulk 吸附缓冲)、DIF1111[1×NL×nS]
(NH4/NO3 排布)、DISP、K_mineralization、Vmax/Km_nit、Vmax/Km_denit、K_nh4_dissolve、Kurea、pH、MSOC_volume、
Theta_f(field capacity)。**不碰** N 浓度 state / Retard/Dsh/Smin 工作数组 / tS 诊断(留 allocateWorkspace)。
gate-A 取径 A:真调 setInitialNitrogen 只比第①块(bd_source='soil' 对齐 standalone SP),`test_buildNitrogenSoilParams` **15/15**。
已接入 buildState(parameters.soil 之后)。

## buildScope (S1 完成 2026-06-25;含消费点审计)

Explore 代理做了 SCOPE **消费点审计**(run_stemmus [12]-[16] + COMPASS.m 主循环 + ebal_clm45* + RTMo/RTMf/RTMz/RTMt + select_input + loadTimeSeries)。结论:
- **6 子组**:options / leaf_optical / soil_optical / geometry / spectral / canopy_rt。
- **静态 vs 动态**:buildScope 只冻静态;**运行时被重写的只存基线**——Cab/Cca(LUNA 每日覆写)、rss/rbs/rs_thermal(ebal 按 SMC 每步重算)、tts(每小时自算)、Vcmo(LUNA)、LAI_layer(veg)。
- **veg/scope 边界**:LAI/LAI_layer/hc 归 veg(已做);layer_weight(垂直廓线形状,静态)归 scope.canopy_rt,用 veg 的 LAI_profile_* 参数。
- **来源**:4 节在 site.toml(S1);spectral=硬编码 io.define_bands、canopy_rt=run_stemmus 硬编码常数(nlayers=30/litab/...)、optipar/soil 谱/atmosphere=[model_library] 外部文件 → 都 S2。

**S1 ✓**:`init.buildScope(raw.config)` → `parameters.scope.{options, leaf_optical, soil_optical, geometry}`。
4 节纯抽取 + friendly→legacy 名(N_mesophyll→N、leaf_rho/tau_thermal→rho/tau_thermal、fluor_*→kNPQs/qLs/stressfactor/beta、
cca_function_of_cab→Cca_function_of_Cab、observer_zenith→tto、relative_azimuth→psi)+ bool→0/1;Cab/rss/rs_thermal 等存基线。
gate-A `test_buildScope` **23/23**(passthrough 对 site.toml);已接入 buildState(test_buildState 90)。
- **S2 ✓ (2026-06-25)**:`buildScope(raw, built)` 加 `scope.{spectral, canopy_rt, optical_library}`。
  spectral=复用 `io.define_bands()`;canopy_rt=run_stemmus 硬编码常数(nlayers=30/nlincl=13/nlazi=36/litab/lazitab/x/xl)
  + 静态 layer_weight/shape/z(`vegetation.calc_lai_profile` 用 reference LAI=LAI_constant,LAI-无关,运行时 LAI_layer=layer_weight×LAI);
  optical_library=[model_library] 文件路径解析(data/model_library/)。**LAI/LAI_layer 当前值不进 S2**(留 veg/主循环)。
  gate-A `test_buildScope` **32/32**(S1 23 + S2 9:spectral==define_bands、canopy_rt 常数、sum(layer_weight)==1、layer_weight==calc_lai_profile、路径构造)。已接入 buildState(scope 在 veg 之后,test_buildState 90)。
  > **遗留(数据迁移)**:`data/model_library/` 目录还不存在——optipar 等实际在 `data/input/fluspect_parameters/` 等处,
  > 且 site.toml 写的 `FLEX-S3_std.atm` 现仅有 `FLEX-S3_Dry.atm`。optical_library 现是**意向路径 spec**(未验存在),
  > 待建 model_library 目录 + 迁移/对名后再验存在或加载。spectral/canopy_rt 不受影响。
- **S3 待做**:legacy ScopeParameters/options/canopy/leafbio adapter(若老主循环仍需老 structs)。

## 尚未替换的 run_stemmus_initialize 部分

- `buildVeg`: cfg.Veg 等价(物候/LAI/冠层/根系)。
- `buildNitrogenSoilParams` (D-soil): 很薄，主要固化 `RHOKG` 和 N 求解器需要的排布。
- `allocateWorkspace`: Solver/State/Flux/OLD/cache 分配。
- `buildForcing`: 替换 `MeteoData`。
- `buildVeg`: 物候、LAI、冠层/根系相关初始化。
- `buildInitialPlantPools`: 植物 C/N pool 初始化。
- `buildScope`: SCOPE 参数、光谱库、冠层结构、波段、方向性设置。
- `initialize`: 顶层编排，最终替换老初始化入口。

## 建议下一步

`[7b]`/`reconcileBulkDensity` 已完成。下一步顺序：

```text
1. buildConfig          §8 名字映射 + 开关联动 + 非网格派生量(喂老主循环必需)
2. buildNitrogenSoilParams (D-soil)   固化 RHOKG + N 求解器排布
3. buildForcing         替换 MeteoData
4. buildVeg             物候/LAI/冠层
5. buildInitialPlantPools (D-plant)   植物 C/N 池(需 veg 的 CN)
6. buildScope           SCOPE 参数/光谱/冠层/波段
7. allocateWorkspace    Solver/State/Flux/OLD/cache 分配
8. initialize           顶层编排, 替换老入口
```

`buildConfig` 排最前，因为 forcing/veg/scope/E 都要吃统一后的配置名与开关联动结果。

## 链条完成 (2026-06-25)

`loadSite → checkSite → buildState → allocateWorkspace` 全部落地并验证，`initialize` 门面收口。

- **allocateWorkspace** ✓ — 单写 `built.workspace`(step.trial/previous · derived · solver.{coeff,water,thermal,gas,nitrogen,SAVE} · flux) + `built.diagnostics`(6 个松散遗留 buffer: plantN_diag/nsc_diag/NH4_total/NO3_total/CHK2/CHK3)，**不碰** `built.Y/parameters/p`。关键决策: legacy `h`≡`Y.h`(accepted)，trial 暖启=Y，`step.old` 不单列；granule/冰相 N 池补进 Y 单写者 `buildInitialState`(zeros)。审计+决策全档: `docs/allocate_workspace_audit.md`。gate-A `test_allocateWorkspace` 52/52。
- **initialize** ✓ — `model = init.initialize(siteDir)`，返回单个 model 容器(siteDir/raw/check_report/built/status/created_at)，不返回裸 built、不 `assignin base`、不吞 checkSite 错误。`test_initialize` 21/21。
- 全套 **599 检查 / 20 文件 / 0 失败**，checkcode 干净。

S3 legacy adapter 跳过(长期重构不做旧入口短期桥)。剩余: 新主循环 `runOneStep(model, kt)` 消费 `model.built`(workspace/Y 双状态推进)，以及把老 `COMPASS.m` 逐步迁到新容器。

## B 阶段 runOneStep 迁移台账 (2026-06-28)

逐段核对老 `COMPASS.m` pre-loop 段（[2]-[6] + checkpoint），按 §1 state/aux/param 归位。
**init 侧能定的已落地；运行时演化态留 runOneStep。**

### 已落地（init 侧，runOneStep 直接读）

| 老 COMPASS.m | 新落点 | 测试 |
|---|---|---|
| [6] `debug_substep_*` 策略(growth/shrink/grow_after_success/hold_after_reject/grow_iter_threshold/enabled) | `control.numerics.substep`(buildNumerics，可选 config 覆盖，默认=老硬编码值) | test_numericsSubstep 9/9 |
| [6] `debug_substep_min_dt` | `control.time.Delt_t_min`(不在 substep 重复) | — |
| [3] setBoundaryCondition DSTOR/DSTOR0/RS | `Y.surface_water`(prognostic，=0) | test_surfaceWaterInit 9/9 |
| [3] DSTMAX/NBChh/IRPT1/IRPT2 | `workspace.surface`(scratch，NBChh=NBCh_top≠1?2:1) | 同上 |
| [3] fmax/Ks0 | `parameters.soil.fmax(=0)/ks0(=ks(end))` | 同上 |
| [3] BtmPg/UseHydrostaticAirBottom | `boundary`(buildBoundary) | 同上 |
| [2] Grid/Time/SoilConstants/Groundwater | domain / control.time / soil.constants / runtime.groundwater | (前期) |
| [2b] Registry | runtime.diagnostics(Registry) | (前期) |
| [2b] checkpoint/diagnostic/partial interval + output 目录(COMPASS.m:148-171) | runtime.output(output_dir/checkpoint_dir/**nc_dir/csv_dir/legacy_csv** 比老更细) + runtime.checkpoint(checkpoint_interval/diagnostic_interval/partial_output_interval，默认 24/360 与老一致；env STEMMUS_OUTPUT/PARTIAL_OUTPUT_INTERVAL 都接) | 逐字段核对 ✓ 2026-06-28 |
| [4]/[4b]/[5] plantHydro/W储/LAI 日状态 | Y.plant_hydraulic / p.plant_hydraulic / vegetation_state | (前期) |
| [6] TransportCoefficient | workspace.derived(结构替代，容器名没了内容在) | — |

### 留给 runOneStep（运行时演化态 / 主循环逻辑，B 阶段建）

- **变步长运行时态**：`Delt_t`(当前步长)、`Delt_old`(上步长)、substep `retry_count`(重试计数) →
  runOneStep stepping state（消费 `control.numerics.substep` 策略 + `+solver` 的
  accumulate/restore/finalize/initHourlySubstep 执行函数）。
- **地表水文每步演化**：`Y.surface_water` + `workspace.surface` 的逐步更新（入渗/积水/产流/NBChh 动态切换）。
- **checkpoint 写盘**：`mkdir` **已归 ③ beginRun**（2026-07-01 落地 `ensureOutputDirs`）；prepareRuntime
  只记路径(pure allocator)。剩下的**写文件**——每步 checkpoint save / partial output 归 runOneStep 步末 hook、
  最终 NetCDF/CSV flush 归 finalizeRun。（早期"mkdir 归 runOneStep"作废，**以 beginRun 负责 mkdir 为准**。）
- **STEMMUS_MAX_KT** 短窗口 env(老 [2])：smoke 便利，未迁，runOneStep 可补。
- `debug_scope_output_csv`：SCOPE 调试 CSV 开关，runOneStep 局部。
- **主循环本体**：select_input→`prepareScopeStep` 喂 SCOPE；逐步 Picard(water/thermal/gas/nitrogen)
  + ebal + 植物 C/N + 物候/LAI 日更 + LUNA + VOD。
- **SCOPE 叶片/土壤光学谱(老 COMPASS.m [7])**：`fluspect`(fversion) 算 leafopt(rho/tau/reflZ/tranZ)
  + lidf/emis/rs_thermal/rt_thermal 分支。老版是"循环前**静态预计算一次**"(假设 leafbio 每步不变)，
  但 **Cab 现已每日从 LUNA 更新，静态前提失效**(注释自承"vary 就 move back into loop") →
  runOneStep 里应**每日(LUNA 触发时)重算 fluspect**，不再循环前算一次。
  `buildScope` 只冻**输入参数**(Cab/Cw/Cs/N/V2Z/rho_thermal + spectral 波段 IwlP/IwlT)，RT 计算归 runOneStep。
  red-edge SIF 诊断索引(re_wl 590-850 / IwlP_re / sif_wl / IwlF_in_re) → runOneStep + diagnostics(radiation_spectral)。
  **注意：SCOPE/VOD/正演原代码只读，要改先 copy 到独立工作区**(见 feedback_no_edit_vod_code)。
- **地下水耦合(老 COMPASS.m [9])** → **✅ 全部在 ③ beginRun 完成(2026-07-01)**：gw_Dep / indexBotmLayer /
  tempBotm=Y.T(idx) / **NBCh_bottom=1**(新字段) / 土柱底检查 + **Dunnian 饱和径流(forcing adapter)**——见下「③ beginRun 收口」段。changwu Coupling=0 不触发。
- **resume 续算加载(老 COMPASS.m [9])** → **✅ 已在 ③ beginRun 完成(2026-07-01)**：`applyResumeCheckpoint`
  return-value(替掉老 assignin/evalin)，从新 Y 重建 workspace + 恢复 runtime memory + 设 next_kt。见下「③ beginRun 收口」段。

### ③ beginRun 收口完成 (2026-07-01)

老 COMPASS.m 的 pre-loop 一次性副作用（mkdir / resume / 地下水）全迁进 `src/beginRun.m`（③ 层）；
`runCOMPASS` 的 `~stopAfterInit` 分支跑 `initialize → prepareRuntime → beginRun → for kt runOneStep → finalizeRun`。

- **3.2 输出目录 + 输入校验**：`ensureOutputDirs`（幂等 mkdir output_dir/checkpoints/nc/csv，失败报错）+
  `validateBeginRun`（runtime.output 齐全）。`test_beginRun` 7/7。
- **3.1 resume**（return-value，**无 assignin/evalin**）：`applyResumeCheckpoint(model, io.loadCompassCheckpoint(file))`——
  `built.Y=ckpt.Y` → `allocateWorkspace`(从新 Y warm-start trial) + `initializeSoilAux`(重算 h_frez/derived) →
  恢复 runtime prognostic memory → `run_control.next_kt=ckpt.kt+1`；`runCOMPASS` 从 `next_kt` 起跑。
  **checkpoint schema v1**（`io.saveCompassCheckpoint`）：存 `Y`（含 snow/surface/canopy_water 状态,
  2026-07-05 增量1-3 搬入 Y）+ runtime `step/substep/calendar` +
  **prognostic memory**(phenology / lai_state / luna / leaf_resp / groundwater / management)；**不存** workspace/aux
  (从 Y 重算)、Registry(续算重填)。**luna handle 用 `Climate10dBuffer.toStruct/fromStruct` 显式序列化**，
  不赌 MATLAB save(handle)。测试：climate10dBufferSnapshot 11 / compassCheckpoint 10 / resume 14。
- **3.3 groundwater 主体**（仅 Coupling==1；changwu no-op）：gw_Dep / indexBotmLayer / tempBotm=Y.T(idx) /
  **NBCh_bottom=1**(新字段，非旧 NBChB) / 土柱底检查(gw_Dep>tot_depth → error)。`test_groundwaterSetup` 9/9。
- **✅ Dunnian 饱和地表径流(forcing adapter，2026-07-01 完成)**：老 COMPASS.m:325-326 在循环前改
  `MeteoData.Precip_msr`(日 mm)，因老链日→时降尺度在循环内、改 raw 能传到；新链 `buildForcing` 已在 init
  就把 `raw.Precip` 展开成 `hourly.Precip`(cm/s)——改 raw 白改、改 hourly 又不匹配 `updateDunnianRunoff` 的
  日 mm 接口。**方案 = forcing adapter**（不在 beginRun 裸改 raw、不在 runOneStep 每步扣）：
  `forcing = init.applyGroundwaterForcingAdjustment(forcing, gwDep, built)` → `groundwater.updateDunnianRunoff(
  raw.Precip, gwDep)` 得 `[runoffDunnian, corrected]` → 存 `forcing.runoffDunnian`(诊断) → **若 precip 真变了**
  (gw_Dep≤1cm 饱和才全归零)才 `raw.Precip=corrected` + 重建 `init.buildHourlyForcing`→`init.buildDerivedForcing`
  →`init.buildAtmosphericGasDrivers`；gw_Dep>1cm 是纯 no-op(hourly 保持 bit 一致)。
  - **前置重构 = `buildForcing` 拆 public**：nested `buildHourly`→`init.buildHourlyForcing`(F2)、
    `buildDerived`→`init.buildDerivedForcing`(F3)(adapter 可重建)；共享 `getForcingConfig` 提到
    `+init/private/`；`buildRawLayer`/`buildSolar`/`getcol` 留 buildForcing。`gas` 本就 public。
    **golden A/B `isequaln` byte-identical + test_buildForcing 81/81** 证拆分零功能漂移(纯搬迁)。
  - **接线**：beginRun 3.3 在 `applyGroundwaterSetup` 之后调 adapter(消费它算出的 gw_Dep)。
  - **测试**：`test_dunnianForcingAdjustment` 16/16(非饱和 no-op / 饱和归零 + 只动 precip 不动热力学 /
    beginRun 浅水位循环前降水清零)；回归 groundwaterSetup 9 / beginRun 7 / runCOMPASS 24 / resume 14 全绿。
  - **后续(未做)**：若接动态 MODFLOW(每步水位变)再升级成 runOneStep 每步 forcing correction；当前静态
    gw_Dep，forcing adapter 更干净。changwu Coupling=0 不触发。

**runCOMPASS 五层现状**：①initialize ②prepareRuntime ③beginRun **全实（3.1 resume / 3.2 output-dir /
3.3 groundwater 主体 + Dunnian adapter 全收口）**；④runOneStep = substep 框架（物理占位 `attemptSolve`）；
⑤finalizeRun 待（等 ④ 有诊断数据再 flush）。

### 初始化迁移审计 (2026-07-01, 45-agent workflow 交叉验证)

对 `COMPASS.m` 25-343（[1]-[9] + resume）逐小节 vs 新链做了 22 单元 × (审计→对抗验证) + 覆盖批判。
结论：**16 已迁(干净/表面差异) · 5 部分迁移属设计内延后到 ④/⑤ · 1 真缺口(s4)**，无整块漏搬。

**s4 真缺口处理（Init_Plant）**：
- **✅ 4 个校验门已补**（`test_plantValidationGates` 9/9，changwu 干净通过）：
  `buildVeg:badPsilMin`(psil_min<0) / `buildVeg:invalidRootDistMode`({0,1}) /
  `buildVeg:badRootContact`(enabled 时阈值有限负 + start>end) /
  `buildPlantHydraulicState:badInitialWaterPotential`(init ψ 有限 ≤0)。均引 Init_Plant.m 出处。
- **🟢 低优先级留置**：`dyn_K_sat/dyn_psi_63/dyn_c_weibull`（Weibull 脆弱性曲线）buildVeg 未搬（其余 5 个
  dyn_* 已搬）。**只影响 off-by-default 的 V2 预报方案**（`prognosticPlantWater`，仅 demo/test/旁挂 wrapper
  调用，不在主链）；**live VOD 走 ψ_leaf→PV 曲线→RWC→含水量，不用 Weibull**（见 [[project_compass_vod_water_path]]）。
  等真激活 V2 再顺手补 3 行。

**两处"已迁但行为变了"（保留新行为 + 此处标注）**：
- **冻结门**：旧 `StartInit` 无条件算初始冻结头；新 `initializeFreezingHead` 仅 `lFreezing==1` 才算
  （关冻融+初始零下时分叉；changwu 开冻融→无差；新行为更对）。
- **Dunnian**：旧 `COMPASS.m:325-326` 只改日 `Precip_msr`，主循环读构造时定死的小时 `MeteoData.Precip`
  → **旧 Dunnian 是隐性 no-op**；新 adapter 重建 hourly→**真扣降水**（更对；仅耦合场景+饱和才触发）。

**次要差异（可接受，已记）**：[3] buildBoundary 丢 `ismember([1 2 3])` BC 范围校验；[2b] MAX_KT 冒烟 Registry
按全 tS 开 + `COMPASS_ONLINE_GS_DIAG` 前缀翻标志不注册组；[6] `write_periodic_partial` 码内默认 true→false
(site.toml 抵消)；resume `leaf_resp` 存了没恢复(当前空)；`COMPASS.m:147 debug_scope_output_csv` 未迁(纯 debug)。

### ④ runOneStep 迁移前必读（审计标出的延后项 + 接线雷）

> 主循环逐块迁移进度见 [compass_complete_architecture.md](../architecture/compass_complete_architecture.md) §13。

**设计内延后到 ④ 的每步项（接物理时一起迁）**：
- [5] LAI Method-B 每步状态：`wt_day_min` / `lai_day_prev_idx` / `nTmax`（`buildVegState` 已存 init 槽，
  每步累加器/日检测/forcing 钳位在 ④）。
- [6] `Ta_history` 数组 + 每步 `Tgrowth` 滚动均值（init 槽在 `prepareCalendar/prepareLeafResp`；
  legacy `luna_clim.push` **每步无条件**推，与 luna_enabled 无关——④ 里照抄别加门）。
- [6] `solveCtl`(EB 迭代) / `TransportCoefficient`（每步 solve scratch，`runOneStep` 已有 `gatherSolveControls`）。
- [7] SCOPE **RT 计算**（fluspect/lidf/emis/thermal 填/红边）——因 LUNA 每日改 Cab，静态假设失效，
  **每日重算**（输入已冻结在 `buildScope/scope_static`；`rss/rbs` 也每步重算，别冻）。
- [8] 每步输出数组（EVAP/Infiltration/QET/QEB/…→ `reg.write` Registry）+ `SWE`/`Tb`/`Tbh`/监测深度。
- [1] `ISFT` 门（`COMPASS.m:434`）。

**接线雷（④ wiring 时务必）**：
- **phen 年重置**：用 `Y.vegetation.phenology.year`(=0 哨兵；2026-07-05 增量#5 从 runtime.phenology 搬来) 做 `cur_year~=phen_year` 首日重置门，
  **不能**用 `runtime.calendar.current_year`(=year(start)，每步覆写) —— 否则丢"首日重置"保证。
- **steps_per_day**：legacy 变步长下 `COMPASS.m:1288` 每步重算 `round(86400/Delt_t)`；新链存 init 时单值。
  Tgrowth_window/LUNA 窗长在变外步下会漂——④ 若要严格一致需按当前 Delt_t 重算。
- **Dunnian**：新 adapter 已在 ③ 循环前一次性处理，④ **不要**再每步扣（除非将来接动态 MODFLOW 水位）。

### 命名债

老 `debug_substep_*` 前缀误导（非调试，是生产变步长控制器）。新 `control.numerics.substep.*`
已用干净名；runOneStep 接入时老前缀一并退役。

**env 前缀统一（2026-06-28 起）**：项目从 STEMMUS-SCOPE 演化来，环境变量沿用 `STEMMUS_*`（共 23 个），
新代码已零散用 `COMPASS_*`。统一机制 = `io.getEnv(suffix)`：**`COMPASS_<suffix>` 优先、`STEMMUS_<suffix>`
兜底**（老脚本不断），`test_getEnv` 5/5。
- **已统一**：`prepareRuntime` 的 OUTPUT / PARTIAL_OUTPUT_INTERVAL / WORKER_ROOT / RESUME_FILE /
  ONLINE_GS_DIAG / PS_MAX。
- **待统一**（随 runOneStep / 各文件迁移时逐个改用 `io.getEnv` 传 suffix）：DIAG_INTERVAL /
  DISABLE_ROOT_RETURN / END_TIME / LFREEZING / LFUSION_OVERRIDE / LSOILAIREFC / MAX_KT / NMSRMN /
  SOILPROPERTY_FILE / SPIKE_DIAG / START_TIME / THREEWAY_COMPARE / THREEWAY_DAYS / TOL_WATER /
  WATER_RESIDUAL_SNAPSHOT_{DIR,KT,MAX,STAGE}（散在 solver / 主循环 / IO / 诊断）。
- **退役条件**：所有读取点都走 `io.getEnv` 后、过渡期满，删 `STEMMUS_` 兜底分支。

**入口改名(2026-06-29)**：新链顶层入口 `COMPASS_main` → **`runCOMPASS`**(camelCase 动宾，与 buildState/runOneStep
同风格；COMPASS 保持大写专名；与 runOneStep 配套=跑整个 vs 跑一步)。老脚本 `COMPASS.m` 过渡期**并行、逐步废弃**；
等其退役再议是否让 runCOMPASS 接管 `COMPASS` 名(注意 Windows 文件名大小写不敏感，不能用 Compass/compass 变体)。
配套：`test_COMPASS_main`→`test_runCOMPASS`；error id `runCOMPASS:runOneStepPending`；引用更新
prepareScopeStep/test_scope_step_parity/run_load_site。docs/superpowers 历史快照未改(保留当时命名)。

## runCOMPASS 五层执行流 (2026-06-29)

老 COMPASS.m 的 pre-loop [1]-[9] + 主循环，在新链拆成 5 层，`runCOMPASS` 串起来：

```text
runCOMPASS(siteDir):
  ① init.initialize(siteDir)      构建状态  buildState/allocateWorkspace/initializeSoilAux
  ② init.prepareRuntime(model)    pre-loop 服务+状态槽（pure allocator：不碰 Y / 不碰盘）
  ③ beginRun(model)               进循环前·做一次（有副作用）：resume 加载 / mkdir / GW 预设
  ④ for kt: runOneStep(model,kt)  每步物理推进
  ⑤ finalizeRun(model)            跑完·做一次：诊断 flush(NetCDF/CSV)
```

**为什么 ③ 不并进 ②**：② 是纯分配器（doc 契约：no prognostic write / no file I/O，可测可重入）；
③ 三件事都有副作用——resume **覆盖 Y/KT/TIME**、mkdir **写盘**、GW 改 forcing+控制流——并进 ②
会毁掉它的纯度（测试要 mock 磁盘、不能安全重调、resume 时序错乱）。所以 ③ 独立。

**resume 只在 ③ 做一次**（不在每个 runOneStep 反复做）：进循环前加载恢复点 → 从恢复后的 KT 续跑；
checkpoint **写入**才是 ④ 每步末尾的 hook，⑤ 收尾再 flush 全量诊断。

**进度**：①② 完成；③ `beginRun` / ⑤ `finalizeRun` = 骨架壳（noop + TODO，2026-06-29 建，接进
runCOMPASS 的 ~stopAfterInit 分支）；④ `runOneStep` 下一步。③④⑤ 的 TODO 明细见上面「留给 runOneStep」清单。
