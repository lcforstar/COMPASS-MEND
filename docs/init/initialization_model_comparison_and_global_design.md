# COMPASS 初始化完成口径与区域/全球化设计稿

本文档给后续 Claude/Codex 接力使用。目标不是泛泛比较陆面模型，而是回答一个具体问题：

> COMPASS 的初始化怎样才算完全做完，并且以后能从单点站点自然扩展到区域/全球网格，同时和中间求解模块分开？

## 1. 当前结论

按新边界，`init.initialize(siteDir)` 作为“站点输入 -> 可运行状态容器”的初始化门面已经打通：

```matlab
model = init.initialize(siteDir)
  raw    = io.loadSite(siteDir)
  report = io.checkSite(raw)
  built  = init.buildState(raw)
  built  = init.allocateWorkspace(built)
```

当前回归验证：22 个初始化/runtime/区域合同测试文件，695 checks，0 failures。

如果按老 `src/COMPASS.m` 的口径，也就是“从 `run_stemmus_initialize` 到进入
`KT` 主循环前的全部准备”，现在已经拆成两段：

- `init.initialize` / `init.initializeRegional`：输入、参数、初值、forcing、workspace。
- `init.prepareRuntime`：diagnostics registry、checkpoint/resume policy、groundwater
  预处理、LUNA buffer、phenology state、Picard diagnostics、
  hysteresis/dual-porosity 状态记录、SCOPE static library 记录。
  （canopy water 不再是 runtime 容器：state→`Y.canopy_water`、history→`built.diagnostics.canopy_water`、config→`params.veg.canopy`，2026-07-05 增量#3。）

因此推荐边界：

```text
init.initialize      = 输入、参数、初值、forcing、workspace
init.prepareRuntime  = 进入主循环前的运行服务
runOneStep           = 真正推进一个时间步
```

不要把这些 runtime 服务塞回 `buildState`，否则初始化又会变成老 `COMPASS.m` 那种混杂入口。

## 2. 外部模型怎么做

### 2.1 ClimaLand

ClimaLand 是最值得学的初始化结构。

它的核心流程是：

```julia
domain -> model(parameters, forcing, components)
Y, p, coords = initialize(model)
set_ic!(Y, p, t0, model)
set_initial_cache!(p, Y, t0)
ODEProblem(tendency!, Y, (t0, tf), p)
```

关键点：

- `domain` 先独立创建，可以是单点 `Column`，也可以是全球/区域 `SphericalShell`。
- `Y` 是 prognostic state，只放会随时间推进的量。
- `p` 是 auxiliary/cache/driver，不是主状态，但会在每步更新。
- `set_ic!` 只设初值，不承担读文件、输出、checkpoint。
- Simulation 层才接 SciML 积分器、driver update callback、diagnostics callback。

对 COMPASS 的启发：

- 你现在的 `built.Y` 对应 ClimaLand 的 `Y`。
- 你现在的 `built.p` / `workspace` 对应 ClimaLand 的 `p/cache`。
- `init.initialize` 不应该管求解器推进；`runOneStep` 或未来 integrator 才管推进。
- 区域/全球能力应该从 `domain` 和数组形状开始，而不是从主循环里硬改。

### 2.2 Noah-MP

Noah-MP 是传统 Fortran 陆面模型路线。它很强在区域/全球格点运行。

典型结构：

```fortran
read static NetCDF fields: lon_lat, landmask, veg_soil, tbot, VG parameters ...
CALL NOAHMP_INIT(...)
time loop:
  read forcing
  per grid cell:
    CALL NOAHMP_SFLX(...)
```

关键点：

- 静态地理场和参数场来自 NetCDF。
- 初始化一次性把很多格点数组分配/填好。
- 单步核 `NOAHMP_SFLX` 是大接口，输入/输出状态都在数组里。

对 COMPASS 的启发：

- 区域/全球必须支持共享静态 NetCDF 数据，不应每个格点一个站点文件夹。
- `domain` 里要有 mask、lat/lon、cell id。
- forcing 要能是 `time x cell` 或 `cell x time` 的数组，而不是单点 Excel/CSV。
- 求解核应消费已经准备好的数组，不应该自己读文件。

### 2.3 DSSAT/CSM

DSSAT/CSM 是阶段机路线：

```fortran
RUNINIT  -> run-level initialization
SEASINIT -> season initialization
RATE     -> rate calculation
INTEGR   -> integration
OUTPUT   -> output
SEASEND  -> season end
```

关键点：

- 初始化不是一个函数，而是按 `DYNAMIC` 阶段分散到各模块。
- 优点是农业管理、季节、多作物试验流程清楚。
- 缺点是状态流不如 ClimaLand 那种 `Y/p/cache` 明确。

对 COMPASS 的启发：

- 植被物候、管理事件、年度重置以后可以借鉴阶段机。
- 但当前初始化重构不建议回到 `DYNAMIC` 全局阶段变量，否则会降低可测试性。

### 2.4 DAISY

DAISY 是对象式农业土壤-植物模型。它有 `Daisy::initialize` 和 `Column::initialize`，并且 `Column` 持有 soil、vegetation、geometry、weather、management 操作。

对 COMPASS 的启发：

- `Column` 概念清楚：一个土壤-植被-气象耦合单元。
- 区域/全球可以理解成很多 column 的集合，但不能把每个 column 都变成一个独立文件夹和独立脚本。

### 2.5 HYDRUS

本地 HYDRUS 资料主要是手册/PDF，适合作为土壤水分和溶质数值精度参照，不适合作为 COMPASS 全耦合初始化架构的主要参照。

对 COMPASS 的启发：

- Richards 方程、滞回、双孔、溶质边界这些局部数值逻辑应保持严谨。
- 但 HYDRUS 的主要优势是单点土壤剖面，不是区域/全球生态-辐射-植物水力初始化。

## 3. 和 COMPASS 当前初始化的优劣对比

| 维度 | COMPASS 当前新初始化 | ClimaLand | Noah-MP | DSSAT/DAISY/HYDRUS |
|---|---|---|---|---|
| 数值稳定 | `Y/trial/previous/workspace/runtime` 已分开，比老 `COMPASS.m` 稳；求解器迁移仍待接 | 最强，`Y/p/cache/ODEProblem` 边界清楚 | Fortran 生产型稳定，格点数组成熟 | HYDRUS 土壤水盐很稳；DSSAT/DAISY 阶段机稳但耦合粒度较粗 |
| 使用难度 | 单点站点包清楚；区域/全球已有 NetCDF static 入口和 land-mask 打包合同，但真实产品映射仍待配置 | 抽象重，学习难 | 输入表和 NetCDF 多，Fortran 接口大 | DSSAT 试验管理友好，HYDRUS 单点 GUI 友好 |
| 运行速度 | MATLAB 单点可接受；区域/全球若不数组化会慢 | Julia/ClimaCore 更适合大规模并行 | Fortran 快，格点批处理成熟 | 多为单点/日尺度优势，不适合直接承载 COMPASS 全耦合 |
| 物理精确 | 对长武最细：水热气氮、植物水力、SCOPE、C/N 都进入初始化 | 通用陆面耦合正规 | 标准陆面能水过程成熟 | HYDRUS 土壤最细；DSSAT/DAISY 农业 C/N 管理强 |

判断：

- COMPASS 当前科学过程细节强，尤其适合长武果园单点。
- 当前主要短板不是物理项少，而是初始化/runtime/solver 边界还没完全像 ClimaLand 那样干净。
- 若要区域/全球化，必须先把 `domain`、`Y`、`forcing`、`workspace` 的数组契约定好。

更具体地说，COMPASS 新初始化和这些成熟模型相比，优劣不在同一个层面：

### 3.1 数值稳定

优势：

- 已经把 accepted state、trial state、previous state 分开：
  `built.Y` 是已接受状态，`workspace.step.trial` 是本步试算，`workspace.step.previous`
  是回退依据。这比老 `h/hh/hOLD/backup` 散在主循环里稳定。
- `parameters`、`boundary`、`forcing`、`Y`、`workspace` 分家后，初始化不会在求解中途被隐式重建。
- 单点和区域路径共用同一套初始化合同，区域只是多了 `nCell` 维度，不是另写一套物理。

短板：

- 还没有像 ClimaLand/SciML 那样把残差函数、状态向量和积分器完全统一起来。当前只是为手写
  Picard/Newton 求解准备好了状态边界。
- `runtime.diagnostics.registry` 仍是进入主循环前的运行服务，不是严格的数值状态。
- 如果 NetCDF 只覆盖某些派生土壤字段，例如只改 `theta_s` 不同步更新相关热参/特征含水量，
  物理一致性要靠上游产品映射层保证。

结论：初始化层的稳定性已经明显强于老 COMPASS，但还没有达到 ClimaLand 那种“求解器也被状态合同约束”的程度。

### 3.2 使用难度

优势：

- 单点用户仍然只需要站点文件夹：
  `model = init.initialize("data/sites/changwu")`。
- 区域/全球用户可以走共享 NetCDF：
  `model = init.initializeRegional(siteDir, staticNc)`，不用给每个格点建站点目录。
- `io.loadSite -> checkSite -> buildState -> allocateWorkspace -> prepareRuntime` 的顺序清楚，
  比老 `run_stemmus_initialize` 更容易定位输入错误。

短板：

- 配置项多，因为 COMPASS 同时初始化水、热、气、氮、植物水力、C/N、SCOPE，比 Noah-MP 或 HYDRUS
  单一过程入口更重。
- 区域/全球真实产品还缺一个映射层：变量名、单位、缺测、PFT/soil class 转参数表都需要配置化。
- 对新用户来说，`Y/workspace/runtime/p` 的区别需要先讲清楚；DSSAT 那种阶段机对农业用户更直观。

结论：单点使用难度下降了；区域/全球接口已经成形，但产品映射层没做完前，仍不能说“开箱即跑全球”。

### 3.3 运行速度

优势：

- 单点初始化已经是一次性构造，不再依赖 base workspace 和反复脚本副作用。
- 区域初始化采用 land-only columns，只给陆地格点分配数组；`active_index` / `max_cells`
  支持分块初始化。
- `scatterRegionalColumns` 提供从 land columns 回二维/全球网格的输出路线，避免求解器内部带海洋格点空跑。

短板：

- MATLAB 对全球大规模耦合求解不如 Noah-MP Fortran 或 ClimaLand Julia/ClimaCore 天然适合。
- 当前 NetCDF reader 支持 block 合同，但 still 是先读全变量再选择 active cells；超大产品后续要改成真正 chunked `ncread`。
- 主循环 kernel 还没数组化；如果每个 cell 用 MATLAB for-loop 跑完整水热气氮耦合，区域速度会是主要瓶颈。

结论：初始化已经按区域/全球速度方向设计了，但真正性能胜负取决于下一步 `runOneStep` 是否数组化/分块化。

### 3.4 物理精确

优势：

- COMPASS 初始化物理内容比多数通用陆面初始化更细：土壤 VG/热参/N 参数、初始水热气氮剖面、
  植物水力、植物 C/N 初池、物候驱动、SCOPE 静态光谱/冠层 RT 都进了初始化合同。
- 对长武苹果园这种单点研究，实测根系、树结构、冠层/水力参数能被保留下来，精度潜力高于只按 PFT 表的通用模式。
- HYDRUS 只强在土壤水盐，DSSAT/DAISY 强在农业管理；COMPASS 的优势是把 soil-plant-radiation-C/N 放进同一状态入口。

短板：

- 越精细越依赖输入质量。区域/全球时，如果土壤、根系、树结构、叶光学都只能用粗产品或 PFT 默认值，
  物理精确度会退化为“结构精细但参数粗”。
- `scope.optical_library`、真实全球 soil/veg/product 映射、动态 LAI/水力胁迫主循环还要继续接，初始化只是把入口准备好。
- 对一般全球陆面能水通量，Noah-MP 这类模式有更成熟的参数标定和大尺度产品经验；COMPASS 需要更多区域验证。

结论：COMPASS 的物理精细度是优势，也是负担。单点/小区域会很强；全球运行必须降低参数自由度并建立产品映射和验证体系。

### 3.5 总体取舍

当前初始化路线应继续坚持：

- 学 ClimaLand 的状态边界：`domain/Y/p/workspace/runtime/diagnostics` 分清。
- 学 Noah-MP 的区域数据入口：共享 NetCDF + land mask + packed columns。
- 借 DSSAT/DAISY 的农业阶段思想，但不要把全局 `DYNAMIC` 式副作用带回主循环。
- 借 HYDRUS 的土壤数值严谨性，但不要让 HYDRUS 式单剖面入口限制区域化。

所以，COMPASS 初始化的定位不是“最简单”或“最快”，而是：

```text
单点物理细节强 + 区域/全球数组合同已具备 + 求解器尚待迁移/数组化
```

## 4. 初始化完全完成的判据

“初始化完成”不能只看 `init.initialize` 测试绿。应分两层验收。

### 4.1 Stage A: 输入状态初始化

已经基本完成：

- 站点文件读取和 schema 检查
- grid/domain
- control/time/switches/numerics
- boundary/sources
- soil/soil_n
- forcing
- veg/vegetation_state
- plant hydraulic
- plant C/N pools
- scope static config/library
- `Y`
- `workspace`
- loose diagnostics buffer

证据：当前 22 个初始化/runtime/区域合同测试文件 695/695 通过。

### 4.2 Stage B: runtime 初始化

已新增：

```matlab
model = init.prepareRuntime(model)
```

当前已输出：

```matlab
model.runtime.diagnostics.registry
model.runtime.checkpoint
model.runtime.resume
model.runtime.groundwater
model.runtime.picard
model.runtime.soil_modes
model.runtime.scope    (SCOPE runtime workspace; the static optics cache is built.p.scope_static, built at init)
```

对应老 `COMPASS.m` 迁移项：

| 老代码内容 | 新家 |
|---|---|
| `diagnostics.Registry(tS,NL,NN)` + `registerAllDiagnostics` | `runtime.diagnostics` |
| checkpoint interval / partial output / resume env | `runtime.checkpoint` + `runtime.resume` |
| `io.readGroundwaterSettings` 和地下水初始计算 | `runtime.groundwater` |
| LUNA `state`/`last_*` + `leaf_resp.Ta_history` + `core.Climate10dBuffer`(handle) | `Y.plant_acclimation.{luna_state, leaf_resp.Ta_history, luna_climate_buffer}`(2026-07-05~07 增量#6a+6b);config -> `built.p.plant_acclimation`;buffer 存序列化 struct(Option A,`runtime.luna` 整删) |
| `vegetation.integrated.phenology.PhenologyState` + lai_state | `Y.vegetation.{phenology,lai_state}`(2026-07-05 增量#5) |
| `state_h2ocan*`(state) | `Y.canopy_water`(2026-07-05 增量#3) |
| canopy `*_history/residual`(诊断) | `built.diagnostics.canopy_water`(row-write) |
| `Ps_max/Ps_min_diurnal`(config) | `built.parameters.veg.canopy` |
| `solver.initPicardDiagnostics` | `runtime.picard` |
| `hysterIni` / `dualPorIni` | `runtime.soil_modes` |
| SCOPE 光学静态预计算 | `built.p.scope_static`(init cache;2026-07-05 增量#4,原 `runtime.scope_static`) |

当前验证：`test_prepareRuntime` 47/47，全套初始化 + runtime + 区域合同回归 22 个测试文件、695 checks、0 failures。

注意：`prepareRuntime` 目前是 runtime 服务容器落地，不是老主循环 kernel adapter 完成。尤其：

- `runtime.groundwater` 当前按老默认值建立 disabled 容器，未接通地下水开启时的动态预处理。
- `runtime.soil_modes.dual_porosity.status` 在开关打开时标为 `pending_adapter`，未调用旧 `dualPorIni`。
- `built.p.scope_static`（2026-07-05 增量#4 从 `runtime.scope_static` 搬来的 init 静态缓存）记录 spectral/canopy_rt/optical_library，并已做 `precomputeScopeOptics`（lidf_static/emis_static/leafopt_s/rs_static）光学预计算。
- `runtime.resume` 只验证 resume 文件请求存在性，未真正恢复工作区变量；真正恢复应放在 run facade 或 runtime resume adapter。

所以现在可以说“运行服务初始化骨架和主要容器已经完成”，但不能说“老 `COMPASS.m` 进入 KT 前所有副作用已经逐行等价迁移”。

## 5. 区域/全球化初始化设计

### 5.1 domain 先升级

当前单点：

```matlab
built.domain.kind = "column"
built.domain.site.latitude
built.domain.site.longitude
built.domain.NL
built.domain.NN
```

目标区域/全球：

```matlab
built.domain.kind = "column" | "regional" | "global"
built.domain.space.lat      % [nCell x 1] 或 [nY x nX]
built.domain.space.lon
built.domain.space.mask
built.domain.space.cell_id
built.domain.vertical.dz
built.domain.vertical.z_node
built.domain.NL
built.domain.NN
built.domain.nCell
```

设计原则：

- 垂向网格可先全域共享，后续再支持 per-cell soil depth。
- 区域/全球不是多个 `siteDir` 循环，而是一个 `domain` 管多个 cell。
- `site.toml` 只适合单点；区域/全球应有 `grid.toml + static.nc + forcing.nc`。

### 5.2 Y 的数组形状

当前单点：

```matlab
Y.h      % [NN x 1]
Y.T      % [NN x 1]
Y.theta  % [NN x 1]
Y.Theta_L % [NL x nD]
```

区域/全球建议：

```matlab
Y.h        % [NN x nCell]
Y.T        % [NN x nCell]
Y.p_gas    % [NN x nCell]
Y.nh4      % [NN x nCell]
Y.no3      % [NN x nCell]
Y.Theta_L  % [NL x nD x nCell]
```

MATLAB 推荐 `NN x nCell`，因为同一个变量的垂向列是连续逻辑，后续可以按 cell 批量或分块。

不要让 solver 读文件；solver 只读：

```matlab
Y, parameters, forcing_t, boundary_t, workspace
```

### 5.3 parameters 的区域/全球化

单点参数可以是标量或 `[NL x 1]`。

区域/全球参数应允许：

```matlab
parameters.soil.theta_s  % [NL x nCell]
parameters.veg.PFT_type  % [nCell x 1]
parameters.scope.options % scalar/global config
```

规则：

- 标量配置可以广播到所有 cell。
- 土壤、植被、地形、经纬度、海陆 mask 应来自 NetCDF 静态场。
- `data/model_library_nc` 是共享静态库，不放到站点文件夹。

### 5.4 forcing 的区域/全球化

当前单点 forcing 来自站点 CSV/Excel，已转换为小时序列。

区域/全球目标：

```matlab
built.forcing.time_utc
built.forcing.Ta      % [nTime x nCell]
built.forcing.RH
built.forcing.Precip
built.forcing.Rin
built.forcing.Pg
built.forcing.Ca
built.forcing.solar.tts % [nTime x nCell]
```

原则：

- 内部时间一律 UTC。
- 太阳几何用 UTC + lat/lon 算，不用 `round(lon/15)`。
- 全球日尺度 Tmax/Tmin 重建小时温度风险大，优先使用小时 forcing。

### 5.5 workspace 的区域/全球化

当前：

```matlab
workspace.step.trial.h % [NN x 1]
workspace.derived.Se   % [NL x nD]
```

目标：

```matlab
workspace.step.trial.h % [NN x nCell]
workspace.derived.Se   % [NL x nD x nCell]
workspace.solver.water.C7 % [NN x nCell]
```

原则：

- workspace 形状必须镜像 `Y` 和 per-cell 参数。
- 可以先实现单点和 tiny 2-cell 合同测试，不急着全局跑完。
- 后续运行速度优化再做 block/chunk，比如每次跑 100-1000 个 cell。

当前状态：

- `init.allocateWorkspace` 已支持 `built.domain.nCell > 1`。
- 单点 `nCell=1` 仍保持旧形状：node `[NN x 1]`，element `[NL x nD]`。
- 区域/全球合同形状已通过两条测试路径：
  `space struct` 的 tiny 2-cell mock，以及二维 NetCDF grid + land-mask 打包。

```matlab
Y.h                         % [NN x nCell]
Y.Theta_L                   % [NL x nD x nCell]
parameters.soil.theta_s     % [NL x nCell]
forcing.hourly.Ta           % [nTime x nCell]
workspace.step.trial.h      % [NN x nCell]
workspace.derived.Se        % [NL x nD x nCell]
workspace.solver.water.C7   % [NN x nCell]
workspace.flux.SinkS        % [NL x nD x nS x nCell]
diagnostics.NH4_total       % [NN x tS x nCell]
built.diagnostics.canopy_water.*  % [tS x nCell] (history; state = Y.canopy_water.state_* [1 x nCell])
scatterRegionalColumns      % land columns -> original grid_shape
```

公开入口：

```matlab
space = struct("lat", lat(:), "lon", lon(:), "mask", mask(:), "cell_id", cell_id(:));
model = init.initializeRegional("data/sites/changwu", space);

model = init.initializeRegional("data/sites/changwu", "path/to/regional_static.nc");

model = init.initializeRegional("data/sites/changwu", "path/to/regional_static.nc", ...
                                struct("active_index", landBlock));
```

当前入口仍是过渡合同，但已经不再只是口头 mock：

- `space struct` 路径：把已校验的单点站点包广播到多个 land cell。
- `NetCDF static` 路径：`init.readRegionalStatic` 读取 `lat/lon/mask/cell_id`，
  并按命名约定读取覆盖字段；`init.applyRegionalStatic` 把这些字段覆盖到
  `Y/parameters/forcing` 后再分配 workspace。
- NetCDF 的 `lat/lon/mask/cell_id` 可以是二维网格。当前实现默认按 `mask==true`
  打包成 land-only columns，`domain.space.grid_shape` 和
  `domain.space.active_index` 会保留原网格形状和 land cell 位置，供后续输出层
  scatter 回二维/全球网格。
- `init.scatterRegionalColumns(columnData, domain, fillValue)` 已提供反向映射：
  把 land-only columns 放回 `grid_shape`，非陆地格点填 `fillValue`。
- NetCDF 路径支持 `options.active_index` 或 `options.max_cells` 做 block/chunk
  初始化；这允许全球静态文件按 land-cell 子集分批进入 `initializeRegional`，
  避免一次性给全部陆地格点分配求解数组。

当前 NetCDF 命名约定：

```text
soil_<field>             -> built.parameters.soil.<field>
soil_n_<field>           -> built.parameters.soil_n.<field>
veg_root_RAI             -> built.parameters.veg.root.RAI
veg_PFT_type             -> built.parameters.veg.structure.PFT_type
Y_<field>                -> built.Y.<field>
forcing_hourly_<field>   -> built.forcing.hourly.<field>
```

这证明初始化层已经能接真实 `.nc` 网格/静态输入。但它仍不是完整全球产品读取器：
真实全球土壤/PFT/forcing 数据集的变量映射、单位换算、缺测填补、分块读取还需要下一阶段实现。

验证：`test_regionalInitializationContract` 49/49。

## 6. 推荐下一步给 Claude

### Step 1: `init.prepareRuntime` 已写

已完成 runtime 容器，不动主循环求解。

当前调用：

```matlab
model = init.initialize("data/sites/changwu");
model = init.prepareRuntime(model);
```

应存在：

```matlab
model.runtime.diagnostics
model.runtime.checkpoint
model.runtime.groundwater
model.runtime.picard
model.runtime.soil_modes
model.runtime.scope    (SCOPE runtime workspace; the static optics cache is built.p.scope_static, built at init)
```

测试：

```matlab
test_prepareRuntime
```

检查：

- 不改 `built.Y`
- 不改 `built.parameters`
- registry 尺寸等于 `tS, NL, NN`
- canopy water 初值(`Y.canopy_water.state_*`)和 history 长度(`built.diagnostics.canopy_water`)正确
- LUNA 开关关闭时不建重对象，开启时建 buffer
- dual/hysteresis 关闭时是 no-op，开启时有明确错误或完成初始化

当前结果：47/47 通过；完整初始化链 + runtime + 区域初始化合同为 695/695 通过。

### Step 2: 区域/全球初始化合同测试已加

已完成 tiny 2-cell mock，不要求真的跑全球，只证明合同成立。

新增公开入口 `init.initializeRegional(siteDir, regionalInput)`，当前用 tiny 2-cell mock
和临时 NetCDF 静态文件两条路径验证：

```matlab
built.domain.kind = "regional"
built.domain.nCell = 2
Y.h = zeros(NN, 2)
parameters.soil.theta_s = zeros(NL, 2)
forcing.Ta = zeros(nTime, 2)
workspace.step.trial.h = zeros(NN, 2)
```

测试名：

```matlab
test_regionalInitializationContract
```

当前结果：49/49 通过。

### Step 3: 再接 `runOneStep`

`runOneStep` 不读文件，不初始化参数，只做：

```matlab
[model, step_report] = runOneStep(model, kt)
```

推进规则：

```text
accepted state: model.built.Y
trial state:    model.built.workspace.step.trial
previous state: model.built.workspace.step.previous

step start: previous <- Y, trial <- Y
solve:      write trial/workspace
accept:     Y <- trial
reject:     restore previous/backup
```

## 7. 今天任务的真实完成标准

今天如果要说“初始化部分完全做完”，至少要满足：

1. `init.initialize` / `init.initializeRegional` 当前完整回归 695/695 继续全绿。
2. `init.prepareRuntime` 落地，并覆盖老 `COMPASS.m` 进入 KT 前的 runtime 初始化项。
3. `test_prepareRuntime` 全绿。
4. 单点长武可以得到：

```matlab
model = init.initialize("data/sites/changwu");
model = init.prepareRuntime(model);
```

5. 区域/全球化不是口头设计，而是至少有 tiny regional contract test，证明 `domain/Y/parameters/forcing/workspace` 的数组维度能表达多 cell。当前 `test_regionalInitializationContract` 已满足这一条。

当前状态：第 1、2、3、4、5 条均已有当前测试证据。仍需谨慎表述的是：这证明“初始化容器和 workspace 合同可以区域/全球化”，不证明主循环求解器已经能区域/全球并行推进。
