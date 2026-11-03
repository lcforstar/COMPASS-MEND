# COMPASS Initialization Handoff For Claude

本文件是给 Claude 接力用的短版说明。长版设计和外部模型对比见：
`docs/initialization_model_comparison_and_global_design.md`。

## 1. 已学习的外部模型结论

最该学的是 ClimaLand：把 `domain`、prognostic state `Y`、aux/cache `p`、
driver/forcing、diagnostics、solver/integrator 分开。初始化只构造状态和缓存，
不推进时间。

Noah-MP 给区域/全球化的关键启发是：静态地理场和参数场必须来自共享 NetCDF，
不能每个格点一个站点目录；求解核只消费已经准备好的数组。

DSSAT/CSM、DAISY 的阶段机/Column 思路对农业物候和管理事件有参考价值，但
COMPASS 当前不应回到全局 `DYNAMIC` 阶段变量。HYDRUS 适合作为 Richards 方程、
溶质和边界条件精度参考，不适合作为全耦合区域初始化架构。

## 2. 和外部模型相比的优劣重点

数值稳定：

- 新 COMPASS 已把 `Y`、`workspace.step.trial`、`workspace.step.previous` 分开，
  比老 `h/hh/hOLD` 散落在主循环里稳定。
- 仍弱于 ClimaLand 的完整 `Y/p/cache/integrator` 合同；下一步 `runOneStep`
  必须继续守住这个边界。

使用难度：

- 单点入口简单：`init.initialize(siteDir)`。
- 区域/全球入口已成形：`init.initializeRegional(siteDir, staticNc, options)`。
- 真实产品映射层还没完成，所以现在是“初始化合同可用”，不是“全球产品开箱即跑”。

运行速度：

- land-only columns、`active_index`、`max_cells` 已支持区域/全球分块初始化。
- 目前 NetCDF 覆盖仍是先读变量再取 active cells；超大产品后续要做真正 chunked 读取。
- 主循环 kernel 还没数组化，区域/全球速度最终取决于 `runOneStep` 迁移质量。

物理精确：

- COMPASS 初始化内容比通用陆面入口更细：水、热、气、氮、植物水力、植物 C/N、
  物候、SCOPE 静态 RT 都进了合同。
- 单点/小区域精度潜力强；全球运行会受土壤、PFT、根系、冠层、叶光学产品质量限制。

一句话：当前 COMPASS 初始化不是最简单、也不是最快，但物理入口最细；现在已经具备
区域/全球数组合同，下一步瓶颈转到 solver/runtime/output。

## 3. 当前初始化边界

当前边界已经固定为：

```matlab
model = init.initialize(siteDir)
model = init.initializeRegional(siteDir, regionalInput)
model = init.prepareRuntime(model)
```

为方便从主程序视角检查，当前新增了一个不碰老 `src/COMPASS.m` 的新门面：

```matlab
model = COMPASS_main("data/sites/changwu")
```

它只执行：

```matlab
model = init.initialize(siteDir);
model = init.prepareRuntime(model);
```

然后停在第一个 `KT` 前。`model.main_loop` 会记录：

```text
accepted state: model.built.Y
trial state:    model.built.workspace.step.trial
previous state: model.built.workspace.step.previous
next step:      runOneStep(model, kt)
```

如果设置 `stop_after_initialization=false`，当前会明确报
`COMPASS_main:runOneStepPending`，因为新 `runOneStep` 还没有迁移。

`init.initialize` / `init.initializeRegional` 负责：

- `raw = io.loadSite(siteDir)`
- `io.checkSite(raw)`
- `init.buildState(raw)`
- regional 扩展和 NetCDF static 覆盖
- `init.allocateWorkspace(built)`

`init.prepareRuntime` 负责进入 KT 主循环前的运行服务：

- diagnostics registry
- checkpoint/resume policy
- groundwater runtime container
- LUNA 10-day climate buffer
- integrated phenology runtime state
- canopy interception water runtime state
- Picard diagnostics
- soil mode status
- SCOPE static library record

主循环/求解器后续应只接：

```matlab
[model, step_report] = runOneStep(model, kt)
```

不要在 `runOneStep` 里读站点文件、重建参数或重新分配 initialization containers。

## 4. 当前容器结构

核心容器：

```text
model.raw
model.check_report
model.built.domain
model.built.control
model.built.boundary
model.built.sources
model.built.forcing
model.built.parameters
model.built.vegetation_state
model.built.Y
model.built.p
model.built.workspace
model.built.diagnostics
model.runtime
```

重要状态约定：

```text
accepted state: model.built.Y
trial state:    model.built.workspace.step.trial
previous state: model.built.workspace.step.previous
aux/cache:      model.built.p + model.built.workspace.derived
diagnostics:    model.built.diagnostics + model.runtime.diagnostics.registry
```

时间步推进建议：

```text
step start: previous <- Y, trial <- Y
solve:      write trial/workspace
accept:     Y <- trial
reject:     restore previous/backup
```

## 5. 区域/全球初始化入口

公开入口：

```matlab
model = init.initializeRegional("data/sites/changwu", space);
model = init.initializeRegional("data/sites/changwu", "path/to/regional_static.nc");
model = init.initializeRegional("data/sites/changwu", "path/to/regional_static.nc", ...
                                struct("active_index", landBlock));
```

`space` 路径用于直接给 land columns：

```matlab
space = struct("lat", lat(:), "lon", lon(:), ...
               "mask", mask(:), "cell_id", cell_id(:));
```

NetCDF static 路径支持二维网格：

```text
lat      [y x x] or any grid shape, required
lon      same shape as lat, required
mask     same shape, optional, default all true
cell_id  same shape, optional, default 1:nAllCell
```

当前实现默认按 `mask==true` 打包成 land-only columns：

```text
domain.nCell              = number of active land cells
domain.space.grid_shape   = original grid shape
domain.space.active_index = linear indices of land cells
domain.space.pack         = "land_only"
```

这样求解器只给陆地格点分配数组，后续输出层可用 `active_index` scatter 回二维/全球网格。
当前已有工具：

```matlab
gridData = init.scatterRegionalColumns(columnData, model.built.domain, NaN);
```

它把 land-only columns 反向散回原 `grid_shape`，非陆地格点填 `NaN` 或指定
`fillValue`。这让初始化合同同时具备 gather 和 scatter 两个方向。

为避免全球网格一次性分配过大，NetCDF 路径还支持分块：

```matlab
options = struct("active_index", activeLandIndexSubset);
model = init.initializeRegional(siteDir, ncPath, options);

options = struct("max_cells", 1000);
model = init.initializeRegional(siteDir, ncPath, options);
```

`active_index` 是原始 full grid 的线性索引，必须属于 `mask==true` 的陆地格点。

## 6. NetCDF 覆盖字段约定

当前约定：

```text
soil_<field>             -> built.parameters.soil.<field>
soil_n_<field>           -> built.parameters.soil_n.<field>
veg_root_RAI             -> built.parameters.veg.root.RAI
veg_PFT_type             -> built.parameters.veg.structure.PFT_type
Y_<field>                -> built.Y.<field>
forcing_hourly_<field>   -> built.forcing.hourly.<field>
```

支持形状包括：

```text
node fields:     [NN x nCell] or [NN x full_grid]
layer fields:    [NL x nCell] or [NL x full_grid]
layer state:     [NL x nD x nCell] or [NL x nD x full_grid]
forcing fields:  [nTime x nCell] or [nTime x full_grid]
```

对二维网格，`full_grid` 可以是 `[y x x]`；reader 会按 column-major `active_index`
取出 land columns。

## 7. 当前验证证据

最后一次完整初始化回归：

```text
22 test files
695 checks
0 failures
checkcode CLEAN for modified init/runtime/regional files
```

重点测试：

```text
test_initialize                     21/21
test_prepareRuntime                 47/47
test_regionalInitializationContract 49/49
```

`test_regionalInitializationContract` 覆盖：

- `space struct` 的 2-cell 初始化
- 二维 NetCDF grid + land-mask 打包
- `Y.h` 覆盖
- `soil.theta_s` 覆盖
- `veg.PFT_type` 覆盖
- `veg.root.RAI` 覆盖
- `forcing.hourly.Ta` 覆盖
- workspace trial 从覆盖后的 `Y` 暖启
- runtime 的 `nCell`、canopy water、groundwater buffer 形状
- land-only columns scatter 回二维 grid，并填充 inactive cells
- NetCDF `active_index` block 初始化，只分配选中的 land cells

## 8. Claude 下一步不要做什么

不要把 `run_stemmus_initialize` 式的逻辑重新塞回主循环。

不要让 `runOneStep` 读文件。

不要把 `workspace.step.trial` 当 accepted state；accepted state 只能是 `built.Y`。

不要把 diagnostics 当物理输入读回来。

不要为了区域/全球化给每个格点建一个站点文件夹；静态场应该走共享 NetCDF。

## 9. Claude 下一步应做什么

优先顺序：

1. 定 `runOneStep(model, kt)` 的读写边界。
2. 从老 `COMPASS.m` 中剥离一个最小 water/heat/gas/N step kernel，只消费
   `model.built.Y/workspace/parameters/forcing/boundary/runtime`。
3. 接受/拒绝机制按 `Y <-> workspace.step.trial/previous` 做，不要再用散落的
   `h/hh/hOLD` 全局变量。
4. 输出层后续需要支持 `domain.space.active_index`，把 land-only columns scatter 回
   原二维网格。
5. 真实全球产品映射另建配置层：变量名映射、单位换算、缺测填补、分块读取。

当前初始化层已经可以作为 solver-neutral contract 使用；后续重点是主循环迁移和
真实产品映射，不是再重写初始化。
