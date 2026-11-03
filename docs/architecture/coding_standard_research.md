# COMPASS 编码规范调研报告 — 5 仓库对照

调研日期：2026-05-13
调研对象：ClimaLand.jl、DAISY 7.1.0、NoahMP_Dual、STEMMUS-SCOPE、HYDRUS（仅手册）
调研目的：为 COMPASS MATLAB/Julia 双语言编码规范 (`compass-coding-standard` skill) 提供基础

---

## 0. 一句话定位

| 仓库 | 语言 | COMPASS 应学的核心 |
|---|---|---|
| **ClimaLand.jl** | Julia | **现代陆面模型工程范本**：standalone/integrated 分层、prognostic-vars 与 aux-vars 分离、守恒积分量作为 state、统一 diagnostics 注册 |
| **DAISY 7.1.0** | C++ | **农业 C-N-水耦合范本**：反应/输运分离、参数 `.dai` 配置带单位、写回式松耦合接口可单元测试 |
| **NoahMP_Dual** | Fortran | **传统陆面规范**：参数集中 struct 显式 intent(in/out) 形参，**反例**：12000 行单文件 + 多版本 `.WL/.niu` 同名共存 |
| **STEMMUS-SCOPE** | MATLAB | **COMPASS 的基准**：`+package/` 命名空间，**反例**：22 参数巨签名、状态变量大小写混写 |
| **HYDRUS** | Fortran | 手册公式编号 ⇄ 源码行号 1:1 可追溯；反应-传输算子分裂 |

---

## 1. 顶层目录结构 — 各家怎么切

### ClimaLand.jl ⭐ 最优范本
```
src/
  standalone/        # 单体模型，独立可跑独立可测
    Soil/  Vegetation/  Snow/  Bucket/  SurfaceWater/
  integrated/        # 耦合系统组装层
    land.jl  soil_canopy_model.jl  soil_canopy_root_interactions.jl
  shared_utilities/  # 横切关注点
    Domains.jl  Parameters.jl  drivers.jl  boundary_conditions.jl
    sources.jl  implicit_timestepping.jl  models.jl  checkpoints.jl
  diagnostics/       # 输出层完全独立
    Diagnostics.jl  define_diagnostics.jl  default_diagnostics.jl
test/                # 镜像 src 结构
```
**核心思想**：按"耦合度"切，而非按"物理过程"。耦合代码不写进 standalone，避免单体相互依赖。

### DAISY
```
src/
  daisy/             # 科学模型
    soil/  crop/  chemicals/  organic_matter/  upper_boundary/  manager/
  object_model/      # 反射/工厂/参数解析框架
  ui/  util/  programs/
include/             # 公共头，与 src 镜像
lib/                 # .dai 参数库（63 个文件，作物/化学品/管理参数全外置）
test/                # 三层测试：cxx-unit-tests / dai-unit-tests / dai-system-tests
```
**亮点**：通用基础设施（object_model/util）与科学模型（daisy/）严格分开。

### NoahMP_Dual ⚠️ 反例
```
Noah_code/           # 12000 行单文件 module_Noahlsm.F.WL
IO_code/             # 驱动+NetCDF
Run/                 # MPTABLE.TBL / SOILPARM.TBL / GENPARM.TBL（参数全外置）
```
**亮点**：物理与 I/O 严格隔离 + 参数全外置 `.TBL`
**反例**：`.WL/.niu/.v01` 多版本同名文件共存 — **COMPASS 已踩坑** (ebal.m / ebal_clm45.m / ebal_clm45_profile.m)

### STEMMUS-SCOPE
```
src/
  +soilmoisture/  +energy/  +conductivity/  +equations/
  +groundwater/  +init/  +io/  +parameters/  +helpers/  +plot/  +dryair/
  STEMMUS_SCOPE.m  RTMo.m  RTMf.m  biochemical*.m  fluspect_*.m
test/  docs/  utils/  example_data/
Dockerfile  mkdocs.yml  CHANGELOG.md  miss_hit.cfg
```
**亮点**：`+package/` 切分 + 顶层 `miss_hit.cfg` MATLAB linter 配置
**反例**：根目录还混杂 RTMo/biochemical 等 SCOPE 遗产文件，没全部进 package

### → COMPASS 应学
1. ⭐ **建 `src/+integrated/`**（学 ClimaLand）：把当前散在 `+plant/+PlantHydraulic/calSPAC.m` 的 soil-root-plant 耦合代码搬出来
2. ⭐ **建 `src/+diagnostics/`**（学 ClimaLand）：输出注册统一一处
3. **参数从源码外置到 toml/yaml 配置**（学 DAISY 的 `.dai` 风格）
4. **避免** `ebal_clm45.m` / `ebal_clm45_profile.m` 多版本同名（学 NoahMP 的教训）

---

## 2. 模块边界 — 物理过程怎么分

### ClimaLand: 按"耦合度"切
- **standalone/Soil/** 单独分文件：`rre.jl` (Richards), `energy_hydrology.jl` (能量+水耦合), `retention_models.jl` (参数化), `Biogeochemistry/` (生化子模块)
- **integrated/** 写耦合：`soil_canopy_root_interactions.jl`

### DAISY: N 作为"化学品框架的应用"，**不是独立模块**
- `chemicals/` 下 `Reaction` 接口的子类（`reaction_nit.C` / `reaction_denit.C` / `nitrification.C` / `adsorption_*.C`）
- `organic_matter/` 4 池（SOM/SMB/AOM/DOM），矿化产物作为源项喂入无机 N 池
- `soil/transport/movement_solute.C` 水和溶质共用一套数值框架
- `crop/root/solupt.C` 根吸收 `Solupt` 接口抽象

### STEMMUS-SCOPE: "物理过程 × 算法步骤"双维
每个 PDE 5 件套：assemble 系数 → 边界 → 时间导 → solve 三对角 → 收敛检查
- `+soilmoisture/solveSoilMoistureBalance.m`
- `+dryair/solveDryAirEquations.m`
- `+energy/solveEnergyBalanceEquations.m`

### NoahMP: 单文件内调用树
`NOAHMP_SFLX` (L390 主入口) → `ATM/PHENOLOGY/PRECIP_HEAT/ENERGY/WATER/CARBON`
12000 行单文件 + private 声明 — 现代角度看是反模式

### → COMPASS 应学
1. **N 模块按 DAISY 分子模块**：`+chemistry/` (NH4/NO3 反应+吸附) + `+organic/` (有机质池) + `+transport/` (与水共用) + `+root/+uptake/` (被植物调) 四个对等子模块
2. **每个 PDE 求解都做 5 件套** (学 STEMMUS-SCOPE)：assemble → BC → time-deriv → solve → check_conv
3. **耦合代码搬到 `+integrated/`** (学 ClimaLand)

---

## 3. State / Parameters / Forcing / Auxiliary / Diagnostics 容器

### ⭐ ClimaLand 严格四分离
```julia
struct EnergyHydrology{FT, PS, D, BCRH, S} <: AbstractSoilModel{FT}
    parameters::PS              # 静态参数（不可变）
    domain::D                   # 网格
    boundary_conditions::BCRH   # 强迫/BC
    sources::S                  # 源汇
    lateral_flow::Bool
end
```
- **Parameters**: `EnergyHydrologyParameters` `@kwdef struct`，每个字段一行字符串注释+单位
- **Prognostic state (Y)**: `prognostic_vars(model) = (:ϑ_l, :∫F_vol_liq_water_dt)` 返回 Symbol
- **Auxiliary (p / cache)**: `auxiliary_vars(model) = (:K, :ψ, :total_water, ...)`
- **Diagnostics**: 完全独立输出层，只读 Y/p 不回写

⭐ **关键设计**：state 名声明与 state 内容更新分离，框架自动分配 FieldVector

### DAISY: 状态分布在对象里
- `Chemistry` 持有 `auto_vector<Chemical*> chemicals` + `auto_vector<Reaction*> reactions`
- 每个 `Chemical` 自己管自己的 `M`（总量）、`C`（液相浓度）vector
- 不集中在大 struct

### NoahMP
- `TYPE noahmp_parameters` 一个 struct 装 ~200 参数（L184），显式 `intent(in)` 传
- 瞬态状态散列 `REAL, DIMENSION(...), INTENT(INOUT)` 一长列在 SUBROUTINE 形参里

### STEMMUS-SCOPE: 多 struct 串联
22 参数返回的巨签名（`SoilVariables, HeatMatrices, HeatVariables, ForcingData, ModelSettings, SoilProperties, TimeProperties, InitialValues, GroundwaterSettings, ScopeParameters, Constants, BoundaryCondition`）— **反例**

### → COMPASS 应学
1. ⭐⭐ **prognostic / auxiliary / diagnostic 严格分容器**（学 ClimaLand）
   - 现 `SoilVariables` 把 prognostic (hh/TT/Theta_LL) 和 auxiliary (DTheta_LLh / Khh / Se) 混 — 必须拆
2. ⭐ **参数集中到一个 `core.SoilParameters` struct 显式传**（学 NoahMP），不要散
3. **避免** 22 参数巨签名（学 STEMMUS-SCOPE 反面教训），用对象传

---

## 4. 过程模块函数签名风格

### ⭐ ClimaLand: 工厂函数返回闭包
```julia
function ClimaLand.make_compute_imp_tendency(model::RichardsModel)
    function compute_imp_tendency!(dY, Y, p, t)
        # 写 dY.soil.ϑ_l in place
    end
    return compute_imp_tendency!
end
```
- 接口稳定：`(dY, Y, p, t)` 或 `(p, Y, t)`
- model 在闭包里被捕获，运行时无需再传
- 多分派让不同模型共享时间步进框架

### DAISY: 直接函数调用 + 抽象接口
```cpp
NH4Upt = NH4_uptake->value(geo, soil, soil_water, EffectiveDensity,
                           H2OExtraction, Rad, soil_NH4, PotNUpt,
                           NH4Extraction, MxNH4Up, NH4_root_min);
```
- 接口签名 const 引用，没有全局状态，可单元测试
- 写回式松耦合（写 `NH4Extraction` 向量，由 Chemistry 下一步扣除）
- `chemistry.find("NH4")` 用 symbol 查找 — 土壤未跟踪时自动走 else 降级

### NoahMP: SUBROUTINE 形参表显式 intent
```fortran
SUBROUTINE SOILWATER(parameters, ...)
  TYPE(noahmp_parameters), INTENT(IN) :: parameters
  REAL, INTENT(INOUT) :: SH2O(NSOIL)  ! [m3/m3]
  ...
```
- `!in/!out/!inout` 显式声明
- 行末 `[unit]` 注释

### → COMPASS 应学
1. ⭐ **MATLAB 没有多分派**，但可约定**每个模块固定签名**：
   ```matlab
   function [state, aux] = update_aux(state, aux, params, forcing, t)
   function [dstate] = compute_tendency(state, aux, params, forcing, t)
   ```
2. ⭐ **函数签名注释强制写明** `% in/inout/out` (学 NoahMP)
3. ⭐ **接口测试**：每个公共函数能不依赖全局状态独立调用 (学 DAISY 接口)

---

## 5. 参数对象设计

### ⭐ ClimaLand: 不可变 struct + 强类型 + toml 双构造
```julia
Base.@kwdef struct EnergyHydrologyParameters
    "The saturated hydraulic conductivity (m/s)"  # 单位写在 docstring
    K_sat::F
    ...
end
```
- 直接构造 + toml 构造两条路径
- 用 `ClimaParams.jl` 一套 toml 共享物理常数
- 模型只装"模型特有"参数

### ⭐ DAISY: `.dai` 配置带单位解析
```
(defchemical Atrazine herbicide
  (decompose rate 7e-5 [h^-1])
  (adsorption linear (K_clay 7 [cm^3/g])))
```
- 单位写在方括号里参与解析
- 由 `object_model/units.C` 校验
- 嵌套子模型选择（`(adsorption linear ...)`）一体化

### NoahMP: `.TBL` 表
`MPTABLE.TBL` 植被参数，`SOILPARM.TBL` 土壤参数 — 全外置

### → COMPASS 应学
1. ⭐⭐ **配置文件支持单位标注**（学 DAISY）— 当前 `config_*.txt` 没有
2. ⭐ **参数 struct 不可变**（学 ClimaLand）— 当前 `core.ModelSettings` 等是 handle class，可写
3. **物理常数集中一个 `core.PhysicalConstants`**（已有 ✓）

---

## 6. doc string / 注释规范

### ⭐ ClimaLand 三件套
```julia
"""
    c3_rubisco_assimilation_farquhar(...) 

Computes the Rubisco limiting rate of photosynthesis for C3 plants (`Ac`),
in units of moles CO2/m^2/s, as a function of ...
The units of ci, Γstar, and Kmm may be Pa or mol/mol, but they must be consistent.
See Table 11.5 of G. Bonan's textbook, Climate Change and Terrestrial Ecosystem Modeling (2019).
"""
```
**三件套**：(1) 每个输入参数解释一句；(2) **单位明确**（强调一致性）；(3) **文献引用具体到表/章节**

### ⭐ DAISY: 简化但够用
- 文件头 1 行用途
- 物理量后注 `[unit]` 单位
- 变量名贴论文符号（`NH4Upt` / `WLeaf` / `PotNUpt`）
- **不用 Doxygen 标签**

### NoahMP
- `! in/!out/!inout` 形参方向
- 行末 `[m3/m3]` 单位
- 无 Doxygen/FORD

### STEMMUS-SCOPE
- 多行块 `%{ ... %}` 引论文公式编号
- 例 `solveSoilMoistureBalance.m` L3-L8: "STEMMUS Technical Notes section 4, Eq 4.32"

### HYDRUS
- 手册公式编号 ⇄ 源码行号 1:1 可追溯
- 变量名沿论文符号（`Th=θ、Con=K、Cap=C`）

### → COMPASS 应学
1. ⭐⭐ **函数 doc 三件套**（学 ClimaLand）：参数解释 + 单位 + 文献页码
2. ⭐ **变量行末注 `[unit]`**（学 NoahMP/DAISY）
3. ⭐ **公式函数引论文章节号**（学 STEMMUS-SCOPE/HYDRUS）
4. **变量名贴论文符号**（学 DAISY 的 `WLeaf`、`PotNUpt`）

---

## 7. 测试组织

### ⭐ ClimaLand: 镜像 src + 每模块多文件
```
test/standalone/Soil/
  soiltest.jl              # 单元测试主流程
  soil_parameterizations.jl
  soil_bc.jl
  climate_drivers.jl
  conservation.jl          # ⭐ 守恒性专门一个文件
  parameters.jl
test/integrated/
  soil_canopy_lsm.jl
```
- 单元测试 + 集成测试分目录
- **守恒检查是 CI 一部分**

### DAISY: 三层测试
`test/cxx-unit-tests/`（C++ 单元）+ `test/dai-unit-tests/`（DSL 单元）+ `test/dai-system-tests/`（系统级）

### → COMPASS 应学
1. ⭐⭐ **建 `test/+soil/`、`test/+plant/`、`test/+integrated/`**（镜像 src）
2. ⭐⭐ **每个模块至少 `conservation.m` + `parameterizations.m` + `bc.m` 三件套**
3. **当前 `tools/` 下的脚本式 verify_* / diagnose_*** 应该规范化为 `test/`

---

## 8. 守恒检查 / Diagnostics（COMPASS 最缺）

### ⭐⭐ ClimaLand: 累积通量提升为 prognostic state
```julia
# rre.jl:282
prognostic_vars(soil::RichardsModel) = (:ϑ_l, :∫F_vol_liq_water_dt)

# rre.jl:166 在 compute_imp_tendency! 里
dY.soil.∫F_vol_liq_water_dt = -(top_bc - bottom_bc)
```
- ⭐ **累积通量本身是状态变量**
- 守恒残差 = `(total_water[end] - total_water[0]) - ∫F_vol_liq_water_dt` 在后处理直接算
- **积分量是模型内置 state，不靠后处理重新积分**

### ⭐ ClimaLand: 统一诊断注册
`src/diagnostics/define_diagnostics.jl:9-47`:
```julia
add_diagnostic_variable!(
    short_name = "wvpa",
    long_name = "water_volume_per_area",
    units = "m",
    comments = "...",
    compute! = ...
)
```
- 全局注册 short_name / long_name / units / comments / compute! 五件套
- 物理模块只算不写文件

### 通用接口
`shared_utilities/models.jl:498` 的 `total_liq_water_vol_per_area!(cache, model, Y, p, t)` 是默认 fallback (写 0)，各 standalone 模型 dispatch 重写

### → COMPASS 应学
1. ⭐⭐⭐ **把 cum_E / cum_T / cum_drainage / cum_R_n 提升为 state**（学 ClimaLand）
   - 物理模块每步 in-place 累加
   - 不要在 `COMPASS.m` 末尾从输出文件后处理积分
2. ⭐⭐ **`src/+diagnostics/registerDiagnostic.m`**：所有输出在一处注册
3. ⭐⭐ **`test/+conservation/`**：每个 module 守恒 CI 测试

---

# 总结：COMPASS 应当吸收的 9 大原则

按 ROI 排序：

| # | 原则 | 主要参考 | COMPASS 现状 |
|---|---|---|---|
| 1 | **累积通量提升为 prognostic state** | ClimaLand | ❌ 现在后处理积分 |
| 2 | **prognostic / auxiliary / diagnostic 严格分容器** | ClimaLand | ⚠️ SoilVariables 混着 |
| 3 | **diagnostics 注册统一一处** | ClimaLand | ❌ 散在 COMPASS.m |
| 4 | **守恒 CI 测试** | ClimaLand, DAISY | ❌ 无 |
| 5 | **函数 doc 三件套** (参数+单位+文献) | ClimaLand | ⚠️ 部分有 |
| 6 | **N 模块按"反应/输运/吸收"切而非一个文件** | DAISY | ❌ 现 STEMMUS_N.m 一锅 |
| 7 | **参数集中 struct 显式 intent** | NoahMP | ⚠️ 部分用 cfg.* 但散 |
| 8 | **避免多版本同名文件** | NoahMP 反例 | ⚠️ ebal_*.m 有多版本 |
| 9 | **配置文件带单位标注** | DAISY | ❌ config_*.txt 裸数值 |
| 10 | **模块独立可跑性**：standalone 不依赖耦合 | ClimaLand, DAISY | ❌ calSPAC.m 必须有 SoilVariables |
| 11 | **代码 / 参数 / 运行 三分离** | ClimaLand, DAISY | ✅ 已初步分（src/data/run_*）|

---

## 9. 渐进迁移路径（COMPASS 现实约束）

ClimaLand 不是一夜重构出来的（GitHub history 看至少 4 年演化）。COMPASS 必须**渐进迁移**，规则是：

### 不破不立：现有代码不动，新代码强制规范
- `src/COMPASS.m` 不重写（103 处 Sim_output 赋值 + 6 处 NBChh 切换都依赖现状）
- `src/ebal_clm45*.m` 4 个版本暂保留，**禁止再新增 v3/v4**
- `src/STEMMUS_N.m` 1387 行老 N 入口保留，**新 N 模块写在 `+soil/+nitrogen/` 包下**

### 三类工作的规范级别
| 工作类型 | 规范严格度 |
|---|---|
| 新增 src/ 文件 | ⭐⭐⭐ 严格按 skill 规范 |
| 重构现有模块（如 N） | ⭐⭐⭐ 严格按 skill 规范（重构本身是新代码） |
| 修复 bug（如 Thomas SOR） | ⭐ 改最小行数，不强求规范化整个文件 |
| 一次性诊断脚本 | 不要求 |

### 优先重构顺序（基于 ROI + 当前痛点）
1. **N 模块**（你已经定的下一步）— 1387 行 STEMMUS_N.m 拆成 `+chemistry/` + `+organic/` + `+transport/` + `+root/+uptake/`
2. **诊断输出层** — 把 COMPASS.m 里 103 处 Sim_output.* 赋值搬到 `src/+diagnostics/registerXxx.m`
3. **ebal 整合** — `ebal*.m` 4 版本合并为 `+canopy/+energy/solveEbal.m`，子方案通过 `cfg.ebal_scheme` 切（`'clm45'/'clm45_profile'/'classic'`）
4. **SoilVariables 拆分** — 长期目标，拆成 `+core/SoilProgState.m`（hh/TT/Theta_LL）+ `+core/SoilAuxState.m`（Se/Khh/DTheta_LLh）

---

## 下一步

基于本调研，下一文档将是：**`compass-coding-standard` skill**（双语言严格规范），包含：
- 每条原则的"为什么 + ✅/❌示例 + 当前 COMPASS 哪里要改"
- 模块创建脚手架模板
- doc string / 注释模板
- 命名约定决策树
- CLAUDE.md 接通最高警示

预计字数 5000-8000，写完前用户 review 一次。
