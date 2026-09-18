# Pre-soil-loop 变量归属审计 / Variable Ownership Audit (进 ④ 前钉边界)

> 进入真正土壤求解循环(④, `attemptSolve` stub 的真实实现)之前,把 `runOneStep` ③ 块涉及的每个变量
> 按 ClimaLand/Julia 风格分类,钉死 **state / aux / cache / runtime-service / diagnostics** 的边界。
> **本文档是审计,不改物理。** 结论分三档:**保留(现状可)** · **须迁(标目标)** · **随 ④(冻融/水解一体)**。
> 配套 Julia 迁移：`.claude/skills/compass-julia-port` + `docs/architecture/compass_complete_architecture.md`。
> Editor: Jizixing Zhu; Date: 2026-07-05.

## 0. 七类归属定义(判据)

| # | 类 | 判据 | checkpoint | rollback | Julia target |
|---|---|---|---|---|---|
| 1 | **State / Y** | 跨步物理/生物状态,推进下一步 | 是 | 是(失败步) | `ModelState` (prognostic) |
| 2 | **MemoryState** | 跨步"记忆"(非 PDE 但影响未来物理) | 是 | 是 | `ModelState.memory` / `MemoryState` |
| 3 | **parameters** | 静态参数,运行不变 | 否(config 重建) | 否 | `ModelParameters` (immutable) |
| 4 | **p / cache** | init 后预计算缓存,可重建 | 否(rebuild) | 否 | `ModelCache` / `PrecomputedParameters` |
| 5 | **forcing/boundary/sources** | 外部驱动/边界/管理输入 | 否(重读) | 否 | `Drivers` |
| 6 | **workspace / aux** | 当步临时/派生量/物理函数输入包 | 否(rebuild) | (步内 snapshot) | `Aux` (每步 `update_aux!`) |
| 7 | **runtime.services** | IO/checkpoint/calendar/debug/run_control | 部分(计数器/指针) | 否 | `RuntimeServices` |

**核心原则**:`会影响未来物理推进、需 rollback/resume 的 → State/MemoryState` · `当步临时 → Aux` ·
`init 预算 → Cache/p` · `纯运行服务 → RuntimeServices` · `只解释不推进 → Diagnostics`。

---

## 1. 主归属表(按容器)

### 1.1 `built.Y.*` — prognostic State(已对,整块进 checkpoint + snapshot rollback)
| field | 类 | 说明 | Julia target |
|---|---|---|---|
| `h, T, Theta_L, Theta_U, p_gas, theta, bd` | State | 土壤水/热/气 prognostic + 时间层 | `State.soil` |
| `nh4, no3, urea, C_NH4_granule` | State | 土壤 N prognostic(node 底上) | `State.soil_n` |
| `soil_gas` | State | 气相 prognostic | `State.soil` |
| `canopy_air, scope_soil` | State | canopy-air OLD 层 / SCOPE 土壤热史 | `State.canopy` / `State.soil` |
| `plant_cn` | State | 植物 C/N 池 + accum + state(PR-D) | `State.plant_cn` |
| `plant_hydraulic` | State | psi_* + W_{root,stem,leaf}(PR-A V2) | `State.plant_hydraulic` |
| `surface_water` | State | DSTOR/DSTOR0/RS 积水(2026-06-28) | `State.surface` |
| `soil_modes` | State | 滞回 Lenhard RHSW/RASW/…(PR-F 审计) | `State.soil.hysteresis` |
✅ **Y 分类干净**,整块 snapshot(rollback)+ 整块 checkpoint。**保留。**

### 1.2 `built.parameters.*` — parameters(静态)
`soil / soil_n / veg / scope / plant_cn`。✅ **保留**(不变、不 rollback、config 重建)。→ Julia `ModelParameters`。
⚠️ `veg` 是灰区(见 §2.10):内含**参数 + 驱动 + init 模板**混杂,须拆(现在标注,暂不搬)。

### 1.3 `built.p.*` — cache(预算,可重建)
| field | 类 | 说明 | 判定 |
|---|---|---|---|
| `plant_hydraulic` | cache | PSIL_RANGE/VC/KS/gpmax 脆性查表 | ✅ 保留 |
| `plantHydro` | **legacy adapter/cache** | 喂 ebal/V2 老函数的旧形态包 | ⚠️ **legacy adapter**(已 doc 标),勿长期化(见 §2.11) |
→ Julia `ModelCache`。

### 1.4 `built.workspace.*` — aux(当步)
| field | 类 | 说明 | 判定 |
|---|---|---|---|
| `derived` | aux/cache | Se/KL_h/Theta_LL/**Theta_II**/… 每步 refreshSoilAuxiliary 重算 | ✅ 保留(可重建,不 checkpoint) |
| `flux` | aux/source | Srt/SinkS(→④ source)+ PR-E Throughfall/EVAP/CanopyInterception | ✅ 保留(当步 source 项) |
| `surface` | aux/scratch | DSTMAX/NBChh/h_target/hN 顶 BC scratch | ✅ 保留(每步 reset,snapshot 覆盖) |
| `solver` | aux | PDE 系数 bank + assembled systems | ✅ 保留 |
| `scope_step` | aux/step-scratch | 每步 prepareCanopyStep 重建的 SCOPE/ebal 输入+输出包 | ✅ 保留(**非 accepted state**,见 §2.12) |
| `step` | aux+carry | previous(t_{n-1} 锚)/ trial(子步解)/ TIME | ⚠️ previous/trial 是 carry-commit 三元组,概念上属 State 的时间层 |
→ Julia `Aux`(+ `step.previous/trial` → `State` 的 OLD/Trial 三元组,skill trap #5)。

### 1.5 `runtime.*` — **混了 service + memory-state(最大灰区,见 §2 + §3)**
| service(纯运行) | memory-state(影响未来物理,restart-critical) |
|---|---|
| output, diagnostics, checkpoint, resume, calendar, run_control, debug, substep, status, dimensions, created_at, picard | ~~canopy_water(state)→Y, snow→Y, thermal→Y, phenology→Y, lai_state→Y, luna(state/last_*/climate_buffer)/leaf_resp→Y.plant_acclimation~~(增量1-6b **全搬完**), **management(once-only), groundwater**(仅剩这两个 restart-state) |
⚠️ **这是本次审计的重点**:物理记忆现塞在 `runtime` 里,Julia 化会分不清 service vs state。§2/§3 逐个判。

### 1.6 `built.{forcing, sources, boundary, domain, vegetation_state}`
| field | 类 | 判定 |
|---|---|---|
| `forcing.hourly.*`(Ta/Precip/RH/…) | forcing | ✅ 保留 → `Drivers.meteo` |
| `sources.*`(施肥事件等) | sources | ✅ 保留 → `Drivers.management` |
| `boundary.*`(thermal_bottom/…) | boundary | ✅ 保留 → `Drivers.boundary` |
| `domain.*`(grid/site/dz/NN/…) | parameters(几何) | ✅ 保留 → `ModelParameters.domain` |
| `vegetation_state` | **混杂** | ⚠️ **须拆**(§2.10) |

### 1.7 `pending.*` — 本步事务(A' 已建)
`prognostic / runtime / workspace_flux / scope_flux / diagnostics / management`。✅ 唯一提交边界(§7)。→ Julia `StepPending`。

### 1.8 checkpoint 存/不存边界(已对)
- **存**:`Y`(整块;snow/surface/canopy_water/vegetation(phenology+lai_state)/plant_acclimation(luna_state+leaf_resp+luna_climate_buffer) 状态已随 Y 整块走,2026-07-05~07 增量1-6b)+ `runtime.{step,substep,calendar,groundwater,management}`(**runtime 不再存任何 LUNA/物理记忆**)。
- **不存(重建)**:workspace/aux(refreshSoilAuxiliary 从 Y 重算)+ diagnostics Registry(从 resume KT 重填)+ canopy_water histories。
✅ 边界清楚。**但注意**:checkpoint 现在**手抄了一份 runtime restart 白名单** —— 正好印证 §3(这些 runtime 字段其实是 restart-state,该系统化)。

---

## 2. 灰区逐个判定(结论明确)

| # | 变量 | 现位置 | 应属 | checkpoint | rollback | 判定 |
|---|---|---|---|---|---|---|
| 1 | `canopy_water.state_h2ocan*` | runtime | **MemoryState** | 是 | 是 | ✅ **已搬 → Y.canopy_water(2026-07-05,增量#3)**;随 ckpt.Y;regional [1×nCell] |
| 1b | `canopy_water.*_history/residual` | runtime | **diagnostics** | 否 | 否 | ✅ **已拆 → built.diagnostics.canopy_water(2026-07-05,增量#3)**;row-write accept-only,不进 Y |
| 2 | `snow.SWE` | runtime | **MemoryState** | 是 | 是 | **须迁** → Y.snow.SWE(**风险#1,单标量,最先**) |
| 3 | `thermal.Tss` | runtime | **MemoryState** | (现未单存) | 是 | ✅ **已搬 → Y.surface.Tss(2026-07-05,增量#2)**;随 ckpt.Y 补上 checkpoint(修了漏存 gap) |
| 4 | `phenology.now/state` | runtime | **State(veg)** | 是 | 是 | ✅ **已搬 → Y.vegetation.phenology(2026-07-05,增量#5)**;整 struct,随 ckpt.Y |
| 5 | `lai_state`(含 wt_day_min) | runtime | **State(veg)** | 是 | 是 | ✅ **已搬 → Y.vegetation.lai_state(2026-07-05,增量#5)**;全字段,随 ckpt.Y |
| 6 | `luna.state`(+ 6 `last_*`) | runtime | **MemoryState** | 是 | 是 | ✅ **已搬 → Y.plant_acclimation.luna_state(2026-07-05,增量#6a)**;5 base+6 last_* 固定字段,last_* 初值 NaN,随 ckpt.Y |
| 7 | `luna.climate_buffer` | runtime | **MemoryState(handle)** | 是(toStruct) | 是 | ✅ **已搬 → Y.plant_acclimation.luna_climate_buffer(2026-07-06,增量#6b)**;struct-in-Y(Option A:accept-path fromStruct/push/toStruct,不改 Climate10dBuffer 类);随 ckpt.Y;**runtime.luna 整删** |
| 8 | `leaf_resp.Ta_history` | runtime | **MemoryState** | 是 | 是 | ✅ **已搬 → Y.plant_acclimation.leaf_resp.Ta_history(2026-07-05,增量#6a)**;config(atkin/Tgrowth_window)→ built.p.plant_acclimation;随 ckpt.Y |
| 9 | `scope_static` | runtime | **cache/p** | 否 | 否 | ✅ **已搬 → built.p.scope_static(2026-07-05,增量#4)**;init 静态缓存,不推进/不 rollback/不 checkpoint |
| 10 | `vegetation_state` | built | **混杂** | — | — | **须拆** → parameters.veg / forcing.veg_driver / Y.vegetation / p.veg_cache(§5) |
| 11 | `p.plantHydro` | built.p | **legacy adapter** | 否 | 否 | **保留但标死**:adapter/cache,Julia 拆 typed 参数+state(§6) |
| 12 | `workspace.scope_step` | workspace | **aux(step-scratch)** | 否 | (步内) | **保留**:每步重建,**非 accepted state**(别 checkpoint) |
| 13 | `workspace.flux` | workspace | **aux/source** | 否 | (步内) | **保留**:④ source 项 |
| 14 | `workspace.derived` | workspace | **aux/cache** | 否 | (步内) | **保留**:refreshSoilAuxiliary 重算 |
| 15 | `management`(fertilizer.applied) | runtime | **MemoryState(once-only)** | 是 | 是 | **保留在 restart**;已走 pending accept-only |
| 16 | `groundwater` | runtime | State/memory | 是 | ? | 待查(本次未展开,标 TODO) |

**关键澄清**(用户点名的三点):
- `canopy_water.histories` = **输出历史,不是 state** → diagnostics/workspace,**不进 prognostic Y**。
- `luna.climate_buffer` = **restart-critical memory**(影响未来 Vcmax/Cab),必须 checkpoint/rollback。
- `phenology / lai_state` = **state**(影响未来 LAI),不是普通 runtime service。

---

## 3. runtime → services / restart_state 语义分层(§2 用户建议)

**目标命名**(先建分类,字段搬迁按风险分批;checkpoint 白名单 = restart_state,天然对齐):
```
runtime.services.{output, checkpoint, resume, diagnostics, calendar, run_control, debug, substep, status, dimensions}
runtime.restart_state.{canopy_water(state), snow, thermal, phenology, lai_state, luna, leaf_resp, plant_acclimation, management, groundwater}
  (2026-07-05~07 增量1-6b 落地: canopy_water/snow/thermal/phenology/lai_state -> Y; luna(state/last_*/climate_buffer)+leaf_resp
   -> Y.plant_acclimation; runtime 只剩 management + groundwater 两个 restart-state,其余全进 ckpt.Y)
```
- `services` 不碰物理,只管理运行。`restart_state` = 暂存命名,**最终大部分进 Julia `State`/`MemoryState`**。
- checkpoint/resume 必须能区分两类(现在 checkpoint 手抄白名单 = 已隐式区分,系统化即可)。
- **落地方式**:可先**只在文档 + 一个 `classifyRuntimeField` helper** 建立分类,不搬字段(低风险);Julia 版严禁把物理 memory 混进 service。

---

## 4-9. 逐项迁移方向(结论,按风险排序;实际搬迁待用户逐个放行)

- **§4 `scope_static` → `built.p.scope_static`** ✅ **已搬(2026-07-05,增量#4)**:init 静态 SCOPE/RTM 缓存,不推进/不 rollback/不 checkpoint。抽 `init.buildScopeStatic`(逐字搬旧 prepareScopeStatic)在 buildState 建;prepareRuntime 变纯 runtime allocator(不建 cache);6 读点(prepareScopeRuntime/prepareScopeStep/runOneStep×5)改读 built.p。Julia = `Cache`/`PrecomputedParameters`。gate 18/0(test_rtmo_parity 行为兜底 + grep 0 live reader)。
- **§5 `vegetation_state` 拆**:`parameters.veg`(静态)/ `forcing.veg_driver`(LAI/hc/物候时间序列)/ `Y.vegetation`(跨步变)/ `p.veg_cache`(预算)。暂标"待拆分"。
- **§6 `p.plantHydro` = legacy adapter**:doc 已标(本会话)。静态参数应来自 parameters/p;`.init.*` 只作旧函数输入,**不能反向成状态来源**。Julia 不照搬大结构。
- **§7 `pending` = 唯一提交边界**(A' 已建 6 层 + `assertRequiredPending` fail-fast)。下一步细化:`commitStepPending` 可再拆 `commitPrognostic/commitRuntimeRestart/commitWorkspaceFlux/commitScopeFlux/commitDiagnostics/commitManagement`(**可选**,当前单函数分 8 步已清晰;拆分收益递减,视 Julia 需要再定)。原则已满足:失败步零提交,accept 才提交,所有跨步 memory 走 pending。
- **§8 单位转换集中**:另建 `docs/unit_boundaries.md` 登记(`U_N×1e-2`/`EVAP×0.1`/`Resp×12/1e9`/Throughfall cm·s⁻¹/SWE cm/GPP kgC↔umol)。**不改物理**,Julia 集中到 `Units`/boundary adapter。
- **§9 Julia shadow port**:`julia/CompassPreSoilStep/` 只转结构(typed structs + ③ 骨架),物理 stub/parity adapter,**逼类型分类清楚**。这是**验证本审计分类是否合理**的手段,不跑完整模型。

**迁移风险序(搬字段时逐个,每个都改 prepareRuntime/buildInitialState + pending + commitStepPending + checkpoint + beginRun + 测试,并保旧 checkpoint 映射)**:
1. ✅ **`snow.SWE → Y.snow.SWE`(2026-07-05 已搬,模板落地)** 2. ✅ **`thermal.Tss → Y.surface.Tss`(2026-07-05 已搬,增量#2)** 3. ✅ **`canopy_water` 3路拆(state→Y.canopy_water / histories→built.diagnostics / config→params.veg.canopy;2026-07-05 增量#3)** 4. ✅ **`scope_static → built.p.scope_static`(2026-07-05 增量#4;静态 cache,非 state)** 5. ✅ **`phenology+lai_state → Y.vegetation`(2026-07-05 增量#5;两 plain-struct 一个域容器)** 6a. ✅ **`luna.state/last_* + leaf_resp → Y.plant_acclimation`(2026-07-05 增量#6a;config→built.p.plant_acclimation;runtime.leaf_resp 删,runtime.luna 减到 climate_buffer)** 6b. ✅ **`luna.climate_buffer(handle)→Y.plant_acclimation.luna_climate_buffer`(2026-07-06 增量#6b;struct-in-Y,accept-path fromStruct/push/toStruct;删 runtime.luna + snapshotLuna)**。

> **✅ pre-④ 状态归属收口 COMPLETE(2026-07-06):** 增量#1–#6b 全部落地。**所有跨步物理记忆**现在都住在 `Y.<domain>`(soil/snow/surface/canopy_water/vegetation/plant_acclimation/plant_cn/plant_hydraulic 等)、静态派生缓存住 `built.p.*`、只写诊断住 `built.diagnostics.*`;`runtime.*` 只剩**服务**(output/diagnostics/checkpoint/resume/calendar/run_control/debug/substep/dimensions/picard/soil_modes/scope)+ 两个 once-only/边界 restart-state(management/groundwater),**不再拥有任何 LUNA 或物理记忆态**。checkpoint 手抄白名单相应收缩(只 step/substep/calendar/groundwater/management 随 ckpt.runtime,其余随 ckpt.Y)。**下一笔债 = ④ substep solve(attemptSolve 空桩),不是 ownership。**

> **增量 #6a 落点(plant_acclimation,2026-07-05):** 最后一个 ownership 块,拆 6a(plain-state)+6b(handle)。**6a**:`leaf_resp.Ta_history`(逐步生长温度记忆)+ `luna.state`(Vcmax25_z/Jmax25_z/PNlc_z/enzs_z/first_day)+ 6 个 `last_*` 缓存(原懒创建,现提升为**固定字段集,初值 NaN**,`~isfinite(last_Cab_nlc)` 为门)搬进 `Y.plant_acclimation`;4 个派生 config(`luna_enabled/steps_per_day/leaf_resp_atkin_enabled/leaf_resp_Tgrowth_window`)→ 新 `init.buildPlantAcclimationConfig` 建的 `built.p.plant_acclimation`(**不 checkpoint**,init 重建)。**`runtime.leaf_resp` 整删**;**`runtime.luna` 减到只剩 `climate_buffer` handle**(临时 adapter,6b 删)。读点:`prepareScopeStep.buildStepLeafBio`/`applyDynamicCabLeafOptics` 读 lagged `Y.plant_acclimation.luna_state`,`applyTgrowth` 读 `Y.plant_acclimation.leaf_resp.Ta_history`(config 从 built.p);`runDailyLuna`(6a 签名)读 config(built.p)+ buffer(runtime handle)+ state(Y),**只写** `Y.plant_acclimation.luna_state`。commit:`leaf_resp_Ta` 走 `pending.prognostic`,**row-update** `Y.plant_acclimation.leaf_resp.Ta_history(kt)`(非整块覆盖);buffer push + `runDailyLuna` 仍在 accept 路径。`saveCompassCheckpoint`:`snapshotLuna` 减到只存 `climate_buffer`,删 `stateOnlyLeafResp`(state/leaf_resp 随 ckpt.Y)。`beginRun` compat:捕获 fresh-init `Y.plant_acclimation`(固定字段+last_*=NaN)作默认,`Y=ckpt.Y` 后若缺则 overlay 旧 `ckpt.runtime.{leaf_resp.Ta_history, luna.state, luna.last_*}`(缺的保 fresh,尤其 last_*=NaN);**仍保留** `runtime.luna.climate_buffer = fromStruct(...)` handle 恢复(6b 才移入 Y)。**顺带修**:`last_*` 现进 Y → 随 ckpt.Y checkpoint,补掉旧 `snapshotLuna` 不存 last_* 的 resume gap(连续跑零变化,deliberate improvement)。**regional 不广播**(单列)。gate 11/11(含 failstep 双钉 luna_state+Ta_history(kt) 不变 + 宽 grep 抓到 `test_prepareRuntimeSlots`/`test_luna_chain`/`test_runOneStep` 存量读点并 retarget)。**6b 待**:`climate_buffer` handle → `Y.plant_acclimation.luna_climate_buffer`(struct-in-Y,Option A accept-path 临时重建),删 `runtime.luna`/`snapshotLuna`。

> **增量 #6b 落点(luna_climate_buffer,2026-07-06):** ownership 收口最后一步 —— 把 `Climate10dBuffer` **handle** 搬进 Y。Y 不能放 handle,故存**序列化 plain struct**:`buildInitialState` 生 `Y.plant_acclimation.luna_climate_buffer = Climate10dBuffer(steps_per_day,10).toStruct()`。**Option A**(struct-in-Y + accept-path 临时重建,**不改 Climate10dBuffer 类**):`commitStepPending` step-6 只在 **accept 分支** `buf=fromStruct(Y...luna_climate_buffer); buf.push(kt,payload); [ls,bufNext]=runDailyLuna(model,kt,buf); Y...luna_state=ls; Y...luna_climate_buffer=bufNext`。`runDailyLuna` 改 6b 签名 `[luna_state_next,buffer_struct_next]=runDailyLuna(model,kt,buf)`:显式收 buf,**总返** `buf.toStruct()`(即使早退,给 caller 单一写源),自己只读 buf(isWarm/assembleLunaInput)不 push。**`runtime.luna` 整删**(prepareRuntime 删 prepareLuna+行);**`snapshotLuna` 删**(buffer 随 ckpt.Y)。`beginRun` **两层 compat**:(a) 容器缺失分支(pre-6a ckpt,无 Y.plant_acclimation)内 overlay 旧 `ckpt.runtime.luna.climate_buffer`;(b) **子字段缺失分支**(6a-era ckpt,Y.plant_acclimation 有但缺 luna_climate_buffer,buffer 仍在 ckpt.runtime.luna)—— 容器分支后新增,overlay 旧 buffer,缺则回落 fresh 空 buffer。**用户抓的 blocker**:只处理容器缺失会漏 6a-era ckpt → 断链;两层修好。gate 11/11(含 failstep `isequaln(buffer,bufBefore)` 完全不变 + checkpoint_compat 双分支各测 filled==2/==3 + getMeans bit-exact,旧 buffer 用真实 push() 构造保内部一致)。**6b cleanup tail(2026-07-07 已收)**:`pending.runtime.luna_climate` → **`pending.prognostic.luna_climate_entry`**(它是写进 `Y.plant_acclimation.luna_climate_buffer` 的一条 entry payload,非 runtime destination);`pending.runtime` 现为空层(保留给未来真正的 runtime-service payload);`assertRequiredPending` 删 runtime required loop,改判 `pending.prognostic.luna_climate_entry` 非空。纯结构重命名,零物理,gate 11/11。**至此"pending 按目标容器分层"无语义例外。**

> **增量 #5 落点(vegetation,2026-07-05):** 两个 plain-struct 跨步物理记忆(`phenology`={state(PhenologyState 工厂,非 handle),year_mtemp_max,year,now};`lai_state`=Method-B EPIC/frPHU + wt_day_min)搬进**一个** `Y.vegetation` 容器。`buildInitialState` 生 `Y.vegetation`(bit-equiv 旧 preparePhenology/prepareLaiState,整 struct);`prepareRuntime` 删两 local+两行;4 读点分清 **write-pending**(advanceDailyPhenology/applyMethodBLai 产 phenNext/laiStateNext → `pending.prognostic.vegetation`,不写 Y)vs **lagged read-old-committed**(applyCanopyLAI/stepPlantCN 读 `Y.vegetation.phenology` 上一 accepted,保 1-step lag);commit step-5e 原子写 `Y.vegetation`;`saveCompassCheckpoint` 停存 `runtime.{phenology,lai_state}`(随 ckpt.Y);`beginRun` compat **捕获 fresh-init `Y.vegetation` 作默认**(复杂嵌套 struct,不手拼)再 overlay 旧 `ckpt.runtime.{phenology,lai_state}`。**regional 不广播**(phenology/lai_state 今天非 per-cell;多-PFT 物候独立设计)。gate 25/0(含 failstep 两边钉 + 宽 grep 抓到并 retarget 8 处存量测试读点)。

> **增量 #3 落点(canopy_water,2026-07-05):** 前两个增量是单一容器整搬;canopy_water 是 **3 路拆分**。read-back 分析定性:唯一读回的是 `state_h2ocan/_sun/_sha`(→ `Y.canopy_water`,required prognostic);8 个 `*_history/residual` 只写不读(→ `built.diagnostics.canopy_water`,commit 按 kt **row-write accept-only**,复用 N-diag 惯例);静态 `Ps_max/Ps_min_diurnal`(→ `built.parameters.veg.canopy`,`Ps_max` 经新 `io.getEnvDouble` 保 `COMPASS_>STEMMUS_>0.7` 语义),`scheme` 不复制(单源)。**`runtime.canopy_water` 整删**。2 阶段执行:stage1 config+history 搬出(state 暂留 runtime),stage2 state→Y + `materializeCheckpointStateDefaults(Y,rt,nCell)` **逐字段** compat(旧 ckpt state 本就进 checkpoint → resume bit-exact)。`computeCanopyInterception` 保 8 参签名,caller 组 **cw 参数视图**(state 从 Y、config 从 params),函数体 bit-exact。**regional**:`expandBuiltToRegional.broadcastState` 的 nested 白名单补 `snow/surface/canopy_water`(顺带修增量1-2 遗留:三者 regional 下须 `[1×nCell]`)。gate 25/0 绿(含新 `test_psMaxEnv` + `test_regionalInitializationContract` 扩)。

> **增量 #1 落点(snow.SWE,2026-07-05):** `buildInitialState` 生 `Y.snow.SWE=0`;`prepareRuntime` 删 `runtime.snow`;`stepCanopyWaterSnowEvap` 读 `built.Y.snow.SWE`;`emptyStepPending` 把 snow 从 runtime 层移到 prognostic 层;③-block 存 `pending.prognostic.snow`;`commitStepPending` step-5b 写 `Y.snow`(required);`assertRequiredPending` 查 `pending.prognostic.snow`;`saveCompassCheckpoint` 停存 `runtime.snow`(随 `ckpt.Y` 整块走);`beginRun.applyResumeCheckpoint` 经 `materializeCheckpointStateDefaults(Y, ckpt.runtime)` 兼容旧 checkpoint(`rt.snow.SWE`→`Y.snow.SWE`,否则 init 默认 0)。**零物理变更**,gate 9/9 绿。`materializeCheckpointStateDefaults` 是每个后续增量的**唯一扩展点**(一 field 一 clause)。spec/plan:`docs/superpowers/{specs,plans}/2026-07-05-*`。

> **增量 #2 落点(thermal.Tss,2026-07-05):** 比 #1 多一步 **guard 抽纯 helper**——把 `runOneStep.guardSurfaceTemperature` 的 body 逐字搬进纯函数 `vegetation.energy.guardSurfaceTemperature(scopeStep, prevTss, kt)`(无 model/state 读写,Julia-friendly `(state,prev,t)`),local 缩成薄 wrapper。owner 搬:`buildInitialState` 生 `Y.surface.Tss=NaN`(哨兵→kt=1 fallback 用 Ta);`prepareRuntime` 删 `runtime.thermal`;wrapper 的 fallback 源翻到 `built.Y.surface.Tss`(**正常路径 `Tss=th.Tsave` 不变**,只动 reject fallback);`emptyStepPending`/③-block/`commitStepPending`(step-5c)/`assertRequiredPending`(required,**NaN 合法**)全搬 `pending.prognostic.surface`;`materializeCheckpointStateDefaults` 加 surface clause(旧 ckpt 无 rt.thermal→默认 NaN)。**关键发现**:`runtime.thermal` 历史上**从没进 checkpoint**(Tss resume 后丢成 NaN),搬到 `Y.surface` 随 `ckpt.Y` 顺带**修了这个漏存 gap**——非 resume bit-exact,resume 行为改善(用户批准 Option A)。gate 11/11 绿(含新 `test_guardSurfaceTemperature` 12 检查:reject/accept 全测)。

---

## 10. 测试门槛 + 两个新测试(§10 用户要求)

**每次搬字段必跑**:`test_runOneStep` · `test_compassCheckpoint` · `test_resume` · `test_luna_chain` ·
`test_methodB_water_stress` · `test_stepPlantHydraulicsV2` · `test_computeCanopyWaterEvap` ·
`test_stepCanopyWaterSnowEvap_failstep` · `test_stepPlantCN_golden` · `test_stepPlantCN_failstep` ·
`test_assembleNitrogenVariables` · `test_rtmo_parity` · `test_sif_step`。

**新增 2 测试(2026-07-05 已建,随增量 #1 上线)**:
1. ✅ **`tests/test_ownership_allowlist.m`**:断言 `~isfield(runtime,'snow')` 且 `Y.snow.SWE==0`(runtime 服务层不再出现已搬的物理状态字段;后续增量按同款逐字段扩断言)。
2. ✅ **`tests/test_checkpoint_compat.m`**:造一个旧路径 checkpoint(`ckpt.Y` 无 snow + `ckpt.runtime.snow.SWE=1.9`),resume 后断言 `Y.snow.SWE==1.9`(`materializeCheckpointStateDefaults` 映射成功)。

---

## 结论 / Verdict

**现在做**(本轮,零物理风险):**本审计表**(边界钉死)。✅ **极低风险第一搬已完成** = `snow.SWE → Y.snow.SWE`(单标量,带旧 checkpoint 映射 + 13 测试 + 2 新测试,gate 9/9 绿,2026-07-05)——这是**状态归属迁移模板**,后续 5 个字段照此逐个搬(待用户逐个放行)。

**随 ④**:`Theta_II_prev`(冰时间层)、`DSTOR0` 积水维护、滞回反向点(PR-F ③.10)—— 消费 ④ 冻融/水解产物,一体落地。

**判断**:不宜现在继续堆 ④。**先按本表把 state/aux/cache/runtime-service/diagnostics 边界钉死**,再逐个(风险序)搬迁 + 验证,最后 Julia shadow port 验分类。边界定死后进土壤循环,转 Julia 才不二次大拆。
