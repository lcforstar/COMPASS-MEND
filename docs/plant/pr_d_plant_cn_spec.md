# PR-D spec -- assemble/scatter bridge for the plant C/N + root-N-uptake block

> **MERGED into the migration ledger.** The PR-D migration record (what/why/how + golden results) now
> lives in [compass_complete_architecture.md](../architecture/compass_complete_architecture.md) §13 (the single source of
> truth for the COMPASS.m -> runOneStep migration). This file is kept as the DETAILED APPENDIX: the
> per-field assemble/scatter maps (§1/§2), the 2-slot decision (§0), and the forcing-index contract (§4.1).
> COMPLETE (all steps DONE 2026-07-04).
>
> **[owner 更新 2026-07-05 增量#5]** 本文档下方所有 `runtime.phenology` 说法是**当时**的 owner;phenology
> 已搬入 `Y.vegetation.phenology`(state-ownership 迁移增量#5)。phen_alloc 现读 LAGGED
> `model.built.Y.vegetation.phenology`;env_turn 现读 fresh pending phenNext(1 步 lag 不变)。历史文字保留脉络,当前 owner 以 `Y.vegetation.phenology` 为准。

> Migrate legacy `COMPASS.m:1195-1421` (root N uptake + plant C/N + NEE) into the new container
> chain. Same pattern as the ebal migration's `assemblePlantHydro.m`: the 4 physics functions
> (`runPlantNDemand`, `applyRootUptake`, `runPlantCNStep`, `computeNEE`) are already in the package
> tree and reused verbatim; the new chain lacks a per-step CALLER and a `NitrogenVariables`
> assemble/scatter bridge. Source of this spec: the `map-plant-cn-assemble-scatter` workflow
> (7 agents, per-field maps). Date: 2026-07-04.

## 0. Architecture decision (resolve before wiring)

Legacy `NitrogenVariables` plant pools are `tS x 1` time series, read at `idx_prev=max(1,KT-1)` and
written at `KT`. New-chain `Y.plant_cn.{C,N,transfer}.*` are SCALARS (`buildInitialPlantPools`). The
4 functions do `pool(KT-1)` reads / `pool(KT)` writes with the REAL KT -- a naive assemble of a
2-element array breaks (`pool(17000)` on length-2). Note the asymmetry: the DIAGNOSTICS
(`plantN_diag`, `nsc_diag`) are already `tS x 1` struct arrays (`allocateWorkspace.m:199-200`), so
`nsc_diag(KT)` writes at real KT work; only the prognostic POOLS are scalar.

**Chosen: 2-slot local remap (keeps the 4 functions verbatim -> golden-safe).**
- `assembleNitrogenVariables(built, kt)` builds pool fields as 2-element columns `[prev; prev]`
  (index 1 = accepted scalar from `Y.plant_cn`, index 2 = write slot).
- Call the 4 functions with a LOCAL `kt2=2` (so `idx_prev=max(1,1)=1` reads the accepted scalar,
  writes land at index 2). Forcing that legacy indexes by real KT (`MeteoData.Ta(KT)`,
  `current_date(KT)`, `Precip`) is passed as already-resolved values by the glue, NOT re-indexed
  inside -> real KT and pool KT decouple.
- `nsc_diag`/`plantN_diag`/`turnover_diag`/`frPHU_diag`: **in-flight 2-slot diag scratch** (written at
  kt2=2 by the reused functions), mapped back to the REAL-kt `built.diagnostics(kt)` row by scatter.
  (Implemented this way, NOT full tS arrays passed in, because runPlantCNStep uses ONE kt for pools AND
  diag rows; see assembleNitrogenVariables doc.)
- `scatterNitrogenVariables` writes index-2 pool values back to `Y.plant_cn` scalars on accept.

**Commit discipline (mirrors phenNext/laiStateNext):** all 4 runStep functions produce PENDING
output; `Y.plant_cn`, accumulators, flags, `flux.Srt`, `flux.SinkS`, diag rows commit ONLY on
accept. A failed outer step must not advance pools or fire one-shot events.

## 1. assembleNitrogenVariables -- READ fields, grouped by container

### A. parameters.soil_n (static, read-only, never scatter back)
- `RHO_bulk` -> `parameters.soil_n.RHO_bulk` (buildNitrogenSoilParams.m:32), NL x 1 [g/cm3]
- `RHOKG` -> `parameters.soil_n.RHOKG` (:33), NL x 1 [-]
- `DIF1111` -> `parameters.soil_n.DIF1111` (:35-37), **1 x NL x nS** (singleton lead dim -- trap)

### B. soil state for SoilStateN
- `SoilVariables.Theta_LL` -> `workspace.derived.Theta_LL` (allocateWorkspace.m:110, per step), NL x 2 [-]
- `SoilVariables.Theta_s`  -> **`parameters.soil.theta_s`** (NOT soil_n -- trap), NL x 1 [-]
- `SoilVariables.TT`       -> reconstruct from **`Y.T`** (bottom-up node temp) via `init.assembleSoilStateN` -- **DONE (step 2): `TT = Y.T(:)`, NO flip, KT=1 golden**, NN x 1 [degC]
- `SoilVariables.C_NH4`    -> `Y.nh4` (node, bottom-up), NN x 1 [ug/cm3 liq]
- `SoilVariables.C_NO3`    -> `Y.no3` (node, bottom-up), NN x 1 [ug/cm3 liq]

### C. Y.plant_cn prognostic pools (READ prev; sec 2 scatters back)
`Y.plant_cn.<grp>.<organ>` (buildInitialPlantPools.m:75-86). Name remap (drop `C_`/`N_`/`_xfer`):
- C pools -> `Y.plant_cn.C.{labile,leaf,root_fine,livestem,deadstem,livecroot,deadcroot,fruit,root_storage,stem_storage,xsmr}`
- N pools -> `Y.plant_cn.N.{... same minus xsmr}`
- transfer -> `C_leaf_xfer->transfer.C.leaf`, `C_root_xfer->transfer.C.root`, `C_livestem_xfer->transfer.C.livestem` (+ N)
- baseline (read-only) -> `Y.plant_cn.baseline.{C_pool_init_const,N_pool_init_const}`

### D. Persistent accumulators + flags -- **DONE in step 1a** (`Y.plant_cn.accum/state`)
- 15 accumulators incl. `cum_U_N` -> `Y.plant_cn.accum.*` (init 0)
- 6 flags/scalars (`turnover_state,onset_mobilized_flag,resp_night_hours,C_deadstem_soy,N_deadstem_soy,pruned_this_year`) -> `Y.plant_cn.state.*`

### E. built.diagnostics (READ: only `nsc_diag(kt).R_a_actual_kg`, by computeNEE)
`built.diagnostics.nsc_diag` (allocateWorkspace.m:200). But NEE must read the row from the PENDING
runStepPlantCN output (not yet committed) -- see sec 4.

### F. non-NitrogenVariables intermediates (assemble from other containers -- sec 3)
`plantHydro.*` -> `built.p.plantHydro` (assemblePlantHydro, cached at init); `domain.{NL,NN,dz,tot_depth,sumdz}`;
`Delt_t` -> `control.time.Delt_t`; switches `lRwuef/lNefc` -> `control.switches`; `phen_state/phen_now` ->
`runtime.phenology`; `canopy.LAI` -> `workspace.scope_step.canopy.LAI`; forcing -> `built.forcing.hourly.*`;
`current_date` -> `runtime.calendar`.

## 2. scatterNitrogenVariables -- WRITE fields, committed ONLY on accept

### Prognostic (carry -> Y / workspace.flux)
- all C/N/transfer pools (index-2) -> `Y.plant_cn.{C,N,transfer}.*` (last write of the step)
- `C_xsmr` -> `Y.plant_cn.C.xsmr` (may be < 0 = debt)
- `cum_U_N` (glue, COMPASS.m:1280) -> `Y.plant_cn.accum.cum_U_N` -- plant-N conservation inflow; MUST be prognostic
- 14 other `cum_*` -> `Y.plant_cn.accum.*`
- `resp_night_hours,onset_mobilized_flag,pruned_this_year,C_deadstem_soy,N_deadstem_soy,turnover_state` -> `Y.plant_cn.state.*`
- `SinkS` (glue, :1272) -> `workspace.flux.SinkS` (allocateWorkspace.m:159, NL x nD x nS) -- **load-bearing: drives cvSinkSnh4/no3; 2026-05-13 bug**
- `Srt` (refined by applyRootUptake, :1265) -> `workspace.flux.Srt` (:157, NL x 2) -- **overwrites the PR-C.b scope_step.Srt; the refined one is committed**

### Diagnostic (scatter maps the in-flight 2-slot scratch -> real-kt built.diagnostics row)
The reused functions write diag rows at kt2=2 (2-slot scratch seeded by assembleNitrogenVariables);
scatter copies slot 2 -> `built.diagnostics.<name>(kt)`. (NOT full tS arrays passed in -- runPlantCNStep
uses one kt for pools AND diag; see §0 / assembleNitrogenVariables doc.)
- `nsc_diag(2)` -> `built.diagnostics.nsc_diag(kt)` (has `R_a_actual_kg`, consumed by NEE same step). **DONE.
  NOTE: nscDiagTemplate + setInitialNitrogen were missing `NPP_actual_kg` (applyNSCCoupling.m:304 writes it);
  the whole-struct `nsc_diag(KT)=nsc_diag` then hit MATLAB:heterogeneousStrucAssignment. Fixed 2026-07-04.**
- `plantN_diag(2)` -> `built.diagnostics.plantN_diag(kt)` -- `A_croot` added to plantNDiagTemplate (step 1b). **DONE.**
- `turnover_diag(2)` -> `built.diagnostics.turnover_diag(kt)` (step 1b). **DONE.**
- `frPHU_diag(2)` -> `built.diagnostics.frPHU_diag(kt)` (step 1b). **DONE.**
- `demand_NH4_KT,demand_NO3_KT` (glue) -> `built.diagnostics.root_demand.{nh4,no3}(kt)` (step 1b, guarded). **DONE.**
- `root_uptake_M1_diag` (glue) -> `built.diagnostics.root_uptake_M1_diag` (step 1b, guarded). **DONE.**

computeNEE writes nothing to NitrogenVariables -- NEE/Resp -> `scope_step.fluxes.{NEE,Resp}` (transient,
committed accept-only; rides the scope_step snapshot, so a failed step keeps neither).

## 3. Intermediate assembly (build inside the runStep functions)

- **RLD** (COMPASS.m:1211-1224): preferred `plantHydro.rootProfile.root_length_density_m_m3`
  (-> `built.p.plantHydro.rootProfile`); fallback `RAI(1:NL) ./ (2*pi*r_root_m*(DeltZ*1e-2))`
  (RAI truncate to 1:NL, r_root_m -> `built.p.plantHydro.r_root_m`, DeltZ -> `domain.dz`).
- **SoilStateN** (:1204-1209): `{Theta_LL->derived, Theta_s->parameters.soil, RHO_bulk/RHOKG/DIF1111->soil_n, TT->from Y.T}`.
  Frame trap: `Y.nh4/no3/T` bottom-up node vs `derived.Theta_LL` working frame vs `theta_s` top-down.
- **phen_alloc** (:1255-1256, feeds runPlantNDemand) -- from the LAGGED committed phenology
  (`runtime.phenology.state`, NOT pending phenNext): `gdd_budburst=state.gdd5-state.gdd5_at_budburst`
  (NaN pre-budburst); `DS_fruit=state.DS_fruit` (0 if NaN); `phen=runtime.phenology.now`;
  `harvested=double(logical(state.harvested_flag))`.
- **env_turn** (:1374-1382, feeds runPlantCNStep) -- from the FRESH pending phenology (phenNext):
  `frPHU=phen_HUI(day_idx)` with **`day_idx=ceil(kt/steps_per_day)`, `steps_per_day=round(86400/Delt_t)`**
  (day-index trap, not kt); `phen=phenNext.now`; `phen_state=phenNext.state` (whole object);
  `doy_now`->runtime.calendar; `LAI=scope_step.canopy.LAI`; `DLAI`->veg.phenology.DLAI (0.45);
  `use_leafc_sla=strcmpi(veg.canopy.LAI_engine,'leafc_sla_phen')`; `year`->runtime.calendar.current_year.
- **U_N_total_g_m2_s** (:1277): `(SPUptake(1)+SPUptake(2))*1e-2` -- keep the 1e-2 in the glue.
- **Phenology already migrated** as `advanceDailyPhenology` (runOneStep.m:117). env_turn reads pending
  `phenNext`; phen_alloc reads committed `runtime.phenology` (lagged). Do NOT re-run dispatch.

## 4. Wiring order in runOneStep (3) block

Insert AFTER `phenNext = advanceDailyPhenology(...)` (runOneStep.m:117), pre-substep. Gate: outer
`if switches.lRwuef==1`, inner `if switches.lNefc==1` (else water-only applyRootUptake + degenerate NEE).

```
NV = init.assembleNitrogenVariables(model.built, kt);       % gather sec-1 reads (2-slot pools)
phenLag = model.runtime.phenology;                          % LAGGED for phen_alloc
[demand_ctx, root_ctx, NV] = runStepNDemand(model, NV, kt, phenLag);
[NV, sinkS, srtRefined, U_N_total] = runStepRootUptake(model, NV, root_ctx, kt);  % OVERWRITES scope_step.Srt
[NV, R_a_plant_kg] = runStepPlantCN(model, NV, demand_ctx, U_N_total, phenNext, kt);
neeOut = runStepNEE(model, NV, kt);                         % pending workspace.fluxes.{NEE,Resp}
plantCnNext = NV;                                           % PENDING; scattered on accept
```

**Accept branch** (after `model.runtime.phenology = phenNext;`):
```
model.built = init.scatterNitrogenVariables(model.built, plantCnNext, kt);
model.built.workspace.fluxes = neeOut;
```

- **Srt interaction:** PR-C.b `applyRwuToSrt` staged `scope_step.Srt` from ebal RWU. `runStepRootUptake`
  reads it, `applyRootUptake` refines it, and the REFINED Srt is committed to `flux.Srt` (supersedes PR-C.b).
- **NEE ordering:** `nsc_diag(kt).R_a_actual_kg` is produced by runStepPlantCN and consumed by
  runStepNEE the same step. Read it from the PENDING `NV.nsc_diag(kt)`, NOT `built.diagnostics`.
  Use `nsc_diag(kt).R_a_actual_kg` (actual), NOT `R_a_plant_kg` (potential -- double-counts ~272 gC/m2/yr).

### 4.1 Forcing-index contract -- kt2 vs real kt (**the step-4 trap**)

The reused functions are called with a LOCAL `kt2=2` so their pool reads (`idx_prev=max(1,kt2-1)=1`)
hit the accepted scalar. But some of them ALSO index REAL forcing by the SAME argument. Audit of every
`KT`-driven read across the 4 functions:

| function | real-forcing / real-kt read | fix |
|---|---|---|
| `runPlantNDemand` | **`MeteoData.Ta(KT)`** (line 135) -- the ONLY in-function forcing-index | pass a **2-slot meteo view** so `localMeteo.Ta(2) = built.forcing.hourly.Ta(realkt)`; `current_date_KT` is already a pre-resolved SCALAR (glue passes `current_date(KT)`), so the wrapper resolves it once = `start_time + seconds(readStepTime + dtOuterActual)` (SAME ts as `advanceDailyPhenology`, keeps demand DOY == phenology DOY) |
| `runPlantCNStep` | none -- `KT` indexes ONLY pools + diag rows | safe to call with `kt2=2` |
| `applyRootUptake` | none -- no `KT` argument at all | call as-is (glue-resolved args) |
| `computeNEE` | none -- takes resolved scalars `(GPP,R_h,R_a_plant)` | glue computes them at real kt |

Second real-kt read lives in the **NEE glue** (COMPASS.m:1396-1398), NOT inside `computeNEE`:
`R_h = equations.soil_respiration(meteo.Ta, P30_cm)` with `P30_cm = sum(MeteoData.Precip(i0:KT))*Delt_t`,
`i0=max(1,KT-719)` -- a 30-day precip window on REAL kt + current-step Ta. The `runStepNEE` wrapper must
build this with the REAL kt (`built.forcing.hourly.{Ta,Precip}`), never `kt2`. `computeNEE` itself stays a
pure scalar function -> no kt2 confusion.

`built.forcing.hourly.Ta(kt) == built.workspace.scope_step.meteo.Ta` (prepareScopeStep.m:258 direct slice,
no transform) == legacy `MeteoData.Ta(kt)`. The 2-slot meteo view is built by
`vegetation.integrated.buildDemandMeteoView(built, kt)` (pinned by `tests/test_ndemand_forcing_view.m`).
**HARD RULE for step 4: never call `runPlantNDemand(2, <real MeteoData>, ...)` directly.**

## 5. Golden strategy -- KT=1 component parity per function

Reuse each legacy function verbatim on assembled inputs; assert bit-parity vs a legacy `COMPASS.m`
run at the same KT. The stub soil solve freezes `Theta_LL/T`, so multi-step trajectories diverge ->
parity is defined on **KT=1** (assembled inputs identical to legacy: `idx_prev=max(1,0)=1`=init).

| Function | compare |
|---|---|
| runPlantNDemand | `demand_ctx`, `root_uptake_ctx` field-by-field |
| applyRootUptake | `Srt, SinkS, SPUptake, uptakeDiag` (tolerate empty uptakeDiag when Nefc!=1) |
| runPlantCNStep  | updated NV pools + `plantN_diag(kt)`,`nsc_diag(kt)` + `R_a_plant_kg` |
| computeNEE      | `fluxes.NEE` scalar |

## 6. BLOCKERS + traps

**BLOCKERS (resolve before PR-D runs):**
1. Container for accumulators/flags -- **DONE (step 1a: Y.plant_cn.accum/state)**.
2. `TT` reconstruction from `Y.T` (bottom-up node) -- **DONE (step 2)**: `init.assembleSoilStateN`, `TT = Y.T(:)` no flip, pinned by `tests/test_assembleSoilStateN.m` (6/0, KT=1).
3. Diagnostics gaps: add `built.diagnostics.{turnover_diag,frPHU_diag,root_demand,root_uptake_M1_diag}` +
   `A_croot` to plantNDiagTemplate -- step 1b.
4. tS x 1 vs scalar -- the 2-slot remap (sec 0).

**TOP TRAPS:**
- `Theta_s` is `parameters.soil.theta_s`, NOT soil_n (soil_n has `Theta_f` only).
- `DIF1111` is 1 x NL x nS (singleton lead dim).
- `SinkS` scatter is load-bearing -> `workspace.flux.SinkS`; skipping = cvSinkS*=0 bug.
- `applyRootUptake` OVERWRITES the PR-C.b Srt -- commit the refined one.
- `SPUptake x 1e-2` unit conversion in the glue.
- `day_idx = ceil(kt/steps_per_day)` for `phen_HUI` (not kt).
- 1-step phenology lag: phen_alloc -> committed runtime.phenology (lagged); env_turn -> pending phenNext (fresh).
- NEE uses `nsc_diag(kt).R_a_actual_kg` (actual), never `R_a_plant_kg` (potential).
- `RL`/`Lmax` passed but unused inside applyRootUptake (keep for signature compat).
- frame flip: `Y.nh4/no3/T` bottom-up vs `derived.Theta_LL` working frame vs `theta_s` top-down.
- `LAI_engine` default mismatch (`epic_phen_mix` legacy fallback vs `leafc_sla_phen` buildVeg default) --
  verify config sets it explicitly (changwu: leafc_sla_phen).
- `fluxes.Resp` written in umol but computeNEE consumes kg.
- all 4 runStep outputs commit only on accept.

**Step order:** 1a accum/state (DONE) -> 1b diagnostics slots (DONE) -> 2 TT reconstructor (DONE) ->
3 assemble/scatter (**DONE**: `init.assembleNitrogenVariables` + `init.scatterNitrogenVariables`, 2-slot
pools + 2-slot diag scratch, scatter accept-only; `tests/test_assembleNitrogenVariables.m` 41/0) ->
4 wire (**DONE**: `vegetation.integrated.stepPlantCN` orchestrator + `vegetation.integrated.buildDemandMeteoView`
+ runOneStep (3)-block `stepPlantCarbonNitrogen` [lRwuef gate] + accept-branch scatter; `test_runOneStep`
103/0. En route fixed the nsc_diag NPP_actual_kg template drift, §2) ->
5 golden (**DONE**: `tests/test_stepPlantCN_golden.m` 17/0 = per-function KT=1 parity vs an independent
legacy glue [demand_ctx/pftP, applyRootUptake Srt/SinkS/SPUptake bit-exact, runPlantCNStep pools/nsc_diag/
accum/R_a, computeNEE]; `tests/test_stepPlantCN_failstep.m` 7/0 = accept-only commit discipline).

**PR-D COMPLETE.** All hold/uncommitted.
