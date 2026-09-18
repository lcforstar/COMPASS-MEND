# Plan: Soil Interface Temperature Outer Closure (Nested Scheme C)

**Date:** 2026-07-19  
**Audience:** Codex / implementers  
**Status:** Design only — implement after reading this plan end-to-end  
**Editor context:** Changwu thermal BC; KT7464 numerical recovery via pure Dirichlet; energy closure still fails  

---

## 0. One-sentence goal

Keep a **stable Dirichlet soil top** for the water–heat–gas Picard, but add a **scalar outer loop on interface temperature \(T^*\)** so that

\[
G_{\mathrm{soil}}(T^*) = G_{\mathrm{sfc}}(T^*)
\]

**without re-running ebal** and **without putting aggregate Robin into the \(h/U/P_g\) matrix**.

---

## 1. Background (do not re-litigate)

### 1.1 What already works (keep)

| Item | Evidence / location |
|------|---------------------|
| `surface_temperature` is **strict Dirichlet** (no silent upgrade to Robin when `top_heat_response` exists) | `src/+soil/+thermal/+boundary/materializeRequest.m` |
| Changwu production kind | `data/sites/changwu/site.toml`: `heat_top_kind = "surface_temperature"` |
| KT7464 fair replay | `dt = 3600 s`, 13 iters, one attempt |
| KT1440→1960 incl. meltouts | snow water/energy residuals ~ roundoff |
| Meltout state machine, single snow advance, ledgers | **do not touch** unless a test proves a direct conflict |

### 1.2 What fails (pure Scheme A)

Pinning ebal `Tss` once:

- Soil **reaction heat flux** ≠ surface **available energy** from ebal/snow budget.
- Fixture gaps ~70–140 W m⁻²; production trajectory peak interface mismatch ~**378 W m⁻²**; large step energy residual.
- **Cannot** ship as production by loosening thresholds or rewriting ledgers.
- Year-long run **not** started on pure A for that reason.

### 1.3 What is forbidden

1. Re-run `stageCanopyEnergy` / `ebal` inside the \(T^*\) loop (except existing **event** restage paths for meltout/large cover change).
2. Aggregate partial-snow Robin  
   \(q(T_N)=(1-f)q_b(T_N)+f G_{si}(T_N)\)  
   inside the coupled soil Picard (KT7464 collapse).
3. Silent remap of `surface_temperature` → `surface_energy_response` in `materializeRequest`.
4. “Closing” energy only in diagnostics (`buildSurfaceSoilEnergyBalance`) without changing \(T^*\).

---

## 2. Classification

| Name | Role |
|------|------|
| **Inner** | Existing soil Picard \(h/U/P_g\) with **Dirichlet** \(T_N = T_{\mathrm{trial}}\) (Scheme A kernel) |
| **Outer** | Scalar root find on \(T^*\) so fluxes match (Scheme **C** target, nested implementation) |
| **Not** | Classical one-shot Neumann (Scheme B); not matrix Robin |

Call it in code/docs: **`soil_interface_temperature_closure`** (nested Dirichlet + flux match).

---

## 3. Physics (frozen surface maps)

After §3 of `runOneStep` (ebal + postEbal + snow advance), **freeze** surface maps once per outer atmospheric step (or per soil substep start if substepping—see §5).

### 3.1 Target flux \(G_{\mathrm{sfc}}(T)\) [W m⁻²], into soil positive

Use **affine evaluations only** (no ebal):

**Bare contribution** (from staged SCOPE fluxes / `top_heat_response` fields already on workspace):

\[
G_b(T) = (1-f)\Big[
  q_b + q_b'(T - T_{b,\mathrm{ref}})
\Big]
\]

where \(q_b\), \(q_b'\), \(T_{b,\mathrm{ref}}\) come from  
`bare_ground_external_energy_flux_*` and linearization temperature  
(or equivalent fields already assembled into `workspace.surface.top_heat_response`).

**Snow contribution** (from `workspace.snow.ground_heat_response`):

\[
G_{si}(T) = f \cdot C_{si}\,(T_{\mathrm{snow}} - T)
\]

plus any **temperature-independent** redistribution / meltout deposit terms already defined in the snow–soil contract (add as constant offset, not re-derived in ebal).

**Grid target:**

\[
G_{\mathrm{sfc}}(T) = G_b(T) + G_{si}(T) + G_{\mathrm{redistrib,const}}
\]

with \(f\), \(C_{si}\), \(T_{\mathrm{snow}}\) frozen from the staged snow transaction.

### 3.2 Soil reaction flux \(G_{\mathrm{soil}}(T)\)

For trial Dirichlet \(T_N = T\):

1. Materialize heat top as `surface_temperature` with `temperature_C = T`  
   (and `surface_skin_temperature_C` / drivers consistent so materialize does not fall back incorrectly).
2. Run the **existing** soil attempt / heat solve path.
3. Read top heat flux into domain from the authoritative boundary result  
   (same quantity used by energy diagnostics / `boundary_flux_summary` /  
   `top_response_evaluated_flux` path—**one** definition, SI W m⁻²).

### 3.3 Residual

\[
F(T) = G_{\mathrm{soil}}(T) - G_{\mathrm{sfc}}(T)
\]

Solve \(F(T^*) = 0\) (or \(|F| < \varepsilon_G\)).

**Initial guess:** guarded `Tss` from ebal (`thermal.Tss`).

**Suggested \(\varepsilon_G\):** start with `max(5, 0.01 * max(1, |G_sfc|))` W m⁻² configurable under numerics; report final \(|F|\) always.

---

## 4. Where to implement

### 4.1 New modules (preferred)

| File | Responsibility |
|------|----------------|
| `src/+surface/+energy/evaluateFrozenSurfaceSoilHeat.m` | Given frozen workspace maps + scalar \(T\), return \(G_{\mathrm{sfc}}\) [W m⁻²] and breakdown (bare/snow/const). **No ebal.** |
| `src/+surface/+energy/closeSoilInterfaceTemperature.m` | Scalar root (secant / guarded Newton). Calls soil flux evaluator callback. Returns \(T^*\), \(G_{\mathrm{soil}}\), \(G_{\mathrm{sfc}}\), iters, converged flag. |
| Optional helper | `readSoilTopHeatFlux_W_m2(built, soilAdvanceOrBoundary)` — single SI conversion point |

Editor signature per coding standard on new/changed model functions.

### 4.2 Call site

**Primary (recommended for v1):**  
Inside `solver.advanceSoilEquationOuterStep` (or a thin wrapper it calls), **once per accepted soil outer interval / first trial of each substep window**, **before** locking the heat BC for that interval:

```text
T* = closeSoilInterfaceTemperature(...)
→ write T* into heat drivers / surface_skin override for materialize
→ then existing materialize + attemptSolve loop for that dt
```

**Acceptable simpler v1:**  
In `runOneStep` §4, **once per atmospheric hour** before `advanceSoilEquationOuterStep`:

- Close \(T^*\) with a representative `dt` (outer step length)  
- Override staged skin temperature used by heat materialize to \(T^*\)  
- Then advance soil as today  

Document that substep-level closure is a follow-up if hour-level residual remains large under adaptive dt.

### 4.3 Do **not** put outer loop in

- `ebal.m` / `stageCanopyEnergy` (except existing event restage)
- `buildSurfaceSoilEnergyBalance.m` (diagnostics only)
- `materializeRequest` as Robin

### 4.4 Interaction with `reconcileSurfaceSoilHeatTransaction`

| Module | Role after this plan |
|--------|----------------------|
| Outer \(T^*\) closure | Match **soil reaction** to **frozen surface available energy** |
| `reconcileSnowGroundHeat` | Replay snow at **accepted** soil interface \(T\); keep meltout/event restage rules |
| Event restage ebal | **Only** cover/domain/meltout events (existing); not used as routine energy fixer |

Order after soil accept:

1. Soil advanced with \(T^*\) (or closed per substep)  
2. Existing snow ground reconciliation  
3. Diagnostics `buildSurfaceSoilEnergyBalance` **must improve** vs pure A  

---

## 5. Algorithm (outer loop)

```text
T0 = thermal.Tss   % guarded
Freeze maps from workspace (f, qb, qb', Tbref, Csi, Tsnow, redistrib)

if |F(T0)| < eps_G
  T* = T0; return   % skip outer work (common bare / mild cases)
end

Secant (or hybrid):
  T1 = T0 + small perturbation or one Picard heat-only probe
  for k = 1:max_outer (e.g. 6)
    evaluate G_soil(Tk) via Dirichlet soil heat (see cost policy)
    evaluate G_sfc(Tk) from frozen maps
    if |F(Tk)| < eps_G: success
    update Tk by secant; clamp to physical band e.g. [-60, 90] °C
  end
  on failure: keep best |F| trial, set flag, do not silent-success
```

### 5.1 Cost policy (required)

1. **Skip outer** when `|F(Tss)| < eps_G` (or snow inactive and bare residual small).  
2. Cap `max_outer` (default **6**).  
3. Prefer **secant** (2 flux evaluations minimum).  
4. v1 may re-run full `attemptSolveSoilEquations` per trial (simplest, slower).  
5. Optional v1.1: freeze water/gas state, re-solve **heat-only** for outer trials (faster)—only if existing APIs allow without refactor explosion.  
6. Log: `outer_iters`, `|F|` initial/final, \(Tss\), \(T^*\), \(G_{\mathrm{soil}}\), \(G_{\mathrm{sfc}}\).

**Do not** re-enter ebal in this loop.

---

## 6. Driver / materialize wiring

For each outer trial and for the final soil advance:

1. Keep `heat_top_kind = surface_temperature` in site config.  
2. Ensure heat drivers pass  
   `surface_temperature_C = T_trial`  
   (today often from `scope_step.thermal.Tss` in `materializeSoilBoundaryTrial`—either update a **trial override** field on workspace.surface, e.g.  
   `workspace.surface.soil_interface_temperature_C = T_trial`,  
   read with priority over Tss in materializeSoilBoundaryTrial,  
   or temporarily set the skin used by materialize).  
3. **Do not** attach Robin in `surface_temperature` branch.  
4. After success, publish \(T^*\) to diagnostics and keep snow reconcile using the **accepted** soil surface temperature (time-mean if that is already the project convention).

---

## 7. Configuration

Add optional numerics (names indicative):

```toml
# e.g. under numerics or heat section
soil_interface_closure_enable = true
soil_interface_closure_tol_W_m2 = 5.0
soil_interface_closure_max_iters = 6
soil_interface_closure_skip_if_abs_F_below_W_m2 = 5.0
```

Default: **enabled** for Changwu when heat + snow paths active; allow disable for pure A A/B experiments.

---

## 8. Acceptance tests (must pass before year run)

### 8.1 Regression (must not break)

| Test | Expect |
|------|--------|
| Unit: materialize `surface_temperature` never becomes response | existing / extend |
| Snow coupling unit tests | pass |
| KT1951-style meltout | pass; single integrate; correct fail reasons |
| KT7464 fair step | still ~3600 s class convergence (outer iters extra OK if not collapsing dt) |

### 8.2 Energy (why this plan exists)

| Check | Pure A (baseline) | Target with closure |
|-------|-------------------|---------------------|
| Full-snow fixture interface | ~136 W m⁻² open | within \(\varepsilon_G\) (or ≪ 136) |
| Strong meltout fixture | ~71 W m⁻² open | within \(\varepsilon_G\) (or ≪ 71) |
| Production short trajectory max interface mismatch | ~378 W m⁻² | much reduced; document max |
| Step total energy residual | large | within agreed ledger tolerance |

Use the **same** diagnostic definitions that failed pure A (`buildSurfaceSoilEnergyBalance` / fixture harness)—do not redefine success by softening metrics.

### 8.3 Replay order

1. Fixture clamps (full snow, meltout)  
2. KT7464 step (+ optional short prepare)  
3. KT1440→1960 or KT1951 window  
4. Only then consider year run  

---

## 9. Out of scope (this PR / this task)

- Full CLM snow layers in soil thermal matrix  
- Dual soil columns for partial cover  
- Replacing ebal skin solver  
- Re-enabling production `surface_energy_response` aggregate Robin for partial cover  
- “Fix” by ledger rescaling without \(T^*\)  

Follow-up (later): resistance-based \(T_{\mathrm{interface}}\) inside `Tsave` formula; heat-only outer solves.

---

## 10. Implementation checklist for Codex

1. Read this plan + `materializeRequest.m` surface_temperature branch + `materializeSoilBoundaryTrial.m` + `runOneStep` §3–4 + `reconcileSurfaceSoilHeatTransaction`.  
2. Invoke **compass-coding-standard** before editing model code under `src/`.  
3. Implement `evaluateFrozenSurfaceSoilHeat` + `closeSoilInterfaceTemperature`.  
4. Wire trial \(T\) override into heat materialize path.  
5. Call outer closure from agreed site (§4.2).  
6. Diagnostics/logging fields for outer loop.  
7. Config flags + defaults.  
8. Run §8 tests; paste numbers vs pure A.  
9. **Do not commit / do not start year run** until energy fixtures pass and user approves.  
10. Keep meltout / single-advance snow logic unchanged unless tests force a minimal fix (document why).

---

## 11. Success definition

- Numerical: no return of KT7464-style subsecond collapse under partial snow.  
- Physical: interface \(G_{\mathrm{soil}} \approx G_{\mathrm{sfc}}\) under the same audit that rejected pure A.  
- Architectural: ebal once per normal hour; outer loop is **scalar \(T^*\) only**; soil Picard remains Dirichlet-top.  

---

## 12. Short message to paste if Codex only gets a blurb

> Implement nested soil-interface temperature closure: outer scalar \(T^*\) so \(G_soil(T^*)=G_sfc(T^*)\) using frozen post-ebal affine bare flux + snow conductance only (no re-ebal). Inner soil Picard keeps strict Dirichlet (already fixed in materializeRequest; changwu uses surface_temperature). Do not put aggregate Robin back into the Picard matrix. Call from soil advance / runOneStep §4. Skip outer when residual at Tss is small; max ~6 secant iters. Prove energy fixtures that failed pure A; keep meltout tests green; no year run until energy passes.
