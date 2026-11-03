# runOneStep read/write contract (proposal, 2026-06-25)

Priority #1 from `docs/claude_initialization_handoff.md` §8: define the read/write
boundary of `[model, step_report] = runOneStep(model, kt)` BEFORE writing it.
Grounded in a read-only audit of the legacy `src/COMPASS.m` main loop (line anchors below),
not guessed.

## Scope (proposed)

- `runOneStep` advances **one outer forcing interval** (`Delt_t_outer`, one `kt`). It OWNS the
  inner adaptive-substep loop and the Picard loop; the caller just loops `kt = 1..tS`.
- **First kernel = soil only** (water + heat + gas + nitrogen), per handoff §8.2. Canopy /
  ebal / plant-hydraulics are NOT migrated yet; their two outputs the soil step needs
  (`Srt`, `Tss`) are **injected as step-frozen inputs** (see Coupling).

## The coupling finding (central constraint)

The soil solve cannot run from `{Y, forcing, boundary, parameters}` alone. Legacy computes
two coupling fields ONCE before the Picard loop and holds them constant across it:

| field | legacy origin | consumed by | COMPASS.m |
|---|---|---|---|
| `Srt` (root water uptake sink) | ebal/plant hydraulics (RWU) | Richards RHS sink | produced :1054, used :1731 |
| `Tss` (surface skin temperature) | ebal | heat top BC (`BCT`) | produced :1821, used by thermal solve |
| top water BC (Throughfall/EVAP/infil-pot) | canopy water + ebal | Richards top BC | used :1731 |

→ The contract must expose a **step-coupling input** the soil kernel reads but never writes.
Proposed home: `model.runtime.coupling.{Srt, Tss, throughfall, evap, ...}` (one writer =
the future canopy/ebal step; for now seeded by a stub — see Decision D1).

## READS (per container)

| container | what | role |
|---|---|---|
| `built.Y` | h, T, p_gas, nh4, no3, Theta_L, C_NH4_granule, C_*_ice | accepted t_n (old time level) |
| `built.forcing` | slice `[kt]`: Ta/RH/wind/Precip/Rin/Rli/Ca/... | atmospheric driver |
| `built.boundary` | BC codes/values (top/bottom water/heat/gas/N) | BC, + `Tss` patched into heat top BC |
| `built.parameters` | soil (VG/thermal), soil_n (rates/RHOKG), grid dims | static params |
| `runtime.coupling` | `Srt`, `Tss`, throughfall, evap (step-frozen) | injected ebal/plant outputs |
| `runtime` | groundwater container, substep dt state, NIT/Tol_* | run services |
| `control` | switches (lThmrlefc/lSoilairefc/lNefc/lFreezing/lDualPor), time | gating |

## WRITES (per container)

| container | what | when |
|---|---|---|
| `workspace.step.trial` | hh/TT/P_gg/C_*_new (TDMA outputs) | every Picard iter |
| `workspace.step.previous` | h^{n-1}/backups | step entry + on accept |
| `workspace.derived` | Se/KL_h/KfL_h/Theta_LL/RHOV/L/TransportCoefficient/*BAR | every Picard iter |
| `workspace.solver` | coeff bank + assembled C1..C9/RHS + SAVE (per equation) | every Picard iter |
| `workspace.flux` | Srt copy, reaction rates Smin/Snit/Sden/Svol, V_L/V_A | per step |
| `built.Y` | h←hh, T←TT, p_gas←P_gg, nh4/no3←C_*_new, granule/ice advance | **on accept only** |
| `built.diagnostics` | CHK2/CHK3, NH4_total/NO3_total, plantN/nsc (N step) | post-substep |
| `runtime` | substep dt for next step, resume checkpoint on failure | per step |

Single-writer holds: the soil kernel only writes Y (on accept), workspace, diagnostics,
runtime. It does NOT mutate parameters/forcing/boundary (Tss is injected via runtime.coupling,
not written back). Audit confirmed no legacy single-writer violations.

## Accept/reject choreography (legacy → new containers)

Legacy 3-level structure: outer KT loop (:349) ⊃ adaptive substep loop (:1566) ⊃ Picard (:1657).

```
runOneStep(model, kt):
  fslice          = slice(forcing, kt)
  outer_backup    = Y                              # :367  (rollback whole step)
  ws.step.previous = Y                             # extrapolation + rollback anchor
  sub_elapsed = 0;  dt = runtime.substep.dt
  while sub_elapsed < Delt_t_outer:               # :1566
    sub_backup = {Y, ws}                           # :1573
    ws.step.trial = extrapolate(Y, ws.step.previous)   # :1594  hh=h+α(h−hOLD)
    for kit = 1..NIT:                              # :1657 Picard
      updateDerived(ws, Y, trial)                  # :1684 UpdateSoilWaterContent
      assemble(ws.solver, ...)                     # :1696 calculateMatricCoefficients
      trial.h     = solveWater(ws, Y, bc, Srt)     # :1731 → CHKWATER, QMT/QMB
      trial.p_gas = solveGas(ws, Y, ...)           # :1801 (if lSoilairefc)
      trial.T     = solveThermal(ws, Y, bc.Tss)    # :1823 → CHK, QET/QEB (if lThmrlefc)
      if Error_Water<Tol_Water && Error_Energy<Tol_Energy: break   # :1849
    if converged:
      ws.step.previous = Y;  Y = trial             # :2023 ACCEPT  (immediate, clean)
      sub_elapsed += dt;  dt = grow(dt, iters)     # :2046 ×2 after 4 low-iter successes
    elseif dt > dt_min:
      {Y, ws} = sub_backup;  dt = dt*0.5           # :1897 SUBSTEP REJECT, retry
    else:
      Y = outer_backup;  step_report.failed = true # :1922 OUTER ABANDON
      runtime.checkpoint = resumeState;  break      # :1999 emit resume
  if ~failed:
    solveNitrogen(Y, ws, ...)                      # :2235 post-substep → Y.nh4/no3 + diag
  step_report = finalize(QMT,QMB,QET,QEB,errors,substeps,iters)   # :2076
```

Note: legacy commits the previous KT's trial at the NEXT KT's top (:468-505, a lagged commit).
The new design accepts immediately at substep convergence (`Y←trial`), per handoff §3 — cleaner,
and diagnostics read `Y` directly. This is an intentional deviation (the lag is an artifact, not
physics) — flagged as Decision D4.

## step_report (proposed)

`converged`(bool), `failed`(bool), `n_substeps`, `picard_iters`(last), `dt_used`,
`QMT/QMB`(water top/bottom), `QET/QEB`(energy top/bottom), `Error_Water/Error_Energy`,
`reject_count`. (The registry still owns the ~185 logged fields; step_report is the caller's
lightweight per-step summary.)

## Decisions for ratification

- **D1 — Srt/Tss source for the first kernel.** Soil kernel reads `runtime.coupling.{Srt,Tss,...}`.
  For the first milestone, how is it filled (ebal/plant not migrated yet)?
  (a) **stub**: `Srt=0`, `Tss=forcing.Ta` — validates the container plumbing + TDMA + accept/reject
      against a controlled case (recommended first); then
  (c) migrate ebal+plant to write `runtime.coupling` for the faithful coupled step.
  vs (b) capture Srt/Tss from a legacy run for a faithful 1-KT gate-A.
- **D2 — substep ownership.** `runOneStep` owns the adaptive substep + Picard internally
  (one call = one outer forcing step) — confirm, vs exposing substeps to the caller.
- **D3 — coupling/runtime homes.** `Srt/Tss/top-BC` live in `runtime.coupling` (single writer =
  future canopy/ebal step); adaptive-dt state + resume checkpoint live in `runtime` — confirm.
- **D4 — immediate accept.** Drop the legacy lagged commit; accept `Y←trial` at convergence —
  confirm (cleaner; deviates from legacy line ordering but not physics).
- **D5 — gate-A strategy.** Without ebal/plant, a bit-faithful 1-KT match needs injected
  legacy Srt/Tss. Accept a **plumbing gate-A** first (mass/energy balance closes, accept/reject
  transitions fire, Y advances monotonically on a stub) and defer bit-faithful comparison to
  after ebal/plant migration — confirm.
```
