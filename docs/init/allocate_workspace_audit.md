# allocateWorkspace audit (2026-06-25)

Read-only audit of the legacy workspace/diagnostic allocation, to design
`init.allocateWorkspace` (single-writer for `built.workspace` + `built.diagnostics`).
Driven by the contract: `built.Y` = accepted prognostic state; `built.workspace` =
time-step scratch (old/trial/derived/solver/flux); `built.diagnostics` = output buffers;
`built.p` stays small (plant_hydraulic VC/KS/gpmax only).

Dims: `NN = NL+1` (nodes), `NL` (layers/elements), `nD = 2` (time-level second dim of
element arrays), `nS` (sub-phases), `Nmsrmn = tS` (output steps).

## Sources audited

| # | routine | file |
|---|---|---|
| 1 | SolverCoefficients | `src/+core/SolverCoefficients.m` |
| 2 | StateVariables | `src/+core/StateVariables.m` |
| 3 | FluxVariables | `src/+core/FluxVariables.m` |
| 4 | StartInit / defineSoilVariables | `src/StartInit.m`, `src/+init/defineSoilVariables.m` |
| 5 | UpdateSoilWaterContent | `src/UpdateSoilWaterContent.m` |
| 6 | solveSoilMoistureBalance (+ assemble/calcTimeDeriv) | `src/+soil/+water/…` |
| 7 | solveThermalEquations | `src/+soil/+thermal/…` |
| 8 | solveDryAirEquations | `src/+soil/+gas/…` |
| 9 | setInitialNitrogen blocks ②③ | `src/+init/setInitialNitrogen.m` |
| 10 | diagnostics Registry / flushToNetCDF | `src/+diagnostics/…` |

## Headline findings

1. **`h` IS the accepted state** (`= h^n`, fixed across the Picard loop): the temporal
   terms `(hh−h)/dt`, `(TT−T)/dt`, `(P_gg−P_g)/dt`, `(C_*_new−C_*)/dt`,
   `(Theta_LL−Theta_L)/...` prove `h/T/P_g/C_NH4/C_NO3/Theta_L` are the previous-time-level
   accepted values, and `hh/TT/P_gg/C_*_new/Theta_LL` are the trial. So
   **`workspace.step.old ≡ built.Y`** (no separate copy). OLD/TRIAL pairs confirmed:
   `h↔hh` (water/calculateTimeDerivatives), `T↔TT` (thermal/assemble:60-79),
   `P_g↔P_gg` (gas/assemble:56,65,72), `Theta_L↔Theta_LL` (conductivity/calculateDTheta:48,80),
   `C_NH4↔C_NH4_new` / `C_NO3↔C_NO3_new` (nitrogen transport).
2. **`hOLD/TOLD/P_gOLD`** = `h^{n-1}`, only for the extrapolated first guess
   (`COMPASS.m:1597-1605`) + rollback → `workspace.step.previous`.
3. **`FluxVariables` is ~95% diagnostic** (QL*/QV*/QA*/QN*, boundary cvTop/cvBot, Nmsrmn
   series) → `built.diagnostics`; only `*BAR`/`D*DZ` edge averages are derived-cache.
4. **`SolverCoefficients` is ~entirely solver scratch** (Chh…Cag PDE-coeff bank + C1…C9 +
   RHS + N coeffs) → `workspace.solver`; the `*BAR` edge averages are derived-cache.
5. **P_g cadence differs**: air pressure `P_g` is held fixed across adaptive substeps while
   h/T advance (gas solver audit). Container home unchanged (Y.soil.P_g); flagged for the
   main-loop wiring, not allocateWorkspace.

---

## Master table — by target home

allocateWorkspace allocates everything below EXCEPT the "(already in built.Y)" and
"(already in built.parameters)" rows, which are listed only to show provenance.

### A. built.Y — accepted prognostic state  (NOT allocated here; owned by buildInitialState/buildPlant*)

| var | shape | note |
|---|---|---|
| h, T, P_g | NN×1 | accepted soil head / temp / air pressure (`h^n`) |
| Theta_L, Theta_U, Theta_I | NL×nD | accepted liquid / total / ice content (paired with trial below) |
| C_NH4, C_NO3, C_Urea | NN×1 | accepted N concentrations (`C^n`) |
| C_NH4_granule | NN×1 | **prognostic** slow-release granule pool (has step-to-step memory) — flag (audit split diag/prog; mass-balance ⇒ Y) |
| C_NH4_ice, C_NO3_ice | NL×2 | **prognostic** freeze-entrapped N (memory across thaw) — flag (default Y) |
| h_frez, hh_frez, T_CRIT, TT_CRIT | NN×1 | freeze-curve state |

### B. built.workspace.step.trial — t_{n+1} candidate the solver is solving

| var | shape | init |
|---|---|---|
| hh, TT, P_gg | NN×1 | copy of Y.{h,T,P_g} (warm start) |
| C_NH4_new, C_NO3_new, C_Urea_new | NN×1 | copy of Y.C_* |
| Theta_LL, Theta_UU, Theta_II | NL×nD | derived from trial head (recomputed) |
| hh_frez | NN×1 | 0 |

### C. built.workspace.step.previous — t_{n-1} (extrapolation guess + rollback)

| var | shape | init |
|---|---|---|
| hOLD, TOLD, P_gOLD | NN×1 | copy of Y.{h,T,P_g} |
| hOLD_frez | NN×1 | copy |
| C_NH4_OLD, C_NO3_OLD, C_Urea_OLD | NN×1 | copy of Y.C_* |
| h_backup, T_backup, P_g_backup, C_NH4_backup, C_NO3_backup, C_Urea_backup | NN×1 | [] (filled on substep entry, for rollback) |

### D. built.workspace.derived — recomputable each step from state

| group | vars | shape |
|---|---|---|
| moisture | Se, KL_h, KfL_h, KfL_T, KL_T, DTheta_LLh, DTheta_UUh, Ratio_ice, Gamma_hh, COR, CORh | NL×nD (Gamma_hh/COR/CORh NN×1) |
| thermal | L, LL, W, WW, Lambda_eff, c_unsat, EHCAP, ETCON, ZETA, DhT, MU_W, f0, L_WT, D_Ta | NN×1 / NL×nD |
| vapor/gas | RHODA, RHOV, RHOV_s, DRHOV_sT, DRHOVh, DRHOVT, DRHODAt, DRHODAz, Xaa/XaT/Xah, HR, D_A, D_V, Theta_g, Theta_V, k_g, Sa, Alpha_Lg, POR_C, Eta, Kg | NN×1 / NL×nD |
| edge averages | KLhBAR, KLTBAR, RHOVBAR, EtaBAR, Lambda_effBAR, c_unsatBAR, Theta_LLBAR, DTDBAR, DETBAR, DEhBAR, DRHOVhDz, DRHOVTDz, DIFBAR_NH4, DIFBAR_NO3, E1BAR, E2BAR, RetardBAR | NL×1 |
| gradients | DhDZ, DTDZ, DPgDZ, DNH4DZ, DNO3DZ | NL×1 |
| nitrogen transport | Retard, fw, DshNh4, DshNo3, DifNh4, DifNo3, Vpore, WFPS, CRootNh4, CRootNo3, cc | NL×nD (cc NN×1) |

### E. built.workspace.solver — PDE-coefficient bank + assembled tridiagonal (per equation)

| sub | vars | shape |
|---|---|---|
| `.coeff` (SolverCoefficients PDE bank) | Chh,ChT,Khh,KhT,Kha,Vvh,VvT,Chg / CTh,CTT,CTa,KTh,KTT,KTa,VTh,VTT,VTa,CTg,Kcvh,KcvT,Kcva,Ccvh,CcvT / Cah,CaT,Caa,Kah,KaT,Kaa,Vah,VaT,Vaa,Cag,Kcah,KcaT,Kcaa,Ccah,CcaT,Ccaa | NL×nD |
| `.water` | C1,C2,C4,C5,C6 (NN×2), C4_a,C5_a,C7,C9 (NN×1), RHS (NN×1) | — |
| `.thermal` | C1,C2,C3,C4,C5,C6 (NN×2), C4_a,C5_a (NL×1), C7 (NN×1), RHS (NN×1) | — |
| `.gas` | C1,C2,C3,C4,C5,C6 (NN×2), C4_a,C5_a (NL×1), C7 (NN×1), RHS (NN×1) | — |
| `.nitrogen` | R1,R1O,R2,R2O,E1,E2,B1,B2,F1,F2,G1,G2 (NL×nD), C1,C2,C3,C4,C5,C6 (NL×nD), C3_a,C4_a,C5_a,C6_a,C7,C9 (NN×1), C1O,C2O (NL×nD), RHS,RHS1 (NN×1), N_SAVE (3×3×5) | — |
| `.SAVE` | SAVE | 3×3×5 (shared; slots 1/2/3 = water/thermal/gas boundary flux) |

### F. built.workspace.flux — per-step sink / velocity scratch consumed in assembly

| var | shape | note |
|---|---|---|
| Srt | NL×nD | root water uptake sink (→ water C9) |
| SinkNh4, SinkNo3, AUptakeA, SinkS | NL×nD (SinkS NL×nD×2) | N uptake sinks |
| Smin, Snit, SVnit, Sden, SVden, Svol, SVvol, Shys | NL×nD | N reaction rates |
| V_L, V_A | NL×nD | pore-water / air velocities |
| Dsh_h, Dsh_T, Dsh_A | NL×nD×nS | dispersion sensitivities |

### G. built.diagnostics — output buffers (only written, never read by physics)

| group | vars | shape |
|---|---|---|
| water/vapor/air fluxes (FluxVariables) | QL,QL_h,QL_T,QL_A,QL_D,QL_disp,QL_Dts,QL_dispts,QLts / QV*,QVA,QV_D,QV_A,QV_disp,QV_*ts / QA_*,Qa,QA_*ts | NL×1 (QV* NL×nD) |
| N fluxes | QNH4, QNO3 | NL×1 |
| boundary series (Nmsrmn×1) | cvTop,cvBot,cvTop1,cvBot1, cTop_NH4,cBot_NH4,cTop_NO3,cBot_NO3, cvSnit,cvSmin,cvSden,cvSvol, cvSinkSnh4,cvSinkSno3, Cvnit_N2o,Cvdenit_N2o,Cvdenit_N2, cvnh4,cvno3,cvnh4OLD,cvno3OLD | Nmsrmn×1 |
| boundary scalars/step | QMT,QMB,QET,QEB,AVAIL0 | scalar/KT |
| state series (Nmsrmn×1) | VInfil, SPUptakeNh4,SPUptakeNo3, SAUptakeP,SAUptake,SAUptakeAN,SAUptakeA, SPOT | Nmsrmn×1 |
| node diag | Root_Stress, WFPS, Vpore | NN×1 / NL×nD |
| N totals | NH4_total, NO3_total | NN×tS (only if lNtotal) |
| convergence | CHK2, CHK3 | NN×1 (inf init; per-Picard residual) |
| plant C/N (struct-array, tS rows) | plantN_diag (12 fields), nsc_diag (15 fields) | tS×1 struct |

NOTE: most of group G is **already registry-managed** (`src/+diagnostics`, ~185 fields across
~22 groups, auto-allocated `nan([tS dims])` per `register()`). The LOOSE legacy arrays needing
an explicit home in `built.diagnostics` are: `plantN_diag`, `nsc_diag`, `NH4_total`, `NO3_total`,
`CHK2/CHK3`. (`N_SAVE` is ephemeral solver scratch → workspace.solver.nitrogen, not logged.)

### H. built.parameters / dropped legacy fields — static provenance only

Static fields such as `Imped`, `Ss`, `h_trans`, `Disp` and process rates are already
owned by `buildSoil*` / `buildNitrogenSoilParams` / `buildVeg`. `RHOL` is a physical
constant. The 2026-06-26 completeness pass moved the genuinely recomputed runtime
scratch fields `D_Vg` and `eT/eTheta/ePH` into `workspace.derived`; see the final
section below. Dead legacy fields (`ALPHA`, `alpha_h`, `bx`, `b1`, `b11`) are not
allocated.

---

## Proposed container skeleton

```
built.workspace
  .step
    .trial      {hh, TT, P_gg, C_*_new, Theta_LL/UU/II, hh_frez}
    .previous   {hOLD, TOLD, P_gOLD, hOLD_frez, C_*_OLD, *_backup}
  .derived      {Se, KL_h, KfL_h, RHOV, Theta_g/V, *BAR, D*DZ, Retard, ...}
  .solver
    .coeff      {Chh..Cag PDE bank}
    .water / .thermal / .gas / .nitrogen   {C1..C9, RHS, N coeffs, N_SAVE}
    .SAVE       (3×3×5 shared)
  .flux         {Srt, SinkNh4/No3, Smin/Snit/Sden/Svol, V_L, V_A, Dsh_*}

built.diagnostics   (loose legacy buffers; registry stays separate)
  {plantN_diag, nsc_diag, NH4_total, NO3_total, CHK2, CHK3}
```

## allocateWorkspace contract (proposed)

- `built = init.allocateWorkspace(built)` — single-writer for `built.workspace` +
  `built.diagnostics`; **does not touch** `built.Y`/`built.parameters`/`built.p`.
- Reads dims from `built.domain.grid` (NN/NL/nD/nS) + `built.control.time` (Nmsrmn=tS).
- trial/previous initialized as warm-start **copies of `built.Y`** (so the first Picard
  guess + first temporal term are well-posed); derived/solver/flux as `zeros`; diagnostics
  as `nan([tS dims])`.
- gate-A: shapes + init values vs the legacy `core.*Variables`/`setInitialNitrogen` blocks
  ②③ constructors (element-for-element where they overlap).

## Decisions — RATIFIED + IMPLEMENTED (2026-06-25)

1. **`workspace.step.old` collapses into `built.Y`** — DONE. legacy `h`≡`Y.h`; trial
   warm-starts from Y; on accept `Y ← trial`, old Y → `step.previous`.
2. **Granule + freeze-ice N pools** (`C_NH4_granule`, `C_*_ice`) → `built.Y` as prognostic —
   DONE. Added to the Y single-writer `buildInitialState` (zeros at IC), NOT to workspace.
3. **solver granularity**: per-equation `.solver.{water,thermal,gas,nitrogen}` + shared
   `.coeff`/`.SAVE` — DONE.
4. **diagnostics scope**: allocateWorkspace seeds only the 6 LOOSE legacy buffers
   (plantN_diag/nsc_diag/NH4_total/NO3_total/CHK2/CHK3); registry keeps the ~185 — DONE.

## Status

`init.allocateWorkspace(built) -> built` implemented (`src/+init/allocateWorkspace.m`),
single-writer for `built.workspace` + `built.diagnostics`, does not touch `built.Y`.
gate-A `tests/test_allocateWorkspace.m` 52/52; full init suite 578/578; checkcode clean.
Next: `init.initialize` (top-level orchestration loadSite→checkSite→buildState→allocateWorkspace).

## Legacy-coverage completeness check (2026-06-26)

Programmatic diff of legacy `core.{StateVariables,SolverCoefficients,FluxVariables}`
properties vs the new container leaf-fields (rename map P_g→p_gas, C_NH4→nh4,
hh→trial.h, …):

- **SolverCoefficients 81/81** → `workspace.solver`(+`*BAR`→derived). 100% ported.
- **FluxVariables 8/59 by name** — the 8 are gradients/`*BAR`→`workspace.derived`; the
  other 51 are OUTPUT fluxes (QL*/QV*/QA*/QN*/cv*) owned by the diagnostics REGISTRY
  (snake_case, ~325 fields), by design NOT stored in `built.*`. Accounted for.
- **StateVariables 100/122**; the 22 not-by-name split, after checking real consumers
  in the compute path:
  - registry diagnostics (VInfil/SPUptake*/SAUptake*/SPOT/Root_Stress), constant RHOL
    (→ PhysicalConstants), `Disp`→`parameters.soil_n.DISP` (case): not gaps.
  - **MIGRATED to `workspace.derived` (genuinely recomputed + consumed each step):**
    `T_CRIT`/`TT_CRIT` (NN×1, HT_frez → KfL_T frozen-flux), `eT`/`eTheta`/`ePH` (NL×1,
    solveNitrogen rate modifiers), `D_Vg` (NL×1, vapor/gas term in matric coeffs).
    Added 2026-06-26; zero-init (recomputed before use). Physical-equivalence check
    deferred to when the freeze/N/gas runtime lands in `runOneStep`.
  - **NOT migrated (dead — no compute consumer; confirmed by grep):** `alpha_h`
    (CalculateLAIAndStressFactors uses a function-LOCAL var, never reads the field),
    `bx`, `b1`, `b11`, `ALPHA`. Left in legacy `StateVariables`; will vanish when the
    old entry retires. Migrating them would pollute the new structure with phantom state.
- **NitrogenVariables** (struct): split across `Y.plant_cn` (pools) + `Y`
  (concentrations/granule/ice) + `diagnostics` (plantN/nsc/NH4_total/NO3_total/CHK2/CHK3)
  + `workspace.flux` (Smin/Snit/Sden/Svol/Sink) + `parameters.soil_n` (rates).

gate-A after migration: `test_allocateWorkspace` green (adds 6-field presence/shape/zero
+ dead-field exclusion checks); input suite 12/12 green; checkcode clean.
