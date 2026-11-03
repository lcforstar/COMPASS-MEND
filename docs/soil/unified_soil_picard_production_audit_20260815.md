# Unified Soil Water–Heat–Gas Modified-Picard Production Audit

Date: 2026-08-15
Editor: Jizixing Zhu

## Scope and decision boundary

This audit covers the production integration of the reviewed soil
water–heat–dry-air changes. It does not introduce Newton iteration, disable
`KLa`, or reinterpret a conservation diagnostic as a physical equation. The
accepted solver remains Modified Picard. Ordinary sequential Picard is used
first; a block Modified-Picard map is activated only after the coupled
water–heat–air iteration stalls.

The local `data/sites/changwu (2)` package is used only as the annual validation
fixture. That folder and all local ZIP archives remain untracked and are not
part of the commit or tag.

## Implemented contracts

### 1. One element-flux owner

`soil.coupling.evaluateCoupledElementFlux` now evaluates liquid-phase Darcy
flow, gas-phase Darcy flow, liquid water, water vapor, dry air, and their energy
carriers from one candidate state and one coefficient epoch. The active liquid
pressure term is

```text
q_l,Pg = -(K_l / Gamma_w) grad(P_g)
```

and the same untruncated `grad(P_g)` is used by water, dry-air, and energy
consumers. The former gas/heat-only `±1e5` gradient clip remains diagnostic only.
The current arithmetic endpoint FEM rule is preserved; this patch changes
ownership and cross-equation consistency, not the spatial discretization.

Water vapor is owned once:

```text
j_v = j_v,h + j_v,T + j_v,Pg
```

The vapor advection term is no longer reconstructed and added a second time in
another equation bank. Heat consumes the same component fluxes:

```text
j_E = j_conduction + c_l T j_l + (c_l T + L_v) j_v + c_a T j_a
```

### 2. One exact candidate residual and ledger owner

`soil.coupling.evaluateResidualLedger` detaches one candidate, refreshes its
complete constitutive state, assembles all enabled equations from the shared
flux family, and returns:

- layer residual blocks for water, dry air, and energy;
- exact column ledgers from the same candidate;
- exact applied top boundary fluxes;
- dry-air admissibility and exact surface state;
- a candidate-state hash and a caller-unchanged assertion.

The same exact dry-air and heat boundary fluxes used for acceptance are
published to the outer accepted ledger. This removes the previous mismatch in
which acceptance could use the refreshed candidate while the outer ledger used
the preceding frozen Picard map.

### 3. Surface-water complementarity

The surface boundary is an active-set problem with retained surface-water depth
`d`, surface pressure head `h`, and overflow reaction `r`:

```text
d1 - d0 = dt (P_liq - E_actual - q_in - r)
0 <= d1  _|_  d1 - h >= 0
0 <= r   _|_  dmax - d1 >= 0
```

The dry boundary uses its own head-scale tolerance and requires a solved dry
endpoint plus a crossing of the configured dry-head threshold. Candidate head,
water flux, surface ledger, and active-set decision are inseparable. A rejected
branch is retried from the same accepted state and cannot enter the ledger.

### 4. Dry-air admissibility and production conservation gate

Dry-air density and storage are constrained to their nonnegative physical
domain with a coupled step fraction and explicit active-set boundary reaction.
The production column-mass gate approved for this integration is

```text
tol_air = max(1e-8 g cm-2,
              1e-6 * accepted_column_dry_air_mass)
```

Both the absolute and relative terms and the effective tolerance are reported.
This is a conservation acceptance gate, not a change to the dry-air equation.

### 5. Modified-Picard globalization

The ordinary scalar Modified-Picard map remains the default. At deterministic
iteration 20, only when both heat and dry air are enabled and the ordinary map
has stalled, Anderson history is reset and a frozen `h/U/P_g` block map is
formed. A common residual line search accepts only a physically admissible
candidate that reduces the exact scaled merit. No Newton step is executed.

### 6. Surface–snow–soil energy budget

The scalar surface/soil interface solver now reserves the configured combined
energy tolerance for already-consumed soil and snow equation residuals. The
interface mismatch is accepted only inside the remaining triangle-inequality
budget. This prevents an interface residual just below its own tolerance from
pushing the combined surface–snow–soil ledger above that same tolerance.

### 7. Supporting production corrections

- The static vertical grid uses a 1 cm top band and reviewed depth bands below.
- Hourly shortwave alignment follows the forcing clock.
- Winter Wheat is selected through PFT 25.
- Vegetation candidate materialization preserves the accepted/trial contract.

## Current verification evidence

### Focused and contract tests

The clean staged snapshot passed all 23 direct MATLAB contract gates and all 16
`matlab.unittest` cases exercised by the focused production set. This includes
shared fluxes, exact residual repeatability, water-vapor non-duplication,
dry-air admissibility and mass correction, all four equation-switch
combinations, block Modified Picard, surface active sets, MODFLOW boundary
transactions, and the full snow–surface–soil transaction.

The two failures found during the first whole-suite audit were real staged
regressions relative to the old baseline and were closed before publication:

- MODFLOW test: the assertion required `1e-6 mm` but inherited the site's
  `1e-4 mm` production gate. The regression now explicitly solves at the
  asserted `1e-6 mm` gate; the production tolerance was not weakened.
- Warm snow meltout: `0.0105886 W m-2` exceeded the `0.01 W m-2` combined gate
  because the interface alone consumed the full budget. The remaining-budget
  criterion reduced the final combined residual to `3.14392 J m-2` over the
  hour (`8.73e-4 W m-2`).

The repository-wide runner is not a clean release gate in this MATLAB
environment: tests mutate the working directory/path, several legacy tests
require `run_stemmus_initialize`, ONNX tests require an unavailable toolbox and
an `F:\\...` artifact, one test requires unavailable `struct2array`, and MATLAB
R2025a can terminate after results with `free(): chunks in smallbin corrupted`.
These infrastructure failures are reported, not hidden. Focused tests were
therefore rerun from isolated clean snapshots with their documented entry
points.

### Difficult and stable windows

Production-default results before the final annual rerun:

| Window | Result | Trials / accepted | Iterations | Minimum step | Wall time | Max water residual | Max dry-air residual | Max energy residual |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| KT6014 | PASS | 2 / 2 | 24 | 3600 s | 3.74 s | `7.76e-6 mm` | `8.69e-8 g cm-2` | `2.21e-4 J cm-2` |
| KT8418 | PASS | 6 / 6 | 172 | 1800 s | 10.60 s | `3.99e-7 mm` | `1.15e-7 g cm-2` | `1.83e-4 J cm-2` |
| KT8602–8603 | PASS | 48 / 37 | 1763 | 225 s | 84.91 s | `5.04e-7 mm` | `1.69e-7 g cm-2` | `2.50e-4 J cm-2` |
| Orchard KT1–168 | PASS | 502 / 502 | 6771 | 3600 s | 197.18 s isolated | `1.74e-6 mm` | `2.91e-7 g cm-2` | `7.58e-4 J cm-2` |

All values were finite; `KLa` fallback count was zero and the coupled-system
Newton count was zero. The existing safeguarded local `U -> T` state inversion
may still use its scalar Newton/bisection implementation; that is not a coupled
soil nonlinear solve.
The difficult KT8602–8603 window used a minimum accepted line-search factor of
`0.5`; the other listed windows retained `alpha=1`.

For the orchard reference, the strict absolute-only dry-air gate required
58,718 nonlinear iterations, 1,540 trials, a 56.25 s minimum step, and
1,585.71 s. The approved mixed gate required 6,771 iterations and full-hour
steps: an 8.67-fold reduction in nonlinear iterations and an 8.04-fold isolated
wall-time speedup.

### Final annual Changwu Winter Wheat run

The corrected production source completed one uninterrupted cold-start
KT1–8760 Winter Wheat run with `status=PASS_SCIENTIFIC`:

| Metric | Annual result |
|---|---:|
| Outer hours accepted | 8760 / 8760 |
| Nonlinear attempts | 15,281 |
| Accepted / rejected attempts | 15,133 / 148 |
| Total Modified-Picard iterations | 231,774 |
| Minimum accepted substep | 225 s |
| Wall time | 6,729.57 s (112.16 min) |
| Maximum water residual | `7.19174e-5 mm` |
| Maximum dry-air residual | `1.73412e-7 g cm-2` |
| Maximum effective dry-air tolerance | `1.74184e-7 g cm-2` |
| Maximum dry-air residual / accepted inventory | `9.99992e-7` |
| Maximum energy residual | `5.56230e-3 J cm-2` |
| Maximum accepted `abs(grad(P_g))` | `3471.02` internal gradient units |
| State finite gate | PASS |
| Hidden `KLa` fallback | 0 |
| Coupled Newton use | 0 |

The final 59-node state remained finite. Its main ranges were
`h=[-123803.92,-1213.59] cm`, `T=[11.525,33.640] degC`,
`U=[21.716,49.172]` in the model's internal-energy state units, and
`P_g=[877402.57,877990.94]` in the active internal pressure units. Final
endpoint unfrozen water content was `[0.06287,0.23453]`; ice content was zero.

The exact accepted-attempt residual sums were `1.33127e-4 mm` signed water
(`4.57921e-3 mm` absolute), `-2.08784e-5 g cm-2` signed dry air
(`1.65767e-3 g cm-2` absolute), and `5.99324e-2 J cm-2` signed energy
(`2.39199 J cm-2` absolute). These are sums of local nonlinear acceptance
residuals, not a replacement for the model's cumulative physical ledger.

The block path was used by 1,947 accepted attempts. Its accepted line-search
distribution was `alpha=1` for 1,870 attempts, `alpha=0.5` for 76, and
`alpha=0.25` for one. Thus 96.05% of accepted block steps used the full update;
only 3.95% required damping. The remaining 13,186 accepted attempts stayed on
the ordinary scalar pressure-coordinate path.

Hot bare/sparse-canopy hours emitted the existing `EBAL45-DIAG` because its
temperature diagnostic triggers at `Ts>50 degC`. The printed EBAL component
residuals remained under the configured `5 W m-2` gate and the annual run did
not fail or invoke a soil fallback. This debug verbosity is a remaining
engineering cleanup item, not evidence of a soil-solver rejection.

Artifacts:

```text
/tmp/compass_production_unified_tests/
  wheat_00001_08760_default_production_final_v2/
```

SHA-256:

```text
summary.txt        ceb077548396932e0328e0cc60e664ae8adc8416844d30c91440f6f438047c4b
attempts.csv       616c0da287d9102e2d80050100878876e7e7e9f295263d9a40756e2742d10481
accepted_state.mat a581af7d994862d0c82eda34e33fd2e1a1527fa20fc753918a4cff617d364179
```

## Publication gate

The production commit/tag gate is closed by the following evidence:

1. the corrected KT1–8760 run completed with finite accepted states;
2. annual water, energy, and dry-air residuals satisfied their configured gates;
3. no hidden `KLa` fallback or coupled Newton step was reported;
4. the staged diff passed `git diff --cached --check`;
5. local site packages and archives remained excluded.
