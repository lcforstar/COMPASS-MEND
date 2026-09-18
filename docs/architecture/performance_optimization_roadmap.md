# COMPASS runtime optimization roadmap

Date: 2026-06-16

Scope: one-year Changwu apple orchard run with the coupled physics kept on:

- `lFreezing=1`
- `lThmrlefc=1`
- `lSoilairefc=1`
- `lNefc=1`
- hourly outer forcing step, adaptive inner substeps

Goal: reduce the current coupled run toward less than 1 h without disabling freeze-thaw, thermal conduction, dry-air coupling, or nitrogen cycling.

## Current Evidence

### Original Pre-P2 Evidence

The original performance report indicated total profiled runtime of about `37777 s` (`10.49 h`). The dominant path was:

| Component | Profile time | Interpretation |
|---|---:|---|
| `COMPASS` | `37777 s` | full profiled run |
| `UpdateSoilWaterContent` | `20756 s` | soil water state and hydraulic-property refresh |
| `calculateHydraulicConductivity` | `20705 s` | main hydraulic/freeze-thaw property kernel |
| `getSoilConstants` | `7801 s` | high-frequency constant lookup inside scalar kernels |
| Dall'Amico/theta/dtheta helpers | several thousand seconds | repeated scalar helper calls |
| `ebal_clm45` | `3465 s` | canopy energy/photosynthesis/hydraulic loop |
| `solveNitrogen` | `2050 s` | soil N transport and reactions |
| `Registry.write` | `2020 s` | per-hour diagnostic writes through `containers.Map` buffers |
| periodic/full flush | `~1183 s + 810 s + 471 s` | partial/full NetCDF and CSV output |

The current code already contains important stability and acceleration mechanisms:

- Soil Picard Anderson acceleration exists in `src/COMPASS.m`.
- ebal Anderson acceleration exists in `src/ebal_clm45.m`.
- Adaptive substepping and rollback exist in `src/COMPASS.m`.
- Current `checkpointInterval=720`, so periodic partial output is monthly, not daily.

Therefore the remaining large gain must come from reducing the cost per Picard/substep solve, not from simply adding another convergence wrapper.

## Profile-Derived Opportunity Classification

The profiler PDF is a screenshot-style MATLAB profile report. The root page
shows that the first-order cost is not generic MATLAB overhead; it is repeated
scalar physics work inside the water/thermal/freeze-thaw solve. The following
groups should guide the next optimization passes.

| Class | Profile evidence | Code path | Priority | Action |
|---|---:|---|---|---|
| Hydraulic theta kernel | `calculateHydraulicConductivity` `20705 s`; `calculateTheta_LL` `2770 s`; `calculateTheta_UU` `2373 s`; `calculateTheta` `2151 s`; `calculateGamma_hh` `1112 s`; `calculateSe` `1102 s` | `src/+soil/+water/+conductivity/*`, called by `UpdateSoilWaterContent` inside the Picard loop | Already active | P2 MATLAB array kernel, scalar reference retained |
| DTheta and K residual loop | `calculateDTheta_UUh` `2704 s`; `calculateDTheta_LLh` `2508 s`; `calculateDTheta` `2333 s`; `calculateKL_h` `171 s`; `heaviside1` `584 s` | still scalar inside `calculateHydraulicConductivityArray.m` after P2.2 | Immediate | P2.3: vectorize derivatives and hydraulic/thermal conductivities before touching broader modules |
| Dall'Amico phase-change scalar kernel | `calculateDallAmicoPhaseChange` `3446 s`; children include `vgSeFromHead` `264 s`, `clamp` `208 s`, `boundedPhaseStep` `166 s`, `vgHeadFromSe` `120 s`, `getThermalCapacity` `195 s`, `limiterCode` `140 s` | currently called per node in both scalar and array hydraulic kernels | Immediate after P2.3 | P3/P2.4: add `calculateDallAmicoPhaseChangeArray.m` with array outputs and scalar-vs-array tests |
| Soil water matrix assembly | `calculateMatricCoefficients` `1869 s`; `solveSoilMoistureBalance` `1100 s`; `calculateBottomTempBoundary` `769 s`; tridiagonal solves are much smaller | `src/COMPASS.m` Picard loop and soil water/thermal coefficient functions | High after hydraulic kernel | Hoist constants, precompute geometry factors, reduce repeated coefficient assembly if state is unchanged |
| Thermal and dry-air parameter refresh | `solveThermalEquations` `656 s`; `calculateEnergyParameters` `354 s`; `calculateThermalConductivityCapacity` `161 s`; `calculateThermalLiquidDispersion` `152 s`; `solveDryAirEquations` `240 s`; `calculateDryAirParameters` `238 s`; gas diffusivity/conductivity `~158 s` combined | `src/+soil/+thermal/*`, `src/+soil/+gas/*` | Medium-high | Vectorize parameter refresh and cache geometric/constant terms; preserve conservation and coupling |
| Canopy energy and biochemistry | `ebal_clm45` `3465 s`; `biochemical` `1868 s`; `compute_assimilation` `614 s`; `solveEbalHydraulics*`/plant hydraulics `~470 s`; `Plantrootconductivity` `374 s` | `src/ebal_clm45.m`, `src/+plant/*`, `src/+PlantHydraulic/*` | Medium | Warm starts and caching after hydraulic kernel is cheaper; do not reorder photosynthesis/hydraulic coupling without A/B |
| Nitrogen and organic matter | `solveNitrogen` `2050 s`; `computeOrganicSources` `1865 s`; `OrganicStandard.tick` `1269 s`; `OM.turnover_pool` `414 s`; `turnover` `122 s` | `src/+soil/+nitrogen/*`, organic matter classes/functions | Medium | Hoist invariant C/N terms, cache temperature/moisture response factors, keep NH4/NO3 residual tests strict |
| Output and diagnostics | `Registry.write` `2020 s`; `saveCompassPeriodicPartialOutput` `1183 s`; `flushToCsv` `810 s`; `flushToNetCDF` `471 s`; `writematrix`/`writetable` `~1590 s` | `src/+diagnostics/*`, `src/+io/saveCompassPeriodicPartialOutput.m` | Mostly done | P1 registry/output work; keep CSV optional and partial output cheap |
| Micro-hotspots from repeated utilities | `Geometry1D.total_soil` `777 s`; `Geometry1D.cell_volume` `569 s`; `Geometry1D.dz` `353 s`; `Geometry1D.zminus` `126 s`; `PhysicalConstants.instance` `431 s`; `mean` `862 s`; `firstNonSingletonDim` `120 s`; `daisy_assert` `128 s`; `Calc_Esat_Flatau` `1062 s` | many modules | Low until big kernels are fixed | Cache geometry/constant properties locally in hot loops; avoid optimizing generic MATLAB utilities before module-level work |
| Radiation transfer | `RTMf` `1820 s`; `RTMz` `788 s`; `RTMo` `336 s`; `RTMt_sb` `58 s`; `RTMo.calcJ*` `~111 s` | canopy radiation/SCOPE path | Medium-low | Cache or reuse radiation only if forcing/canopy state permits; high scientific risk if angular/spectral state is stale |

Original interpretation at that time: P2/P3 were the right next moves. That
decision has now been executed and verified. The later post-P2/P3 profile below
supersedes this original priority order.

### Post-P2/P3 Two-Year Profile Reset

Latest profile:

- Folder: `data/output/_profile_2yr/`
- Raw profile: `profile_info.mat`
- Run window: `2021-01-01` to `2023-01-01` (`730 d`, `17520` hourly steps)
- Config: `lFreezing=1`, `lNefc=1`, `lSoilairefc=1`, `HydraulicKernel=array`
- Profile total for `COMPASS`: `14152.589 s`
- Rollbacks: `3763`
- Fatal/fail count: `0`
- Final NetCDF output only: `flushToNetCDF = 20.938 s`; periodic partial output was not active

The post-P2/P3 ranking is different from the original profile. The hydraulic
freeze-thaw kernel is now a residual cost, not the leading bottleneck.

| Rank | Function/block | Profile time [s] | Calls | Interpretation |
|---:|---|---:|---:|---|
| 1 | `ebal_clm45` | `3041.165` | `17520` | canopy energy, biochemistry, plant hydraulics |
| 2 | `soil.water.calculateMatricCoefficients` | `2948.191` | `1614466` | new soil-water coefficient bottleneck |
| 3 | `solveSoilMoistureBalance` | `1759.974` | `807233` | dominated by coefficient assembly |
| 4 | `vegetation.biochemical` | `1559.995` | `302708` | Ci/stomata/FvCB and output packing |
| 5 | `calculateBottomTempBoundary` | `1167.447` | `807233` | repeated bottom-temperature/date-series work |
| 6 | `solveThermalEquations` | `1073.495` | `807233` | thermal coefficient/matrix solve |
| 7 | `Calc_Esat_Flatau` | `1029.887` | `756771` | saturation vapor utility, repeated polynomial/mask work |
| 8 | `calculateHydraulicConductivityArray` | `701.712` | `824754` | P2/P3 residual; no longer first priority |
| 9 | plant hydraulics chain | `~1400` combined | `151354+` | root-side/root-soil conductance inside ebal |
| 10 | diagnostics/output | `~141` registry + `21` NetCDF | millions | not a current bottleneck |

Line-level evidence for the next work:

- `Calc_Esat_Flatau`: coefficients are rebuilt on every call and the 9-term
  polynomial repeatedly indexes `mask_water` / `mask_ice`.
- `calculateBottomTempBoundary`: every Picard call recomputes layer means and
  regenerates the full Julian-day series before writing one hourly bottom
  temperature value.
- `soil.water.calculateMatricCoefficients`: the current `NL x nD` nested loop is
  the new soil-water hotspot after the hydraulic property kernel was vectorized.
- `Plantrootconductivity`: about `38.7M` endpoint visits through
  `rootSideConductance`, `rootSoilContactFactor`, and related scalar helpers.
- `ebal_clm45`: the largest lines are repeated calls to `vegetation.biochemical`,
  `heatfluxes_clm45`, and `solveEbalHydraulicsLayered`; the P5 warm-start
  measurement showed the ebal iteration count is not controlled by cold-start
  distance.

Status update after P4-a/P4-b:

- P4-a was implemented as a conservative Flatau utility cleanup: skip unused
  `qs`/`dqs_dT` outputs at hot call sites and remove repeated masked
  gather/scatter in the polynomial loop. Horner evaluation and persistent
  constants were not adopted because the accepted implementation preserves the
  existing polynomial accumulation order and avoids persistent model state.
- P4-b was implemented as bottom-temperature date-series caching plus
  `mean(X,2)` replacement for the layer-mean loop. The layer-mean loop was the
  real cost; date-series caching was a small secondary gain. This also removes
  the dominant source of the global `mean` hotspot from the two-year profile.
- The active next bottleneck is therefore P4-c, not P4-a/P4-b follow-up work.

Current conclusion:

- P2/P3 achieved their purpose: `calculateHydraulicConductivityArray` is now
  about 5% of profiled COMPASS time.
- Do not spend the next effort on MEX or further hydraulic-kernel polishing until
  the new soil matrix, vapor utility, bottom-temperature, plant-hydraulic, and
  biochemical hotspots are reduced.
- The next sequence should be low-risk utility/cache work first, then the
  larger soil-water matrix vectorization.

New immediate order:

1. Commit the verified P4-a/P4-b changes and the unrelated `n` header fix as
   separate rollback points.
2. P4-c: vectorize `soil.water.calculateMatricCoefficients`.
3. Fresh short profile after P4-c to confirm whether soil matrix assembly or
   canopy/thermal work is now dominant.
4. P4-d: vectorize or cache endpoint logic in `Plantrootconductivity`.
5. P4-e/P5b: reduce repeated `biochemical`/`ebal`/radiation work without changing nonlinear
   solver semantics.

Defer for now:

- `Registry.write` and final NetCDF output: too small in the current profile.
- ebal previous-hour warm starts: measured and rejected.
- MEX: wait until the new MATLAB hotspots are reduced.

Status update after P4-c (2026-06-16, fresh 1-year array-kernel profile):

- P4-c landed: `soil.water.calculateMatricCoefficientsArray`, selected by
  `cfg.MatricKernel` (default flipped to `'array'`). Verified bit-identical
  (`max_diff=0`) function-level across all 6 lhThmrl/gate combos, and end-to-end
  on a 3-day winter freeze-thaw window and a summer spin-up window; ~10x on the
  matric-kernel micro-benchmark. Commits `4b040f8` (opt-in) + `889bb80` (default).
- Fresh 1-year profile (array matric + array hydraulic), normalized x2 vs the
  two-year baseline:
  - water matric: `2948 s` (#6) -> `~322 s` (`...Array`), **~9x**.
  - `calculateBottomTempBoundary`: `1167 s` (#9) -> dropped out of the top 50.
  - `mean`: `1098 s` / 230M calls (#10) -> `~300 s` / ~26M calls (#34); the
    bottom-temperature layer-mean loop was ~90% of the global `mean` calls.
  - `Calc_Esat_Flatau`: ~10% off only (real saving was always small; the
    profiled `1030 s` was call-overhead inflated).
  - Net: profiled COMPASS 1-year `6693 s` vs two-year baseline normalized
    (`7076 s`) = **~5.4% drop**. Cumulative pre-P2 -> now: `~17834 -> ~5590 s/yr`
    (**~3.2x**, of which P2 is the bulk and P4 is ~5-6%).
- New ranking (1-year, post-P4): Canopy/biochem/radiation `5700 s` (class #1) >
  Thermal/dry-air `2241 s` > P2 hydraulic `590 s` > Soil matrix `299 s` (was the
  big scalar-matric cost, now small). Top leaves: `ebal_clm45 1876 s` (aggregator),
  `biochemical 921 s`, `RTMf 404 s`, `solveEbalHydraulics* 337+332 s`,
  `compute_assimilation 318 s`, `Plant_hydraulics`/`calSPAC 315+313 s`,
  `Plantrootconductivity 271 s`, `rootSideConductance 167 s` (18.8M calls).
- Next targets: P4-d plant-hydraulics endpoint vectorization
  (`Plantrootconductivity` + `rootSideConductance` + `sel_root`) is the remaining
  pure-vectorization win inside the canopy class; the thermal/dry-air cluster
  (`solveThermalEquations 710 s`, `calculateEnergyParameters 395 s`, ~384k calls
  each) is a second vectorization candidate. The ebal/biochemical iteration count
  itself is not reducible by warm-start (P5 measured + rejected).

## Code-Grounded Levers (scout audit of classes 4-7, 2026-06-15)

Static read of the actual source for the four classes that become visible after
the hydraulic kernel (P2/P3). Candidate levers only — NOT run-verified; each must
be A/B'd per the Verification Policy before adoption. Two incidental findings
surfaced (one fixed).

**Class 4 — Soil water matrix assembly**
- `calculateBottomTempBoundary` (called every Picard iter, `COMPASS.m:1670`)
  regenerates the entire year's Julian-day series (~365 values) via
  `generateJDSeries` on every iteration, only to produce one scalar deep-temperature
  boundary that barely changes hour-to-hour (~44k-88k calls/yr). Lever: cache
  `jd_series` and damp depth `De` once per outer step (KT), index inside the
  Picard loop. ~20-50x on that function, zero physics risk. **Highest-ROI single
  lever found.**
- `calculateMatricCoefficients:19-21` re-fetches `RHOL/Gamma0/Gamma_w` + `cfg.Grid`
  256x per Picard iter → cache once (5-10%). `calculateTimeDerivatives` repeats
  `/Delt_t` → hoist `inv_Delt_t` (free, zero risk).

**Class 5 — Thermal and dry-air parameter refresh**
- Nine thermal/gas functions called unconditionally every Picard iter
  (`COMPASS.m:1656-1668`), each an NL×2 per-node loop re-fetching
  `core.PhysicalConstants.instance()`. Dominant driver is iteration count
  (~260k/yr), not loop size — same disease the hydraulic kernel had.
- Levers: cache physical constants + `cfg` flags once per KT (15-25% + 5-10%, low
  risk); vectorize the de Vries per-layer thermal-properties loop (20-30% of that
  function, medium risk — conditional logic). Redundancy: dry-air diffusivity
  `D_A` is computed independently in both `calculateVaporVariables` and
  `calculateSoilThermalProperites` — compute once, share.

**Class 6 — Canopy energy and biochemistry (iteration-count, NOT vectorization)**
- `ebal_clm45:136-137` cold-starts `Tch/Tcu` from forcing every hour (no warm
  start). Nesting: ebal `while` (30-50 iters, stress cases 300-900) × hydraulic
  loop (6) × Ci iteration (50) → ~360-600 `biochemical` calls per hour.
- Levers: warm-start `Ci` from the previous accepted ebal iteration (40-50% fewer
  Ci iters); cache `fDYL`/sun position per hour (now recomputed 2-12x); cache
  pre-`sfactor` `Vcmax25`/temperature-response across the hydraulic loop.
- CAVEAT: scout estimates `biochemical` is only ~1-2% of de-inflated annual
  runtime; its value is stress-case reliability, not raw time. True rank pending
  the post-P2/P3 1-year profile.

**Class 7 — Nitrogen and organic matter**
- `solveNitrogen:183` Picard loop (`max_iter=100`, ~5 typical) recomputes the
  response factors `eT/eTheta/ePH/Retard` (per-node loops in `mineralization`,
  `calculateRetardationFactors`, `soluteDispersion`) every iteration, although
  `TT/Theta_LL/PH1` are frozen during the N loop (only concentrations change).
- Lever: hoist the response factors out of the iteration (compute once per KT) →
  50-70% on the N solve, low-medium risk; concentration-dependent
  nitrification/denitrification/transport MUST stay inside. Then vectorize the
  quadruple-nested `soluteDispersion` Dsh assembly.

**Two incidental findings (not optimizations):**
1. The `calculateBottomTempBoundary` year-long date-series regeneration above — a
   ~20-50x waste, not a correctness bug.
2. `soluteDispersion.m` `lArtD=1` branch used a stale `MN` (= `NN`) instead of
   `ML+ND-1` for the artificial-dispersion `DD` term → every element used the
   bottom node's flux. **Fixed 2026-06-15 (commit `cd87b77`); dormant (`lArtD=0`
   in production).** Also: the `eT` block in the same file is a no-op because
   `EA=0` is hardcoded — dispersion temperature modulation is effectively off
   (open science question for the N-module owner).

## Priority Order

### P0 - Measurement Harness

Objective: create reproducible timing and numerical-difference evidence before changing kernels.

Work:

- Add a short-window performance harness for 3-day and 10-day runs.
- Keep all target physics switches on.
- Time major blocks with low-overhead timers:
  - `UpdateSoilWaterContent`
  - `calculateHydraulicConductivity`
  - Dall'Amico phase change
  - `solveThermalEquations`
  - `solveDryAirEquations`
  - `ebal_clm45`
  - `solveNitrogen`
  - `Registry.write`
  - partial/final flush
- Save benchmark metadata beside temp output, not in active long-run output folders.

Files likely touched:

- `tests/benchmark_compass_runtime.m` or `tools/benchmark_compass_runtime.m`
- optionally `src/+diagnostics/PerfTimer.m`

Estimated implementation time: 0.5-1 day.

Expected speed gain: none directly.

Acceptance:

- Same benchmark can be rerun with a single command.
- Benchmark output reports per-block wall time, call count, and total runtime.
- Output path is temp-isolated.

### P1 - Output and Registry Quick Wins

Objective: remove avoidable overhead that does not affect physics.

Work:

- Refactor `diagnostics.Registry` internals:
  - assign integer ids at registration time;
  - store buffers in cell arrays or direct properties;
  - avoid repeated `containers.Map` get/set in `write`;
  - fast-path `transform='none'`.
- Add runtime output policy:
  - keep NetCDF as primary output;
  - make legacy CSV optional or final-only;
  - stop writing full CSV during periodic partial snapshots by default.
- Change periodic partial output:
  - write resume/checkpoint state cheaply;
  - optionally write NetCDF partial only;
  - avoid full CSV duplication unless explicitly requested.

Files likely touched:

- `src/+diagnostics/Registry.m`
- `src/+diagnostics/flushToCsv.m`
- `src/+diagnostics/flushToNetCDF.m`
- `src/+io/saveCompassPeriodicPartialOutput.m`
- `src/+io/saveCompassResumeCheckpoint.m`
- `src/+core/ModelSettings.m`
- `data/config_main.txt`

Estimated implementation time: 1-2 days.

Expected speed gain: about 15-30% of current runtime in output-heavy profiled runs; likely 10.5 h to 7-8.5 h depending on partial flush frequency and profiler overhead.

Acceptance:

- NetCDF output unchanged for final runs.
- CSV can still be generated when requested.
- Existing plot scripts either keep working in CSV-enabled mode or are given an offline converter path.
- Conservation diagnostics match before/after within roundoff.

Risk:

- Medium for backward compatibility, low for physics.

### P2 - MATLAB Array Hydraulic Kernel

Objective: remove the scalar-node execution pattern in the main bottleneck while keeping formulas in MATLAB.

Current bottleneck:

`src/+soil/+water/+conductivity/calculateHydraulicConductivity.m` loops over `NL x nD`, constructs a per-node `SV` struct, and calls many scalar helper functions for every node and Picard iteration.

Work:

- Add a new array implementation:
  - `src/+soil/+water/+conductivity/calculateHydraulicConductivityArray.m`
- Represent endpoint/node fields as `NL x 2` arrays:
  - `hh`, `hh_prev`, `hh_frez`, `h_frez`, `TT`, `TT_CRIT`
  - `Theta_L`, `Theta_U`, `Theta_II`
  - `theta_s`, `theta_r`, `alpha`, `n`, `m`, `Ks`
- Inline or batch-compute:
  - `calculateGamma_hh`
  - `calculateTheta_m`
  - `calculateTheta_UU`
  - `calculateTheta_II`
  - `calculateTheta_LL`
  - `calculateDTheta_UUh`
  - `calculateDTheta_LLh`
  - `calculateSe`
  - `calculateKL_h`
- Remove inner-loop calls to:
  - `io.getSoilConstants`
  - `core.PhysicalConstants.instance`
  - per-node `SV` struct construction
- Keep the old scalar implementation available as a reference.
- Add a config or environment selector:
  - `HydraulicKernel=scalar`
  - `HydraulicKernel=array`

Files likely touched:

- `src/+soil/+water/+conductivity/calculateHydraulicConductivityArray.m`
- `src/+soil/+water/+conductivity/calculateHydraulicConductivity.m`
- `src/UpdateSoilWaterContent.m`
- `src/+core/ModelSettings.m`
- `data/config_main.txt`
- new A/B test script under `tests/`

Estimated implementation time: 3-5 days for a conservative first version.

Expected speed gain: large. This is the main MATLAB-layer optimization. Target is 2-4x faster for `calculateHydraulicConductivity`; whole-run runtime may drop from 7-8 h toward 3-5 h after P1.

Acceptance:

- Short-window A/B with `scalar` vs `array` reports:
  - `Theta_LL`
  - `Theta_UU`
  - `Theta_II`
  - `DTheta_LLh`
  - `DTheta_UUh`
  - `Se`
  - `KL_h`
  - `KfL_h`
  - `Ratio_ice`
- Suggested tolerances:
  - water contents: max absolute difference <= `1e-10` to `1e-8`
  - conductivities: absolute difference <= `1e-10` where values are not near zero; use relative tolerance outside tiny-K zones
  - hourly fluxes and conservation residuals unchanged within numerical noise
- Old scalar path remains callable for debugging.

Risk:

- High for numerical equivalence because the scalar helper functions contain multiple branch paths.
- Low for scientific intent if all formulas are copied exactly and A/B tests pass.

### P3 - Dall'Amico Array Kernel

Objective: remove per-node struct-returning phase-change calls.

Current bottleneck:

`calculateDallAmicoPhaseChange` takes scalar inputs and returns a struct. It is called inside the hydraulic endpoint loop whenever freeze-thaw and thermal time inputs are active.

Meaning:

`calculateDallAmicoPhaseChangeArray.m` would not be a new physical process. It
is the array form of the existing Dall'Amico/ClimaLand-style freeze-thaw
relaxation already implemented in `calculateDallAmicoPhaseChange.m`. The scalar
function currently computes one endpoint at a time and returns a struct:

- next liquid water `theta_l_next`
- next ice water `theta_i_next`
- equilibrium liquid target `theta_l_star`
- liquid phase-change source `source_liquid_s`
- thermal relaxation time `thermal_time_s`
- limiter label (`freezing-limited`, `melting-limited`, or none)

The array version should take `NL x nD` numeric arrays and return numeric arrays
for those same values. If implemented correctly, the only intended effect is
speed. The physical outputs, conservation behavior, limiter decisions, and
diagnostic fields must match the scalar path within the same tolerance used for
P2.

Other impacts to control:

- It changes only the implementation path used by `calculateHydraulicConductivityArray.m`; the default scalar reference path must remain available.
- It must fall back cleanly when `HydraulicKernel=scalar` or when array prerequisites are not met.
- It must preserve the Dall'Amico diagnostic fields written into `SoilVariables`: `DallAmicoThermalTime`, `DallAmicoThetaLStar`, `DallAmicoPhaseSource`, and `DallAmicoLimiterCode`.
- It must be tested on limiter edge cases because branch equality, not just arithmetic formula equality, is the main risk.

Work:

- Add:
  - `src/+soil/+water/+conductivity/calculateDallAmicoPhaseChangeArray.m`
- Inputs are `NL x 2` arrays:
  - `theta_l`, `theta_i`, `TT`, `theta_s`, `theta_r`, `alpha`, `n`, `m`, `rhoC`, `lambda`, `dz`, `dt`
- Outputs are arrays, not structs:
  - `theta_l_next`
  - `theta_i_next`
  - `theta_l_star`
  - `source_liquid_s`
  - `thermal_time_s`
  - `limiter_code`
- Inline helper logic:
  - `vgHeadFromSe`
  - `vgSeFromHead`
  - `boundedPhaseStep`
  - `clamp`

Files likely touched:

- `src/+soil/+water/+conductivity/calculateDallAmicoPhaseChangeArray.m`
- `src/+soil/+water/+conductivity/calculateHydraulicConductivityArray.m`
- tests from P2

Estimated implementation time: 2-4 days.

Expected speed gain: significant in freeze-thaw periods. Combined with P2, target whole-run runtime is roughly 2.5-4 h before MEX.

Acceptance:

- Scalar vs array Dall'Amico output equality on synthetic edge cases:
  - warm no-change
  - freezing-limited
  - melting-limited
  - dry near residual
  - saturated near porosity
- Full short-window A/B stays within P2 tolerances.

Risk:

- Medium-high around limiter branches and exact boundary behavior.

### P3b - Soil Matrix and Geometry Cache Cleanup (superseded by P4-c)

Objective: reduce repeated coefficient and geometry work that remains visible
after the hydraulic array kernel.

Current evidence:

The profile shows `calculateMatricCoefficients` at about `1869 s`,
`solveSoilMoistureBalance` at about `1100 s`, and repeated `Geometry1D` property
accessors (`total_soil`, `cell_volume`, `dz`, `zminus`) together above `1800 s`.
Those are not as dominant as the hydraulic kernel, but they will become more
visible after P2/P3.

Status: the post-P2/P3 two-year profile confirmed this class as the new
soil-water bottleneck (`calculateMatricCoefficients = 2948.191 s`). The active
implementation plan is now P4-c.

Work:

- Cache `cfg.Grid` geometry vectors and scalar constants once per run or per
  outer step instead of resolving object properties inside hot loops.
- Hoist terms in `calculateMatricCoefficients` that depend only on geometry,
  boundary type, or static soil parameters.
- Keep state-dependent matrix terms (`h`, `Theta`, `K`, freeze-thaw coupling)
  recomputed exactly where needed.

Files likely touched:

- `src/COMPASS.m`
- soil water matrix/coefficient assembly functions
- `src/+core/Geometry1D.m` only if call-site caching is insufficient

Estimated implementation time: 1-3 days after P2/P3.

Expected speed gain: moderate, likely several hundred seconds to about one
thousand seconds in the original one-year profiled run if repeated geometry
lookups are a real runtime cost outside profiler overhead.

Acceptance:

- Matrix coefficients match baseline for a frozen state snapshot.
- Water balance residuals and Picard convergence statistics match baseline.
- No change in accepted/rejected substep behavior except roundoff-level changes.

Risk:

- Medium. This touches numerical assembly, so it must be snapshot-tested before
  a full run.

### P4-a - Saturation Vapor Utility Fast Path

Objective: remove repeated scalar utility overhead from `Calc_Esat_Flatau`
without changing the Flatau polynomial or phase switch.

Current evidence:

The two-year profile shows `Calc_Esat_Flatau` at `1029.887 s` across `756771`
calls. Line-level profiling shows most time in the 9-term polynomial loop and
masked assignment lines. The current implementation also recreates four
coefficient vectors on every call.

Work:

- Skip `qs` and `dqs_dT` when callers only request saturation vapor pressure
  and its temperature derivative.
- Replace repeated water/ice masked gather/scatter in the polynomial loop with
  one coefficient-row selection per input point.
- Keep the original 9-term polynomial accumulation order; do not use Horner
  evaluation unless a future equivalence test accepts roundoff differences.
- Preserve current vector/scalar input behavior and output shape.
- Keep the exact same phase rule: water for `T >= 0`, ice for `T < 0`.

Files likely touched:

- `src/+equations/Calc_Esat_Flatau.m`
- focused equivalence test under `tests/` or `tools/`

Status: implemented and function-level verified.

Expected speed gain: high for a utility change; microbenchmark showed about
2.4x on the hot path, implying hundreds of seconds in the profiled two-year
run if call mix is similar.

Acceptance:

- Dense scalar/vector temperature grid from below freezing to hot canopy states.
- Old vs new `es_Pa`, `des_dT`, `qs`, `dqs_dT` match within roundoff.
- 3-day model A/B reports zero or machine-roundoff output differences.

Risk:

- Low. Formula and branch rule are unchanged.

### P4-b - Bottom Temperature Boundary Cache

Objective: remove repeated date-series and layer-mean work from
`calculateBottomTempBoundary`.

Current evidence:

The two-year profile shows `calculateBottomTempBoundary` at `1167.447 s` across
`807233` Picard calls. The function regenerates the full Julian-day series and
loops over all layers for two means on every call, although only one hourly
bottom-boundary value is written.

Work:

- Cache `jd_series`, `JD0`, and hour-to-day index mapping for the run.
- Replace the layer loop with `mean(Lambda_eff,2)` and `mean(c_unsat,2)`.
- Avoid recomputing invariant date terms inside the Picard loop.
- Keep the damping-depth calculation tied to current effective thermal
  properties; do not freeze state-dependent `Lambda_effAvg` or `c_unsatAvg`
  unless A/B proves it is safe.

Files likely touched:

- `src/+soil/+thermal/calculateBottomTempBoundary.m`
- possibly a small cache helper or fields in `cfg.Time`

Status: implemented and model A/B verified.

Expected speed gain: moderate in true wall time but high for the targeted loop.
The date cache is small; the main gain is replacing the repeated layer-mean
loop. This also removes most of the global `mean` hotspot seen in the
two-year profile.

Acceptance:

- Bottom temperature time series `Tbh` matches baseline within roundoff.
- `QET/QEB`, soil temperature profile, and energy residual match A/B tolerances.
- No change in rollback count beyond roundoff-level solver-path noise.

Risk:

- Low to medium. Date caching is safe; caching state-dependent thermal means is
  not safe unless explicitly verified.

### P4-c - Soil-Water Matric Coefficient Array Pass

Objective: remove the new post-P2 soil-water coefficient bottleneck.

Current evidence:

The two-year profile shows `soil.water.calculateMatricCoefficients` at
`2948.191 s` and `solveSoilMoistureBalance` at `1759.974 s`. The hot path is the
`NL x nD` nested loop that recomputes water/thermal coupling coefficients every
Picard call.

What P4-c mainly accelerates:

- The coefficient assembly before the water-balance linear solve, not the
  tridiagonal solve itself.
- Per-Picard endpoint work over `NL x nD` (`128 x 2` in the current grid):
  freeze/thermal derivative correction, derivative limiting, vapor/liquid
  coupling terms, and optional dry-air coupling terms.
- Repeated scalar writes to the coefficient fields:
  `Chh`, `ChT`, `Khh`, `KhT`, `Kha`, `Vvh`, `VvT`, and `Chg`.
- Repeated updates to stateful derivative fields:
  `SoilVariables.DTheta_UUh`, `SoilVariables.DTheta_LLh`, and
  `TransportCoefficient.DTheta_LLT`.
- Repeated scalar calls to `limitDerivative`, which appears as a large global
  profiler hotspot because it is invoked for every endpoint in every Picard
  call.

Expected implementation shape:

- Keep the current function as the authoritative physics location unless the
  diff becomes clearer with a local helper; do not create a parallel physics
  fork unless a selector/reference path is needed for A/B.
- Extract node-aligned vectors once, form `NL x nD` endpoint arrays, and apply
  mutually exclusive masks for the `lhThmrl` and non-`lhThmrl` branches.
- Implement derivative limiting as array clamps with the same non-finite to
  zero behavior as `limitDerivative`.
- Preserve the subtle indexing relation `MN = i + j - 1`: endpoint columns share
  neighboring nodes, so array reshaping must be audited against the scalar loop.
- Do not freeze any state-dependent terms. `hh`, `h`, `TT`, `T`, vapor density
  derivatives, freeze heads, water contents, and conductivities remain current
  Picard-state inputs.

Work:

- Convert the endpoint loop to array operations where branch masks are
  mutually exclusive and formula order can be preserved.
- Hoist constants and `cfg.Grid` values outside the endpoint loop.
- Keep state-dependent quantities (`hh`, `h`, `TT`, `T`, freeze-thaw
  derivatives, `Theta_*`, `K*`) recomputed each Picard iteration.
- Preserve `limitDerivative` semantics, including exact cap values.

Files likely touched:

- `src/+soil/+water/calculateMatricCoefficients.m`
- snapshot/equivalence harness for coefficient arrays

Estimated implementation time: 2-4 days.

Expected speed gain: large. This is the next real kernel-level optimization
after P2/P3. In the two-year profile it targets about `2948 s` directly and may
also reduce part of the `solveSoilMoistureBalance` apparent cost because that
solve is fed by these coefficient arrays.

Acceptance:

- Frozen-state coefficient snapshot A/B for `Chh`, `ChT`, `Khh`, `KhT`, `Kha`,
  `Vvh`, `VvT`, `Chg`, and updated derivative fields.
- 3-day and 10-day model A/B within the P2/P3 tolerance policy.
- Picard iteration and rollback statistics are unchanged or explainably
  roundoff-level different.
- Fresh post-P4-c profile confirms whether the next dominant target is canopy
  work, thermal/dry-air refresh, or plant-root conductance.

Risk:

- Medium-high. This touches numerical assembly, so it needs the same discipline
  as the hydraulic-array work.

Decision after P4-c:

- If `calculateMatricCoefficients` is no longer dominant, move to P4-d
  (`Plantrootconductivity`) or P4-e/P5b (canopy biochemistry/radiation),
  depending on the fresh profile.
- If `solveThermalEquations`, `calculateEnergyParameters`, or dry-air parameter
  refresh rises to the top, add a thermal/dry-air array pass before touching
  canopy physics.
- If the profile becomes dominated by one compact numerical kernel after P4-c,
  consider P7 MEX only then, with the verified MATLAB array implementation as
  the reference.

### P4-d - Plant Root Conductance Endpoint Vectorization

Objective: reduce scalar endpoint work inside ebal plant hydraulics.

Current evidence:

The two-year profile shows the plant-hydraulic chain around `1400 s` combined,
with `Plantrootconductivity` at `360.926 s` and about `38.7M` endpoint visits
through `rootSideConductance`, `rootSoilContactFactor`, and related helpers.

Work:

- Vectorize root-soil endpoint conductance calculations across layers/endpoints.
- Precompute static root geometry terms such as layer thickness, root density,
  root area, and radial geometry factors.
- Keep the same conductance mode behavior (`clm5`, `radial_darcy`,
  `radial_darcy_series`) and diagnostics.

Files likely touched:

- `src/+plant/+PlantHydraulic/Plantrootconductivity.m`
- related helper tests or snapshot harness

Estimated implementation time: 1-3 days after P4-c.

Expected speed gain: moderate; useful because it sits inside the ebal loop.

Acceptance:

- Root conductance vector, contact diagnostics, and plant hydraulic outputs
  match baseline for representative wet/dry/rootless states.
- 3-day model A/B keeps `Tp_t`, `psil`, `wt`, `GPP`, and ET within tolerance.

Risk:

- Medium. It affects hydraulic stress, which feeds stomata and GPP.

### P4-e - Biochemical and ebal Call Reduction

Objective: reduce repeated biochemical work without relying on previous-hour
warm starts.

Current evidence:

The two-year profile shows `biochemical` at `1559.995 s`, plus
`compute_assimilation`, `medlyn_gs`, and output-packing costs. The largest ebal
lines are repeated sunlit/shaded calls to `vegetation.biochemical` before and inside
the hydraulic inner loop.

Work:

- Identify fields in `biochem_in` that are unchanged within a hydraulic inner
  loop and avoid repacking them.
- Add a strict no-fluorescence or reduced-output fast path only if the requested
  diagnostics do not use fluorescence fields.
- Avoid changing Ci iteration tolerances or the photosynthesis/stomatal
  coupling equations unless a separate science decision is made.
- Do not implement previous-hour ebal warm starts; P5-0 measured no useful
  iteration reduction.

Files likely touched:

- `src/ebal_clm45.m`
- `src/+plant/biochemical.m`

Estimated implementation time: 2-5 days after P4-c/P4-d.

Expected speed gain: moderate to high, but riskier than P4-a/P4-b.

Acceptance:

- GPP, leaf temperature, `Tp_t`, stomatal resistance, and fluorescence
  diagnostics match baseline when fluorescence output is enabled.
- If a reduced-output mode is added, it must be opt-in and clearly documented.

Risk:

- Medium-high. This path directly controls GPP, transpiration, and energy
  balance.

### P4-z - Adaptive Substep Initial Guess and Rollback Reduction

Objective: reduce rejected substeps only after the cost per Picard iteration has
been reduced.

Current evidence:

The two-year profile has `3763` rollbacks and no fatal failures. This is worth
studying, but the current top costs are per-Picard function costs, not only
failed 3600 s trial steps.

Work:

- Carry a rolling `preferred_sub_dt` across outer hours.
- If recent hours had thermal or water non-convergence, start the next outer
  hour at the last stable substep size.
- Relax back to 3600 s only after several low-iteration accepted substeps.
- Add diagnostics for rejected substep count, accepted substep count, initial
  trial dt, and final accepted dt.

Files likely touched:

- `src/COMPASS.m`
- `src/+solver/initHourlySubstepAccum.m`
- `src/+solver/accumulateAcceptedSubstep.m`
- diagnostics registration if new fields are exported

Estimated implementation time: 1-3 days after P4-c.

Expected speed gain: profile-window dependent.

Acceptance:

- Fewer rejected substeps in the same difficult-window benchmark.
- No increase in outer-step fatal rollback.
- Water/energy conservation residuals unchanged or improved.

Risk:

- Medium. It changes numerical path but not equations.

### P5 - ebal and Biochemistry Warm Starts

Objective: reduce repeated canopy energy/photosynthesis iterations.

Current behavior:

`ebal_clm45` initializes leaf/canopy states from current forcing each hour and repeatedly calls `vegetation.biochemical`. Hydraulic inner loops may re-run biochemistry when water stress changes.

Measurement result (P5-0, 2026-06-16):

P5-0 was implemented as opt-in, write-only diagnostics and then used only to
measure whether previous-hour canopy states predict ebal iteration count. The
diagnostic scaffold was removed after the decision because no production path
uses it.

Summary:

| Window | Daytime hours | Mean ebal iter | corr(cold-distance, iter) | Interpretation |
|---|---:|---:|---:|---|
| April | 386 | 4.73 | 0.305 | weak positive relation during leaf-out |
| May | 434 | 5.38 | 0.027 | no useful relation |
| June | 420 | 5.92 | -0.075 | no useful relation |
| July 1-10 | 140 | 5.54 | -0.072 | no useful relation |
| All daytime valid-warm hours | 1736 | - | -0.013 | no useful relation |

Warm-start states were often much closer to the converged state than the cold
initial guess (`Tcu` distance about 10.6 -> 2.7 degC, `Tch` about 3.0 -> 1.5
degC, `Ci` warm distance about 6 ppm), but ebal still required about 5-6
iterations. Hours where warm was more than 5 degC closer still averaged about
6.2 iterations. Soil surface temperature already carries previous-hour state
through `soil.Ts`, so it had little remaining warm-start opportunity.

Decision:

- Do not proceed to P5-1/P5-2.
- Do not add `prev_canopy`, reset gates, or an ebal warm-start config switch.
- The iteration floor appears controlled by the coupled energy/photosynthesis/
  stomatal/hydraulic solve sequence, not by the hour-level initial distance.
- Continue with P5b-style cache work or the next post-profile hotspot instead.

Risk avoided:

- Warm-starting ebal would change nonlinear solver initial conditions and add
  cross-hour state management, with near-zero measured speed benefit.

### P5b - Radiation Cache Work

Objective: reduce repeated canopy-radiation utility work only after the new
soil-water and canopy-hydraulic hotspots are handled.

Current evidence:

The post-P2/P3 two-year profile shows `RTMf` at `558.406 s`, `RTMo` at
`268.163 s`, `RTMz` at `254.020 s`, and `RTMt_sb` at `54.195 s`. These costs are
real but are smaller than `ebal_clm45`, water matrix assembly, biochemistry,
bottom-temperature boundary, and `Calc_Esat_Flatau`.

Work:

- Cache radiation results only when forcing, solar geometry, canopy structure,
  and spectral inputs are unchanged for the needed interval.
- Cache angle/spectral helper products that are purely geometric.
- Avoid caching across iterations when leaf/canopy/soil temperature, LAI,
  wetness, or optical state has changed.

Files likely touched:

- `src/ebal_clm45.m`
- radiation transfer functions

Estimated implementation time: 2-5 days after P4-c/P4-e.

Expected speed gain: moderate and profile-dependent.

Acceptance:

- Radiation diagnostics, GPP, ET, leaf temperature, and energy balance residual
  remain within configured tolerances.
- Cache hit diagnostics prove the cache is actually used.

Risk:

- Medium-high. Incorrect reuse can silently stale radiation or vapor-pressure
  state.

### P6 - Nitrogen Loop Cleanup

Objective: reduce repeated N reaction/transport work after water/thermal bottlenecks are addressed.

Current status:

In the original profile `solveNitrogen` was about `2050 s`. In the post-P2/P3
two-year profile it is `141.143 s` across `17520` calls, so it is not a near-term
runtime priority. Keep this section as a later cleanup item unless a future
profile after P4 shows N rising again.

Work:

- Identify terms constant over the N Picard loop and hoist them out.
- Cache temperature response factors and pH terms when state changes are negligible.
- Avoid repeated O2 transport where fallback or slow-changing state is valid.
- Keep conservation checks strict.

Files likely touched:

- `src/+soil/+nitrogen/solveNitrogen.m`
- `src/+soil/+nitrogen/+chemistry/*`
- `src/+soil/+nitrogen/+transport/*`
- `src/+soil/+gas/solveO2Transport.m`

Estimated implementation time: 3-5 days after P2-P4.

Expected speed gain: moderate; likely 10-25% of N module time.

Acceptance:

- NH4/NO3 cumulative residuals unchanged.
- N2O/N2, mineralization, nitrification, denitrification remain numerically consistent.
- Plant N uptake and plant N residual unchanged within numerical tolerance.

Risk:

- Medium-high because conservation accounting is sensitive.

### P6b - Organic Matter Tick and Turnover Cleanup

Objective: reduce organic source/turnover overhead after hydraulic and canopy
work no longer dominate.

Current evidence:

The profile shows `computeOrganicSources` around `1865 s`,
`OrganicStandard.tick` around `1269 s`, and organic/soil C/N pool methods such
as `turnover_pool`, `total_C`, `total_N`, and `soil_C` contributing hundreds of
seconds. This is related to nitrogen/C coupling but is not the first bottleneck.

Work:

- Cache pool totals within a tick when no pool update has occurred between
  reads.
- Hoist temperature/moisture modifiers that are reused across organic matter
  and nitrogen source calculations.
- Keep process-level cumulative diagnostics unchanged.

Files likely touched:

- organic matter pool classes and source functions
- `src/+soil/+nitrogen/solveNitrogen.m`

Estimated implementation time: 2-4 days after P6.

Expected speed gain: moderate.

Acceptance:

- Total C/N pool sums and process cumulative fluxes match baseline.
- Plant C/N residuals and NH4/NO3 residuals remain unchanged within tolerance.

Risk:

- Medium-high because hidden ordering assumptions in pool updates can affect
  conservation.

### P7 - MEX Compilation of the Remaining Kernel

Objective: reach the less-than-1-hour target only after the post-P4 MATLAB
profile identifies a remaining compact numerical kernel worth compiling.

Work:

- Use the verified MATLAB implementation as the MEX source target.
- Candidate targets are no longer automatically hydraulic/freeze-thaw. Re-rank
  after P4; likely candidates may include water matrix coefficients, plant-root
  conductance, or biochemical inner kernels.
- Keep runtime fallback:
  - use MEX only when present and enabled;
  - otherwise use the MATLAB reference path;
  - scalar/debug paths remain callable where they exist.
- Avoid passing large MATLAB objects to MEX:
  - pass numeric arrays and scalar config values only.

Files likely touched:

- whichever post-P4 function remains dominant
- MEX build script under `tools/`

Estimated implementation time: 4-8 days after P4 and a fresh profile.

Expected speed gain: profile-dependent. Do not compile the hydraulic kernel just
because it used to be the bottleneck; the two-year profile shows it has already
fallen behind other modules.

Acceptance:

- MEX vs MATLAB short-window output differences pass the relevant module
  tolerances.
- MEX fallback works if the compiled file is missing.
- Build instructions are reproducible on the target Windows/MATLAB installation.

Risk:

- High engineering risk, low scientific risk if array MATLAB reference remains authoritative.

## Suggested Timeline

| Stage | Work | Time | Runtime target after stage |
|---|---|---:|---:|
| 0 | P0 benchmark harness | done | completed |
| 1 | P1/P1.5 registry, NetCDF, CSV/partial-output quick wins | done | completed |
| 2 | P2 hydraulic MATLAB array kernel | done | completed; hydraulic kernel no longer dominant |
| 3 | P3 Dall'Amico array kernel | done | completed; P2/P3 verified by 3-day, 10-day, and 1-year A/B |
| 4a | `Calc_Esat_Flatau` fast path | done; pending commit grouping | utility hotspot reduced |
| 4b | `calculateBottomTempBoundary` cache/vectorization | done; pending commit grouping | bottom-boundary loop and global mean hotspot reduced |
| 4c | water `calculateMatricCoefficients` array pass | 2-4 days | next major kernel speedup; targets coefficient assembly |
| 4d | `Plantrootconductivity` endpoint vectorization | 1-3 days | moderate canopy-hydraulic speedup |
| 4e | biochemical/ebal call reduction | 2-5 days | moderate/high but higher risk |
| 4z | adaptive substep/rollback strategy | 1-3 days | profile-window dependent, after P4c |
| 5 | ebal previous-hour warm starts | measured and rejected | do not implement |
| 5b | radiation cache work | 2-5 days | defer until after P4 |
| 6 | N cleanup | 3-5 days if future profile justifies it | currently low priority |
| 6b | organic matter tick and turnover cleanup | 2-4 days if future profile justifies it | currently low priority |
| 7 | MEX remaining dominant kernel | 4-8 days after fresh profile | target <1 h if MATLAB P4 is insufficient |

These are engineering estimates for careful implementation with A/B verification, not just coding time.

After P4-c, do not choose the next stage from the old two-year profile alone.
Run a fresh short profile because P4-a/P4-b/P4-c remove several of the current
top entries. The expected remaining candidates are:

- Canopy block: `ebal_clm45`, `vegetation.biochemical`, `heatfluxes_clm45`, and
  radiation transfer. This is likely the biggest remaining science-sensitive
  class, but previous-hour ebal warm starts have already been rejected.
- Plant hydraulics: `Plantrootconductivity`, `rootSideConductance`, and related
  endpoint helpers. This is a contained vectorization target inside the canopy
  loop.
- Thermal and dry-air refresh: `solveThermalEquations`,
  `calculateEnergyParameters`, `calculateDryAirParameters`, and gas transport
  helpers. These are good candidates if the fresh profile shows soil-side
  parameter refresh rising after P4-c.
- Residual hydraulic/freeze-thaw work: currently no longer first priority, but
  it can become a MEX candidate only if it reappears as a compact post-P4
  hotspot.
- Rollbacks/adaptive substeps: optimize only after per-Picard costs are reduced,
  because the present top costs are still repeated work inside accepted Picard
  solves.

## Verification Policy

Every performance change must report both speed and numerical equivalence.

Required checks:

- 3-day smoke benchmark
- 10-day difficult-window benchmark, preferably including freeze-thaw or thermal non-convergence windows
- 1-year final benchmark after major stages

Required numerical outputs:

- `Theta_LL`, `Theta_UU`, `Theta_II`
- `DTheta_LLh`, `DTheta_UUh`
- `KL_h`, `KfL_h`, `Ratio_ice`
- `QMT`, `QMB`, `QET`, `QEB`
- `GPP`, `Reco`, `NEE`
- `water_balance_residual`
- `energy_balance_residual`
- plant C/N residuals
- NH4/NO3 cumulative residuals

Required runtime outputs:

- total runtime
- time in hydraulic kernel
- time in Dall'Amico kernel
- Picard iteration count distribution
- accepted/rejected substep counts
- ebal iteration count distribution
- output flush time

## Implementation Rule

Do not jump directly from the current scalar implementation to MEX. The safe path is:

1. scalar baseline
2. MATLAB array implementation
3. scalar-vs-array A/B equivalence
4. MEX generated from the array implementation
5. MATLAB-array-vs-MEX A/B equivalence
6. one-year performance run

This keeps the scientific reference path readable in MATLAB while allowing the production path to be compiled.
