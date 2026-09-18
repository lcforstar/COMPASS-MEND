# Water-only Richards capability roadmap

Date: 2026-07-12
Scope: one-dimensional water-only Richards simulations under natural and
laboratory boundary conditions  
Current milestone: W1 boundary contract plus conservative variably saturated
core are complete; W2 benchmark matrix is in progress

Coupled-physics hold: before any I6 water-vapor-air-heat wiring, complete the
source-equation audit and staged validation protocol in
`coupled_water_vapor_heat_equation_audit.md` and
`coupled_water_vapor_heat_validation_protocol.md`. W2 remains a deliberately
water-only numerical and boundary benchmark; coupled coefficient repairs must
not be tuned inside W2 cases.

## 1. Target capability

The water-only configuration should become a standalone soil-water model with
the following properties:

1. the same Richards kernel is used by natural land-surface simulations and
   controlled soil-column experiments;
2. top and bottom boundaries support natural forcing, prescribed heads,
   prescribed fluxes, closed boundaries, free drainage, and explicit time
   series;
3. nonlinear convergence, adaptive substeps, rollback, and water accounting
   remain valid under every supported boundary;
4. root uptake and internal sources are physically limited and accounted for;
5. standard cases can be compared against analytical solutions, HYDRUS, or
   another high-resolution reference;
6. the architecture remains directly portable to plain-data Julia/SciML
   components.

The comparison has two distinct benchmarks:

- ClimaLand is the architecture benchmark: independent policy, driver, state,
  request, solver, and diagnostic ownership.
- HYDRUS is the numerical and physical benchmark: Richards convergence,
  experimental boundaries, atmospheric switching, root uptake limitation,
  water conservation, and standard benchmark behaviour.

The goal is capability parity where it is scientifically relevant, not a
line-by-line copy of either model.

### 1.1 Spatial execution and point-context hierarchy

COMPASS does not treat "ecosystem" and "land surface" as separate physical
models. The same coupled column physics is used at a field site and inside a
regional or global calculation. The first distinction is spatial execution;
the second distinction, for a point column, is the input context:

```text
spatial execution
|-- point
|   |-- field: natural atmosphere, canopy, snow, management, and root uptake
|   `-- laboratory: prescribed boundary protocol and optional internal sources
`-- regional/global
    `-- many instances of the same accepted point-column kernel
```

The intended user-facing selectors are:

```toml
[simulation]
spatial_mode = "point"   # point | regional
point_context = "field"  # field | laboratory; used only for point
```

These selectors describe execution and driver ownership. They do not replace
the existing process switches or create a second Richards implementation.
Their runtime wiring is a later execution-layer increment; W1 keeps the
laboratory and field entry points explicit while sharing the scheduler and
kernel.

### 1.2 Reference-model synthesis

The architecture combines complementary lessons from four local reference
models:

- CoLM: the run loop reads regional forcing, dispatches patches through
  `CoLMDRIVER`, calls one `CoLMMAIN` column kernel per patch, and applies
  lateral flow after vertical-column physics;
- Noah-MP: the offline driver loops over grid cells and calls one
  `NOAHMP_SFLX` kernel per cell, demonstrating a clear grid-executor versus
  point-kernel boundary;
- CLM5: the natural driver orders begin-balance, canopy hydrology, surface
  fluxes, soil temperature, hydrology/drainage, ecosystem processes, and final
  balance checks;
- ClimaLand: prescribed laboratory boundary types and atmosphere-driven land
  boundaries use the same component/timestep framework, while coupled
  infiltration is owned by the soil-surface composition rather than duplicated
  in each component.

COMPASS retains its own stronger accepted/rejected trial transaction. The
resulting target is:

```text
PointExecutor
|-- FieldPointDriver ---------\
`-- LaboratoryDriver ----------+--> PointColumnKernel

RegionalExecutor --> cell/block views --> PointColumnKernel --> regional assembly
```

The regional executor must call a plain-data point kernel, not the current
full `runOneStep` with global clock, output, checkpoint, and diagnostic
registry ownership.

## 2. Accepted baseline that must not be redesigned

Current-step root-water conservation is closed. Uptake ownership was then
separated before W2 so water, nitrogen, and PlantCN no longer share one commit
adapter.

Verified on the current new chain:

- real water-only trajectory reached growing-season `kt=2335`;
- heat-accumulation LAI was `0.3301`, above `ALAI_MIN=0.293`;
- `Rnctot=2.48107 W m-2`;
- `Tp_t=3.70364217206e-8 cm s-1`;
- finite-element `integral(Srt)=Tp_t`;
- WBC transpiration equals both `Tp_t*dt*10` and
  `integral(Srt)*dt*10`;
- one-tile, 0.4/0.4/0.2 multi-tile and reject/no-leak paths pass.

The accepted single-source uptake policy is:

1. `root_water_uptake` maps to `lRwuef` and solely controls current `Srt`;
2. `lRwuef=1` maps current ebal `RWU` directly into
   `workspace.scope_step.Srt`;
3. `lRwuef=0` stages and commits exact-zero `RWU/Srt`, clearing stale uptake;
4. `root_nitrogen_uptake` maps to `lRsuef` and solely controls current
   `SinkS` and `SPUptake`;
5. `plant_cn_cycle` maps to `lPlantCN`; PlantCN may provide a demand cap and
   consume `U_N_total`, but it cannot calculate, change, or commit `Srt` or
   `SinkS`;
6. with PlantCN demand available, root N uses the demand-capped DAISY+Nacry
   path; without PlantCN, it falls back to the original passive
   `max(Srt,0)*max(C,0)` mass-flow path, with no demand limit.

`lActRSU` is a deprecated legacy compatibility field. Active membrane uptake
is intrinsic to the unified M1 DAISY+Nacry kernel; there is no second active
uptake switch in the new chain.

Cold, low-LAI ebal non-convergence warnings are a separate canopy energy
balance issue. They are not evidence of an Srt or Richards wiring failure.

## 3. Capability gaps relative to HYDRUS-class use

| ID | Gap | Current state | Planned milestone |
|---|---|---|---|
| H1 | Clear head/flux input contract | Completed in W1 | W1 |
| H2 | Natural and prescribed providers | Completed in W1 | W1 |
| H3 | Static laboratory boundary validation | Bonan fixed-head/free-drainage case passes; remaining cases are open | W2 |
| H4 | HYDRUS-like nonlinear convergence contract | Completed for the new water path: mixed storage, TolTh/TolH, active-set stability, and an independent mass gate | W1c |
| H5 | Actual-flux water ledger | Completed for the new water path using QMT, QMB, integrated Srt, and Theta_UU plus vapor | W1c |
| H6 | Atmospheric complementarity and surface water | Completed: finite DSTOR, two-pair active set, accepted-only runoff/storage, and boundary residual gate | W5 |
| H7 | Root/source availability limiting | Current Srt wiring is correct; capacity limits are missing | W6 |
| H8 | Advanced lower boundaries and benchmark matrix | Not yet implemented | W7-W8 |

Important current limitations:

- atmospheric surface water is complete for liquid 1-D ponding; retained
  frozen ponding still fails explicitly until surface ice/enthalpy is modeled;
- seepage face, dynamic groundwater, and lysimeter-style lower boundaries are
  not yet complete;
- the first W2 reference benchmark is complete; the remaining static boundary
  and analytical cases are still open.

### 3.1 W1c conservative variably saturated increment

The new transactional water path now retains signed positive pressure head and
uses one active-set Richards solve across h=0. For h>=0, condensed storage is
`theta_s + Ss*h`, capacity is `Ss`, and warm conductivity is `Ks`. Mixed-form
storage and WBC use `Theta_UU + Theta_V*rho_v/rho_l`; accepted ice is owned by
`Y.Theta_I`. Positive-head `Ss*h` remains in total storage but outside the
physical liquid/ice pore-volume partition, so the later heat equation cannot
mistake elastic storage for sensible or latent phase volume. The actual-flux ledger is
`-QMT - integral(Srt) - (-QMB) - dStorage` after unit conversion.

The laboratory driver now also accepts isothermal freeze-thaw Richards cases:
temperature is prescribed, `L_f` and phase storage are active, while the heat
and gas equations remain disabled. A real frozen 60 s flux step closes below
`1e-6 mm`; the observed regression residual is about `-1.16e-9 mm`.

All h=0 crossings are diagnostic outputs, including multiple crossings from
perched saturation. See `saturated_flow_route_B_design.md` for the implemented
algorithm and validation evidence.

## 4. Ownership model

### 4.1 BoundarySpec

`built.boundary.water` owns immutable configuration only:

- endpoint kind;
- source kind;
- constant value or series metadata;
- units and public sign convention.

It does not own forcing values, observations, trial state, actual flux, or
surface-water storage.

### 4.2 Drivers

Natural and prescribed data are separate from the boundary policy:

- natural surface payload: current `Throughfall` and `EVAP`;
- prescribed laboratory data: `built.forcing.water_boundary` loaded from
  optional `water_boundary.csv`;
- irrigation schedules: management sources, not boundary policy;
- `h_sur_cm`: observation compatibility field, not an automatic control.

### 4.3 BoundaryRequest

The scheduler materializes one immutable request for each exact trial
interval. It contains explicit start time, end time, duration, top request,
bottom request, and sign convention.

The provider is pure:

- no cursor is advanced;
- no model clock is updated;
- no prognostic state is changed;
- retrying the same physical interval produces the same request;
- changing `dt` re-evaluates the new interval exactly.

### 4.4 BoundaryResult

BoundaryResult is audit output, not state. It separates:

- requested kind and value;
- applied mathematical kind;
- actual solved top and bottom flux;
- accepted or rejected trial status.

For an atmospheric request, the applied top kind may be `fixed_flux`,
`surface_storage`, or `fixed_head`. W5 resolves those active sets inside the
nonlinear trial and commits only a mass- and complementarity-consistent result.

### 4.5 Prognostic state

`Y.surface_water` owns surface-water storage. It must change only on accepted
substeps. Boundary policy and observations must never own or overwrite this
state.

### 4.6 Internal sources

Root uptake and future subsurface irrigation are internal source terms. They
must not be inserted into top or bottom boundary fluxes.

## 5. Single user-facing site schema

New sites keep all water settings in the existing `[boundary]` table:

```toml
[boundary]
water_top_kind    = "atmospheric"
water_bottom_kind = "free_drainage"
```

Historical integer meanings remain only in the internal mapping that generates
solver aliases. Every prescribed mode reads the fixed top_value or bottom_value
CSV column. Constant modes require an unchanged column, variable_head uses
linear, and variable_flux uses previous. Historical integer fields, the split
kind/source/column/interpolation schema, and the earlier nested draft must be
migrated and are not accepted as alternative schemas.

The public flux sign is always positive into the soil domain. Legacy top and
bottom signs are converted only at the solver adapter.

See `docs/water_boundary_site_input.md` for all supported inputs.

## 6. Milestones

### W0: preserve the accepted natural baseline

Keep current Srt, independent SinkS, state-layer, multi-tile, reject, and natural-boundary
regressions green before and after every later milestone.

### W1: boundary contract and experiment driver — complete

Completed work:

- canonical typed top and bottom boundary schema;
- strict separation of head and flux fields;
- single public `flux_into_domain_cm_s` sign convention;
- generated legacy solver aliases from the single typed site schema;
- pure materialization for constants, natural surface input, and prescribed
  series;
- exact interval integration for piecewise-constant flux series;
- endpoint evaluation for previous/linear head series;
- optional `water_boundary.csv` input path;
- per-tile and per-retry BoundaryRequest materialization;
- shared `solver.advanceWaterOnlyOuterStep` for tiling, retry, rollback,
  mid-tile carry, and accepted-flux accumulation;
- natural `runOneStep` and standalone experiment driver use the same
  `solver.attemptSolveWaterOnly`;
- standalone experiments bypass RTMo, ebal, and PlantCN without copying the
  Richards solver or scheduler;
- solver adapter support for atmospheric, fixed head, fixed flux, zero flux,
  and free drainage;
- separate requested, applied, and actual BoundaryResult fields;
- Changwu migrated to the single typed site schema while preserving old
  `water_top=3/water_bottom=3` physics.

W1 responsibility refinement:

- the field entry owns natural `Throughfall` and `EVAP` production;
- the laboratory entry accepts prescribed boundary values only from the single
  configured `water_boundary.csv` contract;
- the laboratory entry rejects an atmospheric top policy instead of injecting
  fake zero natural forcing;
- both entries share `solver.advanceWaterOnlyOuterStep` and
  `solver.attemptSolveWaterOnly`;
- a prescribed negative top flux is an experimental outflow, not a natural
  evaporation demand.

W1 intentionally does not include atmospheric complementarity redesign,
surface-water evolution, seepage face, convergence-policy changes, mixed-form
default selection, root uptake limiting, or dry-air I6 work.

W1 initialization and execution closure completed before W2:

- `[simulation] spatial_mode` selects `point` or `regional`, and
  `point_context` selects `field` or `laboratory` for a point calculation;
- `init.resolveSimulationPlan` resolves active inputs from the explicit context
  and module switches before any optional CSV or SCOPE library is read;
- `data/sites/lab_water` is a minimal standalone example containing only
  `site.toml`, `soil.csv`, `soil_initial_profile.csv`, and `water_boundary.csv`;
- the laboratory path does not read or build meteorological forcing, roots,
  plant pools, vegetation observations, canopy/SCOPE state, field diagnostics,
  checkpoints, or management services;
- false switches permit the corresponding laboratory input groups to be
  absent, while a true unsupported laboratory component fails explicitly;
- `run_site("changwu")` selects the field loop and `run_site("lab_water")`
  selects the standalone Richards loop from one repository-root entry;
- both loops retain the accepted W1 boundary contract and call the same
  Richards trial kernel. No W2 convergence or boundary physics is included in
  this initialization increment.
- field uptake ownership is independent of PlantCN: `lRwuef` owns `Srt`,
  `lRsuef` owns `SinkS`, and both fluxes use the same accepted-step transaction
  so rejected steps cannot leak them;
- `lRsuef=0` returns exact-zero root-N sinks; with PlantCN disabled and
  `Y.plant_cn` absent, `lRsuef=1` runs passive water-linked NH4/NO3 uptake;
  with PlantCN enabled, its CN-based demand selects the DAISY+Nacry path.

### W2: static experimental boundary correctness

Validate every supported static endpoint with an analytical solution,
conservation identity, or high-resolution reference:

- fixed top and bottom head;
- fixed top and bottom flux;
- zero flux;
- free drainage;
- head and flux time series;
- closed column;
- saturated Darcy flow;
- constant-flux infiltration;
- horizontal infiltration;
- capillary rise;
- layered material interface.

Required gates:

- grid convergence;
- time-step convergence;
- request and applied condition match the intended experiment;
- actual solved flux is recorded independently;
- no retry or tiling dependence in prescribed water amount.

Current W2 evidence (updated 2026-07-12):

- the homogeneous Bonan clay fixed-state/free-drainage case now has fresh
  COMPASS, ClimaLand, and HYDRUS runs;
- COMPASS uses positive named phase derivatives in the separate conservative
  Richards assembler, mixed-form storage, no Anderson,
  `MaxIt=20`, `TolTh=1e-4`, and `TolH=0.1 cm`;
- the negative-capacity source, public-code masking path, superseded `abs`
  compatibility patch, and final conservative boundary are recorded in
  `docs/w2_vg_capacity_sign_audit.md` and
  `docs/water_heat_freeze_conservative_derivation.md`;
- all 1,000 COMPASS steps were accepted; after the unified total-water ledger,
  the rerun cumulative balance residual is `-0.000365 mm`; the previously
  reported final COMPASS-HYDRUS RMSE was `0.001442 m3/m3`;
- the current ClimaLand source test against its Bonan artifact fails at RMSE
  0.027383, so that artifact is retained as historical evidence rather than a
  golden gate;
- the sourced SV62 `+5 cm` head/free-drainage case now preserves positive
  pressure head, all `h=0` crossings, augmented storage, actual boundary
  fluxes, and full pressure-head profiles. The ClimaLand sampled-depth policy
  is now a canonical site option rather than a benchmark override:
  `[soil_input].vertical_parameter_mode="sample_linear"` interpolates primitive
  `soil.csv` parameters to element centers and holds endpoint values. The fresh
  60 min COMPASS ledger is: initial storage `38.877015 mm`, final storage
  `350.152487 mm`, top input `311.283648 mm`, bottom output `0.008112 mm`, and
  residual `6.32e-5 mm`;
- a strict sensitivity run (`TolTh=1e-5`, mass tolerance `1e-7 mm`, MaxIt=80)
  accepted all 120 steps without retry and reduced the ledger residual to
  `5.63e-6 mm`, but the 60 min COMPASS-ClimaLand augmented-profile RMSE stayed
  `0.03021 m3/m3`. The remaining mismatch is therefore not a convergence-
  tolerance or water-loss defect; it remains a W2 spatial/interface-flux
  discretization audit;
- the accepted growing-season Srt chain was rerun after the storage changes:
  at `kt=2335`, `LAI=0.3301`, `Tp_t=3.68385204645e-8 cm/s`, and the FE
  integral of Srt is identical. Single-tile, three-tile, reject/no-leak, and
  `lRwuef=0` paths still pass;
- the complete sourced-case and execution matrix is maintained in
  `docs/w2_boundary_benchmark_matrix.md`;
- W2 remains open for closed-column, saturated Darcy, prescribed flux,
  capillary-rise, horizontal, variable-series, evaporation, and layered cases.

### W3: nonlinear convergence contract

- use damped `delta thetaEq` in unsaturated regions;
- use damped `delta h` near saturation;
- keep active boundary mode fixed within one residual/Jacobian evaluation;
- keep conservative mixed storage fixed as the only time formulation, and
  expose `TolTheta` and `TolH` as production controls;
- complete mass-led acceptance auditing;
- decide the mixed-form default only after W2 passes.

### W4: actual-flux water accounting

Unify the ledger around:

- actual top flux;
- actual bottom flux;
- actual `integral(Srt)`;
- surface-water storage change;
- runoff;
- soil storage change;
- absolute and relative residual.

Requested precipitation, evaporation, or infiltration must not replace actual
solved boundary flux in the water balance.

### W5: atmospheric complementarity and surface water

**Implemented 2026-07-17.** The 1-D column keeps the HYDRUS/STEMMUS
finite-depth surface reservoir; it does not adopt a mixed finite-element or
Nitsche discretization. For every coupled Picard map it now enforces

\[
d^{n+1}-d^n=\Delta t\left(P_l-E_s+Q_{MT}-R\right),
\]

with `QMT` positive outward. The complete finite-storage constraints are

\[
0\le d^{n+1}\ \perp\ d^{n+1}-h_{top}\ge0,
\qquad
0\le R\ \perp\ DSTMAX-d^{n+1}\ge0.
\]

The first pair is contact between retained water and top pressure; the second
is overflow against depression-storage capacity. Together with the mass
identity they define `0 <= d <= DSTMAX` and `R >= 0`. The active
sets are `supply_limited` (`d=R=0`, `h<=0`), `storage`
(`0<d<DSTMAX`, `h=d`, `R=0`), and `capacity_limited`
(`d=DSTMAX`, `h=DSTMAX`, `R>=0`). `DSTMAX=0` is the no-ponding Signorini
special case.

For each trial/Picard cycle:

1. materialize requested liquid supply and evaporation demand;
2. solve with the current active mode;
3. compute actual top flux and surface head;
4. test dry-limit, ponding, and runoff constraints;
5. switch mode inside Picard, clear Anderson history, and reassemble the
   complete water/heat/dry-air map at the same `dt`;
6. independently require nonlinear convergence, a stable active set, the
   surface-reservoir residual, and the existing soil/air mass gates;
7. on an unresolved active-set cycle, return a controlled failure so the
   outer scheduler shrinks `dt` rather than terminating with an assertion;
8. commit surface storage and runoff only after acceptance.

The accepted outer-step diagnostics are `surface_pressure_head`,
`surface_water_available_flux`, `soil_top_water_flux`,
`surface_water_storage_change`, `surface_runoff_flux`, and the final active-set
code. The internal storage-corrected reaction used by the active set is

\[
R^*=P_l-E_s-q_{in}-\frac{d^{n+1}-d^n}{\Delta t},
\]

and equals runoff within the accepted surface-mass tolerance. This
storage-corrected reaction is the
Gatti-like KKT quantity; the raw `q_available-q_in` is retained as a forcing
diagnostic but is not mislabeled as runoff because it includes releasable
surface storage. Maximum mass and both complementarity-product residuals are
also written for every accepted outer step.

This is the node-level engineering analogue of the pressure/flux
complementarity in Gatti et al. (2024), Sect. 2.2 Eq. (2), generalized from
their explicitly no-ponding slope to finite `DSTOR`. Their RT0-P0,
Nitsche, and hybridized formulations remain a future 2-D/Julia direction,
not part of the MATLAB tridiagonal solver. The implementation also follows
the 2025 CATHY result that boundary switching needs its own convergence gate;
adaptive time-step shrink remains the numerical fallback consistent with the
2025 advanced Noah-MP hydrology experiments.

References:

- Gatti et al. (2024), *Computer Methods in Applied Mechanics and
  Engineering*, arXiv:2407.07865.
- [CATHY boundary-switch convergence study (2025)](https://www.frontiersin.org/journals/water/articles/10.3389/frwa.2025.1553578/full).
- [Advanced Noah-MP hydrology experiments (2025)](https://hess.copernicus.org/articles/29/547/2025/hess-29-547-2025.html).

Required cases include light rain, ponding, excess runoff, drainage of a pond,
dry-limit evaporation, mode-switch retry, reject/no-leak, a real
`0.01 cm s^-1` pulse, and controlled active-set-cycle rejection. Noah-MP's
50 mm experiment is not a default parameter: `DSTMAX` remains a site-scale
microtopographic/depression-storage input.

W5 changes the nonlinear boundary solve, not rainfall information content.
Changwu daily totals are still disaggregated to hourly forcing upstream, and
soil substeps hold that hourly rate fixed. Event-aligned minute forcing is a
separate forcing-interface task; adaptive soil `dt` cannot reconstruct peaks
that are absent from the input data.

### W6: root uptake and source capacity limits

- positive extraction cannot exceed available water divided by `dt`;
- hydraulic redistribution cannot exceed receiving pore capacity divided by
  `dt`;
- record requested, applied, actual, and unmet transpiration;
- preserve SPAC ownership of root-water physics and independent root-N
  ownership of `SinkS`; PlantCN remains an optional demand provider/consumer;
- maintain `integral(Srt)=Tp_actual`.

### W7: advanced lower boundaries and heterogeneity

Add only with benchmark coverage:

- seepage face;
- dynamic groundwater head;
- capillary-rise lower boundary;
- lysimeter/free-outflow variants;
- layered and strongly heterogeneous columns.

### W8: HYDRUS-class benchmark matrix

Build a reproducible matrix across boundary types, soil textures, saturation
states, grid resolutions, time steps, and uptake modes. Compare profiles,
cumulative fluxes, storage, convergence history, and mass residuals.

## 7. Time-series rules

- `time_s` is elapsed seconds from simulation start.
- Flux series are piecewise constant on
  `[time_s(i), time_s(i+1))`.
- Trial-average flux is computed from the exact interval integral.
- The final row is a coverage endpoint; no extrapolation is allowed.
- Head series are evaluated at the trial endpoint.
- Head series use automatic linear interpolation at the trial endpoint.
- Times must be finite and strictly increasing.

Signed fluxes follow the single public convention: positive is into the soil
domain. Therefore:

- positive top flux is downward input;
- negative top flux is upward removal, including a controlled evaporation
  experiment;
- positive bottom flux is upward recharge;
- negative bottom flux is downward outflow.

The transitional Richards kernel retains its historical top `BCh/QMT`
upward-positive convention. The adapter applies
`legacy_top_flux = -public_top_flux` on input and the inverse conversion on
output. Consequently, an old prescribed evaporation value `+E` migrates to
CSV `top_value=-E`; copying the old sign directly would reverse the experiment.

For example, a laboratory top boundary with `top_value=-1e-6` under
`constant_flux` prescribes upward removal at `1e-6 cm s-1`. Natural evaporation
is different: it remains a nonnegative demand generated by the energy-balance
chain, and the surface transaction determines the actual signed soil flux.

These rules guarantee that retry and alternative substep tiling do not change
the prescribed cumulative water amount.

## 8. `h_sur_cm`, irrigation, and surface water

These quantities have different owners:

| Quantity | Meaning | Owner | Default role |
|---|---|---|---|
| `h_sur_cm` | observed surface pressure/matric head | observation forcing | diagnostics and calibration |
| irrigation | managed water application | management source | routed by application method |
| surface water | ponded water depth | `Y.surface_water` | accepted-only prognostic state |

Irrigation is not itself a boundary kind. Planned routes are:

- `overhead`: above canopy, then interception and snow partition;
- `surface`: direct surface liquid supply;
- `subsurface`: internal soil source;
- `ponded_head`: explicit prescribed surface pressure head.

No route may be counted twice. A W2 fixed-flux laboratory input is not yet a
complete agricultural irrigation implementation.

## 9. ClimaLand architecture mapping

| ClimaLand concept | COMPASS counterpart |
|---|---|
| concrete water BC types | `built.boundary.water.top/bottom.kind` |
| water flux BC | `fixed_flux` with explicit public sign |
| state BC | `fixed_head` with explicit pressure/matric head |
| free drainage | bottom `free_drainage` |
| time-varying input | pure natural or prescribed materializer |
| prescribed precipitation | natural liquid-supply driver |
| atmosphere-driven BC | atmospheric request plus surface policy |
| runoff strategy | future W5 surface-water policy |
| pond model | `Y.surface_water` |
| model sources | Srt and future subsurface irrigation |

Primary local references:

- `D:/OneDrive/应用/Model/ClimaLand/src/standalone/Soil/boundary_conditions.jl:25-98`;
- `D:/OneDrive/应用/Model/ClimaLand/src/standalone/Soil/boundary_conditions.jl:164-211`;
- `D:/OneDrive/应用/Model/ClimaLand/src/standalone/Soil/boundary_conditions.jl:916-939`;
- `D:/OneDrive/应用/Model/ClimaLand/src/shared_utilities/drivers.jl:44-83,175-187,1383-1391`;
- `D:/OneDrive/应用/Model/ClimaLand/src/standalone/SurfaceWater/Pond.jl:26-93`;
- `D:/OneDrive/应用/Model/ClimaLand/docs/src/tutorials/integrated/handling_soil_fluxes.jl:1-58`.

Architecture parity means correct ownership and composition, not matching type
names.

## 10. Invariants

1. Policy, driver, state, request, and result each have one owner.
2. Each unit of precipitation, snowmelt, or irrigation enters one accounting
   path only.
3. Component composition is explicit; it is not inferred from empty fields.
4. Observations never silently become control values.
5. Requested, applied, and actual values remain separate.
6. Providers are stateless; only accepted trials commit state and time.
7. Public units and signs are unique; legacy signs remain inside adapters.
8. Unsupported combinations fail explicitly.
9. Natural Changwu Throughfall, EVAP, Srt, and water-only results remain
   regression-compatible.
10. Model and checkpoint state remain plain data for future Julia migration.

## 11. Execution order

`W0 -> W1 -> W2 -> W3 -> W4 -> W5 -> W6 -> W7 -> W8 -> I6/I7`

W1 is closed. The next implementation increment is W2, beginning with closed
column, saturated Darcy, constant-head, and constant-flux benchmarks before
changing convergence or atmospheric physics.

Before regional execution is enabled, the point pathway must pass an explicit
cross-context equivalence gate: given identical accepted initial state, root
sink, and effective top/bottom requests, the laboratory and field adapters must
produce the same Richards trial result and actual boundary fluxes. This gate is
scheduled after the common W2-W4 numerical and ledger work and before the W5
natural surface transaction is declared complete.
