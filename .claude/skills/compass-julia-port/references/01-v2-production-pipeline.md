# COMPASS V2.0 production pipeline baseline

Audit date: 2026-07-22. Re-verify against source before each migration task.

## Runtime lifecycle

The canonical current chain is:

```text
src/runCOMPASS.m
  -> init.initialize(siteDir) or init.initializeRegional(...)
  -> init.prepareRuntime(model)
  -> beginRun(model)
  -> runOneStep(model, kt)
  -> finalizeRun(model)
```

`runCOMPASS` stops on the first rejected forcing interval, preserves the last accepted state, records
the failure and attempts an accepted-state checkpoint. It does not silently skip failed forcing.

`beginRun` performs one-time side effects after runtime allocation: resume materialization, output
directory creation and configured groundwater setup/forcing adjustment. Resume rebuilds disposable
soil auxiliary data from accepted state rather than accepting serialized workspace as truth.

`finalizeRun` flushes registered diagnostics and optional reports. It is not a TODO shell.

## Initialization chain

`init.initialize` executes:

```text
io.loadSiteConfig
  -> init.resolveSimulationPlan
  -> io.loadSite(siteDir, plan, config)
  -> io.checkSite(raw, plan)
  -> init.buildState(raw) or init.buildWaterOnlyState(raw, plan)
  -> init.allocateWorkspace
  -> init.initializeSoilAux
  -> init.alignInitialThermalBoundary
```

The returned model retains `model.raw` for in-memory provenance. This is not the same as writing an
immutable raw-input archive or hashing every source file; no such production guarantee should be
claimed without new implementation and tests.

### Stage 1: raw read

`io.loadSite` is a read-only assembly stage. It parses `site.toml` and enabled tables at native
granularity. It performs no grid construction, interpolation, name mapping, van-Genuchten inversion
or unit conversion.

The simulation plan determines which inputs are active. Laboratory water-only mode can omit
vegetation/SCOPE/full forcing inputs that a field run requires. Optional values may be absent only
when the plan and schema permit it; when a value is written, the parser must validate and consume it
or reject/declare it unsupported.

### Validation

`io.checkSite` validates raw inputs before Stage 2. Validation is plan-aware. A schema requirement
does not prove downstream consumption, so config and column consumers still need call-chain tracing.

### Stage 2: build

`init.buildState` is a transform-only stage with this current order:

```text
domain
-> control (switches, numerics, time)
-> boundary and sources
-> environment-owned groups
-> forcing
-> soil parameters and bulk density
-> moisture/thermal properties and soil modes
-> soil nitrogen parameters
-> vegetation parameters/state
-> SCOPE parameters/static cache
-> initial Y
-> plant hydraulic and C/N state
-> canopy/SCOPE state
```

Its principal containers are `domain`, `control`, `boundary`, `sources`, `forcing`, `parameters`,
accepted `Y`, auxiliary `p`, and initialization-time vegetation state. `allocateWorkspace` then adds
trial/scratch and diagnostic buffers. These MATLAB containers are evidence for semantics; they are
not a Julia type design to copy verbatim.

## One outer-step transaction

Current `runOneStep` performs six major stages:

1. validate the requested interval from accepted model time;
2. stage solver, surface-water and management auxiliaries from accepted state;
3. stage the coupled atmosphere–vegetation–surface candidate graph;
4. advance coupled soil `h/U/P_g` with adaptive, boundary-aware substeps and bounded surface–soil
   heat reconciliation;
5. after a complete soil endpoint, create the soil N candidate when enabled and atomically commit all
   candidates, ledgers and time; otherwise restore the outer snapshot;
6. publish clock, diagnostics and periodic artifacts with accept-only semantics.

This sequence is part of the numerical method. Julia may improve types and package ownership but may
not reorder these operations without a separately justified physics change and new oracle.

## State and transaction observations

- `model.built.Y` is accepted prognostic state.
- trial/previous state and scratch live in workspace and are disposable on rejection.
- parameters and control are not accepted state.
- runtime contains counters and non-recomputable run memory.
- diagnostics consume accepted state/report and must not drive the solve.
- request, applied result and ledger are distinct concepts even when MATLAB currently packages them
  in adjacent structs.
- checkpoint content must be derived from exact resume requirements; auxiliary workspace is rebuilt.

Do not reduce this to a universal `Y/p` ODE convention. The production step includes nested nonlinear
solves, active boundaries, event restaging and an atomic multi-domain transaction.

## Orientation and time

Current V2 convention is depth positive with internal soil arrays ordered bottom-to-surface: node 1
is deepest and the last node is the surface. Preserve this convention inside Julia. Any reversal must
be an explicit import/export adapter with a tested contract.

Keep configured forcing time, forcing row, trial time, accepted model time and diagnostic timestamp
distinct. In the audited forcing path, configured `forcing_start` plus row index constructs the time
axis; a required `date` column is not by itself proof that date values drive the model. Re-check this
before implementing Julia input time semantics.

## Scale status

- Point/site: current production target.
- Laboratory: point + laboratory water-only path is tested.
- Regional: `init.initializeRegional` and regional static overrides are tested, but the coupled
  `h/U/P_g` production solve explicitly rejects `nCell > 1`; this is initialization capability, not
  runnable regional physics.
- Global: no current V2 production global runtime was found in this audit.

Julia architecture should make later spatial execution possible without claiming capabilities that
do not yet run.

## Known audit-sensitive input facts

- A TOML key present in a site file may still be ignored by current MATLAB. Trace consumers.
- The 2026-07-22 Changwu audit found written debug/advanced keys without production initialization
  consumers, while some transport switches were determined by `buildSwitches` cascades. Re-check the
  current code before relying on this snapshot.
- No production caller was confirmed for old input/output helper paths merely because their files
  still exist. Static existence is not activity.
