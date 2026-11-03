# COMPASS V2.0 to Julia refactor contract

## Goal

Re-express verified V2.0 mathematics, physics, discretization and transaction semantics as typed,
testable Julia modules. This is a semantic refactor, not line-by-line transliteration and not a
replacement by another land model.

Every design must also conform to the binding migration plan and architecture documents listed in
`00-v2-source-authority.md`. This reference refines those documents; it does not replace them.

## Design order

For every capability, design in this order:

1. source map and activity status;
2. domain/control volume/orientation;
3. prognostic state and time levels;
4. governing balances and constitutive closures;
5. parameters, forcing, boundary and sources;
6. discretization and numerical solver;
7. request, applied result and conservation ledger;
8. candidate, commit, rollback and failure report;
9. accepted-state diagnostics and I/O;
10. parity and invariant tests.

Directory or type design before items 1–6 is provisional only.

## Ownership

| Concern | Julia owner | Constraint |
|---|---|---|
| atmospheric drivers | Atmosphere | read-only current-step forcing |
| canopy, plant hydraulics, roots, phenology, plant C/N | Vegetation | owns its prognostic state |
| radiation, surface water, snow and surface energy memory | Surface | owns surface candidates/state |
| water, heat/internal energy, ice, gases and soil N | Soil | owns soil state and local solvers |
| external groundwater driver | Groundwater boundary | no invented aquifer state |
| cross-domain exchange | typed Interfaces/Coupling | no long-lived state ownership |
| local numerical algorithms | Solver | may mutate trial workspace, never commit |
| atomic acceptance | Step | sole writer of accepted multi-domain state |
| closure checks | Conservation | read-only pre-commit decision support |
| outputs/restart | Diagnostics/IO | accepted state/results/ledgers only |

Single-domain modules must not import integrated orchestration to obtain hidden state.

## Core categories

- `Parameters{FT}`: immutable and concrete.
- domain-owned accepted state: only Step replaces/commits it.
- `ForcingAtTime{FT}`: immutable slice frozen for a trial/retry interval.
- `Aux{FT}` and `Workspace{FT}`: recomputable or disposable; preallocated where hot.
- `Request{FT}`: demand and boundary intent.
- `Result{FT}`: actual applied flux/state response.
- `Ledger{FT}`: conservation quantities in the same control volume and time discretization as state.
- `StepReport`: convergence, iteration, acceptance, reason and artifacts.
- runtime memory: explicit fields needed for exact resume but not part of a physical field.

Do not create a generic mutable Core/model bag. Avoid untyped `Dict{Symbol,Any}` in production
physics and avoid one universal flux struct shared by unrelated domains.

## Input-first contract

Implement input as visible stages:

```text
site/experiment/spatial source
  -> raw parse with provenance
  -> raw validation
  -> active component plan
  -> typed normalized configuration
  -> grid/time/forcing/parameter/initial-state derivation
  -> complete model inputs
```

Optional TOML means a key can be absent only when the schema supplies a documented default or the
component is inactive. If present, the key must be consumed and validated; silently ignoring written
keys is forbidden. Preserve raw values and provenance separately from normalized/derived values.

Use the same typed component input contract for laboratory, point, regional and future global runs.
Spatial adapters supply per-column values and execution partitions; they do not alter column physics.

## Output and restart contract

Define each process's in-memory result/ledger/diagnostics contract during P3–P10, but implement the
complete persistence layer only in P11, after model layers and before P12 integrated execution. This
keeps every process testable without violating the binding migration order.

- Diagnostics register stable variables with units, dimensions, orientation and source owner.
- Writers consume accepted state, applied results and committed ledgers only.
- NetCDF/other files are projections; they are not model state.
- Checkpoints store all and only information required for exact resume: accepted prognostic state,
  committed ledgers and non-recomputable runtime memory.
- Aux/workspace/registries are rebuilt unless exact evidence proves they carry non-recomputable state.
- Schema version and downgrade behavior are explicit and tested.

## Spatial contract

`Column` is the minimum physical and atomic transaction unit. Tile/patch/grid/global layers may map
parameters, forcing, area weights, routing and I/O, but cannot change the column equation, orientation,
solver order or acceptance rule.

Preserve the V2 vertical convention: depth positive; arrays bottom-to-surface. Do not reverse Julia
internals for stylistic alignment with another model. Adapt external datasets explicitly.

Regional initialization support in MATLAB is not evidence of regional solve support. Add a spatial
capability only after an executable multi-column path, failure isolation, aggregation and tests exist.

## Reference-model use

- CLM5/CTSM: ask about control volumes, state ownership, active sets, area aggregation and closure.
- ClimaLand: ask about concrete `FT` types, standalone/integrated boundaries and lifecycle APIs.
- CoLM: ask about invariant/state/forcing/flux/accumulated-flux separation and patch/worker/I/O roles.

Never copy their governing equations, coordinate convention or integrator over a verified COMPASS
choice merely because their software architecture is mature.

## Process completion gate

A process is complete only when all are present:

- complete typed input;
- candidate state and applied output;
- exact ledger and closure calculation;
- explicit failure and rollback behavior;
- accepted-state diagnostics;
- V2 source map and activity status;
- parity tests for active paths;
- invariant/disabled tests for inactive retained paths;
- type-stability and allocation checks appropriate to the hot path.

Stubs, all-zero returns, interfaces without physics, or types without execution are scaffolding and
must be labeled `not_implemented`.
