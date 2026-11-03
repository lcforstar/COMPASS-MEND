# COMPASS V2.0 source authority and audit protocol

## Purpose

Use this file before every Julia migration task. It defines how to establish facts from current
MATLAB V2.0 production code. It deliberately contains no V1.0 process dossier.

Resolve path placeholders using SKILL.md §1. Reuse unchanged evidence for local follow-ups.

## Authority

The current production dependency closure under `<MATLAB_ROOT>/src` is the
highest authority for current MATLAB behavior. Tests and active site inputs establish expected
behavior. V2 physics documents explain intent, but code wins for the description of what currently
executes.

The following documents are binding for the Julia target architecture, ownership, spatial hierarchy,
implementation order and completion gates:

- `<MIGRATION_PLAN>`
- `<JULIA_ROOT>/docs/architecture/compass_complete_architecture.md`
- `<JULIA_ROOT>/docs/architecture/REFERENCE_MODEL_ARCHITECTURE_REVIEW.md`

They do not rewrite MATLAB facts, but Julia implementation may not depart from them. A proposed
departure is an `ARCHITECTURE_CONFLICT` and blocks coding until the design or formal document is
explicitly revised.

Reference models provide design questions only. V1.0 and unreachable legacy code are excluded unless
the user explicitly requests historical provenance.

## Required audit sequence for initial process migration

For a local follow-up, update only affected evidence and contracts. N/A with a reason is valid.
A conflict blocks only dependent implementation, not unrelated authorized work.

1. Identify the current migration phase and exact deliverable.
2. Start from `runCOMPASS`, an independent production entry, or a test that invokes the target.
3. Trace callers and callees until every input, state mutation, candidate, result, ledger, commit,
   rollback, diagnostic and output is accounted for.
4. Read every relevant file completely, including local functions. Search results locate files; they
   do not replace reading.
5. Trace each TOML key and input column to its consumer. A written key/column is not necessarily used.
6. Trace every switch and enum through defaults, cascades and runtime selection.
7. Mark unreachable alternatives and experimental paths separately.
8. Compare source behavior with tests and physics docs. Record conflicts rather than resolving them
   silently.
9. Produce the source map and physics contract below.
10. Design the Julia refactor and verification fixture.
11. Check that design against all three binding migration/architecture documents before editing `.jl`.

## Source map template

```markdown
### Capability

- V2 production entry:
- Upstream caller:
- Direct callees read completely:
- Tests/fixtures:
- Config switches and actual consumers:
- Input files/columns and actual consumers:
- Activity status:
- Known unreachable/legacy alternatives:

### Mathematical and physical contract

- Domain/control volume/orientation:
- Prognostic unknowns and time levels:
- Parameters and units:
- Forcing and frozen retry interval:
- Governing balance:
- Constitutive closures:
- Boundary/source request:
- Applied result:
- Discretization and solve order:
- Acceptance/rejection rule:
- Conservation ledger:
- Diagnostics:

### Julia design before code

- Physical owner:
- Parameters type:
- State type:
- Aux/workspace type:
- Request/result/ledger types:
- Candidate/commit/rollback boundary:
- Disabled capability status:
- Parity and invariant tests:
- UNRESOLVED:
```

Any non-empty `UNRESOLVED` that affects mathematical meaning, units, orientation, activity or state
ownership blocks implementation.

## Fact labels

Use these labels in audits and plans:

- `V2_SOURCE_CONFIRMED`: read in current production closure.
- `V2_TEST_CONFIRMED`: exercised by a current test or fixture.
- `V2_DOC_INTENT`: stated in V2 documentation but not yet proved active.
- `JULIA_DESIGN_PROPOSAL`: target architecture, not current MATLAB fact.
- `REFERENCE_MODEL_PATTERN`: borrowed design question or pattern only.
- `LEGACY_ONLY`: V1 or unreachable code; not a migration oracle.
- `UNRESOLVED`: V2 evidence incomplete or conflicting.
- `ARCHITECTURE_CONFLICT`: proposed Julia design violates a binding migration/architecture decision.

## Conflict gate

Pause dependent implementation when any of these affects the current change and remains unresolved:

- source and test expect different behavior;
- source and physics document disagree on boundary or sign;
- a config key is written but has no confirmed consumer;
- an output is diagnosed from request rather than applied result;
- a state has more than one writer;
- an enabled path lacks a valid input or failure contract;
- regional/global support is inferred from array shape without a runnable solver path.

Resolve the conflict in the plan or V2 source first. Do not hide it with a Julia default.
