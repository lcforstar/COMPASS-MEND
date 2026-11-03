# COMPASS.jl coding standard

Read this before creating or editing model `.jl` files. The source and physics rules in
`compass-coding-standard` remain mandatory.

## File organization

Organize a physics file so a mathematical reader can follow it top to bottom:

1. public types and invariants;
2. governing contract/docstrings;
3. public entry functions;
4. assembly and solve functions in physical order;
5. small private helpers near the logic they support;
6. no unrelated utility dump at the bottom.

Julia does not require “variables first, main loop second, helpers last” as a universal layout.
Types, constructors and methods should be grouped by the physical abstraction they implement.

## Types

- Parameterize numerical structs by concrete `FT<:AbstractFloat` where useful.
- Prefer immutable structs for parameters, domains, forcing and requests.
- Use mutable storage only for accepted state under Step control, runtime memory, or preallocated hot
  workspace with a clear owner.
- Keep field types concrete. Do not place abstract arrays/functions in hot structs without a justified
  parametric type.
- Encode boundary/process choices with types or validated enums, not magic integers spread through
  kernels.
- Preserve array location and dimensions in names, types or docstrings.

## Functions

- Use `!` only when arguments are mutated.
- Make ownership visible in signatures; do not fetch hidden global model state.
- Separate materialization (`request -> coefficients/applied boundary`) from solve and commit.
- Pure constitutive functions may return values; hot assembly/solve functions should reuse workspace.
- Provide `(dY, Y, aux, t)` only where the V2 process genuinely behaves as a tendency. Do not force
  nested Picard, active-boundary or transaction logic into an ODE facade.
- Return structured results and failure reasons rather than boolean success plus hidden mutations.

## Numerical fidelity

- Preserve V2 units, sign convention, orientation, coefficient timing, solve order and acceptance rule.
- Centralize conversions at adapters/boundaries.
- Guard mathematical domains explicitly; do not add clamps unless V2 behavior or physical analysis
  justifies them and conservation accounts for them.
- Snapshot only the state/ledger required for rollback; a rejected step must not alter accepted values.
- Keep forcing and external boundary snapshots fixed across retries when V2 does.

## Documentation

Public model functions must document:

- argument/return types, shapes, units, orientation and time levels;
- physical equation or algorithm;
- boundary/source sign convention;
- V2 source path and active condition;
- failure/acceptance behavior;
- verified literature or `[SOURCE UNVERIFIED YYYY-MM-DD]`.

Follow compass-coding-standard §10: new model functions use `# Editor: Wenjie Liu; Date: YYYY-MM-DD.` Preserve historical attribution. Append concise `# Change:` notes
for later substantive changes.

## Tests

Mirror source ownership in `test/`. For each implemented capability include the applicable set:

- constructor/schema/default/fail-fast tests;
- MATLAB raw and derived-input parity;
- pure numerical parity over normal and edge domains;
- conservation and sign/orientation tests;
- accepted/rejected transaction and rollback tests;
- disabled-path tests;
- checkpoint exact-resume parity;
- `@inferred`/`@code_warntype` checks and hot-loop allocation budgets.

Use tolerances established before looking at the Julia result. Test failure reasons and activity state,
not just final floating-point values.

## Review checklist

- [ ] V2 source map is complete and current.
- [ ] Current migration phase permits this file.
- [ ] Physical owner and control volume are unambiguous.
- [ ] State/aux/forcing/request/result/ledger/diagnostics are not mixed.
- [ ] No `model.built` clone, generic flux bucket or hidden global state exists.
- [ ] Orientation remains bottom-to-surface internally.
- [ ] Request differs from actual applied result where physics requires it.
- [ ] Only Step commits accepted state.
- [ ] Active and disabled capabilities are labeled honestly.
- [ ] Tests prove parity, conservation and rollback as applicable.
