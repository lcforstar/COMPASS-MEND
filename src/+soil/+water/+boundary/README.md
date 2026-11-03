# Water boundary contract

This package owns the canonical water-boundary policy contract. It separates:

- immutable boundary definitions (`kind`, source, units and sign convention);
- external forcing or management sources;
- per-trial boundary requests;
- actual solver boundary fluxes;
- prognostic surface-water state.

`normalizeSpec` is the configuration adapter. Current site files expose one
flat typed schema inside `[boundary]` using `water_top_*` and
`water_bottom_*`. It returns one canonical plain-data spec plus temporary
aliases required internally by the current Richards kernel. Historical
integer fields and the earlier nested draft fail with explicit migration
messages. It performs no forcing evaluation and no physical state update.

`materializeRequest` is the pure trial provider. It evaluates constants, natural
surface requests, and explicit boundary series for one exact trial interval.
Piecewise-constant flux series are integrated before conversion to an interval
average, so retry and substep tiling preserve the prescribed water amount.

W1b routes each materialized request through
`solver.advanceWaterOnlyOuterStep` into the existing
`solver.attemptSolveWaterOnly` / `soil.water.solveWaterTrial` chain. The same
scheduler is used by natural `runOneStep` and
`solver.runWaterOnlyExperimentStep`; the standalone driver bypasses canopy,
energy-balance, and plant production without copying retry or Richards logic.

`buildResult` packages the immutable request together with the applied static
kind and actual solved top/bottom flux. Requests and results are audit outputs,
not prognostic state. Conservation accounting based on these actual fluxes is a
separate W4 increment.

The public flux convention is positive into the soil domain. Legacy `QMT/QMB`
and positional signs remain confined to compatibility adapters.

Architecture reference: ClimaLand.jl
`src/standalone/Soil/boundary_conditions.jl:25-98` and
`src/shared_utilities/drivers.jl:44-83`.
