# COMPASS-MEND 1.0 project context

This package is the MATLAB COMPASS-MEND implementation with the accepted
Changwu winter-wheat input baseline. Preserve original author/Change records
and the upstream rights statements in NOTICE.md. Open describes the original
modular architecture; it does not grant a software license.

- Read src/compass_skill_routing.md and the applicable local coding skill
  before changes affecting model semantics.
- Canonical run flow: run_site -> runCOMPASS -> init.initialize ->
  init.prepareRuntime -> beginRun -> accepted runOneStep -> finalizeRun.
- src contains production code. data/sites/changwu contains scientific inputs;
  data/model_library_nc supplies required optical/reference data.
- Current example: winter wheat, 2000-06-15 through2001-06-15; source_hourly
  MEND and mend_dynamic sorption, FAST10, passive root water/N sinks, no active
  prescribed root-residue input and no plant-C/N or automatic litter input.
- Retain the native units, orientation, accepted/trial transaction boundaries,
  water/heat/C/N/O2 checks and original tolerances. No state clamping or relaxed
  checks to make tests pass.
- For water/heat/gas/freezing/top-boundary work, read the relevant sections of
  docs/soil/compass_soil_water_heat_physics.md. Changwu atmospheric heat top is
  the matched fixed G* Neumann boundary, not an imposed skin temperature.
- Run the focused checks appropriate to a change. tests/verifyStandaloneRelease.m
  provides a bounded public-entry smoke check; it is not an annual validation.
- Generated output and local credentials do not belong in version control.
- Source of the accepted input values and their active/compatibility status:
  docs/release/ACCEPTED_BASELINE.md and docs/parameters/.

Current coupling/site editor: Wenjie Liu. Preserve the original COMPASS and
third-party attribution; do not infer ownership or a new license from this file.
