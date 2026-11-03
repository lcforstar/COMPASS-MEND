# Standalone release smoke check

For transport regressions, add `src` and `tests` to the MATLAB path. Run
`verifyDomWaterTiming()` and `verifyDomDriverCompatibility()` for conservative
solute redistribution, accepted water drivers and boundary cases. Run
`verifyDinTransportNumerics(pwd)` to exercise the DIN numerical options
through independent copies of the Changwu site and check the small-Péclet limit.

`verifyDomWaterTiming()` also exercises direction reversal and a dry donor at
an unequal-width interface. `verifyOxygenSinkComponents()` checks the small
nitrification sink alongside much larger respiratory oxygen consumption. The
real storm, old checkpoint rejection, and continuous vs. restarted column
comparison are recorded in [the face audit](../data/verification/cand30_local_face_audit_20261103_01/非均匀网格与20cm通量核查.md).

From a fresh MATLAB session at the project root:

```matlab
addpath('tests');
receipt = verifyStandaloneRelease(pwd, [], 4);
```

The test runs the public `run_site` entry for four accepted hourly steps with
the normal full diagnostic registry. It checks source resolution within this
copy, native C/N/DIN/O2 conservation,29 MEND profile/column flux outputs and
their NetCDF records, and all three fertilizer-N cohort diagnostics.

Output is isolated under `data/output/release_smoke_<timestamp>/`. The test
does not modify scientific inputs and does not require the original working
directory, verification archive, server helper or a previous checkpoint.

This bounded test is not an annual or a complete historical regression-suite
validation. The accepted annual example and reference metrics are documented
in `docs/release/ACCEPTED_BASELINE.md`.

For the Stage 0 dissolved-oxygen unit repair, provide the actual frozen
pre-repair kernel and a new output directory:

```matlab
addpath('tests');
report = verifyOxygenDissolvedUnitsStage0(outputRoot, baselineKernelPath);
```

This checks pressure-form Henry concentrations, oxygen-response limits,
phase inventories and dissolved outflow. With identical fixed inputs, the
old and repaired kernels must produce bitwise identical dry mole fractions
and complete ledgers; only dissolved mol/L output changes by a factor of
1000. It does not advance the model clock and does not establish coupled
annual-trajectory parity. Both inputs and the receipt are kept in the
isolated Stage 0 verification directory for this repair.

For the layer-capable Qmax/Kba/Kdes interface and source sorption algebra:

```matlab
addpath('tests/diagnostics');
report = verifyMendDynamicSorption(fullfile(pwd, 'data', 'output', 'sorption_check'));
```

This check initializes the Changwu site and uses detached reaction probes;
it does not advance its accepted model clock or change its input files. It
checks scalar/constant-profile parity, heterogeneous profiles, single-layer
responses, paired C/N transfers and the original temperature laws. It is
distinct from an independent full Fortran-trajectory comparison.
