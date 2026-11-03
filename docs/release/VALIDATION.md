# Independent release verification

Editor: Wenjie Liu; Date: 2026-10-31.

Status: PASS on Windows, MATLAB R2025a25.1.0.2943329, with JVM enabled.

The current Qmax/Kba/Kdes layer-interface revision passed a fresh four-hour
standalone run and a separate old-code/old-scalar, new-code/old-scalar,
new-code/uniform-profile comparison. Cold and accepted states,29 MEND
fluxes, C/N/DIN/O2 budgets, registry values and NetCDF numerical values were
identical by `isequaln` and serialized-value SHA256. File-container bytes
were not used as a substitute for physical-value comparison.

The native18-to-58-layer mapping, sole ownership, physical-domain rejection,
zero Kba/Kdes limits and RK4 initialization without inactive sorption inputs
also passed. Detached heterogeneous and single-layer probes passed the
original source algebra and C/N checks. No tolerance was relaxed.

The launcher, GUI and profile-rendering checks listed below are the initial
package checks; these unchanged tools were not rerun for this migration.

- A fresh MATLAB process restored its default path and resolved8 production
  entry points exclusively inside this release.
- The public run_site entry accepted4 hourly steps under the normal full
  registry. Original C/N/DIN/O2 conservation checks passed.
- All29 MEND profile and column-rate diagnostics were finite in memory and
  read back from native NetCDF. Three fertilizer-N cohort diagnostics retained
  their original tolerance and zero numerical correction.
- Root run.bat separately started and completed1 requested step, returned0
  in wait mode and wrote this run's actual MATLAB PID/completion sidecar.
- A deliberately invalid MATLAB path returned a nonzero exit without hanging
  at pause; GUI logic reports failure even if all KT rows preceded a flush
  error, and distinguishes short completion from a full annual run.
- Bash launch scripts passed syntax/LF checks. The optional GUI's Tk import
  and receipt logic passed without opening a window.
- The native profile viewer rendered PNG/PDF and its five panels were inspected.
- All817 production src files are synchronized with the current working
  code. Five interface files changed from the initial package to support
  soil-owned sorption profiles and elementwise rate construction. The
  physical default values and initial pools retain the accepted values;
  the release default output path remains data/output/changwu.

This is a bounded independent-copy smoke test, not a new annual simulation,
Linux runtime validation or execution of the complete historical test suite.
See validation.json and ACCEPTED_BASELINE.md for the distinct evidence scopes.
Generated smoke output is archived outside the distributable package.
Legacy scalar-hourly checkpoints are not silently promoted to the new
profile configuration: restart signatures still compare these shapes
strictly. The current profile configuration was validated from cold start.

CSV display cleanup: soil.csv and organic_initial_pools.csv now have UTF-8 BOM and one leading comment each. Machine headers and data bytes are unchanged. All four before/after tables (working and standalone copies) compare exactly through the native MATLAB reader. No model timestep or source edit was required. Original header history is retained in docs/parameters/INPUT_HEADER_HISTORY.md.

The single leading CSV comment now aligns a Chinese meaning and unit with every column (40 soil and70 organic-pool fields). All four native-reader tables remain exact; no parameter value or source file changed.

Input layout cleanup: named soil values and all399 TOML values are unchanged. The soil front20 columns retain their order, empirical mineralization follows the inorganic-N group, and inactive linear DOM-Kd/microbial-transport columns are last. Native reading, checkSite and buildState parameters/initial states compare exactly in both copies; no source edit or model timestep was required.

Applicable original-author TOML explanations were restored with source attribution. Crop-specific source notes were excluded, and outdated unit/boundary statements were adapted to current consumers. All non-comment code lines and all399 values remain exact, including native MATLAB serialized parsing; soil/input pools and production source did not change.
