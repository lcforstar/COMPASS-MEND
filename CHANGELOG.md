# Changelog

## 2026-11-03 — DOC/DON flux at unequal cell interfaces and O₂ sink accounting

- Interpolate liquid flux from the two cell centres to their shared interface by centre-to-face distance. This flux drives advection, mechanical dispersion and donor checks. Active nonuniform-grid DOM runs now bind cell widths and the face reconstruction to the checkpoint.
- Retain the nitrification and respiration components of the linearized O₂ sink during replay. This preserves a small physical sink when subtracting two rounded total rates would give zero. The convergence and budget tolerances retain their original values.
- Tested the June storm and a continuous 168-hour run, a KT83–96 restart, rejection of the previous nonuniform-grid operator, and the historical August 2001 O₂ failure. [Paired hydrology and carbon transfers](data/verification/cand30_local_face_audit_20261103_01/非均匀网格与20cm通量核查.md)

## 2026-11-03 — DOM liquid-water timing and inorganic nitrogen transport numerics

- Advance free DOC and DON on each accepted water-flux interval. Explicit solute steps use the liquid volume at their own start. The existing concentration boundaries, immobile dynamic pools and source-label mass ledgers follow the accepted transfers.
- Pass the site's upwind and artificial-dispersion options to the inorganic nitrogen transport solver. The exponential upwind weight now uses its analytic limit near a zero signed Péclet number.
- Bind the active DOM water driver and nonzero DIN transport options to checkpoint compatibility. Active upwinding also records the weight algorithm. Older trajectories that used a different operator require an explicitly documented research branch or a new run.
- Add focused DOM and DIN regression checks. Earlier calibration output describes the transport implementation used by its original run.

## 2026-11-01 — Dissolved oxygen concentration units (Stage 0)

- Remove the extra factor of 1000 from the gas-species kernel's dissolved
  concentration output, which is already in mol/L water.
- Preserve the separate /1000 conversion for phase storage, all transport
  equations and constants, site inputs and native conservation thresholds.
- Add independent Henry, oxygen-response and phase-inventory checks, including
  bitwise fixed-input dry-mole-fraction and ledger parity with the frozen old kernel.
- Re-evaluate coupled MEND runs separately: earlier annual fits describe the
  pre-repair oxygen response and are not post-repair calibration evidence.

## 2026-10-31 — Layer-capable MEND sorption interface

- Move Qmax, Kba and Kdes_ref into 18-row soil profiles while preserving their default values and legacy TOML-scalar compatibility.
- Update the site to 11 MEND soil columns / 35 microbe TOML keys and regenerate all current parameter catalogues.
- Preserve initial pools, management inputs, source tracing and native conservation thresholds.
- Record this interface regression separately from initial-package checks and the historically adopted annual run; no new annual evaluation is claimed.

## 1.0 — 2026-10-31 (initial package assembly)

- Assemble a standalone Changwu winter-wheat COMPASS-MEND MATLAB package.
- Adopt the accepted K1 mineral-fertilizer0-20cm inputs and B_ref organic-N
  partition, FAST10 and three fertilizer-origin cohorts.
- Keep all production src files byte-identical to the validated working code.
- Provide portable launchers, explicit run receipts, bounded smoke/integrity
  checks, and a complete parameter catalogue.
- Include the required NetCDF model libraries and original input provenance.

This is the package version, not a renaming of an upstream COMPASS release.
