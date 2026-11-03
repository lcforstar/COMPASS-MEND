# Three-case water benchmark figure contract

Core conclusion: A common mass-conservation and input-parity gate separates a
validated homogeneous Richards case from an unresolved layered case and a
real experiment whose measured conductivity law is not yet available in
COMPASS.

- Figure archetype: quantitative grid.
- Backend: Python and matplotlib only.
- Target output: double-column, 183 mm wide, editable SVG/PDF and 600 dpi TIFF.
- Panel a: Bonan homogeneous-clay final profile and model references.
- Panel b: SV62 60 min layered profile with the current diagnostic COMPASS run.
- Panel c: TV86 open-column 0.40 h digitized observations and initial condition.
- Panel d: absolute accepted-step water-balance residual against the 0.05 mm gate.
- Hero evidence: panel d, because profile agreement is not admissible without
  conservation and exact physical-input parity.
- Reviewer risk: SV62 completes numerically but fails the conservation gate;
  TV86 has digitized rather than tabulated observations and no exact COMPASS
  K(theta) run. Both limitations must remain visible in the figure and caption.
