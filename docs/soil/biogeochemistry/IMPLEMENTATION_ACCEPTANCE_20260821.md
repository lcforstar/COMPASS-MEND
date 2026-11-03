# DOM–MEMC–CoLM implementation acceptance

Date: 2026-08-21  
Implemented repository: `C:\Users\Administrator\Desktop\COMPASS-DAISY-DOC`

## Result

The optional DOM, explicit MEMC microbial-carbon, PlantCN litter, CoLM passive vertical redistribution, 3600 s outer transaction, 300 s microbial reaction clock, and schema-2 restart path are implemented. The original Daisy path remains the strict fallback when DOM, explicit microbes, and litter are all disabled.

The mistakenly targeted `COMPASS-COMPASS_MATLAB` tree was restored to Git `HEAD` content. `git diff --quiet -- COMPASS-COMPASS_MATLAB` returns zero. NTFS/index timestamp metadata still causes `git status` to show several tracked files as modified, but their content diff is empty; no new coupling files remain in that tree.

## Source-equation contract

- MEMC SHA-256: `3F30A1B214FC4A4FF2DB337C7DEBB8BD3E1DE71A7113FC04AEA0BB08F06EC481`.
- CoLM SHA-256: `70E59F3BB69B40099517FED4F4B482887AA3A0A7977F973087A510122D9DC731`.
- MEMC `F1`–`F10` and derivatives are frozen, including upstream `dEP` and `dTot` expressions.
- CoLM Patankar/harmonic/tridiagonal equation structure is frozen.
- Coupling adaptations are restricted to state ownership, grid order, units, numerical integration, and transaction scheduling.
- DON mineralization has an NH4 source and an identically zero direct NO3 source.

MEMC has one `MB` pool and mortality, but no active/dormant split, microbial-N state, or active motility equation. These processes were not invented under the frozen-equation requirement. The CoLM operator supplies passive profile redistribution only.

## Verification evidence

| Test | Result |
|---|---:|
| Biogeochemistry mode and ownership | 4/4 pass |
| DOM equilibrium/transport, DON routing, MEMC equation and clock locks | 5/5 pass |
| CoLM vertical transport | 2/2 pass |
| Checkpoint topology compatibility | 3/3 pass |
| Site input validation | 34/34 pass |
| Checkpoint round trip | 24/24 pass |
| Resume continuity | 19/19 pass |
| Accepted nitrogen transaction | 3/3 pass |
| Five-domain contract | 5/5 pass |
| Full Changwu configuration, all optional modules on, one 3600 s step | accepted |

Full-coupling first-hour values:

- `sum(DOC) = 3.00854629964e-6 g C cm-3` across elements;
- `sum(DON) = 0` from a blank DON cold start with no first-hour DON input;
- `sum(MB): 11.6856706499 -> 11.6808009094 mg C g-1` across elements;
- minimum accepted NH4 concentration: `1.30730673574 ug cm-3 water`;
- minimum accepted NO3 concentration: `11.2566944573 ug cm-3 water`;
- all inspected DOM, MB, NH4, and NO3 values were finite and nonnegative.

## Initialization and production use

Blank DOC/DON profiles are accepted only through `runCOMPASSBiogeochemistrySpinup`. Production runs must resume a compatible schema-2 checkpoint. `soil.csv` retains `organic_turnover_factor`; fixed DOM, MEMC, litter, and coupling parameters live in `site.toml` between `[organic]` and `[boundary]`.
