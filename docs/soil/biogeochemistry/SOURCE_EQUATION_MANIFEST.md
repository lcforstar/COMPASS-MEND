# DOM–MEMC–CoLM source-equation manifest

Date: 2026-08-21  
Implementation repository: `C:\Users\Administrator\Desktop\COMPASS-DAISY-DOC`

## Frozen MEMC reaction equations

- Upstream file: `C:\Users\Administrator\Desktop\MEMC-main\MEMC-main\R\model_component.R`
- SHA-256: `3F30A1B214FC4A4FF2DB337C7DEBB8BD3E1DE71A7113FC04AEA0BB08F06EC481`
- Local transcription: `src/+soil/+microbe/memcFluxes.m` and `memcDerivatives.m`
- Selected upstream dynamics: `F1="MM"`, `F2="MM"`, `F8="LM"`.
- `F1`–`F10` and all nine state derivatives are transcribed without algebraic correction. In particular, the upstream expressions `dEP = F9_em - F10_ep` and `dTot = Input_POC + Input_DOC - (F4 - F5)` are retained exactly.
- Coupling changes only state ownership, units, numerical integration, and call order. The reaction equations are not altered.
- MEMC software time is days and kinetic rates are `d^-1`. COMPASS converts
  accepted seconds to model days with `dt_s/86400`; see `MEMC_TIME_BASIS_AUDIT.md`.
- MEMC adsorption/desorption equations remain present, but `K_ads=K_des=0` because equilibrium sorption is owned by the DOM module. This prevents double counting.

MEMC has one microbial biomass pool (`MB`) and an explicit mortality flux (`F8`); it does **not** contain an active/dormant biomass split or an active-motility equation. Those processes have therefore not been invented locally under the frozen-equation constraint. CoLM mixing provides passive profile redistribution only.

## Frozen CoLM vertical-transport equation

- Upstream file: `C:\Users\Administrator\Desktop\CoLM\CoLM-master\main\BGC\MOD_BGC_Soil_BiogeochemLittVertTransp.F90`
- SHA-256: `70E59F3BB69B40099517FED4F4B482887AA3A0A7977F973087A510122D9DC731`
- Local transcription: `src/+soil/+litter/colmVerticalTransport.m`
- Preserved operations: harmonic interface diffusivity, Patankar power-law `A(Pe)=max[0,(1-0.1|Pe|)^5]`, coefficient assembly, source placement, ghost equations, and tridiagonal solution.
- The operator is applied to Daisy pools or MEMC solid/microbial pools according to the active owner. It is passive redistribution, not a new microbial reaction.

## COMPASS-owned coupling equations

- `src/+soil/+dom/equilibrate.m`: linear equilibrium sorption, with total DOM as the prognostic owner.
- `src/+soil/+dom/transportConservative.m`: conservative finite-volume DOC/DON transport on the accepted water-flux field.
- `src/+soil/+dom/advanceAcceptedCandidate.m`: Michaelis–Menten DON mineralization whose sole immediate inorganic-N product is NH4; direct DON-to-NO3 source is fixed to zero.
- Outer transaction interval: 3600 s. MEMC reaction interval is configurable and defaults to 300 s; accepted inner intervals sum exactly to the outer interval. These seconds are divided by 86400 before evaluating day-based MEMC derivatives.

## Automated equation locks

- `tests/soil/nitrogen/test_dom_memc_equations.m`
- `tests/soil/nitrogen/test_colm_litter_transport.m`
- `tests/init/test_resolveBiogeochemistryMode.m`
- `tests/io/test_biogeochemistry_checkpoint_compatibility.m`

These tests lock equation identities, mass closure, non-negativity, switch ownership, restart topology, and the 3600/300 s clock contract.
