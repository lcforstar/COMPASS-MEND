# MEND production-core migration contract

Date: 2026-09-12

## Decision

The production organic-carbon reactor is MEND, not MEMC. COMPASS retains the
soil grid, accepted-state transaction, dissolved-organic-matter transport,
plant-first DIN reservation, inorganic-N solver, litter/AOM delivery, and
diagnostics. No production reaction path may call `memcFluxes`,
`memcDerivatives`, or the former detached dormancy correction.

The source baseline is `wanggangsheng/MEND`, commit
`92323c7326e9a30ae44e794d114612b90a9e3b1c`, principally
`src/MOD_MEND.F90` and `src/MOD_MEND_TYPE.F90`. Site parameters are explicit;
the source example parameter vector is not silently treated as a Changwu
calibration.

The Changwu microbial DIN settings (`VNup_MB`, `KsNH4_MB`, `KsNO3_MB`) remain
`0.001/0.001/0.001`, matching the upstream namelist's *initial* entries
28–30. The final assignment in that file is `0.1/0.00018/0.00041`; it is
neither the current Changwu setting nor a Changwu calibration. Published and
Fortran parameter labels give the two Ks in mg N cm-3, while the upstream
active formula divides mineral inventories annotated mg N g-1 by those Ks
without an explicit conversion. COMPASS consistently evaluates the adapter
in mg N g-1 dry soil. Numerical formula parity at equal coordinates therefore
does not establish physical-unit or site-parameter parity.

The one-layer check `tools/diagnostics/verifyMendUpstreamDinFormula.m` compares
the active [upstream immobilization expression](https://github.com/wanggangsheng/MEND/blob/92323c7326e9a30ae44e794d114612b90a9e3b1c/src/MOD_MEND.F90)
with the COMPASS microbial-only request. At 20 degC, MBA-C/N = 0.8/0.1 and
NH4/NO3 = 0.002/0.005 mg N g-1, both the upstream initial and final vectors
match the independent expression to 1e-14 mg N g-1 h-1. The resulting NH4/NO3
requests are 5.0822592e-8/3.1764120e-7 (initial) and
5.25876265e-6/3.11963591e-5 (final). This is **pointwise formula parity**;
the upstream executable bundled in the repository is for macOS and this
Windows environment has no Fortran compiler, so a 1-10-day Fortran-versus-
MATLAB trajectory benchmark has not been run. It must also control the
different DOM sorption and plant/DIN competition policies.

## State and ownership

| Production meaning | Authoritative COMPASS state | Compatibility aggregate |
|---|---|---|
| oxidative/lignin POM | `POM1` | `POC=POM1+POM2` |
| hydrolytic/cellulose POM | `POM2` | `POC=POM1+POM2` |
| mineral-associated OM | `MOC` | direct MEND MOM |
| active microbial biomass | `MB` plus physical transport phases | direct MEND MBA |
| dormant microbial biomass | `BD` plus physical transport phases | direct MEND MBD |
| two POM enzymes | `ENZP1`, `ENZP2` | `EP=ENZP1+ENZP2` |
| MOM enzyme | `EM` | direct MEND ENZM |
| total dissolved C/N | `soil_dom.doc_total_g_cm3`, `don_total_g_cm3` | no second accepted QOC/QON owner |
| sorbed dissolved phase | equilibrium diagnostic from total DOM | `QOC=0`, `QON=0` in accepted solid ledgers |
| respiration | cumulative `IC` and accepted CO2 flux | no negative-pool repair |

All MEND solid and microbial C pools use top-down mg C g-1 dry soil. Paired N
pools use bottom-up g N cm-3 in accepted COMPASS state and are converted once
at the reaction boundary. DOM remains bottom-up g element cm-3 bulk soil.

## Carbon equations

The MEND MM configuration (`iKinetics=0`) is used:

- each POM pool is decomposed by its own ENZP pool; fraction `frPOM2DOM`
  enters DOM and the complement enters MOM;
- MOM is decomposed by ENZM and enters DOM;
- aqueous DOM is taken up by MBA with `(Vg+Vm)/Yg`, while growth and active
  maintenance respiration retain the source split;
- dormant maintenance is `VmD*MBD`;
- MBA mortality is split to DOM and POM, with non-DOM residue partitioned
  0.1/0.9 between POM1/POM2;
- ENZP production follows the current POM fractions; ENZM is independent;
- enzyme turnover enters DOM;
- MBA-to-MBD and MBD-to-MBA use the source substrate and soil-water-potential
  response functions;
- dynamic MEND DOM sorption is disabled because COMPASS total-DOM equilibrium
  partition plus transport is the sole sorption owner. This is an explicit
  coupling adaptation, not a MEND source claim.

During each carbon RK stage, the integrated DOC coordinate is total DOC.
The aqueous concentration supplied to MEND uptake and dormancy is recomputed
from that stage's total using `theta/(theta+rho_b Kd_DOC)`. Thus sorbed DOC
changes with each reaction stage without introducing a second accepted pool.
The paired nitrogen RK stages likewise repartition total DON using
`theta/(theta+rho_b Kd_DON)` before aqueous DON uptake. Water content, bulk
density, and Kd remain fixed over the reaction candidate and are updated by
their respective COMPASS owners outside this reactor.

Source rates are hourly. COMPASS evaluates them as derivatives per hour and
integrates them with adaptive RK stages. The source's one-hour inventory cap is
replaced by candidate rejection/substep reduction, so 300/150/75 s convergence
is meaningful and no final-state clipping hides donor overdraw.
This is a material source-fidelity boundary for DOM uptake, not just a small
numerical tolerance: the source `fV_MM` cap was active in 8421 of 8760 hours
in the 2000-06-15 to 2001-06-14 aligned mass-coordinate case. At source hour
721, the same-state/same-parameter MATLAB raw DOM-to-MBA rate is 7.00 times
the capped Fortran transfer. The [annual audit](../../model_comparison/results/compass_mend_20000615_20010614/physical_carbon_alignment/README.md)
contains the active-microbe C ledger. A source-exact comparison would need
the full hourly discrete reaction and sorption sequence; adding a cap to each
300-second RK stage is not equivalent.

## Environmental response

The source Arrhenius factors are retained at 20 degC reference: LIG 53,
CEL 36.3, ENZ/DOM 47, Km 30, and maintenance 20 kJ mol-1. The source Gaussian
pH responses and the `fSWP`, `fSWP_OPT`, `fSWP_A2D`, and `fSWP_D2A` functions
are retained. COMPASS matric head is converted by
`SWP_MPa = h_cm * 9.80665e-5`. Changwu must explicitly select the MEND biome
and SOM response categories.

## Nitrogen and transaction order

For each accepted coupled step, using the matching carbon RK substep ledger:

1. reserve the COMPASS root NH4/NO3 request in a shadow inventory;
2. repartition total DON using COMPASS's DON Kd at every matching carbon RK
   stage; aqueous DOC uptake consumes only aqueous DON, while enzyme synthesis
   and turnover transfer N at the configured fixed enzyme C:N. At each stage,
   the active-microbe C:N sets
   intrinsic NUE and releases `(1-YgN)` of DON uptake to mineral N;
3. evaluate the first microbial DIN request without MEND vegetation terms,
   using free NH4 and liquid NO3 after root priority as the source-formula
   substrates; cap its actual uptake by the separate post-root total NH4
   (liquid plus equilibrium-sorbed) and liquid NO3 inventories;
4. return uptake-associated mineralization to NH4 only after first uptake;
5. assemble the full C/N candidate; enzyme N is removed from MBA-N at the same
   stage and an overdraw rejects the candidate;
6. calculate MBA and MBD C:N corrections separately; only MBA requests the
   remaining DIN, while each pool can release excess N independently;
7. route MEND overflow from the corrected MBA or MBD candidate to CO2 and
   update that same microbial state;
8. constrain all inorganic-N reactions by the remaining species-specific
   inventories and scale products with actual reactant consumption;
9. commit DOM, MEND C/N, DIN, plant uptake, cumulative fluxes, and management
   events atomically.

Plant priority is a COMPASS-MEND coupling assumption. MEND vegetation uptake
and its `KsNH4_VG/KsNO3_VG` denominator terms are excluded.

The organic-N candidate reports gross `don_uptake_mg_g`, retained
`don_assimilated_mg_g`, intrinsic `nue_mineralization_mg_g`, separate MBA/MBD
C:N releases, and their summed `din_release_mg_g`. The latter is delivered
once to the NH4 ledger; the same mineralized N is not kept in MBA-N.

## Acceptance gates

- no production references to MEMC derivative/flux functions;
- exact unit and top-down/bottom-up conversion tests;
- independent scalar/source-vector tests for every MEND carbon flux and
  environmental response;
- zero-DOM, zero-DIN, zero-MBA/MBD, zero-POM and zero-enzyme boundaries;
- same-stage enzyme N turnover test with initially zero enzyme;
- root reservation, first/second immobilization, overflow, and inorganic-N
  actual/request/unmet diagnostics from one candidate;
- C and N nonnegativity without final clipping;
- per-layer and column C/N closure plus rejected-step rollback;
- 300/150/75 s short coupled convergence for no fertilizer, fertilizer, low
  DIN, high root uptake, and high microbial immobilization;
- short full-physics run before any annual run.
