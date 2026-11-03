# Source organic C/N exchange reference

Diagnostic-only, scalar one-box hourly kernel; source MOD_MEND.F90 subMEND
lines 1367–1578 plus intracellular enzyme ledger at 1605–1650. No production
owner changes. Inputs are pre-hour C, organic N, free NH4/NO3/NO2, parameters,
and external C/N additions. Fixed source enzyme CN=3, microbial CN bounds
2–14 for the compared case; reject nonzero source plant uptake.

Order: carbon core without final N enzymes; donor-ratio organic N transfers;
first DIN kinetics and mineralization; dormant CN correction; active CN
correction and second uptake/overflow; intracellular enzyme C/N accounting.
Outputs include paired C/N states, pre-inorganic-chemistry NH4/NO3 and actual
exchange ledgers. Enzyme factor derives from the supplied pre-hour organic
and mineral N owners, not from archived reaction fluxes. Verify sum C+CO2
and organic N+NH4+NO3 changes equal external additions.

Excluded: inorganic nitrification/denitrification, N gas transport, NH4
Langmuir redistribution, leaching, vegetation and soil transport. Thus compare
final organic states and exchange fluxes directly with source, but compare
returned DIN at the pre-inorganic-chemistry boundary. Do not compare it with
source final DIN without accounting for the excluded processes. No autonomous
annual DIN-equivalence claim. Tests: 8760 same-state hourly checks for each
iHR; relative material p95 <1e-4, absolute state error <1e-10 mg/g; manufactured
excess-N release, second-uptake and zero/partial-DIN overflow cases. Zero DIN
uses a finite zero-uptake extension to source's undefined 0/0 split.
