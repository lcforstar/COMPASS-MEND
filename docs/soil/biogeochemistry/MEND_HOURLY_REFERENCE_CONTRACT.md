# Isolated hourly source-carbon reference contract

Scope: diagnostic-only Python implementation under tools/diagnostics; no
production owner or process switch changes. Source: MOD_MEND.F90 subMEND,
fV_MM and fAds, commit 92323c7, iKinetics=0, iHR=0 or 1. Full C core includes
POM1/POM2/MOM, free DOM, dynamic QOM, active/dormant biomass, three enzymes,
and respiratory carbon. The final implementation also includes aggregate
intracellular nitrogen-enzyme carbon, whose production/turnover changes MBA
after the main C section. Supply its allocation factor from pre-hour N
stocks; no measured reaction flux is supplied. Nitrogen evolution and C:N
overflow remain outside this kernel and must be checked inactive for parity.

Inputs: same pre-hour pools [mg C/g dry soil], hourly source parameter rows,
and POM1/POM2/DOM external additions [mg C/g per hour]. No reaction flux or
post-hour state may enter the kernel. Read archive reaction outputs only in
the independent comparator. One call is exactly one source hour, not an ODE
stage. Apply independent MM caps, post-uptake sorption, source biomass
redistribution, and final simultaneous pool balances. Include respiratory
carbon in total mass closure. Fail on materially negative states or bad mass.

Verification: independent known-stock conservation/zero-reaction fixture;
8760 reset-per-hour replay against full archived C states and fluxes; separately
report first initialization transition; continuous replay from source accepted
hour 1 through hour 8760 with no state resets, using the same archived hourly
parameters/additions. Printed-state comparison gate after initialization:
95th-percentile relative error <1e-4 for stocks >1e-9 mg/g on reset replay;
continuous replay target <1e-3. Record any exceeded gate as failure rather
than relaxing it. Continuous replay is conditional on supplied parameter
history, not an independent recreation of environmental or nitrogen kinetics.
