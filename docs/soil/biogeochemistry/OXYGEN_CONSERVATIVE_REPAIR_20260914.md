# Oxygen and C/N interface repair contract

Production chain: runOneStep → soil.nitrogen.advanceAcceptedCandidate →
MEND candidate / solveNitrogen → solveO2Transport → soil.commitCandidate.
Enabled by the existing nitrogen and explicit-microbial switches.

O2 owns gas-phase node concentration in mol/cm³ gas and an accepted phase
snapshot (bottom-up endpoint liquid/gas fractions and node temperature).
The snapshot describes its old storage capacity, not a second oxygen stock.
New capacity uses the actual liquid fraction, excluding ice. The interval
balance is M_new C_new − M_old C_old = dt(boundary − actual consumption).
Boundary reactions are recovered from the original rows before elimination.

Use lumped storage, existing molecular diffusion, and donor-upwind gas
advection to obtain a positive transport operator. A nonnegative active-set
solve limits prescribed oxygen removal; it does not clip a negative state.
Return requested/actual/unmet layer removal and the full column balance.
Organic respiration is reserved before additional nitrification. If the
organic proposal exceeds oxygen supply, recompute the detached MEND C/N
proposal from the identical accepted owners with a common per-layer reaction
multiplier. Scale the complete reaction flux ledger, preserving its internal
C/N stoichiometry, rather than editing respiration after the state update.
Surface AOM respiration belongs to the atmospheric surface, not soil-pore
oxygen; incorporated AOM has no direct respiration in the existing kernel.
Transport, external additions, and roots are not multiplied. Nitrification
is limited by the remaining oxygen supply and its nitrogen products follow
the actual rate. An unresolved proposal fails before atomic commit.

This is an explicit oxygen-limited COMPASS coupling policy, not a claim of
upstream hourly MEND equivalence. No anaerobic carbon pathway is introduced.

Verification: heterogeneous-density zero-reaction invariance; detection of
AOM component/bulk N and urea changes; variable-capacity and frozen-liquid
O2 storage; no-O2/no-input reaction rejection through returned actual demand;
Dirichlet/closed-bottom mass balance; large advective Peclet positivity;
full production smoke and accepted diagnostic checks. Failed candidates may
not update accepted owners, management events, or cumulative diagnostics.
