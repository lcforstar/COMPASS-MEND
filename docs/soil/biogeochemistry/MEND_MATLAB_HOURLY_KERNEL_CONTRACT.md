# MATLAB source-hourly C/N kernel

## Scope and source evidence (2026-09-14)

| Contract | Evidence / decision |
|---|---|
| Entry | `soil.microbe.mendSourceHour`; detached value-input/value-output kernel. Production `runOneStep -> soil.nitrogen.advanceAcceptedCandidate` currently still uses the RK4 path. |
| Source | MEND `92323c7326e9a30ae44e794d114612b90a9e3b1c`, `MOD_MEND.F90:1147-1578,1605-1650`; precise executable only changes output precision relative to the documented aligned build. |
| Active subset | MM kinetics 0, iHR 0 or 1, fDOM=1, plant uptake disabled, enzyme C:N=3, microbial C:N bounds 2–14. Other source branches are outside this kernel's contract. |
| Inputs | Struct fields named like the archived source variables/parameters; scalar or aligned column states in mg element/g dry soil. Each column entry is an independent box, with no grid permutation. Additions are amounts per hour, not rates per second. |
| State | Separate QOM and free DOMT C/N; separate aggregate intracellular ENZNmt C/N; cumulative CO2. No accepted COMPASS owner is modified. |
| Outputs | New C/N, pre-inorganic-chemistry NH4/NO3, actual C/N transfers, mineral uptake/release, overflow, and C/N residuals. NO2 is an unchanged boundary input to enzyme allocation. |
| Math | Source simultaneous capped hourly C update; source donor-C:N transfers; first DIN uptake; dormant correction; active correction and second uptake; intracellular enzyme production/turnover. |
| Transaction | Pure function, no I/O or global/cache writes. Invalid proposals throw before returning. Commit, root priority, oxygen and transport are outside this kernel. |
| Verification | Independent Fortran final organic pools and actual exchange for both full-year iHR cases; continuous organic-state evolution with externally supplied source DIN; manufactured DIN/overflow/turnover cases; invalid input and no-mutation checks. |

This is a fixed **3600 s discrete map**, not an ODE derivative. Running it
twelve times with the same hourly caps does not implement a 300 s RK4 step.
The API rejects any other duration. The iHR=0 respiratory caps remain separate
from the DOM uptake cap exactly as in source.

Conservation per box is `sum(delta C pools including CO2) = organic C input`
and `sum(delta organic N) + delta NH4 + delta NO3 = organic N input + deposition`.
An invalid negative result is rejected, not repaired by clipping. The source
negative-active-C correction is retained as a signed MBA-to-MBD transfer.

Explicit finite extensions: zero DIN produces zero first uptake; an empty
C/N donor carries zero N with zero C flux; zero growth and maintenance under
iHR=1 produces zero respiration. N-only organic donor states and undefined
POM enzyme allocation are rejected. Intracellular enzyme C/N must initially
be consistent with the fixed source C:N, otherwise the source overwrite would
manufacture nitrogen.

## Production integration gate

The source kernel alone does **not** validate the coupled model. Before making
it a selectable production mode, the following interfaces must be implemented
and tested together:

- `soil_dom` owns QOM/DON sorbed inventories and transports only the free
  fraction; diagnostics must not re-equilibrate them with Kd.
- Intracellular enzyme C/N must appear in initialization, accepted storage,
  microbial transport, stock diagnostics and spin-up convergence. Current
  BIOPARTICLE transport only carries MBA/MBD C/N.
- Root-reserved free DIN drives kinetics; the total mineral inventory limits
  actual uptake. Oxygen retries must rebuild paired C/N proposals from the same
  accepted state, including overflow and enzyme costs.
- Checkpoint restoration and rejected outer steps must retain every new store.

No site switch is added until those consumers exist; no annual production
alignment claim follows from the standalone tests below.

## Implemented DOM interface and verification boundary

`soil.dom.phaseState` now reads explicit schema-2 dynamic total/sorbed DOC and
DON inventories. `soil.dom.advanceAcceptedCandidate` transports the difference
with Kd=0 and retains the sorbed state unchanged, including its contribution to
total storage diagnostics. This is a tested production DOM function, but the
site initializer does not yet construct dynamic states. The current microbial
RK4 entry explicitly rejects dynamic phase states, preventing a partially
migrated caller from silently re-equilibrating them.

The dynamic transport test covers changed liquid water, independence from Kd,
closed-column C/N conservation, invalid-owner failures and MAT serialization
followed by identical continuation. **The MAT roundtrip is not a full COMPASS
checkpoint/outer-transaction acceptance test.** Existing equilibrium-mode
production O2 replay and its checkpoint test also pass; they do not establish
oxygen coupling of the new source-hour kernel.

Remaining production work: dynamic initialization/old-checkpoint migration;
paired source-hour owner adapter and parameter mapping (including intracellular
enzyme C/N); explicit root-first/free-vs-total DIN coupling and common C/N/O2
retry semantics; intracellular enzyme transport; accepted-state diagnostics
and spin-up coverage; then full production short-window and annual comparison.
