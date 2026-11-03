# Mineral-fertilizer nitrogen source tracer

Change: 2026-10-28 Wenjie Liu — `source_hourly + mend_dynamic` is supported
with tracer schema 5. The explicit `dom.don_sorbed_g_cm3` label is a subpool
of `dom.don_total_g_cm3`, never another inventory added to column storage.
Linear-equilibrium tracing retains schema 4. Existing dynamic runs with
tracing disabled contain no historical fertilizer-source information and
cannot retrospectively acquire labels by switching this flag at resume.

Change: 2026-09-29 Wenjie Liu — current schema 4 also records organic top
terminal losses. Prescribed external DON/particle inputs remain unlabelled;
they dilute the source fraction and never create fertilizer labels. See
[organic boundary types](ORGANIC_BOUNDARY_TYPES_20260929.md). DIN equations
and boundary interfaces are unchanged by this organic-only extension.

Change: 2026-09-29 Wenjie Liu — tracer state schema 3 additionally accounts for
DON and disjoint MB/BD/intracellular-ENZN bottom terminal losses under their
respective open boundaries. These are not internal organic/mineral recycling.
See [particle boundary acceptance](BIOPARTICLE_BOTTOM_BOUNDARY_IMPLEMENTATION_20260929.md)
and [DOM boundary acceptance](DOM_BOTTOM_BOUNDARY_IMPLEMENTATION_20260929.md).

## Purpose and invariant

`fertilizer_source_tracer` is a passive mass label for nitrogen introduced by
the mineral-fertilizer schedule. It answers how much fertilizer-source N is in
each mineral, dissolved-organic, particulate-organic, microbial, and enzyme
pool, and records its accepted gross entry into and release from organic pools.

The tracer is downstream of the total-N solution. It cannot supply a rate,
coefficient, boundary condition, convergence criterion, or timestep to the
water, heat, gas, carbon, microbial, or total-nitrogen equations. The production
A/B regression requires every non-tracer dynamic field to be bitwise identical
with the tracer switched on and off.

Enable it in the site configuration:

```toml
[switches]
fertilizer_source_tracer = true
```

The current schema requires `solve_nitrogen`, native soil C, DOM, and explicit
microbes, with `microbe.reaction_method = "source_hourly"` and
`dom.sorption_model = "linear_equilibrium"` or `"mend_dynamic"`.
The RK4 reaction path remains unsupported by this source-hour replay.

## Source definition

One aggregate cohort labels every event in the mineral-fertilizer schedule,
including its exact liquid NH4, liquid NO3, urea, and solid NH4-granule
increments. Native soil N, atmospheric deposition, groundwater N, biological
fixation, litter, straw, manure, organic-fertilizer organic N, and
organic-fertilizer TAN are unlabelled.

This is a source label, not a simulated isotope. It assumes perfect mixing in
each donor pool and no source preference or isotope fractionation. It therefore
must not be interpreted as a mechanistic 15N fractionation model.

## Accepted calculation order

For each accepted outer step, the label replays the already accepted total
trajectory in the same physical order:

1. Stage the exact nodal mineral-fertilizer increments.
2. Attribute accepted DON transport using its frozen total endpoints and
   directed face-exchange masses, with bounded implicit label fractions.
3. Attribute root-first DIN reservation from the event-preceding accepted pool.
4. Replay source-hour MEND gross transfers through PON, MON, DON, active and
   dormant microbes, extracellular enzymes, intracellular enzyme phases, and
   the two DIN uptake stages.
5. Attribute accepted BIOPARTICLE mobile transport in the same way; interleave
   the exact total-N attachment, detachment, and straining transfers in order.
6. Attribute urea hydrolysis and granule dissolution from their donor-label
   fractions.
7. Solve NH4 and NO3 labels with the final frozen total-N ADE matrices,
   boundary types, transport coefficients, and accepted gross rates; apply a
   bounded conservative co-projection so that `0 <= label <= total`.
8. Commit the label, total state, management flag, and route ledgers in one
   accepted transaction.

For a gross accepted transfer `G` from donor pool `D`, its labelled transfer is

```text
G_label = G * label(D) / total(D)
```

using the donor fraction at the corresponding accepted gross-route stage. Net
pool changes are never used to infer two simultaneous opposing routes.

## State and balance

### Dynamic DON phase attribution (2026-10-28)

For dynamic sorption the model owns total DON and its explicit sorbed
subpool. Tracing uses the same ownership: `free = total - sorbed`.
Transport replays the accepted **free-DON** exchange ledger, leaving the
immobile sorbed label untouched. Source-hour uptake reads the free donor
label; it does not infer its fraction from total DON.

The independently accepted adsorption amount `A` and desorption amount `D`
exchange labels between the free pool `F` after its accepted uptake and
all same-hour gross organic inputs, and the starting sorbed pool `S`.
This includes same-hour inputs because the physical source map can borrow
from them and requires nonnegative final phases, without requiring a
nonnegative sorption intermediate before those inputs. Bounded mixing solves

```text
[F + D,   -D] [f_free  ] = [label_free_after_non_sorption_routes]
[   -A, S + A] [f_sorbed]   [label_sorbed_start     ]
adsorption_label = A * f_free
desorption_label = D * f_sorbed
```

This is a passive attribution discretization on the accepted gross
transfers, matching the implicit-mixing policy used for transport; it does
not alter the physical sorption algorithm. It permits simultaneous gross
recycling without overdraw and preserves a uniform source fraction. Both
gross routes are exported separately. Their net transfer is internal to
total DON and is not a new source or terminal loss. Free and sorbed label
bounds, reconstructed total endpoints, and per-layer label closure are
checked independently. Schema-4 checkpoints are rejected by the existing
compatibility-signature check when the requested dynamic tracer is schema 5;
no historical sorbed labels are invented.

If the source map emits a signed sorption route (for example, negative
adsorption above its `Qmax`), tracing preserves that sign. The mixing solve
uses directed transfers `A+ = max(A,0) + max(-D,0)` and
`D+ = max(D,0) + max(-A,0)`; each original signed route receives the
fraction of its corresponding donor. No physical gross amount is clipped.
This keeps source-supported final phases traceable even when nominal
adsorption reverses direction.

Nitrification tracing consumes the accepted frozen total nitrification
amount. It automatically inherits the physical moisture and oxygen limits:
98% enters nitrate and 2% enters nitrification N2O. A second moisture factor
must not be applied to the label.

### Transport attribution update (2026-09-22)

The explicit total transport can accept a nonnegative total endpoint without
being positivity-preserving for arbitrary source-label distributions. A sampled
Changwu case exposed this at KT=84. Merely replaying the total's accepted time
grid for labels was therefore insufficient. The total solver, its accepted time
grid and physical states remain unchanged by this repair.

Each accepted tile now exports before/after total concentrations and forward/
backward face masses. The face masses preserve the original accepted net flux
and its opposing gross advection/dispersion components, including dispersion
mixing at zero net total flux. In layer-mass units, label fractions satisfy:

```text
(diag(m_new + outgoing) - incoming_exchange_matrix) f_new = label_mass_old
label_mass_new = m_new .* f_new
```

This matrix has nonpositive off-diagonals, row sums `m_old`, and column sums
`m_new`. Thus it preserves uniform fractions, boundedness and label mass without
altering totals. Implementation uses the algebraically equivalent diagonal
`m_old + incoming` to retain row dominance, checks total-ledger roundoff and
independently audits label conservation. Only fraction excursions within a
floating-point roundoff bound may be snapped to [0,1], with the resulting mass
recorded; material excursions fail. No clipping of materially negative labelled
concentrations or global mass redistribution is used.

This is an implicit source-attribution discretization, not an exact replay of
the old explicit label matrix. It can change labelled spatial mixing and the
resulting route metrics. It assumes within-tile implicit mixing and retains
the same perfect-mixing/no-source-preference interpretation. Refinement studies
are still needed for quantitative path sensitivity; a short regression does
not establish full-year convergence. Historical tracer results and campaigns
must not be mixed with results from this changed tracer algorithm.

Microbial phase exchanges use frozen accepted **nitrogen** transfers, avoiding
a synthetic-carbon reconstruction and keeping C/N ownership explicit.

The extended sampled run also exposed a separate reaction-phase attribution
defect: scaling phase labels independently to a new aggregate can exceed a
phase's resized total. Phase attribution now first carries each label with
the total phase resize. A decrease in the aggregate labelled target removes
labels proportionally; an increase is assigned proportionally to each phase's
remaining unlabelled capacity. The complementary label receives the
complementary allocation. This preserves the gross-derived aggregate target,
zero/full-label limits and every phase bound. It is a numerical attribution
remap, not an added physical interphase transfer. The allocation and residual
are explicitly retained under the source-replay `phase_attribution` ledger.
The old independent label-rescaling rule is no longer used.

Change: 2026-09-22 Wenjie Liu — bounded transport attribution contract.

The accepted owner is `Y.soil_fertilizer_n_tracer`. Mineral pools retain the
total solver's bottom-up node representation in `ug N cm-3`; DON and organic
pools retain bottom-up layer values in `g N cm-3`. Microbial mobile, attached,
and retained phases are explicit. Checkpoints bind the enabled flag, schema,
and immutable source definition and carry the complete accepted label state.

The column balance is

```text
fertilizer input + labelled boundary inputs
  = current labelled storage + labelled terminal losses
    + mirrored numerical correction + closure residual.
```

Terminal losses include root uptake, NH3 volatilization, nitrification N2O,
denitrification N2O/N2, and labelled boundary losses. Organic entry and release
are internal gross crossings and therefore are reported but do not appear as
terminal balance terms.

Production diagnostics add 58 switch-gated fields to `soil_nitrogen.nc`: four
mineral node profiles, ten non-overlapping DON/organic layer profiles, five
column storages, step and cumulative route ledgers, and MEND, DIN, and column
closure residuals. In particular:

- `fertilizer_tracer_organic_entry_gross_{step,cumulative}`
- `fertilizer_tracer_organic_release_gross_{step,cumulative}`
- `fertilizer_tracer_organic_storage`
- `fertilizer_tracer_total_storage`
- `fertilizer_tracer_cumulative_closure_residual`

Output profiles are written surface-to-bottom; the accepted MATLAB state stays
in its owner-native bottom-up order.

## Interpretation limits

- Entry and release are gross crossings. The same atom can be fixed, released,
  and fixed again, so cumulative gross values are not mutually exclusive fate
  fractions and can exceed a one-pass allocation.
- The aggregate cohort does not separate application years or individual
  fertilizer events.
- After labelled organic N is released and remixes with other fertilizer-source
  DIN, this schema records the release route but does not retain a separate
  “ever organic” history tag for its later fate.
- Total-N freeze/thaw liquid-to-ice partitioning is currently disabled in the
  upstream production kernel. The tracer explicitly fails if a nonzero
  partition route becomes active; liquid/ice label replay must be implemented
  before that total route is restored.
- A short production window proves route execution and conservation, not an
  agronomically important release magnitude. Long-term release inference still
  requires calibrated parameters and observation-based evaluation.

## Verification

Run from the repository root in MATLAB:

```matlab
addpath('src');
addpath('tools/diagnostics');
verifyFertilizerTracerEventStaging;
verifyFrozenFertilizerTracerADE;
verifyFertilizerTracerSourceReplay;
verifyFertilizerSourceTracerProduction;
verifyChangwuSubstepConvergence;
```

The production acceptance spans the Changwu mineral event, subsequent organic
entry and release, tracer-on/off bitwise parity, an uncommitted event candidate,
and a checkpoint/resume continuation. The convergence regression uses the same
3600 s forcing with maximum accepted soil tiles of 3600, 1800, 900, and 450 s
and includes mineral-label, organic-label, and cumulative-route trajectories.
The dated numerical record is
`docs/model_comparison/FERTILIZER_SOURCE_TRACER_ACCEPTANCE_20260921.md`.
