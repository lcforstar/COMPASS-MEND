# Labile sugar floor for COMPASS NSC allocation

## Purpose

This note defines a scientifically defensible way to stop the COMPASS carbon
allocator from emptying the soluble-sugar pool (`C_labile`) into starch-like
storage every time step.

The intended pool semantics are:

- `C_labile` = soluble sugars and immediately mobile carbon, used for
  respiration, transport, osmotic adjustment, signalling, and near-term growth.
- `C_root_storage` and `C_stem_storage` = starch-dominated reserve carbon,
  seasonally buffered and remobilisable.

The allocation rule should therefore protect a labile sugar floor first, then
send only true excess carbon to storage. When labile sugar falls below the
floor, storage should remobilise back to labile carbon before any hard
starvation/debt guard is used.

This mechanism does not, by itself, fix false substrate down-regulation of
respiration. Respiration down-regulation and labile-pool realism must be treated
as two separate controls:

- Respiration substrate control: read total remobilisable NSC
  (`C_labile + C_root_storage + C_stem_storage`) so starch reserves can support
  maintenance respiration after remobilisation.
- Labile sugar floor: maintain a small, realistic soluble-sugar working pool so
  `C_labile` is not emptied into starch-like storage every step.

With this separation, respiration down-regulation becomes a starvation safety
valve controlled by total available carbohydrate, while the sugar floor controls
pool partitioning and future sugar-dependent processes such as osmotic
adjustment, cold hardening, and microwave/VOD interpretation.

## Literature basis

1. Martinez-Vilalta et al. (2016) compiled seasonal NSC measurements from
   121 studies and 177 species. Their synthesis reports common starch depletion
   but rare depletion of total NSC or soluble sugars. Minimum total NSC was
   still high, averaging about 46% of seasonal maxima. The authors interpret
   this as dual NSC function: starch is mostly a reservoir, while soluble sugars
   support immediate functions and remain above a critical threshold.
   Verification status: local PDF present at
   `data/references/fruit_tree/Martinez-Vilalta_et_al-2016-Ecological_Monographs.pdf`.
   DOI: https://doi.org/10.1002/ecm.1231

2. Hartmann and Trumbore (2016) emphasise that NSC pools should not be treated
   as one undifferentiated storage pool. Sugars and starch have distinct roles,
   and NSC storage buffers diel, seasonal, and inter-organ mismatches between
   carbon supply and demand.
   Verification status: DOI/Crossref metadata verified only; full text was not
   found in the local reference library during this note. Do not cite this in
   code comments until the full text is downloaded or otherwise verified.
   DOI: https://doi.org/10.1111/nph.13955

3. Thalmann and Santelia (2017) review starch remobilisation under abiotic
   stress. Starch degradation releases sugars and derived metabolites that help
   plants maintain function under stress.
   Verification status: DOI/Crossref metadata verified only; full text was not
   found in the local reference library during this note. Do not cite this in
   code comments until the full text is downloaded or otherwise verified.
   DOI: https://doi.org/10.1111/nph.14491

4. Reyes et al. (2016) is the apple-specific source/sink study to keep in the
   calibration evidence set. Its Crossref metadata confirms the paper and
   apple resource-limitation context; apple-specific soluble sugar/starch
   parameterisation should still be tied to measured orchard data or full-text
   extraction before hard numeric defaults are claimed.
   Verification status: local PDF present at
   `data/references/fruit_tree/MaxGrowth_2016_AppleReserve_FrontiersPlantSci.pdf`.
   DOI: https://doi.org/10.3389/fpls.2016.00233

## Do not use a fixed sugar:starch ratio as the first implementation

A fixed soluble-sugar:starch ratio is not the right first closure because:

- soluble sugars and starch have different turnover times;
- the ratio is organ-specific and seasonal;
- apple winter cold acclimation can raise soluble sugars while starch falls;
- the Martinez-Vilalta 46% value applies to minimum total NSC relative to
  seasonal maximum, not directly to soluble sugar alone.

The 46% result should be used as a sanity check on total NSC depletion, not as
`C_labile = 0.46 * NSC`.

## Preferred hierarchy for setting the labile floor

### Level 1: observation-constrained floor

Use this whenever measured soluble sugar data exist.

For organ `o`:

```text
C_labile_target_obs(t) = sum_o qSS_o(season, organ) * B_live_o(t)
```

where:

- `qSS_o` is the low-end soluble sugar concentration, preferably the 5th or
  10th percentile of measured soluble sugar carbon per live organ carbon;
- `B_live_o(t)` is live organ carbon in COMPASS units `[g C m-2]`;
- the season bins should be broad enough to avoid overfitting, for example
  dormant, early growth, fruit expansion, post-harvest, and winter hardening.

This is the most defensible option because the threshold is tied to measured
soluble sugar, not inferred from starch or total NSC.

### Level 2: respiration-buffer floor

Use this as the default when no soluble sugar observations are available.

```text
C_labile_target(t) = max(C_labile_abs_min,
                         labile_resp_days * Rm_potential_day(t))
```

where:

- `C_labile_target` is `[g C m-2]`;
- `C_labile_abs_min` is a small nonzero numerical/physiological floor
  `[g C m-2]`;
- `labile_resp_days` is the number of days of potential maintenance
  respiration protected as soluble sugar;
- `Rm_potential_day = R_m_potential_kg * 1000 * 86400` is `[g C m-2 d-1]`.

`R_m_potential_kg` must be the temperature/N-driven maintenance respiration
before substrate down-regulation. Do not derive the sugar floor from the
already-throttled `R_a_cpool_kg` passed into `computeCarbonAllocation`; doing so
creates a circular dependency in which a false down-regulation shrinks the
target sugar pool.

Recommended first sensitivity set:

```text
labile_resp_days = [0.25, 0.5, 1.0, 2.0]
```

Do not lock a single value before checking model behaviour. The smallest value
that prevents persistent zero `C_labile` while preserving realistic storage
seasonality is preferable.

### Level 3: optional seasonal modifier

Add only if observations or validation show that the respiration-buffer floor
is insufficient.

```text
C_labile_target = C_labile_target_base * f_season
```

with a simple bounded modifier, for example:

- `f_season > 1` during winter hardening or cold dormant periods;
- `f_season = 1` otherwise.

This should remain optional because it adds a parameter and can hide allocation
errors if introduced too early.

## Proposed code logic

This change has two code parts. They solve different problems and should be
implemented together.

### Part A: respiration down-regulation reads total NSC

Modify `src/+plant/+integrated/applyNSCCoupling.m`, Step 1b.

The current substrate term uses `C_labile` for the Michaelis-Menten
down-regulation:

```matlab
Cl_prev = max(C_pools_prev.C_labile, 0);
R_m_g = R_m_g * (Cl_prev / (resp_Km + Cl_prev));
```

This makes the down-regulation fire even when starch-like storage is abundant.
The substrate for respiration down-regulation should instead be total
remobilisable NSC:

```matlab
C_resp_substrate = max(C_pools_prev.C_labile, 0) + ...
                   max(C_pools_prev.C_root_storage, 0) + ...
                   max(C_pools_prev.C_stem_storage, 0);
R_m_g = R_m_g * (C_resp_substrate / (resp_Km + C_resp_substrate));
```

The same substrate concept must be used for any substrate-proportional
respiration option that currently treats `C_labile` as the only fast pool,
unless that option is explicitly reinterpreted as sugar-only respiration.

### Part B: allocation protects the labile sugar floor

Modify `src/+plant/+PlantCarbon/computeCarbonAllocation.m`, Step 3c.

Current logic:

```matlab
avail   = C_pools_prev.C_labile + GPP_in_KT - R_a_out_KT - display_total;
surplus = max(avail, 0);
```

Proposed logic:

```matlab
C_labile_target = computeLabileSugarTarget(...);       % [g C/m2]
avail = C_pools_prev.C_labile + GPP_in_KT - R_a_out_KT - display_total;
surplus = max(avail - C_labile_target, 0);
```

After storage refill:

```matlab
C_labile_pre = avail - C_to_storage_tot;
labile_deficit = max(C_labile_target - C_labile_pre, 0);
C_storage_to_labile = min(labile_deficit, C_root_storage + C_stem_storage);
```

Then update:

```matlab
C_labile_new = C_labile_pre + C_storage_to_labile;
C_storage_new = C_storage_prev + C_to_storage_tot - C_storage_to_labile;
```

Root/stem storage drawdown should be split by available storage, not by a fixed
preference that empties one pool first, unless there is a physiological reason
to impose priority.

### Coordination with `applyNSCCoupling` Step 5

`applyNSCCoupling.m` already has a hard post-allocation guard:

- if `C_labile < 0`, it raises `C_labile` to zero;
- in `resp_scheme == 4`, it draws from `C_root_storage` and `C_stem_storage`;
- any residual deficit becomes `C_xsmr` debt.

The new sugar-floor remobilisation should occur inside
`computeCarbonAllocation` before that guard and should target
`C_labile_target`, not merely zero. Step 5 should remain as a last-resort
conservation guard for residual negative `C_labile` after allocation and
remobilisation. To prevent double drawdown:

1. `computeCarbonAllocation` should return storage pools after any
   `C_storage_to_labile` remobilisation.
2. `applyNSCCoupling` Step 5 should compute `cpool_deficit_g` from that returned
   post-remobilisation `C_labile`.
3. If the sugar-floor remobilisation has already kept `C_labile >= 0`, Step 5
   draws nothing.
4. If storage cannot maintain even zero `C_labile`, Step 5 handles only the
   residual negative amount.

This gives a clear order:

```text
respiration substrate check on total NSC
-> allocation and storage refill with labile floor protected
-> storage-to-labile remobilisation toward labile floor
-> hard zero/debt guard only for remaining deficit
```

## Diagnostics required

Add diagnostics so the mechanism can be audited:

- `C_labile_target` `[g C m-2]`
- `C_labile_pre_remob` `[g C m-2]`
- `C_labile_deficit_before_remob` `[g C m-2 timestep-1]`
- `C_storage_to_labile` `[g C m-2 timestep-1]`
- `C_root_storage_to_labile` `[g C m-2 timestep-1]`
- `C_stem_storage_to_labile` `[g C m-2 timestep-1]`
- `resp_substrate_g` `[g C m-2]`
- `resp_beta_sub` `[-]`, after switching the substrate to total NSC

These are diagnostics only; physical code should not read them back.

## Acceptance checks

Before accepting this change:

1. Carbon conservation remains at machine precision:

```text
dC_labile + dC_xsmr + dC_organs + dC_storage = GPP - R_a
```

2. `C_labile` no longer spends long periods exactly at zero in ordinary
   growing-season conditions.

3. Storage remains seasonal and bounded rather than monotonically exploding.

4. Starch-like storage is allowed to decline during night, stress, or source
   limitation, while soluble sugar is buffered.

5. The respiration down-regulation signal should no longer fire just because
   `C_labile` is intentionally small while starch-like storage is abundant.
   This is checked through `resp_substrate_g` and `resp_beta_sub`, not by
   increasing `C_labile_target` until the Michaelis-Menten term stops firing.

6. Sensitivity over `labile_resp_days = [0.25, 0.5, 1.0, 2.0]` is reported
   before choosing a default.

7. Annual mean autotrophic respiration is not expected to move strongly from
   this change alone. This change corrects NSC pool semantics and false
   substrate down-regulation; it is not the fix for an overly high baseline
   respiration formulation such as leaf `Rd` or `r_m10`.

## Initial recommendation

Implement both code parts together:

1. change respiration substrate down-regulation to total remobilisable NSC;
2. add the labile sugar floor and storage-to-labile remobilisation for pool
   realism.

For the sugar floor, implement Level 2 first because it is parsimonious and
uses existing model fluxes. Treat `labile_resp_days` as a calibration parameter,
not a known constant. Use Level 1 as soon as measured soluble sugar data are
available.
