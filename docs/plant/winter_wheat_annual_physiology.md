# Winter-wheat annual leaf and mobile-carbon physiology

> **版本状态：2026-09-30用户缩小范围后的历史设计。** 本文保留回退前的设计、文献来源和验证记录。按用户最新范围，APSIM茎WSC、站立枯叶、实际叶氮LUNA反馈及新增冠层fail-fast属于撤回项；地上保留冬小麦物候、CoLM植被C/N分配和原CoLM凋落物，其它沿用COMPASS。本文的历史测试不能解释为当前生产验证；当前回退状态和后续检查以[本轮范围回退报告](winter_wheat_scope_rollback_20260930.md)为准。

Editor: Wenjie Liu; Date: 2026-09-30.

## Evidence and scope

This design adds annual green-leaf / standing-senescent-leaf bookkeeping and
annual stem water-soluble carbohydrate (WSC) supply to grain. It does not
restore the disabled perennial wood-carbon reserves or next-year nitrogen
storage. Changwu observations are unchanged; reference parameters below are
not claimed as site measurements or cultivar calibration.

Primary mechanistic evidence is Takahashi, Chevalier and Rupp (2001),
*Storage and Remobilization of Soluble Carbohydrates after Heading in Different
Plant Parts of a Winter Wheat Cultivar*, Plant Production Science 4:160–165,
DOI [10.1626/pps.4.160](https://doi.org/10.1626/pps.4.160). Their winter-wheat
measurements identify stem fructan/sucrose accumulation during early filling
and depletion during late filling. This establishes an annual mobile-carbon
process, rather than perennial woody storage.

Reference mechanisms inspected on 2026-09-30:

- [APSIM Classic Wheat.xml](https://github.com/APSIMInitiative/APSIM-Classic-Wheat-Doc/blob/b2dfe824c1f2c7fb22e13be2cd4ae3b851427615/Wheat.xml):
  leaf dry-matter senescence is proportional to canopy area senescence
  (`x_dm_sen_frac_leaf = 0 1`, `y_dm_sen_frac_leaf = 0 1`); standing leaf
  detachment defaults to zero. `stemGrowthStructuralFraction = .70 0.0`
  supplies an explicitly identifiable reference partition. COMPASS uses
  0.70 structural / 0.30 mobile in vegetative new-stem growth, without
  claiming APSIM phenological-stage parity.
- [APSIM Next Generation Wheat.json](https://github.com/APSIMInitiative/ApsimX/blob/865099771e8b88981b350139bd7c57824e88c5e7/Models/Resources/Wheat.json):
  `Stem/DMRetranslocationFactor` is zero before the reproductive interval and
  0.5 in its reproductive interval. This supplies the reference daily
  mobilizable proportion, not a measured Changwu rate.
- [APSIM RetranslocateNonStructural.cs](https://github.com/APSIMInitiative/ApsimX/blob/865099771e8b88981b350139bd7c57824e88c5e7/Models/PMF/RetranslocateNonStructural.cs):
  `CalculateBiomass` bounds supply by start-of-round storage, after any
  reallocation, times the retranslocation factor. Actual transfer is
  separately debited from the donor.
- [APSIM GenericOrgan.cs](https://github.com/APSIMInitiative/ApsimX/blob/865099771e8b88981b350139bd7c57824e88c5e7/Models/PMF/Organs/GenericOrgan.cs):
  green and dead pools are distinct; senescent mass moves internally to dead
  material, while harvesting/detachment are separate external transactions.
- Existing COMPASS grain-fill retained leaf C:N = 65 originates from
  [CoLM MOD_Const_PFT.F90](https://github.com/CoLM-SYSU/CoLM/blob/626347a9156ab97db4ca846f47cd1b7f2728f705/main/MOD_Const_PFT.F90).
  Reusing it for the residual N in standing senescent material is a declared
  COMPASS surrogate, not a published dead-leaf C:N measurement.

### Local CoLM winter-wheat cross-check

The user's local CoLM source was read from
`C:/Users/Administrator/Desktop/CoLM/CoLM-master/main/BGC/MOD_BGC_Veg_CNPhenology.F90`,
`MOD_BGC_Veg_NutrientCompetition.F90` in that BGC directory, and
`main/MOD_Const_PFT.F90`. Unmodified copies and SHA-256 hashes are retained in
`outputs/review_aboveground_implementation_20260930/annual_colm_local_evidence/manifest.json`.
This local audit does not establish an upstream revision for that checkout.

`CropPhenology` lines 1063–1093 calculates winter-wheat maturity heat units as
`0.42*gdd1020 + 440`, with `HUI=gddplant/gddmaturity`. `gdd1020` is the
cross-year climate quantity updated in `CNPhenologyClimate` lines 317–345;
for the northern hemisphere it aggregates April–September heat above 10°C.
It is not the current sowing-to-date crop heat accumulation. Crop heat
`gddplant` uses the crop base temperature and a vernalization multiplier
in phase 2. Thus the maturity formula cannot be transplanted by substituting
the current crop's accumulated heat for the climatology input.

The 0-based winter-wheat entry 21 has `baset=0`°C, `grnfill=0.67` HUI,
`mxmat=270` days after planting, `arooti=0.1`, `arootf=0`, `fleafi=0.9`,
`bfact=0.1`, and `astemf=0.05`. These are CoLM reference parameters, not
observations for Changwu. `CropPhenology` distinguishes planting/emergence,
vegetative growth, and grain filling; its original `HUI>=1` or `idpp>=mxmat`
branch directly harvests. COMPASS intentionally separates physiological
maturity from harvest at the observed LAI-zero event. The reference's
planting seed is 3 g C m^-2, rather than a woody reserve.

`CNAllocation` in `MOD_BGC_Veg_NutrientCompetition.F90` lines 381–466
provides stage-dependent allocation and a one-time leaf/stem/root N
retranslocation to the mobile N pool at filling. The current COMPASS retained
N limits have the corresponding `fleafcn`, `fstemcn`, and `ffrootcn` sources.
Its `cpool_to_grainc_storage` at line 190 receives new allocation into the
grain organ's storage compartment; it is not an N-free stem WSC donor
supplying grain. No corresponding independent stem WSC owner or stem-to-grain
carbon retranslocation was found in the local BGC call path. CoLM therefore
does not supply replacements for the reference 0.30 WSC allocation or
0.50/day WSC mobilization used here.

`CNBackgroundLitterfall` lines 1403–1423 transfers green leaf carbon to soil
litter and partitions leaf N between litter and retranslocation. It has no
independent standing-dead leaf owner. The winter-wheat `lflitcn=25` controls
that litter flux, whereas `fleafcn=65` is the live-leaf N retention limit
during grain filling. Neither is an independently measured standing-dead leaf
C:N. Using 25 in the standing pool would change the retained N and grain N
supply; it cannot be presented as an equivalent source-driven replacement.
The 65 surrogate remains explicit and requires independent senescent-leaf
N measurements or sensitivity analysis before a quantitative site claim.

## State, forcing and parameters

All plant pools are ground-area inventories: C in g C m^-2, N in g N m^-2.
Existing `C.leaf` and `N.leaf` represent green tissue. New
`C.leaf_senescent` and `N.leaf_senescent` represent standing senescent leaf
material. New `C.stem_wsc` is nitrogen-free mobile carbohydrate stored in
the annual stem. It is excluded from structural stem nitrogen demand and
from the perennial `root_storage` / `stem_storage` pools.

`crop_green_lai_previous` is accepted runtime memory of the canopy-area
forcing used by this process, with an initialization marker. It is physical
event memory included in snapshots and checkpoints, not a diagnostic driving
physics. A failed trial never advances it.

Compiled parameters are the sole process input. The current reference
defaults are `stem_wsc_alloc_fraction = 0.30`,
`stem_wsc_retrans_daily_fraction = 0.50`, and the existing retained leaf
C:N through `senescent_leaf_cn`. These are explicitly uncalibrated. The
existing growth respiration fraction supplies the grain conversion yield;
no separate unexplained transport-cost or grain-number parameter is added.

## Conservative closure

Let L0 be the previous accepted positive green LAI and L1 the current
consumed green LAI. Prescribed canopy area decline gives

    s = max(0, 1 - L1/L0), bounded by [0,1].

This applies only after the canopy has had a positive area. Zero measured
LAI before canopy appearance does not destroy the sub-observation seed leaf.
On a physical harvest event, shoot removal takes precedence. Increasing or
unchanged area does not senesce tissue. This is an area-constrained canopy
senescence closure, not a claim to model leaf cohorts, drought-triggered
senescence or detached litter independently.

Green leaf carbon dC = s*C_leaf moves to standing senescent leaf carbon.
Green leaf nitrogen dN = s*N_leaf_start is removed from the same accepted
tissue fraction; start-of-step N excludes any explicit grain-fill donor
transfer already applied. Newly allocated C/N is ineligible in this round.
Standing tissue receives min(dN,dC/CN_senescent); the balance enters mobile
N. Thus total plant C and N are invariant under senescence. No litter or
CO2 is invented by the internal change. Already depleted leaf N cannot be
exported a second time.

New vegetative stem allocation is split into structural carbon and N-free
WSC before the structural C:N ceiling. The matching N-demand equation uses
only the structural share. WSC allocation receives its own applied result,
and is included in the same carbon budget as structural growth. Existing
stem material is not retrospectively relabelled or fabricated as WSC.

During active filling, the hourly/subdaily mobilizable fraction is

    f(dt) = 1 - (1-f_day)^(dt/86400).

This integrates the reference daily proportion and avoids applying 0.5 on
every hourly step. Supply is bounded by the accepted pre-step WSC inventory;
new WSC is not loaded and discharged in the same step. The stoichiometric
grain sink is bounded by both existing grain-N carbon capacity and physical
mobile N. Carbon transferred to newly formed grain incurs the existing
growth-respiration conversion loss. Donor debit equals grain-C credit plus
this actual atmospheric C flux. Mobile-N debit equals grain-N credit.

COMPASS lacks grain number and potential individual grain size in the
current state. Therefore the sink is a declared stoichiometric adaptation,
not APSIM's independent per-grain demand closure. It cannot support a
validated absolute potential-yield claim. A future grain-number/sink module
must be supported by cultivar/site evidence and independently validated.

Harvest removes standing senescent leaf C/N and residual stem WSC with
the existing Changwu standing-shoot export. Root residue and retained seed
remain under their existing contracts. The mobile-sugar cleanup flux remains
explicit. New internal-transfer accumulators are diagnostic accepted-state
ledgers only; they never enter external mass-input/output sums.

## Validation gates

Focused checks cover: unchanged/increasing LAI, known 50% area loss, no
seed-leaf loss before a positive canopy, nonnegative bounded donor debits,
leaf N already below residual requirement, C/N closure, N-free WSC loading
with consistent structural N demand, zero GPP with WSC and N continuing
grain growth, no grain growth without C or N, explicit respiratory loss,
daily/subdaily retranslocation equivalence for unconstrained supply,
harvest once with new pools, deterministic retry, accepted-state adapter
roundtrip, initialization inventory, and carbon/nitrogen budget inclusion.
Existing lifecycle/funding checks retain their tolerances. Short coupled
checks and restart checks follow integration; no annual run is required to
claim only these local mechanisms are verified.

The accepted soil N amount also constrains the final crop candidate. A
requested-uptake plant calculation is replayed from its original accepted
input when the applied soil uptake differs. The replay uses the same frozen
demand, time interval and events, and asserts that already-consumed external
litter/harvest ledgers are unchanged. Merely subtracting the shortfall from
mobile N fails when the new WSC/grain sink has already consumed that N. Replay
inputs are transient transaction material, not checkpoint or physical-state
owners; non-crop candidate materialization retains its original behavior.

Verification artifacts after these changes:

- `annual_physiology_gate_06/verification.json`: 47 pure process/state gates;
  initialization, adapters, C/N/CO2 closure, N limits, maturity, harvest and retry.
- `actual_root_N_gate_03/verification.json`: 29 actual-uptake replay gates;
  an explicit bounded soil-interface shortfall reproduces the negative mobile-N
  failure of the previous correction, while the replay closes mass. It also
  compares non-crop materialization with the frozen original source bitwise.
- `original_lifecycle_gate_02/verification.json`: 73 original lifecycle gates;
  original pre-edit inputs, references and tolerances are preserved.
- `original_funding_gate_03/verification.json`: 106 original carbon-funding
  gates, preserving the original reference cases and tolerances.

All paths above are under
`outputs/review_aboveground_implementation_20260930/`. Each named run completed
with MATLAB exit 0. The two new suites supply synthetic accepted stocks and
explicit interfaces and execute zero coupled soil-model steps. They establish
conservation and transaction behavior, not site-specific rates, yield or
field validation. Separate coupled accepted-step and restart evidence must
be read from the final implementation audit.
