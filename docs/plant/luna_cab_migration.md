# LUNA -> Vcmax/Cab dynamic coupling: new-chain migration

> Migrate the Route-A LUNA -> dynamic Vcmax + dynamic leaf Cab feedback from the legacy
> `src/COMPASS.m` monolith into the new container chain (`runOneStep` + `+init`).
> Legacy `COMPASS.m` unchanged (frozen). Date: 2026-07-04.
>
> **[owner 更新 2026-07-05~07 增量#6a+6b]** 本文档下方所有 `runtime.luna.state` / `runtime.luna.last_*` /
> `runtime.luna.climate_buffer` 说法是**当时**的 owner;现已**全部**搬入 `Y.plant_acclimation`:LUNA 日优化状态 +
> `last_*` → `Y.plant_acclimation.luna_state`(增量#6a),10 天气候 buffer → `Y.plant_acclimation.luna_climate_buffer`
> (增量#6b,序列化 struct)。`runDailyLuna` 现签名 `[luna_state_next,buffer_struct_next]=runDailyLuna(model,kt,buf)`,读/写
> `Y.plant_acclimation.luna_state`(config 从 `built.p.plant_acclimation`,buffer 由 caller 在 accept 路径 `fromStruct→push→toStruct`);
> 消费者(link 3/4)读 lagged `Y.plant_acclimation.luna_state`。**`runtime.luna` 已整删**。历史文字保留脉络,当前 owner 以 `Y.plant_acclimation.*` 为准。

## 1. Why this needed a real migration (not a 15-line move)

The dynamic-Cab re-fluspect block a reviewer might point at (`COMPASS.m:852-866`) is only the
LAST link. It reads `NitrogenVariables.luna_last_Cab_nlc`, a value **produced** by the daily
LUNA optimization at `COMPASS.m:982-1046` (which runs after ebal). The new chain's
`buildLunaClimate` only *pushed to the climate buffer*; it never ran the optimization, so
`luna_last_Cab_nlc` had no producer. Porting the consumer without the producer = dead or wrong
code (the "croot_stem three-stage broken-link" hazard). So the migration is a 4-link chain,
producer-first.

## 2. The lagged loop (no circular dependency)

```
GPP(day N)  --ebal-->  LUNA climate buffer  --daily opt(day N)-->  luna.state + last_Cab_nlc
                                                                          |
                                        (read at step top, day N+1)       v
   biochemical Vcmax  <--use_luna/Vcmax25_LUNA--  leafbio      RTMo <-- leafopt (re-fluspect)
```

- The daily update is a **lagged** producer: day N's optimization feeds day N+1's optics/biochem.
- **First LUNA day is static**: `computeLuna` seeds (Vcmax25=50) on `first_day`, and the whole
  first day the consumers are gated off (`~first_day`). Dynamic behaviour starts day 2+. This
  matches the intuition "first step static, dynamic from the second".

## 3. The four links (legacy -> new)

| # | New-chain site | Legacy | What it does |
|---|---|---|---|
| 1 | `runOneStep.buildLunaClimate` (moved AFTER `runStepEbal`) | `COMPASS.m:624` push + `987-995` backfill | Build the COMPLETE climate-buffer entry in one shot (forcing Ta/RH/CO2/PAR + ebal GPP/Tleaf/rb, **rb=`thermal.rbc`**, see 5.1); pushed ONLY on accept |
| 2 | `runOneStep.runDailyLuna` (accept branch, after the push) | `COMPASS.m:982-1046` | Once/day when the buffer `isWarm`: `assembleLunaInput -> computeLuna`, commit `runtime.luna.state` (+ `first_day=false`) and `last_*` incl. `last_Cab_nlc` (`nLeafToCab`). Accept-only |
| 3 | `init.prepareScopeStep/buildStepLeafBio` (now takes `rt.luna`) | `COMPASS.m:682-686` | When enabled and past the first LUNA day: `use_luna=1`, `Vcmax25_LUNA/Jmax25_LUNA`, `Vcmo` override. `biochemical.m` reads `Vcmax25_LUNA` when `use_luna==1` (and skips fDYL/FN on it) |
| 4 | `runOneStep.applyDynamicCabLeafOptics` (after `applyCanopyLAI`, before `runStepRTMo`) | `COMPASS.m:688-709` + `852-866` | Gated on `cab_from_luna`: override `leafbio.Cab` (co-vary `Cca`) with the lagged `last_Cab_nlc`, re-run fluspect, overwrite the OPTICAL band (`spectral.IwlP`) of `scope_step.leafopt` |

### Accept-only / rollback discipline
Links 1 and 2 run in the ACCEPT branch (climate push, then daily LUNA), so a failed outer step
advances neither the 10-day buffer nor a daily update. Consumers (3, 4) read `runtime.luna.state`
at step top = the last accepted daily value (lagged), so they are naturally safe.

## 4. Design decision "A" (clean, NOT bit-exact on daily GPP)

Legacy pushes climate BEFORE ebal (GPP=0) then backfills GPP after ebal (`COMPASS.m:995`). The
buffer's daily rollover (`Climate10dBuffer.push`, `gpp_yesterday = gpp_today`) fires inside that
pre-ebal push, so on the LAST KT of a day, that hour's GPP lands in the NEXT day
(`gpp_yesterday` = hours 1..23).

The new chain builds the entry ONCE after ebal with real GPP -> the rollover includes the last
hour -> `gpp_yesterday` = hours 1..24. **Chosen "A" (clean one-shot).** Consequence:
`gpp_yesterday` differs from legacy by exactly the last hour of each day (near-sunset, low GPP,
< ~1%), so `Vcmax25`/`Cab` differ slightly. This is a **deliberate fix** of the legacy
pre-ebal-rollover leak, accepted as NOT bit-exact. `Climate10dBuffer` is shared with legacy, so
the buffer's rollover cannot be changed without touching legacy.

## 5. Input mapping notes (verified against legacy)

### 5.1 `rb` maps to `thermal.rbc` (fixes a legacy field-name bug; diverges from legacy)
`ebalClm45` outputs `thermal.rbc` (canopy boundary-layer resistance) and `thermal.rbc_wv` (the
water-vapor Wallace-Verhoef variant) -- there is **no `thermal.rb`**. Legacy `COMPASS.m:987` reads
`getOr(thermal,'rb',60)`, which never finds `thermal.rb`, so the LUNA buffer's `rb` (and `rb10`) is
a **constant 60** in every legacy run -- a latent wrong-field-name bug.
`buildLunaClimate` maps `rb = getField(thermal,'rbc',60)` -> the real canopy boundary-layer
resistance. LUNA uses `rb` as `gb_mol = cf/rb` for **CO2** diffusion (`photosynthesisLuna`), so
`rbc` is correct (the water-vapor `rbc_wv` stays available as a diagnostic/alternative). This is a
**deliberate physics fix that DIVERGES from legacy-as-run** (which used 60) -- the SECOND intended
divergence, alongside the design-"A" GPP shift (Section 4). Pinned by `tests/test_luna_chain.m`
(`buffer rb == scope_step.thermal.rbc`).

### 5.2 Other inputs
- `assembleLunaInput` IGNORES its 2nd arg (N_leaf): LUNA uses `lnc = 1/(slatop*leafcn)`, a PFT
  constant (`Climate10dBuffer.m:248`, CLM5 `lnc_opt=.false.`). So no prognostic N pool needed.
- `canopy.dayl` is never assigned anywhere in `src` -> legacy silently uses `dayl_sec = 12*3600`
  (`COMPASS.m:1008`); the new chain uses the same constant.
- `rabsorb = 1 - leafbio.rho_thermal - leafbio.tau_thermal`; `slatop=SLA0`, `leafcn=CN_L` from
  the PFT table; `nrad=1` (`COMPASS.m:1016-1019`).
- `optipar` from the static cache: `built.p.scope_static.library.optipar`; `spectral` / `options`
  (`calc_PSI`) from `built.p.scope_static` (moved from runtime.scope_static 2026-07-05 increment 4).
  `fversion` selection mirrors `precomputeScopeOptics`.

## 6. Validation status

| Check | Result |
|---|---|
| Regression `test_runOneStep` (LUNA off/first day -> links inert) | 103/0 |
| Regression `test_rtmo_parity` (RTMo leaf-optics path) | 114/0 bit-exact |
| Activation smoke, natural 24-step run | link 2 fires (`first_day 1->0`), link 3 `use_luna=1` |
| Activation smoke, forced valid LUNA state (1 step) | link 3 `Vcmax25_LUNA=60`; link 4 `leafopt.refl` optical band changes 0.104 |

**Pending (rigorous):** per-link COMPONENT parity vs legacy (feed matched buffer/state/leafbio,
assert `Vcmax25`/`Cab`/`leafopt` match within design-A rel-tol). Links 2 & 4 REUSE the exact
legacy functions (`computeLuna`, `nLeafToCab`, fluspect) with 1:1-mapped inputs, so given matched
inputs the outputs are bit-identical; the parity test would lock the INPUT-assembly wiring.

**Why no live end-to-end golden yet:** the substep soil/heat/N solve is still a stub
(`attemptSolve` no-op), so the new chain's soil state is frozen and its multi-day GPP diverges
from a full legacy run independent of the LUNA code. End-to-end LUNA/Cab parity is only possible
once the real solve lands; until then, component parity is the correct golden.

## 7. Intentional numerical divergences from legacy-as-run

Everything else in this migration reuses the exact legacy functions (`computeLuna`, `nLeafToCab`,
fluspect, `ebalClm45`) and is bit-exact given matched inputs. The following are the **only two
places** the new chain intentionally produces DIFFERENT numbers from a legacy run. Both are
deliberate improvements (fixing latent legacy defects); list them when explaining LUNA/Cab
differences in GMD / parity work.

| # | Quantity | legacy-as-run | new chain | Why (intent) | Magnitude |
|---|---|---|---|---|---|
| **D1** | daily GPP total `gpp_yesterday` (LUNA `gpp_day`) | hours 1..23 (last hour of each day leaks to the NEXT day) | hours 1..24 (full day) | Design "A": `buildLunaClimate` builds the entry once AFTER ebal, so the buffer's daily rollover includes the current hour. Legacy pushed before ebal (GPP=0) then backfilled, so the rollover fired before the last hour's GPP was added. Fixes the pre-ebal-rollover leak. | < ~1% on `Vcmax25`/`Cab` (leaked hour is near-sunset, low GPP) |
| **D2** | LUNA buffer `rb` (-> `rb10`) | constant **60** [s/m] | real `thermal.rbc` (canopy boundary-layer resistance) | Legacy `getOr(thermal,'rb',60)` reads `thermal.rb`, a field ebal never sets -> silent fallback to 60 (wrong-field-name bug). LUNA uses `rb` as `gb_mol=cf/rb` for CO2 diffusion, so the real `rbc` is correct (not the water-vapor `rbc_wv`). | `rb` now varies with canopy state instead of 60 -> shifts the LUNA CO2-supply term |

Both diverge from **legacy-as-run** (what the legacy code actually computes), not from legacy
intent. Neither is a bug in the new chain. Pinned by `tests/test_luna_chain.m` (D1 indirectly via
the GPP-from-fluxes check + the accept-only rollover; D2 directly via `rb == thermal.rbc`).
