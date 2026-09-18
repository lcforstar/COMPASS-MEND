# Stomatal–photosynthesis decoupling under heat: formal plan (v3)

> **Status (2026-07-16):** Formal decision after ebal I/O cleanup, output review, and literature
> (Jones GCB 2026, Sicangco New Phyt 2026, Marchin GCB 2023, Duursma 2019, Eller SOX 2020).
> **This revision supersedes v2 recommendations that packaged
> `gmin(T) + kmax(T) + 2-D three-objective optimizer` as the default implementation path.**
>
> **Work order now:** record plan here → **do not implement stomatal physics yet**.
> ebal/ebalProfile migration is already done; do **not** reopen a parallel “migration”
> solely for this topic. Implementation starts only after Stage 0 diagnostics pass.

---

## 0. One-line decision

> **Separate hot soil from hot leaf → independent 2×2 tests of passive H₂O conductance and plant xylem viscosity → only with leaf-level evidence add hydraulic stomata (LTO nesting, hydraulic cost only) → thermal-damage / three-objective optimisation last (research branch).**

### Do not do now

| Forbidden | Why |
|-----------|-----|
| Fix the seven ~50–51.5 °C points with leaf stomata | Those are **sunlit soil temperature**, not leaf temperature (low LAI, weak wind, soil absorbing radiation). |
| `gs_total = gs_Medlyn + gmin(T)` into CO₂/`Ci` | Doubles residual with `g0`; cuticle is mainly **H₂O**; Duursma 2019: `g0`, `gmin`, `gcuti` are not interchangeable. |
| Immediate sun/shade 2-D three-objective optimizer | Expensive research path; Sicangco: instantaneous thermal cost does not clearly improve heatwave fits. |
| Edit deleted `ebal_clm45*.m` or fork bulk vs profile strategies | Use shared evaluator on cleaned `ebal` / `ebalProfile` only. |
| Treat prior TH-SOX / Medlyn_TH metrics as current baseline | Re-run same-definition regression on **current** kernels before any claim. |

### Implementation timing

| Now | Later |
|-----|--------|
| This document + Stage 0 analysis scripts / diagnostics only | Stage 1–2 model flags and physics |
| Keep ebal I/O as-is (`assembleEbalInput` → kernel → `applyEbalOutput`) | Shared `plant` stomatal evaluator when Stage 1 starts |

---

## 1. Corrected scientific judgments

### 1.1 Hot soil ≠ hot leaf

- Seven 50–51.5 °C points: **local sunlit soil**, not canopy leaf. Small energy residual + low LAI + weak wind → strong soil absorption, **not** Medlyn high-T stomatal failure.
- Offline SOX shadow and thermal-acclimation summaries for the 2-yr profile: **zero samples** in both `Ta > 30 °C` bins (`docs/figures/R_nature/_profile_2yr/20_sox_offline_shadow_summary.csv`, `21_thermal_acclimation_shadow_summary.csv`). **Current annual outputs cannot test heatwave decoupling.**

### 1.2 Passive conductance: H₂O only

Split:

\[
g_{\mathrm{leaf,H_2O}}
=
g_{\mathrm{stomatal,H_2O}}
+
g_{\mathrm{passive,H_2O}}(T)
\]

| Path | Conductance used |
|------|------------------|
| Photosynthesis / `Ci` | Active stomatal **CO₂** conductance only |
| Dry-leaf transpiration & energy balance | Active stomatal H₂O + **passive** H₂O\((T)\) |
| Wet-leaf evaporation | Existing wet path only; **no** passive add-on |

Re-define or re-calibrate Medlyn `g0` (“minimum/residual stomatal”) so it does not double-count passive cuticular/residual leak (Duursma et al. 2019, doi:10.1111/nph.15395).

### 1.3 Active vs passive decoupling undecided

Marchin et al. 2023 (GCB, doi:10.1111/gcb.16929): heatwave `gs–A` decoupling in woody species may be **cuticular passive leak and/or active stomatal regulation**. Treat as **competing hypotheses**, not a default dual enable.

### 1.4 SOX is hydraulic/drought first, not thermal protection

- Eller et al. 2020 (SOX, doi:10.1111/nph.16419): trade-off between photosynthetic gain and xylem conductance loss → natural replacement for “`sfactor` first on Vcmax → Medlyn lowers gs”.
- Sicangco et al. 2026 (doi:10.1111/nph.71020): hydraulic-cost models reproduce decoupling better than USO; **extra respiratory / instantaneous thermal costs do not clearly improve fits**. “Survival-value thermal cost” remains a **research hypothesis**, not the main package.

### 1.5 Jones et al. 2026 (GCB): structure of \(T_{\mathrm{leaf}}\) nesting

Jones et al. 2026, *Stomatal Decoupling From Photosynthesis Under High Temperatures Is Consistent With Stomatal Optimisation* (GCB 32:e70972, doi:10.1111/gcb.70972):

| Scheme | Meaning | Decoupling |
|--------|---------|-------------|
| **AT** | \(T_{\mathrm{leaf}}=T_a\) | No |
| **LT** | Outer loop: fix \(T_l\) → optimal \(g_s\) → update \(T_l\) (common LSM style) | **No** (EB may close, cooling benefit **outside** objective) |
| **LTO** | Leaf temperature **inside** optimisation (parallel with \(g_s\)) | **Yes** |

COMPASS today (Picard: fixed \(T_{cu}/T_{ch}\) → Medlyn → heatfluxes → update \(T\)) is closer to **LT**. Stage 2 hydraulic stomata must specify **LTO nesting**, not only a new formula with the same outer loop.

### 1.6 Hydraulics topology is not “two free \(\psi\) only”

`solveEbalHydraulicsLayered`: shared root–stem supply, then layer LAI / demand → sun/shade layer \(\psi\) and `wt`. Any future decision variables must keep:

- shared plant water supply constraint;
- layer transpiration allocation;
- layer \(\psi_{\mathrm{leaf}}/wt\) diagnostics;
- sun/shade + layer water-mass closure.

---

## 2. Staged roadmap (authoritative)

### Stage 0 — Mechanism diagnosis (no stomatal equation change)

**Two windows separately:**

1. **Soil-radiation heat:** e.g. KT 5052–5054, 5099–5102 (and analogues).
2. **True leaf heatwave windows:** require layer \(T_{\mathrm{leaf}}\), \(A\), \(g_s\), VPD, LE, \(\psi_{\mathrm{leaf}}\), supply ratio.

**Minimum aligned diagnostics:**

`Tleaf / Tsoil / Ta`, `Rn`, `H`, `LE`, `A`, `gs`, `VPDleaf`, `sfactor`, `supply_ratio`, \(\psi_{\mathrm{leaf}}\), wet fraction, LAI, wind.

**Gate to Stage 1 is not “soil \(T>50\)°C”.** Require leaf-side behaviour such as:

\[
\frac{\partial A}{\partial T_{\mathrm{leaf}}} < 0,
\qquad
\frac{\partial g_{\mathrm{leaf,H_2O}}}{\partial T_{\mathrm{leaf}}} \ge 0
\]

and that the **current** model fails to jointly explain observed leaf temperature, transpiration, and leaf water potential (when observations exist). Always list \(T_{\mathrm{leaf}}\) and \(T_{\mathrm{soil}}\) separately.

---

### Stage 1 — Two independent low-complexity factors (2×2)

| Passive H₂O | Plant xylem viscosity | Role |
|-------------|----------------------|------|
| off | off | Baseline (current) |
| on | off | Test \(g_{\mathrm{passive,H_2O}}(T)\) alone |
| off | on | Test plant-path viscosity alone |
| on | on | Interaction |

**Rules:**

- \(g_{\mathrm{passive}}(T)\) only in **dry-leaf vapour** exchange; not in CO₂ diffusion / `Ci` iteration.
- Revisit physical meaning of `gs_g0` / Medlyn residual vs passive.
- Viscosity factor only on **plant xylem** path; **never** multiply soil conductivity that already has viscosity.
- Do **not** hard-code universal `(T/298)^7`; wrap as `plant_viscosity_factor(T, …)` with documented temperature variable and validity range (Sicangco-style \(k(T)\) is optional form, not law).
- Both flags **default off**; off path bit-identical (or regression-identical) to current.

Illustrative passive form (H₂O only; \(T_p\) = empirical breakpoint from \(g_{\min}\)–\(T\) data, not \(T_{\mathrm{crit}}\)):

```text
g_passive(T) = gmin25 * Q10^((T-25)/10)           if T <= Tp
g_passive(T) = g_passive(Tp) * exp(q*(T-Tp))      if T > Tp
```

---

### Stage 2 — Active hydraulic stomata (only if Stage 1 insufficient)

1. SOX / ProfitMax **shadow diagnostic** first.
2. Then online active stomatal branch.
3. **Hydraulic cost only** in v1; no thermal-damage cost in the objective.
4. Decision variable is **\(g_{\mathrm{stomatal}}\)** (or \(\psi\) control that maps to it). Passive \(g_{\mathrm{passive}}\) is **not** optimised but **must** enter E and whole-plant water accounting.
5. Nesting: **LTO** (Jones): candidate stomatal opening updates \(T_{\mathrm{leaf}}\) **inside** the score, not only outer Picard LT.

**`sfactor` double-penalty trials** (current `biochemical` multiplies same hydraulic `sfactor` onto Vcmax/Jmax/TPU, then Medlyn lowers gs via \(A\)):

| Arm | Meaning |
|-----|---------|
| Medlyn + instantaneous `sfactor` on capacity | Current |
| Hydraulic \(g_s\) + **no** same instantaneous `sfactor` on capacity | Primary online trial |
| Both on | Diagnostic for double penalty only |

Do **not** permanently delete all drought biochemical suppression. Later evidence may add a **slow** damage/acclimation state separate from instantaneous xylem vulnerability.

---

### Stage 3 — Optional research branch

Only if Stage 2 still fails heatwave data:

- cumulative PSII heat damage and recovery timescales;
- survival value / future carbon;
- three-objective C–hydraulic–thermal optimisation.

Needs real apple heat tolerance, co-located \(T_{\mathrm{leaf}}\), \(g_s\), \(A\), \(\psi_{\mathrm{leaf}}\), P50/P88 — **not** a quick fix for current numeric anomalies.

---

## 3. Placement in cleaned energy architecture

```text
buildVeg / site.toml
    → explicit stomata parameter slice (flags default off)
    → assembleEbalInput
    → shared plant stomatal evaluator   ← new when implementing
    → ebal / ebalProfile (same call)
    → out.stomata (+ heat) diagnostics
```

| Fact | Implication |
|------|-------------|
| New chain currently passes mainly `medlyn_g0/g1` (`prepareScopeStep`) | Need full stomata config slice in `buildVeg` + assembly |
| `biochemical` errors on non-Medlyn `gs_model_type` | Type-2 path needs shared evaluator **or** biochemical extension; not a one-line config toggle |
| `runStepEbal` assemble / Anderson retry / apply | **Unchanged** for this work |
| Bulk vs profile | One shared evaluator; no duplicated policy |

---

## 4. Literature anchor (verified DOI list)

| Role | Reference | DOI |
|------|-----------|-----|
| LTO vs LT structure | Jones et al. 2026 GCB | 10.1111/gcb.70972 |
| Hydraulic cost vs thermal cost | Sicangco et al. 2026 New Phytol. | 10.1111/nph.71020 |
| Decoupling, water access | Marchin et al. 2023 GCB | 10.1111/gcb.16929 |
| Passive gmin concepts | Duursma et al. 2019 New Phytol. | 10.1111/nph.15395 |
| SOX LSM | Eller et al. 2020 New Phytol. | 10.1111/nph.16419 |
| USO / Medlyn | Medlyn et al. 2011 GCB | 10.1111/j.1365-2486.2010.02375.x |
| ProfitMax | Sperry et al. 2017 PCE | 10.1111/pce.12852 |
| Heatwave gold data | Drake et al. 2018 GCB | 10.1111/gcb.14037 |

Local copies: user OneDrive plant-hydraulics methods folder; desktop Jones PDF; repo `data/references/` where present.

---

## 5. Superseded content (v2 — do not implement as default)

The following v2 ideas remain **historical / research notes only**:

- Default package “gmin + kmax + three-way carbon–hydraulic–thermal 2-D optimizer” as main path.
- Thermal cost as primary driver of stomatal opening (Sicangco does not support as main fix).
- Treating sun/shade as only two free \(\psi\) without shared supply and layer closure.
- Wiring plans against `ebal_clm45.m` / root golden ebal (deleted; use `vegetation.energy.ebal` / `ebalProfile`).
- Claiming TH-SOX metrics without re-running on cleaned kernels.

Optional long-form equations for ProfitMax / thermal logistic remain in git history of this file (v2) if needed; they are **not** the implementation checklist.

---

## 6. Next concrete actions (when work resumes)

1. **Stage 0 only:** diagnostic notebook/script for soil-heat KT windows vs candidate leaf windows; export the aligned table in §2.
2. If leaf gate fails → stop; do not open Stage 1 flags.
3. If leaf gate passes → Stage 1 config flags + shared evaluator skeleton (passive H₂O + viscosity), defaults off, regression tests.
4. Stage 2 only after Stage 1 failure on co-located \(g_s\)–\(T_{\mathrm{leaf}}\)–\(A\)–\(\psi\).

**Editor note:** Model-code edits under `src/+plant/` require `compass-coding-standard` when implementation starts. This markdown-only update does not.

---

## 7. Change log

| Date | Change |
|------|--------|
| 2026-06-13 | v2: Sicangco-driven demotion of pure thermal-cost opening; gmin competition; three-way prototype language |
| 2026-07-16 | **v3 formal decision:** hot soil vs leaf; H₂O/CO₂ split; 2×2 Stage 1; hydraulic Stage 2 + LTO; three-objective demoted; no code now; ebal migration already done |
