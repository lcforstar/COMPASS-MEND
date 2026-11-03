# COMPASS skill routing

Single task-to-skill map for this repository. Updated 2026-09-15.
Instruction priority, authorization, reading scope and completion follow [AGENTS.md](../AGENTS.md).

## 1. Select by deliverable

Use the smallest applicable skill set. Announce a skill and its purpose once when first used, read
its body with an available file-reading tool, and load supporting references only for the selected mode.
Audit targets are material being reviewed; reading a Skill for an audit does not invoke its workflow.
Do not select a specialist workflow from an isolated keyword such as paper, citation or diagram.

## 2. Resolve paths

COMPASS-owned canonical bodies:

| Skill | Repository-relative path |
|---|---|
| compass-coding-standard | `.claude/skills/compass-coding-standard/SKILL.md` |
| compass-julia-port | `.claude/skills/compass-julia-port/SKILL.md` |

Resolve these against the current repository root, independent of username or drive. An optional
Grok/Claude/Codex pointer is an entrypoint, not a competing authority. For other skills, use the current
host's registered name and exact path. Do not search other machines' directories or install a library
merely because a historical skill name is missing. Use an equivalent available workflow or ordinary tools
when sufficient; a missing irreplaceable scientific contract blocks only work that depends on it.

## 3. Task routes

| Requested work | Applicable workflow |
|---|---|
| Change/review COMPASS model semantics: physics, I/O contracts, state, boundaries, solver, diagnostics | `compass-coding-standard` |
| Design, implement or review COMPASS Julia migration | coding-standard, then `compass-julia-port` |
| Local bug or routine software feature | Focused diagnosis/design; add coding-standard only for model semantics |
| Scientific hypotheses, mechanisms or candidate research directions | `scientific-brainstorming` |
| Quantitative scientific figure | Available scientific visualization skill; `nature-figure` for publication figures, using the requested journal's actual requirements |
| Simple flowchart or model diagram | One suitable diagram skill; use native vector geometry when sufficient |
| Complex illustrated/editable scientific reconstruction | `scientific-figure-workflow` or `nature-image2ppt`, according to the source and deliverable |
| Model-observation statistics | Appropriate analysis tools; `nature-statistics` only for manuscript statistical reporting |
| Literature discovery or claim support | `nature-academic-search` for coordinated literature work; ordinary primary-source search for a small lookup |
| Explicit Nature/CNS-only citation search | `nature-citation`; do not impose its journal whitelist on general references |
| Local Zotero lookup/import/export | `zotero:Zotero` |
| Source-grounded paper question or summary | Read relevant source sections and answer the requested question |
| Full-paper bilingual reader/translation artifact | `nature-reader` |
| Manuscript drafting or polishing | `nature-writing` / `nature-polishing` when their actual mode fits; journal-specific rules only for the target journal |
| Proposal composition/revision | `researchwrite`, scaled to the requested section or project |
| Academic Word document with citation fields | `academic-manuscript` plus the available documents skill |
| General Word/PDF/spreadsheet/presentation artifact | Current host's documents / pdf / spreadsheets / presentations skill |
| Paper-based research slide deck | `nature-paper2ppt` |
| Data Availability, peer review or reviewer response | `nature-data`, `nature-reviewer` or `nature-response` respectively |
| Create or revise Skill instructions | `skill-creator`; no model-code workflow for instruction-only changes |

Use a visualization skill for user-facing scientific figures when available. If none is available,
use the relevant plotting tools while retaining data integrity, units, legibility and export verification.
Pure formatting or non-model Python post-processing may skip coding-standard; no line-count exemption
applies to changes in equations, units, signs, state ownership, boundary or acceptance semantics.

For literature-backed scientific claims, prefer relevant local material in `data/references/` when
available, then verify primary sources as needed. Do not fabricate sources or infer implementation
from a reference model. This does not delay an explicit request for current online documentation.

## 4. Optional full research project route

For an end-to-end research project, preparation, scientific design, execution/analysis, writing and
presentation form a useful sequence. Execute only stages within the current request. A local edit,
figure or review does not require the entire research pipeline, a manuscript or a slide deck.
