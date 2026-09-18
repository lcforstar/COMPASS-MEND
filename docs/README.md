# COMPASS 文档目录（按 `src/` 域对齐）

文档按模型包结构存放，便于「改代码时找对应说明」。

**近期主线（雪已收尾 → 优先流 → 迟滞 → 组合）：**  
[`architecture/near_term_physics_roadmap.md`](architecture/near_term_physics_roadmap.md)

权威 skill 路由仍在代码树：`src/compass_skill_routing.md`（`docs/architecture/compass_skill_routing.md` 为副本，以 `src/` 为准）。

---

## 与 `src/` 对应

| 目录 | 对应 `src/` | 内容 |
|------|-------------|------|
| [`architecture/`](architecture/) | 总览 / 编排 | 五域与路线图、全球蓝图、编码规范调研、版本清单 |
| [`step/`](step/) | `runOneStep.m`、`+step` | 主循环契约、状态清单与迁移 |
| [`snow/`](snow/) | `+surface/+snow` | 部分雪盖 / 两 tile 计划 |
| [`surface/`](surface/) | `+surface` | 地表–土壤热耦合、界面温闭合方案 |
| [`soil/`](soil/) | `+soil` | 土壤水热气总物理与审计 |
| [`soil/water/`](soil/water/) | `+soil/+water` | Richards、water-only、W2、边界 |
| [`soil/thermal/`](soil/thermal/) | `+soil/+thermal` | 冻融能量推导等 |
| [`soil/nitrogen/`](soil/nitrogen/) | `+soil/+nitrogen` | C–N 耦合进度与审计 |
| [`plant/`](plant/) | `+vegetation` | 光合、水力、物候、四象限、叶水势 |
| [`radiation/`](radiation/) | `+surface/+radiation` | VOD / 微波路线 |
| [`solver/`](solver/) | `+solver` | 内外层时间步等 |
| [`init/`](init/) | `+init` | 初始化、workspace、handoff |
| [`io/`](io/) | `+io` | NetCDF 输出设计 |
| [`diagnostics/`](diagnostics/) | `+diagnostics` | 输出变量目录 |
| [`groundwater/`](groundwater/) | `+groundwater` | （预留） |
| [`equations/`](equations/) | `+equations` | （预留） |
| [`deep_learning/`](deep_learning/) | `+deep_learning` | PFEM 暖启动等 |
| [`model_comparison/`](model_comparison/) | — | 与外模型对比 |
| [`figures/`](figures/) | — | 图件（非 md 过程文档） |
| [`archive/`](archive/) | — | 历史计划、Codex 原文备份、superpowers |

---

## 优先阅读（当前工作）

1. [`architecture/compass_complete_architecture.md`](architecture/compass_complete_architecture.md) — **五区架构、迁移台账、Julia 参考窗口与全球扩展的唯一权威总纲**  
2. [`soil/compass_soil_water_heat_physics.md`](soil/compass_soil_water_heat_physics.md) — **土壤水–热–气–冻融权威书**（Grok/代理**每次新对话必读**；雪与热边界 §2.8.1、§5.7）  
3. [`architecture/grok_compass_expert_context.md`](architecture/grok_compass_expert_context.md) — **Grok 专家情境**（含必读规则，可贴系统提示）  
4. [`surface/2026-07-20_atmospheric_matched_flux_guide.md`](surface/2026-07-20_atmospheric_matched_flux_guide.md) — 大气匹配通量短文（原理/公式/代码表）  
5. [`architecture/near_term_physics_roadmap.md`](architecture/near_term_physics_roadmap.md) — **雪收尾 → 优先流 → 迟滞 → 组合** 的分块架构与验收顺序  
6. [`surface/surface_soil_thermal_coupling_progress.md`](surface/surface_soil_thermal_coupling_progress.md) — 热耦合进度与全年验收  

---

## Codex 记忆原文

已从 `C:\Users\zjzx\.codex\memories\extensions\ad_hoc\notes\` 复制到  
[`archive/codex_memory/`](archive/codex_memory/)。  
以仓库内整理版路线图为准；原文仅作溯源。

---

## 测试与工具（同样按域）

| 树 | 说明 |
|----|------|
| [`tests/`](../tests/README.md) | 回归/单元，目录与上表同构；先 `run('tests/setupPath.m')` |
| [`tools/`](../tools/README.md) | 诊断/标定/基准；先 `run('tools/setupPath.m')`（MATLAB） |

## 整理说明

- 2026-07-19：根目录散落 `*.md` 按上表迁入；`superpowers/` → `archive/superpowers/`；`notes/` → `archive/notes/`。  
- 2026-07-19：`tests/`、`tools/` 按同一域分类；旧 `tests/test_xxx.m` 扁平路径可能失效。  
- 非 md 的 pptx/py/pdf 进 `archive/presentations/` 或仍在 `figures/`。  
- 旧路径链接可能失效；请用本 README 或按包名搜索。
