# MEND耦合接口的代码查阅补充

本页用于定位函数和核对实现细节。接口规则请先阅读[中文耦合说明](MEND_SOURCE_PRODUCTION_CONTRACT.md)，测试结果见[验收表](../../model_comparison/MEND_ACCEPTANCE.md)。以下映射依据2026-09-15核对的正式代码。

当前长武选择 `microbe.reaction_method="source_hourly"` 和 `mend_ihr=0`。氧重试最多24次；额外供氧比例受限时乘0.8裕量，累计因子小于1e-12时置零。这些是数值策略，不是生态速率参数。碳状态从浅到深、单位mg/g干土；有机氮与DOM库存从深到浅、单位g/cm³土体。旧QOC/QON兼容字段保持零，胞内酶C:N=3。

已核对文件哈希见[代码证据记录](../../model_comparison/mend_interface_contract_evidence.json)。这份记录不是最终发布证书。

## 函数与实现位置

| Item | Active source evidence and implementation contract |
|---|---|
| Entry | `runOneStep -> soil.nitrogen.advanceAcceptedCandidate`; organic candidate construction precedes the retained inorganic N solve and atomic commit. |
| Selection | `microbe.reaction_method = source_hourly` selects the paired hourly reaction. Explicit `dom.sorption_model = linear_equilibrium` selects the current Changwu Kd closure; `mend_dynamic` retains the historical source sorption experiment. Source-hour enzyme states are retained in either closure. |
| State | `soil_dom`: bottom-up total and sorbed C/N, g/cm3. `soil_microbes`: top-down solid/microbial/enzyme C, mg/g. `soil_organic`: bottom-up paired N, g/cm3. QOC/QON compatibility stores remain zero. Intracellular ENZN is independently accounted, CN=3. |
| Initialization | Preserve total DOC/DON; linear equilibrium stores totals only and derives phases at every read. Dynamic mode partitions once using initial Kd then stores explicit sorbed mass. Zero ENZN adds no C/N. Mode schemas reject incompatible checkpoints; equilibrium signatures also bind both Kd profiles. |
| Parameters | `mendEnvironmentalParameters` retains existing source responses; source MOD_MEND.F90:2008-2027 and 2329-2332 give Qmax unchanged, Kads=Kba*Kdes_ref with Ea=5 kJ/mol, Kdes with Ea=20. Aligned namelist final Pinitial: Qmax=1.5, Kba=6, Kdes=.006, fpENZN=1. |
| Root/DIN | Root-first shadow inventory is reserved once. Free NH4/NO3 drives first kinetics; total post-root NH4 and liquid NO3 bound actual uptake. Source plant denominator terms are omitted because COMPASS root allocation owns plant competition. First uptake is capped to each available species; subsequent release/second uptake follows source order. |
| Inorganic boundary | COMPASS owns inorganic chemistry, NH4 equilibrium sorption, gas and transport. It has no NO2 state; enzyme allocation explicitly uses resolved free NH4+NO3 only. This interface difference is not full MEND inorganic equivalence. |
| Oxygen | Propose the complete paired source hour, then apply one layer factor to *all* its C/N inventory increments and actual exchange amounts, including overflow/enzyme costs. This is a conservative interpolation of a discrete map, not a rate-rescaled source hour. O2 retries rebuild from identical accepted owners; factor zero means zero reaction increment. |
| Microbial phases | Signed source dormancy compensation is carried as the corresponding reverse transfer separately for C and N. Enzyme synthesis follows active-biomass phase shares; intracellular enzyme turnover removes the same fraction from each enzyme phase and returns to MBA. ENZN phases follow the same BIOPARTICLE transport/exchange operator as cell biomass, preserving C/N and CN=3. This is an explicit COMPASS transport closure absent from the source one-box model. |
| Ordering | Management/litter/AOM candidate additions, free DOM transport, paired hourly reaction, microbial transport, inorganic chemistry/transport; one accepted commit. MEND external plant/soil input functions are not additionally applied. |
| Clock | Source discrete map accepts exactly 3600 s. Other outer durations fail before commit; soil internal substeps do not call the hourly map repeatedly. |
| Restart and output | Whole-owner commits/checkpoints preserve added fields; topology signatures reject incompatible modes. Diagnostics use the DOM phase reader and include ENZN; periodic convergence includes free/sorbed DOM and enzyme phases. |
| Acceptance | Source Fortran reset/continuous parity unchanged when coupling options absent; manufactured free-vs-total DIN and species caps; zero/partial-O2 paired closure; phase transport conservation; actual runOneStep, checkpoint continuation, failure-state immutability, short-window production run before background annual launch. |

The biological core is source-identical in its tested subset. Root priority,
resolved mineral species, oxygen interpolation, phase transport and ordering
are explicitly coupled-model differences; no whole-model identity claim follows.
