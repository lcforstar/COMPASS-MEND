# 长武肥料来源氮 Morris 敏感性分析方案

日期：2026-09-21

## 1. 当前结论

肥料来源氮示踪已经具备开展全局筛查的数值基础。本方案采用 Morris
elementary-effects 方法，先回答“哪些过程控制肥料氮进入有机库、从有机库释放及其后续
去向”，再决定哪些参数进入率定或 Sobol 二阶段分析。

本阶段是**敏感性筛查，不是率定**。参数范围分为 MEND 文献先验、PFLOTRAN
源码示例、长武氯离子弥散代理和无直接观测支持的工程包络。后两类不能解释为
长武可信区间。

Morris 轨迹采用偶数层级 \(p\)，归一化步长

\[
\Delta=\frac{p}{2(p-1)}
\]

每条轨迹恰好逐个改变每个参数一次。对输出 \(y\) 的 elementary effect 为

\[
EE_i=\frac{y(\mathbf{x}+\Delta_i)-y(\mathbf{x})}{\Delta_i}.
\]

报告 \(\mu\)、\(\mu^\*\)、\(\sigma\) 和有效轨迹数；其中 \(\mu^\*\) 表示总体影响
强度，\(\mu\) 的符号表示平均方向，较大的 \(\sigma\) 表示非线性或参数交互，但不能
单独区分二者。方法依据为
[Morris (1991)](https://doi.org/10.1080/00401706.1991.10484804)。

## 2. 可执行参数集

生产年默认包含 30 个自由度，完整机器可读规范由
`tools/sensitivity/changwuFertilizerMorrisSpec.m` 从当前站点输入实时生成。
基线不从旧文档抄录，以避免以下已确认漂移：

- 当前 `active_fraction_initial=0.01`，不是旧审计中的 0.5；
- 当前 `doc_leach_1d=0.01 d^-1`，不是旧表中的 0.0001；
- 当前 DOC Kd 表层样点为 0.5、其余为 3.98 cm³ g⁻¹，DON Kd 为 3.98。

| 参数组 | 自由度数 | 参数化原则 |
|---|---:|---|
| 无机氮转化/输运 | 5 | 尿素水解、硝化容量、反硝化容量、NH₄ 吸附用剖面乘子；DIN 弥散用全剖面标量 |
| MEND C/N 动力学 | 16 | 源小时独立配置；成对 Ks 共同缩放；派生的 \(V_m\)、死亡率和休眠维持率不重复抽样 |
| DOM 输运 | 2 | DOC/DON Kd 成组缩放；DOM 弥散为全剖面标量 |
| 微生物运移 | 3 | 附着/脱附共同缩放，截留保持正值，微生物弥散为全剖面标量 |
| 水力背景 | 2 | Ks 与 VG-\(\alpha\) 仅用全剖面乘子，保持剖面形状 |
| 有机肥背景 | 1 | `k_manure_1d` 仅在保留有机肥管理事件的生产年情景中使用 |
| AOM 周转 | 1 | `doc_leach_1d` 同时作用于有机肥和自然 PlantCN 凋落 AOM，两个情景均保留 |

主要 MEND 边界来自 Wang et al. (2020) Table 4
([Journal of Hydrology 585, 124777](https://doi.org/10.1016/j.jhydrol.2020.124777))。
以下范围仍明确标为未验证工程包络：无机氮速率、NH₄ Kd、DOM Kd、微生物截留、
水力剖面乘子、肥料背景分解，以及存在体积基/干土质量基歧义的 MEND Ks。

所有土层参数只释放一个剖面乘子；原有结构零值保持为零。首轮不逐层释放参数，
也不把同一方程中的 Vmax 与 Km 同时全部释放。

## 3. 明确排除项

- `reaction_substep_s`、DOM/BIOPARTICLE transport substep 等数值控制只用于
  收敛验证，不参与率定或 Morris 排名。
- 当前 `source_hourly` 要求微生物 C:N 下/上界固定为 2/14，酶 C:N 固定为 3；
  `mend_cn_mb_min/max`、`mend_cn_enzp/enzm` 不能作为生产自由参数。
- `linear_equilibrium` 下 MEND kinetic sorption 的 Qmax/Kba/Kdes 被显式关闭，
  不进入筛查。
- `active_fraction_initial`、`mend_lci_initial` 和
  `mobile_fraction_initial` 是初始化先验，只有在一致预热的独立设计中研究。
- `mend_input_fraction` 和初始 SOC 四库比例是 simplex，不能逐项无约束采样。
- 尿素 Q10/WFPS 系数和 NH₃ 挥发特征长度目前硬编码，尚无配置接口，本轮不能抽样。

## 4. 两个科学情景

1. **production_year（主分析）**：保持当前 NM 管理，2000-06-15 10:00 的矿质尿素
   事件和 11:00 的有机肥背景均存在。该情景回答真实当前配置下的净过程控制。
2. **mineral_pulse_only（机制隔离）**：只删除有机肥管理事件，保留自然 PlantCN
   凋落和唯一矿质肥脉冲。该情景只排除失活的 manure 分解参数；`doc_leach_1d`
   仍控制自然凋落来源 AOM 的 DOC/DON 淋溶，因此继续参与筛查。该情景用于区分
   有机肥背景对肥料来源氮固定/释放的间接影响。

两个情景必须分别排名，不能把管理情景差异混入同一 Morris 参数效应。

## 5. 输出与解释

每个样本至少保存以下接受态指标：

- 肥料来源氮输入；
- 毛有机进入和毛有机释放；
- 最终 NH₄、NO₃、尿素、有机库和总示踪储量；
- 最终 DON、PON1、PON2、颗粒有机氮合计、矿物结合有机氮（MON）、微生物氮和酶氮分库；
- 根系 NH₄/NO₃ 吸收；
- NH₃、硝化 N₂O、反硝化 N₂O/N₂；
- NH₄/NO₃ 顶部与底部分项边界损失、上下边界合计、数值修正和累积闭合残差；
- 上述关键量相对于肥料输入的比例。

毛进入/释放是可重复循环的内部通量，可能超过一次输入，不能解释为互补的最终去向。
最终有机储量、矿质储量、根吸收、气体损失和边界损失才用于终端质量分配。

24 h 仅用于框架验收；此前该窗口的释放量极小，不能据此形成科学排名。正式筛查使用
完整 8760 h 生产年。若研究多年有机库再释放，应在首年筛查后延长至多年，并保持同一
参数集、输入来源定义和预热合同。

## 6. 默认规模与运行方式

生产年默认 \(k=30\)、\(p=6\)、\(r=10\)，共
\(r(k+1)=310\) 个年度样本。先生成设计并审计：

~~~matlab
addpath('src');
addpath('tools/sensitivity');
opt = struct('scenario',"production_year", ...
    'n_trajectories',10,'n_levels',6,'seed',20260921, ...
    'n_steps',8760,'design_only',true);
runChangwuFertilizerMorris( ...
    fullfile('outputs','changwu_fertilizer_morris_year1'),opt);
~~~

审计 `parameter_spec.csv` 和 `design.csv` 后，将
`design_only=false` 运行。首版串行执行；MATLAB 进程环境变量是全局状态，不能在
同一进程内用 `parfor` 修改输出/续跑变量。批次级续跑会严格比较站点文件 SHA-256、
源码 SHA-256、`data/model_library_nc` 静态物理库 SHA-256、冻结的物理环境变量、设计、
情景和参数向量；每个实际执行样本在开始前和接受结果前再次核对三类输入哈希及物理环境。
已通过样本跳过，失败或中断样本从 KT1 在新的 attempt 目录重跑；身份漂移直接终止 campaign。

完整年度筛查结束后，按多个机制输出的 \(\mu^\*\) 联合选出约 8–12 个参数进入 Sobol
分析。选择时必须同时满足质量闭合、非负、无步失败和独立率定数据可约束性，不能只按
某一个毛通量的单指标排名。只要某个输出存在任一参数的缺失轨迹，该输出仍保留
`mu/mu_star/sigma/n_valid` 供故障诊断，但整个输出的 `rank_mu_star` 均置为 `NaN`；
完全相同的 `mu_star` 使用并列平均名次。

## 7. 验收入口

~~~matlab
addpath('src');
addpath('tools/sensitivity');
verifyChangwuFertilizerMorrisFramework( ...
    fullfile('outputs','changwu_fertilizer_morris_smoke_20260921'));
~~~

验收包括固定种子复现、每步单参数变化、linear/log 映射往返、线性合成函数的
\(\mu/\mu^\*/\sigma\)、30 参数基线完整重建与生产初始化逐字段一致、隔离情景保留
`doc_leach_1d` 并只排除 `k_manure_1d`、调用者物理环境精确恢复、两参数一条轨迹的
24 h 三样本生产试验，以及第二次运行精确跳过已完成样本。

## 8. 2026-09-21 工程验收记录

- MATLAB R2025a 静态检查：9 个新增/受影响 MATLAB 文件均为 0 条 Code Analyzer 告警。
- 最终框架验收：30 参数生产年规范、29 参数矿质脉冲隔离规范、固定种子设计、linear/log
  映射、合成函数 \(\mu/\mu^\*/\sigma\)、不完整排名抑制、并列名次、基线完整重建、
  调用者物理环境恢复、3/3 生产烟雾样本和精确续跑跳过全部通过；最大柱示踪闭合残差为
  \(2.64\times10^{-19}\ \mathrm{g\ N\ cm^{-2}}\)。
- 最终源码身份下的全参数边界轨迹：\(r=1,p=6,k=30\) 的 31/31 个 24 h 样本通过；
  最大柱示踪闭合残差为 \(4.6354\times10^{-19}\ \mathrm{g\ N\ cm^{-2}}\)。同目录精确续跑时
  attempt 数保持 31→31。该结果仅证明参数包络在短窗内可执行，不构成参数重要性排序。
- 示踪 A/B：24 h tracer-on/off 的所有非示踪总量状态与运行态逐字段位级一致；事件原子性、
  checkpoint 重启和来源示踪闭合均通过。
- 边界轨迹曾暴露 MBA 直接释放重构的 1 ULP 负余量。修复仅把不超过 8 ULP 的负 surplus
  规范为零并保持接受态总释放不变；64 ULP 负控仍失败，因此没有放宽真实负通量。
- 样本身份现绑定并逐样本复核源码、站点、实际 `model_library_nc`、六个物理环境变量及
  MATLAB release/architecture；通过/失败结果、sample manifest 和 42 项指标 schema 在汇总时
  再次校验。PON1/PON2 与 MON 分开报告，NH₄/NO₃ 的上下边界损失也分项报告。

验收产物：

- `outputs/changwu_fertilizer_morris_smoke_frozen_final_20260921`
- `outputs/changwu_fertilizer_morris_pilot_r1_24h_frozen_final_20260921`
- `outputs/fertilizer_tracer_post_ulp_regression_20260921`
- `outputs/changwu_fertilizer_morris_year1_final_20260921`

正式 310 样本设计已经冻结但**尚未执行全年样本**，campaign 签名为
`3aaa1b5926fa2ee9cd749ee0a36266d70d795aee16edb6359782cf535325d44f`。冻结身份为：

- 源码 SHA-256：`0a7aa6de7ac7b2732bd37dfeaa81d48f8189d965d1dd3884102e7cc82e63cc3e`；
- 长武站点 SHA-256：`9ee407a12181f6efc3c7340886838af17b0adeb4e57b48f010eadd57d5cad1b8`；
- 静态物理库 SHA-256：`1cf207a8d46ce1210eaca195ff6f7331140e661bdcb78f22e9bee298408f9184`；
- MATLAB `R2025a`，`win64`。

该目录、三组输入哈希、运行环境合同和 campaign 签名共同构成后续全年运行的唯一续跑身份。
旧 `outputs/changwu_fertilizer_morris_year1_frozen_20260921` 的 `28348…` 签名缺少新增身份字段，
仅保留为预冻结历史记录，不得继续运行或引用为正式设计。24 h 的 `morris_summary.csv` 不得用于
科学解释；只有 310 个 8760 h 样本完整通过后才可形成正式 Morris 排名。
