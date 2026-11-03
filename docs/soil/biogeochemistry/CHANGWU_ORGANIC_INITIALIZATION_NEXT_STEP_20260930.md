# 长武冬小麦：年度新基准之后的未观测有机库预热

Editor: Wenjie Liu; Date: 2026-09-30.

本文件为当前活动源码的只读审阅与实施准备，不是新的预热结果。地上范围按用户决定固定：冬小麦 C/N 分配采用已核实的 CoLM 系数，保留 COMPASS 其余设置、实测 LAI 到零的一次收获和实际根 N 耦合；本文件不引入其他植被物理。未改模型、站点 CSV、观测或历史实验，也未启动或停止 MATLAB。

## 结论

现有代码已经支持“自然循环预热未观测有机库 → 在新 t0 条件化到实测水分、NH4/NO3/SOC → 独立初态制品”的技术链。下一步可以复用 `runChangwuBiogeochemistryPilot` 运行新的最多三个年度自然周期，先获得逐库漂移、完整 C/N 收支和条件化调整量；三周期只是首批预算，当前 pilot 明确不会认证收敛，也不会放行正式初态。

正式率定前仍需补的是**科学验收和诊断驱动**：带单位的库稳定政策、水热与全年通量稳定、不同起始猜测的依赖、条件化后首年响应。无需增加新的植被或土壤过程方程。任意指定的 DOC/DON 是试算起点，不能当作缺失的 2000 年实测数据去寻找，也不能将预热值改名为实测。

## 已支持的边界与证据

| 事项 | 当前活动实现 | 可作出的声明 |
|---|---|---|
| 实测 t0 水分和总 DIN | `src/+init/initialize.m:55–56`；`conditionMeasuredSoilInitialState.m`，新鲜单柱、未冻土、零时间；NH4/NO3 用 FE 容量除干土质量的观测算子回算 | 实测约束在该网格映射上得到重建；外推层不升级为逐层测量 |
| 周期间续接 | `captureBiogeochemistryCycleCarry.m:15–46`、`applyBiogeochemistryCycleCarry.m:11–52` | 携带完整 accepted Y、根水记忆和接受的温度历史；时钟重启，管理重新布置；累积账本不清零，周期通量由端点差计算 |
| 正式 t0 有机库条件化 | `materializeMeasuredBiogeochemistryInitialState.m:59–143` | 固定 fresh IOM，逐层按实测 SOC 缩放独立响应 C 库、对应 N、DOC/DON、微生物相；保持 fresh 水热、DIN、植被原位不变 |
| AOM 去向 | 同上 `:84–86,229–254`；`biogeochemistryBudget.m:35–104` | 未反应 litter/straw/manure AOM 单独带入并计入完整库存；不把它重复加进现行 SOC 观测算子 |
| 初始化调整账 | 同上 `:121–140`；`exportMeasuredBiogeochemistryInitialArtifact.m:61–74` | 显式记录相对供体的初始化 C/N 库存变化、逐层 lambda 和保护状态；不冒充物理输入、呼吸或矿化；正式期累计账本与示踪标签 fresh |
| 有单位的逐库比较 | `spinupStateVector.m:3–18,58–173`、`compareSpinupCycles.m:16–69` | 可按独立库和空间位置比较 g C/N/cm²，检查相态与别名，不让层间正负漂移抵消；库门通过本身不认证周期平衡 |
| 初态制品绑定/拒绝 | `exportMeasuredBiogeochemistryInitialArtifact.m:29–79`、`applyMeasuredBiogeochemistryInitialArtifact.m:21–83` | source、参数、网格、forcing、管理、模型库、观测目标和供体 checkpoint 必须精确匹配；正式导入拒绝未放行制品；仅显式诊断 opt-in 能做短验证 |

新鲜水/DIN 条件化与 SOC 有机库条件化是两个明确阶段。前者不改变有机库；后者不借用实测 DIN 补 DON。SOC 缩放会改变有机 N 总量，当前没有逐层 TN 实测约束可以将这些 N 分库唯一识别；因此它们仍是条件先验，必须报告初始化调整量。

## 数据、SOC 算子与收支保持一致

用户来源澄清在 `data/sites/changwu/PARAMETER_PROVENANCE.md:3–29`：NH4/NO3/SOC 初态有实测依据，DOC/DON 为试算猜测。历史 DON 数值与 2006 年剖面的对应记录保留，但不是 2000 年 DOC/DON 实测初值。

当前 SOC 算子为：POM1 + POM2 + MOC + MB + BD + ENZP1 + ENZP2 + EM + ENZN + DOC_total + IOM + NEC_C。`scoreChangwuCarbon.m:159–184` 与条件化器 `materializeMeasuredBiogeochemistryInitialState.m:200–252` 一致。POC/EP 等只是合计别名，不能重复计入；DOC_total 已包含线性平衡下的溶解与吸附总量，不能再次作为新碳加入。`reserveInitialDOM.m:8–21` 首次初始化分别从 MOC/MON 预留 DOC/DON；应用暖态条件化不再次预留，避免双扣款。

完整 C/N 控制体还包括土柱、地表 AOM 和植物。`biogeochemistryBudget.m:7–13,51–104,114–135` 使用实际累积输入/输出和实体库存；MEND 呼吸、AOM 呼吸、植物呼吸与收获均作为实际输出。根吸 N、植物归还、矿化和尿素水解是内部转移，不能再计一份外部输入。负碳债务及 writeoff 另列，账面闭合不等于物理碳闭合。

后续日期 DOC/DON 的提取方法和相态口径仍需匹配观测算子；这与用户指定 t0 猜测值是不同问题。实测 LAI 文件和已建立的水/DIN/SOC 校准日期继续沿用；不将初始日期本身重复计为动态拟合能力。

## 已完成的技术检查与不能使用的旧预热

当前回退后 `outputs/vegetation_scope_rollback_20260930/spinup_schema/verification.json` 为 PASS，29 项库映射/单位/相态/逐库阈值合同检查，`periodic_equilibrium_certified=false`。地上范围检查及 120 个耦合接受步已记录在同目录的 `cn_rollback_verification.json`，它们验证生产边界和守恒，不代表年度初态收敛。

历史技术报告可用于说明已有接口被实际行使过，但不替代当前模型身份的年度验收：

- `outputs/bgc_pilot_check_20260930/verification.json`：22 项、5 个接受小时，PASS，明确没有年度预热与科学放行。
- `outputs/measured_initial_artifact/20260930_095136/verification.json`：19 项、2 个接受小时，PASS，未放行制品拒绝、篡改拒绝及诊断导入可运行。
- `outputs/measured_initialization_20260929/20260929_203618/verification.json`：29 组代数/保护状态/短窗重启检查，PASS，`scientific_release=false`。
- `outputs/bgc_cycle_carry_check_20260929_204002/verification.json`：10 个接受小时，PASS；包括回放接受态、管理、温度历史和 restart 合同。

旧 `changwu_bgc_pilot_20260930` 已因作物 C/N 合同修复停机。停止记录 `STOPPED_FOR_CROP_CN_REPAIR_20260930.json` 明确：末次落盘进度 KT8490，完成 **0** 个年度周期；最后实际保存的 checkpoint 为 KT7920；无最终初态制品、无收敛放行，禁止在修复机制下续算。旧 progress 的 `running` 字符串是被保留的历史状态。

此外，旧 bridge 年度审计 `outputs/spinup_annual_budget_20260930/20260930_095359/verification.json` 的物理 C 残差 0.0029605862332376942 g/cm²、`numerical_budget_passed=false`，不能拿其账面闭合或旧年末状态作为当前预热通过的证据。新年度基准须独立解决这一验收。

## 新年度基准之后的最小执行顺序

1. **先收新年度基准结果。** 确认全年 8760 接受步、LAI 约定的一次收获、根 N 与根残体归还、xsmr/writeoff 均为零、独立物理 C/N 收支通过及当前水/DIN/SOC 匹配报告。基准使用当前固定地上范围；它属于相同试算初态的基准，不宣称已预热。
2. **冻结新预热身份。** 在新目录复制当前 src、Changwu 原始输入/forcing/管理、模型库、校准工具、pilot driver 与既有 44 行 design。记录全部 hash；sample41 必须仍等于 7 个基准参数。旧实验保持历史只读，首次周期重新由实测 t0 初始化。
3. **行使当前 pilot 的短接口。** 在新证据目录运行 `verifyChangwuBiogeochemistryPilot` 的 5 接受小时技术例，确认回退后的新鲜 t0、两次续接、条件化导入、实际收支和未放行拒绝仍成立。若仅出现身份拒绝，重新创建本身份制品；不放宽拒绝条件。
4. **执行基准起点最多 3 个自然年度周期。** 用同一冻结身份，禁止加速和中途测量重置；每周期终点保留 checkpoint、carry、逐库原始漂移表、水热快照、实际全年通量、作物事件和独立收支。当前 driver 上限为 3，不自动续年。
5. **审阅漂移，再决定预算和政策。** 先列最坏库/层/相及其绝对变化与对称相对变化，分开列快库和慢库。若仍有漂移，报告未评估/未收敛；不能由“跑完三年”认定平衡。
6. **再做起点依赖和条件化后响应。** 既有 2026-09-29 方案提出原猜测、0.5 倍及 2 倍 DOC/DON/微生物快库三种确定性探针，差量在 MOC/MON 内守恒重分配。它们是数值探针，非实测置信区间。逐层供体不足或 C:N 不可行应拒绝，不能裁剪或借 DIN；当前 pilot 没有这种分支接口，需独立诊断驱动。各分支在其自己的冻结身份内回放，比较条件化后的 t0 与校准期预测。
7. **初态协议通过后才接正式率定。** 先水，再 DIN，再低维碳过程；有机初态/碳参数改变后复查 DIN。参数变化须重新建立本参数的暖态，精确身份命中才能直接复用；验证数据不用于选择分配、预热预算或阈值。

第 4 步现成接口示例如下，仅用于实施准备，本次没有执行。冻结根目录必须以 `frozen` 命名，输出目录必须尚不存在，设计表另存于新 campaign；不可把新年度基准 checkpoint 塞给本 pilot：

```matlab
cd('C:/Users/Administrator/Desktop/COMPASS-MEND/outputs/changwu_bgc_colm_cn_20260930/frozen');
addpath(fullfile(pwd,'src'),fullfile(pwd,'tools','calibration'), ...
    fullfile(pwd,'tools','diagnostics'));
result=runChangwuBiogeochemistryPilot( ...
    'C:/Users/Administrator/Desktop/COMPASS-MEND/outputs/changwu_bgc_colm_cn_20260930/pilot', ...
    struct('cycles',3,'design_file', ...
    'C:/Users/Administrator/Desktop/COMPASS-MEND/outputs/changwu_bgc_colm_cn_20260930/design.csv'));
```

## 收敛/放行门：现成检查与待补政策

| 门 | 最小合同 | 当前支持程度 |
|---|---|---|
| 数值与守恒 | 每接受步 DIN 残差 ≤1e-6 µg N/cm²；独立 C ≤1e-10 g/cm²、N ≤1e-12 g/cm²；内部植物→土壤和根 N 转移闭合；xsmr/writeoff 精确零 | pilot 已实现；这些是既定数值闭合门，不是库稳定阈值 |
| 同相位与事件 | 完整365天/小时步，同日同钟端点；播种/出苗/收获各一次，实际施肥/有机管理与回放合同一致；作物终点到收获相位 | pilot 已实现与报告；`periodic_state` 名称仅指相位，不代表平衡 |
| 逐库稳定 | 每位置 `abs(end-start) <= abs_element_g_cm2 + rtol*max(abs(start),abs(end))`；C、N 分设有单位绝对门，报告最坏位置，连续同相位周期通过 | 纯比较器已实现；本身份科学政策及放行评估尚未冻结 |
| 水热与通量稳定 | 同相位水/冰、h、T/U、气压；全年 MEND/AOM/植物呼吸、归还、收获、气态 N 与各边界净量/gross 分量；同时报告内部矿化/固定 | pilot 保存水热端点和外部年度通量但 `stability_assessed=false`；还需明确单位/阈值与诊断评估，内部矿化/固定需另读实际诊断 |
| 起点记忆 | 同一参数、管理下的低/原/高快库分支保持实测 SOC/DIN 和总 C/N 重分配；比较暖态及条件化后校准期响应 | 已有设计，无当前 pilot 生产分支消费者/实施证据 |
| 条件化可用性 | 实测水/DIN 位保持、逐层 SOC 回算、lambda/调整量完整，无负库、无重复 DOM/AOM；显式新 t0，不导入旧累积量 | 条件化器、制品导出/导入已实现；需新年度/分支证据判断首年启动响应 |
| 证据与放行 | 同一 identity + 同一 donor checkpoint SHA 的 `source_audit_verified` 与 `convergence_verified` PASS 证据，完整周期供体 | 制品接口已实现；当前 report-only pilot 不产生 convergence attestation，最终制品仍 NOT_READY |

一般 `runCOMPASSBiogeochemistrySpinup.m:37–45,144–201` 有自动停止与旧 native-unit pool-change 门：默认 rtol=1e-3、absolute_floor=1e-12、最少3周期/连续2次。该函数比较 C/N，但统一的不同 native units 绝对 floor 不能被当作统一物理容差；这些默认也不是长武科学阈值。本次最小方案用 report-only pilot 和有单位比较器，不据一般驱动的 stop 字符串签署科学放行。

可从既有零库存/小库存和短窗算术检查估计数值分辨门，再结合预测敏感性预先冻结科学阈值；不能事后以最坏漂移倒推刚好过门的阈值。是否可报告“有限预热/条件先验”或“周期平衡”，取决于实际证据，两个声明要分别记录。三周期预算结束仍保留 `not_assessed`；正式放行不能靠更改 JSON 的布尔值实现。

## 本次交付范围

本次只增加此审阅文件与 `outputs/changwu_organic_initialization_readiness_20260930/readiness.json`。源码行号、当前/历史技术证据和下一步依赖已核实；没有声称当前新年度基准已经完成，没有生成任何已预热初态。新年度评分由主任务独立产出。
