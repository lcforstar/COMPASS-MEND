# 长武冬小麦地上范围收敛与地下率定衔接

Editor: Wenjie Liu; Date: 2026-09-30.

本次按用户最新指令收敛范围：地上保留冬小麦物候与碳氮过程，根、叶、茎、籽粒的分配以用户本地 CoLM 冬小麦代码为参照，凋落物沿用已有方案；其它设置保持 COMPASS。地下水、碳、氮是后续率定重点。用户明确要求撤回本轮新增的站立枯叶库。

## 回退和保留的边界

| 处理 | 内容 | 科学边界 |
| --- | --- | --- |
| 回退 | 本轮实际绿叶氮驱动 LUNA、功能氮投影、季节记忆重置、日长与 PAR 适配等扩展 | 恢复审阅前 COMPASS LUNA 的输入和数值实现；站点原有 `luna` / `cab_from_luna` 开关与光合、光学参数不改变 |
| 回退 | 新增株高/SAI 模式拒绝检查 | 恢复 COMPASS 原消费者；不修改当前站点株高、SAI、LAI 或水力数值 |
| 回退 | APSIM 茎 WSC 分配与向籽粒再转运、LAI 下降转站立枯叶、对应新增库存及生产适配器 | 不保留生产启用开关或通过新增库改变收获与总库存；历史设计及输出保留为旧版本记录 |
| 保留 | 冬小麦播种、出苗、春化、积温、灌浆和成熟/收获分离 | 实测 LAI 零值边沿触发一次收获；成熟后不继续器官生长。原观测不平移、不修改 |
| 保留 | CoLM 冬小麦器官分配及已核参数、灌浆器官氮再转运 | 分配采用 CoLM；COMPASS 的 FlexibleCN、呼吸和实际库存限制仍使用自身合同，不声称完整移植 CoLM BGC |
| 保留 | 冬小麦不沿用多年生树木储氮/木质储备、实际根吸氮、一次收获及剩余碳氮去向 | 保证植物输入与土壤实际扣除一致、收获不重复、凋落物去向不丢失 |
| 保留 | 接受态 AOM/MEND 实际呼吸与生态系统碳诊断 | 用户此前明确要求接入；不改土壤分解动力学，不再用经验 Rh 代替实际所有者输出 |
| 保留 | checkpoint 的 Throughfall/EVAP 因果记忆 | 修复现有地下边界预测器的重启输入，未改水热算法。独立旧证据见 `outputs/wheat_restart_cache_probe_20260930_r2/harvest/history_causality.json` |

地下预热库存适配器 `spinupStateVector` 中专为新增植物库增加的可选 schema 分支也恢复到原版本。这只撤回新增植物库支持，保留原有地下所有者、单位转换及未知库拒绝检查；土壤物理、分解和运输计算未改变。

CoLM 分配来源核对记录见 `outputs/review_aboveground_implementation_20260930/annual_colm_local_evidence/`。外部 CoLM 文件只作参照，不修改。SLA/Cdm 不在本次工作范围。AOM 底物含凋落物、根残体和管理有机肥；活根自养呼吸与根残体异养分解应分别记账。

## 恢复来源、备份及兼容性

选择性恢复来源为 `outputs/changwu_crop_lifecycle_20260930/frozen/`。它保有本轮植被审阅前的地下修改；不使用 Git HEAD 整体覆盖当前工作区，不改冻结副本。

所有被回退的源文件先逐字节备份至 `outputs/vegetation_scope_rollback_20260930/pre_rollback/`。文件清单记录 before/frozen/after SHA256；混合文件按功能块撤回。`stageDailyVegetation` 与 `buildVegState` 的物候候选、日历和实际吸氮重算输入保留。

新生产 checkpoint 签名记录 CoLM 分配/原凋落物范围、编译参数、收获策略、实际根吸氮与碳诊断身份；不允许带新增 WSC/站立枯叶库的旧轨迹静默接入。旧结果及 checkpoint 保留，后续基准使用新目录与冻结版本。

## 验证记录

本次聚焦验证均通过，未放宽原容差：

另有24项独立静态审计通过，核对恢复来源、备份哈希、生产调用、站点/配置范围和395个受保护地下文件。最终证据为 `outputs/vegetation_scope_rollback_20260930/independent_final_scope_audit.json`；唯一地下目录内的回退是上述植物库列表适配器。

| 验证 | 结果 | 证据目录（`outputs/vegetation_scope_rollback_20260930/` 下） |
| --- | --- | --- |
| COMPASS LUNA 与冻结源的数值比较 | 两个隔离 MATLAB 进程；4场景各480小时更新、20日日记录，完整状态与 buffer `isequaln=true` | `luna_reference/frozen_comparison.json` |
| CoLM C/N 生产范围 | 8项：原库 schema、完整茎分配/吸氮需求、成熟停止生长等 | `cn_scope/verification.json` |
| 冬小麦生命周期 | 73项，包括其它PFT历史一致性、灌浆N与一次收获 | `cn_lifecycle/verification.json` |
| 实际碳供给 | 106项，无未融资呼吸债务/收获核销 | `cn_funding/verification.json` |
| 实际根吸氮重算 | 29项，输入冻结、非负库存、守恒及幂等性 | `cn_actual_N/verification.json` |
| AOM/MEND/植物实际碳交换 | 12项，单位、层序、所有者与NEE恒等式 | `actual_carbon/verification.json` |
| 预热库存适配器 | 29项、135库存行；原地下库存与拒绝检查保留 | `spinup_schema/verification.json` |
| 耦合灌浆/收获及真实重启 | 120个接受步；灌浆正同化与实际根吸氮；收获1次；完整Y及8个碳通量字段重启逐位一致；5项坏checkpoint拒绝 | `cn_coupled/verification.json` |

耦合事件窗最大 C 残差为 `1.985e-15 g/cm²`、N 为 `2.140e-16 g/cm²`、DIN 为 `3.683e-12 µg/cm²`；原门限分别为 `1e-10`、`1e-12 g/cm²` 和 `1e-6 µg/cm²`。NEE分量恒等误差为0；实际根N的层厚积分最大误差 `1.271e-21 g/m²/s`。

合成事件窗使用实际长武强迫和生产求解器，并明确设置植物初始库存及阶段以触及灌浆/收获。它验证实现与重启，不能代替长武田间观测验证、全年轨迹或有机分库平衡初始化。历史扩展的47/39/38项验证只属于回退前版本；对应报告和脚本已注明状态。

## 下一阶段

地下观测及参数入口核查见 [地下率定范围报告](../../outputs/underground_calibration_scope_20260930/README.md)。已有 LAI 和 `NM冬小麦实测.xlsx` 继续使用。

1. 固定本次地上方案，另建新版本年度基准，覆盖已有水/N观测日期、灌浆和一次收获，检查完整水与碳氮账本。
2. 先率定现有 COMPASS 土壤 Ks/α；保持冠层、水力、根系和凋落物设置固定。
3. 对未观测有机分库建立预热及实测起点守恒条件化。DOC/DON 的 2000 试算初值不是实测，也不是必须补找的验证数据。
4. 在同一初态方案下，进行地下 DIN 与 MEND 碳参数的敏感性筛选及低维联合率定，避免碳初态改变后沿用旧 DIN 最优参数。
5. 单列水/N留出日期及长窗 SOC、DOC/DON、MBC/MBN 验证；文献不同年份/处理的呼吸用作独立量级和季节型证据，不冒充当前 NM 年观测。

原年度基准只记录到 KT720/8760，无最终结果且本轮未查到活动 MATLAB；停止原因未记录。它未覆盖首个非初态动态观测，更不构成本次限定版本的全年验证。本次不覆盖或自动续接该历史轨迹。
