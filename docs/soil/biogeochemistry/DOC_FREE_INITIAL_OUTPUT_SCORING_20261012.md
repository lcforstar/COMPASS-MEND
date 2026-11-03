# DOC 初值、输出、评分修正：doc_free_d_v2

Editor: Wenjie Liu  
Date: 2026-10-12

## 定义与范围

用户确认的 DOC 初值是独立溶液态库 D，不含 QOC；MAOC = MOC + QOC。
2006-09-15 观测仍是暂代 2000-06-15 的初态代理，不伪称同期动态观测。
本次更正接口、诊断和选优合同，不修改吸附/反应方程，不启动全年或六年模拟。
历史冻结代码、输入、结果及评分全部保留。

内部 `soil_dom.doc_total_g_cm3` 仍是 D+Q 的唯一状态所有者。这一数值实现不等于
将观测 DOC 定义成 D+Q，也不要求为 D/Q 增加重复的可写状态。

## 初值

`applyChangwuProxyDomInitialState` 将代理 DOC 插值至 0–300 cm 的计算层：

\[
D_0 = DOC_{obs}\,\rho_b\,10^{-6},\qquad
Q_0 = D_0\rho_bK_d/\theta_l,\qquad T_0=D_0+Q_0.
\]

观测单位 mg C/kg 干土；D、Q、T 单位 g C/cm³ 土壤；容重 g/cm³、Kd cm³/g。
含水量读取已初始化接受态的液态水缓存，不使用含冰总水量或任意水分下限。
正 D 与零液态水、非有限值、负值或不足的初始 SOC 预算均拒绝。
300 cm 以下原总库保持不变，不把其旧总量再次解释成 D。

为避免夹带新的科学参数化，本次保留原有的初始分配规则：从响应有机碳库按原比例
预留 D+Q，配对有机氮同比例调整，DON 继续使用原先 **DON+QON 总量**输入口径。
逐层总 SOC 保持不变；这不是物理时间步或外部碳通量。
该既有规则**不保证总氮不变**，新增逐层有机 N 调整及列积分审计，明确
`total_N_conservation_claimed=false`、`physical_flux_added=false`。
水、热、无机氮、植被、IOM及未反应外加有机物状态不改。

重要：在现有 Kd 下，真实初态测试的最低响应库保留比例约 0.4183，
有机 N 初始化调整为 −0.0582492 g N/cm²。它们是原比例分配规则的结果，
不是验证了本站合理的初始微生物量或 QOC。正式新筛选前仍需审查 D/Q/MOC 分配，
不能据“接口测试通过”宣布初始生态状态已得到观测验证。

通用 `soil_initial_profile.csv` 新增明确字段 `doc_free_mgkg`，含义为 D；
原 `doc_mgkg` 为历史 D+Q 接口，仅为兼容保留；`doc_mgL` 仍是水相浓度。
三者互斥，禁止相同数值混用单位。自由 D 输入预留的是反算后的 D+Q 总碳。
通用冷启动的自由 D 字段目前只支持非冻结初态；启用冻结且温度非正时明确拒绝，
不悄悄声称初始相分配正确。当前长武代理路径在液态含水量重建后执行，不受此限制。

## 输出与评分

| 输出 | 含义 |
|---|---|
| `doc_carbon_profile` | D，g C/cm³ 土壤，不含 QOC |
| `doc_aqueous_concentration_profile` | 水相浓度，mg C/L 水，与 mg/kg 干土不同 |
| `qoc_carbon_profile` | Q |
| `moc_carbon_profile` | MOC，不含 Q |
| `maoc_carbon_profile` | MOC+Q，仅诊断汇总，不重复加进 SOC |
| `doc_plus_qoc_carbon_profile` | D+Q，仅明确标名的总量诊断，不当 DOC 评分 |

`diagnostics.docFreeLayerMean` 使用目标深度区间与计算层的真实重叠厚度，按干土质量
加权 D。0–20 cm 年末 ≥100 mg/kg 门槛原数值暂不改变，但判据改为 D，
它仍是用户指定的可行性约束，不是伪造的 2001 年观测。

设计、运行规格、初始化审计、结果和评分均携带 `doc_contract_version=doc_free_d_v2`。
仅接收 `proxy_doc_free_initial_0_20_mg_kg`、`proxy_doc_free_final_0_20_mg_kg`、
`proxy_doc_free_persistence_pass`；无版本或旧 D+Q 字段不能回退补位。
高 Q 不能掩盖低 D；数值达标与通过标志必须同时成立，原水氮/SOC损失及物理检查保留。
无机氮筛选只能继承已完成并通过新版 D 约束的实际碳候选，不再硬编码历史碳参数。

## 参数范围

详见 [边界与来源清单](../../../tools/remote/doc_proxy_calibration/BOUNDS_AND_DOC_CONTRACT.md)。
当前可搜索三轴：Vd 0.1–100 h⁻¹、Vg 0.001–0.1 h⁻¹、KsDOM 0.05–0.5 mg C/g 干土。
KsDOM 的质量基准来自2015论文；不直接搬用2020论文的体积基准数字。
8+2候选中9项有论文先验；ratio 的2–20只查到作者代码依据，所以固定使用，禁止独立搜索。
吸附 Kd 不偷用动力学吸附 Kba/Qmax/Kdes 的范围；没有适配依据的轴不释放。
这些是来源明确的搜索先验，不是本站后验或普适物理极限。

在候选生成、站点写入、启动前及评分处检查版本/上下限；篡改、NaN/Inf、未知轴和
活动轴的土层覆盖被拒绝。边界/轮数耗尽报告未达标，不扩界、不改门槛。
新参考 KsDOM=0.05 与历史0.01不同，设计记录这项变化，不冒充同一参考组。

## 验证与部署边界

`verifyDocFreeInitialContract`：代数反算、旧接口兼容、互斥/冻结/预算拒绝、层序、
干土质量加权、真实预热初态 D 回算、SOC、保护态、N调整及幂等性；零模型步。
真实剖面最大 D 回算误差 1.14e-13 mg/kg，SOC误差 5.21e-18 g/cm³。

`verifyDocFreeCandidateInitialization`：真实站点副本中新参考 KsDOM=0.05 与候选0.1，
Vg减半；二者初始化状态位级一致，只有声明的基础/派生参数改变，共用初态文件及原始
观测/输入哈希不变。直接MATLAB入口的同版本边界篡改也在初始化前拒绝；零模型步。

`verifyDOMAqueousDiagnostics`：线性平衡和动力学吸附制造态、非零 MOC、D/Q/MAOC
恒等式、无重复 SOC、拒绝步不输出及 NetCDF读回；零模型步。

Python检查覆盖已下载的真实冻结运行模板 → 新双轴 → 三轴/N11 的补丁链、合同转发、
旧口径拒绝、非法边界、文件篡改、源码依赖清单及有限搜索耗尽。
所有测试使用副本，不覆写旧运行。测试通过不等于全年模型通过。

未来准备新运行必须携带 `SOURCE_OVERLAY` 源文件和边界模块、重算源码哈希，
验证 `require_v2_runtime`；仅改 design 版本号不能使旧冻结实现成为新版。
服务器补丁包是待部署代码，不更改 current、不终止进程、不自动开始新筛选。
