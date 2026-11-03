# DOC/DON 底部排水边界：核查与建议

Editor: Wenjie Liu; Date: 2026-09-29.

状态：以下保留的是实施前核查与设计快照。随后经用户放行，已在五组试验完成后
实现并启用长武 DOC/DON 出流边界；实现范围及验收见
[实施记录](DOM_BOTTOM_BOUNDARY_IMPLEMENTATION_20260929.md)。旧 campaign/检查点未改写。
采用 `compass-coding-standard` 的边界、接受子步、状态所有权和守恒要求；物理背景
复核 §2.8.1、§5.7。下文“当前”“建议”等词以实施前状态为准。

## 1. 结论与当前事实

建议长武田间自由排水场景让 DOC/DON 的水相随底部液水排出；边界采用零弥散梯度。
这是与无机氮一致的物理原则，但不能原样复制无机氮有限元边界行。
MEND 在此提供局地 C/N 反应核；参考模型不包含淋溶，不是 COMPASS 必须封住 DOM
底边界的理由。边界由宿主的水文场景确定，局地反应守恒仍应保持。

| 对象 | 当前可达代码/配置 | 已核实行为 |
|---|---|---|
| 水 | `data/sites/changwu/site.toml`，`water_bottom_kind="free_drainage"` | 水允许底部排出 |
| NH4/NO3 | 同文件底部代码均为 3；`soil.nitrogen.transport.applyBoundaryConditions` | 底行施加 C1=C2；不是总溶质零通量 |
| 无机氮底通量 | `soil.nitrogen.transport.calculateBoundaryFlux` | 类型 3 从保留的原始边界方程恢复实际通量，不简单等于 q×C |
| DOC/DON | `soil.dom.transportConservative` | 仅 n−1 个内部面，没有上下端面传质；末尾明确 closed_zero_mass_flux |
| DOM 水驱动 | `soil.dom.advanceAcceptedCandidate`，`solver.recordAcceptedTransportTile` | 接受子步 dt、单元中心 QL 和配对液水 θ；缺独立底面液水通量 |
| 肥料 DON 标签 | `soil.nitrogen.tracer.replayTransportLedger` | 目前假定闭端面和标签柱质量不变 |
| 微生物颗粒 | `soil.microbe.bioparticleTransport` | 共用 DOM 输运核，但有独立闭合柱合同；不能随 DOM 改动而默认开放 |

当前没有把 DOC/DON 底边界改开放的有效配置开关。只在 TOML 添加一个新名字并不会
使其生效，必须连通配置读取、输运、标签和账本。

## 2. 数学合同

物种 X 分别为 DOC-C 和 DON-N。以 q_b>0 表示液水向下离开土柱，
C_X,b 为底层**水溶液浓度**，则出流分支：

    F_X,out = q_b * C_X,b

其中 q [cm/s]、C [g 元素/cm³ 水]、F [g 元素/cm²/s]；
质量出口为各接受子步的 F*dt 之和，不用小时平均 q 乘小时末 C 替代。
边界弥散/扩散项为零：零浓度梯度不等于零总溶质通量。

当前线性平衡吸附合同为：

    S_X = (theta_liquid + rho_bulk * Kd_X) * C_X,aq

S 为水相与吸附相总量 [g/cm³ 土]。只能用 C_aq 算水携通量，不能把 S 或 mg/kg
干土值直接乘排水量。吸附态不直接随液水平移，但线性平衡下可在后续子步解吸补给
水相；因此不是“吸附态永远不能淋失”。动态吸附模式须保留其显式吸附态，只运移
实际水相，不改为瞬时再平衡。

如果未来允许底部向上入流，则须给定外部 DOC 和 DON 水浓度，采用流入浓度驱动的
质量通量边界：F_X,in=q_up*C_X,external。没有外部浓度时不能自动用 0、复制底层
浓度或封住溶质却继续进水。当前长武建议先实施“仅出流模式”；出现实质性反向液流
时明确拒绝，提示需要双向模式与外部数据。零流时无边界平流交换。

出流时外部浓度不用于把内部浓度钉死为零；不采用 C_bottom=0 的吸收性 Dirichlet。
零弥散梯度是深处截断边界的近似，不是所有土壤场景的精确边界；若深层浓度梯度或
浅地下水重要，应检查加深域或独立的外部浓度/交换条件。

“流出 5 m”只表示离开当前计算域。COMPASS 没有本地地下水库状态，不能把它直接
命名为已到达真实地下水的负荷，也不能假称模拟了 5 m 以下的去向。

## 3. 水驱动中必须先补齐的接口

`QL(1)` 是最深单元中心液流，当前无机氮经 `interpolateElementToNode` 把它复制给
底节点。它不是自动成立的底面实际液流证据。

另一方面，`calculateCoupledBoundaryFlux` 恢复的是完整水方程反力，涉及凝结水和
水汽；不能不分相就把总水 QMB 称为 DOC/DON 的液态载流。

实现应在接受水子步层提供具名底面液流合同（方向、单位、物理面、dt、来源），说明
其与总水反力以及水汽通量的关系。若当前自由排水物理明确规定底部无水汽交换，需在
该模式下显式定义并验证，而不是在 DOM 层偷偷把 QMB 当液流。拒绝子步不得进入
这个序列；不能把最后一个水子步复制到整个小时。此项是上线前置验收，不以猜测代替。

## 4. 配置所有权及不动项

- 边界模式、是否允许回灌、外部入流水浓度属于 `site.toml` 的边界配置，不属于
  `soil.csv`，也不作为拟合参数。
- Kd、弥散长度等剖面参数继续放 `soil.csv`。DOC-C 与 DON-N 各用自己的相分配。
- 建议保留显式 closed 模式给闭柱/历史回归；长武新版本显式选 advective_outflow。
  这些是建议的模式名，当前没有对应实现。
- 顶部本轮保持无额外 DOM 边界交换：降水本身不凭空添加 DOC/DON，蒸发不携走
  非挥发 DOM；已有残体/有机肥/AOM 输入按原通道计入，避免重复加源。
  若以后加入灌溉 DOM 或溶解态径流，应独立定义。
- 不改 MEND 反应方程，不让 POC/MOC 或附着微生物随溶质出口流走。
  移动微生物颗粒的底部出口是另一项物理选择，本次不连带开放。

## 5. 账本和标签必须同步

DOM 输运独立守恒改为：

    M_after - M_before = boundary_input - boundary_output

DOC 与 DON 各自产生逐子步及累计的 bottom_input / bottom_output（g C 或 N/cm²），
输出可另换算 kg/ha。内部层间交换柱积分抵消；淋溶不是呼吸、矿化或数值损失。
完整 C/N 库存账也须加入这两项。现有“after−before=0”断言不能简单删除或放宽，
应改为扣除已施加边界量后的同源闭合式。

DON 肥料来源标签必须同步拥有外部有向交换、流出标签量和累计 terminal 项，维持
0≤tag≤total 与肥料总量账；当前 groundwater 为未标记来源，未来入流需沿用明确
的来源定义。不得通过关闭 tracer 来放行新出口。现有隐式有界混合方法消费冻结的
总量交换图，扩展时应包含端面流出/流入，不能再要求标签柱库存恒定。

检查点签名绑定 DOM 边界模式和外部输入定义；不能拿闭边界检查点无说明地恢复成
开边界运行。新版本另立实验身份，保存旧结果，不删除原样本。

## 6. 最小实施和验证顺序

1. 先完善实际底面液水接受序列和分相说明，核对同源水量账。
2. DOM 有限体积核增加显式端面输入，旧 closed 调用逐位不变；共享核的
   BIOPARTICLE 调用仍明确保持 closed。
3. 每个非负性自适应子步把边界通量与内部通量一起施加。不能事后扣库、截为零，
   或让标签另选时间步。DOC/DON 独立相平衡与各自质量账均闭合。
4. 同步 DOM 总账、DON 来源标签、公开诊断、配置校验和 checkpoint 签名。
5. 必须通过：单层手算出流；q=0；零水相浓度；Kd=0/强吸附；动态模式吸附态保持；
   小库存缩步而非 clamp；反向流缺浓度拒绝；显式双向浓度模式（若本轮实现）；
   非等长水子步及流向反转；拒绝/重试不重复累计；标签上下界与总账；
   closed 模式回归；连续与重启一致；时间步缩小的边界负荷收敛。
6. 做同一新源码下 closed/open 的成对基线，先短窗后延伸到首个观测日，分别报告
   累计 C/N 出口和 DOC/DON、NH4/NO3、MBC/CO2 的变化。边界正确性独立于拟合优劣。
   只有发现原敏感性排序实质改变，才决定是否补跑相关样本；不先重跑全部 44 组。

正在运行的五组试验保持原边界，继续用来辨别初态、摄取、补给机制。
新底边界是增加真实的质量去向，不能预设它会解决“DOC 已经偏低”的问题。

## 7. 外部参照及边界

- [HYDRUS 官方 Third-Type Solute](https://www.pc-progress.com/en/OnlineHelp/HYDRUS3/ThirdTypeSolute.html)：
  说明按输入溶质通量施加 Cauchy 边界，水出流/零流时切换边界形式；大气蒸发边界
  不默认带走溶质。
- [PFLOTRAN 官方 TRANSPORT_CONDITION](https://documentation.pflotran.org/user_guide/cards/subsurface/transport_condition_card.html)：
  ZERO_GRADIENT 给出入流 Robin、出流零扩散梯度；DIRICHLET_ZERO_GRADIENT 也明确
  出流只考虑平流。其 MEMBRANE_FILTER 则是水可过而溶质不可过，说明水边界和
  溶质边界是不同的物理选择。

这些支持本建议的边界原则，不是 COMPASS 已经实现这些选项的证据。
