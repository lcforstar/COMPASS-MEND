# DOC、DON 与微生物边界类型

Editor: Wenjie Liu; Date: 2026-09-29.

本次只改有机物输运及其来源标签/诊断，不修改 NH4/NO3 方程、接口、输入或边界配置。
配置与物理含义参照当前无机氮的可解释分支；不把有限元节点矩阵直接复制给有机物有限体积单元。

## 1. 当前长武设置

DOC、DON、BIOPARTICLE 均为 **top=2，bottom=3**；外部浓度为零。
含义：没有外部有机物随入渗输入，不因蒸发损失非挥发性有机物；底部零浓度梯度，
水向下排出时携带底层水相/移动态物质。零梯度不是零通量。

DOC/DON 的总量与游离水相仍由各自 Kd 和液态水量决定；只运走水相，不直接搬走吸附库存。
微生物只运移 MB_mobile、BD_mobile 和 ENZN_mobile 的 C/N。attached、retained 不直接跨边界；
原有附着、脱附、截留算子不变。胞内酶仍与细胞共享输运类型和单调时间步，不变为独立自由酶。

唯一生产配置是 `data/sites/changwu/site.toml` 的 `[boundary]`，没有改变历史站点副本和 campaign。

## 2. 可选类型与限制

| 边界 | 类型 | 本次有机物实现 |
|---|---|---|
| 顶部 | 1 | 指定外部水相浓度；半单元扩散/弥散梯度与迎风对流计算实际交换 |
| 顶部 | 2 | 入渗时按外部浓度输入；蒸发不携走有机物；零外部弥散通量 |
| 底部 | 1 | 指定外部水相浓度；可有向内或向外的浓度梯度交换 |
| 底部 | 2 | 水流进入时使用指定外部浓度；水流排出时零梯度随流输出 |
| 底部 | 3 | 零梯度、只出不进；反向水流报错，必须选有明确外部浓度的边界 |

**不是“全部与无机氮逐行相同”。** 不复制的差异/限制：

- 无机氮类型 1 钉节点浓度；这里指定边界面浓度，不覆盖整层总库存，实际交换必须入账。
- 无机氮顶部类型 3 的装配与通量恢复并不采用同一明确语义。因此有机物顶部只开放 1/2，
  拒绝 3；无机氮 schema 曾列出的 4 没有对应生产分支，也不添加有机物 4。
- 底部有机物保留已经验收的 **接受物理边界 QMB + 显式零水汽 liquid_only**。
  不为仿照 DIN 而把它悄悄换成底单元 QL。顶部按 DIN 的元素到节点适配方式采用接受的顶单元 QL，
  不把含水汽的总 QMT 当成液态水。它不是一个新增、独立求解的顶界面液态水反力。
- 目前生产水边界联动仍限定 `water_bottom_kind="free_drainage"`；定水头、地下水反补等
  需要另行定义水相通量合同，不在本次启用。底部类型 2 的正向入流在独立算子测试中覆盖，
  **不代表已经完成地下水反补的整模型验收**。
- 类型 3 的反向拒绝比 DIN 无条件施加梯度约束更严格；不凭空假设柱外 DOC、DON 或细胞组成。

## 3. 配置和单位

```toml
[boundary]
doc_top = 2
doc_bottom = 3
doc_top_value_mg_L = 0.0
doc_bottom_value_mg_L = 0.0
don_top = 2
don_bottom = 3
don_top_value_mg_L = 0.0
don_bottom_value_mg_L = 0.0
dom_bottom_water_phase = "liquid_only"

bioparticle_top = 2
bioparticle_bottom = 3
bioparticle_top_carbon_mg_L = [0.0, 0.0, 0.0]
bioparticle_top_nitrogen_mg_L = [0.0, 0.0, 0.0]
bioparticle_bottom_carbon_mg_L = [0.0, 0.0, 0.0]
bioparticle_bottom_nitrogen_mg_L = [0.0, 0.0, 0.0]
bioparticle_bottom_water_phase = "liquid_only"
```

单位是 **mg C 或 N/L 水**，不是 mg/kg 土，也不是细胞个数；解析时一次乘 `1e-6` 成 g/cm³ 水。
微生物向量顺序为 `[MB, BD, ENZN]`，每一项 C/N 必须配对给出，不能用缺省 C:N 猜测。
ENZN 外部输入必须同时有细胞输入；输入组成还应符合启用的生物反应模型的化学计量适用域。
类型 3 的外部浓度不参与计算，要求显式填零，避免非零配置被静默忽略。
未知字段、负值、NaN、缺值、未实现类型均拒绝。

旧 `dom_bottom_kind` / `bioparticle_bottom_kind` 的 closed、advective_outflow 配置仍可独立使用，
但**不得与同一 owner 的新编号键并存**。需闭柱对照时删去该 owner 的新编号/浓度键，
改用旧 `*_bottom_kind="closed"`，从独立新运行开始。旧封闭时间步策略不改。

## 4. 数学、所有权与账本

COMPASS 单元 bottom-up；q 正向上。对任一边界定义 `q_in` 进入土柱为正：
底部为 QMB，顶部为 -QL_top。单元水相浓度为 C，外部浓度为 Ce。
类型 1 的入量率为 `(max(q_in,0)+a) Ce`，出量率为 `(max(-q_in,0)+a) C`，
其中 `a=(theta D + alpha abs(q_in))/(dz/2)` [cm/s]。常规内部面算子不变。
类型 2 无外界弥散；顶部不输出非挥发物，底部按对流方向区分输入与输出。

每个接受输运子步显式保留四路质量：top_inflow、top_outflow、bottom_inflow、bottom_outflow。

`storage_after - storage_before + top_outflow + bottom_outflow - top_inflow - bottom_inflow = roundoff`

DOC/DON 可自适应细分，微生物所有 C/N/ENZN 共用按出流系数计算的单调时钟。
只有接受的水子步进入驱动表；累计量先留在候选状态，外层 KT 成功后一起提交。
底部累计排出沿用原有 bottom_export_cumulative；其他三路放在 boundary_exchange_cumulative，
不建立底部累计量的第二可写副本。

来源标签只消费总量接受台账。外部 DON/细胞 N 是**未标记**输入，只稀释矿质肥料来源比例；
有机顶部输出是新增终端损失，不能归到内部循环。肥料标签 schema 从 3 升到 4。
边界类型、所有载体外部浓度及合同版本进入已有 checkpoint compatibility signature，
不兼容断点拒绝，不静默迁移旧 checkpoint。

## 5. 生产路径与输出

`io.checkSite → normalizeBoundary/normalizeTransportBoundary → init.buildSoilDOM/buildSoilMicrobes`
初始化只读边界与累计状态。`solver.advanceSoilEquationOuterStep` 保留接受水驱动，
`acceptedBottomDrivers → materializeTransportBoundary → transportConservative`
计算 DOM 或 BIOPARTICLE 的实际输入/输出。无机氮 `advanceAcceptedCandidate` 的既有调用点不改。

`soil_organic.nc` 保留原来的底部 step/flux/cumulative，并新增 DOC、DON 和
`bioparticle_{mb,bd,enzn,total}_{carbon,nitrogen}` 的三个其他方向，每方向 step/flux/cumulative，
共 90 个新增 float64 变量。单位 g element/cm² 或 g element/cm²/s。
标签增加 `fertilizer_tracer_{don,bioparticle_mb,bioparticle_bd,bioparticle_enzn}_top_loss_{step,cumulative}`。

## 6. 验证证据

测试入口 `tools/diagnostics/verifyOrganicBoundaryTypes.m`；证据目录
`outputs/organic_boundary_types_20260929/`。修改前源码、站点和模型库冻结在 `frozen/`，
不使用含历史干预的六年 reference 当基线。

- `unit`：96 个可接受的单元手算组合（入渗/蒸发/静止、底部正负零流、两种 Kd），
  输入与双向扩散、闭合、未标记输入稀释、非法配置/逆流拒绝；默认路径精确回归。
- `capture` 与 `production`：修改前后各 48 h，比较全部接受态，只有新增零累计字段和
  标签 schema 元数据例外；独立全有机 C/N 账、KT24→25 断点一致、配置不兼容拒绝、输出核对。
- `configured`：顶部 1/2 × 底部 1/2/3 六组生产短窗，显式外部浓度只在测试 raw 配置中注入；
  合成肥料标签探针验证非零顶部有机 N 终端损失，没有改实际站点初态或观测。
- `restart_nonzero`：三个新增方向的累计量均非零，KT1→2 保存/续跑接受态逐位一致。
- DIN 49 个非 tracer 源文件及站点 NH4/NO3 边界字段与本次冻结快照逐字节/逐字段对比。

以上五个模式均已 PASS。默认 48 h 独立全有机 C/N 最大闭合残差分别为
`8.43309622028742e-16` / `1.2402525950500066e-16` g/cm²，肥料标签最大残差
`7.665648172651418e-19` g N/cm²；新增 90 个诊断变量通过输出检查。
六组非零配置最大 C/N 残差分别为 `3.3294325381213764e-16` /
`7.721390964232693e-17` g/cm²。96 个手算组合最大核残差 `2.0868112493487782e-20` g/cm²。
强流颗粒测试通过共同时间步、附着/截留不直接输出、重复调用幂等和有界标签检查。
旧接受子步重放、80 组随机有界来源输运、肥料输运诊断回归也通过。

首轮 production 的积分和守恒完成后，测试脚本因本机无 `struct2array` 中止；
已改用基础 MATLAB 的 `cell2mat(struct2cell(...))` 并重跑通过，未更改模型或测试容差。

具体结果以同名 `*_verification.json` 为准；旧底部实现记录及其证据
保留历史归属。本次不启动年度筛查、参数率定或六年基线，也不把工程验收称作科学率定完成。
