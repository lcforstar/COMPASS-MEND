---
name: compass-coding-standard
description: "COMPASS 模型语义的数学、物理与工程规范。修改或审查方程、状态、单位、I/O 合同、边界、求解器或诊断时使用；纯格式、指令维护和非模型后处理不触发。"
---

# COMPASS V2.0 Coding Standard

本技能只以 **COMPASS V2.0 当前生产代码** 为模型事实源。它规定先如何读懂模型，再
规定如何写代码。不得从 COMPASS V1.0 档案、旧类、旧主循环或记忆直接推导当前实现。

## 1. 事实权威顺序

发生冲突时按以下顺序裁决：

1. `<当前 MATLAB 工作仓库>/src` 中从当前生产入口可达的代码；
2. 同一仓库中覆盖该活动路径的测试、站点配置和输入数据；
3. V2.0 数学物理文档和架构文档；文档与代码不一致时记录差异，以活动代码为当前事实；
4. 仅在 Julia 迁移任务中，按 compass-julia-port §1 定位正式迁移计划和目标架构文档；
   它们约束 Julia 设计，不改写 MATLAB 当前事实，也不是普通 MATLAB 修改的前置条件；

5. 参考模型用于有明确目的的科学设计和架构比较，不自动替代 COMPASS 当前物理；
6. V1.0 或其他 legacy 材料仅在用户明确要求历史溯源时读取，不得作为迁移依据。

禁止把文件名、注释、旧 dossier、搜索摘要或曾经的运行经验当成已确认事实。必须检查
调用者、被调用者、配置开关、数据形状和测试。

## 2. 按语义影响理解源码

用户当轮明确要求优先于本技能的默认流程；安全边界和更高层要求仍适用。
读取目标函数的完整相关逻辑，核对受影响调用链、配置消费者、输入输出和相邻合同。
涉及状态或耦合时检查提交/回滚与守恒；涉及初始化或诊断时检查对应时间层。
只在发现未解释依赖、跨域影响或源码矛盾时扩大阅读范围。

阅读完成条件是能解释本次改动的活动条件、物理含义、数据流、受影响合同和验证依据。
局部修改不要求从主入口重建整个模型依赖图，也不要求读取无关局部函数、迁移计划或类型设计。
完整过程迁移或跨域重构需要覆盖该过程的活动依赖闭包，并区分活动、条件活动、仅初始化、
实验关闭和遗留无调用路径。可复用已经核实且未变化的证据，补读本次变更的依赖。

实现前定义与任务相称的完成标准；新物理过程明确闭合，迁移明确 parity，局部修复明确
可判别的失败证据和验证。所有权、类型或接口设计仅在本次改动涉及它们时展开。

### 源码证据表

完整过程迁移或跨域重构记录下表；局部修改可在说明中记录受影响项，无需新建整套材料：

| 项目 | 必填内容 |
|---|---|
| 生产入口 | 文件、函数和上游调用者 |
| 活动条件 | TOML 键、switch、枚举、实验条件 |
| 输入 | 名称、类型、形状、单位、方向、时间层、来源 |
| 状态 | owner、accepted/trial/candidate、是否可重建 |
| 输出 | applied result、状态候选、ledger、诊断、失败报告 |
| 数学 | 控制体、守恒方程、本构、边界、离散、求解顺序 |
| 事务 | snapshot、commit、rollback、拒绝后不变量 |
| 验证 | V2.0 测试、reference window、parity 和守恒门 |

不适用项标记 `N/A` 并说明原因。无法确认的事实标为 `UNRESOLVED`，不得编造。
只有影响本次数学含义、单位、方向、边界、状态所有权、接受规则或目标架构的缺口才阻塞
依赖它的实现；继续完成独立工作。Julia 的 `ARCHITECTURE_CONFLICT` 先在授权范围内修订设计；
确需改变用户指定的目标架构且尚未获授权时，提交具体冲突与方案供用户决定。

## 3. V2.0 当前生产边界

以每次工作时的源码复核为准。当前基线入口是：

```text
runCOMPASS
  -> init.initialize / init.initializeRegional
  -> init.prepareRuntime
  -> beginRun
  -> runOneStep (每个 KT 的原子事务)
  -> finalizeRun
```

当前 `runOneStep` 和 `finalizeRun` 是实际生产实现，不是 V1 档案所称的 no-op/TODO。
`src/COMPASS.m`、`run_stemmus_initialize.m`、`SoilVariables`、`StateVariables` 等名称不得在
未证明属于当前活动闭包时被描述为 V2.0 生产架构。

## 4. 数学物理先于软件结构

每个过程代码必须能按下列顺序被复核：

1. domain、坐标、控制体、正方向和变量位置；
2. prognostic unknowns 与 accepted/trial/candidate 时间层；
3. governing balance 及每项单位；
4. constitutive closure、参数、假设和有效范围；
5. boundary、source/sink、request/result；
6. 空间和时间离散；
7. residual、矩阵、RHS 或 tendency 的装配；
8. 数值求解、收敛判据、限制和失败原因；
9. accept/reject、ledger、storage residual；
10. 仅从接受态生成的 diagnostics。

不能以“结构更优雅”“接口更通用”“代码更短”为理由改变方程、单位、符号、边界、
离散、求解顺序或接受判据。物理改变必须与纯重构分开提交、分开验证。

## 5. 所有权和变量分类

每个量必须先归类，再决定容器：

| 类别 | 规则 |
|---|---|
| parameter | 运行时只读；修改意味着重建模型 |
| prognostic state | 不能由其他量重建；只有事务提交层可写入接受态 |
| forcing | 外部驱动；每步切片，在同一 trial/retry 中冻结 |
| auxiliary/cache | 可由接受态、参数和 forcing 重建；不得充当第二套状态 |
| request | 相邻域提出的需求，不等于实际施加量 |
| result | 过程实际施加或求解得到的量 |
| ledger | 与状态使用同一离散控制体的守恒账本 |
| diagnostics | 只消费接受态/result/ledger，不得反向驱动物理 |
| runtime memory | 无法由状态重建但又不属于物理场的运行记忆；checkpoint 合同须明确 |

禁止：

- 通用 `workspace_flux`、`scope_flux` 或可跨域随意写入的 flux 桶；
- 把 request 当作 applied result；
- 让 diagnostics 或输出文件反向决定物理；
- 让 auxiliary 在拒绝步后被当作接受态；
- 建立多个可写的同一物理状态副本；
- 为 Julia 复制一个通用可变 `model.built` 仓库。

## 6. 事务、边界和守恒

- 外层 KT 是原子事务：候选过程全部成功后才 commit。
- 拒绝必须保持接受态、接受时间、累计 ledger 和事件状态位级不变，除明确允许的审计信息。
- trial 可以更新预分配 workspace，但不得泄漏到 accepted state。
- 边界合同必须携带类型、单位、正方向、参考位置、时间区间和需要的导数/切线。
- request 经过 solver/materializer 后产生 applied result；诊断记录二者差异。
- 累计通量只能在接受后进入 ledger。
- 禁止用状态 clamp 掩盖质量/能量误差；限制器必须有物理含义并进入 ledger。

涉及土壤水、热、冻融、气体和顶边界时，强制阅读
`docs/soil/compass_soil_water_heat_physics.md` 的相关章节，至少复核 §2.8.1 和 §5.7。
生产 `heat_top_kind = atmospheric` 是 matched fixed `G*` Neumann；不得写成 SCOPE skin
temperature Dirichlet，除非目标站点明确选择其他活动边界。

## 7. 空间、单位和时间层

- V2.0 土柱内部方向保持现有约定：深度为正，数组 bottom-to-surface，node 1 最深，末节点为表面。
- Julia 不得擅自反转该方向；倒序只在明确的外部适配/输出边界按合同发生。
- 每个数组说明 node/layer/endpoint/tile/cell 位置和形状。
- 单位转换只允许出现在命名明确的 adapter/boundary 函数中。
- 变量名存在单位歧义时添加后缀，例如 `_s`、`_K`、`_cm`、`_W_m2`。
- 明确区分 configured forcing time、row index、accepted model time 和 diagnostic timestamp。
- 不能假定输入文件中存在的列或 TOML 键已被模型消费；必须追踪到读取点。

## 8. MATLAB V2.0 实现规则

- 新物理核使用显式值输入/输出；禁止新增 `global`、`assignin`、`evalin` 或隐式工作区依赖。
- package 边界反映物理 owner，不以方便为由跨域写状态。
- 新接口优先小而 typed 的 request/result/ledger；旧形状兼容仅留在 adapter。
- 保持接受态和 trial 分离；快照与恢复由明确的事务层负责。
- 配置中写了但没有生产消费者的键，必须 fail-fast、明确弃用或实现真实消费，不能静默忽略。
- 旧函数只有在当前生产闭包或明确历史对照中才可修改；不得把清理 legacy 当成迁移进度。

## 9. Julia 实现规则

编写 `.jl` 时还必须使用 `compass-julia-port`。

- 采用具体、参数化的不可变参数类型；浮点类型统一由 `FT` 传播。
- 预报状态由物理域拥有；aux/workspace 可重建；Step 是唯一 commit/rollback 层。
- 热循环使用就地函数和预分配；不得为了零分配牺牲数学可读性或事务正确性。
- 不要求把所有过程塞入统一 `(dY,Y,p,t)`；先复现 V2.0 的离散和求解语义，再为合适的纯 tendency 提供 SciML 接口。
- 不得用 DifferentialEquations.jl、ClimaCore 或参考模型默认算法替换 COMPASS stepper，除非单独证明数值等价。
- MATLAB 未启用但存在的过程可完整迁移为 disabled 状态，但不得自动启用或伪称通过生产 parity。

## 10. 文档、来源和署名

公共模型函数/docstring 至少说明：

- 输入和返回值的类型、形状、单位、方向和时间层；
- 方程/算法及关键假设；
- 改编或迁移时记录源实现路径和函数；
- 有文献性主张时给出可核实来源；未核实的已有来源标 `[SOURCE UNVERIFIED YYYY-MM-DD]`，不编造 DOI 或公式编号；
- 错误和拒绝条件。

新增模型函数的文档末尾使用当前编辑者署名（MATLAB 加 `%`，Julia 加 `#`）：

```text
Editor: Wenjie Liu; Date: YYYY-MM-DD.
```

修改已有模型函数时保留历史作者与 Change 记录；本次实质性修改追加：

```text
Change: YYYY-MM-DD Wenjie Liu — 说明本次物理、合同或文档含义变化。
```

没有既有署名的函数，在本次实质修改时补 Editor。纯格式调整不批量补署名或追加 Change。
这些是未来编辑的模板，不授权批量替换历史作者。

## 11. 验证与完成门

先按改动影响选择必要检查，复用未失效的验证证据：

- 局部 MATLAB 修复：可判别的回归检查和受影响合同；不要求 Julia 类型稳定性或迁移 parity。
- 新物理或跨域修改：适用的单位、极限、守恒、失败/rollback 和耦合短窗检查。
- Julia 过程迁移：源证据、完整输入输出、相关 parity 与接受/拒绝轨迹；关键路径检查类型与分配。
- 全年、区域或资源测试：仅在交付范围、阶段完成门或变更风险需要时运行。
- 只按实际验证范围声明能力；初始化容器或单元测试通过不等于区域/全球生产运行通过。

必要检查通过后，只因新修改、失败或未解决疑点扩大或重跑。完成指当前请求的交付物完成，
并不要求每次把整个物理过程重新迁移一遍。既有无关失败单独报告，不扩展为无关修复。

如果测试失败，修复原因；不得放宽容差、跳过测试、改基准或用 NaN/零值占位来宣布完成。
