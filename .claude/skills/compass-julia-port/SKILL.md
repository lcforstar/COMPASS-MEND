---
name: compass-julia-port
description: "从当前 COMPASS MATLAB 实现迁移、设计或审查 COMPASS Julia 模型过程。保留物理与数值语义；普通 MATLAB 修改、无关 Julia 文件和纯格式调整不触发完整迁移流程。"
---

# COMPASS V2.0 → Julia Port

本技能不提供“凭旧档案直接翻译”的捷径。它强制执行：**先完整理解 V2.0 生产实现，
再形成数学物理合同和 Julia 重构设计，最后才编码。**

同时使用 `compass-coding-standard`；数学物理和源码证据服从前者，Julia 迁移流程服从本技能。

## 1. 源码与目标路径

`<MATLAB_ROOT>` 是本任务确认的当前 MATLAB checkout；通常为当前工作仓库。
`<JULIA_ROOT>` 是用户指定或已有项目配置明确的 COMPASS.jl checkout。`<MIGRATION_PLAN>`
使用用户指定路径，或项目中实际存在的 COMPASS_JULIA_MIGRATION_PLAN.md。不得假设旧用户名、
Linux 路径或相邻同名目录就是正确版本。只在已知工作区内定位相关文件；必要依据缺失时，
明确阻塞的设计部分并继续独立工作，不开展全盘搜索。

首次迁移一个过程需要下列证据；后续局部修改只更新受影响的证据、合同与测试。
用户当轮明确要求优先于默认流程，已有授权和答案持续有效。

### 过程迁移源码门

首次实现或实质重构一个 Julia 模型过程前：

1. 完整阅读 COMPASS.jl `AGENTS.md` 和以下三份强制目标依据中与任务有关的完整章节：
   - `<MIGRATION_PLAN>`
   - `<JULIA_ROOT>/docs/architecture/compass_complete_architecture.md`
   - `<JULIA_ROOT>/docs/architecture/REFERENCE_MODEL_ARCHITECTURE_REVIEW.md`
2. 阅读 `references/00-v2-source-authority.md`，建立本任务的 V2.0 源码证据表。
3. 从 `<MATLAB_ROOT>/src/runCOMPASS.m` 或正式独立入口追踪目标过程的
   活动依赖闭包。
4. 完整阅读相关 MATLAB 文件及局部函数、调用者、直接被调用函数、配置、输入适配器、
   状态构造、事务提交/回滚、诊断和测试。
5. 标记每条路径为 `active`、`conditional_active`、`initialization_only`、
   `implemented_disabled` 或 `legacy_unreachable`。
6. 写出控制体、状态、方程、本构、边界、离散、求解顺序、接受规则和守恒账本。
7. 设计 Julia owner、typed contracts、文件归属、迁移门和 parity fixture。
8. 对照三份正式文档复核目录、所有权、空间层、接口、实施顺序和完成门。
9. 影响本次数学、状态、边界、接受规则或架构的缺口必须先解决；不适用项标 N/A，其他缺口仅阻塞依赖它的部分。

不得根据函数名、注释、旧 dossier 或 V1.0 类结构猜测。不得把“代码存在”等同于“生产启用”。

## 2. 唯一事实源和禁止项

- MATLAB V2.0 当前生产代码和测试是方程、离散、单位、方向、求解顺序和失败行为的 oracle。
- 上述迁移总纲和两份架构 MD 共同规定 Julia 目标架构、状态所有权、空间层、实施顺序和
  完成门，是实现的约束而非建议；不能覆盖 MATLAB 当前事实，也不能被代码自行背离。
- CLM5/CTSM、ClimaLand、CoLM 用于架构对标，不得替换 COMPASS 物理或 stepper。
- V1.0 只允许在用户明确要求历史溯源时读取；不得从 V1 dossier 迁移任何 Julia 代码。
- 不得复制 V1 的 `SoilVariables`、`StateVariables`、旧 `COMPASS.m` 主循环或 handle 类架构。
- 不得把 V2 的 `model.built` 原样复制成 Julia 的通用可变仓库。
- 不得预先决定所有过程都适合统一 ODE tendency；先重现 V2 数值语义。

源码、测试或正式文档的矛盾影响本次合同且无法在既有授权内解决时，暂停受影响实现并记录冲突；继续独立工作。
不得在代码里静默选择其中一方，也不得先写代码再让文档追认。

## 3. 带日期的 V2.0 基线（核对受影响路径）

截至 2026-07-22，生产生命周期是：

```text
runCOMPASS
  -> init.initialize / init.initializeRegional
  -> init.prepareRuntime
  -> beginRun
  -> runOneStep
  -> finalizeRun
```

`runOneStep` 已是实际的原子外步事务：分阶段构造地表—植被候选、推进耦合土壤
`h/U/P_g`、在土壤接受后推进 N candidate、统一 commit，失败则恢复 snapshot。
`finalizeRun` 已实际写 diagnostics。旧档案中“runOneStep no-op、finalizeRun TODO、
`src/COMPASS.m` 才是真实主循环”的结论已经失效。

初始化和运行事实见 `references/01-v2-production-pipeline.md`。这些内容是带审计日期的基线，
不是替代重新读源码的永久真理。

## 4. 迁移顺序

服从 `<MIGRATION_PLAN>` 的当前阶段门。原则顺序是：

1. P0/P1 架构、reference windows 和公共 contracts；
2. P2 完整输入层：raw/provenance、验证、派生、typed config、forcing、parameters 和初态；
3. P3 大气；
4. P4 植被；
5. P5 地表；
6. P6 土壤水；
7. P7 水—热；
8. P8 水—热—冰；
9. P9 水—热—气；
10. P10 土壤/植被 C/N；
11. P11 完整输出持久化层：registry/history、NetCDF/CSV、checkpoint/resume 和 finalize；
12. P12 耦合、求解器、Step、Conservation 和完整集成运行；
13. P13 多站点、区域和全球执行扩展。

不能因为某个公式容易而越过输入或前置模型阶段。P3–P10 必须定义并测试内存中的
result/ledger/diagnostics contract，但正式 CSV、NetCDF、checkpoint 和其他持久化只能由
P11 实现。代码中存在但当前未启用的过程，在所属阶段完整迁移并保持 disabled；不能提前
接入生产 driver。

## 5. 每个过程必须迁移完整输入和完整输出

### 完整输入

- immutable parameters/configuration；
- accepted prognostic state 和必要旧时间层；
- immutable current-step forcing；
- boundary/source/management request；
- auxiliary/cache 和预分配 workspace；
- 单位、符号、控制体、时间区间、数组位置、vertical orientation；
- V2.0 来源函数、活动调用路径、开关和实验条件。

### 完整输出

- trial/candidate prognostic state；
- 实际 applied result，而不是 demand/request；
- 唯一 conservation ledger 和累计积分；
- convergence/iteration/accept/reject/failure report；
- accepted-state diagnostics；
- parity fixture、守恒、rollback、失败和 disabled-path 测试。

只迁移公式、只建类型、只返回恒零值或只复制输出字段都不算完成。

## 6. Julia 重构边界

具体合同见 `references/02-v2-to-julia-refactor-contract.md`。核心规则：

- 每个 prognostic state 只有一个物理域 owner；Column 是最小物理、守恒和原子事务单元。
- Parameters 运行时不可变；ForcingAtTime 在 retry 内冻结；Aux 可重建；Step 独占 commit。
- request、applied result、ledger、diagnostics 分开建 typed records。
- Coupling 只编排交换，不拥有长期状态；Solver 不 commit；Diagnostics/IO 不反向驱动物理。
- 保持 V2.0 纵向约定：深度为正、数组 bottom-to-surface；Julia 内部不擅自反转。
- 用具体 `FT` 类型、明确形状和预分配；避免抽象字段、动态字典和通用 flux bucket。
- 只有适合纯 tendency 的过程才提供 `(dY,Y,aux,t)`；事务式/嵌套求解过程保留显式 solver API。

## 7. 空间尺度必须如实声明

分开声明：

- laboratory：V2.0 当前有 point + laboratory 的 water-only 初始化/运行合同；
- point/site：当前生产闭合目标；
- regional：V2.0 有区域容器初始化和 override 测试，但生产 `h/U/P_g` 求解仍对 `nCell>1`
  fail-fast；不能称为区域模型已可运行；
- global：当前 V2.0 没有全球生产 runtime；Julia 只可先设计空间/执行合同，不得宣称实现。

站点、实验、区域和全球共享同一 Column 物理语义；空间层只负责索引、参数分布、forcing、
聚合、调度和 I/O，不能改变柱内方程和接受判据。

## 8. 未启用与遗留代码

采用以下状态词：

| 状态 | 含义 |
|---|---|
| `enabled_validated` | V2 活动 oracle、parity、守恒和 rollback 已通过 |
| `implemented_disabled` | 完整实现和测试，默认关闭 |
| `implemented_unverified_disabled` | 已实现但没有活动 oracle，禁止生产启用 |
| `not_implemented` | 尚未迁移 |

`legacy_unreachable` MATLAB 代码默认不迁移。若用户要求保留，先证明其物理价值和依赖，
再作为 disabled capability 单独设计，不能污染活动路径。

## 9. 验证顺序

1. 输入 raw fidelity：MATLAB 与 Julia 读取值、类型、缺省和 provenance 对齐；
2. derived input parity：网格、时间轴、forcing、参数和初态对齐；
3. pure relation parity：覆盖有效域、边界、极端值和单位转换；
4. 单过程短窗：state/result/ledger/failure 对齐；
5. 事件窗：降雨、冻融、积雪、施肥、物候和边界切换；
6. accepted/rejected trajectory 和 rollback parity；
7. 仅在当前阶段交付或风险需要时运行全年/资源测试；
8. 当前阶段适用的门通过后才接入下一层；局部修改不重跑无关阶段或未失效的检查。

容差必须来自数值分析和 V2 reference window，不能事后为通过测试而放宽。

## 10. 按需读取的参考档案

| 场景 | 必读文件 |
|---|---|
| 任意 Julia 迁移或设计 | `references/00-v2-source-authority.md` |
| 输入、初始化、runtime、主循环、输出、checkpoint | `references/01-v2-production-pipeline.md` |
| 类型、owner、模块、空间层、迁移完成门 | `references/02-v2-to-julia-refactor-contract.md` |
| 编写或审查任意 `.jl` | `references/03-julia-coding-standard.md` |

物理细节不再存放静态 V1 dossier。迁移某个过程时必须从 V2.0 当前源码重新生成该过程的
source map，并将稳定结论写入项目正式架构/迁移文档，而不是把旧快照重新塞回技能。
