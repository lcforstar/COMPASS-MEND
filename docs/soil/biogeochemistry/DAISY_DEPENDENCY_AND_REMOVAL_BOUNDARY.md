# DAISY 依赖清单与删除边界

日期：2026-08-29

## 1. 结论与边界

本次清理删除的是旧 DAISY 有机质生产对象、状态拓扑、配置入口和专属预热路径，
不是删除所有以 DAISY 为科学来源、对照模型或验证资产的实现。

以下内容明确保留：

- DAISY 来源的根氮吸收方法及其 Nacry 耦合；
- DAISY 来源或 DAISY-style 的硝化、反硝化动力学；
- 植物碳分配、NSC、周转和物候中的 DAISY 文献或源码参考；
- `tools/verification/` 下的 DAISY 对照与复现工具；
- `daisy_om_matlab/` MATLAB 对照实现；
- `daisy-7.1.0/` 上游参考源码。

这些保留项不得拥有 COMPASS 原生有机碳生产状态，也不得重新成为
`runCOMPASS` 或 `runCOMPASSBiogeochemistrySpinup` 的有机质反应内核。

## 2. 已删除的生产依赖

| 类别 | 已删除或替代内容 | 当前所有者 |
|---|---|---|
| 生产入口 | `runCOMPASSOrganicSpinup.m` | `runCOMPASSBiogeochemistrySpinup.m` |
| 旧动态预热 | `soil.organic.runDynamicSpinup` | 原生 BGC 循环控制器 |
| 旧反应候选 | `soil.organic.proposeStep` 及 DAISY SMB/SOM/AM 对象 | AOM、DOM、MEMC/MEND 原生过程 |
| 主开关 | `soil_organic_cycle` | `soil_carbon_model` |
| 初始化模式 | `daisy_equilibrium` 和 DAISY 池分配 | 测量 SOC + POC/MOC/MB/IOM 先验 |
| 状态 | SMB、SOM、DAISY AM/AOM 容器和 buffer | 原生 AOM/POC/MOC/MB/DOM/MEND 状态 |
| 参数 | DAISY 专属慢/快 SMB、SOM、AOM 周转配置 | `soil.csv` 与原生 `[organic]` 参数 |

## 3. 保留但必须隔离的来源依赖

| 范围 | 状态 | 允许用途 | 禁止用途 |
|---|---|---|---|
| 根氮吸收 | 生产保留 | DAISY 单根供给与 Nacry 膜通量耦合 | 拥有土壤有机碳状态 |
| 硝化/反硝化 | 生产保留 | 来源明确的反应形式 | 恢复 DAISY 有机质反应器 |
| 植物分配/物候 | 生产保留 | 文献、公式和行为参考 | 构造 DAISY 土壤池 |
| `tools/verification/` | 验证保留 | 离线 parity、历史复现 | 加入生产 MATLAB path 作为反应 owner |
| `daisy_om_matlab/` | 参考保留 | 对照实现 | 从生产入口调用 |
| `daisy-7.1.0/` | 上游源码保留 | 科学溯源、离线验证 | 编译或运行时生产依赖 |

## 4. 原生术语不得误删

`AOM` 是 COMPASS 当前凋落物、秸秆和粪肥添加有机质的原生科学概念；
不能仅因名称也见于 DAISY 而删除。`POC`、`MOC`、`MB`、`DOM`、`MEMC`、
`MEND`、`BIOPARTICLE` 及其 C/N 状态均由 COMPASS 当前生产路径拥有。

## 5. 验收扫描规则

搜索 `DAISY` 只能生成候选清单，不能直接生成删除列表。每个命中必须归入：

1. 生产保留的科学来源；
2. 离线验证或上游参考资产；
3. 历史文档；
4. 错误所有权文字，需改为 native/COMPASS；
5. 旧有机质生产依赖，必须删除。

完成门为：当前生产有机质闭包中不存在第 5 类；第 1、2 类符合本文件隔离边界；
原生 AOM 等同名科学概念没有被误删。
