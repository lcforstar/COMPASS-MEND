# COMPASS–MEND正式耦合候选基线

**后续变更：** 用户于2026-09-16指定DOC/DON采用平衡吸附，当前工作区已据此修改。本文归档仍为动力学吸附RC1历史基线，其文件核验报告当前工作区差异属于预期。新设置和验证见[平衡吸附说明](MEND_EQUILIBRIUM_SORPTION.md)，不可将本页全年结论直接移用于新模式。

版本：`compass-mend-rc1-20260916`。本页是截至2026年9月16日的当前状态入口。旧文档中的“等待恢复输入”和失败记录属于阶段历史。

## 固定了什么

这是一份通过本站点数值验收的候选基线：MEND负责层内有机碳氮、微生物及酶反应；COMPASS负责环境、无机氮、根竞争、供氧、分层、运输及外部输入。它不等于完整复制独立MEND，也不要求耦合后的年末库存与独立MEND相同。

冻结包位于[候选目录](baselines/compass-mend-rc1-20260916/)，包括`candidate.zip`和`manifest.json`。793个文件覆盖716个正式源码文件、完整长武站点目录、选定验收驱动和本轮结果证据。文件逐项使用SHA-256标识；当前工作区源码与已通过的全年源码完全一致。归档创建后已重新验证压缩包、全部文件及当前工作区，一致性检查通过。

仓库存在此前工作留下的未提交改动，因此版本身份以文件哈希为准，不能只用Git提交号代替。此次未改模型公式或站点参数。

## 适用的运行设置

| 项目 | 本次验收配置 |
|---|---|
| 站点 | 长武冬小麦，`data/sites/changwu/site.toml` |
| 时段 | 2000-06-15至2001-06-15，共8760小时 |
| 土柱 | 500 cm、58个非均匀土层 |
| 有机反应 | MEND小时源核；耦合全年采用iHR=0 |
| 初始有机库 | COMPASS站点分层初态，含先验分库假设 |
| 外部输入 | 现有有机肥、AOM和植物凋落物均启用 |
| 其他耦合 | 实际环境、根系氮预留、无机氮、氧约束、DOM和微生物运输 |
| 重启 | 必须保留根系吸水历史；使用当前兼容断点 |

## 已经通过什么

| 检查 | 结论与证据 |
|---|---|
| 独立源核复现 | iHR=0/1分别完成8760小时；这是源核证据，不是两种模式的耦合全年证据，见[MEND验收记录](MEND_ACCEPTANCE.md) |
| 暂缓新输入的耦合全年 | 收支和运输通过，见[初始库全年验收](FORMAL_ANNUAL_ACCEPTANCE.md)；该历史运行曾暴露重启历史缺失问题 |
| 重启历史修复 | 保存和恢复上一小时根系吸水通量；修复后短窗及359小时对照通过，见[修复记录](ROOT_WATER_RESTART_FIX.md) |
| 有机肥/AOM全年 | 8760小时输入数量、放置、收支及施肥后重启通过，见[有机肥验收](MANURE_INPUT_ACCEPTANCE.md) |
| 有机肥+AOM+凋落物全年 | 8760小时通过；真实凋落物431.7798 kg C/ha、10.28047 kg N/ha，见[凋落物验收](LITTER_INPUT_ACCEPTANCE.md) |
| 跨凋落物事件重启 | 8233–8760共528小时；事件后立即重启一致，3个完整断点一致，229个输出变量逐小时一致（含108个剖面变量），见同上 |

主要通过结果保存在`annual_runs/manure-input-annual-20260915/`、`annual_runs/litter-input-annual-20260915/`和`annual_runs/litter-event-restart8232-20260916/`。前两次后期重启驱动失败是检查文件和事件窗口选择错误，保留原记录，不作为通过证据。

## 如何使用与核对

正式调用顺序是`runCOMPASS` → `init.initialize` → `init.prepareRuntime` → `beginRun` → 每小时`runOneStep` → `finalizeRun`。正式入口在`src/runCOMPASS.m`；各模块的责任和输入顺序见[中文接口说明](../soil/biogeochemistry/MEND_SOURCE_PRODUCTION_CONTRACT.md)。`tools/diagnostics/`内脚本是验收工具，不需要加入生产物理模块。

先在PowerShell核对当前工作区是否仍对应此候选版本：

```powershell
python 'C:\Users\Administrator\Desktop\COMPASS-MEND\tools\diagnostics\verifyMendCandidate.py'
```

输出要求`archive_verified=true`，两个工作区差异列表为空。此命令只检查文件，不启动模型。

要重复已通过的全年验收，可使用现有入口并指定新的目录名：

```powershell
powershell -NoProfile -ExecutionPolicy Bypass -File 'C:\Users\Administrator\Desktop\COMPASS-MEND\tools\diagnostics\runMendLitterAnnual.ps1' -RunName 'litter-annual-repeat-20260916'
```

该入口清理相关进程环境变量，采用独立源码副本和输出目录，保留实际施肥与凋落物。它是带收支检查的验收运行；日常模型运行则调用正式`runCOMPASS`。更改输入、参数或过程开关后应建立新的实验目录，不再把结果称为此候选配置的重复运行。

## 尚未覆盖的范围

- 本站点全年未触发氧反应缩放。零氧/部分供氧构造检查不能替代长期真实缺氧情景验证。
- 尚未独立覆盖非零修剪、木质死亡等输入事件，也没有完成植物自身全部C/N预算的独立审计。
- 原全年没有8760完整Y断点，最终小时验证限于已保存输出；三个中间断点比较的是完整内部状态。
- 此次才完整归档站点输入。此前运行仅对源码、驱动和`site.toml`做了哈希核验，不能追溯证明全部气象和其他输入文件在运行期间未变化。归档不包含MATLAB运行时和站点目录外的全部依赖，不能视作独立可迁移运行包。
- 数值守恒及重启一致不等于参数观测验证；初始分库、肥料参数与材料分解参数仍需科学校验。其他站点、不同管理情景和多年稳定性不在本轮结论内。

这份候选基线可用于后续有明确控制变量的科学实验。若目标是跨机器长期复现，下一项是补齐运行实际读取的数据依赖清单与运行前哈希校验；若目标是扩展科学适用范围，则按新增情景验证，避免无改动地重复已有全年试验。
