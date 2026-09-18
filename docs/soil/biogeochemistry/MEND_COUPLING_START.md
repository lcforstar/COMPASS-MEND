# 按接口约定固定正式耦合路径

**当前状态（2026-09-16）：** 真实有机肥/AOM、凋落物全年和跨事件528小时重启对照均已通过，候选基线已冻结。使用说明与验收范围以[正式耦合候选基线](../../model_comparison/MEND_COUPLED_CANDIDATE.md)为准。下文保留最初暂缓输入阶段的运行记录，不代表当前仍需关闭输入。

更新：2026年9月15日。

当前仓库已经具备MEND小时反应与COMPASS主循环之间的调用。这个阶段的工作是验证这些接口在同一份正式源码中一起工作，不再使用临时源码中的输入替换、休眠补回或100厘米隔离。

## 本阶段运行什么

使用长武站点的分层初始有机库和实际气象、水热、无机氮、根系、氧气及层间运输。保留公共氧气反应比例。为遵守分阶段恢复约定，运行入口暂时屏蔽新的有机肥事件和植物凋落物输入，初始AOM必须为空。不注入独立MEND的均匀源输入，也不修改原站点文件。

这是“已有土壤有机库在实际环境中的耦合检查”。它不代表已经验证站点全年真实有机物补给。

## 如何重复检查

在仓库根目录的PowerShell执行：

```powershell
powershell -NoProfile -ExecutionPolicy Bypass -File .\tools\diagnostics\runMendContractAcceptance.ps1 -Hours 168
```

入口先检查冻结参考数据，再运行源有机反应和氮接口回归，最后使用正式`runOneStep`进行组合检查。结果写到独立的`docs/model_comparison/contract_runs/时间戳/`目录。

- `status.json`：是否完成、退出码、运行范围。
- `source_baseline.json`、`din_interface.json`：源反应和氮接口结果。
- `production.json`：主循环、完整状态重启一致性、零氧/部分供氧和失败候选检查结果。
- `code-hashes.json`：本次源码与配置记录；运行中发生变更会使本次验收失败。
- `matlab.log`：实时KT进度和错误。

组合检查包含多个独立分支，模型输出目录只用于这些检查，不应当作一条完整科学试验轨迹。检查结论以JSON结果和日志为准。

## 通过之后还需要什么

本轮已完成：最终入口运行至168/168小时，退出码0。结果位于`docs/model_comparison/contract_runs/20260915-201910-501/`。有机碳、氮反应守恒残差最大分别为1.8319e-15和2.1075e-16 mg/g；完整状态第二小时断点续跑逐位一致，旧模式断点拒绝、零氧/部分供氧、错误候选输入不变、酶碳氮及各相合计检查通过，逐小时确认AOM保持为空。iHR=0/1源反应各8760小时回归和氮接口构造检查也通过。本次未修改正式模型方程或站点文件，修改的是组合验收入口。

短窗通过后，在固定源码和配置上做单独的全年连续运行，并检查更长时间的缺氧反馈、库存与运输预算。随后再恢复真实有机输入、AOM和凋落物，分别验证输入量、放置深度、时序及碳氮守恒。短窗或一次重启通过，都不替代这些工作。

接口规则见[中文接口说明](MEND_SOURCE_PRODUCTION_CONTRACT.md)，已有证据见[验收表](../../model_comparison/MEND_ACCEPTANCE.md)。

## 全年连续运行已启动

2026-09-15使用`startMendContractAnnual.ps1`启动独立连续运行，目录为`docs/model_comparison/annual_runs/compass-mend-contract-annual-20260915/`。运行使用已验收模型源码的原样副本；驱动仅暂缓有机肥事件和凋落物，逐小时检查AOM为空、微生物与酶相态、反应C/N、根预留、氧气预算及整柱有机运输闭合。初态保存在`initial.mat`，终态和完整预算仍需运行结束后审计。

新全年驱动先通过24小时预检，结果位于相邻的`compass-mend-contract-preflight-20260915/`目录，退出码0，C/N反应残差最大1.638e-15/1.912e-16 mg/g，氧预算残差最大8.257e-19 mol/cm²。

全年运行每24小时更新`progress.json`并打印`CONTRACT ANNUAL`进度；每168小时保存模型断点及`audit_KTxxxxx.mat`配套预算。续跑必须使用相同源码与配置，并提供同一KT的模型断点和预算记录，不能把已有累计量从零开始。只有`runner-status.json`显示`completed`且`annual-result.json`确认8760小时完成，才能进行最终全年判定。启动本身不代表全年通过。

在仓库根目录查看实时进度：

```powershell
Get-Content '.\docs\model_comparison\annual_runs\compass-mend-contract-annual-20260915\matlab.log' -Encoding UTF8 -Tail 20 -Wait
```
