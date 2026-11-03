# MEND 分解与 DOC 摄取参数层化合同

编辑：Wenjie Liu；2026-09-29。

> 2026-10-12 接口更新：本文保留四项参数的历史迁移与逐位验证记录；当前正式接口、
> 8+2 候选清单、独立 MEND 层组映射及释放顺序见
> `CHANGWU_MEND_LAYER_CANDIDATES_20261012.md`。新增接口不表示十项参数已同时进入优化。
>
> 2026-10-11 长武配置记录：四项动力学在所有采样层原本完全相同，且没有
> 分层微生物动力学观测支持。它们已恢复为 `site.toml [microbe]` 的全站标量，
> 数值和方程不变。本文其余内容保留为逐层输入能力及2026-09-29历史迁移记录；
> `resolveMendLayerParameters` 仍支持有独立分层证据的站点使用完整四列配置。

## 范围与唯一输入来源

本次只迁移四项碳动力学参数的输入所有权和数组形状，不改变先验数值、方程、
初态、施肥、边界、时间步或守恒阈值。没有启动率定。

| soil.csv 列 | 原先验（各行相同） | 单位 | 作用 |
|---|---:|---|---|
| `mend_vd_pom_h` | 50 | h⁻¹ | 共用的 POM1/POM2/MOM 参考分解容量 |
| `mend_ks_pom1_mg_g` | 60 | mg C/g 干土 | POM1 参考半饱和常数 |
| `mend_vg_h` | 0.00425194 | h⁻¹ | 微生物 DOC 生长摄取尺度 |
| `mend_ks_dom_mg_g` | 0.00010034 | mg C/g 干土 | DOC 摄取参考半饱和常数 |

这些数值从迁移前 Changwu `site.toml` 原样移动，不声称是长武最优值，也不是由
文献浓度单位重新换算出来的另一套参数。新列追加到原表末尾，原有 29 列保持。

固定参数仍在 site，例如 `mend_ks_pom_ratio`、CUE、alpha、休眠响应、酶分配。
微生物 DIN 摄取参数本次不层化。POC/MOC/微生物初始分库仍在原初态文件。

同一组输入必须满足下列之一：

1. soil 包含全部四列，site 不再包含这四个键；这是当前 Changwu 模式。
2. soil 不包含任何一列，site 包含全部四个正有限标量；仅供未迁移旧站点兼容。

部分列、NaN/Inf、零/负值、两处同时定义（即使值相同）均报错。不存在静默
回退或 soil/site 优先级覆盖。公共列名与单位由 `io.adapters.mendLayerParameterSpec`
列出，`resolveMendLayerParameters` 同时服务检查器、初始化和采样器。

## 深度、形状和单位

`soil.csv` 的深度含义沿用站点的 `soil_input.vertical_parameter_mode`：

- `sample_linear`：各行是采样坐标，线性插值到计算单元中心，范围外保持端点值。
- `horizon_bounds`：各行是层底边界，按既有分层规则映射。

本次不改变该选择。单独改一行采样点时，sample_linear 模式下相邻插值区间
也会变化；这与“直接修改一个已映射计算层参数”的单层测试不同。

`buildSoilParameters` 给出 bottom-up 的 NL×1 参数，不对这四项做日/秒换算。
`buildSoilMicrobes` 在明确的 MEND 适配边界翻转为 top-down。POM 双分量
使用 NL×2（列 1=POM1，列 2=POM2），不是把两个层误当两种底物。
环境响应只在逐层参数上施加原温度、水势和 pH 因子；旧标量在环境适配器中广播。

以下派生关系保持不变并按层计算：

```text
Vd_POM1 = Vd_POM2 = Vd_MOM = mend_vd_pom_h
Ks_POM2 = Ks_POM1 / mend_ks_pom_ratio
Ks_MOM  = Ks_POM1 * mend_ks_pom_ratio
Vm      = Vg * alpha / (1 - alpha)
rMORT   = min(0.99, Vm * mortality_multiplier)
VmD     = Vm * dormant_maintenance_ratio
```

因此本次允许共用容量随深度变化，**没有增加独立 MOC 分解参数**，也没有将
派生的维持/死亡速率作为新的独立优化变量。

## 采样及重启

`changwuCarbonCalibrationSpec` 和 `materializeChangwuCarbonSample` 读取真实输入
所有者。旧的 `mend_vd_pom_h` 优化坐标现在表示第一个土壤采样行的参考值，
通过 `候选值/原首行值` 缩放整条剖面；不把已有非均匀剖面压平。
`mend_ks_dom_scale` 缩放 soil 中的 KsDOM。基准样本比例为 1，不改变原剖面。
本次未增加 BO 维度，正式深度分组和参数范围仍需后续设计。

新 checkpoint 签名记录规范化的 top-down 动力学剖面及相应派生量：改变某层
参数后不能接着用原参数的 checkpoint。旧 checkpoint 没有这份证据，会拒绝
自动续跑；它们保留作历史结果，不伪造新签名。未迁移旧站点可以重新初始化运行。

原水氮 campaign 的代码/站点身份亦已改变，不能追加到旧 campaign。已完成的
44 个结果不覆盖、不重标记；迁移前源码、站点和率定工具保存在
`outputs/mend_layer_parameters_20260929/pre_migration/`。

## 验证依据

`tools/diagnostics/verifyMendLayerParameters.m` 的 capture 阶段在修改前运行 168 h，
保存初态、环境响应、KT25 和 KT168 接受态；verify 阶段执行以下检查：

1. 旧 soil 列和非迁移 TOML 值不变，迁移后的常数剖面与旧环境响应一致。
2. 旧标量输入、两种深度映射、单一采样深度和乱序列名。
3. 重复所有者、缺列、非法值必须失败。
4. 分别扰动四个参数的一个已映射层，检查真实 source-hour C/N 反应候选只在
   该层变化，排除错误广播/层序反转；测试不包含本来会跨层传播的输运。
5. 原参数 168 h 的逐步 C/N 断言、初态/KT25/KT168 与迁移前基准逐位比较。
6. KT24 保存后官方重启至 KT25，与不间断结果逐位比较；更改参数或缺少
   动力学签名的 checkpoint 必须拒绝。
7. 现有碳采样器的基准和剖面缩放。

实际结果以 `outputs/mend_layer_parameters_20260929/verification.json` 为准。
测试通过不等于这四项参数已率定，也不替代初态/浸提口径的科学核查。

本次实际运行结果：168 h 对照通过，初态、KT25、KT168 的完整接受态逐位相同；
官方重启 KT24→25 也逐位相同。最大逐小时反应 C/N 残差为
`1.7208456881689926e-15 / 2.055969059203575e-16 mg/g`。
`supplemental.json` 记录关闭微生物时的签名、NL=1/17 和错误层序拒绝检查通过。
`input_migration_audit.json` 验证 18 行原有 29 列和全部非迁移 TOML 值未改；
包括观测工作簿、初态在内的其他站点文件逐字节相同。原有率定单位/解析/聚合/
合成敏感性检查也为 0 失败，日志为 `supplemental.log`（其中 failed sample
属于测试自建的合成夹具，不是新生产运行失败）。

已知无关问题：现行 TOML 编辑器 schema 中 microbe 节仍保留旧 MEMC 定义；
它不是生产 `io.checkSite` 的验证来源，本次未借层化任务重写整个 schema。
