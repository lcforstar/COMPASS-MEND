# MEND逐层动态吸附参数接口（2026-10-31）

Editor: Wenjie Liu; Date: 2026-10-31.

当前接受基线仍为 `K1_mineral_0_20cm`。本次把三个动态吸附参考参数由
`site.toml [microbe]` 迁到 `soil.csv`，使其能按深度配置。18个原始参数行
均保留原参考数值；这次迁移没有为参数添加深度梯度，也不是重新率定。

| soil.csv键 | 18行当前值 | 公共单位 | 有效含义 |
|---|---:|---|---|
| `mend_qmax_mg_g` | 1.5 | mg C/g dry soil | 吸附碳容量Qmax |
| `mend_kba_g_mg` | 6 | g dry soil/mg C | 与Kdes_ref共同构造吸附系数的亲和参数 |
| `mend_kdes_mg_g_h` | 0.006 | mg C/g dry soil/h | 20°C解吸参考容量速率，同时参与吸附系数构造 |

每个键只能有一个所有者：soil逐层列或旧版TOML标量。当前长武站选择soil所有权，
TOML中不再保留这三键；同时给两处赋值应报重复所有权错误。旧版TOML标量接口保留，
不得以接口扩展为由强迫旧 `rk4_equilibrium` 站点新增不生效的动力学吸附参数。

`mend_vertical_parameter_mode="sample_linear"` 把18个原始采样点映射至当前58个物理计算层，
然后按微生物模块的顶向下顺序使用；原始参数行不是18个常值土层。
raw700/1000 cm行保留输入，但在当前500 cm柱内不贡献插值权重。
Qmax必须大于0，Kba和Kdes_ref不得小于0；配置值还须为有限实数。

当前 `source_hourly` / `mend_dynamic` 路径逐层计算：

```text
Arr(T,Ea) = exp[Ea*1000/R * (1/293.15 - 1/(T+273.15))]
Kads(T) = Kba * Kdes_ref * Arr(T, 5 kJ/mol)
Kdes(T) = Kdes_ref * Arr(T,20 kJ/mol)
ads_C = Kads(T) * D_remaining_C * (1-Q_C/Qmax)
des_C = Kdes(T) * Q_C/Qmax
```

D_remaining_C为当步DOM摄取后的可用碳，Q_C为反应前吸附碳；两者基准单位为mg C/g干土。
上述为未施加原生供体上限前的参考通量（mg C/g干土/h）；实际转移保留原来的供体限制和时间映射。
Kads为h^-1；Kdes_ref及Kdes是容量速率，不能只标成h^-1。
DON交换沿用既有配对C/N路由与真实供体组成，不额外引入独立的DON吸附参数。
`linear_equilibrium` 模式不再叠加该动力学吸附/解吸。

当前 `soil.csv` 为18行×40列，含11个MEND列；`site.toml [microbe]` 留35键，
全部TOML为399叶键。MEND仍是46个唯一配置键，文件位置改变没有增加新的率定轴。
参数目录列出198条soil深度值和35条TOML值，共233行。

初始有机C/N仍由 `organic_initial_pools.csv` 单独拥有；本次不改58层×70列初值、
施肥日期和氮量、0–20 cm入土配置、FAST10/TAN先验、水力参数或守恒阈值。
Qmax/Kba/Kdes的既有值源于上游示例/先验沿用，不能因此声称是长武直接测量。

此前8760步结果是已采用的历史基线证据。本次接口迁移须另核短步、逐层与旧标量
兼容性及同有效参数回归；是否通过以本轮实际验证回执为准，不宣称重新完成全年。
原优化器保持暂停；本次不自动启动率定。

相关实现：`mendLayerParameterSpec` → `resolveMendLayerParameters` →
`buildSoilMicrobes` / `buildMendHourlyState` → `advanceSourceCandidate` → `mendSourceHour`。

## 本次验证回执

Windows MATLAB R2025a：旧/新代码三种配置各4小时的接受态、29通量、原生预算及诊断数值严格字节一致；异质分层、所有权、物理域、零速率极限、旧RK4缺三键初始化和独立桌面版本4步检查均PASS。原C/N/DIN/O2门槛保持。未重新运行全年。新soil向量配置不静默接管旧标量hourly检查点。
