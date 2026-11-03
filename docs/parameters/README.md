# 当前站点参数目录

Editor: Wenjie Liu; Date: 2026-10-31.

这些表重新读取迁移后的实际长武输入，保留源数值字面精度；历史目录仅复用已经审计的单位、作用和来源元数据，未单独审计的项目明确标注。

| 文件 | 数据行数 | 内容 |
|---|---:|---|
| `site_toml_parameters.csv` | 399 | 全部TOML叶键，含35个microbe键 |
| `soil_layer_parameters.csv` | 720 | 18个参数采样点×40列，含11个MEND列 |
| `MEND_parameters.csv` | 233 | 198条soil分层值与35条microbe值，共46个唯一参数 |
| `organic_initial_CN_allocation_58_layers.csv` | 58 | 未改变的初始C/N主池分配摘要 |
| `organic_initial_N_changes_22_layers_66_cells.csv` | 66 | 先前采用B_ref的22层×3个N库历史改动；不是本次吸附迁移的新改动 |

当前Qmax/Kba/Kdes_ref由soil拥有，TOML不重复定义；旧TOML标量配置接口保留。
18行是采样点，sample_linear映射到58个物理计算层。原始700/1000 cm行保留但不参与当前500 cm柱内插值权重。
Kdes_ref公共单位为mg C/g干土/h，既构造吸附系数也构造解吸容量速率。

来源可追踪不等于参数已有长武直接测量。旧全年结果和本次接口回归分开记录；这些CSV不代替运行守恒回执。

输入CSV现采用UTF-8 BOM，开头各保留一行与参数列逐一对应的中文含义和单位；原始注释与署名见 [输入表历史说明](INPUT_HEADER_HISTORY.md)。

当前soil分组顺序与TOML整理说明见 [输入文件排列](INPUT_LAYOUT.md)；TOML历史注释见 [归档](SITE_TOML_COMMENT_HISTORY.md)。

原作者通用注释的保留与当前实现适配见 [TOML注释来源](UPSTREAM_TOML_COMMENTS.md)。
