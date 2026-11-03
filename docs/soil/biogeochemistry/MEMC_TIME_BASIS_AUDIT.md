# MEMC 模型时钟审计

日期：2026-08-29

## 结论

COMPASS 中的 MEMC 反应时钟固定为 **day**，动力学速率统一使用 `d^-1`。
`advanceAcceptedCandidate.m` 把秒制子步除以 86400 后送入 MEMC 常微分方程是正确的；
当前 `soil.csv` 的日尺度参数不得再乘或除以 24。

原始 Wang et al. (2013) 论文中的速率单位是 `h^-1`。只有直接从该论文重新引入
参数值时，才需要先执行 `k_d = 24 k_h`，然后写入 COMPASS 的日尺度接口。

## 证据链

1. MEMC 官方 `R/model_component.R` 把 `memc_solve` 的输入定义为 daily time steps，
   输出时间单位为 days。
2. 官方当前 `data-raw/memc_params.csv` 把 `V_p`、`V_m`、`V_d`、`r_ep`、
   `r_em`、`K_ads`、`K_des` 和 `m_r` 全部标为 `d^-1`。
3. 上游提交 `3bf984e1e9ce8133a0f25a6722731264410f19a3`（2025-01-08，
   “correct the time label”）同时把求解图轴由 hour 改为 days、把参数标签由
   `h^-1` 改为 `d^-1`，没有对数值做 24 倍缩放。对当前 MEMC 软件包而言，
   这是明确的日尺度语义。
4. 上游提交 `63936c567bad881a04bad2fac35ec9b5f42b439c`（2025-12-08，
   “Canonical Wang 2013”）把 raw 参数更新为当前的 42、0.6、0.5168 等值。

## 上游仓库内部不一致

当前上游 README 和已提交的 `data/*.rda` 没有与 `data-raw/memc_params.csv`
同步：它们仍显示旧的 14、0.25、3.0 等数值；README 还残留 `h^-1` 文本，
而 `.rda` 已显示 `d^-1`。因此：

- 方程以 `R/model_component.R` 为准；
- 当前 intended raw 参数以 `data-raw/memc_params.csv` 为准；
- README 的旧单位和旧 `.rda` 不作为 COMPASS 参数来源；
- `mr_1d=0.01` 不是当前 raw 表默认值（raw 为 0），只能标记为待率定的
  非零工程先验，不能称为 MEMC 默认值。

## COMPASS 锁定措施

- 初始化写入 `rate_time_basis='day'` 和 `rate_unit_seconds=86400`；
- MEMC/MEND 反应器拒绝任何不是 86400 s 的速率时钟；
- 自动测试比较 24 个连续小时与一个完整模型日，要求状态一致；
- 当前时钟问题已关闭，但这不等同于参数已经率定。

## 主要来源

- MEMC source repository: <https://github.com/Microbial-Explicit-Model/MEMC>
- Current raw parameter table:
  <https://raw.githubusercontent.com/Microbial-Explicit-Model/MEMC/main/data-raw/memc_params.csv>
- Wang, G., Post, W. M. & Mayes, M. A. (2013),
  DOI: <https://doi.org/10.1890/12-0681.1>
