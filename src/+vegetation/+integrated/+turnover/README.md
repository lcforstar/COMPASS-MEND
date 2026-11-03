# +vegetation/+integrated/+turnover — N pool turnover + harvest + retransl + storage (M2-B.4)

## 物理问题

让植物 N 池有"出口" — M2-B.3 之前 N_leaf/N_fruit 只进不出, 物理不真实.

## 4 个 N 转移事件 (按 frPHU 触发)

```
Event            | 触发条件                      | 文献
-----------------|-------------------------------|----------------------------
(a) Apple harvest | frPHU > harvest_threshold     | EPIC (Williams 1989)
                  | (一次性, 每年最多 1 次)
(b) Leaf senesce  | frPHU > DLAI_threshold (0.45) | Aerts 1996, Niinemets 2010
                  | + LAI > 0.1
(c) Root turnover | 每日, 全年                    | Wells & Eissenstat 2002
(d) Spring release| 0 < frPHU < spring_release_endphu (0.20)
                  | (storage → labile)            | Tagliavini 2003, Aerts 1996
```

## 模块组织

```
+vegetation/+integrated/
  ├── turnoverAndHarvest.m         主调度 (caller-facing)
  └── +turnover/                    4 子函数 (私有)
      ├── appleFruitHarvest.m
      ├── leafSenescenceRetransl.m
      ├── fineRootTurnover.m
      └── springStorageRelease.m
```

## 输入 / 输出

**Inputs** (caller, COMPASS.m M2-B.4 hook):
- `pools_in.{N_labile, N_leaf, N_root_fine, N_stem_wood, N_fruit, N_storage}` (6 池)
- `state_in.{harvested_this_year, year_now}`  (跨 KT 状态)
- `pftParams.{harvest_threshold, harvest_fraction, root_turnover_rate, spring_release_*, leaf_long_days, retransl_fraction, is_perennial, has_fruit, has_storage_pool}`
- `env.{frPHU, LAI, DLAI, year}`
- `dt`

**Outputs**:
- `pools_out` (6 池更新, 守恒: Σ pools_out + Σ diag.* × dt × cum = Σ pools_in)
- `diag.{harvest_today, leaf_litter_today, root_litter_today, storage_release_today}` (g N/m²)
- `state_out` (写回 caller, 持续跨 KT/year)

## 跟外部耦合

```
COMPASS.m 主循环:
  ...
  plant N allocation (M2-B.1b → pools_new) ...
  ↓
  turnoverAndHarvest(pools_new, state, pftP, env, dt)
  ↓
  N_leaf 跌 → litter (cum_N_leaf_litter)   → M3 SOM 入口
  N_root  跌 → litter (cum_N_root_litter)  → M3 SOM 入口
  N_fruit 跌 → 出地 (cum_N_harvest)         → 不进 SOM (商业带走)
  N_storage ↑ (秋季) / N_labile ↑ (春萌)
```

## 守恒检查

整株 N 守恒方程 (每 KT):
```
Σ(N_labile, N_leaf, N_root_fine, N_stem_wood, N_fruit, N_storage)_KT
 + Σ(cum_N_harvest, cum_N_leaf_litter, cum_N_root_litter)_KT
 = Σ_init + Σ U_N × dt
```

## 参考文献

- Williams JR (1989) EPIC. *Trans.ASAE* 32:497 — harvest threshold
- Aerts R (1996) *Oikos* 67:469 — N resorption from senescing leaves
- Niinemets Ü (2010) — leaf longevity
- Wells CE, Eissenstat DM (2002) *JPGR* 21:324 — Apple fine root turnover
- Tagliavini M et al (2003) — Apple stem N as winter storage

## 实施进度

- [x] S1: PFT 表 [M] 5 字段 × 24 PFT
- [x] S2: 4 子函数 + turnoverAndHarvest 主调度
- [ ] S3: COMPASS.m 调度 hook + cum_* 累积 + N_storage prognostic
- [ ] S4: 9 诊断字段
- [ ] S5: 单元测试 + 7天 + 全年 + §12.25 存档
