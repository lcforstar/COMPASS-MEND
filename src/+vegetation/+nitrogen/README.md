# +vegetation/+nitrogen/

**Phase C M2-B.1b 新建包** (2026-05-19). 详 docs/cn_coupling_path_and_literature_audit.md §12.17.

## 物理范围

植物 N 池动态: organ allocation + LUNA-style leaf N → Vcmax 反馈 (M2-B.3) + 多年生季节储备 (M2-B.4).

**5 池 prognostic state** (cn_coupling §12.13.1):
- `N_labile`     mobile N pool (M2-A 已有, M2-B.1b 继续用)
- `N_leaf`       叶 N (M2-B.0 占位为 0, M2-B.1b 物理填)
- `N_root_fine`  细根 N (吸收主体)
- `N_stem_wood`  木质部 N (multi-year PFT 累积, annual=0)
- `N_fruit`      果实 N (收获带走, M2-B.5)
- `N_storage`    季节储备 (multi-year PFT, 春释/秋收, M2-B.4)

**Phase C M2-B.1b 范围** (本阶段):
- ✓ M2-B.1a: ebal 抽 `fluxes.Rd_leaf` → LPJ R_m,leaf
- ⏳ M2-B.1b: 组织 allocation (N_labile → 4 organ pools) + LPJ NPP 计算

**留待后续阶段**:
- M2-B.2: DAISY fixed_sink C_root 反算 (替代 M2-A C_min_dynamic)
- M2-B.3: LUNA 叶 N → Vcmax25_dyn 反馈
- M2-B.4: N_storage 季节储备 + retransl
- M2-B.5: N_fruit + harvest 事件

## 输入 / 输出

**输入** (parameters):
- 通过 `vegetation.getPFTValue(PFT_type, '*')` 从 PFT 表读:
  - C 分配比例: `f_alloc_leaf/root/stem/fruit`
  - 器官目标 C:N: `CN_L` (= CN_target_leaf), `CN_target_root/stem/fruit`
  - Sitch 2003 呼吸: `r_m10` (Table 3), `R_g_fraction` (p.169 = 0.25)

**输入** (forcing): `T_air, T_soil, fluxes.GPP, fluxes.Rd_leaf` (来自 ebal/SCOPE)

**输入** (state): `N_labile, N_leaf, N_root_fine, N_stem_wood, N_fruit`

**输出** (state updated):
- `N_labile, N_leaf, N_root_fine, N_stem_wood, N_fruit`

**输出** (diagnostic per KT):
- `NPP, R_m_leaf, R_m_root, R_m_stem, R_g, R_a_plant`
- `A_leaf, A_root, A_stem, A_fruit` (实际 allocation)
- `f_N_growth` (N 限制因子)

## 跟其他模块耦合点

| 耦合点 | 哪里来 | 哪里去 | 阶段 |
|---|---|---|---|
| 蒸腾流带 N | `Plantrootconductivity → Srt` | `+nitrogen/+root_uptake/nitrogenUptake` | M1 |
| 水胁迫调 demand | `hydraulic.supply_ratio` | `computeOrganDemand` (M2-B.2 用) | M2-B.2 |
| 叶 N → Vcmax | `computeLeafN_Vcmax (M2-B.3)` | `cfg.Veg.Photo_Vcmax_manual → SCOPE` | M2-B.3 |
| 新 LAI 反馈 | `N_leaf + CN_target_leaf` | `phenology` (M3) | M3 |

## ClimaLand 对齐 (§12.16b)

每个 .m 函数:
- 输入: `(state, params, forcing, dt)` 或子集
- 输出: `(state_new, diag)`
- 4 节 doc string: Arguments / Returns / Method / Reference
- 公式注释 specific 到方程号

未来 Julia 迁移: 直接对应 `prognostic_vars(::PlantNState)` + `compute_tendency!(dY, Y, p, t)` 风格.

## 文件清单

- `computeNPP.m`           LPJ Sitch 2003 完整 NPP (M2-B.1b)
- `computeOrganDemand.m`   D_N,i = G_pot,i / CN_target (M2-B.1b)
- `allocatePlantN.m`       N_labile → 4 organ pools 更新 (M2-B.1b)
- `computeLeafN_Vcmax.m`   LUNA Ali 2016 (M2-B.3 加)
- `seasonalStorage.m`      N_storage 春释/秋收 (M2-B.4)
- `harvestEvent.m`         N_fruit 收获 (M2-B.5)
