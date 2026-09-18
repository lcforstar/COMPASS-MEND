# COMPASS 求解器调试现状与待办

> 维护：朱寄子星 · 最后整理：2026-05-13（v33 log-SOR + Tol_Water=50 修复）
> 老版本（按时间堆叠的工作日志）已被本版替换，git 历史保留原文。
> 状态：历史试验档案。旧 Fourier 梯度混入 ebal 的支路及配置入口已于 2026-07-18 删除；当前权威热接口见 `surface_soil_thermal_coupling_progress.md`。

本文档围绕三个用途组织：
1. **现状（§1-§3）**——新会话或新人接手时第一时间需要看的内容。
2. **未解决 / 待办（§4-§5）**——下一步该做什么。
3. **历史与档案（§7-§9）**——为什么改成今天这样、改过哪几轮。

---

## §1. 当前生产配置

**模型主跑路径**：`COMPASS.m → ebal_clm45.m → resistances_clm45.m`（CLM4.5 风格）。
**配置时段**：长武苹果塬 2021-01-01 ~ 2022-01-01，8760 KT，NL=128（动态调整 NL=200）。

| 文件 | 关键参数 | 含义 |
|---|---|---|
| `data/config_environment.txt` | `Z_measure=10`, `Wind_Height_Convert=1`, `Wind_Obs_Height=2.0` | U2→U10 对数廓线转换 |
| `data/config_stemmus_scope.txt` | `lFreezeLatentCoupling=1` (V8 Fourier-G), `lFreezing=1`, `lSFCC=1`, `lSWCC=1`, `lThmrlefc=1` | 冻融潜热 V8 注入 + 严格冻融 |
| `data/config_vegetation.txt` | `hydraulic_p50=-3.46`, `hydraulic_aa=3.5`, `hydraulic_gpmax=2e-08`, `root_contact_enabled=0`, `root_dist_mode=0` (Jackson) | 植物水力 Case D + Jackson 根分布 |
| `data/InitialConditions.xlsx` | 当前=10/2 状态（KT=36 复现现场） | 1/1 备份在 `InitialConditions_jan01_backup.xlsx`，10/2 在 `InitialConditions_oct02_for_step6612.xlsx` |
| `src/+soil/+water/solveTridiagonalMatrixEquations.m:48-62` | 对数空间自适应 SOR (ω=0.5 触发阈 CHKWATER>200cm) | v33 修 KT=36 数值病态 |
| `src/COMPASS.m:1339-1346` | `Tol_Water=50` cm + `STEMMUS_TOL_WATER` 环境变量 | v33 严收敛容差 |

**Richards 边界**：`NBCh_top=3` (atmospheric), `NBCh_bottom=3` (free drainage), `DSTMAX=0` (no ponding), `hcrit=-15000 cm`.

**v33 全年残差**: −39.4 mm = **3.9% of P**（修复前 v32: 8.0%，OK 阈值 <5%）。详 `docs/figures/v33_vs_v32_comparison.md`。

**主要诊断输出**（`data/output_full_year/`）：
- 物理：SoilMoisture / SoilMatricHead / SoilTemperature / Throughfall / TopWaterInput / Runoff
- 能量：Rntot / Htot / lEtot / Gtot / lE_CanopyIntercept_W / fwet / h2ocan_mm
- 根系：RWU / RWU_pos / RWU_neg / HydraulicLift_intensity

---

## §2. 当前未解决偏差（按物理重要性）

### a. 100-200 cm 暴雨期不响应（最大物理偏差，详 §5 / B-2）
2021-07-15 暴雨期 v16 实测 100 cm Δθ=0.008 / 200 cm Δθ=0.006，实测 0.150 / 0.162。**水流根本没传到这个深度**——dual-porosity 已确认救不了，需要 dual-permeability 或显式 macropore。

### b. 500/700 cm 深层水势相位异常（详 §4 / U7）
v25 全年实测 KGE：500 cm = -0.04（PBIAS +29.9% 偏湿），700 cm = -0.56（反相位）。可能跟 hydraulic lift 简化版（B-1）+ 剖面只到 1000 cm 缺深层水库（U7）相关。

### c. 冬季 H 偏低 -21 W/m²（详 §5 / B-4）
ebal 残差里缺 `S_freeze`（冻融潜热源项）和 `S_storage`（表层储能项）。V8 Fourier-G 注入只解决"G 量级"，没修源项缺失。

### d. 8/9 月叶水势季节偏差（论文残差，case scan 不能解）
- 9 月 midday 偏不负 0.4 MPa（秋末水分胁迫低估）
- 8 月 predawn 偏负 0.27 MPa（夜间恢复不足，疑似 B-1 hydraulic lift 缺位）

---

## §3. 与 STEMMUS-SCOPE 公开版的差异（合并对比表）

参考公开版：`F:/STEMMUS_SCOPE-main/`（含冻融）+ `F:/STEMMUS_Zhujizixing -马铃薯/`（朴素无冻融）。

### 3.1 求解器与数值

| 维度 | 公开版 | COMPASS |
|---|---|---|
| Picard 迭代 | 朴素 | 带 Anderson 加速 |
| 子步策略 | 固定整步 | 自适应（dt 缩放） |
| 步失败处理 | 整步重试 | StepFailed(KT) 标记继续 |
| 输出聚合 | 单步直接 | 子步加权平均到 KT |
| 冻融路径 | 无 / lSFCC=0 | lSFCC=1 严格（hh+hh_frez） |
| **Thomas SOR** | omega=1（无阻尼） | **对数空间 SOR（CHKWATER>200cm 时 ω=0.5，hh 几何混合）** ⭐ v33 |
| **Tol_Water** | 1 cm（朴素严格） | **50 cm（默认）+ STEMMUS_TOL_WATER 环境变量覆盖** ⭐ v33 |

⭐ v33 修复（详 §9 已解决档案 §9.X）：
- 公开版大雨步 KT=36 Picard 第一次解出 hh 残差 **3926 cm（39 米水头）** 是数值病态而非物理。
- 对数空间 SOR + Tol=50 联合修复后：**KT=36 首次 Water=528 cm（87% ↓）**，72h 残差 8.5% → 0.28%，全年残差 8.5% → 3.9%。
- 0 FATAL，冻融期 (1-2 月) 完全跑通。

### 3.2 物理过程

| 物理过程 | 公开版 | COMPASS | 说明 |
|---|---|---|---|
| **截留模型** | 无（Throughfall=Precip） | Gash 1979 + h2ocan 状态变量 | v29 加入 |
| **fwet 干湿叶分离** | 无（lE_old 公式） | CLM5 风格干湿拆分 | v29 加入；阴阳叶共用单 fwet（待 B-5 升级） |
| **截留蒸发潜热 lE_CI** | 0（隐式漏算）| `lE_wet`（峰 96 W/m²，年均 ~20）| v29 修复，详 §9 已解决档案 |
| **SPAC 蒸腾输入** | `lEctot` 含截留 ❌ | `lE_dry` 仅气孔蒸腾 ✅ | v29 修复 |
| **u 下限钳位** | 无 | MeteoData 入口 `max(1.0, U)` | v28 对齐 CLM5/NoahMP |
| **ustar 公式** | W&V Eq 30 + Paulson 1970 (2 区间) | resistances_clm45 (CLM4.5 4 区间) | v28 统一口径 |
| **RWU 通量含正负** | line 545 含正负，line 553 钳 1e-20 | 拆 RWU_pos/RWU_neg，仅 net 进 sink | v27 |
| **Σ RWU = Tp_t** | 否（SPAC 自然闭合） | 否（discount 自动闭合） | — |
| **Hydraulic redistribution** | 半允许（数值不稳） | 简化版（不真改 θ） | 待 B-1 升级 |
| **根分布** | empirical | Jackson 1996（0-50cm 60%） | v27 |
| **GPP 单位** | kg C/m²/s | 同 | — |
| **能量闭合诊断** | 内置 | 后处理 | — |

### 3.3 Switches 行为差异（不是 bug，仅记录）

| 维度 | 当前 COMPASS | 朴素 STEMMUS |
|---|---|---|
| 冻融 | lSFCC=1 严格路径 | 无 |
| 双孔 | lDualPor=0（默认） | 无 |
| 滞回 | lHystrs=0 / 1 仅 5% K 修正 | 无 |
| 大气边界 | NBCh_top=3 | 同 |

### 3.4 已确认非 bug 的细节

| 项 | 状态 | 说明 |
|---|---|---|
| `calculateMatricCoefficients` 差分守卫 4 条件 | 与公开版一致 | `+soilmoisture/calculateMatricCoefficients.m:43-48` |
| `HT_frez` 在 Picard 内用旧 TT | 与公开版一致 | `STEMMUS_SCOPE.m:625-631` |
| `hh_frez` 加进水方程梯度 DhUU | 与公开版一致 | `+soilmoisture/calculateMatricCoefficients.m:40-41` |
| `Gamma_hh=0 → DTheta_LLh=0` | 与公开版一致；暴雨期不触发 | — |
| `DTheta_UUh/LLh` 标签 | **公开版有 bug，COMPASS 反而是对的** | 不要往回改 |

---

## §4. 已知未解决问题（U 系列）

### U2. Theta_LL 不能超过孔隙度 ⭐⭐⭐
COMPASS 状态范式（hh 主控）下 `Theta_LL ≤ θ_s` 是数学约束。要允许 `ϑ_l > ν` 必须切到 ClimaLand 的"ϑ_l 主控"范式，等于改造整个状态变量体系（Route B，预计 1-2 周）。
**决策**：当前不动，配套留 Route B 阶段实施。详见 §6 Limitation。

### U5. SUPG 人工弥散两难 ⭐⭐⭐
- **副作用**：抬升干层 Khh 至 1e-8（远大于 VG 真实 K(h)<1e-11）→ 形成 50 cm 抽水链瓶颈。
- **必需性**：2026-05-08 实测删 SUPG 后 KT=1 矩阵奇异（RCOND=2.92e-17，Water 残差 9.25e+6），dt 砍到 112.5 s 也跑不动第一步。
- **修复方向**：条件 SUPG，仅在 `Pe>1` 或 `hh<-1e4` 时启用。3 处共 ~10 行代码。

```matlab
% 仅在 Pe>1 时启用
Pe_local = abs(Vvh(i,j)) * DeltZ(i) / max(Khh(i,j), 1e-15);
if Pe_local > 1
    Khh(i,j) = max(Khh(i,j), abs(Vvh(i,j)) * DeltZ(i) / 2);
end
```

### U7. 剖面深度仅 1000 cm（缺深层水库）⭐⭐
Shao 2024 WRR 长武苹果 25 龄根系达 2000 cm，深层 200-2000 cm 贡献年 ET 22.9%（极旱年 26%）。当前 `Tot_Depth=1000`, `root_max_depth=800`，缺 1000-2000 cm 段水库。
**修复**：扩展 grid 到 2000 cm + 加 IC 表深层数据 + 更新 root_max_depth=2000。

### U10. 顶层 N1 hh 撞 hcrit 触发阈值过敏感 ⭐
- 旧跑 6/30~7/4 N1 最低 hh = -28706 ~ -33804 cm（远低于 hcrit=-15000）
- 根因：`solveSoilMoistureBalance.m:72` `isRaining = Throughfall > 1e-10` 即 3.6e-6 mm/h，微量穿透雨触发 isRaining=true → 钉死被旁路
- 修复：阈值改成 `> 0.1/36000`（0.1 mm/h）
- 当前不改原因：v25 hd=-3e4 钉死会先于 hcrit 触发，hd 接管后 U10 不再是问题。**待全年验证 hd 是否真接管完，再决定是否补 U10 修复**。

### U12 / B-2 共因。详 §5 / B-2。

### U13. N 模块物理修正路线（2026-05-13 启动）⭐⭐⭐
完整路线见 [docs/n_module_refactor_progress.md](n_module_refactor_progress.md)。

阶段 0 完成：N 模块按 DAISY 风格拆 5 子包 (`+chemistry/+organic/+transport/+root_uptake/+fertilizer`)，删除 8 个冲突/dead 文件，N 守恒检查改成 ClimaLand-style 累积通量。

阶段 1 待执行（按顺序）：
1. **P0a** unfrozen water 替代 — N 输运用未冻液态水，不用总 θ_LL（冬季浓缩效应）
2. **P0b** O2 运移冻土期 fallback（避免 Theta_g≈0 数值病态）
3. **P1**  N 沉降从 Excel 时间序列读（顶部 BC 时间变化）
4. **P1.5** 施肥开关单测（守恒账闭合验证）
5. **P2**  1 个月测速跑
6. **P3**  全年跑 + N 守恒 < 5%
7. **P4**  详细规范化重写（Skill §2 + ClimaLand 风格）

---

## §5. 待办路线图（B 系列，按优先级）

### B-1. CLM5 完整版 hydraulic redistribution ⭐⭐⭐ 最高优先级
**问题**：v27 简化版只让湿层 `RWU_pos × discount` 等效"少吸"，**不动 θ 做层间搬运**。夏季夜间 / 冬季落叶期 `Tp_t=0` 时 sink=0 直接放过 hydraulic lift（实测 v27 KT=8 凌晨 7 点 `F_pos=1.16, F_neg=-1.16, F_net=0`，整夜没水移动）。**全年累积可能损失 70+ mm 水力 lift**。

**实测疑似关联**：500 cm PBIAS +29.9% 偏湿（深层应被 lift 抽水补浅层，模型没做）。

**实现路径**（参考 `D:/OneDrive/CTSM-release-clm5.0.37/src/biogeophys/SoilWaterPlantSinkMod.F90:305-336`）：
1. v27 已 expose `RWU_pos_node / RWU_neg_node` 到 `hydraulic` 结构体
2. 新建 `applyHydraulicRedistribution.m`：post-Richards 把 |F_neg| 按 |RWU_neg| 加权分配回 Theta_LL（湿层 -= / 干层 +=），Σ Δθ × DZ × dt × 10 = 0
3. 数值约束：干层加水 ≤ θ_s，湿层减水 ≥ θ_r

**工作量**：1-2 天（搞清主循环 + 写函数 + 插入 + 单元测试）。

**期望**：浅层夜间 θ +0.005~0.02、ψ_leaf 跌落缓冲、500 cm 偏湿改善、年 ET/P 不变。

---

### B-2. Macropore / 优先流（解决 100-200 cm 不响应）⭐⭐⭐ 第二优先级

**实测对照**（2021-07-15 暴雨期 Δθ）：

| 深度 | 实测 Δθ | v16 (dual-pora) Δθ | 当前/实测 |
|---|---|---|---|
| 50 cm | 0.180 | 0.070 | 39% |
| 70 cm | 0.216 | 0.046 | 21% |
| **100 cm** | **0.150** | **0.008** | **5%** ❌ |
| **200 cm** | **0.162** | **0.006** | **4%** ❌ |

**已确认是结构性物理缺失，不是 bug**：
- v16 hh 时间序列：72h 暴雨期 100 cm hh 仅波动 -470→-453 cm（约 40 cm 水头）。**水根本没传到此深度**
- 表层 K_s = 3e-4 cm/s，但 100 cm K(h) 暴雨峰仅 1e-5 cm/s，平均水穿透 100 cm 需 116 天
- 实测 50/70/100/200 cm 同步响应（峰值都在 7/16）= **优先流通道直接把水送下去**
- HYDRUS dual-porosity (iDualPor=1) 不是 macropore：只解 mobile Richards，immobile 靠 SinkIm 被动接收，没有大空隙独立 K(h) 快速通道

**升级路线**：

| 方案 | 工程量 | 适用场景 | 优先级 |
|---|---|---|---|
| A: dual-permeability 顺序求解 | 1-2 天，~500 行 | 基质-裂缝**连续**双介质 | ⭐⭐⭐ 起步首选 |
| B: dual-permeability 块三对角耦合 | 1 周，~1500 行 | 同 A，更稳但工程量大 | ⭐ A 不够再升 |
| C: Daisy biopore (MIM/kinematic wave) | 2-3 周 | **离散根孔/裂缝**优先流 | ⭐⭐ 看响应形状 |

**A vs C 判别**（依据 100 cm 暴雨期小时分辨率响应形状）：
- 平滑爬升（几小时到一天）→ 走 A
- 突然脉冲（暴雨当下出现，50-70 cm 还没湿）→ 走 C

黄土塬 + 10 龄苹果根穿孔/裂缝结构 → **C 更可能正确**。Daisy §4.7 几何参数：管半径 r_b、间距 D_b、深度 z_top/z_bot；Poiseuille 重力流 Eq. 4.34；基质↔大孔交换 Eq. 4.36。

**论文价值**：黄土塬苹果园 macropore 显式建模 + STEMMUS-SCOPE 框架，论文核心创新点之一。

---

### B-3. 双孔 (lDualPor=1) 冬季冻融失稳 ⭐⭐
**功能性 bug**（开关打不开），不是物理偏差。当前 lDualPor=0 也能跑，但缺溶质运移能力。

**根因（2026-05-08 v23 实测）**：v16 双孔修复版 `dualPorIni.m` 把 `VG.Theta_s` 砍掉 `ts_im`，冬季 Theta_LL 离 `θ_s_eff` 太近，加上冰相 `Theta_II × ρI/ρL` 撑爆容量，VG 雅可比 dθ/dh 在 Se→1 奇异 → Water 残差爆 1.5e+10。

**修复清单**：
1. ✅ `dualPorIni.m` HYDRUS InitW 标准（v16 已做）
2. 🔴 **待修**：`Theta_V = POR - Theta_LL - Theta_II×ρI/ρL - Theta_im`（冰占 mobile，不重复扣 immobile 部分）
3. 🟡 **待验证**：诊断断言 `Theta_LL + Theta_II*0.917 ≤ Theta_s_m_eff + 1e-6`

**关键物理共识**（HYDRUS WATFLOW.f line 1037 一致）：冰只在 mobile 域形成，immobile 不冻。

**注意**：B-3 修好也救不了 100-200 cm 不响应（dual-porosity 仍受 K(h) 卡死，要 B-2 才行）。

---

### B-4. ebal 加 S_freeze + S_storage 源项 ⭐⭐
解决冬季 H -21 W/m² 偏低的根本。

**当前 ebal**：`Rn = LE + H + G`（经验闭合）
**CLM5/NoahMP 完整版**：`Rn = LE + H + G + S_freeze + S_storage`

**S_freeze 冻融潜热**：
- ebal 不知道 `L_f × ∂θ_ice/∂t`（土壤热扩散方程算了但没回灌 ebal 残差）
- V8 Fourier-G 注入只修"G 量级"，不修源项缺失
- 真正修法：在 EBers 残差里加 `S_freeze = L_f × ρ_I × ΔTheta_II / dt`，~100 行大改
- 量级：1 月 ~30 W/m²，年均 ~5

**S_storage 表层储能**：
- ebal 假定 `∂(c_v·T)/∂t = 0`，但表层 0-N cm 升降温有储能
- CLM5: `c_v(0) × dT/dt × dz` 入 EBers
- 量级：小时 ±50，年均 ~0

**其他能量平衡缺失项（与 CLM5/NoahMP 对比，按优先级）**：

| 缺失项 | 量级 | CLM5 名 | NoahMP 名 | 优先级 |
|---|---|---|---|---|
| ~~lE_CI 截留蒸发潜热~~ | ~~~20 W/m² 年均~~ | `fwet × eflx_lh_vege` | FCEV | ✅ v29 已修 |
| **S_freeze 冰相变** | 1 月 ~30, 年均 ~5 | `xmf` | (内嵌) | ⭐⭐ |
| **S_storage 表层储能** | 小时 ±50, 年均 ~0 | `Σ Δt_soisno_j/fact_j` | (内嵌) | ⭐ |
| CANHS 冠层显热储热 | 小时 ±50, 年均 ~0 | (内嵌 cgrnd) | CANHS | ⭐ |
| PAH 降水带来的热 | ~0-2 | (无) | PAH | ⭐ |
| eflx_snomelt 雪融 | 长武≈0 | `eflx_snomelt` | QMELT | — |

---

### B-5. 干湿叶 × 阴阳叶四象限混合 ⭐⭐⭐ 论文核心创新点

**起源**：v29 fwet 实施后发现 CLM5/NoahMP/SCOPE 公开版**都用单一 fwet**（全冠层平均），不分阴阳叶。我们提出**四象限混合**作为新物理。

**单一 fwet 的物理 bug**：暴雨刚停，阳叶顶层接受了大部分雨被全湿（fwet_sun=1），阴叶底层几乎没湿（fwet_sha=0.2），但模型用 fwet_avg=0.6 当成两者都是 0.6 → 阳叶蒸发被低估 / 阴叶蒸发被高估。

**四象限独立 fwet**：

|  | fwet_dry | fwet_wet |
|---|---|---|
| 阳叶 sunlit | 1-fwet_sun | fwet_sun |
| 阴叶 shaded | 1-fwet_sha | fwet_sha |

4 个 lE 通量分别算：`lE_sun_dry / lE_sun_wet / lE_sha_dry / lE_sha_wet`。

**雨水 routing（顶部冠层先接收，Beer-Lambert）**：
```matlab
P_to_sun = (1 - exp(-k_rain * LAI*Ps)) * Precip
P_to_sha = Precip - P_to_sun
```

**为什么是创新**：

| 方面 | CLM5/NoahMP | SCOPE 公开版 | **v30 (我们)** |
|---|---|---|---|
| 阴阳叶辐射 | 简化（big-leaf） | ✅ 严格分 | ✅ 严格分 |
| 干湿叶分离 | ✅ 全冠层 fwet | ❌ | ✅ 阴阳叶独立 fwet |
| 雨水冠层 routing | 单层 | 不算 | **顶→底 Beer-Lambert** |
| 阳叶 vs 阴叶 lE | 平均 | 分阴阳 | **4 象限独立** |

**实施步骤**（3-5 天）：
1. 状态变量加倍 `state_h2ocan_sun` + `state_h2ocan_sha`
2. 雨水 routing：`calculateRainfallPartitioning` 输出 `P_to_sun / P_to_sha`
3. fwet 拆双：`fwet_sun = (h2ocan_sun/max_sun)^(2/3)`，阴叶同
4. `interception` 结构体扩展两路
5. heatfluxes_clm45 调用：阳叶传 fwet_sun，阴叶传 fwet_sha
6. 守恒检查双路

**潜在数值挑战**：`Ps`（sunlit fraction）随太阳角变化 → max_storage 随 KT 变 → 需要白天/夜间水库再分配逻辑。

**文献依据**：Choudhury & Monteith 1988（fwet 二层分）、Brutsaert 1982（surface wetting heterogeneity）、Sellers et al 1986（SiB wet/dry sunlit/shaded）。

---

### B-6. O2 模块完善 — R_O2_sink 补根呼吸 + V_A 瞬态检查 ⭐⭐ N 论文配套

**起源**（2026-05-16 N 量级排查）：发现 N 浓度量级偏低、硝化年通量偏高时，对 O2 模块做了一次审计。模块核心算法（Galerkin FEM + Millington-Quirk 扩散 + Henry buffer + Michaelis-Menten f_aerobic）写得很好，但有 3 处可改：

**a. R_O2_sink 漏根呼吸**（`src/+soil/+gas/calcO2Sink.m`）
- 当前：`R_O2 = R_O2_nit (硝化) + R_O2_min (异养矿化)`
- 真实土壤 O2 sink 中**根自养呼吸占 50-70%**（生长季），微生物异养占 20-40%，硝化只占 5-10%
- **后果**：R_O2_sink 系统性偏小 → C_O2 偏高 → f_aerobic 偏高 → 硝化加速 → NH4 池被掏空更快
- **修法**：从 `carbon_fluxes.nc` 取 `R_root` (autotrophic root respiration) → 按 O2:CO2 = 1:1 化学计量加进 R_O2_sink
  ```matlab
  R_O2_root = root_respiration_per_layer * 1e-6 / 12 * 1;  % g C → mol O2
  R_O2_sink = R_O2_nit + R_O2_min + R_O2_root;
  ```

**b. V_A（干空气 Darcy 速度）瞬态过程缺失**
- 当前驱动 O2 入土对流的 V_A 由 STEMMUS dry-air 方程给（lSoilairefc=1），量级 < 0.001 cm/s
- 但实际土壤 O2 大量交换在**瞬态过程**：
  - 雨水入渗的"活塞效应"（雨水挤走表层空气）
  - 夜间土壤冷却 → 气体收缩 → 大气 O2 被吸入（thermal pumping）
  - 风压扰动 (wind pumping)
- **后果**：深层 O2 补给慢 → 根区 f_aerobic 偏低 → 根区硝化偏慢（**跟 a 的方向相反，可能互相抵消**）
- **首次检查**：先输出 V_A 全年时序，看雨后/温变是否有脉冲信号；如果 V_A 全年≈恒值，说明瞬态没被捕获，再补
- 修法可能要在 calculateGasConductivity 加 source term，或在 solveO2Transport 加显式雨脉冲项

**c. 底边界 Neumann zero-flux 对 1000 cm 深剖面是否合理**
- 当前 `use_bot_dirichlet=false` → 深 1000 cm 处零通量
- 长武果园地下水位 ~50m，物理上 1000-5000 cm 段还是 vadose zone，O2 应该跟更深层交换 → 严格说 Dirichlet @ atmospheric 更对
- **但**：KG (NH4 吸附) 在 600 cm 以下 = 0，深层硝化贡献本来就极小 → 实际影响可以忽略
- **优先级最低**，论文可写成"边界条件假设"

**论文角度**：
- (a) 是必修——审稿人会问 "Did you include root respiration in O2 demand?"
- (b) 是 "interesting limitation"，可以放 discussion，论文不必修
- (c) 一句 model description 解释清楚就过

**触发条件**：等 §B-N（N 化学 MM 改造，2026-05-16 启动）跑通全年、NH4 量级落到 5-10 mg/kg 后，再回头看 f_aerobic 季节图。如果生长季 f_aerobic 还是 > 0.8（应该是 0.5-0.7），就先做 (a)。

---

## §6. 已知 Limitation：暴雨饱和期 27 mm 水量残差

**现象**（KT=372 雨高峰）：TopWaterInput=76 mm，StorageChange=49 mm，**残差 -27 mm**。

**物理诊断**：KT=372 顶 14 节点 hh > 0（最高 +26 cm），土壤顶 30 cm 进入饱和流。VG lSFCC=1 路径下 `hh ≥ -1` 时 Theta_LL = θ_s 钉死，水进 RHS 但 θ 不变 → storage 不增。DTheta_LLh/UUh 在饱和区 = SV.Ss ≈ 1e-5（v25 hh 钳位修复后），时间项接收量约 0.0006 mm，远不够吸收 27 mm。

**与 ClimaLand 公开版对照**：

| 维度 | ClimaLand | COMPASS lSFCC=1 |
|---|---|---|
| 状态变量 | `ϑ_l, θ_i` 独立 | `hh` 独立，θ 反算 |
| 饱和判据 | `ϑ_l > ν_eff` | `hh ≥ -1` |
| Theta 钉死 | 不存在（ϑ_l 可超 ν） | 钉死 = θ_s |
| 饱和区 ψ | `(ϑ_l - ν_eff)/S_s` | hh 求解器自由解 |
| dψ/dϑ 饱和区 | 1/S_s | SV.Ss |

**核心冲突**：公开版基于"hh 主导"模型（θ 是 hh 的函数，θ_s 是上限），ClimaLand 基于"ϑ_l 主导"模型（ϑ_l 可超 ν，ψ 是 ϑ_l 的函数）。**两者求解架构根本不同**——不是参数调整能切换。

**多相系统增加的复杂度**：ClimaLand 没有显式气相方程，可以简单允许 ϑ_l > ν。COMPASS 有 `lSoilairefc=1` 干空气方程，孔隙 `ν = θ_液 + θ_冰 + θ_气`，若允许 Theta_LL > POR → Theta_V < 0 物理破产。必须配套：饱和时 Theta_V = 0 + Air_sub 跳过求解 + 气-水交界面切换逻辑。

**完整变饱和流 + 三相运动改造计划**（1-2 周量级，Route B）：
1. Theta_LL/UU 允许 > POR（calculateTheta_LL.m, calculateTheta_UU.m）
2. DTheta 饱和区 = SV.Ss（已修）
3. Theta_V 饱和钉零 + 气相方程切换
4. storage 公式包含积水：`storage_mm = Σ(θ × DeltZ) + max(hh(NN),0) × 10`
5. 蒸发模型饱和修正：hh(NN)>0 时 rss→0，蒸发优先消耗"积水"（lake-evap 模式）

**当前决策**：接受 27 mm 残差作为已知 limitation，先跑全年验证整体物理。Route B 列下一阶段研究 work。

---

## §7. 关键参考

### 7.1 关键文件

```
src/COMPASS.m                                                  # 主循环
src/+core/MeteoData.m                                          # 风速 U2→U10 转换 convert_wind_height
src/ebal_clm45.m                                               # 能量平衡 (CLM4.5)，V8 Fourier-G 注入第 624 行附近
src/ebal_clm45_profile.m                                       # 能量平衡 (剖面版)
src/heatfluxes_clm45.m                                         # 干湿叶分离 lE_dry/lE_wet
src/+soil/+water/solveSoilMoistureBalance.m                    # 水方程入口
src/+soil/+water/solveTridiagonalMatrixEquations.m             # Tridiag 求解
src/+soil/+water/calculateMatricCoefficients.m                 # 矩阵系数
src/+soil/+water/calculateBoundaryConditions.m                 # BC 组装
src/+soil/+water/+conductivity/calculateDTheta_LLh.m           # DTheta 导数
src/+soil/+water/+conductivity/calculateDTheta_UUh.m
data/output_full_year/                                         # 当前跑输出
PlotResults/paper/                                             # 论文图脚本
```

### 7.2 参考代码（外部仓库）

```
F:/STEMMUS_SCOPE-main/                                          # GitHub 公开版（含冻融）
F:/STEMMUS_Zhujizixing -马铃薯 - 变步长/                       # 朴素 STEMMUS（无冻融）
D:/OneDrive/CTSM-release-clm5.0.37/                             # CLM5 源码（B-1 / B-4 参考）
D:/OneDrive/应用/Model/ClimaLand/                              # ClimaLand 公开版（§6 参考）
```

### 7.3 关键诊断脚本

| 脚本 | 用途 |
|---|---|
| `tools/download_era5_flux_soiltemp_gee.py` | ERA5 GEE 下载 (Colab) |
| `PlotResults/paper/plot_era5_vs_compass.py` | ERA5 vs COMPASS 全年对比图 |
| `PlotResults/paper/plot_psil_paper.py` | 叶水势 vs 实测主图 |
| `tools/calibrate_psil_to_obs.py` | 植物水力 case scan |
| `tools/calibrate_wind_z0.py` | 风速 z0 率定 |
| `tools/make_IC_from_yearfull.py` | 从全年 KT 取剖面生成新 IC |
| `tools/verify_freeze_coupling.py` | 冻融耦合诊断 |
| `tools/diagnose_h_low_r.py` | H 相关性低诊断 |
| `tools/diagnose_freeze_thaw_latent.py` | 冻融潜热量级估算 |
| `tools/compare_G_v4_vs_orig.py` | G_ebal vs QET vs ERA5 对比 |
| `tools/diagnose_rss_winter.py` | rss 冬季合理性诊断 |

### 7.4 节点编号约定（容易踩坑）

```
CSV 输出：N1=顶层(0cm), N129=底层(1000cm)
MATLAB 内部：m=1=底, m=128=顶
QLiquid/RWU/CondHydraulicAvg 文件用 MATLAB 内部 m
SoilMoisture/SoilMatricHead CSV 用 m_csv = 130 - m_matlab
```

### 7.5 验证清单（每次跑完检查）

```bash
# 1. 暴雨期 hh(NN) 是否反弹（应 > -1000 cm，不是 -15000）
awk -F',' 'NR>=370 && NR<=378 {print $1, $2}' data/output_full_year/Sim_output_SoilMatricHead.csv

# 2. TopWaterInput 是否接近 Throughfall（应 80%+）
paste <(awk -F',' 'NR>=370 && NR<=378' data/output_full_year/Sim_output_Throughfall_mm.csv) \
      <(awk -F',' 'NR>=370 && NR<=378' data/output_full_year/Sim_output_TopWaterInput_mm.csv)

# 3. Runoff 是否归零或 < 5 mm 总和
awk -F',' 'NR>=370 && NR<=378' data/output_full_year/Sim_output_SurfaceRunoff_mm.csv

# 4. 10/30/50/70/100 cm θ 是否上升（10cm 必须 > 0.40）
awk -F',' 'NR>=370 && NR<=378 {print $1, $16, $22, $26, $32, $39}' data/output_full_year/Sim_output_SoilMoisture.csv

# 5. 重画对比图
python tools/compare_rain_after_4bug_fix.py
```

---

## §8. 历史改动日志（倒序）

### v33 (2026-05-13) — Thomas TDMA 数值病态修复 ⭐⭐⭐

**问题**: 大雨步 KT=36 (10/3 11:00, 39 mm/h) Picard 第一次迭代解出 hh 残差 3926 cm (39 米水头)，数值病态。Substep 缩 dt 到 28-56s 跑 39+ 子步，但每个 sub-step 都"擦边收敛"(Water ≈ 95-99 cm，Tol_Water=100 cm)，让 Theta_LL 慢慢偏高，全年累积单步 dS spike +23 mm，年绝对残差 749 mm = 8.5% 降雨。

**根因**: Thomas TDMA solver 在大雨边界条件下数值病态，自 v25 起就存在但被年总残差稀释看不出 step-level spike。**非 four_quadrant / Anderson / Tol_Water 100 直接造成** (A/B 实验排除)。

**修复 (2 处)**:
1. **对数空间自适应 SOR** (`src/+soil/+water/solveTridiagonalMatrixEquations.m:50-68`)
   - 当 `max(CHKWATER) > 200 cm` (数值病态信号), ω=0.5 在 log|hh| 空间混合 (几何平均)；其余 ω=1
   - hh 跨 1e-6 to 1e7 cm 量级，log 空间更符合 VG 曲线分布
2. **Tol_Water 100 → 50 cm + 环境变量钩子** (`src/COMPASS.m:1339-1346`)
   - 默认严收敛，冻融期出 FATAL 时 `STEMMUS_TOL_WATER=100` 兜底
   - Tol 扫描 (100/50/20/10/5/1) 证明 50 是甜蜜点 — 浅层不会"假收敛"也不会跑满 maxiter

**验证**:

| 指标 | 修复前 (v32_full) | 修复后 (v33) | 改善 |
|---|---|---|---|
| 72h 雨季 KT=36 残差 | −9.35 mm (8.5%) | **−0.31 mm (0.28%)** | 30× |
| 全年残差 | −749 mm (8.5%) | **−39.4 mm (3.9%) ✓ OK** | 19× |
| KT=36 首次 Picard Water | 3926 cm | **528 cm** | 87% ↓ |
| KT=36 子步数 | 39+ | 20 | 2× ↓ |
| 冻融期 FATAL | 0 | **0** | - |
| 全年耗时 | ~1 h | 4.3 h | 4× ↑ |

**物理输出几乎不变** — VWC RMSE 各深度 ±0.001、叶水势 −2.658→−2.655 MPa、ERA5 LE/H/Rn/G 完全一致。**纯数值层修复，不改变科学结论**。

**vs STEMMUS-SCOPE 公开版差异**:
- 公开版用朴素 Tol_Water=1 cm + omega=1 (无 SOR)，冻融期容易 FATAL 退出
- COMPASS v33 用 50 + log-SOR + Anderson + adaptive substep 四件套，雨季和冻融期都稳

**详见**:
- `docs/figures/v33_vs_v32_comparison.md` 七件套对比报告
- 记忆 `project_compass_thomas_sor_fix.md` (skill 调用)
- 记忆 `project_compass_kt36_rain_bug.md` (诊断历程，A/B 排除)

---

### v30.1 (2026-05-11) — 全年跑通 + 出图流水线 + ERA5 对比

**全年模拟**（`output_v30_yearfull_v8`，lFreezeLatentCoupling=0）
- ✅ 全年 8760 KT 零 FATAL 跑通，StepFailed 全 0
- ✅ OUTER-STEP-SKIP bug 修复（子步耗尽后 `continue` 跳过输出装配，消除死锁）
- ✅ persistent fid 污染清除（删 `calculateHydraulicConductivity` / `calculateMatricCoefficients` 诊断块）
- ✅ checkpoint 目录 mkdir 兜底（防止首次运行因目录不存在保存失败）

**一键出图流水线**（`PlotResults/run_six_panel_analysis.m`）
- ✅ 从 `config_stemmus_scope.txt` 自动读 OutputPath 推断 tag，零手动配置
- ✅ 六项标准图（01_soil ~ 06_conservation）移至 `PlotResults/`，与 `tests/` 彻底分离
- ✅ `overwrite=false` 默认保护已有图（二次调用不覆盖）
- ✅ 无 log 文件时 conservation 图 Picard 残差子图显示 `(no log)` 不崩溃

**ERA5-Land 对比**（`PlotResults/plot_era5_vs_compass.py`）
- ✅ 自动推断 tag，图存 `docs/figures/<tag>/07_era5/` 不覆盖其他 run
- 对比结果（v8 vs v30_single_yearfull，日均）：

| 变量 | v8 RMSE | single RMSE | v8 R | single R | 结论 |
|---|---|---|---|---|---|
| LE | 18.3 | 17.9 | **0.915** | 0.912 | 相当 |
| H  | 15.7 | 15.3 | 0.611 | 0.595 | 相当；H 与 ERA5 本身差异大 |
| Rn | 10.2 | 10.2 | **0.983** | 0.982 | 无差 |
| G  | 12.3 | 12.3 | 0.564 | 0.576 | 无差 |

> **H 偏差分析**：v8 与 ERA5 日均 H 的 RMSE=15.7 W/m²、R=0.611，与 single 几乎相同。ERA5-Land 空间分辨率 0.1°（~10 km），而 COMPASS 对应单点苹果塬；H 对局地地表粗糙度/LAI/土壤质地极敏感，大尺度再分析与站点模型之间存在结构性偏差，不属于 COMPASS 本身的误差。LE 的 R=0.915 已显示能量分配整体合理。下一步优先改善**根系分配**（RWU 深度分布），而非追 ERA5-H。

**土壤水分 vs 实测（v8 vs single，全年 KGE/RMSE）**

| 深度 | v8 RMSE | single RMSE | v8 KGE | single KGE |
|---|---|---|---|---|
| 10 cm | **0.0356** | 0.0362 | **0.746** | 0.734 |
| 30 cm | **0.0393** | 0.0409 | **0.803** | 0.796 |
| 50 cm | **0.0455** | 0.0466 | 0.708 | 0.724 |
| 70 cm | **0.0437** | 0.0442 | 0.688 | 0.721 |
| 100 cm | **0.0787** | 0.0798 | 0.205 | 0.251 |

> 浅层（10-30 cm）RMSE 小幅改善，深层（70-100 cm）KGE 退步——与风速转换和 lFreezeLatentCoupling 关闭的影响一致，整体属于同量级波动，非系统性提升。

### v30 (2026-05-10) — 论文准备阶段
- ✅ ERA5-Land vs COMPASS 对比基础设施（GEE 下载 + 时区/单位/G 残差对比脚本）
- ✅ 风速 U2→U10 对数廓线转换（`MeteoData.convert_wind_height`），长武 0.91→1.26 m/s
- ✅ 冻融潜热耦合 V1-V8 实验，**V8 (Fourier-G + ×100 单位修复) 最佳**：bias -16.7 / R 0.84 / 完美稳
  - V1-V7 失败原因分别为：崩溃 / 不稳 / H 相位反转 / 永不触发 / 死循环 / SCOPE-GUARD / 仍欠精度。详见 git log。
- ✅ 植物水力 case scan 4 case，**Case D (P50=-3.46, aa=3.5, Braeburn)** 最佳：RMSE 0.238 / Bias +0.055 / 100% in 2σ
- ✅ 论文图整理 PlotResults/paper/，叶水势主图 RMSE 0.256 / Bias +0.029 / 75% in 2σ
- ✅ data/output 清理释放 109 MB

### v29 (2026-05-09) — CLM5 fwet 截留蒸发物理（13 处文件，~150 行）
- ✅ 统一截留参数命名（`canopy_storage_per_LAI=0.1`，CLM5/NoahMP 全球标准）
- ✅ 冠层水库状态变量 `h2ocan`，主循环重排：Gash → 收水 → drip → fwet → ebal → 扣 lE_wet → 守恒
- ✅ fwet 拆分干湿叶面（`heatfluxes_clm45.m`），CLM5 Tech Note Eq 5.97
- ✅ ebal_clm45 + ebal_clm45_profile 同步加 interception 入参，**SPAC 输入 Tp 必须用 lE_dry 不用 lEctot**
- ✅ Sim_output 新增 4 个诊断字段（h2ocan_mm / fwet / lE_CanopyIntercept_W / canopy_water_residual_mm）
- ✅ 断点续跑加 5 个 v29 状态字段
- 24h + 120h 暴雨双验证：0 nonconv，0 守恒报警；暴雨期叶温被压低 3.1°C，截留蒸发潜热峰值 96 W/m²

### v28 (2026-05-09) — ustar 统一 CLM5 标准（4 处）
- ✅ MeteoData 入口 `max(1.0, U)`（实测+日→小时两个分支），对齐 CLM5 BareGroundFluxes:256 / NoahMP Noahlsm:2077 / ClimaLand gustiness=1.0
- ✅ resistances_clm45.m 删二次钳位 `Va = u`（u 已在入口钳）
- ✅ resistances.m 标记 dead code（W&V 2 区间稳定度已不被主流程调用）
- ✅ ustar 公式确认对齐 CLM5 FrictionVelocityMod.F90:473-487（4 区间稳定度）
- 长武 65.3% 时段 u<1 m/s，统一钳位让两套 ustar 路径用同一 u_eff

### v27 (2026-05-09) — CLM5 风格 RWU 重写（3 处）
- ✅ 根分布 empirical → Jackson 1996（`root_dist_mode=1→0`），24h 实测 0-100 cm 吸水占比 4%→56%
- ✅ CLM5 风格 hydraulic redistribution 简化版（`ebal_clm45.m:1003-1062`）：拆 RWU_pos/RWU_neg，`flux × discount` 仅净流入进 sink，避免负 sink 让 hh 飙到 +1e11
- ✅ 诊断输出 RWU_pos/RWU_neg/HydraulicLift_intensity
- 24h 水量残差 v25→v27：8.78→0.52 mm（17 倍改善）；hydraulic lift 强度 70.8 mm/24h

### v25 (2026-05-08) — hh 钳位 + LAI 物候对齐（多处）
- ✅ VG dry 端 `hd: -1e7→-3e4`（避开雅可比病态区，4 个量级 dry tail）
- ✅ lSFCC=1 主路径 4 文件 dry 钉死分支补全（calculateTheta_LL/UU + calculateDTheta_LLh/UUh）
- ✅ hh 钳位 `if hh(i) > -1e-6, hh(i) = -1e-6`（对齐公开版 fixHeat），全年 0 FAIL Iter=100
- ✅ LAI 物候 14 维 GA+fmincon 拟合（RMSE=0.348）：sowing_day 104→79, GMHU 300→72, PHU 4555→3921
- ⚠️ 三方程 SUPG 删除尝试**已回退**（KT=1 矩阵奇异，详 §4 / U5）
- ✅ 数值安全网保留：`DPgDZ` ±1e5 钳位、`KTT ≥ Lambda_eff`

### v16 (2026-05-08 之前) — 双孔修复
- ✅ Bug B: SinkIm 时序+守恒（dualPorExchange 移到 substep 开始前）
- ✅ Bug A+C: HYDRUS InitW 拆分（`VG.Theta_s -= ts_im`，θ_obs 拆 mobile + immobile）
- 72h 累积水量残差 v8→v16：12.57→3.37 mm（73% 改善）
- ⚠️ 100/200 cm 仍未响应（v16 反而比 v15g 更钉死，详 §5 / B-2）

---

## §9. 已解决问题档案（按问题编号）

<details>
<summary><b>U1. 饱和区水量不闭合（旧 466 mm/年残差）</b> — v25 hh 钳位修复</summary>

- **修复**：v25 起 hh 强制 ≤ -1e-6，承压物理改用 DSTOR 表面积水库表达
- **效果**：v27 24h 水量残差 0.52 mm（v25 同口径 8.78 mm/24h），全年残差预计 < 30 mm
- **仍未解决的小尾巴**：暴雨期 DSTOR 蒸发路径未接入 SCOPE rss（设计见 `docs/saturated_flow_route_B_design.md`，搁置）
</details>

<details>
<summary><b>U9. 全年水量不闭合（旧 780 mm 残差）</b> — v25 hh 钳位 + v27 RWU 修复</summary>

- **v25 全年实测**：全系统残差 -63.28 mm（占 P 6.4%），改善 92%
- **v27 24h**：0.52 mm（v25 同口径 8.78 mm/24h）
- **顶层 P=CI+Runoff+Infil 完美闭合**（残差 0.000 mm）
- **仍未解决（< 30 mm/年）**：vapor flux 与 Es 统计口径不完全一致；冰相年净增 166 mm 的边界精度
</details>

<details>
<summary><b>B0. 顶层 N1 撞 hcrit 矩阵奇异</b> — v25 hh 钳位修复</summary>

- **关闭依据**：v25 全年跑 0 次 FAIL Iter=100（旧 1014 次），0 个 nonconv_diag_KT 文件（旧 200+ 个）
- **真因**：不是 hcrit 撞底，是 `solveTridiagonalMatrixEquations.m` 没钳 hh > 0 → 冻土条件 `KfL_h × 10^(-Imped·Ratio_ice) → 0` 让 Thomas 主元 → 0，hh 飙到 +5e7
- **修复**：`if hh(i) > -1e-6, hh(i) = -1e-6`（对齐公开版 fixHeat）
</details>

<details>
<summary><b>S2. ustar 低风速失真（v25 实测 3.79× 偏差）</b> — v28 ustar 统一 CLM5 标准</summary>

- **v25 现象**：长武 65.3% 时段 u<1 m/s，CLM-style ustar / W&V ustar = 3.79×（p95 11×）
- **v28 修复**：MeteoData 入口统一 `max(1.0, U)`，下游两套 ustar 用同一 u_eff
- **物理依据**：CLM5/NoahMP/ClimaLand 三家陆面模型一致——阵风 (gustiness) 修正，零风时大气仍有 ~1 m/s 量级湍流脉动
- **结论**："3.79×"不是 bug，是物理；v25 失真是 W&V 路径未做钳位导致两套不一致；v28 统一后自然消除
</details>

<details>
<summary><b>U4. 18 处硬编码 -1e7 残留</b> — v25 主路径已修，余项不阻塞</summary>

- **当前状态**：剩余 -1e7 都在 lSFCC=0 / Clapp / 反算工具等主跑不走的路径
- **主跑路径** lSFCC=1 + lFreezing=1 已通过 `getSoilConstants.hd=-3e4` 生效
- **优先级降为 P1**：保持代码一致性，但不阻塞
</details>

<details>
<summary><b>U8. ET 时间分配异常（4-6 月 E/T=2.1 反常）</b> — v25 LAI 物候对齐修复</summary>

- **旧跑**：4-6 月 Es 1.70 mm/d / Tr 0.81 mm/d / E/T=2.10 ❌
- **根因**：早生长季 LAI 卡在 0.20（已通过 LAI 14 维拟合修复 → 4 月 LAI 0.61）
- **预期**：4-6 月 Es 降到 ~1.0，Tr 升到 ~1.5，E/T 接近 1.0
</details>

<details>
<summary><b>B-4 (lE_CI 部分). 截留蒸发潜热漏算</b> — v29 修复</summary>

- **旧 bug**：SCOPE lEtot 只算干叶气孔/边界层，没算湿叶面 fwet 蒸发；Gash 模型扣水量 ~260 mm/y 但对应 ~20 W/m² 年均潜热没出现在 lE 里
- **v29 实测**：暴雨期 lE_wet 峰值 96 W/m²，叶温被压低 3.1°C，守恒残差 < 1e-6 mm
- **仍未解决（→ B-4 剩余）**：S_freeze / S_storage / CANHS / PAH 详 §5 / B-4
</details>

<details>
<summary><b>2026-05-07 5 bug 修复成果</b> — 暴雨入渗物理建立</summary>

- 暴雨期 Picard 收敛（KT=372 仅 3 次迭代）
- TopWaterInput ≈ Throughfall（83.4/85.4 mm）
- 表层 10 cm θ Δ48h = +0.260（实测 +0.276）
- Runoff < 2 mm
- **未达标的部分**：30-100 cm 响应偏小 → 演变为后续 U12 / B-2 macropore 议题
</details>
