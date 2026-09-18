# 呼吸 / NSC 改进完整路线 (Respiration & NSC Improvement Roadmap)

> 2026-06-12 起。整合"自养呼吸 2× 高估 + xsmr 负债 −90"的全程诊断,落成分阶段路线。
> 关联: data/references/nsc_respiration/README.md (5 篇文献 + 公式),
> docs/superpowers/specs/2026-06-11-xsmrpool-carbon-buffer-design.md (xsmr 原 spec),
> 记忆 project_compass_clm5_xsmrpool_gap / project_compass_nsc_reserve_n。

---

## 0. 诊断总结 (已确证)

| 问题 | 数据 | 根因 (已核源码/文献) |
|---|---|---|
| 自养呼吸高估 2× | Ra=738 vs 实测~366;**真 CUE=(GPP−Ra)/GPP=0.34** vs 正常 0.4-0.6 | 与 CLM5-FruitTree 母本同病 (Dombrowski 2022);R_m 基准偏高 |
| xsmr 负债 −90 (苹果"死") | 全年末 C_xsmr −90,夏季 Δ−73 | **非池坏/非抄错**:CLM5 港代码逐行核对一致。真因=CLM5 `dayscrecover=30` 的回血是**季节尺度**,对付**日夜循环**太慢(夜攒 1.86 > 日还限速 0.77),债爬到 −90 |
| C_labile 堆 32 g 闲置 | none=0.05,tcmm=32(节流盈余碳) | COMPASS **缺** CLM5 的余碳下调 (excess_cflux→downreg,NutrientCompetition...:314-320),盈余碳无处去只能堆 |

**注意 (诚实):** `net_primary_productivity` 诊断口径=cpool 资助的显示生长,**不含 xsmr 负债**,会把 CUE 报成 0.62 (虚高);真 CUE 必须用 (GPP−r_a_actual)/GPP=0.34。

---

## 1. 核心决定:C_labile 的身份 (文献定调)

**C_labile = 快周转可溶性糖底物**,给**一切**供燃料(生长 + 全部呼吸,含夜),呼吸受其**底物限制**;
**C_storage(根/茎)= 淀粉慢储**,季节缓冲。这是标准两池 NSC,与 Jones 2024 快/慢池一一对应。

依据 (四篇一致,均已读 PDF):
- **Thornley & Cannell 2000** Ann Bot 85:55 — C 底物池供 growth+phloem+N+维持全部;Eq 17 在池上做维持 turn-down。**且第 59 页点名批负债那套**:"maintenance...unrelated to assimilate supply...carbohydrate concentrations can be driven negative, leading to nonsense results" = 正是 xsmr −90。
- **Jones 2024** Commun Earth Environ 5:148 — 快池 S_f 直接喂呼吸 R=(1−Yg)·kf·S_f。
- **O'Leary 2017** Plant Physiol 174:2261 — 夜呼吸 ∝ 碳水储 (starch r²=0.73)。
- **Bruhn 2022** Nat Commun 13:5650 — 夜呼吸下降=底物耗尽。

**推论:维持呼吸应从 C_labile(糖)付、被 β turn-down 限制;退役会变负的 xsmr 债池。** = nodebt 路线。

---

## 2. 分阶段路线

### P0 — 验证 (进行中)
- 并行跑 9 方案: none/bruhn/tc_mm(km2,km6)/jones + **nodebt(km2,km6)** + **fastrec(km2,km6)**。
- 切换靠 env `COMPASS_RESP_SCHEME` + `COMPASS_RESP_KM`(临时实验开关,P4 收掉)。
- **验收(预登记):** ① 真 CUE 进 0.4-0.6;② 负债压到 ~0 / 有界(苹果活);③ GPP/LAI/产量(C_fruit)不掉;④ C/N 守恒 <1e-8;⑤ C_labile≥0 且日内潮汐;⑥ 能复现 Bruhn 夜降 ~25%/8h。
- **产出:** 9 方案对比表 + 选定 (预期 nodebt) + K_m 定值。

### P1 — 定身份 + 退债池 (核心改动)
- 采纳 §1 定义:**维持呼吸从 C_labile 付,β turn-down 限制,无 xsmr 债池**(nodebt scheme=4 转正)。
- **维持优先级保留**:先付 R_m(GPP+C_labile),`NPP=availc−R_g` 给生长——CLM5 的"先保维持"逻辑不丢。
- **两池缓冲细化:** nodebt 可付上限从 `GPP+C_labile` 放宽到 `GPP+C_labile+f·C_storage`,让淀粉慢储在持续缺糖时 remob 兜底(对齐文献完整两池;FATES storage-cushion 同理)。
- xsmr 债池代码:**保留为可选分支**(scheme 开关),默认关;或整体退役——看 P0 实证 + 你拍板。

### P2 — 余碳下调 (可选,补 CLM5 缺口)
- COMPASS 缺 CLM5 的 excess_cflux→downreg。节流盈余碳现在堆 C_labile;nodebt 已用它付呼吸(消化一部分)。
- 若 P1 后 C_labile 仍异常堆积 / GPP 虚高:补一个余碳下调(碳用不完→调小光合),对标 CLM5 NutrientCompetitionCLM45defaultMod.F90:314-320。
- **优先级低**:先看 nodebt 是否已把堆积消掉。

### P3 — 标定 + 机制验证
- **K_m 定值:** 用 P0 的 km2/km6 包络,标到真 CUE 落 0.4-0.6 中段 + 生长季 β≈0.9(不误伤健康期)。标 `[site-calibrated]` (§5.5)。
- **Bruhn 夜降验证:** 检查 nodebt 的夜间 C_labile 耗→R_m 降是否复现 ~25%/8h(把 Bruhn 从候选转成验证靶)。
- **leaf Rd 口径:** 复核 r_a_actual 是否含 turn-down 后叶 Rd;补 resp_beta_sub 进 NetCDF 输出(当前在 nsc_diag 未注册)便于复盘。

### P4 — 清理 (定案后)
- 剪掉落选 scheme 分支(bruhn/jones/with-debt 三选一留),把胜者**提升进 config**(`cfg.resp_scheme` + `cfg.resp_Km`),撤掉临时 env 开关 + getPFTParameterSet 的 env 读取块。
- 删临时诊断 / 标 `[EXP]` 的代码;按 §7.3 走诊断生命周期。
- checkcode 清零。

### P5 — 文档 + 提交
- 更新 README(本路线 + 最终决定 + 文献)、output_catalog(新增字段)。
- commit + tag(如 v34.1-substrate-limited-respiration)。
- 更新记忆:C_labile 身份定论 + xsmr 债池去留。

### P6 — 剩余 limitation (诚实 + future)
- **CUE 到 0.71 的 gap:** nodebt 预期进 0.4-0.6(正常区);**实测 0.71 在正常区顶/之上,属果园特异**(管理/测量/不连续冠层),非标准呼吸物理能给。写成 limitation。
- **Route A (future GMD 卖点):** 不连续果园冠层 (Dombrowski 2022 点名空白) → 内冠层叶氮陡降 → 治叶 Rd 大头。重,留后续。
- **GPP 虚高:** 缺余碳下调可能让 GPP 略乐观 (1208 vs obs 1223 其实接近,但机制上要留意)。

---

## 3. 决策树 (一图流)

```
C_labile 是什么?
  ├─ "可溶性糖,烧给一切(含夜呼吸)"  ← 文献(T&C/Jones/O'Leary/Bruhn)一致 + COMPASS 原生收支
  │     └→ nodebt:R_m 从 C_labile 付,β 限制,退债池,storage 做慢缓冲   ★ 推荐
  └─ "纯生长池,维持走 GPP+缓冲"  ← CLM5 cpool 定义
        ├─ fastrec:留 xsmr 池,但回血不限速(白天剩余 GPP 立刻还清)
        │     ⚠ 风险:还债优先于生长,债大时饿生长(你 Q2)
        └─ (原 xsmr 30天回血)  ← 已证 −90,T&C 2000 点名"nonsense"
```

## 4. 验收标准 (预登记,所有方案同口径)
1. 真 CUE = (GPP − Σr_a_actual)/GPP ∈ [0.40, 0.60]
2. xsmr/负债 end ≈ 0 或有界(不持续变负;苹果活)
3. GPP 与 none 偏差 < ~10%;LAI 峰 ~2.8;C_fruit 不显著掉
4. C/N 守恒残差 < 1e-8
5. C_labile ≥ 0 全程,且日内有潮汐(被用起来)
6. (加分)复现 Bruhn 夜降 ~25%/8h

## 5. 一句话路线
**把 C_labile 定成"糖,烧给一切"(文献一致)→ 维持呼吸从它付、β 底物限制、退役会变负的 xsmr 债池(T&C 2000 点名该退)→ 淀粉慢储做季节缓冲 → 标定 K_m → 清理提交 → 剩余到 0.71 的 gap 归为果园特异/Route A future。**
