# COMPASS C/N 循环统一审计 + 融合设计

> **目标**:终结"LPJ + CLM5 + CLM5-FruitTree + Aerts + EPIC 多套机制拼接"的缝合怪,
> 梳理一套**属于 COMPASS 自己的、自洽的、有出处的** C/N 循环。
>
> **本文档已于 2026-06-06 统一重写**:只保留经"验证 workflow(7 bug 对抗反驳)+ 全 24 PFT 矩阵
> + 全生育期归类 + 参数血统表"验证后的**权威结论**。早期未经验证的初判(曾误判 aleafstor/cn_mode/
> 三重死亡链)已删除,以本版为准。
>
> 来源:2026-06-06 多 agent 审计 + 验证 workflow(wuxghbil4)+ 用户决策
> ("不要拼接 / 全 PFT 查 / 主用一个删另一个 / 先验证再动手")。
> ⚠️ 本文档只审计+设计。实施前每步须经用户批准 + compass-coding-standard skill。

---

## 0. 缝合怪从何而来

COMPASS C/N 循环缝了 **7 个来源**,每加一个补丁就多一处张力:

| 来源 | 贡献机制 | 时间基准 |
|---|---|---|
| LPJ-GUESS (Sitch 2003) | NPP×f_alloc / Sykes 物候 / APHEN_MAX / L&T 呼吸 / nstore 储备目标 | gC/m²/d, base 5/10°C |
| CLM5 (FlexibleCN) | 浮动 C:N / fcur 越冬 / onset-offset ramp / 日长落叶 / xsmrpool | gC/m²/s, Q10 base 20°C |
| CLM5-FruitTree | 4相 GDD 分配 / Cesaraccio 冷激 / peaklai / aleafstor / storage_flag / onset-xfer | GDD 驱动 |
| Aerts 1996 | retransl_fraction=0.7 | — |
| EPIC (Williams) | frPHU S-curve / spring_release / harvest_DS | 积温 |
| SITE (长武反标) | f_alloc / GDD_bloom / TSum_fruit | — |
| COMPASS 自创 | applyNSCCoupling NSC 抽储 | 无出处 |

---

## 1. ⭐ 7 个 bug 最终判定(验证 workflow + 对抗反驳,2026-06-06)

| Bug | 最终判定 | 严重度 | 致 year-2 死? | 要修? |
|---|---|---|---|---|
| **B1 APHEN_MAX 提前杀叶** | **REAL,已修且忠实 CLM5** | HIGH | 是(主因) | ✅ 已修 |
| **B2 aleafstor 未实装** | **NOT_A_BUG** | — | 否 | ❌ 已实装,误判 |
| **B3 NSC 抽 storage** | **存疑(被反驳)** | — | 否(反驳) | ⚠️ 实测验证 |
| **B4 cn_mode 双判定** | **NOT_A_BUG** | — | 否 | ❌ 忠实 CLM5,误判 |
| **B5 peaklai C/N 不对称** | **REAL(唯一存活真 bug)** | MEDIUM | 否 | ✅ 修 |
| **B6 f_alloc 三层覆盖** | NEEDS_RUNTIME_TEST | MEDIUM | 否 | 待测 |
| **B7 gdd5 双累积** | NOT_A_BUG(注释撒谎,逻辑对) | LOW | 否 | ❌ 只改注释 |

### B1 — APHEN_MAX 提前杀叶 ✅ 已修(致 year-2 死的主因)
- **机制**:苹果 budburst DOY79 + APHEN_MAX=210(LPJ growth.cpp:50)= DOY289 强制 phen=0 → LAI 归零。
  早于 CLM5 日长落叶(crit_dayl=10.9h,长武 lat35.25 在 DOY294-295 才触发)。提前 ~6 天杀叶
  → 丢掉采后(DOY285采收)挂叶光合充储备窗口(CLM5-FruitTree storage_flag→fcur=0)→ 储备饿死 → year-2 死。
- **验证铁证**:DOY289 dayl=11.06h(>10.9 日长未触发),DOY295=10.85h。
- **修复**(seasonDecidPhenology.m:323-327):`has_dayl_offset = crit_dayl_h>0`,有日长机制就跳过 APHEN_MAX。
- **为什么忠实 CLM5**:CLM5 季节落叶(CNPhenologyMod.F90:957)**根本没有 APHEN_MAX**——纯日长触发。
  APHEN_MAX 是 LPJ 因为缺日长落叶才需要的兜底。COMPASS 所有 season_decid PFT 都有 crit_dayl_h,
  所以 APHEN_MAX 对它们全是冗余兜底,跳过安全。**(全 PFT 矩阵曾担心"7 PFT 主用",已调和:它们同时有日长机制,APHEN_MAX 是冗余非主用。)**

### B5 — peaklai C/N 不对称 ✅ 唯一存活真 bug(中等,不致 year-2 死)
- **机制**:peaklai 砍叶时,C 侧(computeCarbonAllocation.m:175,221)把 capped 叶碳 100% 进 C_storage,
  但 N 侧(computeOrganDemand.m:113-125)只把 f_leaf 归零、**没对应 N 进 storage**(注释明说"drop without reroute")。
  → storage 进 C 无 N → storage C:N 偏高 → 春萌 onset-xfer 注入偏高 C:N。
- **忠实 CLM5 修法**(NutrientCompetitionFlexibleCNMod.F90:1653-1657):peaklai 时 aleaf→1e-5,
  freed 份额**重分配给真器官(根,phase1)/果/茎**,而非孤立进 storage。c_allometry 和 n_allometry
  用同一组 f 系数 → 接收器官 C 和 N 一起涨,**无 C:N 漂移**。
- **COMPASS 修法**:peaklai 的 freed f_leaf 在 C 侧和 N 侧都转给根(对称),而非 C 进 storage/N 留 labile。

### B3 — NSC 抽 storage ⚠️ 存疑(被反驳,需实测)
- **机制存在**:applyNSCCoupling.m:90-92 把 C_labile+C_storage 当呼吸燃料;:124-131 C_labile<0 时抽 C_storage(只floor=0)。
- **反驳**:"冬季烧穿 storage 48→5"被模型输出推翻——**C_labile 是大缓冲累积器**(森林 34→155,苹果 v54 到 62.7gC),
  正常 MR 先被 C_labile 吸收,轮不到 storage。"C_labile 是过路池"是假前提。
- **状态**:存疑。若 B1 修复后 year-2 仍死,再查 storage 真凶(可能是 E2 双释放,见 §3)。
- **CLM5 对照(若要修)**:xsmrpool 缓冲池(NutrientCompetitionCLM45defaultMod.F90:707-744),
  MR 缺碳记 xsmrpool(可负),春季光合优先恢复,storage 不被 MR 直接抽。

---

## 2. ⭐ 死亡链最终结论(收敛)

**year-2 死 ≈ B1(APHEN_MAX)一个主因,已修。**
- ~~三重混用叠加~~ → 删除(早期误判)。
- ~~aleafstor 未实装~~ → 删除(B2 已实装)。
- ~~NSC 烧穿 storage~~ → B3 被反驳,存疑。

**最快事实判定**:只带已修 B1 跑一次。
- year-2 活 → B1 是唯一主因,B3 反驳对,收工。
- year-2 还死 → 查 storage 真凶(B3/E2)。

---

## 3. ⭐ 参数血统表 + 7 大混用热点(防复发核心)

**完整表见 `data/_taxonomy.txt`(16KB)。** 每参数标:苹果值/来源/驱动/冲突风险。

### 7 大混用热点(按风险)
| # | 热点 | 两来源 | 谁活谁死 | 风险 |
|---|---|---|---|---|
| **E2** | **春季储备双释放(CRITICAL,两个都活!)** | EPIC springStorageRelease(spring_release_days=60) + CLM5-FT springOnsetXfer(fstor2tran=0.5/ndays_stor=50) | **都活,同一 N_storage→生长 两套 remobilize → 可能双释放** | 🔴 最高,storage 异常隐藏推手,之前漏查 |
| B1 | budburst 三来源 | LPJ Sykes k_chillb + CLM5 onset_gdd + FT Cesaraccio crequ/baset | use_cesaraccio_chill=1→Cesaraccio 活,Sykes k_chillb=1150 **死但在表** | 🟡 改 k_chillb 对苹果无效 |
| A1 | 分配两来源 | SITE f_alloc_* + FT 4相 fleafi/arooti | use_fruittree_alloc=1→FT 活,f_alloc_* 对C侧**死**但仍用于 storage-split | 🟡 改 f_alloc_leaf 调叶碳"静默无效" |
| E1 | retransl 重复 | Aerts retransl_fraction(活) + LPJ nrelocfrac(**死,从不读**) | 都=0.7 掩盖 | 🟡 死参数 |
| D1 | 茎 C:N 单/双池 | FT CN_target_stem=60(demand) + 双池 CN_lw=60/CN_dw=500(turnover) | demand 和实际存的背离 | 🟡 latent |
| C1 | 呼吸 LPJ+CLM5 配对 | Sitch r_m10=0.066 + CLM5 R_g_fraction=0.11(原 Sitch 0.25) | g(T)是 Sitch 式但配 CLM5 的 0.11 | 🟡 已审计迁移 |
| B2 | 落叶日长 vs 叶龄 | CLM5 crit_dayl + LPJ APHEN_MAX | guard 已解决(B1) | ✅ 已解 |

### 防复发铁则
1. **开关门控的死参数陷阱**:改 `f_alloc_*` / `k_chill*` / `Vcmax25` 前,先查管它的开关
   (use_fruittree_alloc / use_cesaraccio_chill / lLUNA)——开关关掉时这些参数是死的,改了"静默无效"。
2. **两个真共存的热点(E2 春释放 / C1 呼吸)需物理调和**,不是切开关能解决。
3. **E2 是新发现高危**:理想每跑只留一套春季释放机制(springStorageRelease 或 springOnsetXfer)。

---

## 4. 全 24 PFT 机制主从矩阵

**核心:多数"打架"其实是「同一机制对不同 PFT 主从相反」,不是真冲突。删机制前必查全 PFT。**

### 4.1 PFT 分类
| 类群 | PFT | 物候引擎 | 储备 |
|---|---|---|---|
| 常绿 | 1,2,4,5,9 | evergreenPhenology(phen=1 常数) | 无 |
| 季节落叶 | 3,7,8,11,12,**23**,24 | seasonDecidPhenology(Sykes+CLM5日长,含苹果) | 8个有 |
| 胁迫落叶 | 6,10,13,14 | stressDecidPhenology(4累积器,草无日长) | 灌有草无 |
| 作物 | 15-22 | cropPhenology(DAISY DS,一年生) | 无 |

### 4.2 禁删机制(多 PFT 主用)
APHEN_MAX(7 season_decid,但有日长时是冗余兜底)/ crit_dayl(7)/ Sykes(6,不含苹果)/
静态 f_alloc(23个)/ has_storage_pool(8)/ retransl_fraction(18)。

### 4.3 可门控机制(仅苹果,use_* 标志,其他 PFT 自动跳过)
Cesaraccio 冷激 / 4相 GDD 分配 / aleafstor 采后叶→储备。

### 4.4 苹果(PFT23)= 最复杂(bug 重灾区)
季节落叶基础 + 整条果实线:Cesaraccio + 4相GDD + 开花GDD_bloom + 果实DS时钟 + 采收 + aleafstor + peaklai。
叠加机制最多 → bug 最多。

---

## 5. 融合方法论(用户共识)

**策略 B「按过程选优」+ 四条铁律 + 全 PFT 主从锁。** 不推倒重来,把补丁理顺成正式分层。

### 四条铁律
1. **一过程一主机制**:每过程一个当家,同源机制降级 fallback(门控,非删)。
2. **统一时间尺度/单位**:全 gC/m²/s,诊断转 /d。base 温度文档化。
3. **C/N 严格对称**:任何碳流必有对应氮流,同套拆分(B5 就是违反此条)。
4. **参数血统标签**:每参数标 `[SOURCE: 模型 文件:行]`,同过程禁混不同血统参数(防 §3 热点复发)。

### 全 PFT 主从锁(安全底线)
删机制前查全 24 PFT:主用 PFT >1 → 禁删(只能门控);=1(仅苹果)→ 可门控降级。
**一处条件门控,按 PFT 自动主从**(如 APHEN_MAX 跳过逻辑 = `crit_dayl_h>0`)。

### GMD 反"缝合怪"说法
"COMPASS 非简单拼接,而是对每个 C/N 过程**机制选优**:物候 = LPJ Sykes 冷激(温带成熟)+ CLM5 日长落叶
(大叶标准)+ FruitTree Cesaraccio(果树专用),分配 = CLM5-FruitTree 4相,储备 = CLM5 fcur+onset-xfer,
并通过统一时间尺度、C/N 对称约束、参数血统标注、PFT 主从门控确保有机融合而非冲突。"

---

## 6. 实施优先级(待用户批准)

**P0 — 验证 B1 是否足够(零成本,直接出事实)**:
- 只带已修 B1 跑 2-3 年,看 year-2 活不活 + storage 跨年闭环。

**P0.5 — 若 P0 不够,修确认的真 bug**:
1. B5 peaklai C/N 对称(freed 叶份额 C+N 都转根,忠实 CLM5:1653-1657)。
2. E2 春季储备双释放(springStorageRelease / springOnsetXfer 二选一)。
3. (条件)B3 xsmrpool 缓冲(仅当实测确认 MR 真烧穿 storage)。

**P1 — 理顺血统(不致死)**:
- 死参数清理(nrelocfrac/k_chillb 苹果死/f_alloc_* 苹果死),加血统标注。
- B7 注释纠正,D1 茎 C:N 单双池统一,C1 呼吸配对调和。

---

## 7. ⭐ 实施完成记录(2026-06-06,改完后状态)

本节记录基于上方设计**实际实施**的改动 + 验证状态。以此为"改完后"真相。

### 7.1 双池储备(C/N_root_storage + C/N_stem_storage)
**动机**:CLM5 每器官有 storage 子池(frootc_storage/livestemc_storage...),COMPASS 旧版单一 C_storage
→ peaklai 砍叶碳进单池无对应 N(B5 C:N 漂)+ "转根 vs 转储备"二选一矛盾。
**做法**:单池拆成 **粗根储备 + 茎(边材)储备**(生理:NSC 主存于粗根+干边材,叶/细根可忽略;
Chapin 1990 / Loescher 1990)。
**关键设计(用户两次纠正)**:
- **按来源器官归库,非固定比例**:根分配 (1-fcur) → 根库;茎/叶/果 → 茎库;叶秋季回流 → 茎库(就近)。
  两池此后**各自独立演化,比例自由漂移**(不被锁定)。
- **f_storage_root(=0.55)仅用于初值分配**(KT=1 无历史),运行时不再用。site-calibrated 标注,待实测。
- **融合型 onset**:春萌把 root+stem 合为总储备喂 springOnsetXfer(NSC 韧皮部可流动),按抽前占比回扣两池。
- **overdraft 先抽茎库后根库**(保根库给春萌 bootstrap)。
**改动文件**(12):setInitialNitrogen / computeCarbonAllocation / allocatePlantN / turnoverCarbonMirror /
turnoverAndHarvest / runPlantCNStep / applyNSCCoupling / COMPASS.m(守恒求和+输出) /
registerPlantNitrogen(输出注册) / getPFTParameterSet(f_storage_root) / runPlantNDemand(CONS-1 修) /
springOnsetXfer(单池接口保留,收合并总储备)。

### 7.2 B5 修复(peaklai C/N 对称)— 用户提醒"C/N 必须一致"
**最终做法(忠实 CLM5)**:peaklai 触发时 **C 侧和 N 侧都 `f_root += f_leaf; f_leaf=0`**
(computeCarbonAllocation + computeOrganDemand 两侧),C 和 N 共用同一 f_root → freed 叶碳氮一起转根
生长 + (按来源)根储备,**C:N 严格对称,无孤立 C 流**。出处 CLM5-FruitTree
NutrientCompetitionFlexibleCNMod.F90:1653-1657(aleaf→1e-5,aroot += freed;c/n_allometry 共系数)。
删除了中途的 A_C_peaklai_to_storage 特殊通道(那个造成 C/N 不对称)。

### 7.3 E2 解决 + 春季释放统一
**问题**:springStorageRelease(旧,Tagliavini 2003,储备→labile)+ springOnsetXfer(CLM5,储备→xfer→器官)
两套同窗双抽储备。
**做法**:**删 springStorageRelease 调用**,onset-xfer 成为唯一春季释放;**去掉 onset-xfer 的 use_leafc_sla 闸**
→ 所有 8 储备 PFT(3,6,7,8,10,11,23,24)、两个 LAI 引擎(leafc_sla/EPIC)都走 onset-xfer。
全 PFT 核实安全(8 PFT 都有 budburst_doy;EPIC 下 onset-xfer 走碳流不驱 LAI,无冲突)。
清死变量 use_leafc_sla / spring_signal。springStorageRelease.m 文件保留(无调用,git 历史清晰)。

### 7.4 其余 7 热点处理状态
- **B1 APHEN_MAX**:已修(crit_dayl PFT 跳过叶龄 cap,忠实 CLM5 无 APHEN_MAX)。
- **A1/B1-budburst/E1/D1/C1/B7**:开关门控的死参数(k_chillb 苹果死/f_alloc 苹果死/nrelocfrac 死)
  + 单双池茎 C:N + 呼吸配对 + gdd5 注释 —— **均非致死,标 P1 后续清理**(加血统标注、删死参数)。未在本轮动。

### 7.5 验证状态
- **冒烟(1天,output_v57_dualstorage_1day)**:exit 0,24 KT 不崩;守恒 **C 3.4e-13 / N 3.6e-15**(机精);
  双池初值正确(根 112.7 + 茎 92.2 = 204.9,根占比 0.550)。
- **2 年跑(output_v58_dualstorage_2yr)**:进行中。验 year-1 守恒 + **year-2 复活到 LAI~2.8**
  (v54/v55 死在 year-2,本轮双池+onset-xfer 应解)+ 多年 storage 跨年闭环。

### 7.6 CLM5 对照审查结论(workflow w8gb658ps)
守恒全闭合(对抗验证确认双池所有碳氮流)。CONS-1(runPlantNDemand 漏改)= 唯一致命,已修。
DS-1(融合 onset)/DS-4(overdraft 抽储)非严格 CLM5 但守恒安全、对抗验证判"可接受简化"。

---

---

## 8. ⭐ labile/storage 碳路径诊断(2026-06-06,用户逐层追问逼出的根因)

**问题现象**(v58 2yr,DOY170):茎(边材)储备从 92.2 被抽到 1.6;C_labile 生长季钉在 0。

**用户连环质疑**:茎库为啥0 → labile 为啥空 → 有 labile 为啥还抽 storage → 有降呼吸吗 → 分配有度吗 → 是不是 N 主导没改彻底。

### 8.1 实测碳路径(v58 KT4320,逐日累积量)
| DOY | labile | cumGPP | cumRa | C_leaf | 根库 | 茎库 |
|---|---|---|---|---|---|---|
| 1 | 40.6 | 0 | 0.3 | 0 | 112.7 | 92.2 |
| 60 | 21.0 | 0 | 19.9 | 0 | 112.7 | 92.2 |
| 79(春萌) | 13.5 | 1.4 | 28.7 | 7.3 | 56.4 | 46.1 |
| 130 | 0.0 | 140.7 | 102.7 | 59 | 58.1 | 42.1 |
| 170 | 0.0 | 440.7 | 270 | 94 | 61.7 | **1.6** |

### 8.2 根因链(三者叠加,非单一)
1. **冬季维持呼吸吃 labile**:DOY1-79 GPP=0 但 cumRa 涨 28,labile 40→13(冬季呼吸第一杀手)。
2. **生长季 labile 是"过路池"**:C 主导分配 A_C=NPP×f_alloc(flexible_cn,**这部分是对的 C 主导**)每步把 GPP 进账几乎全分光 → labile≈0,**没留缓冲余额**。
3. **labile 一旦算微负 → Step5 overdraft 先抽茎库**:茎库被当一线缓冲,而其充值((1-fcur)茎分配)跟不上抽取 → 46→1.6。根库因"先保护根库 + 持续充值"而稳(56→61)。

### 8.3 澄清两个澄清点(纠正过程中的误判)
- **分配是 C 主导,改对了**:config=flexible_cn,走 A_C=NPP×f_alloc + clampPoolCN(超额C留labile)。
  审查中一度误看 c_follows_n 分支(A_C=A_N×CN,N主导)——那分支**未生效**(config 不是 c_follows_n)。
- **降呼吸机制存在但被架空**:applyNSCCoupling Step2 cap `NSC_avail = labile + 根库 + 茎库`,
  把 storage 算进"呼吸可用燃料",所以呼吸不降(cum_storage_drawdown_count=0)、改抽 storage。
  即"有降呼吸,但 storage 进了 cap → 降呼吸排在抽光 storage 之后"。

### 8.4 真正修法(待定,治本三选/组合)
1. **cap 不算 storage**(或只留越冬保护底):呼吸超 GPP+labile 即触发降呼吸(碳限制),storage 纯越冬。
   = 用户"有降呼吸 storage 不该0"的直接解。一行改 Step2 NSC_avail。
2. **分配留 labile 缓冲**:不把 GPP 全分光,留一部分在 labile(CLM5 cpool/availc 约束;NutrientCompetitionFlexibleCNMod.F90:338 availc + :478 excess_cflux + :482 downreg)。
3. **overdraft 不抽 storage 越冬底**(或两库等比例,不先榨茎库)。

**对照 CLM5**:CLM5 有 availc(可分配碳上限)+ excess_cflux(分不掉的不透支)+ downreg(碳限制下调)+
xsmrpool(呼吸债池,不碰 storage)。COMPASS 缺 availc 上限 + storage 被算进呼吸/overdraft 燃料。
**注**:用户指出 xsmrpool"扣负等光合补"非严格生理;"缺碳抽 storage"反而真实(储备应急),
故修法重点是**给 storage 留越冬保护底 + cap 不算 storage 让降呼吸接管**,而非照搬 xsmrpool 赊账。

### 8.5 处理决定
- **现状是 COMPASS 既有行为**(单池 v53 也"夏季抽 storage→秋季回填"),year-1 能 work。
- **双池让它现形**(茎库小,抽干扎眼)。
- **待 v58 跑完看 year-2**:若 year-2 靠根库 bootstrap + 回填能活 → 现状可暂留,§8.4 作下一轮治本;
  若 year-2 死或茎库不回填 → 立即上 §8.4(cap 不算 storage + 越冬保护底)。

---

*创建 2026-06-06,统一重写;§7 实施记录 + §8 碳路径诊断 2026-06-06 追加。*
*来源:多 agent 审计 + 全24PFT矩阵 + 验证workflow + 用户逐层追问(labile/storage 碳路径)。*
*("不要拼接/全PFT查/主用删另一个/先验证再动手/按来源归库非固定比例/C-N必须一致/分配要有度")。*
