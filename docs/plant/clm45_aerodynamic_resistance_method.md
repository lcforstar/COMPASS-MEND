# CLM4.5 改进空气动力阻力与土壤热边界模块方法说明

> 本文档整理当前 `ebal_clm45` 分支相对原始 `ebal` 的最新改动，面向论文方法部分、补充材料和代码维护。
> 组织顺序为：物理动机 -> 运行逻辑 -> 公式 -> 参数与输入 -> 与原版对比验证。
> 本方案不是纯 CLM4.5 big-leaf 复刻，而是保留 SCOPE/Wallace-Verhoef 多通道阻力网络，并注入 CLM4.5/CLM5 中更物理的稳定度、粗糙度、边界层和土壤热参数化。

---

## 1 模型总体框架

原始 SCOPE `ebal` 的能量平衡使用 Wallace-Verhoef 型阻力网络，将冠层和土壤分为两条通道：

```text
冠层叶片 -> 冠层内阻力 rawc -> 冠层上方阻力 raa -> 大气
土壤表面 -> 土壤通道阻力 raws -> 冠层上方阻力 raa -> 大气
```

其中

$$
raw_c = r_{wc} + r_{bc}
\tag{1}
$$

$$
raw_s = r_{ws} + r_{bs}
\tag{2}
$$

原始结构的优点是与 SCOPE 的多层辐射传输、sunlit/shaded 冠层和土壤-冠层通量分流一致；缺点是部分空气动力学参数静态化或经验化较强。本研究新增 `ebal_clm45` 链路，在不破坏 SCOPE 多层结构的前提下进行以下改进：

1. 保留 Wallace-Verhoef 多通道网络；
2. 将上方空气动力阻力改为 CLM4.5 风格的稳定度函数和 Stanton 修正；
3. 将冠层叶片边界层阻力 `rbc` 改为 Choudhury-Monteith 形式；
4. 将土壤表面阻力 `rss` 改为土壤含水量动态函数；
5. 将土壤热红外发射率 `epsilon_s` 改为表层饱和度动态函数；
6. 将土壤边界层阻力 `rbs` 改为随当前 LAI 动态缩放；
7. 引入 SAI，并用 $LAI+SAI$ 动态计算动量粗糙度 `z0m` 和零平面位移 `d`；
8. 在 `AMainloop` 中同时运行原始 `ebal` 和新 `ebal_clm45`，输出逐时对比表。

新增链路核心文件为：

| 文件 | 作用 |
|---|---|
| `src/ebal_clm45.m` | 新能量平衡主函数，接入动态阻力、热红外发射率和对比输出 |
| `src/resistances_clm45.m` | 新空气动力阻力计算函数 |
| `src/heatfluxes_clm45.m` | 支持热通道和水汽通道分离的通量计算 |
| `src/AMainloop.m` | 同步运行 `ebal` 与 `ebal_clm45` 并输出 `EBAL_compare.csv` |
| `src/+core/VegetationData.m` | 新增 SAI 配置与 InPlant 第 6 列读取 |
| `data/config_vegetation.txt` | 新增常数 SAI 输入 |

---

## 2 运行逻辑

### 2.1 主循环中的双方案对比

每个时间步 `KT` 中，主循环先运行原始能量平衡：

```matlab
[iter, fluxes, rad, thermal, ...] = ebal(...)
```

随后使用同一组输入运行新方案：

```matlab
[iter_clm45, fluxes_clm45, rad_clm45, thermal_clm45, ...] = ebal_clm45(...)
```

两套结果写入 `EBAL_compare`，并最终输出：

```text
DataOutput/EBAL_compare.csv
```

对比变量包括：

```text
Htot, lEtot, Hctot, lEctot, Hstot, lEstot
raa, rah, raw, rbc, rbs, rss
ustar, ustar_clm
Tsave, Tcave
SAI, L_plus_S, zo, d
epsilon_s, rs_thermal, Sat_surf
```

### 2.2 `ebal_clm45` 内部顺序

`ebal_clm45` 每个时间步进入能量平衡迭代前，先更新慢变量和边界参数：

```text
表层土壤含水量 Theta_surf
  -> 动态 rss
  -> 动态 epsilon_s / rs_thermal
  -> 动态 rbs
  -> SAI + LAI
  -> 动态 z0m 和 d
  -> resistances_clm45
  -> heatfluxes_clm45
  -> 能量平衡迭代
```

能量平衡迭代内仍保留 SCOPE 原始结构：

```text
RTMt_sb 热红外辐射传输
  -> 冠层/土壤净辐射
  -> 冠层/土壤显热与潜热
  -> 土壤温度和叶温更新
  -> 收敛判断
```

---

## 3 CLM4.5 风格空气动力阻力

### 3.1 保留 SCOPE/Wallace-Verhoef 多通道网络

纯 CLM4.5 陆面模型多采用 big-leaf 或较低维度的阻力表达。SCOPE 的辐射、光合和荧光计算是多层结构，并区分 sunlit/shaded 叶片。因此本方案没有将冠层压缩为单层，而是保留：

$$
r_{aa} = r_{ai} + r_{ar} + r_{ac}
\tag{3}
$$

$$
raw_c = r_{wc} + r_{bc}
\tag{4}
$$

$$
raw_s = r_{ws} + r_{bs}
\tag{5}
$$

其中 `rawc` 对应冠层通道，`raws` 对应土壤通道。这样可以与 SCOPE 中的通量分解保持一致：

$$
H_c = LAI \left(F_c H_{ch} + \overline{H_{cu}}\right)
\tag{6}
$$

$$
H_s = F_s H_s
\tag{7}
$$

潜热通量同理。这一结构来自 SCOPE 原始 Wallace-Verhoef 阻力网络。

### 3.2 四段稳定度函数

原始 W&V 实现中部分位置使用两段 Paulson 稳定度修正。本方案在 `resistances_clm45.m` 中引入 CLM4.5 风格四段稳定度函数，用于动量和热量稳定度修正：

$$
\psi_m = \psi_m(\zeta), \quad \psi_h = \psi_h(\zeta)
\tag{8}
$$

其中

$$
\zeta = \frac{z-d}{L}
\tag{9}
$$

$L$ 为 Monin-Obukhov 长度。四段函数区分：

```text
强不稳定
弱不稳定
弱稳定
强稳定
```

这样可以避免两段 Paulson 形式在强稳定或强不稳定条件下的过度外推。

来源：CLM4.5 Technical Note 第 5 章关于 Monin-Obukhov 相似理论的空气动力阻力参数化；Paulson 稳定度函数；Businger-Dyer 型稳定度函数。

### 3.3 CLM4.5 风格 `rah/raw`

新方案将冠层上方热量和水汽空气动力阻力显式分离：

$$
r_{ah} =
\frac{
\left[\ln\left(\frac{z-d}{z_{0m}}\right)-\psi_m(z)+\psi_m(z_{0m})\right]
\left[\ln\left(\frac{z-d}{z_{0h}}\right)-\psi_h(z)+\psi_h(z_{0h})\right]
}
{\kappa^2 U}
\tag{10}
$$

水汽通道同理：

$$
r_{aw} =
\frac{
\left[\ln\left(\frac{z-d}{z_{0m}}\right)-\psi_m(z)+\psi_m(z_{0m})\right]
\left[\ln\left(\frac{z-d}{z_{0w}}\right)-\psi_h(z)+\psi_h(z_{0w})\right]
}
{\kappa^2 U}
\tag{11}
$$

其中 $U$ 为参考高度风速，$\kappa$ 为 von Karman 常数。代码中将 $r_{ah}$ 输出为 `thermal.rah`，将 $r_{aw}$ 输出为 `thermal.raw`。

### 3.4 Stanton 修正的热粗糙度

热粗糙度使用 CLM4.5 风格 Stanton 修正：

$$
z_{0h} = z_{0m}
\exp\left[-0.13 \left(\frac{u_* z_{0m}}{\nu}\right)^{0.45}\right]
\tag{12}
$$

水汽粗糙度当前取：

$$
z_{0w} = z_{0h}
\tag{13}
$$

其中 $\nu$ 为空气运动黏度。该修正表达了动量粗糙度和标量粗糙度不同的物理事实，避免简单令 $z_{0h}=z_{0m}$ 导致标量交换过强。

来源：CLM4.5 Technical Note 第 5 章空气动力阻力与粗糙度参数化。

---

### 3.5 当前仍保留并实际使用的 SCOPE/Wallace-Verhoef 阻力路径

当前 `ebal_clm45` 并没有把 SCOPE/Wallace-Verhoef 阻力体系完全替换掉，而是保留其冠层内部和土壤穿冠层交换部分，再将上方空气动力交换替换为 CLM-style 的 `rah/raw`。因此当前方案应理解为：

```text
SCOPE/W&V 冠层内部阻力
+ CLM-style 上方热量/水汽空气动力阻力
```

在 `resistances_clm45.m` 中，冠层风速衰减系数为：

$$
n = \frac{C_d LAI}{2\kappa^2}
\tag{14a}
$$

其中 $C_d$ 为冠层阻力系数，$LAI$ 为叶面积指数，$\kappa$ 为 von Karman 常数。$n$ 越大，表示冠层内部风速衰减越强。

原 SCOPE/W&V 摩阻速度仍被计算为：

$$
u_* =
\frac{\kappa U}
{\ln \left( \frac{z-d}{z_{0m}} \right) - \psi_m}
\tag{14b}
$$

该 $u_*$ 主要用于冠层内部扩散、冠层高度风速和 W&V 冠层内阻力计算。由此计算热量湍流扩散系数：

$$
K_h =
\frac{\kappa u_* (z_r-d)}{\phi_h}
\tag{14c}
$$

其中 $z_r=2.5h$ 为粗糙子层顶部高度，$h$ 为冠层高度。

W&V 网络中仍保留以下阻力：

$$
r_{ai} =
\frac{1}{\kappa u_*}
\left[
\ln\left(\frac{z-d}{z_r-d}\right)
-\psi_h(z)+\psi_h(z_r)
\right]
\tag{14d}
$$

$$
r_{ar} =
\frac{1}{\kappa u_*}
\left(\frac{z_r-h}{z_r-d}\right)
-\psi^*_h(z_r)+\psi^*_h(h)
\tag{14e}
$$

$$
r_{ac} =
\frac{h\sinh(n)}{nK_h}
\left[
\ln\left(\frac{e^n-1}{e^n+1}\right)
-
\ln\left(
\frac{e^{n(z_{0m}+d)/h}-1}
{e^{n(z_{0m}+d)/h}+1}
\right)
\right]
\tag{14f}
$$

$$
r_{ws} =
\frac{h\sinh(n)}{nK_h}
\left[
\ln\left(
\frac{e^{n(z_{0m}+d)/h}-1}
{e^{n(z_{0m}+d)/h}+1}
\right)
-
\ln\left(
\frac{e^{0.01n/h}-1}
{e^{0.01n/h}+1}
\right)
\right]
\tag{14g}
$$

这些量的物理含义为：

| 阻力 | 含义 | 当前用途 |
|---|---|---|
| `rai` | 测风高度到粗糙子层顶部的惯性子层阻力 | 组成诊断量 `raa` |
| `rar` | 粗糙子层顶部到冠层顶部的阻力 | 组成诊断量 `raa` |
| `rac` | 冠层内部空气阻力 | 组成诊断量 `raa` |
| `rws` | 土壤表面到冠层空气空间的空气动力阻力 | 实际进入 `raws` |
| `rbc` | 叶片边界层阻力 | 实际进入 `rawc` |

其中：

$$
r_{aa}=r_{ai}+r_{ar}+r_{ac}
\tag{14h}
$$

在当前 `ebal_clm45` 中主要作为诊断量保存为 `thermal.raa`，并不再直接进入 `heatfluxes_clm45` 的显热和潜热通量计算。

实际仍参与通量计算的 SCOPE/W&V 阻力是：

$$
raw_c = r_{wc}+r_{bc}
\tag{14i}
$$

$$
raw_s = r_{ws}+r_{bs}
\tag{14j}
$$

其中 `rawc` 表示叶片到冠层空气空间的冠层内部路径阻力，`raws` 表示土壤表面到冠层空气空间的路径阻力。它们随后与 CLM-style 上方空气动力阻力组合：

$$
r_{a,leaf,h} = (LAI+1)(r_{ah}+raw_c)
\tag{14k}
$$

$$
r_{a,leaf,w} = (LAI+1)(r_{aw}+raw_c)
\tag{14l}
$$

$$
r_{a,soil,h} = (LAI+1)(r_{ah}+raw_s)
\tag{14m}
$$

$$
r_{a,soil,w} = (LAI+1)(r_{aw}+raw_s)
\tag{14n}
$$

因此当前热量和水汽路径可概括为：

```text
冠层显热:  leaf temperature -> rawc -> rah -> atmosphere
冠层潜热:  leaf vapor -> stomatal resistance -> rawc -> raw -> atmosphere
土壤显热:  soil surface -> raws -> rah -> atmosphere
土壤潜热:  soil vapor -> soil surface resistance -> raws -> raw -> atmosphere
```

对应 `heatfluxes_clm45` 中的通量表达为：

$$
H = \frac{\rho c_p}{r_{a,h}}(T_c-T_a)
\tag{14o}
$$

$$
\lambda E =
\frac{\rho\lambda}{r_{a,w}+r_s}(q_i-q_a)
\tag{14p}
$$

这里的优点是：SCOPE/W&V 的冠层内部和土壤穿冠层交换结构仍被保留，因而可以继续服务于 sunlit/shaded 叶片、冠层辐射、光合和植物水力耦合；同时，冠层上方与大气之间的交换由 CLM-style `rah/raw` 控制，使热量和水汽空气动力阻力可以分离，并响应 Monin-Obukhov 稳定度。

需要注意的是，该结构不是纯 CLM canopy flux network。它是：

```text
SCOPE two-leaf energy balance framework
+ SCOPE/W&V within-canopy resistance
+ CLM-style atmospheric stability and scalar aerodynamic resistance
```

因此 `raa` 没有在当前通量路径中与 `rah/raw` 重复相加；实际进入通量的是 `rawc/raws + rah/raw` 的混合路径。

### 3.6 当前需要诊断的物理一致性问题：`ustar` 与 `ustar_clm`

与 CTSM/CLM5 `release-clm5.0.37` 源代码对比后，当前 `ebal_clm45` 的能量平衡符号方向没有发现明显错误。当前约定为：

```text
H  > 0: 地表/冠层向大气的显热通量
LE > 0: 地表/冠层向大气的潜热通量
G  > 0: 进入土壤的热通量
```

因此地表能量平衡残差：

$$
R_n - LE - H - G
\tag{14q}
$$

与 CLM 地表能量收支方向基本一致。

目前更需要诊断的是摩阻速度与阻力闭环的一致性。CLM 中 `ustar`、稳定度函数、`rah/raw` 和 Monin-Obukhov 长度属于同一套迭代闭环：

```text
L -> stability functions -> ustar -> rah/raw -> H/LE -> new L
```

当前本地代码中存在两个摩阻速度：

| 变量 | 来源 | 当前用途 |
|---|---|---|
| `ustar` | SCOPE/W&V 风速廓线 | 用于冠层内部阻力、`Kh`、`rawc/raws`，并用于更新 Monin-Obukhov 长度 `L` |
| `ustar_clm` | CLM-style `fm` 公式 | 用于计算 scalar roughness `z0h/z0w`，并参与 `rah/raw` 的 CLM-style 阻力体系 |

本地当前闭环可概括为：

```text
L -> CLM-style stability functions -> ustar_clm -> rah/raw -> H/LE
H/LE -> use SCOPE/W&V ustar -> new L
```

这不是阻力重复计算，也不是符号错误，而是一个混合方案中的物理一致性风险点。若 `ustar` 与 `ustar_clm` 在大部分时段接近，则该混合闭环影响较小；若两者在稳定层结、强不稳定或低风速条件下长期差异较大，则 `L` 的反馈可能与实际用于通量的 `rah/raw` 不完全自洽。

建议将该项作为诊断输出和后续验证重点，而不是立即修改。推荐检查：

```text
thermal.ustar
thermal.ustar_clm
thermal.rah
thermal.raw
thermal.L
fluxes.Htot
fluxes.lEtot
```

推荐诊断指标为：

$$
R_{u_*} = \frac{u_{*,clm}}{u_*}
\tag{14r}
$$

若 $R_{u_*}$ 长期接近 1，说明当前混合阻力闭环基本稳定；若 $R_{u_*}$ 经常大于 2 或小于 0.5，则应优先检查稳定度反馈、低风速下限和 `rah/raw` 对通量的敏感性。

---

## 4 冠层叶片边界层阻力 `rbc`

### 4.1 原始公式

原始 W&V 经验公式为：

$$
r_{bc,old} = \frac{70}{LAI}\sqrt{\frac{w}{u_{z0}}}
\tag{14}
$$

其中 $w$ 为叶宽，$u_{z0}$ 为近冠层粗糙高度风速。该公式简单、稳定，但对冠层高度、叶宽动态和 LAI 上限响应不足。

### 4.2 Choudhury-Monteith 型公式

新方案采用用户 `untitled7.m` 和 `Calc_Resistance_Full_Optimized.m` 中的 Choudhury-Monteith 型叶片边界层阻力。首先计算冠层高度风速：

$$
U_h = \frac{u_*}{\kappa}\ln\left(\frac{h-d}{z_{0m}}\right)
\tag{15}
$$

风速衰减参数按冠层高度分段：

$$
n =
\begin{cases}
2.5, & h \le 1 \\
2.306 + 0.194h, & 1 < h < 10 \\
4.25, & h \ge 10
\end{cases}
\tag{16}
$$

叶片边界层阻力：

$$
r_b =
\frac{100}{n}
\left(\frac{w_{leaf}}{U_h}\right)^{1/2}
\frac{1}{1-\exp(-n/2)}
\tag{17}
$$

再将其按叶面积归一化为冠层通道阻力：

$$
r_{bc} =
\frac{r_b}{D_{LAI}}
\tag{18}
$$

其中

$$
D_{LAI} =
\begin{cases}
2LAI, & LAI \le 0.5LAI_{max} \\
0.5LAI_{max}, & LAI > 0.5LAI_{max}
\end{cases}
\tag{19}
$$

该公式相比原式的主要改进是：显式响应 $h$、$d$、$z_{0m}$、$u_*$、叶宽、$LAI$ 和 $LAI_{max}$。

来源：Choudhury and Monteith canopy boundary-layer resistance framework；用户已有蒸散阻力函数 `Calc_Resistance_Full_Optimized.m`。

---

## 5 动态土壤表面阻力 `rss`

原始 SCOPE 可使用常数 `soil.rss`。新方案在 `ebal_clm45` 每步用表层含水量重算：

$$
r_{ss} = 38113 \exp(-12 \theta_s / \theta_{sat})
\tag{20}
$$

并限制在：

$$
50 \le r_{ss} \le 2000
\tag{21}
$$

其中 $\theta_s$ 为表层土壤体积含水量，$\theta_{sat}$ 为饱和含水量。该形式使土壤变干时水汽扩散阻力快速增大，土壤湿润时允许更强蒸发。

当前代码实现位置：

```text
src/+CalculateEvapotranspiration/Calculate_Soil_Resistance.m
src/ebal_clm45.m
```

来源：Mahfouf and Noilhan 型土壤蒸发阻力参数化；Sellers 陆面模型中类似土壤表面阻力表达。

---

## 6 动态土壤热红外发射率 `epsilon_s`

### 6.1 物理意义

土壤长波辐射发射能力为：

$$
L_{\uparrow} = \epsilon_s \sigma T_s^4
\tag{22}
$$

原始模型中热红外反射率 `soil.rs_thermal` 为固定值，因此：

$$
\epsilon_s = 1-rs_{thermal}
\tag{23}
$$

也是固定的。

### 6.2 表层饱和度动态化

新方案定义表层饱和度：

$$
Sat_s = \frac{\theta_s}{\theta_{sat}}
\tag{24}
$$

并截断到 $[0,1]$。热红外反射率按干湿两端线性插值：

$$
rs_{thermal}
= \rho_{dry} + (\rho_{wet}-\rho_{dry})Sat_s
\tag{25}
$$

其中当前默认：

$$
\rho_{dry}=0.086,\quad \rho_{wet}=0.020
\tag{26}
$$

因此

$$
\epsilon_s = 1-rs_{thermal}
= 1-\left[0.086 + (0.020-0.086)Sat_s\right]
\tag{27}
$$

即：

```text
干土 Sat=0：epsilon_s=0.914
中等湿度 Sat=0.5：epsilon_s=0.947
饱和 Sat=1：epsilon_s=0.980
```

该参数进入热红外辐射传输 `RTMt_sb`，并进入土壤温度更新的长波导数项：

$$
4\epsilon_s \sigma (T_s+273.15)^3
\tag{28}
$$

来源：Idso 型土壤热红外发射率与含水量关系；CLM5 动态土壤光热属性思想。需要说明的是，CLM4.5 Technical Note 中可见地表发射率 $\epsilon_g$，但土壤发射率常取固定值，本研究的式（25）是进一步动态化改进，不是 CLM4.5 原公式。

---

## 7 动态 `rbs`

### 7.1 问题

原始 `select_input.m` 中在初始化阶段计算：

$$
r_{bs} = r_{bs,ref}\frac{LAI}{4.3}
\tag{29}
$$

若之后 LAI 随时间变化，而 `rbs` 不再更新，就会出现幼苗期和盛花期使用同一土壤边界层阻力的问题。

### 7.2 新方案

为避免重复缩放，`ebal_clm45` 首先恢复参考值：

$$
r_{bs,ref} = r_{bs,input}\frac{4.3}{LAI}
\tag{30}
$$

然后每步用当前 LAI 重新计算：

$$
r_{bs,dyn} = r_{bs,ref}\frac{LAI}{4.3}
\tag{31}
$$

这样在当前时间步与原始初始化一致，但在 LAI 随季节变化时可以同步更新。

代码输出：

```text
thermal.rbs
thermal.rbs_static
thermal.rbs_init
```

验证时应满足：

$$
r_{bs,dyn} = r_{bs,init}\frac{LAI}{4.3}
\tag{32}
$$

来源：SCOPE 原始 `equations.calc_rssrbs` 中的 LAI 缩放逻辑。

---

## 8 SAI 输入与 CLM5 动态粗糙度

### 8.1 SAI 输入逻辑

SAI 表示茎面积指数，用于与 LAI 共同描述冠层有效面积：

$$
L+S = LAI + SAI
\tag{33}
$$

当前输入优先级为：

```text
InPlant 第 6 列 SAI > config_vegetation.txt 中的 SAI 常数 > 默认 0.05
```

`config_vegetation.txt` 中可写：

```text
SAI=0.05
```

如果使用 Excel 动态植被输入，`InPlant.xlsx/InPlant.xls` 列顺序为：

```text
DOY | LAI | hc | Vcmax(optional) | Cab(optional) | SAI(optional)
```

如果第 6 列不存在，则全时段使用 `SAI=0.05`。

### 8.2 CLM5 动态粗糙度公式

原始 `canopy.zo` 和 `canopy.d` 可能在初始化后冻结。新方案每步用当前 `LAI`、`SAI` 和冠层高度 `hc` 计算：

$$
z_{0m}
= z_{0g} + (0.13h_c - z_{0g})
\left[1-\exp\left(-0.5(LAI+SAI)\right)\right]
\tag{34}
$$

$$
d =
\left[
1-
\frac{1-\exp\left(-\sqrt{7.5(LAI+SAI)}\right)}
{\sqrt{7.5(LAI+SAI)}}
\right]h_c
\tag{35}
$$

当前默认裸土动量粗糙度：

$$
z_{0g}=0.01\ \text{m}
\tag{36}
$$

公式物理含义：

```text
LAI+SAI 很小时，z0m 接近裸土粗糙度；
LAI+SAI 增大时，z0m 渐近到 0.13hc；
d 随冠层变密变高而增大，但应保持 d < hc。
```

来源：CLM5 Technical Note 粗糙度和零平面位移参数化，用户给定为 Eq. 5.127-5.128。该部分属于 CLM5 动态粗糙度思想，而不是 CLM4.5 Tech Note 中的固定植被粗糙度处理。

---

## 9 热通量与水汽通量分离

原始通量计算中热量和水汽常共用同一个空气动力阻力。新 `heatfluxes_clm45` 支持：

$$
H = \frac{\rho c_p}{r_{a,h}}(T_c-T_a)
\tag{37}
$$

$$
\lambda E =
\frac{\rho}{r_{a,w}+r_s}
\lambda(q_i-q_a)
\tag{38}
$$

当前实现已将 $\rho$ 从固定标准空气密度改为随 forcing 状态变化的空气密度：

$$
q_a = e_a \frac{100M_w}{M_a p}
\tag{38a}
$$

$$
T_v = (T_a+273.15)(1+0.61q_a)
\tag{38b}
$$

$$
\rho_a = \frac{p}{R_d T_v}
\tag{38c}
$$

其中 $p$ 为气压，$e_a$ 为空气水汽压，$q_a$ 为空气比湿近似，$T_v$ 为虚温。该处理对应 CLM 中 `forc_rho` 进入湍流交换系数的思想，使 $H$ 和 $\lambda E$ 能随气压、气温和湿度变化而改变，而不是固定使用标准状态空气密度。代码中该诊断量输出为 `thermal.rhoa`。

其中 $r_{a,h}$ 为热通道阻力，$r_{a,w}$ 为水汽通道阻力。当前仍使用温度差 $T_c-T_a$，尚未切换到位温或虚位温。此前评估表明，该项对当前低高度农田情景影响通常小于阻力、`rss`、`z0m/d` 等主要项。

---

## 10 参数与输出

### 10.1 新增输入参数

| 参数 | 位置 | 默认值 | 含义 |
|---|---|---:|---|
| `SAI` | `data/config_vegetation.txt` | 0.05 | 茎面积指数常数 |
| `SAI` 第 6 列 | `InPlant.xlsx/InPlant.xls` | 可选 | SAI 时间序列，优先于配置常数 |
| `Rho_Thermal_Soil_Dry` | `cfg.Env` | 0.086 | 干土热红外反射率 |
| `Rho_Thermal_Soil_Wet` | `cfg.Env` | 0.020 | 湿土热红外反射率 |

### 10.2 新增诊断输出

| 输出 | 含义 |
|---|---|
| `thermal.rah` | CLM4.5 风格热通道空气动力阻力 |
| `thermal.raw` | CLM4.5 风格水汽通道空气动力阻力 |
| `thermal.rhoa` | forcing 状态空气密度，用于热量和水汽湍流通量 |
| `thermal.rbc` | 新 Choudhury-Monteith 冠层边界层阻力 |
| `thermal.rbc_wv` | 原 W&V 经验 `rbc`，用于对照 |
| `thermal.rbs` | 动态土壤边界层阻力 |
| `thermal.rbs_init` | 恢复出的参考 `rbs` |
| `thermal.rss` | 动态土壤表面水汽阻力 |
| `thermal.Sat_surf` | 表层土壤饱和度 |
| `thermal.epsilon_s` | 动态土壤热红外发射率 |
| `thermal.rs_thermal` | 动态土壤热红外反射率 |
| `thermal.SAI` | 当前 SAI |
| `thermal.L_plus_S` | 当前 $LAI+SAI$ |
| `thermal.zo` | 动态动量粗糙度 |
| `thermal.d` | 动态零平面位移 |
| `thermal.ustar_clm` | CLM4.5 风格摩擦速度 |

---

## 11 与原版的数值验证

### 11.1 最终版本单时刻对比

以下为同一组输入下，原始 `ebal` 与最终 `ebal_clm45` 在 `KT=1` 的对比输出：

| 变量 | 原版 | CLM45 改进版 | 差值 |
|---|---:|---:|---:|
| $H_{tot}$ (W m⁻²) | -8.61 | -6.03 | +2.57 |
| $\lambda E_{tot}$ (W m⁻²) | 7.89 | 10.72 | +2.83 |
| $H_c$ (W m⁻²) | -19.25 | -13.78 | +5.47 |
| $H_s$ (W m⁻²) | 10.64 | 7.75 | -2.89 |
| $\lambda E_c$ (W m⁻²) | 0.84 | 0.23 | -0.61 |
| $\lambda E_s$ (W m⁻²) | 7.05 | 10.49 | +3.44 |
| $T_{s,ave}$ (°C) | 15.84 | 15.80 | -0.04 |
| $T_{c,ave}$ (°C) | 11.51 | 11.01 | -0.50 |

阻力和粗糙度诊断：

| 变量 | 结果 | 物理解释 |
|---|---:|---|
| `raa` 原版 | 20.70 s m⁻¹ | 原 W&V 冠层上方阻力 |
| `raa` CLM45/W&V 网络 | 18.10 s m⁻¹ | 动态 `z0m/d` 后冠层交换增强 |
| `rah/raw` | 68.22 / 68.22 s m⁻¹ | CLM4.5 风格热/水汽通道阻力 |
| `rbc` 新/旧 | 8.19 / 12.38 s m⁻¹ | Choudhury-Monteith 公式降低叶片边界层阻力 |
| `rbs` 动态/初始缩放/参考 | 6.98 / 6.98 / 10.00 s m⁻¹ | 验证 $10\times LAI/4.3$，未重复缩放 |
| `rss` | 200.67 s m⁻¹ | 中等湿度土壤表面阻力 |
| `Sat_surf` | 0.437 | 表层土壤中等偏湿 |
| `epsilon_s` | 0.943 | 介于干土 0.914 和湿土 0.980 之间 |
| `rs_thermal` | 0.057 | 与 $1-\epsilon_s$ 一致 |
| `SAI/L+S` | 0.050 / 3.050 | 使用配置常数 SAI |
| `zo` 动态/旧 | 0.0530 / 0.0292 m | 动态粗糙度增强动量交换 |
| `d` 动态/旧 | 0.3963 / 0.4104 m | 零平面位移仍小于冠层高度 |
| `u*` 原版/CLM45 | 0.1953 / 0.2239 m s⁻¹ | 动态粗糙度使摩擦速度升高 |

### 11.2 公式一致性验证

**1. 动态土壤热红外发射率**

诊断输出：

```text
Sat_surf = 0.437
rs_th = 0.057
epsilon_s = 0.943
```

代入式（25）：

$$
rs_{thermal}
= 0.086 + (0.020-0.086)\times 0.437
= 0.057
\tag{39}
$$

因此：

$$
\epsilon_s = 1-0.057 = 0.943
\tag{40}
$$

与模型输出完全一致，说明动态发射率接入正确。

**2. 动态 `rbs`**

诊断输出：

```text
rbs_dynamic = 6.98
rbs_init = 10.00
```

当前 $LAI=3.0$，由式（31）：

$$
r_{bs,dyn}=10.00\times \frac{3.0}{4.3}=6.98\ \text{s m}^{-1}
\tag{41}
$$

与输出一致，说明没有发生重复缩放。

**3. 动态粗糙度**

诊断输出：

```text
SAI = 0.050
LAI+SAI = 3.050
zo = 0.0530 m
d = 0.3963 m
```

当前 $LAI+SAI=3.05$，动态 `zo` 大于旧值 `0.0292 m`，代表冠层动量粗糙度增强；`d=0.3963 m` 仍小于冠层高度，满足物理约束。对应摩擦速度由约 `0.1953` 增至 `0.2239 m s⁻¹`，空气动力交换增强，`rah/raw` 从先前约 `77.55 s m⁻¹` 降至 `68.22 s m⁻¹`。

### 11.3 能量分配响应验证

动态 `z0m/d` 和阻力改进后，土壤通道表现为：

```text
Hs: 10.64 -> 7.75 W m⁻²
LEs: 7.05 -> 10.49 W m⁻²
```

在当前表层土壤饱和度 `0.437`、`rss=200.67 s m⁻¹` 的条件下，土壤处于中等湿润状态。新方案增强空气动力交换，同时允许较强土壤蒸发，因此土壤可用能量从显热更多转向潜热。这一方向符合物理预期。

冠层通道表现为：

```text
Hc: -19.25 -> -13.78 W m⁻²
LEc: 0.84 -> 0.23 W m⁻²
Tcave: 11.51 -> 11.01 °C
```

冠层平均温度下降约 `0.50 °C`，说明新边界层阻力与上方空气动力阻力增强了叶片与空气的耦合。由于当前时刻光合和蒸腾较弱，冠层潜热变化较小。

土壤平均温度：

```text
Tsave: 15.84 -> 15.80 °C
```

变化仅 `-0.04 °C`，说明动态 `epsilon_s` 和动态空气阻力没有引起土壤温度数值漂移，能量平衡仍稳定。

---

## 12 与原始 SCOPE 和纯 CLM4.5 的关系

本方案不是纯 CLM4.5 big-leaf 替换。原因是 SCOPE 的核心优势是多层辐射传输、sunlit/shaded 光合和荧光计算，若直接使用纯 CLM4.5 单层阻力，会破坏以下对应关系：

| 模块 | SCOPE 需求 | 本方案处理 |
|---|---|---|
| 冠层辐射 | 多层、sunlit/shaded | 保留 W&V 多通道网络 |
| 冠层通量 | 冠层/土壤分流 | 保留 `rawc/raws` |
| 光合与荧光 | 多层叶片输出 | 不降维为 big-leaf |
| 空气动力学 | 需更物理的稳定度和粗糙度 | 注入 CLM4.5/CLM5 参数化 |
| 土壤蒸发 | 需响应表层水分 | 动态 `rss` 与 `epsilon_s` |

因此本方案可概括为：

```text
SCOPE/W&V 多层阻力骨架
+ CLM4.5 四段稳定度与 Stanton 标量粗糙度
+ Choudhury-Monteith 冠层叶片边界层阻力
+ CLM5 LAI+SAI 动态 z0m/d
+ 动态 rss 与 epsilon_s
```

---

## 13 当前仍未改动的部分

以下部分目前保留，作为后续可选工作：

1. **`L` 仍由能量平衡外迭代反算**  
   尚未改为基于虚位温 Richardson 数的一次性内算。该改动风险较高，且当前验证显示主要误差已由 `z0m/d`、`rss`、`rbc` 等项控制，暂不优先。

2. **`rwc` 暂未动态化**  
   `rawc = rwc + rbc` 中 `rbc` 已动态化。`rwc` 缺少明确的本地公式来源，直接按 LAI 缩放可能会误调冠层通道，因此保留常数，后续根据 `rawc` 季节响应再决定。

3. **热通量仍用温度差而非位温差**  
   当前冠层高度和参考高度下，位温修正量级较小。若后续扩展到高塔、山地或强稳定边界层，可进一步引入 $\theta$ 或 $\theta_v$。

---

## 14 结论

当前 `ebal_clm45` 改进链路通过同一输入下与原始 `ebal` 的并行对比证明：

1. 动态 `epsilon_s` 的数值严格落在干湿土物理范围内，并与表层饱和度公式一致；
2. 动态 `rbs` 精确满足原 SCOPE 的 LAI 缩放公式，且避免重复缩放；
3. CLM5 动态 `z0m/d` 使 `ustar` 合理增大、`rah/raw` 合理降低；
4. 冠层和土壤通量响应方向与物理预期一致；
5. 土壤平均温度没有异常漂移，说明新阻力和热红外参数未破坏能量平衡稳定性；
6. 新方案保留 SCOPE 多层结构，避免纯 CLM4.5 big-leaf 降维造成的辐射-能量-光合耦合损失。

最终建议：将 `ebal_clm45` 作为当前候选改进方案，继续用全天或全季节观测的 `LE/H`、土壤温度和叶温进行定量筛选。若观测显示新版 `LEs` 偏高，再优先校准 `rss` 或土壤通道阻力，而不是先改 `L` 内算。

---

## 15 参考来源

1. Oleson, K. W., et al. CLM4.5 Technical Note. Chapter 5: Biogeophysics, aerodynamic resistance, Monin-Obukhov similarity, scalar roughness and surface flux parameterizations.
2. CLM5 Technical Note. Section on vegetation roughness length and displacement height, Eq. 5.127-5.128 for $z_{0m}$ and $d$ as functions of $LAI+SAI$.
3. Wallace, J. S., and Verhoef, A. Two-source canopy-soil resistance network used in SCOPE-style energy balance.
4. Choudhury, B. J., and Monteith, J. L. Canopy boundary-layer resistance formulation for vegetation-atmosphere exchange.
5. Mahfouf, J.-F., and Noilhan, J. Soil resistance parameterization for bare soil evaporation.
6. Sellers, P. J., et al. Land surface model soil evaporation resistance concepts.
7. Idso, S. B. Soil thermal infrared emissivity dependence on soil moisture.
8. SCOPE/STEMMUS-SCOPE local implementation files: `ebal.m`, `resistances.m`, `RTMt_sb.m`, `equations.calc_rssrbs`, and current `ebal_clm45` branch.
