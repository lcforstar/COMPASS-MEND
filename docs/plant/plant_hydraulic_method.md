# 植物水力耦合模块方法说明

> 本文档面向论文方法部分写作，同时兼顾内部理解。
> 全文以"物理动机 → 数学表达 → 参数说明"的顺序组织，不出现实现细节代码。
> 最后一节给出模块级数值验证结果，可直接引用为补充材料图表。

---

## 1 模型总体框架

本模型将土壤-植物-大气连续体（SPAC）简化为**一条共享上游通路、两个并联叶端分支**的拓扑结构（图 1）：

```
土壤各层 ──→ 根系 ──→ 茎端共享节点 ──→ 阳叶
                                    └──→ 阴叶
```

根据此结构，阳叶和阴叶共享同一个根-茎水力通道，不能各自独立地从土壤取用全量水分。SCOPE 能量平衡模块分别计算阳叶和阴叶的**蒸腾需求**，植物水力模块再判断共享的根-茎系统能够实际供应多少水，并将实际供水量按需求比例反分配给阳叶和阴叶。

当前计算逻辑分以下步骤依次完成：

1. 初始化阶段读取根系几何。若存在实测根长密度 $B$，则优先使用 measured-B，并由 $RAI=2\pi r_{root}Bz_{sup}$ 推出端点根面积；若无实测 $B$，则使用模型给定的 $RAI$，并反推等效根长密度。
2. 土壤水分模块逐时步更新端点土壤导水率 $K_{soil}$ 和土壤水势；植物水力模块将 $K_{soil}$ 转为 SI 单位，并用 CONTACT 因子得到 $K_{eff}$。
3. 根据 `gsr_formula` 选择土壤-根导度模式：`katul`、`clm5` 或 `radial_darcy_series`。其中 `clm5` 和 `radial_darcy_series` 都采用 soil-side 与 root-side 串联阻力结构；`radial_darcy_series` 进一步用圆柱根径向 Darcy 积分计算 soil-side conductance。
4. 将各有限元端点导度装配到共享土壤节点，得到节点导度 $g_{sr,n}$、总导度 $G_{sr}$，并用 $g_{sr,n}$ 加权计算根系感受到的有效土水势；该有效土水势包含从土壤节点到根颈/地表方向的重力修正。
5. 利用木质部脆弱性曲线积分，构建**共享供水能力**与候选总蒸腾量的函数关系。
6. 用**二分法**求解使供需平衡的实际总蒸腾量。
7. 按需求比例分配给阳叶/阴叶，并**反求各自叶水势**。
8. 将叶水势代入脆弱性曲线，得到**水力胁迫因子** $w_t$，反馈至气孔阻力和光合计算。

---

## 2 土壤-根导度

### 2.1 端点导度计算

土壤导水率 $K_{soil}$（单位 cm s$^{-1}$）由 STEMMUS 土壤水分模块逐时步更新，其维度为 $N_L \times 2$，对应每个土层的两个有限元端点。模型先将其转换为 SI 单位并施加根-土接触修正：

$$
K_{eff,e}=K_{soil,e}\times 10^{-2}\, f_{contact}(\psi_{soil,e})
\tag{1}
$$

其中 $e=(l,d)$，$l$ 为土层编号，$d \in \{1,2\}$ 为端点编号。随后根据配置的 `gsr_formula` 选择三种 soil-root conductance 之一。

**Katul-style bulk conductance**：

$$
g_{sr,e}^{Katul}
=
\frac{K_{eff,e}\sqrt{RAI_e}}
{\pi \rho_w g Z_r}
10^6
\tag{2a}
$$

该式为经验根区 bulk scaling，$\sqrt{RAI_e}$ 表示根面积对有效吸水界面的几何修正。它适合保留原始 macroscopic root-water uptake 行为，但不显式区分土壤侧径向阻力和根组织阻力。

**CLM5/CTSM-style series conductance**：

$$
dx_{root,e}=(\pi B_e)^{-1/2}
$$

$$
k_{s,e}^{CLM5}=\frac{K_{eff,e}}{dx_{root,e}}
$$

$$
k_{r,e}
=
f_{root,e}\frac{RAI_e K_{root,max}}{L_{path,e}}
$$

$$
g_{sr,e}^{CLM5}
=
\frac{k_{s,e}^{CLM5}k_{r,e}}
{k_{s,e}^{CLM5}+k_{r,e}}
\frac{10^6}{\rho_w g}
\tag{2b}
$$

其中 $B_e$ 为端点根长密度，$K_{root,max}$ 对应可选配置 `root_hydraulic_krmax`，单位按 CLM5 `krmax` 语义为 m s$^{-1}$，表示根组织导度的面积-路径长度尺度参数；$L_{path,e}$ 为根到根颈的路径长度近似，当前为 $0.25+z_e$ m；$f_{root,e}$ 为可选根段 PLC 因子，若无 root-specific vulnerability 参数则取 1。该式采用 CLM5/CTSM 的 soil-side 与 root-side 串联阻力结构。当 `root_hydraulic_krmax=0` 时，root tissue side 在 `gsr` 中视为非限制，避免把整株 `hydraulic_gpmax` 同时用于 soil-root conductance 和后续植物脆弱性曲线。

**Radial Darcy + root-side series conductance**：

$$
r_{outer,e}=\frac{1}{\sqrt{\pi B_e}}, \qquad
\Lambda_e=\ln\left(\frac{r_{outer,e}}{r_{root}}\right)
$$

$$
k_{s,e}^{radial}
=
\frac{2\pi K_{eff,e}B_e z_{sup,e}}{\Lambda_e}
=
\frac{K_{eff,e}RAI_e}{r_{root}\Lambda_e}
$$

$$
g_{sr,e}^{radial}
=
\frac{k_{s,e}^{radial}k_{r,e}}
{k_{s,e}^{radial}+k_{r,e}}
\frac{10^6}{\rho_w g}
\tag{2c}
$$

其中 $z_{sup,e}$ 为端点支撑厚度，当前对每层两个端点取 $z_{sup,e}=\Delta z_l/2$。若有实测根长密度，则 $B_e$ 直接来自观测，并由 $RAI_e=2\pi r_{root}B_ez_{sup,e}$ 推出端点根面积；若无实测根长密度，则以模型给定 $RAI_e$ 为主，并反推等效 $B_e$。这样保证根间距和吸水面积来自同一套根几何。

上述三种计算**仅对根区范围内** $RAI_e > 0$ 或等效根面积大于零的端点进行；非根区端点即使 $K_{soil}$ 非零，对植物吸水也无贡献。节点 $n$ 的总土壤-根导度为

$$
g_{sr,n} = \sum_{e:\,m(e)=n} g_{sr,e}
\tag{2d}
$$

其中 $m(e)$ 表示端点 $e$ 装配到的全局土壤节点。内部节点由相邻两个单元端点共同贡献，边界节点只接收一个相邻端点贡献。

### 2.2 根-土接触修正（CONTACT 机制）

土壤干燥过程中，根表面与土壤颗粒之间可能出现气隙，使根际实际导水能力低于体积土壤导水率所指示的水平。本模型采用 Lei 等（2023）提出的分段线性 CONTACT 因子 $f_{contact}$ 对端点导水率进行修正：

$$
f_{contact}(\psi) =
\begin{cases}
1 & \psi \geq \psi_{start} \\[4pt]
\dfrac{\psi - \psi_{end}}{\psi_{start} - \psi_{end}} & \psi_{end} < \psi < \psi_{start} \\[8pt]
0 & \psi \leq \psi_{end}
\end{cases}
\tag{3}
$$

其中 $\psi$ 为节点土水势（MPa，负值），$\psi_{start}$ 和 $\psi_{end}$ 分别为接触开始退化和完全丧失时对应的土水势阈值。本研究参照 Lei 等（2023）向日葵实验情景中使用的值，取 $\psi_{start} = -0.05$ MPa、$\psi_{end} = -1.5$ MPa，并通过配置文件向用户开放。这里的 $\psi_{end}$ 是当前参数化下的根-土接触完全丧失阈值，不是苹果根系的普适生理常数，应用到不同物种和土壤质地时必须重新校准。

将 $f_{contact}$ 代入 $K_{eff,e}$ 后，三种 `gsr` 分支中的端点等效土壤-根导度都会随土壤变干而降低，从而在不修改体积土壤水力模块的前提下，捕捉干旱过程中根-土界面阻力增大的物理效应。

### 2.3 总导度与根区加权土水势

对所有根区节点求和，得到总土壤-根导度：

$$
G_{sr} = \sum_{n=1}^{N_N} g_{sr,n}
\tag{4}
$$

$G_{sr}$ 越大，根系越容易从土壤取水；相同蒸腾通量下，根水势下降越小。

土壤节点水头 $h_n$（cmH₂O）换算为水势（MPa）：

$$
\psi_{soil,n} = h_n \times 9.807 \times 10^{-5}
\tag{5}
$$

以各节点土壤-根导度为权重，计算植物根系实际"感受到"的根区平均土水势。这里的权重不是单纯根量权重，而是由根几何、土壤导水率、根-土接触状态和所选 `gsr_formula` 共同决定的有效导度权重：

$$
\psi_s = \frac{\sum_n g_{sr,n} \cdot \psi_{soil,n}^{eff}}{G_{sr}}
\tag{6}
$$

其中 $\psi_{soil,n}^{eff}=\psi_{soil,n}-\psi_{gravity,n}$，表示扣除土壤节点到根颈/地表方向重力势后的有效土水势。该加权平均保证了**根量多且水力连通性强**的土层贡献更大；根多但土壤极干、导水率低或 CONTACT 因子接近零的土层不会仅因根量大而获得高权重。这体现的是多层土壤异质供水的**导度加权**图像，而不是单纯的总根量加权。

---

### 2.4 `gsr` 的三种土壤-根导度模式

当前 `Plantrootconductivity` 支持三种 soil-root conductance (`gsr`) 计算模式。三者的共同目标是把每个有限元端点的土壤水势与根水势联系起来：

$$
U_e = g_{sr,e}\left(\psi_{soil,e}-\psi_{root}\right)
$$

其中 $U_e$ 为端点吸水通量（m s$^{-1}$），$g_{sr,e}$ 为土壤-根导度（m s$^{-1}$ MPa$^{-1}$），$\psi_{soil,e}$ 和 $\psi_{root}$ 为 MPa 单位水势。若计入从土壤节点到根颈/地表的重力修正，则使用 $\psi_{soil,e}-\psi_{gravity,e}$ 作为有效土壤水势。

#### 2.4.1 Katul-style bulk conductance

默认经验分支采用 bulk root-zone length scale：

$$
g_{sr,e}^{Katul}
=
\frac{
K_{eff,e}\sqrt{RAI_e}
}{
\pi \rho_w g Z_r
}
10^6
$$

其中 $K_{eff,e}$ 为经 CONTACT 修正后的土壤导水率（m s$^{-1}$），$RAI_e$ 为端点根面积指数，$Z_r$ 为有效根区深度（m），$\rho_w$ 为水密度，$g$ 为重力加速度。该式保留 Lai and Katul (2000) 类型 macroscopic root-water uptake 思路：根区吸水能力由土壤导水率、根系空间分布和根区尺度共同控制。这里的 $\sqrt{RAI_e}$ 是经验几何缩放项，不是显式径向 Darcy 积分。

代码对应：

```matlab
layer_endpoint_gsr = K_eff_m_s * sqrt(root_area_endpoint(ML, ND)) ...
    / depth_factor * pa_per_mpa;
```

来源：Lai and Katul (2000) 提出的动态 root-water uptake 与潜在/实际蒸腾耦合框架；本实现为该类 macroscopic soil-root conductance 在当前有限元端点结构下的工程化形式。

#### 2.4.2 CLM5/CTSM-style series conductance

CLM5/CTSM plant hydraulics 将 soil-to-root pathway 表示为土壤-根界面阻力和根组织阻力串联。首先用根长密度 $B_e$ 估计根间距：

$$
dx_{root,e}=(\pi B_e)^{-1/2}
$$

土壤侧水头单位导度为：

$$
k_{s,e}^{CLM5}=\frac{K_{eff,e}}{dx_{root,e}}
=K_{eff,e}\sqrt{\pi B_e}
$$

根组织侧导度使用 whole-plant 最大导度 $g_{pmax}$ 转换到水头单位后，按端点根面积比例分配：

$$
k_{pmax}=\frac{g_{pmax}}{10^6/(\rho_w g)}
$$

$$
k_{r,e}=k_{pmax}\frac{RAI_e}{\sum_e RAI_e}
$$

然后两段阻力串联：

$$
k_{combined,e}
=
\frac{k_{s,e}k_{r,e}}{k_{s,e}+k_{r,e}}
$$

最后转换到 MPa 单位：

$$
g_{sr,e}^{CLM5}=k_{combined,e}\frac{10^6}{\rho_w g}
$$

代码对应：

```matlab
r_soil_m = sqrt(1.0 / (pi * L_root_density));
k_s = K_eff_m_s / r_soil_m;
k_r = rootSideConductanceFromArea(...);
k_combined = (k_s * k_r) / (k_s + k_r + 1e-30);
layer_endpoint_gsr = k_combined * m_water_per_MPa;
```

注意：$RAI_e$ 已经进入 $k_{r,e}$ 的分配，因此 `layer_endpoint_gsr` 后面不能再额外乘一次 $RAI_e$，否则会重复计入根量。

来源：CLM5/CTSM Technical Note, Plant Hydraulics, Section 2.11.1.2 Root Spacing 与 2.11.2.1.3 Soil-to-root；其中给出 $dx_{root,i}=(\pi L_i)^{-1/2}$、$k_{s,i}=k_{soil,i}/dx_{root,i}$，以及 soil-root interface resistance 与 root tissue resistance 串联的 $k_{3,i}=k_{r,i}k_{s,i}/(k_{r,i}+k_{s,i})$。CTSM 开发者论坛也说明，代码实现中 RAI 应修正 `krmax`，也就是 root conductance term，而不是再乘到最终 soil-to-root conductance 上。

#### 2.4.3 Radial Darcy + root-side series conductance

新 `radial_darcy_series` 分支保留 CLM5/CTSM 的“土壤侧 + 根组织侧串联”结构，但将 soil-side length-scale approximation 替换为显式的 cylindrical radial Darcy integration。

对一个半径为 $r_{root}$、长度为 $L$ 的圆柱根段，径向 Darcy 方程为：

$$
Q=-K A(r)\frac{dh}{dr}
$$

圆柱径向流面积为：

$$
A(r)=2\pi rL
$$

因此：

$$
Q=-K(2\pi rL)\frac{dh}{dr}
$$

整理并从根表面 $r_{root}$ 积分到外边界 $r_{outer}$：

$$
dh=-\frac{Q}{2\pi KL}\frac{dr}{r}
$$

$$
h_{soil}-h_{root}
=
\frac{Q}{2\pi KL}
\ln\left(\frac{r_{outer}}{r_{root}}\right)
$$

所以单根段的径向水头导度为：

$$
G_{soil}
=
\frac{2\pi KL}
{\ln\left(r_{outer}/r_{root}\right)}
$$

其中对数项来自 $\int (1/r)dr=\ln r$，不是经验拟合项。

对有限元端点 $e$，令根长密度为：

$$
B_e=\frac{\text{root length}}{\text{soil volume}}
$$

单位为 m root m$^{-3}$ soil，即 m$^{-2}$。若端点代表的土层支撑厚度为：

$$
z_{sup,e}=\frac{\Delta z}{2}
$$

则单位地表面积上的端点根长为：

$$
L_e=B_e z_{sup,e}
$$

代入圆柱径向 Darcy 导度，得到 soil-side 水头单位导度：

$$
k_{s,e}^{radial}
=
\frac{2\pi K_{eff,e}B_e z_{sup,e}}
{\ln\left(r_{outer,e}/r_{root}\right)}
$$

单位为 s$^{-1}$。

根面积指数与根长密度的几何关系为：

$$
RAI_e=2\pi r_{root}B_e z_{sup,e}
$$

因此上式等价于：

$$
k_{s,e}^{radial}
=
\frac{K_{eff,e}RAI_e}
{r_{root}\ln\left(r_{outer,e}/r_{root}\right)}
$$

根外边界半径由根长密度估计。假设根在土体中近似均匀分布，一单位根长分配到的平均截面积为 $1/B_e$，近似为圆形：

$$
\pi r_{outer,e}^2=\frac{1}{B_e}
$$

所以：

$$
r_{outer,e}=\frac{1}{\sqrt{\pi B_e}}
$$

并且：

$$
\Lambda_e
=
\ln\left(\frac{r_{outer,e}}{r_{root}}\right)
=
\ln\left(
\frac{1}{r_{root}\sqrt{\pi B_e}}
\right)
$$

最终 soil-side radial Darcy conductance 可写为：

$$
k_{s,e}^{radial}
=
\frac{2\pi K_{eff,e}B_e z_{sup,e}}{\Lambda_e}
$$

或等价地：

$$
k_{s,e}^{radial}
=
\frac{K_{eff,e}RAI_e}{r_{root}\Lambda_e}
$$

代码对应：

```matlab
B_m_per_m3 = max(root_density_endpoint(ML, ND), 1e-3);
r_outer_m = 1.0 / sqrt(pi * B_m_per_m3);
r_outer_m = max(r_outer_m, r_root_m * 1.01);
log_term = max(log(r_outer_m / r_root_m), 1e-6);
k_s = K_eff_m_s * root_area_endpoint(ML, ND) / ...
    (r_root_m * log_term);
```

这与以下形式等价：

```matlab
k_s = 2 * pi * K_eff_m_s * B_e * z_support_e / log_term;
```

前提是：

$$
RAI_e=2\pi r_{root}B_ez_{sup,e}
$$

##### Root amount source consistency

该分支必须保证根间距与根表面积来自同一套根几何。

若存在实测根长密度：

$$
B_e=B_{measured,e}
$$

则端点根面积由几何关系推出：

$$
RAI_e^{geom}=2\pi r_{root}B_e z_{sup,e}
$$

此时 root spacing 和 absorbing area 都来自 measured $B$。

若没有实测根长密度，则使用模型给定的端点根面积：

$$
RAI_e=RAI_{model,e}
$$

并反推等效根长密度：

$$
B_e=
\frac{RAI_e}
{2\pi r_{root}z_{sup,e}}
$$

这样避免用 measured $B$ 计算 $r_{outer}$，却用 unrelated modeled $RAI$ 计算吸水面积。

代码对应：

```matlab
if has_measured_B
    root_density_endpoint = [measured_B, measured_B];
    root_area_endpoint = 2 * pi * r_root_m ...
        * root_density_endpoint .* endpoint_support;
else
    root_area_endpoint = RAI;
    root_density_endpoint = root_area_endpoint ./ ...
        max(2 * pi * r_root_m * endpoint_support, 1e-30);
end
```

其中 `endpoint_support = 0.5 * DeltZ_m`。`root_density_endpoint = [measured_B, measured_B]` 表示每层两个局部端点使用同一层实测 $B$；根量守恒来自每个端点只代表半层支撑厚度，而不是把 $B$ 本身除以 2。

##### Root-side series coupling

root-side conductance 与 CLM5-style 分支一致：

$$
k_{r,e}
=
f_{root,e}\frac{RAI_e K_{root,max}}{L_{path,e}}
$$

其中 $RAI_e$ 在代码中为 source-consistent 的 `root_area_endpoint`。这对应 CTSM/CLM5 `root_conductance = fs * rai * krmax / (croot_average_length + z)`。若 `root_hydraulic_krmax=0` 且 `root_side_mode=none`，代码直接令 $k_{combined,e}=k_{s,e}$，即 root-side 不参与串联，$g_{sr}$ 主要表示 soil-side 与根几何/根土接触控制。若 `root_side_mode=stemmus_scope_rrr`，代码移植 STEMMUS-SCOPE 原 `rrr` 根径向阻力：

$$
r_{rrr,e}=4\times 10^{11}
\frac{\theta_{s,e}}{\theta_e}
\frac{1}{B_e z_{sup,e}}
$$

其中 $B_e$ 为端点根长密度，$z_{sup,e}=\Delta z/2$。两端点并联时，$B_e(0.5\Delta z)+B_e(0.5\Delta z)=B_e\Delta z$，因此总根长支持量与原 STEMMUS-SCOPE layer-scale `Rl*DeltZ` 一致。原 `ebal.m` 中水势以水头 m 表示、蒸腾以 m s$^{-1}$ 表示，因此该阻力在离散水头方程中作为秒量纲使用，代码取 $k_{r,e}=1/r_{rrr,e}$ 得到水头导度 [s$^{-1}$]，再与 $k_{s,e}$ 串联。整株植物实际导度仍在 `calSPAC` 中由 `hydraulic_gpmax * vulnerability(psil)` 控制。

CTSM/CLM5 同时允许 signed root-soil exchange：`qflx_rootsoi_col` 为正表示土壤到根，为负表示根到土壤的 hydraulic redistribution，并在土壤水方程中以 `rmx = qin - qout - qflx_rootsoi_col` 进入。因此负 `Srt` 本身不是守恒错误；守恒检查必须把它作为土壤方程内部源汇处理，而不是外部新增水源。然后：

$$
k_{combined,e}
=
\frac{k_{s,e}k_{r,e}}{k_{s,e}+k_{r,e}}
$$

并转换为 MPa 单位：

$$
g_{sr,e}^{radial}
=
k_{combined,e}\frac{10^6}{\rho_w g}
$$

代码对应：

```matlab
k_r = rootSideConductanceFromArea(plantHydro, root_area_endpoint, ...
    ML, ND, k_s, m_water_per_MPa);
k_combined = (k_s * k_r) / (k_s + k_r + 1e-30);
layer_endpoint_gsr = k_combined * m_water_per_MPa;
```

同样不能在最后再乘一次 $RAI_e$，因为根量已经进入 $k_s$，并在启用 `root_hydraulic_krmax` 时也进入 $k_r$。

来源：径向 Darcy 推导对应经典 cylindrical root radial flow / microscopic root water uptake 处理，可追溯到 Gardner (1960) 的根间距思想，后续 HYDRUS 和 de Jong van Lier 类模型也使用圆柱根与 rhizosphere radius 的径向流框架；CLM5/CTSM 提供 root spacing 与 soil/root series resistance 的陆面模型结构；当前 `radial_darcy_series` 是这两类思想在 COMPASS 端点导度装配中的组合实现。

---

## 3 蒸腾供水能力与二分法求解

### 3.1 根-茎水势与共享供水能力

SCOPE 能量平衡给出阳叶和阴叶蒸腾需求 $E_{sun}^{dem}$、$E_{sha}^{dem}$（m s⁻¹），总需求为 $E_{total}^{dem} = E_{sun}^{dem} + E_{sha}^{dem}$。对任意候选总蒸腾量 $E$，根据 Darcy 定律，根水势为

$$
\psi_r = \psi_s - \frac{E}{G_{sr}}
\tag{7}
$$

计入冠层高度 $H$（m）引起的重力水势降，茎端共享节点水势为

$$
\psi_{stem} = \psi_r - \frac{\rho_w g H}{10^6}
\tag{8}
$$

木质部导水能力 $k$（m s⁻¹ MPa⁻¹）随水势降低而下降，其脆弱性曲线取 Weibull 型：

$$
k(\psi) = \frac{k_{max}}{1 + \left(\psi / P_{50}\right)^a}
\tag{9}
$$

其中 $k_{max}$ 为最大整株导水度，$P_{50}$ 为导水损失 50% 时对应的水势（MPa，负值），$a$ 为曲线陡峭程度参数。

对脆弱性曲线在水势轴上积分，定义累计导水函数 $\mathcal{K}$：

$$
\mathcal{K}(\psi) = \int_{\psi_{min}}^{\psi} k(\psi') \, d\psi'
\tag{10}
$$

其中 $\psi_{min}$ 为水力查找表的最低叶水势下限（本研究取 −6 MPa）。在当前苹果参数下，该值表示重度失导但并非 P95 或死亡阈值；若有物种脆弱性曲线观测，可将其重新设为目标失导比例对应的水势。$\mathcal{K}$ 的量纲为 m s⁻¹，与蒸腾通量量纲相同，因此可以直接比较。

从当前茎水势 $\psi_{stem}$ 到 $\psi_{min}$ 之间还能支持的最大叶端水流，即**共享供水能力**：

$$
E_{supply}(E) = \mathcal{K}\!\left(\psi_{stem}(E)\right) - \mathcal{K}(\psi_{min})
\tag{11}
$$

式（11）说明供水能力本身依赖于候选蒸腾量 $E$：$E$ 越大，$\psi_{stem}$ 越低，剩余水势空间越小，供水能力越弱。这一非线性自反馈是整个水力闭合问题的核心。

### 3.2 二分法求解实际蒸腾

实际总蒸腾量 $E^*$ 须同时满足：

$$
E^* \leq E_{total}^{dem}, \quad E^* \leq E_{supply}(E^*)
\tag{12}
$$

即 $E^*$ 不超过需求，也不超过该通量下根-茎系统的供水能力。若需求量处于供水能力允许范围内（$E_{total}^{dem} \leq E_{supply}(E_{total}^{dem})$），则 $E^* = E_{total}^{dem}$；否则，$E^*$ 为方程 $E = E_{supply}(E)$ 的不动点。

由于 $E_{supply}(E)$ 关于 $E$ 单调递减，不动点唯一，采用**二分法**求解：在 $[0,\; E_{total}^{dem}]$ 区间内迭代，每步取中点 $\hat{E}$，若 $\hat{E} \leq E_{supply}(\hat{E})$ 则上移下界，否则下移上界，迭代至区间宽度小于收敛阈值（本研究取 10⁻¹⁰ m s⁻¹），共最多迭代 80 步，通常 30 步内收敛。

---

## 4 阳叶/阴叶分配与叶水势反求

### 4.1 按需求比例分配

二分法仅求解总实际蒸腾 $E^*$，将其按阳叶/阴叶原始需求比例分配：

$$
E_{sun}^* = E^* \cdot \frac{E_{sun}^{dem}}{E_{total}^{dem}}, \quad
E_{sha}^* = E^* \cdot \frac{E_{sha}^{dem}}{E_{total}^{dem}}
\tag{13}
$$

此分配策略保留了 SCOPE 能量平衡中因辐射、温度和气孔阻力差异造成的阳叶/阴叶需求比例，同时保证二者之和不超过共享供水上限。

### 4.2 叶水势反求

确定各叶端实际通量后，重新计算与之一致的根水势和茎水势（式 7、8），再利用累计导水函数 $\mathcal{K}$ 反求叶水势。以阳叶为例：

$$
\mathcal{K}(\psi_{l,sun}) = \mathcal{K}(\psi_{stem}) - E_{sun}^*
\tag{14}
$$

对式（14）求逆，得到阳叶叶水势 $\psi_{l,sun}$；阴叶叶水势 $\psi_{l,sha}$ 同理。由于 $\mathcal{K}$ 单调递增，反查唯一且数值稳定。实现中将 $\mathcal{K}$ 预先列表（水势步长 0.01 MPa），用线性插值完成正反两个方向的查表，避免每时步重复积分。

水势梯度的物理一致性：由式（7）、（8）和（14）可知，两个叶端水势均应低于共享茎水势，茎水势低于根水势，根水势低于根区加权土水势：

$$
\min(\psi_{l,sun}, \psi_{l,sha}) \leq
\max(\psi_{l,sun}, \psi_{l,sha}) \leq
\psi_{stem} \leq \psi_r \leq \psi_s
\tag{15}
$$

在当前白天情景中，阳叶蒸腾需求通常大于阴叶，因此常见关系为 $\psi_{l,sun} \leq \psi_{l,sha}$。若特殊条件下阴叶需求超过阳叶，二者相对顺序可以反转，但二者仍应共同低于共享茎水势。本模块在每个时步运行结束后均检查式（15）所表达的水势方向是否成立，以供调试使用。

---

## 5 水力胁迫因子与能量平衡耦合

### 5.1 水力胁迫因子

将叶水势代入脆弱性曲线（式 9），对最大导水度归一化，得到阳叶和阴叶各自的**水力胁迫因子**：

$$
w_{t,sun} = \frac{k(\psi_{l,sun})}{k_{max}}, \quad
w_{t,sha} = \frac{k(\psi_{l,sha})}{k_{max}}
\tag{16}
$$

$w_t \in [0, 1]$：值为 1 表示无水力胁迫，值越小说明木质部栓塞程度越高、气孔开度越受限。$w_t$ 作为应力乘子传入 SCOPE 的 `biochemical` 子程序，直接缩减最大羧化速率 $V_{cmax}$，从而同时影响光合速率和气孔阻力 $r_{cw}$。

### 5.2 与能量平衡的内迭代耦合

未耦合版本中，本时步的水力胁迫 $w_t$ 只能在**下一次**能量平衡外迭代中才能反映到 $r_{cw}$ 和潜热 $\lambda E$，导致气孔阻力与水力状态在同一温度迭代内不一致。

本研究在能量平衡外迭代的每个温度步内嵌入一个**水力-气孔内迭代环**，使 $w_t$、$r_{cw}$ 和 $\lambda E$ 在同一温度迭代内自洽闭合：

$$
\underbrace{\text{biochemical}(w_t^{(j)}) \to r_{cw}^{(j+1)}}_{\text{气孔更新}}
\xrightarrow{\lambda E^{(j+1)}}
\underbrace{\text{hydraulics} \to w_t^{(j+1)}}_{\text{水力更新}}
\tag{17}
$$

迭代终止条件为 $\max\!\left(|w_{t,sun}^{(j+1)} - w_{t,sun}^{(j)}|,\; |w_{t,sha}^{(j+1)} - w_{t,sha}^{(j)}|\right) < 10^{-3}$，或达到最大内迭代次数（6 次）。实践中多数时步在 2–3 次内收敛。

---

## 6 主要参数说明

| 参数 | 符号 | 默认值 | 单位 | 来源 |
|:--|:--:|:--:|:--:|:--|
| 最大整株导水度 | $k_{max}$ | $8\times10^{-8}$ | m s⁻¹ MPa⁻¹ | 文献 / 校准 |
| 50% 失导水势 | $P_{50}$ | −3.8 | MPa | 苹果茎部实测范围 / 校准 |
| 脆弱性曲线形状 | $a$ | 3 | — | 文献 / 校准 |
| 叶水势下限 | $\psi_{min}$ | −6.0 | MPa | 数值查表下限 / 校准 |
| CONTACT 起始阈值 | $\psi_{start}$ | −0.05 | MPa | Lei 等（2023）|
| CONTACT 终止阈值 | $\psi_{end}$ | −1.5 | MPa | Lei 等（2023）|

当前配置采用苹果参数 $P_{50} = -3.8$ MPa、$a=3$。在该参数下，$\psi_{min}=-6$ MPa 约对应 80% 导水损失，而不是 P95；它在本文中作为水力查找表的数值下限和重度干旱截断边界，不应解释为死亡阈值。若希望使用 P95 作为下限，则应按 $\psi_{P95}=P_{50}(19)^{1/a}$ 重新计算并校准。CONTACT 阈值参照 Lei 等（2023）向日葵实验场景，本文将其作为可校准参数，应用于苹果或其他物种时需结合根系形态、土壤质地和实测水势重新校准。

---

## 7 模块验证

以下验证在固定蒸腾需求（来自白天真实气象强迫下的单步 ebal 结算结果：$E_{sun}^{dem} = 6.0 \times 10^{-6}$ cm s⁻¹，$E_{sha}^{dem} = 1.5 \times 10^{-6}$ cm s⁻¹）的条件下进行，不更新土壤水分，仅扰动植物水力输入，逐项检验各物理过程的单调性与合理性。

### 7.1 KL_h 压力测试：共享供水截断

将根区土壤导水率矩阵整体乘以缩放因子 $s$（$s$ 从 1 递减至 $10^{-8}$），考察实际/需求蒸腾比 $E^*/E^{dem}$、叶水势和水力胁迫因子的响应。

**表 1　KL_h 缩放压力测试结果**

| $s$ | $E^*/E^{dem}$ | $w_{t,sun}$ | $w_{t,sha}$ | $\psi_{stem}$ (MPa) | $\psi_{l,sun}$ (MPa) | $\psi_{l,sha}$ (MPa) |
|:--:|:--:|:--:|:--:|:--:|:--:|:--:|
| $10^{0}$ | 1.000 | 0.763 | 0.939 | −0.769 | −1.625 | −0.966 |
| $10^{-1}$ | 1.000 | 0.763 | 0.939 | −0.769 | −1.625 | −0.966 |
| $10^{-2}$ | 1.000 | 0.763 | 0.939 | −0.769 | −1.625 | −0.966 |
| $10^{-3}$ | 1.000 | 0.763 | 0.939 | −0.769 | −1.625 | −0.966 |
| $10^{-4}$ | 1.000 | 0.762 | 0.938 | −0.771 | −1.627 | −0.968 |
| $10^{-5}$ | 1.000 | 0.755 | 0.935 | −0.789 | −1.649 | −0.986 |
| $10^{-6}$ | 1.000 | 0.674 | 0.896 | −0.965 | −1.884 | −1.169 |
| $10^{-7}$ | 0.798 | 0.138 | 0.428 | −2.330 | −4.421 | −2.643 |
| $10^{-8}$ | 0.179 | 0.076 | 0.131 | −4.271 | −5.517 | −4.509 |

**物理解读：**

- $s \geq 10^{-5}$：基准条件下蒸腾需求较小，$G_{sr}$ 仍能支撑全量需求，$E^*/E^{dem} = 1$。随 $s$ 减小，$G_{sr}$ 等比例缩小，相同通量下 $\psi_r$ 的降幅（式 7：$E/G_{sr}$）加速变大，叶水势随之变负，$w_t$ 逐步下降，但尚未触发供水截断。

- $s = 10^{-7}$：$G_{sr}$ 降至极低，供水能力开始不足以满足需求，二分法求得 $E^* = 0.798 \, E^{dem}$；$\psi_{l,sun}$ 降至 −4.4 MPa，$w_{t,sun}$ 骤降至 0.14，阳叶胁迫显著重于阴叶（阳叶需求占比更大，分配通量更大，茎端以下额外压降更大）。

- $s = 10^{-8}$：供水严重不足，实际蒸腾仅为需求的 18%，$\psi_{l,sun}$ 接近查找表下限 −6 MPa，$w_t$ 逼近 0，模型进入重度干旱保护状态。

上述响应证实：（1）共享供水二分法正常运作，阳叶和阴叶不能各自独立消耗超出根-茎总供水能力的水分；（2）阳叶较阴叶优先受到水力胁迫，符合其需求比例更高的物理预期；（3）水势梯度链 $\psi_{l,sun} \leq \psi_{l,sha} \leq \psi_{stem} \leq \psi_r \leq \psi_s$ 在所有缩放条件下均成立（见式 15）。

### 7.2 CONTACT 响应验证

保持 $K_{soil}$ 不变，将所有土壤节点水势均匀设定为目标值，考察 CONTACT 因子 $f_{contact}$ 和供水响应。

**表 2　CONTACT 机制响应验证结果**

| $\psi_{soil}$ (MPa) | $f_{contact}$ | $E^*/E^{dem}$ | $w_{t,sun}$ | $w_{t,sha}$ | $\psi_r$ (MPa) |
|:--:|:--:|:--:|:--:|:--:|:--:|
| −0.04 | 1.000 | 1.000 | 0.964 | 0.999 | −0.040 |
| −0.50 | 0.690 | 1.000 | 0.863 | 0.977 | −0.500 |
| −1.50 | 0.000 | 0（无供水）| — | — | — |
| −2.00 | 0.000 | 0（无供水）| — | — | — |

**物理解读：**

- $\psi_{soil} = -0.04$ MPa：土壤湿润，$f_{contact} = 1$，根-土导度无损失。根水势等于土水势（−0.04 MPa），水力链极短，$w_t$ 接近 1，蒸腾满足全量需求。

- $\psi_{soil} = -0.50$ MPa：处于两阈值之间（$-1.5 < \psi < -0.05$ MPa），$f_{contact} = 0.69$，有效导水率降低约 31%。$G_{sr}$ 等比例下降，同等需求下 $\psi_r$ 降幅增大，$w_{t,sun}$ 降至 0.86，但总蒸腾量仍满足需求，说明蒸腾需求在此土水势下仍在供水能力范围内。

- $\psi_{soil} \leq -1.50$ MPa：$f_{contact} = 0$，根-土界面导度完全丧失，$G_{sr} \to 0$，植物无法从土壤取水，模型报告"无根导度"异常并将实际蒸腾设为零。这一截断机制忠实地表达了极度干旱下根-土接触断裂的物理过程。

以上验证表明，CONTACT 分段线性机制对土水势变化的响应单调、物理合理，且在导度归零时能够稳定终止计算，不产生数值奇异。

---

## 8 干旱归因：大气干旱与土壤干旱的分离实验

### 8.1 科学问题

植物水分胁迫由两类物理机制驱动，现有研究多单独考虑其中一种，本节围绕三个核心问题展开：

**Q1：VPD 需求增加与土壤供水下降是否产生非加性水力响应？**

当大气干旱（VPD 升高）与土壤干旱（$\psi_{soil}$ 下降）同时发生时，植物响应是否可由两者线性叠加预测？若不能，非线性协同效应的量级有多大？用**协同比**量化：

$$\eta = \frac{\Delta\psi_{l,interact}}{\left|\Delta\psi_{l,atm}\right| + \left|\Delta\psi_{l,soil}\right|} \tag{14a}$$

$\eta \neq 0$ 表示存在非加性协同，$|\eta|$ 越大说明两类干旱的耦合越强，单独建模的误差越大。

**Q2：叶水势下降主要由大气需求端还是土壤供给端控制？**

在不同 VPD × $\psi_{soil}$ 组合下，大气贡献占比 $f_{atm}$ 和土壤贡献占比 $f_{soil}$ 分别为：

---

## 9 2026-04-20 说明：逐层诊断与默认共享路径

为配合非均匀 `LAI_layer`、冠层空气剖面和分层蒸腾需求，`ebal_clm45` 当前会输出逐层：

```text
Tp_sun_layer_demand / Tp_sha_layer_demand
wt_sun_layer / wt_sha_layer
psil_sun_layer / psil_sha_layer
```

但这些量必须按以下原则解释：

### 9.1 默认模式

默认模式保持本文前述主拓扑不变：

```text
soil -> root -> shared stem node -> sunlit leaf / shaded leaf
```

因此默认 `layer_hydraulic_mode = shared_path_scalar_repeated`，含义是：

1. 先按共享根-茎路径求解总供水能力和 `wt_sun / wt_sha`；
2. 再把该聚合解复制到层维度输出；
3. 逐层蒸腾需求仅用于诊断冠层内部蒸腾需求分布，而不改写主水力闭合。

### 9.2 实验模式

若显式设置：

```text
options.use_layer_hydraulic_stress = true
```

则会启用实验性的 `LAI_proportional_v1`：

1. 仍先求聚合总供水；
2. 再按层 `demand / LAI` 比值近似重建逐层叶端 `psi_leaf` 和 `wt`；
3. 用于判断“哪一层先受旱”，而不是重新定义整株主风险。

### 9.3 解释边界

实验模式下若出现：

```text
wt_top < wt_aggregate < wt_bottom
```

并不表示整株新的主风险已经改变，只表示局部叶层相对更受限或更不受限。  
由于脆弱曲线是非线性的，逐层 `wt` 不能反推出新的整株 `wt` 主控制量。

### 9.4 情景验证结论

新增情景矩阵脚本：

```text
tests/check_layered_scenario_transitions.m
```

验证了白天/夜间、湿润/干旱、湿转干、干转湿等情景。结论是：

1. 默认共享路径模式在上述情景下均可稳定运行；
2. 白天干旱显著降低 `supply_ratio` 与 `wt`，夜间需求显著减小；
3. 跃迁情景可跑通，最大能量残差约小于 `0.21 W m-2`；
4. 因此当前推荐的默认主模型配置仍应保持共享路径，而将逐层 `wt` 保留为实验诊断量。

$$f_{atm} = \frac{\left|\Delta\psi_{l,atm}\right|}{\left|\Delta\psi_{l,atm}\right| + \left|\Delta\psi_{l,soil}\right|}, \quad f_{soil} = 1 - f_{atm} \tag{14b}$$

这一比例直接影响调控策略：若 $f_{soil}$ 主导，灌溉是主要手段；若 $f_{atm}$ 显著，遮荫降温更有效。

**Q3：大气干旱和土壤干旱是否在冠层 WUE、波文比和冠层温度上留下可区分信号？**

大气干旱通过气孔导度 $g_s$ 同时影响蒸腾（$E$）和光合（$A$），冠层水分利用效率（Canopy WUE $= A_c/E_c$）的变化方向取决于 $A_c$ 与 $E_c$ 的相对变化；土壤干旱通过根-土导度和水力截断限制实际供水，进而抑制 $A_c$，同时使感热（$H$）升高、冠层温度（$T_c$）升高。两条路径在能量平衡分配上的特征差异可用于**观测诊断**干旱类型。

---

两种干旱机制在模型中天然耦合：$e_a$ 改变通过完整能量平衡（ebal）影响气孔导度和蒸腾需求；土壤水势通过根-茎水力通路（SPAC 二分法）影响实际蒸腾供给。两者共同决定叶水势 $\psi_l$、水力胁迫因子 $w_t$、光合速率 $A$ 和感热通量 $H$。

### 8.2 实验设计

以白天基准状态（1990-09-19 12:00，北京站，$T_a = 18.05$ °C，$e_a = 14.28$ hPa，$\psi_{soil} \approx -0.77$ MPa）为参照点，分三组扰动情景：

| 情景 | 控制变量 | 扰动变量 | 实现方式 |
|---|---|---|---|
| ATM | 土壤保持原状 | $e_a \in [14.28, 4]$ hPa（VPD 从 6.4 升至 16.7 hPa） | 修改 `meteo.ea`，完整重跑 ebal |
| SOIL | $e_a$ 固定为基准值 | $\psi_{soil} \in [-0.04, -2.0]$ MPa | 修改 `SoilVariables.hh`，完整重跑 ebal |
| COMPOUND | 两者均变 | $e_a \in \{14.28, 10, 6\}$ hPa × $\psi_{soil} \in \{-0.04, -0.5, -1.0, -1.5\}$ MPa | 修改 `meteo.ea` 和 `hh`，完整重跑 ebal |

ATM 情景直接扰动实际气象强迫 $e_a$，通过完整能量平衡重新求解叶温、气孔导度和蒸腾，而非仅缩放需求，保证物理一致性。

### 8.3 归因分解方法

#### 8.3.1 分解公式

对阳叶叶水势偏差做加法分解，参照点为观测基准状态（$e_a = 14.28$ hPa，原始土壤剖面 $\psi_s \approx -0.77$ MPa）下的 $\psi_{l,ref}$。因此，SOIL 情景中的 $\psi_{soil}=-0.04$ MPa 是一个人为湿土扰动情景，而不是归因分解的基准点：

$$\Delta\psi_{l,compound}(e_a,\,\psi_s) = \Delta\psi_{l,atm}(e_a) + \Delta\psi_{l,soil}(\psi_s) + \Delta\psi_{l,interact}(e_a,\,\psi_s) \tag{15}$$

各项定义如下：

$$\Delta\psi_{l,atm}(e_a) = \psi_l^{ATM}(e_a) - \psi_{l,ref} \tag{16}$$

$$\Delta\psi_{l,soil}(\psi_s) = \psi_l^{SOIL}(\psi_s) - \psi_{l,ref} \tag{17}$$

$$\Delta\psi_{l,interact}(e_a,\psi_s) \overset{\text{def}}{=} \psi_l^{COMPOUND}(e_a,\psi_s) - \psi_{l,ref} - \Delta\psi_{l,atm}(e_a) - \Delta\psi_{l,soil}(\psi_s) \tag{18}$$

整理后等价于：

$$\psi_l^{COMPOUND} = \underbrace{\psi_{l,ref} + \Delta\psi_{l,atm} + \Delta\psi_{l,soil}}_{\text{线性预测（假设两效应独立）}} + \Delta\psi_{l,interact} \tag{19}$$

#### 8.3.2 计算步骤

分解来自**四个相互独立的数值实验**，各步骤分开运行，互不依赖：

**第一步：基准实验**

使用真实气象强迫（$e_a = 14.28$ hPa）和原始土壤状态完整运行 ebal，记录：

$$\psi_{l,ref} = -1.6297 \text{ MPa}$$

**第二步：ATM 实验（只改 $e_a$，土壤保持原状）**

对每个 $e_a$ 值独立重跑完整 ebal，土壤状态不变：

$$\Delta\psi_{l,atm}(e_a) = \psi_l^{ATM}(e_a) - \psi_{l,ref}$$

例：$e_a = 10$ hPa 时，$\psi_l^{ATM} = -1.6636$ MPa，故 $\Delta\psi_{l,atm} = -0.034$ MPa。

**第三步：SOIL 实验（只改 $\psi_{soil}$，$e_a$ 保持基准）**

对每个 $\psi_{soil}$ 值完整重跑 ebal，$e_a$ 不变：

$$\Delta\psi_{l,soil}(\psi_s) = \psi_l^{SOIL}(\psi_s) - \psi_{l,ref}$$

例：$\psi_s = -0.5$ MPa 时，$\psi_l^{SOIL} = -1.3012$ MPa，故 $\Delta\psi_{l,soil} = +0.3285$ MPa（正值表示该土壤条件比真实土壤更湿，叶水势更高）。

**第四步：线性预测（不运行任何新实验）**

将第二、三步结果代入式 (19) 左侧，计算"若两效应完全独立"的预期叶水势：

$$\psi_{l,\text{linear}} = \psi_{l,ref} + \Delta\psi_{l,atm} + \Delta\psi_{l,soil} = -1.6297 + (-0.034) + 0.3285 = -1.3352 \text{ MPa}$$

此值**从未对应任何真实实验**，仅为假设独立时的理论预测。

**第五步：COMPOUND 实验（同时改 $e_a$ 和 $\psi_{soil}$）**

对 $e_a \times \psi_{soil}$ 网格中每个组合完整重跑 ebal，得到真实物理结果 $\psi_l^{COMPOUND}$。

**第六步：交互项（式 18 的直接计算）**

$$\Delta\psi_{l,interact} = \psi_l^{COMPOUND} - \psi_{l,\text{linear}}$$

例：$e_a = 10$ hPa，$\psi_s = -0.5$ MPa 时，若 $\psi_l^{COMPOUND} = -1.3915$ MPa，则：

$$\Delta\psi_{l,interact} = -1.3915 - (-1.3352) = -0.0563 \text{ MPa}$$

#### 8.3.3 关键说明：为何可以这样分解

交互项由**定义（式 18）**保证成立，是代数恒等式而非假设。其意义在于：若 $\Delta\psi_{interact} = 0$，两效应完全可加，模型行为线性；若 $\Delta\psi_{interact} \neq 0$，则反映了 SPAC 方程的非线性耦合：

$$\psi_{stem} = \psi_{soil} - \frac{E_{actual}}{g_{sr}} - \rho_w g H$$

当 VPD 升高（$E_{actual}$ 需求增大）与土壤变干（$g_{sr}$ 减小）同时发生时，$\psi_{stem}$ 的下降幅度超过两者线性叠加，进入脆弱性曲线加速衰减区，产生**协同加剧**的非线性交互。

交互项为负（$\Delta\psi_{l,interact} < 0$）表示复合效应比线性叠加更强（协同加剧胁迫）；为正表示略有缓解（通常发生在极干土条件下，脆弱性曲线已进入平坦段）。

### 8.4 数值结果与图示

![干旱归因分析六联图](../tests/figures/drought_attribution.png)

**图 8.1** 干旱归因分析（1990-09-19 12:00，北京站）。上排为大气干旱情景，下排左两图为土壤干旱情景，右下图为复合干旱热力图。

---

#### 子图 [1,1]（左上）：大气干旱 — 叶水势 vs VPD

X 轴为 VPD（hPa），Y 轴为叶水势（MPa）。

| $e_a$ (hPa) | VPD (hPa) | $T_p/T_{p,dem}$ | $\psi_{l,sun}$ (MPa) | $w_{t,sun}$ | $A_{sun}$ (μmol m⁻² s⁻¹) | 冠层 WUE | $\beta$ |
|---|---|---|---|---|---|---|---|
| 14.28（基准） | 6.4 | 1.000 | −1.625 | 0.928 | 45.4 | 10.11 | 0.868 |
| 12.00 | 8.7 | 1.000 | −1.658 | 0.923 | 44.9 | 9.54 | 0.785 |
| 10.00 | 10.7 | 1.000 | −1.674 | 0.921 | 44.7 | 9.28 | 0.741 |
| 8.00 | 12.7 | 1.000 | −1.677 | 0.921 | 44.6 | 9.24 | 0.722 |
| 6.00 | 14.7 | 1.000 | −1.691 | 0.919 | 44.4 | 9.12 | 0.696 |
| 4.00 | 16.7 | 1.000 | −1.732 | 0.914 | 43.8 | 8.62 | 0.632 |

**解读**：VPD 从 6.4 升至 16.7 hPa（增幅 160%），$\psi_{l,sun}$ 仅从 −1.625 降至 −1.732 MPa（变化 −0.11 MPa），$w_t$ 从 0.928 降至 0.914（变化仅 0.014）。大气干旱对叶水势的直接影响**极小**，供水始终未截断（$T_p/T_{p,dem} = 1$）。采用苹果实测 $p_{50} = -3.8$ MPa 后，$w_t$ 整体显著高于草本参数（旧值 0.76），反映苹果木质部在此叶水势范围内几乎未发生导水损失。

---

#### 子图 [1,2]（中上）：大气干旱 — 胁迫因子与供水比 vs VPD

双 Y 轴图：左轴为 $w_t$，右轴为 $T_{p,actual}/T_{p,demand}$。

$w_t$ 随 VPD 升高略有下降但幅度极小（0.928 → 0.914），供水比始终 = 1。冠层 WUE 随 VPD 升高从 10.1 下降至 8.6（$E_c$ 增速快于 $A_c$ 下降），波文比 $\beta$ 从 0.87 降至 0.63，表明高 VPD 下更多能量以潜热形式消耗。大气干旱的气孔响应主要通过 VPD 驱动蒸腾增加实现，而不产生液压截断。

---

#### 子图 [1,3]（右上）：归因分解 — 堆叠柱图

X 轴为 VPD 情景（三组），每组三根柱对应三个 $\psi_{soil}$（−0.04 / −0.50 / −1.00 MPa）。柱高由三部分堆叠：蓝色（$\Delta\psi_{atm}$）、橙色（$\Delta\psi_{soil}$）、绿色（$\Delta\psi_{interact}$）。

| $e_a$ (hPa) | $\psi_{soil}$ (MPa) | $\Delta\psi_{atm}$ | $\Delta\psi_{soil}$ | $\Delta\psi_{interact}$ |
|---|---|---|---|---|
| 14.28 | −0.04 | 0.000 | +0.827 | −0.083 |
| 14.28 | −0.50 | 0.000 | +0.329 | −0.049 |
| 14.28 | −1.00 | 0.000 | −0.304 | +0.058 |
| 10.00 | −0.04 | −0.034 | +0.827 | −0.099 |
| 6.00 | −1.00 | −0.063 | −0.304 | +0.056 |

**解读**：橙色柱（土壤贡献）主导总变化；蓝色柱（大气贡献）量级小；绿色交互项为负时表示复合效应强于线性叠加（协同加剧），为正时表示复合效应弱于线性叠加（拮抗或饱和缓冲）。该图的参考点为原始观测土壤剖面，因此 $\psi_{soil}=-0.04$ MPa 表示相对基准更湿的土壤扰动。

---

#### 子图 [2,1]（左下）：土壤干旱 — 叶水势响应

X 轴为 $\psi_{soil}$（MPa，从右往左变干），Y 轴为叶水势（MPa）。

| $\psi_{soil}$ (MPa) | $T_p/T_{p,dem}$ | $\psi_{l,sun}$ (MPa) | $w_{t,sun}$ | $A_{sun}$ (μmol m⁻² s⁻¹) | 冠层 WUE | $\beta$ |
|---|---|---|---|---|---|---|
| −0.04 | 1.000 | −0.890 | 0.987 | 46.7 | 10.19 | 0.825 |
| −0.10 | 1.000 | −0.950 | 0.985 | 46.6 | 10.18 | 0.827 |
| −0.30 | 1.000 | −1.150 | 0.973 | 46.4 | 10.17 | 0.836 |
| −0.50 | 1.000 | −1.352 | 0.957 | 46.1 | 10.14 | 0.847 |
| −0.70 | 1.000 | −1.554 | 0.936 | 45.6 | 10.12 | 0.861 |
| −1.00 | 1.000 | −1.862 | 0.895 | 44.7 | 10.06 | 0.889 |
| −1.50 | 0.000 | −6.000 | 0.203 | 17.3 | 8.82 | 2.475 |
| −2.00 | 0.000 | −6.000 | 0.203 | 17.3 | 8.82 | 2.521 |

$\psi_l$ 随土壤干燥呈**明显非线性加速下降**，在 −0.04 到 −1.0 MPa 区间内下降 0.97 MPa，供水始终满足需求（$T_p/T_{p,dem} = 1$），属于**水势链传导**型胁迫。冠层 WUE 在 −1.0 MPa 前几乎不变（≈10.1），体现苹果在中度胁迫下的水分利用稳定性。$\beta$ 由 0.83 缓升至 0.89，表明胁迫越强越多能量转为感热。

---

#### 子图 [2,2]（中下）：土壤干旱 — 胁迫因子与供水比

$w_t$ 从 0.987（湿土）单调下降至 0.895（−1.0 MPa），供水比在 −1.0 MPa 前保持 1，在 −1.5 MPa 时因根土接触（CONTACT）因子归零而骤降为 0，$w_t$ 跌至 0.203（对应 $\psi_l = -6$ MPa 时苹果残余导水率）。$\beta$ 在 −1.5 MPa 处从 0.889 暴增至 2.475，蒸腾骤停使能量平衡完全转为感热。

**两条线的分离意义**：供水比（橙色）保持 1 而 $w_t$ 下降，说明胁迫来自**叶水势下降驱动气孔部分关闭**，而非供水量不足。供水比骤降（−1.5 MPa 处）是 CONTACT 机制的液压失败信号，而非渐进的供水衰减。与草本参数（$p_{50}=-2.4$ MPa，断崖处 $w_t=0.06$，$\beta=4.6$）相比，苹果参数（$p_{50}=-3.8$ MPa）在极端干旱下仍保留更多残余导水能力，断崖效应更"软"。

---

#### 子图 [2,3]（右下）：复合干旱热力图

X 轴为 $\psi_{soil}$（MPa），Y 轴为 VPD（hPa），颜色为 $\psi_{l,sun}$（MPa）。红色等值线表示供水比 $T_p/T_{p,dem}=0.9$ 的截断边界，黑色虚线表示当前 CONTACT 参数 $\psi_{end}=-1.5$ MPa。

**三个关键观察**：

1. **颜色沿 X 方向（土壤）变化显著**：同一 VPD 下，$\psi_{soil}$ 从 −0.04 变至 −1.0 MPa，$\psi_{l,sun}$ 变化约 1.0 MPa；而同一 $\psi_{soil}$ 下，VPD 从 7 升至 14 hPa，颜色几乎不变（变化 < 0.1 MPa）。**土壤干旱是主控因子**。

2. **CONTACT 黑线近似竖直**：当前液压失败阈值（$\psi_{soil} = -1.5$ MPa，CONTACT 参数 $\psi_{end}$）不随 VPD 移动，表明完全根-土接触截断由土壤水势参数化决定，与大气需求无直接关系。红色供水比等值线用于标识实际蒸腾开始低于需求的区域。

3. **黑线右侧（$\psi_{soil} \leq -1.5$ MPa）**：当前 CONTACT 参数使根-土导度归零，蒸腾骤停，$\psi_{l,sun}$ 跌至下限 −6 MPa；但苹果 $p_{50}=-3.8$ MPa 使 $w_t$ 在此处保留约 0.20，仍有残余光合（$A \approx 17$ μmol）。

### 8.4b 扩展分析：光合、冠层 WUE、波文比与协同比

![扩展分析六联图](../tests/figures/drought_attribution_extended.png)

**图 8.2** 扩展分析（苹果参数 $p_{50}=-3.8$ MPa，同一时步）。上排为大气干旱情景，下排为土壤干旱情景；左列光合与冠层 WUE，中列波文比 $\beta$，右列协同比 $\eta_\psi$ 与贡献占比 $f_{atm}$。

---

![PPT版归因分量热力图](../tests/figures/drought_attribution_components_ppt.png)

**图 8.3** PPT 版归因分量热力图。四个子图分别显示 $\Delta\psi_{atm}$、$\Delta\psi_{soil}$、$\Delta\psi_{interact}$ 和总响应 $\Delta\psi_{compound}$，均相对原始观测基准土壤剖面计算。该图用于直观展示土壤供给端贡献显著大于大气需求端贡献，以及交互项量级较小。

![干旱机制流程图](../tests/figures/drought_mechanism_summary_ppt.png)

**图 8.4** 机制流程图。VPD 主要从需求端提高 $E_{dem}$，土壤干旱主要从供给端降低 $G_{sr}$ 和 CONTACT 因子，两者通过共享 SPAC 供水闭合共同决定叶水势和 $w_t$，再反馈到气孔阻力、潜热、感热和光合。

---

#### 大气干旱扩展（上排）

**左图（$A_{sun}$ 与冠层 WUE vs VPD）**：$A_{sun}$ 从 45.4 μmol 缓降至 43.8 μmol（−3.5%），VPD 增大通过提升蒸腾需求间接抑制光合，但幅度有限。冠层 WUE 从 10.1 降至 8.6（−15%），原因是冠层蒸腾 $E_c$ 随 VPD 增大的速率快于 $A_c$ 的下降，表明**大气干旱以降低冠层水分利用效率为主要代价**，而非大幅压低光合总量。

**中图（$\beta$ vs VPD）**：波文比 $\beta$ 从 0.87 降至 0.63，随 VPD 升高蒸腾驱动的潜热增加，能量平衡向潜热倾斜。大气干旱期间冠层持续散热，不会引发显著升温。

**右上图（$\eta_\psi$ 热图）**：各格数值范围 $|\eta_\psi| < 0.022$，正负混杂，中等 VPD（≈11 hPa）×中等土壤干旱（−0.5 ~ −1.0 MPa）区域出现轻微负协同（协同加剧）；极端土壤干旱（−1.5 MPa）时 $\eta$ 转正（拮抗），因 CONTACT 断裂后 $\psi_l$ 已锁在下限，VPD 附加效应消失。在当前案例和参数范围内，**协同效应微弱**，两类干旱对叶水势的作用近似线性叠加。

---

#### 土壤干旱扩展（下排）

**左图（$A_{sun}$ 与冠层 WUE vs $\psi_{soil}$）**：$A_{sun}$ 在 −1.0 MPa 前从 46.7 缓降至 44.7 μmol（−4%），冠层 WUE 几乎不变（约 10.1），体现苹果在中度土壤胁迫下的**高水分利用稳定性**。$\psi_{soil}=-1.5$ MPa 处出现断崖：$A_{sun}$ 跌至 17.3 μmol（残余 37%）。当冠层潜热接近零时，WUE 的比值会变得不稳定，因此极端截断点只用于说明能量分配和液压失败，不作为稳定的叶片尺度 WUE 结论。

**中图（$\beta$ vs $\psi_{soil}$）**：$\beta$ 在中度胁迫段从 0.83 缓升至 0.89（蒸腾轻微受限），在 −1.5 MPa 处骤升至 2.48（蒸腾骤停、净辐射转向感热）。与旧草本参数（$\beta_{max}=4.6$）相比，苹果的断崖更"软"（$\beta_{max}=2.5$），因 $p_{50}=-3.8$ MPa 在 $\psi_l=-6$ MPa 时仍保有 20% 残余导度、维持部分蒸发散热。

**右下图（$f_{atm}$ 热图）**：大气贡献占比在湿润土壤（$\psi_{soil}=-0.04$ MPa）、高 VPD（15 hPa）时最大，约 8%；随土壤变干 $f_{atm}$ 先升后降，在 $\psi_{soil}=-1.0$ MPa、VPD=15 hPa 时达到最大约 22%；极端干旱（$\psi_{soil}=-1.5$ MPa）时因蒸腾归零、$\psi_l$ 锁定，$f_{atm}$ 骤降至 1.5%。全域 $f_{atm} < 25\%$，**土壤干旱在所有条件下主导叶水势变化**。

---

### 8.5 综合结论

在苹果参数（$p_{50}=-3.8$ MPa，$\psi_{end}=-1.5$ MPa）与本站点气象强迫（1990-09-19 12:00，北京）下：

- **土壤干旱是叶水势和 $w_t$ 下降的主控因子**：$f_{atm}<25\%$，土壤贡献量约为大气干旱的 4 倍以上，在极端胁迫下差距超过 10 倍。
- **大气干旱主要代价是降低冠层 WUE**（10.1→8.6），而非大幅压低 $\psi_l$（变化仅 −0.11 MPa）或截断供水（$T_p/T_{p,dem}=1$）。
- **苹果木质部高耐旱性（$p_{50}=-3.8$ MPa）使模型响应更"软"**：相同叶水势下 $w_t$ 约为草本参数的 1.2 倍，断崖处 $\beta_{max}=2.5$（草本为 4.6），极端干旱下残余光合约 17 μmol（草本仅 6 μmol）。
- **当前参数下的液压失败阈值（$\psi_{soil}=-1.5$ MPa）由 CONTACT 机制决定**，不依赖 VPD，是土壤-根界面参数化给出的供水约束；该阈值应视为可校准参数，而非苹果根系的固定生理常数。
- **在当前苹果参数、该白天时刻和所选扰动范围内，协同效应微弱（$|\eta_\psi|<0.022$）**，两类干旱对叶水势的作用近似线性叠加。该结论不排除在更高蒸腾需求、更浅根系、更低 $g_{p,max}$ 或更脆弱 $P_{50}$ 参数下出现更强非线性协同。

### 8.6 判别指标小结

| 指标 | 大气干旱主导 | 土壤干旱主导 |
|---|---|---|
| $T_p/T_{p,dem}$ | ≈ 1（供水不截断） | < 1（供水截断） |
| $\psi_s \to \psi_r$ 梯度 | 小 | 大 |
| $\psi_r \to \psi_l$ 梯度 | 大（需求拉降） | 小 |
| $w_t$ 下降根本原因 | 叶水势随高 VPD 下降 | 叶水势随供水能力衰减 |
| CONTACT 效应 | 不显著 | 干土末期主导截断 |
| 热力图颜色变化方向 | Y 方向（VPD轴） | X 方向（ψ_soil轴） |
| 冠层 WUE = $A_c/E_c$ 变化方向 | 趋于下降（$E_c$ 增大先于 $A_c$ 下降） | 中度胁迫下可稳定或略升；极端截断时比值不稳定 |
| 冠层温度 $T_c$ 响应 | 轻微升高 | 显著升高（蒸发冷却减弱） |

### 8.7 可扩展分析

基于模型已输出的变量（`Actot`、`lEctot`、`Htot`、`Gtot`、`wt_sun/sha`），可在现有三情景框架下直接扩展以下分析，无需新增实验：

#### 8.7.1 光合速率归因分解

对净光合 $A$（`fluxes.Actot`，μmol m⁻² s⁻¹）做与 $\psi_l$ 完全相同的分解：

$$\Delta A_{compound} = \Delta A_{atm} + \Delta A_{soil} + \Delta A_{interact}$$

大气干旱通过 VPD 升高直接增大气孔导度需求，导致 $c_i/c_a$ 下降、$A$ 受限；土壤干旱通过 $w_t$ 降低气孔导度，间接抑制 $A$。两条路径的贡献比例 $f_{atm}^A$、$f_{soil}^A$ 与叶水势的比例 $f_{atm}^{\psi}$、$f_{soil}^{\psi}$ 一般不相同，比较二者可揭示**气孔调控的非对称性**。

#### 8.7.2 冠层水分利用效率（Canopy WUE）归因

$$\text{Canopy WUE} = \frac{A_c}{E_c} = \frac{\text{fluxes.Actot}}{\text{fluxes.lEctot} / \lambda}$$

- **大气干旱**：VPD 升高使 $E_c$ 增大（VPD 驱动蒸腾）而 $A_c$ 受限（$c_i$ 下降），冠层 WUE 趋于**下降**。
- **土壤干旱**：中度胁迫下气孔调节可使 $E_c$ 与 $A_c$ 同步下降，冠层 WUE 可保持稳定或略升；极端 CONTACT 截断时 $E_c$ 接近零，WUE 比值不再稳定，应结合 $A_c$、$E_c$ 和 $\beta$ 共同解释。
- **复合干旱**：两个方向可能相反，冠层 WUE 的净变化取决于哪类干旱主导，是区分两类干旱的重要观测指标。

#### 8.7.3 能量平衡分配与冠层温度

波文比 $\beta = H/\text{LE}$（`Htot/lEtot`）随干旱强度的变化：

| 情景 | LE 变化 | H 变化 | $\beta$ 变化 | $T_c$ 变化 |
|---|---|---|---|---|
| 大气干旱（VPD↑） | ↑（需求增大） | ↓ | ↓ | 略升（辐射不变） |
| 土壤干旱（$\psi_s$↓） | ↓（供水截断） | ↑ | ↑↑ | 显著升高 |
| 复合干旱 | ↓ | ↑ | ↑ | 最高 |

土壤干旱导致的 $T_c$ 升高会进一步加剧 VPD（$e_s(T_c)$ 升高），形成**正反馈**，这在单时步模型中被 ebal 内部迭代捕捉，在多时步模拟中则会随时间累积。

#### 8.7.4 协同比 $\eta$ 的参数敏感性

改变物种参数（$p_{50}$、$a_a$、$g_{p,max}$、$\psi_{l,min}$）重复上述分析，可回答：

- **耐旱物种**（$p_{50}$ 更负）：$\psi_l$ 下降相同幅度时 $w_t$ 下降更慢，$f_{atm}$ 相对更小。
- **浅根物种**（根系深度浅）：$g_{sr}$ 对干旱更敏感，$f_{soil}$ 更大，协同比 $\eta$ 更大。
- **高 $g_{p,max}$ 物种**：蒸腾需求更易得到满足，大气干旱效应相对土壤干旱更大。

以上分析均可基于现有 `check_drought_attribution.m` 脚本直接扩展，只需在记录 `h.psil_sun` 的同时增加记录 `fluxes_i.Actot`、`fluxes_i.lEctot`、`fluxes_i.Htot` 即可。

---

## 9 与前期方案的对比和创新点

相对于经典 STEMMUS-SCOPE 采用经验土壤水分胁迫因子（SMCsf）的方案，本模块的主要改进包括：

1. **显式水力通路**：用真实的土壤导水率 $K_{soil}$ 和根面积分布 $RAI$ 构建土壤-根导度，代替基于土壤含水量的查找表胁迫因子，使水力响应具有明确的物理单位和传递链。

2. **根区加权土水势**：以导度为权重的加权平均（式 6）代替单层或均匀土水势假设，捕捉多层土壤异质供水的垂向差异。

3. **阳叶/阴叶共享供水**：将 SCOPE 双叶模型与根-茎共享水力通路显式耦合，阳叶和阴叶不再独立竞争供水，保证水量守恒与物理一致性。

4. **脆弱性曲线驱动的叶水势反求**：通过累计导水函数的正反查表（式 10、14），在每个时步直接输出叶水势，而非仅输出一个经验胁迫指数，为后续栓塞回复或 ABA 化学信号机制的引入提供接口。

5. **根-土接触修正（CONTACT）**：引入 Lei 等（2023）的分段线性接触因子，在体积土壤水力参数不变的前提下，描述干旱过程中根际导水率的额外衰减，对干旱末期蒸腾截断有实质性影响。

6. **水力-气孔内迭代耦合**：在能量平衡温度迭代内嵌入水力内迭代环（式 17），使同一时步的气孔阻力、潜热通量和叶水势达到内部一致，消除传统序贯耦合方案中的一步滞后误差。

---

## 10 动静结合扩展：准稳态叶水势闭合与动态植物储水库

### 10.1 扩展动机

当前植物水力模块采用准稳态 SPAC 闭合：在每个时步内，根据土壤水势、土壤-根导度和冠层蒸腾需求，求解共享根-茎供水能力、阳叶/阴叶实际蒸腾和叶水势。这一方案能够保证同一时步内的供需闭合和水势梯度一致，但植物体自身储水没有作为显式状态变量保存，因此对昼夜补水、午后恢复和连续干旱滞后响应的表达仍然有限。

为同时模拟叶水势（pressure state）和植物组织含水量（storage state），可在现有准稳态叶水势闭合外侧引入 ClimaLand 式根-茎-叶动态储水库。该扩展不替代当前的叶水势反求方法，而是将其作为供需闭合核心，并额外追踪植物组织储水量及其换算得到的体积含水量。

### 10.2 状态变量与通量定义

定义单位地表面积上的植物储水量：

$$
W_r(t) = W_{root}, \qquad
W_s(t) = W_{stem}, \qquad
W_l(t) = W_{leaf}
\tag{20}
$$

其中 $W_r$、$W_s$、$W_l$ 分别表示根、茎/木质部和叶片储水量，单位为 $m^3\,water\,m^{-2}\,ground$，等价于等效水柱深度（m）。

定义根-茎-叶水流通量：

$$
q_{sr} = q_{soil \to root}
\tag{21a}
$$

$$
q_{rs} = q_{root \to stem}
\tag{21b}
$$

$$
q_{sl} = q_{stem \to leaf}
\tag{21c}
$$

$$
T = T_{actual}
\tag{21d}
$$

其中 $T$ 为水力限制后的实际蒸腾通量。所有通量均按单位地表面积计，单位为 $m\,s^{-1}$。

### 10.3 动态储水方程

根、茎和叶三个储水库的连续时间方程为：

$$
\frac{dW_r}{dt} = q_{sr} - q_{rs}
\tag{22}
$$

$$
\frac{dW_s}{dt} = q_{rs} - q_{sl}
\tag{23}
$$

$$
\frac{dW_l}{dt} = q_{sl} - T
\tag{24}
$$

对应到模型离散时间步 $\Delta t$：

$$
W_r^{t+\Delta t}
= W_r^t + \Delta t\left(q_{sr} - q_{rs}\right)
\tag{25}
$$

$$
W_s^{t+\Delta t}
= W_s^t + \Delta t\left(q_{rs} - q_{sl}\right)
\tag{26}
$$

$$
W_l^{t+\Delta t}
= W_l^t + \Delta t\left(q_{sl} - T\right)
\tag{27}
$$

因此，当 $q_{rs} > q_{sl}$ 时，$dW_s/dt > 0$，表示茎/木质部储水补给；当 $q_{rs} < q_{sl}$ 时，$dW_s/dt < 0$，表示茎/木质部放水缓冲蒸腾需求。该过程表示储水补给与释放，不等同于木质部栓塞恢复。

### 10.4 与现有叶水势反求方法的耦合

现有模块首先求解准稳态叶水势：

$$
\psi_l^{steady}
=
\mathcal{F}
\left(
\psi_{soil}, K_{soil}, RAI, E_{dem}, g_{p,max}, P_{50}, a
\right)
\tag{28}
$$

其中 $\mathcal{F}$ 表示当前文档第 2--5 节给出的根区加权土水势、共享供水二分法、累计导水函数反查和阳叶/阴叶分配过程。

茎到叶的通量仍可由当前木质部脆弱性曲线积分定义：

$$
q_{sl}
=
\int_{\psi_l}^{\psi_{stem}} k(\psi)\,d\psi
\tag{29}
$$

其中：

$$
k(\psi)
=
\frac{g_{p,max}}
{1+\left(\frac{\psi}{P_{50}}\right)^a}
\tag{30}
$$

因此，新增动态储水库后，木质部在模型中具有两个角色：

1. 作为导水通路：由 $k(\psi)$ 和积分项控制茎-叶水流与叶水势反求；
2. 作为储水库：由 $W_s$ 记录茎/木质部储水的亏损与补给。

这可概括为：

$$
\text{原模型：木质部} = \text{有阻力的导水通路}
\tag{31a}
$$

$$
\text{扩展模型：木质部} = \text{有阻力的导水通路} + \text{动态储水库}
\tag{31b}
$$

### 10.5 植物体积含水量 VWC 的计算

借鉴 ClimaLand 的 compartment water content 表达，单位地表面积储水量可换算为组织体积含水量：

$$
\theta_{leaf}
=
\frac{W_l}{H_l \cdot LAI}
\tag{32}
$$

$$
\theta_{stem}
=
\frac{W_s}{H_s \cdot SAI}
\tag{33}
$$

$$
\theta_{root}
=
\frac{W_r}{H_r \cdot RAI}
\tag{34}
$$

其中 $\theta_{leaf}$、$\theta_{stem}$ 和 $\theta_{root}$ 分别为叶片、茎/木质部和根的体积含水量（VWC，$m^3\,m^{-3}$）；$H_l$、$H_s$、$H_r$ 为叶、茎、根水库的等效厚度或特征长度；$LAI$、$SAI$、$RAI$ 分别为叶面积指数、茎面积指数和根面积指数。

若实测数据为相对含水量（RWC）而非体积含水量，可进一步归一化：

$$
RWC_{leaf}
=
\frac{\theta_{leaf} - \theta_{leaf,min}}
{\theta_{leaf,max} - \theta_{leaf,min}}
\tag{35}
$$

或在缺少 $\theta_{min}$ 与 $\theta_{max}$ 的情况下，采用储水容量归一化的诊断形式：

$$
RWC_{leaf}
=
\frac{W_l}{W_{l,max}}
\tag{36}
$$

茎和根的 RWC 可按同样方式定义。

### 10.6 诊断型动态储水与反馈型动态储水

第一阶段建议采用诊断型动态储水：

$$
\psi_l = \psi_l^{steady}
\tag{37}
$$

即叶水势、水力胁迫和气孔反馈仍由现有准稳态闭合决定，动态储水库只根据 $q_{sr}$、$q_{rs}$、$q_{sl}$ 和 $T$ 更新，并输出 $W$、VWC 和 RWC。该方案不改变当前论文主线和已验证的叶水势结果，适合优先用于叶片含水量实测对比。

第二阶段可采用反馈型动态储水。首先由储水量换算得到储水控制的水势：

$$
\psi_{leaf}^{store}
=
\frac{1}{a_l}
\left(
\frac{\theta_{leaf}}{\nu_l} - 1
\right)
\tag{38}
$$

其中 $a_l$ 为叶片 retention 曲线斜率，$\nu_l$ 为叶片最大或饱和体积含水量。然后将准稳态水势与储水水势耦合：

$$
\psi_{leaf}^{coupled}
=
(1-\lambda)\psi_{leaf}^{steady}
+ \lambda \psi_{leaf}^{store}
\tag{39}
$$

其中 $\lambda$ 为动态储水反馈权重。也可用响应时间尺度定义：

$$
\lambda
=
1-\exp\left(-\frac{\Delta t}{\tau_{plant}}\right)
\tag{40}
$$

反馈后的水势进入原有脆弱性曲线：

$$
w_t
=
\frac{1}
{1+\left(\frac{\psi_{leaf}^{coupled}}{P_{50}}\right)^a}
\tag{41}
$$

随后继续反馈至气孔导度、光合速率、潜热通量和能量平衡。

### 10.7 与栓塞恢复/回水曲线的边界

动态储水库可以模拟茎/木质部储水的亏损与补给：

$$
\frac{dW_s}{dt} > 0
\quad \Rightarrow \quad
\text{stem water storage refilling}
\tag{42}
$$

但这不等同于木质部栓塞恢复。当前脆弱性曲线仍为单值可逆函数：

$$
k = k(\psi)
\tag{43}
$$

同一水势 $\psi$ 对应同一导水率 $k$，不包含历史损伤状态。若未来需要表示栓塞滞回，应额外引入木质部损伤变量：

$$
D_x(t) \in [0,1]
\tag{44}
$$

并令有效导水率为：

$$
k_{eff}
=
g_{p,max}
\min
\left[
1-D_x(t),
\frac{1}{1+\left(\frac{\psi}{P_{50}}\right)^a}
\right]
\tag{45}
$$

由于当前缺少茎水势、PLC 或导水率恢复实测，本文主模型暂不显式模拟 $D_x(t)$，而将 $W_s$ 解释为茎/木质部储水代理变量。

### 10.8 待做清单

1. 在 `plantHydro` 中增加动态储水状态变量：
   `W_root`, `W_stem`, `W_leaf`, `VWC_root`, `VWC_stem`, `VWC_leaf`, `RWC_root`, `RWC_stem`, `RWC_leaf`。

2. 在配置文件中增加储水参数：
   `H_root`, `H_stem`, `H_leaf`, `SAI`, `W_root_max`, `W_stem_max`, `W_leaf_max`, `theta_leaf_min`, `theta_leaf_max`, `theta_stem_min`, `theta_stem_max`。

3. 在现有 `solveEbalHydraulics` 或其 wrapper 后增加诊断型储水更新函数：
   `updatePlantWaterStorageDiagnostic(...)`。

4. 首版只输出动态储水和 VWC/RWC，不反馈 `psi_leaf`，用于和实测叶片含水量、木质部/茎含水量进行独立验证。

5. 若诊断型储水能重现实测 WC 日变化，再增加反馈型开关：
   `plant_storage_feedback_enabled = 0/1`。

6. 反馈型版本中使用式（39）--（41）将储水水势耦合进 `wt`，并重新评估叶水势、蒸腾和能量平衡是否改善。

7. 论文表述中区分三类过程：
   叶水势为 pressure state，VWC/RWC 为 storage state，ET/LE 为 flux state；当前模型模拟储水补给，不声称模拟栓塞恢复滞回。

---

## 11 可选扩展：基于木质部 PLC 连续暴露的水力失败风险

### 11.1 扩展动机

动态储水库描述的是根、茎/木质部和叶片储水量的亏损与补给；木质部脆弱性曲线描述的是水势降低时导水能力的瞬时损失。若要进一步评估持续干旱下的水力失败风险，可借鉴 ORCHIDEE-CAN-NHA 的思路，将木质部导水损失百分比（percentage loss of conductance, PLC）与连续暴露时间相结合。

该模块建议作为诊断型风险指数，而非主模型死亡过程。原因是当前果园应用缺少树体死亡率、茎水势和木质部导水率恢复观测，直接预测死亡率会引入难以约束的经验参数。

### 11.2 由木质部脆弱性曲线计算 PLC

当前模型已有木质部脆弱性曲线：

$$
k(\psi)
=
\frac{g_{p,max}}
{1+\left(\frac{\psi}{P_{50}}\right)^a}
\tag{46}
$$

归一化导水能力为：

$$
\frac{k(\psi)}{g_{p,max}}
=
\frac{1}
{1+\left(\frac{\psi}{P_{50}}\right)^a}
\tag{47}
$$

因此，木质部导水损失百分比可定义为：

$$
PLC(\psi)
=
100
\left[
1
-
\frac{k(\psi)}{g_{p,max}}
\right]
\tag{48}
$$

代入式（47）：

$$
PLC(\psi)
=
100
\left[
1
-
\frac{1}
{1+\left(\frac{\psi}{P_{50}}\right)^a}
\right]
\tag{49}
$$

在当前模型中，水力失败风险建议优先使用茎/木质部水势：

$$
\psi_{stem}
=
\psi_r - \rho_w g H
\tag{50}
$$

实际代码中等价于：

```matlab
pstem = psir - gravity_MPa;
```

因此，茎/木质部 PLC 为：

$$
PLC_{stem}(t)
=
100
\left[
1
-
\frac{1}
{1+\left(\frac{\psi_{stem}(t)}{P_{50,stem}}\right)^{a_{stem}}}
\right]
\tag{51}
$$

若缺少器官特异参数，第一版可令：

$$
P_{50,stem}=P_{50}, \qquad a_{stem}=a
\tag{52}
$$

但论文中应说明这是使用整株/木质部等效脆弱性参数的近似。

### 11.3 日尺度 PLC 指标

ORCHIDEE-CAN-NHA 使用日尺度连续暴露来触发死亡风险。对小时或半小时模型输出，可先聚合为日尺度 PLC。推荐使用每日最大 PLC：

$$
PLC_{stem,d}^{max}
=
\max_{t \in d} PLC_{stem}(t)
\tag{53}
$$

其中 $d$ 表示日期。每日最大 PLC 保留了一天中最危险的水力状态，适合用于水力失败风险诊断。若希望风险更保守，也可使用日均值：

$$
PLC_{stem,d}^{mean}
=
\frac{1}{N_d}
\sum_{t \in d} PLC_{stem}(t)
\tag{54}
$$

### 11.4 连续高 PLC 暴露时间

定义水力失败阈值：

$$
PLC_{crit}
\tag{55}
$$

ORCHIDEE-CAN-NHA 示例中取：

$$
PLC_{crit}=50\%
\tag{56}
$$

定义连续高 PLC 暴露天数：

$$
E_d
\tag{57}
$$

若当日 PLC 超过阈值：

$$
PLC_{stem,d}^{max} > PLC_{crit}
\tag{58}
$$

则暴露天数增加：

$$
E_d = E_{d-1}+1
\tag{59}
$$

若当日 PLC 不超过阈值，则并不立即清零，而是进入湿润中断计数。

### 11.5 湿润中断与暴露清零

定义湿润中断天数：

$$
B_d
\tag{60}
$$

当：

$$
PLC_{stem,d}^{max} \leq PLC_{crit}
\tag{61}
$$

则：

$$
B_d = B_{d-1}+1
\tag{62}
$$

ORCHIDEE-CAN-NHA 设置了一个最小湿润中断阈值：

$$
B_{crit}=5\,days
\tag{63}
$$

只有当低 PLC 状态连续维持足够久：

$$
B_d \geq B_{crit}
\tag{64}
$$

才认为树体从该轮高 PLC 暴露中恢复，并将暴露时间清零：

$$
E_d = 0
\tag{65}
$$

若湿润中断不足 $B_{crit}$ 天，随后又重新超过 $PLC_{crit}$，则继续累积原有暴露时间，而不是从零开始。

完整递推可写为：

$$
E_d =
\begin{cases}
E_{d-1}+1, & PLC_{stem,d}^{max} > PLC_{crit} \\[6pt]
0, & PLC_{stem,d}^{max} \leq PLC_{crit}\ \text{and}\ B_d \geq B_{crit} \\[6pt]
E_{d-1}, & PLC_{stem,d}^{max} \leq PLC_{crit}\ \text{and}\ B_d < B_{crit}
\end{cases}
\tag{66}
$$

湿润中断计数为：

$$
B_d =
\begin{cases}
0, & PLC_{stem,d}^{max} > PLC_{crit} \\[6pt]
B_{d-1}+1, & PLC_{stem,d}^{max} \leq PLC_{crit}
\end{cases}
\tag{67}
$$

### 11.6 水力失败风险指数

定义高风险暴露阈值：

$$
E_{crit}
\tag{68}
$$

ORCHIDEE-CAN-NHA 示例中取：

$$
E_{crit}=15\,days
\tag{69}
$$

可定义无量纲水力失败风险指数：

$$
R_{hyd,d}
=
\frac{E_d}{E_{crit}}
\tag{70}
$$

当：

$$
R_{hyd,d}<1
\tag{71}
$$

表示尚未达到持续高 PLC 暴露阈值；当：

$$
R_{hyd,d}\geq 1
\tag{72}
$$

表示进入高水力失败风险状态。

对于果园应用，建议优先输出 $R_{hyd}$，而不是直接输出死亡率。

### 11.7 ORCHIDEE-CAN-NHA 的死亡率计算

ORCHIDEE-CAN-NHA 在达到连续暴露阈值后，引入每日死亡比例：

$$
m_{frac}
\tag{73}
$$

文献示例中：

$$
m_{frac}=0.006\,day^{-1}
\tag{74}
$$

即每天死亡当前存活树体的 0.6%。当 $E_d \geq E_{crit}$ 时：

$$
M_d = m_{frac}
\tag{75}
$$

否则：

$$
M_d = 0
\tag{76}
$$

若以单棵树的存活概率 $S_d$ 表示，则：

$$
S_d = S_{d-1}(1-M_d)
\tag{77}
$$

累计死亡概率为：

$$
D_d = 1-S_d
\tag{78}
$$

如果连续 $n$ 天处于触发死亡状态，并且每天死亡比例固定为 $m_{frac}$，则累计死亡比例为：

$$
D(n)
=
1-(1-m_{frac})^n
\tag{79}
$$

代入 $m_{frac}=0.006$：

$$
D(n)
=
1-0.994^n
\tag{80}
$$

例如，在已经超过 15 天高 PLC 暴露阈值之后：

$$
D(1)=1-0.994^1=0.006
\tag{81}
$$

即第 1 个触发日死亡 0.6%。若之后又连续触发 15 天：

$$
D(15)=1-0.994^{15}\approx 0.086
\tag{82}
$$

即累计约 8.6% 死亡。若连续触发 30 天：

$$
D(30)=1-0.994^{30}\approx 0.165
\tag{83}
$$

即累计约 16.5% 死亡。

需要注意的是，ORCHIDEE-CAN-NHA 是 cohort 模型，0.6% 表示某一树径 cohort 内树体数量的每日死亡比例，而不是“单棵树死亡 0.6 棵”。对单棵树可解释为每日死亡概率，对果园群体可解释为该类树体的死亡比例。

### 11.8 对果园模型的建议用法

由于果园通常有灌溉、修剪、管理和人为补救措施，且当前缺少死亡观测，本文建议将该模块命名为：

$$
\text{Hydraulic Failure Risk Index}
\tag{84}
$$

而非直接命名为 mortality model。

推荐输出变量：

```text
PLC_stem_hourly
PLC_stem_daily_max
PLC_exposure_days
PLC_wetting_break_days
hydraulic_failure_risk
optional_mortality_fraction
```

在论文中可表述为：

1. $W_s$ / `VWC_stem` 描述木质部储水亏损与补给；
2. $PLC_{stem}$ 描述木质部导水能力损失；
3. $R_{hyd}$ 描述持续水力失败暴露风险；
4. 当前版本不显式模拟栓塞修复滞回，也不将死亡率作为主验证目标。

### 11.9 与当前代码的接入位置

当前植物水力主路径为：

```text
Init_Plant
  -> vegetation.root.calculateRootDistribution
       -> plantHydro.RAI
       -> plantHydro.root_length_density_m_m3   (optional measured-B mode)
       -> plantHydro.gsr_formula
  -> ebal_clm45
       -> vegetation.hydraulics.calSPAC
            -> vegetation.hydraulics.Plantrootconductivity
                 -> rootGeometryEndpointFromPlantHydro
                 -> katul / clm5 / radial_darcy_series gsr
                 -> psis, gsr, sumgsr
            -> shared SPAC supply bisection
            -> psir, psil_sun, psil_sha, wt_sun, wt_sha
```

其中 `rootGeometryEndpointFromPlantHydro` 负责保持根几何来源一致：

```text
measured-B mode:
    B_endpoint   = measured B
    RAI_endpoint = 2*pi*r_root*B_endpoint*z_support

RAI fallback mode:
    RAI_endpoint = model RAI
    B_endpoint   = RAI_endpoint/(2*pi*r_root*z_support)
```

随后 `Plantrootconductivity` 根据 `gsr_formula` 选择：

```text
katul                -> empirical bulk conductance
clm5                 -> CLM5/CTSM root spacing + soil/root series conductance
radial_darcy_series  -> radial Darcy soil side + root-side series conductance
```

当前 `ebal_clm45.m` 已在水力求解后得到：

```matlab
hydraulic.psir
hydraulic.psis
hydraulic.psil
hydraulic.wt_sun
hydraulic.wt_sha
```

并在根水势基础上计算过茎端水势：

```matlab
pstem = hydraulic.psir - gravity_MPa;
```

因此，水力失败风险模块可作为水力诊断后处理接入：

```text
1. 每个时步计算 pstem
2. 由 pstem 和脆弱性曲线计算 PLC_stem
3. 输出 hourly PLC_stem
4. 每日聚合 PLC_stem_daily_max
5. 更新 E_d 和 B_d
6. 输出 hydraulic_failure_risk
```

建议第一版只做诊断输出，不反馈 `wt`、`T_actual` 或碳分配。

---

## 参考文献

CTSM/CLM5 Technical Note. Section 2.11 Plant Hydraulics. Root spacing, soil-to-root conductance, and soil/root series resistance. https://escomp.github.io/CTSM/release-clm5.0/tech_note/Plant_Hydraulics/CLM50_Tech_Note_Plant_Hydraulics.html

CTSM DiscussCESM Forum. Calculation of soil to root conductance. Discussion noting that RAI modifies `krmax` / root conductance term rather than the final soil-to-root conductance. https://bb.cgd.ucar.edu/cesm/threads/calculation-of-soil-to-root-conductance.7679/

Gardner, W. R. (1960). Dynamic aspects of water availability to plants. *Soil Science*, 89, 63-73.

Lai, C. T., & Katul, G. (2000). The dynamic role of root-water uptake in coupling potential to actual transpiration. *Advances in Water Resources*, 23, 427-439. https://doi.org/10.1016/S0309-1708(99)00023-8

de Jong van Lier, Q., van Dam, J. C., Metselaar, K., de Jong, R., & Duijnisveld, W. H. M. (2008). Macroscopic root water uptake distribution using a matric flux potential approach. *Vadose Zone Journal*, 7, 1065-1078. https://doi.org/10.2136/vzj2007.0083

Lei, T., Pang, Z., Wang, X., Li, L., Fu, J., Kan, G., ... & Meng, X. (2023). Relating soil-root hydraulic resistance variation to stomatal regulation in soil-plant water transport modeling. *Journal of Hydrology*, 617, 128879. https://doi.org/10.1016/j.jhydrol.2022.128879
