function biochem_out = biochemical(biochem_in, sfactor, Ci_input)
% biochemical - 叶片尺度 C3/C4 光合–气孔–荧光耦合计算（CLM5 + Medlyn USO）
%               Leaf-scale C3/C4 photosynthesis–stomata–fluorescence (CLM5 + Medlyn USO)
%
% 函数定义 / Function Definition:
%   biochem_out = biochemical(biochem_in, sfactor, Ci_input)
%
% 计算流程 / Algorithm:
%   [1]  解包环境输入 / Unpack env inputs: T, Q, Cs, eb, p[Pa], Type
%   [2]  解包生理参数 / Unpack physio params:
%          Vcmax25, Rdparam, effcon, g_m, Knparams, Medlyn g0/g1, tempcor, po0, Kpep25(C4)
%   [3]  季节日长修正 & 氮胁迫因子 / fDYL & FN:
%          fDYL = (当日日长/夏至日长)²  / (current daylength / summer-solstice daylength)²
%          FN   = 全株氮限制因子（可选）/ whole-plant N stress (optional)
%   [4]  物理常数 & 单位换算 / Constants & unit conversions:
%          p[Pa→bar], O₂=209 mmol/mol (内置常数 / built-in constant)
%   [5]  CLM5 热适应 / Thermal acclimation (Kattge & Knorr 2007):
%          Tgrowth[11-35°C] → R_JV, ΔS_Vcmax, ΔS_Jmax
%   [6]  25°C 参考值 / Reference values at 25°C:
%          Vcmax25×fDYL×FN → Jmax25 / TPU25 / Rd25
%   [7]  温度响应函数 / Temperature response functions:
%          Arrhenius + Peaked 闭包 / Arrhenius + Peaked closures
%   [8]  动力学常数 / Kinetic constants @ Tleaf:
%          Kc, Ko, Γ* (Bernacchi 2001)
%   [9]  Vcmax/Jmax/TPU/Rd 温度修正 / Temperature corrections (Ha + Hd + ΔS)
%   [10] 水力胁迫 / Hydraulic stress:
%          Vcmax/Jmax/TPU × sfactor；Rd 不缩放 / Rd not scaled
%   [11] PSII & 电子传输 / PSII & electron transport:
%          po0=Kp/(Kf+Kd+Kp)，Jpot=0.5·po0·Q，Jact 二次共限制(θ=0.7)，Je=Jpot
%   [12] 湿度 & VPD / Humidity & VPD (Flatau formula)
%   [13] Michaelis-Menten 常数 / MM constants: Km=Kc(1+O/Ko)（C3 only）
%   [14] An–gs–Ci 耦合迭代 / Coupled iteration (≤50步/steps, 阻尼0.5/damp, 阈值0.01Pa):
%          Farquhar An(Ci) ↔ Medlyn gs0(An,VPD) [+ optional g_heat] ↔ Fick Ci=Cs-1.6·An/gs
%          两步平滑共限制 / Two-step smooth co-limitation (0.98/0.95, Collatz 1991/CLM5)
%   [15] 最终 gs & rcw / Final gs & stomatal resistance rcw [s/m]
%          gs = min(gs0 + g_heat, gs_cap); g_heat≡0 when heat_gs_enable=0 (production default)
%   [16] 荧光 / Fluorescence:
%          Ja=Ag/cpe, ps=po0·Ja/Je → van der Tol (2014) → eta/fs/Kn/qE/SIF
%   [17] 输出打包 / Pack outputs
%
% 输入 / Inputs:
%   biochem_in  - 输入参数结构体 / Input parameter struct
%     必需字段 / Required fields:
%       .T        - 叶温 / Leaf temperature [°C or K, scalar/nl×1/13×36×nl]
%       .Q        - 吸收 PAR / Absorbed PAR [μmol photons m⁻² s⁻¹]（已是 APAR / already APAR）
%       .Cs       - 叶面 CO₂ / Leaf surface CO₂ [ppm]
%       .eb       - 边界层水汽压 / Boundary layer vapour pressure [hPa]
%       .p        - 大气压 / Atmospheric pressure [Pa]
%       .Type     - 光合途径 / Photosynthetic pathway ['C3' | 'C4']
%       .Vcmo     - Vcmax25 层基底（含 fVh 梯度）/ Layer Vcmax25 (incl. fVh) [μmol m⁻² s⁻¹]
%       .Rdparam  - Rd/Vcmax25 比例 / Dark respiration fraction [-]
%     可选字段（有默认值）/ Optional fields (with defaults):
%       .Tgrowth    - 10天均温，热适应 / 10-day mean T for acclimation [°C]; default=Tc
%       .g0_medlyn  - Medlyn 最小导度 / min stomatal conductance [mol m⁻² s⁻¹]; default=1e-4
%       .g1_medlyn  - Medlyn 斜率 / Medlyn slope [kPa^0.5]; default=4.86
%       .heat_gs_*  - optional additive heat-induced gs after Medlyn (default enable=0);
%                     see docs/plant/stomata_heat_induced_gs_design.md
%       .Knparams   - NPQ 参数 / NPQ params [Kno, α, β]; default=[2.48, 2.83, 0.114]
%       .tempcor    - 温度响应开关 / T-response switch [0|1]; default=1
%       .latitude/.decl - 纬度/赤纬，用于 fDYL / for fDYL seasonal correction [°/rad]
%       .Nefc/.FN   - 全株氮胁迫开关/因子 / N stress switch/factor
%       .g_m        - 叶肉导度 / Mesophyll conductance [mol m⁻² s⁻¹]; default=Inf（关闭/off）
%   sfactor     - 水力胁迫因子 / Xylem hydraulic stress factor [0-1, scalar or nl×1]
%                 来自脆弱性曲线 / from KS(ψ_leaf)/gpmax
%   Ci_input    - 固定胞间 CO₂（诊断用）/ Fixed Ci for diagnostic mode [ppm]; []= iterate
%
% 输出 / Outputs:
%   biochem_out - 输出结构体 / Output struct
%     .A/.Ag    - 净/总光合 / Net/gross assimilation [μmol CO₂ m⁻² s⁻¹]
%     .Ci       - 胞间 CO₂ / Intercellular CO₂ [ppm]
%     .gs/.rcw  - 气孔导度/阻力 / Stomatal conductance/resistance [mol m⁻² s⁻¹ / s m⁻¹]
%     .gs0_medlyn / .g_heat - pure Medlyn gs0 and optional heat term [mol m⁻² s⁻¹]
%     .Vcmax    - 温度修正 Vcmax / T-corrected Vcmax [μmol m⁻² s⁻¹]
%     .Rd       - 暗呼吸 / Dark respiration [μmol m⁻² s⁻¹]
%     .fs/.eta  - 稳态荧光产额/相对荧光量 / Steady-state fluorescence yield/relative ratio
%     .Kn/.qE   - NPQ 速率常数/能量猝灭 / NPQ rate constant/energy quenching
%     .SIF      - 太阳诱导荧光 / Sun-induced fluorescence [μmol m⁻² s⁻¹]
%
% 使用示例 / Usage Example:
%   biochem_out = vegetation.biochemical(biochem_in, sfactor_layer);
%
% 引用 / References:
%   Farquhar et al. (1980) Planta — C3 三方共限制 / C3 co-limitation
%   Bernacchi et al. (2001) Plant Cell Environ — Kc/Ko/Γ* Arrhenius 参数
%   Collatz et al. (1991/1992) — 平滑共限制曲率 / smooth co-limitation curvature
%   Kattge & Knorr (2007) Plant Cell Environ — CLM5 热适应 / thermal acclimation
%   Medlyn et al. (2011) Global Change Biol — 最优气孔 / optimal stomatal conductance
%   van der Tol et al. (2014) Biogeosciences — 荧光模型 / fluorescence model
%
% 修改者 / Modified by:   朱寄子星 (Zhu Jizixing)
% 修改日期 / Modified date: 2026-04-22
% Editor: Jizixing Zhu; Date: 2026-04-22.
% Change: 2026-07-20 Jizixing Zhu — optional additive g_heat after Medlyn (heat_gs_enable default 0).

% --- [0] 默认参数处理 / Default argument handling ---
if nargin < 2 || isempty(sfactor)
    sfactor = 1;
end
if nargin < 3
    Ci_input = [];
end
if any(isnan(sfactor(:)))
    sfactor(isnan(sfactor)) = 1;
end

% 输出变量预分配（防止 compute_assimilation 嵌套调用前未初始化）
A = [];  Ag = [];  Vc = [];  Vs = [];  Ve = [];
CO2_per_electron = [];
fcount = 0;
warnings = [];

% =========================================================================
% --- [1] 环境输入解包 / Unpack environmental inputs ---
% =========================================================================
if isfield(biochem_in, 'Cs')
    Cs = max(0, real(biochem_in.Cs));   % 叶面 CO₂ [ppm]；强制非负实数
else
    Cs = NaN(size(biochem_in.T));
end
Q  = abs(real(biochem_in.Q));           % 吸收 PAR [μmol m⁻² s⁻¹]；abs 处理遮荫叶近零值
T  = biochem_in.T + 273.15 * (biochem_in.T < 200);  % 叶温 [K]；自动识别 °C 输入
Tc = T - 273.15;                        % 叶温 [°C]，用于饱和水汽压和 Kd
eb = biochem_in.eb;                     % 边界层水汽压 [hPa]
p  = biochem_in.p;                      % 大气压 [Pa]
Type = biochem_in.Type;                 % 'C3' 或 'C4'

% =========================================================================
% --- [2] 生理参数解包 / Unpack physiological parameters ---
% =========================================================================
% Vcmax25 基底值（已含 fVh 叶层垂直梯度，由 ebal 传入）
% [M2-B.3 LUNA 2026-05-20] 若 use_luna=1, 用 LUNA 每日优化输出覆盖 PFT 表常数
use_luna_vc = isfield(biochem_in, 'use_luna') && biochem_in.use_luna == 1 ...
        && isfield(biochem_in, 'Vcmax25_LUNA') && ~isempty(biochem_in.Vcmax25_LUNA);
if use_luna_vc
    Vcmax25_base = biochem_in.Vcmax25_LUNA;     % LUNA 输出 [μmol/m²/s]
else
    Vcmax25_base = biochem_in.Vcmo;             % legacy PFT 表常数
end
RdPerVcmax25 = biochem_in.Rdparam;          % [-]，由 AMainloop 从 photo_Rd_C3/C4 写入

% effcon：每固定一个 CO₂ 消耗的电子数的倒数 [mol CO₂ mol⁻¹ e⁻]
% C3: 1/4（Farquhar 1980 / von Caemmerer 2000 / CLM5，与 Jmax25 参数体系一致）；C4: 1/6
if isfield(biochem_in, 'effcon')
    effcon = biochem_in.effcon;
elseif strcmpi(Type, 'C3')
    effcon = 1 / 4;   % Farquhar 1980 / von Caemmerer 2000 / CLM5 标准：4e⁻ per CO₂
else
    effcon = 1 / 6;
end

% 叶肉导度 g_m [mol m⁻² s⁻¹ bar⁻¹]；默认 Inf（忽略叶肉阻力）
if isfield(biochem_in, 'g_m')
    g_m = biochem_in.g_m * 1e6;    % mol→μmol 单位对齐光合速率
else
    g_m = Inf;
end

fPAR = 1;   % SCOPE 传入 Q 已是 APAR，无需再乘叶片吸收比例

% tempcor：温度响应开关 [0|1]；=0 时所有 Arrhenius 项退化为常数（调试用）
if isfield(biochem_in, 'tempcor')
    tempcor = biochem_in.tempcor;
else
    tempcor = 1;
end

% NPQ 荧光参数 [Kno, alpha, beta]（van der Tol 2014）
% 默认 [2.48, 2.83, 0.114]（一般植物）；干旱模式 [5.01, 1.93, 10]
if isfield(biochem_in, 'Knparams')
    Knparams = biochem_in.Knparams;
elseif isfield(biochem_in, 'Kn0')
    Knparams = [biochem_in.Kn0, biochem_in.Kn_alpha, biochem_in.Kn_beta];
elseif isfield(biochem_in, 'Fluorescence_model') && biochem_in.Fluorescence_model == 0
    Knparams = [5.01, 1.93, 10];    % 干旱参数组
else
    Knparams = [2.48, 2.83, 0.114]; % 标准参数组（棉花数据集标定）
end

% po0：最大暗光化学效率 Fv/Fm [-]，若外部已测定可直接传入
if isfield(biochem_in, 'po0')
    po0 = biochem_in.po0;
else
    po0 = [];   % 空值将在第5节由 Kp/Kf/Kd 计算
end

% C4 PEPcase 活性（25°C） [μmol m⁻² s⁻¹]；默认 Vcmax25/56（Collatz 1992）
if isfield(biochem_in, 'Kpep25')
    Kpep25 = biochem_in.Kpep25;
elseif isfield(biochem_in, 'kp')
    Kpep25 = biochem_in.kp;
else
    Kpep25 = (Vcmax25_base / 56) * 1e6;
end


% Tgrowth：10天滑动均温，用于 CLM5 热适应参数计算 [°C]；缺省=当前叶温（无适应）
if isfield(biochem_in, 'Tgrowth') && ~isempty(biochem_in.Tgrowth)
    Tgrowth = biochem_in.Tgrowth;
else
    Tgrowth = Tc;
end

% Medlyn 气孔参数
if isfield(biochem_in, 'gs_model_type')
    gs_model_type = biochem_in.gs_model_type;   % 1=Medlyn（当前唯一实现）
else
    gs_model_type = 1;
end
if isfield(biochem_in, 'g0_medlyn')
    g0_medlyn = biochem_in.g0_medlyn;   % 最小气孔导度 [mol m⁻² s⁻¹]，默认 1e-4
else
    g0_medlyn = 1e-4;
end
if isfield(biochem_in, 'g1_medlyn')
    g1_medlyn = biochem_in.g1_medlyn;   % Medlyn 斜率 [kPa^0.5]，苹果树≈4.86
else
    g1_medlyn = 4.86;
end

% Optional heat-induced stomatal opening (additive after Medlyn; default OFF).
% Docs: docs/plant/stomata_heat_induced_gs_design.md; offline: tests/plant/offline_heat_gs_demo.m
% [SOURCE UNVERIFIED thresholds] T_on/T_full/wt/D/A gates are provisional until site calibration.
heat_gs = struct();
heat_gs.enable   = local_heat_gs_flag(biochem_in, 'heat_gs_enable', 0);
heat_gs.T_on     = local_heat_gs_num(biochem_in, 'heat_gs_T_on_C', 32.0);
heat_gs.T_full   = local_heat_gs_num(biochem_in, 'heat_gs_T_full_C', 38.0);
heat_gs.g_max    = local_heat_gs_num(biochem_in, 'heat_gs_max', 0.08);
heat_gs.wt_off   = local_heat_gs_num(biochem_in, 'heat_gs_wt_off', 0.4);
heat_gs.wt_on    = local_heat_gs_num(biochem_in, 'heat_gs_wt_on', 0.8);
heat_gs.D_lo     = local_heat_gs_num(biochem_in, 'heat_gs_D_lo_kPa', 1.5);
heat_gs.D_hi     = local_heat_gs_num(biochem_in, 'heat_gs_D_hi_kPa', 3.5);
heat_gs.A_lo     = local_heat_gs_num(biochem_in, 'heat_gs_A_lo', 1.0);
heat_gs.A_hi     = local_heat_gs_num(biochem_in, 'heat_gs_A_hi', 5.0);
heat_gs.gs_cap   = local_heat_gs_num(biochem_in, 'heat_gs_cap', 0.6);

% --- 日长因子 fDYL 与全株氮限制因子 FN ---
% fDYL (Kattge & Knorr 2007 / CLM5 Tech Note §9.2)：
%   fDYL = (当前日长 / 夏至最长日长)²，修正季节性 Vcmax25
%   latitude [°]；decl 来自 surface.radiation.calculatesunpos_ultimate_azi_atan22，单位弧度
if isfield(biochem_in, 'latitude') && isfield(biochem_in, 'decl') && ...
        ~isempty(biochem_in.latitude) && ~isempty(biochem_in.decl)
    lat_rad     = biochem_in.latitude * (pi / 180);   % 纬度 [rad]
    decl_rad    = biochem_in.decl;                    % 赤纬 [rad]，已是弧度
    arg_DYLmax  = -tan(lat_rad) * tan(23.45 * pi / 180);  % 夏至（最长日）
    arg_DYLmax  = max(-1, min(1, arg_DYLmax));
    arg_DYL     = -tan(lat_rad) * tan(decl_rad);          % 当日
    arg_DYL     = max(-1, min(1, arg_DYL));
    fDYL        = (acos(arg_DYL) / acos(arg_DYLmax))^2;
    fDYL        = max(min(fDYL, 1), 0.01);   % 钳位：夏至=1，极端冬季≥0.01
else
    fDYL = 1;   % 无地理信息时不修正
end

% FN：全株氮限制因子 [-] (legacy, photo_Nefc=1 且有时序 FN 数据时启用)
% 注：与 fVh（冠层内N垂直梯度）不重叠——fVh是空间结构分配，FN是整株水平应力
% [M2-B.3 LUNA 后此路径退为 legacy]
%   - lLUNA=1 (默认): N→Vcmax 反馈走 LUNA 机理路径, Current_FN = 1 不影响
%   - lLUNA=0 + photo_Nefc=1 + 外部 FN 时序: 走 legacy 经验缩放
%   - lLUNA=0 + photo_Nefc=0: 无 N 反馈 (静态 PFT 表 Vcmax25)
if isfield(biochem_in, 'Nefc') && biochem_in.Nefc == 1 && ...
        isfield(biochem_in, 'FN') && ~isempty(biochem_in.FN)
    Current_FN = biochem_in.FN;   % 时序氮限制因子，由外部模型传入
else
    Current_FN = 1;   % photo_Nefc=0：默认无氮胁迫
end

% =========================================================================
% --- [3] 物理常数与单位换算 / Physical constants and unit conversions ---
% =========================================================================
Tref = 298.15;      % 参考温度 [K] = 25°C
Rgas = 8.314;       % 理想气体常数 [J mol⁻¹ K⁻¹]
Kf   = 0.05;        % 荧光速率常数 [-]（van der Tol 2014）
Kd   = max(0.8738, 0.0301 .* Tc + 0.0773);  % 热耗散速率常数 [-]，随温度线性增加
Kp   = 4.0;         % 光化学速率常数 [-]
rhoa = 1.2047;      % 干空气密度 [kg m⁻³]（标准状态）
Mair = 28.96;       % 干空气摩尔质量 [g mol⁻¹]

p_bar    = p .* 1e-5;           % 大气压 [Pa→bar]
ppm2bar  = 1e-6 .* p_bar;      % CO₂ ppm→bar 换算因子 [bar ppm⁻¹]
phys = equations.physicalConstants();
if isfield(biochem_in, 'O2a_Pa') && ~isempty(biochem_in.O2a_Pa)
    O2a_Pa = biochem_in.O2a_Pa;
elseif isfield(biochem_in, 'Oa') && ~isempty(biochem_in.Oa)
    O2a_Pa = biochem_in.Oa .* 1e-3 .* p;
elseif isfield(biochem_in, 'O') && ~isempty(biochem_in.O)
    O2a_Pa = biochem_in.O .* 1e-3 .* p;
else
    O2a_Pa = phys.Oa_mmolmol_default .* 1e-3 .* p;
end
if strcmpi(Type, 'C3')
    O_bar = O2a_Pa .* 1e-5;             % O2 partial pressure [bar]
else
    O_bar = zeros(size(T));             % C4 PEP carboxylation uses equivalent O2=0
end

% =========================================================================
% --- [4] CLM5 热适应参数（Kattge & Knorr 2007）/ Thermal acclimation ---
% =========================================================================
% 生长温度钳位到 [11, 35]°C（CLM5 Tech Note Eq.9.6-9.7 有效范围）
Tgrowth_clamped = max(11, min(35, Tgrowth));    % [°C]

% Jmax/Vcmax 比值，随生长温度线性变化（热适应）
R_JV = 2.59 - 0.035 .* Tgrowth_clamped;        % [-]

% Peaked 函数高温抑制的熵参数 ΔS，随生长温度动态调整（热驯化核心）
DeltaS_Vcmax = 668.39 - 1.07 .* Tgrowth_clamped;  % [J mol⁻¹ K⁻¹]
DeltaS_Jmax  = 659.70 - 0.75 .* Tgrowth_clamped;  % [J mol⁻¹ K⁻¹]
DeltaS_TPU   = 490;     % [J mol⁻¹ K⁻¹]（固定，不随 Tgrowth 变）
DeltaS_Rd    = 490;     % [J mol⁻¹ K⁻¹]
Topt_Vcmax = peaked_arrhenius_topt_c(72000, 200000, DeltaS_Vcmax, Rgas); % [degC] diagnostic optimum from current Vcmax parameters
Topt_Jmax  = peaked_arrhenius_topt_c(50000, 200000, DeltaS_Jmax,  Rgas); % [degC] diagnostic optimum from current Jmax parameters

% =========================================================================
% --- [5] 25°C 参考参数 / Reference parameters at 25°C ---
% =========================================================================
% [LUNA fDYL/FN fix 2026-06-14] fDYL (seasonal day-length) and Current_FN (whole-plant N stress)
%   apply ONLY to the NON-LUNA (PFT/manual Vcmo) base. The LUNA path is SKIPPED: computeLuna already
%   ingests day-length (lunaIn.hourpd) and leaf N (lunaIn.lnc), so Vcmax25_LUNA/Jmax25_LUNA are the
%   FINAL 25C capacity -- re-multiplying by fDYL/FN would DOUBLE-COUNT day-length and N. Matches the
%   LUNA README flow (LUNA -> f_T x hydraulic, no fDYL/FN) and CLM5 (dayl_factor only in the non-LUNA
%   vcmax25, PhotosynthesisMod.F90:1378/1385). Vcmax25_base already has fVh (within-canopy N gradient).
if ~use_luna_vc
    Vcmax25_base = Vcmax25_base .* fDYL .* Current_FN;   % [μmol m⁻² s⁻¹] non-LUNA seasonal + N modifiers
end
% [Phase D 2026-05-26 修正] Jmax25 / TPU25 双路径:
%
% 路径 A (use_luna=1, CLM5 strict, LUNA optimization):
%   Source: CLM5 PhotosynthesisMod.F90:1515-1517 (Lawrence 2019 LUNA branch)
%     Vcmax25 = vcmx25_z (LUNA), 已通过 leafbio.Vcmo 覆盖
%     Jmax25  = jmx25_z (LUNA 独立优化, 不走 R_JV 比例)
%     TPU25   = 0.167 × Vcmax25 (CLM5 TPU 始终基于 Vcmax 不是 Jmax!)
%
% 路径 B (legacy, manual/PFT-table Vcmax, Kattge 2007 thermal acclimation):
%   Source: CLM5 PhotosynthesisMod.F90:1404-1405 (non-LUNA default branch)
%     Vcmax25 = leafbio.Vcmo (PFT 表或 config Photo_Vcmax_manual)
%     Jmax25  = R_JV(Tgrowth) × Vcmax25, R_JV=2.59-0.035×Tgrowth (Kattge 2007)
%     TPU25   = 0.167 × Vcmax25 (同样基于 Vcmax)
%
% 两路径区别只在 Jmax25 来源:
%   - 路径 A: Jmax25 独立 LUNA 优化
%   - 路径 B: Jmax25 跟 Vcmax 走 Kattge 经验比例
% TPU25 两路径都是 0.167 × Vcmax25 (CLM5 verbatim).
if isfield(biochem_in, 'use_luna') && biochem_in.use_luna == 1 ...
        && isfield(biochem_in, 'Jmax25_LUNA') && ~isempty(biochem_in.Jmax25_LUNA)
    % 路径 A: LUNA strict (Jmax 独立优化); [LUNA fDYL/FN fix 2026-06-14] NO fDYL/FN -- LUNA output is
    %   the final 25C Jmax (computeLuna already used hourpd + lnc); see the Vcmax25_base block above.
    Jmax25 = biochem_in.Jmax25_LUNA;
else
    % 路径 B: Kattge 2007 thermal acclimation (CLM5 default non-LUNA)
    Jmax25 = R_JV .* Vcmax25_base;     % [μmol m⁻² s⁻¹]
end
% [Phase D 修正] TPU25 = 0.167 × Vcmax25 (NOT × Jmax25)
% Source: CLM5 PhotosynthesisMod.F90:1405 + 1517 verbatim
%   Both LUNA and non-LUNA branches use 0.167 × vcmax25
% 旧 bug: TPU25 = 0.167 × Jmax25 让 TPU25 偏高 ~2 倍 (=0.167×R_JV×Vcmax≈0.33×Vcmax)
TPU25  = 0.167 .* Vcmax25_base;        % [μmol m⁻² s⁻¹]
% --- Rd25 base rate: Vcmax-proportional (default) OR CLM5 Atkin 2015 thermal-acclimated ---
% Default: Rd25 = RdPerVcmax25 × Vcmax25 (CLM5 leaf_mr_vcm method, PhotosynthesisMod.F90:1455;
%   ratio from config_vegetation.txt photo_Rd_C3/C4 via AMainloop→leafbio.Rdparam→biochem_in.Rdparam).
% [Atkin2015 leaf Rd 2026-06-14] Optional growth-temperature-acclimated base rate (CLM5
%   leafresp_mtd_atkin2015, PhotosynthesisMod.F90:1444-1448):
%     lmr25 = lmr_intercept_atkin + 0.2061·lnc − 0.0402·(Tgrowth[°C])
%   lnc = leaf N per leaf area = 1/(SLA·leafcn) [gN m⁻²]; the −0.0402·Tgrowth term is the growth-T
%   acclimation (10-day running mean), so warm-grown canopies down-regulate the basal Rd (homeostasis).
%   Atkin et al. 2015 New Phytol 206:614 (global leaf respiration). CLM5-FruitTree citrus params
%   (clm5_params.c171117__FruitTree.nc, PFT 35): lmr_intercept_atkin=1.756, slatop=0.0280,
%   leafcn=19.7 → lnc=1.81 gN m⁻². Same downstream T-response (f_Arr 46390 / f_Peak 150650) + Kok
%   light inhibition apply to BOTH branches (below), exactly as CLM5. Toggled by biochem_in.leafresp_atkin.
if isfield(biochem_in,'leafresp_atkin') && biochem_in.leafresp_atkin == 1 && ...
        isfield(biochem_in,'lmr_intercept_atkin') && ~isempty(biochem_in.lmr_intercept_atkin) && ...
        isfield(biochem_in,'lnc') && ~isempty(biochem_in.lnc) && biochem_in.lnc > 0
    Rd25 = max(biochem_in.lmr_intercept_atkin + 0.2061 .* biochem_in.lnc ...
               - 0.0402 .* Tgrowth, 0);    % [μmol m⁻² s⁻¹] CLM5 Atkin 2015 (Tgrowth in °C, 10-day mean)
else
    Rd25 = RdPerVcmax25 .* Vcmax25_base;   % [μmol m⁻² s⁻¹]；C3 默认 0.015，C4 默认 0.025
end

% =========================================================================
% --- [6] 温度响应函数 / Temperature response functions ---
% =========================================================================
% Arrhenius 激活能上升段（Bernacchi 2001）
%   f_Arr(Ha, T) = exp(Ha/R/Tref * (1 - Tref/T))，T=[K]，Ha=[J mol⁻¹]
f_Arr  = @(Ha, Tk) exp((Ha ./ (Tref .* Rgas)) .* (1 - Tref ./ Tk));

% Peaked 高温抑制项（CLM5 Eq.9.5）
%   f_Peak(Hd, ΔS, T)：分子=25°C归一化，分母=当前T的去激活因子
%   Hd=[J mol⁻¹]，ΔS=[J mol⁻¹ K⁻¹]
f_Peak = @(Hd, Ds, Tk) ((1 + exp((Tref .* Ds - Hd) ./ (Tref .* Rgas))) ./ ...
                         (1 + exp((Tk  .* Ds - Hd) ./ (Rgas .* Tk))));

% =========================================================================
% --- [7] 动力学常数温度响应 / Kinetic constants at leaf temperature ---
% =========================================================================
% Kc：CO₂的 Michaelis-Menten 常数 [bar]；参考值 404.9 μmol mol⁻¹ @ 25°C
% Ko：O₂的 Michaelis-Menten 常数 [bar]；参考值 278.4 mmol mol⁻¹ @ 25°C
% Gamma_star：CO₂补偿点（无暗呼吸）[bar]；参考值 42.75 μmol mol⁻¹ @ 25°C
% tempcor=0 时 Arrhenius 项=1，常数退化为 25°C 值（用于关闭温度响应调试）
Kc         = 404.9e-6 .* p_bar .* f_Arr(79430, T) .* tempcor  + 404.9e-6 .* p_bar .* (1 - tempcor);
Ko         = 278.4e-3 .* p_bar .* f_Arr(36380, T) .* tempcor  + 278.4e-3 .* p_bar .* (1 - tempcor);
Gamma_star = 42.75e-6 .* p_bar .* f_Arr(37830, T) .* tempcor  + 42.75e-6 .* p_bar .* (1 - tempcor);

% =========================================================================
% --- [8] Vcmax / Jmax / TPU / Rd 温度响应 / Temperature responses ---
% =========================================================================
% 公式：X(T) = X25 × f_Arr(Ha) × f_Peak(Hd, ΔS)
% Ha/Hd 参数来自 CLM5 Table 9.1（Bernacchi 2001 / Kattge & Knorr 2007）
%
%           Ha [J/mol]   Hd [J/mol]   ΔS [J/mol/K]
%  Vcmax     72000        200000       DeltaS_Vcmax（动态）
%  Jmax      50000        200000       DeltaS_Jmax （动态）
%  TPU       72000        200000       490（固定）
%  Rd        46390        150650       490（固定）
Vcmax = Vcmax25_base .* f_Arr(72000, T) .* f_Peak(200000, DeltaS_Vcmax, T);  % [μmol m⁻² s⁻¹]
Jmax  = Jmax25       .* f_Arr(50000, T) .* f_Peak(200000, DeltaS_Jmax,  T);  % [μmol m⁻² s⁻¹]
TPU   = TPU25        .* f_Arr(72000, T) .* f_Peak(200000, DeltaS_TPU,   T);  % [μmol m⁻² s⁻¹]
Rd    = Rd25         .* f_Arr(46390, T) .* f_Peak(150650, DeltaS_Rd,    T);  % [μmol m⁻² s⁻¹]

% [v7 2026-05-28] 白天叶暗呼吸光抑制 (Kok effect).
% Source: CLM5 PhotosynthesisMod.F90:1562-1563 (clm5.0.37, light_inhibit branch; verbatim
%   "if (light_inhibit .and. par_z(p,1) > 0._r8) lmr_z(p,iv) = lmr_z(p,iv) * 0.67_r8"; the
%   sunlit/shaded multilayer form is at :3201-3205). Verified in local clm5.0.37 source.
%   CLM5's own comment attributes 0.67 to Lloyd et al. 2010 Biogeosciences 7:1833-1859
%   (doi:10.5194/bg-7-1833-2010) & Metcalfe et al. 2010 New Phytol 187:608-621 (+ P. Reich
%   pers. comm. 2015). [Old comment's "Funct Ecol" and "Metcalfe 2012" were citation errors;
%   corrected 2026-06-15. The Lloyd/Metcalfe PDFs are not in the local library, so 0.67 is
%   cited as CLM5's documented value, not independently re-verified.]
% COMPASS Rd 温度响应参数 (46390/150650/DeltaS_Rd=490) 与 CLM5 line 1239/1248/1249
% 一致, 仅缺此光抑制步. Q = 吸收 PAR [μmol m⁻² s⁻¹] (line 114), >0 即白天有光.
% element-wise: 白天因子 0.67, 夜间 1.0 (用乘法避免 Q/Rd 维度不一致时逻辑索引错位).
Rd = Rd .* (0.67 + 0.33 .* double(Q <= 0));     % [μmol m⁻² s⁻¹] Q>0→×0.67, Q≤0→×1.0

% C4 PEPcase 活性温度响应（Q10=1.8）
qt        = 0.1 .* (T - Tref) .* tempcor;  % Q10 指数因子
kpepcase  = Kpep25 .* exp(log(1.8) .* qt); % C4 PEP 速率常数 [μmol m⁻² s⁻¹]

% =========================================================================
% --- [9] 水力胁迫施加 / Apply hydraulic stress ---
% =========================================================================
% sfactor = VC(ψ_leaf)/gpmax，来自木质部脆弱性曲线（Plant_hydraulics.m）
% 胁迫等比例缩放 Vcmax/Jmax/TPU（酶和电子传输同步受抑）
% Rd 不缩放：维持呼吸在胁迫下仍维持（生理合理，CLM5 也不对 Rd 施加水力胁迫）
stressfactor = sfactor;
Vcmax = Vcmax .* stressfactor;
Jmax  = Jmax  .* stressfactor;
TPU   = TPU   .* stressfactor;

% =========================================================================
% --- [10] PSII 光化学效率与电子传输 / PSII and electron transport ---
% =========================================================================
% po0：最大暗光化学效率（=Fv/Fm），由 Kp/(Kf+Kd+Kp) 估算或外部测定
if isempty(po0)
    po0 = Kp ./ (Kf + Kd + Kp);    % [-]，典型值 ≈ 0.80
end

% Jpot：潜在电子传输速率 [μmol e⁻ m⁻² s⁻¹]
%   = 0.5 × po0 × Q × fPAR
%   0.5：PSII 和 PSI 平均分配入射光子
%   po0：最大量子效率（暗适应）
%   fPAR：叶片对 PAR 的实际吸收比例（已含 Cab 效率，默认1）
Jpot = 0.5 .* po0 .* Q .* fPAR;    % [μmol e⁻ m⁻² s⁻¹]

% Jact：Jmax 共限制后的实际电子传输速率 [μmol e⁻ m⁻² s⁻¹]
%   求解：θ·Jact² - (Jpot+Jmax)·Jact + Jpot·Jmax = 0，θ=0.7 曲率因子，取小根
thetaJ = 0.7;
Jact   = sel_root(thetaJ, -(Jpot + Jmax), Jpot .* Jmax, -1);
Jact   = real(Jact);    % 消除强光/低温极端条件下的微小虚部

% Je：荧光计算专用电子流 [μmol e⁻ m⁻² s⁻¹]
%   荧光来自 PSII 捕光端（Jpot），不受 Jmax 的下游限制；Jact 才是光合有效值
Je = Jpot;

% =========================================================================
% --- [11] 湿度与 VPD 计算 / Humidity and VPD ---
% =========================================================================
% 使用 Flatau 公式计算饱和水汽压（精度优于 Magnus，分水/冰相处理）
es_hPa = reshape(equations.Calc_Esat_Flatau(Tc, p), size(T)) * 0.01;  % [hPa]
RH     = min(1, eb ./ max(es_hPa, 1e-6));    % 相对湿度 [-]
VPD    = max(0.05, (es_hPa - eb) * 0.1);     % 气孔端 VPD [kPa]；下限 0.05 防零除

% =========================================================================
% --- [12] Michaelis-Menten 常数（C3）/ Michaelis-Menten constants ---
% =========================================================================
if strcmpi(Type, 'C3')
    MM_consts = Kc .* (1 + O_bar ./ max(Ko, 1e-12));   % Km = Kc(1+O/Ko) [bar]
    minCi = 0.3;    % Ci 最小值约束：≥ 0.3×Cs（Rubisco 限制下的物理下限）
else
    MM_consts = zeros(size(T));
    minCi = 0.1;    % C4 的 Ci 可更低（PEP 亲和力强）
end

% =========================================================================
% --- [13] An-gs-Ci 耦合迭代 / Coupled An-gs-Ci iteration ---
% =========================================================================
% 若外部固定 Ci（诊断模式），直接计算一次，不迭代
if ~isempty(Ci_input)
    Ci_ppm = Ci_input;   % 单位约定：ppm（μmol mol⁻¹）
    [A, Ag, Vc, Ve, Vs, CO2_per_electron] = compute_assimilation(Ci_ppm .* ppm2bar);
else
    % 阻尼迭代：初始猜测 Ci = 0.7×Cs（C3 典型 Ci/Ca 比）
    Ci_ppm  = max(minCi .* Cs, 0.7 .* Cs);
    max_iter = 50;
    tol_pa   = 0.01;    % 收敛阈值 [Pa]
    damp     = 0.5;     % 阻尼系数（防止振荡）
    for iter_idx = 1:max_iter
        % a. 光合需求：An(Ci)
        [A, Ag, Vc, Ve, Vs, CO2_per_electron] = compute_assimilation(Ci_ppm .* ppm2bar);
        % b. 气孔供给：Medlyn gs0 + optional heat-induced g_heat (enable default 0)
        gs0 = medlyn_gs(Cs, A, VPD, g0_medlyn, g1_medlyn, gs_model_type);
        g_heat = heat_induced_gs(Tc, VPD, A, sfactor, heat_gs);
        gs = apply_heat_gs_cap(gs0, g_heat, heat_gs);
        % c. 更新 Ci：Ci = Cs - 1.6×An/gs（Fick 扩散，Cs 已含边界层修正）
        Ci_next = max(minCi .* Cs, Cs - 1.6 .* A ./ max(gs, 1e-9));
        Ci_next = min(Cs, Ci_next);     % Ci 不能超过 Cs
        delta_pa = abs((Ci_next - Ci_ppm) .* ppm2bar .* 1e5);   % [Pa]
        Ci_ppm   = Ci_ppm + damp .* (Ci_next - Ci_ppm);
        if all(delta_pa(:) < tol_pa | isnan(delta_pa(:)))
            break
        end
    end
end
fcount = fcount + 1;

% =========================================================================
% --- [14] 最终气孔导度与阻力 / Final stomatal conductance and resistance ---
% =========================================================================
gs0 = medlyn_gs(Cs, A, VPD, g0_medlyn, g1_medlyn, gs_model_type);  % [mol m⁻² s⁻¹] pure Medlyn
g_heat = heat_induced_gs(Tc, VPD, A, sfactor, heat_gs);             % [mol m⁻² s⁻¹] 0 when disabled
gs = apply_heat_gs_cap(gs0, g_heat, heat_gs);                      % [mol m⁻² s⁻¹]
rcw = (rhoa ./ (Mair * 1e-3)) ./ max(gs, 1e-9);    % 气孔阻力 [s m⁻¹]；rhoa/Mair≈41.6 mol m⁻³
gs_from_rcw = (rhoa ./ (Mair * 1e-3)) ./ max(rcw, 1e-9); % [mol m-2 s-1] diagnostic inverse of rcw
% Identity of Medlyn core only (g_heat is intentionally separate and may be >0 when enabled)
gs_medlyn_residual = gs0 - medlyn_gs(Cs, A, VPD, g0_medlyn, g1_medlyn, gs_model_type); % [mol m-2 s-1]

% =========================================================================
% --- [15] 荧光计算 / Fluorescence calculation ---
% =========================================================================
% Ja：实际光合消耗的电子流，由总光合 Ag 反算 [μmol e⁻ m⁻² s⁻¹]
Ja = Ag ./ max(CO2_per_electron, 1e-12);

% ps：PSII 开放反应中心比例（光化学效率）[-]
%   ps = po0 × Ja/Je；当 Je<1 μmol/m²/s（夜间或黎明/黄昏极弱光）强制 ps=po0（暗适应态）
%   防止：弱光下 Ag≈0 → Ja≈0 → ps→0 → ps_rel→1 → 伪最大 NPQ
ps       = po0 .* Ja ./ max(Je, 1e-12);
zero_je  = Je < 1;          % [μmol m⁻² s⁻¹] 弱光阈值
ps(zero_je) = po0(zero_je); % 暗适应态：ps=po0，NPQ=0
ps       = min(max(ps, 0), po0);

% ps_rel：相对光化学效率下降量（NPQ 驱动量）[-]
ps_rel   = max(0, 1 - ps ./ max(po0, 1e-12));

% 荧光速率常数与产额（van der Tol 2014）
[eta, qE, qQ, fs, fo, fm, fo0, fm0, Kn] = Fluorescencemodel(ps, ps_rel, Kp, Kf, Kd, Knparams);
Kpa = ps ./ max(fs, 1e-12) .* Kf;  % 实际光化学速率常数 [-]

% =========================================================================
% --- [16] 叶肉导度修正（可选）/ Mesophyll conductance correction ---
% =========================================================================
Cc = [];
if any(g_m(:) < Inf)
    Cc = (Ci_ppm .* ppm2bar - A ./ g_m) ./ ppm2bar;   % 叶绿体 CO₂ [ppm]
end

% =========================================================================
% --- [17] 输出打包 / Pack outputs ---
% =========================================================================
biochem_out.A      = A;         % 净光合 [μmol CO₂ m⁻² s⁻¹]
biochem_out.Ag     = Ag;        % 总光合 [μmol CO₂ m⁻² s⁻¹]
biochem_out.Ci     = Ci_ppm;    % 胞间 CO₂ [ppm]
if ~isempty(Cc)
    biochem_out.Cc = Cc;        % 叶绿体 CO₂ [ppm]
end
biochem_out.rcw      = rcw;     % 气孔阻力 [s m⁻¹]
biochem_out.gs       = gs;      % 气孔导度 [mol m⁻² s⁻¹] (gs0 + g_heat, capped)
biochem_out.gs0_medlyn = gs0;   % [mol m-2 s-1] pure Medlyn before heat term
biochem_out.g_heat   = g_heat;  % [mol m-2 s-1] heat-induced additive term (0 if disabled)
biochem_out.gs_from_rcw = gs_from_rcw; % [mol m-2 s-1] diagnostic inverse of rcw
biochem_out.gs_medlyn_residual = gs_medlyn_residual; % [mol m-2 s-1] online Medlyn identity check
biochem_out.Cs       = Cs;      % [ppm] leaf-surface CO2 used by Medlyn
biochem_out.O2a_Pa   = O2a_Pa;  % [Pa] atmospheric O2 partial pressure used by C3 kinetics
biochem_out.VPD      = VPD;     % [kPa] leaf-surface vapor pressure deficit used by Medlyn
biochem_out.Tleaf    = Tc;      % [degC] leaf temperature
biochem_out.Q        = Q;       % [umol m-2 s-1] absorbed PAR passed to biochemistry
biochem_out.sfactor  = stressfactor; % [-] hydraulic stress factor applied to capacities
biochem_out.g0_medlyn = g0_medlyn; % [mol m-2 s-1] Medlyn intercept
biochem_out.g1_medlyn = g1_medlyn; % [kPa^0.5] Medlyn slope
biochem_out.RH       = RH;      % 相对湿度 [-]
biochem_out.warnings = warnings;
biochem_out.fcount   = fcount;
biochem_out.Tgrowth = Tgrowth;  % [degC] growth temperature used for thermal acclimation
biochem_out.Tgrowth_clamped = Tgrowth_clamped; % [degC] clamped growth temperature
biochem_out.DeltaS_Vcmax = DeltaS_Vcmax; % [J mol-1 K-1] Vcmax deactivation entropy
biochem_out.DeltaS_Jmax = DeltaS_Jmax;   % [J mol-1 K-1] Jmax deactivation entropy
biochem_out.Topt_Vcmax = Topt_Vcmax;     % [degC] diagnostic Vcmax optimum
biochem_out.Topt_Jmax = Topt_Jmax;       % [degC] diagnostic Jmax optimum
biochem_out.Vcmax25_base = Vcmax25_base; % [umol m-2 s-1] current 25C Vcmax capacity
biochem_out.Jmax25 = Jmax25;             % [umol m-2 s-1] current 25C Jmax capacity
biochem_out.Rd25 = Rd25;                 % [umol m-2 s-1] current 25C dark respiration capacity
biochem_out.Vcmax    = Vcmax;   % 温度修正 Vcmax [μmol m⁻² s⁻¹]
biochem_out.Jmax     = Jmax;    % [umol m-2 s-1] temperature- and stress-adjusted Jmax
biochem_out.Vc       = Vc;      % Rubisco 限制速率 [μmol m⁻² s⁻¹]
biochem_out.Ve       = Ve;      % 光限制速率 [μmol m⁻² s⁻¹]
biochem_out.Vs       = Vs;      % TPU 限制速率 [μmol m⁻² s⁻¹]
biochem_out.Rd       = Rd;      % 暗呼吸 [μmol m⁻² s⁻¹]
biochem_out.Ja       = Ja;      % 实际光合电子流 [μmol e⁻ m⁻² s⁻¹]
biochem_out.ps       = ps;      % PSII 开放中心比例 [-]
biochem_out.ps_rel   = ps_rel;  % 相对光化学下降量 [-]
biochem_out.Kd       = Kd;      % 热耗散速率常数 [-]
biochem_out.Kn       = Kn;      % NPQ 速率常数 [-]
biochem_out.NPQ      = Kn ./ (Kf + Kd);    % 非光化学猝灭 NPQ [-]
biochem_out.Kf       = Kf;      % 荧光速率常数 [-]
biochem_out.Kp0      = Kp;      % 最大光化学速率常数 [-]
biochem_out.Kp       = Kpa;     % 实际光化学速率常数 [-]
biochem_out.eta      = eta;     % 荧光相对量 fs/fo0 [-]
biochem_out.qE       = qE;      % 能量猝灭 qE [-]
biochem_out.fs       = fs;      % 稳态荧光产额 [-]
biochem_out.ft       = fs;      % 同 fs（向后兼容别名）
biochem_out.SIF      = fs .* Q; % 太阳诱导荧光 SIF [μmol m⁻² s⁻¹]
biochem_out.fo0      = fo0;     % 暗适应最小荧光 [-]
biochem_out.fm0      = fm0;     % 暗适应最大荧光 [-]
biochem_out.fo       = fo;      % 光适应最小荧光 [-]
biochem_out.fm       = fm;      % 光适应最大荧光 [-]
biochem_out.Fm_Fo    = fm ./ max(fo,  1e-12);   % Fm/Fo 比值 [-]
biochem_out.Ft_Fo    = fs ./ max(fo,  1e-12);   % Ft/Fo 比值 [-]
biochem_out.qQ       = qQ;      % 光化学猝灭 qQ [-]
return

% =========================================================================
% 嵌套函数：compute_assimilation
% Nested function: compute_assimilation
% =========================================================================
    function [A_loc, Ag_loc, Vc_loc, Ve_loc, Vs_loc, cpe_loc] = compute_assimilation(Ci_bar)
    % 输入：Ci_bar [bar]；输出：A_loc, Ag_loc [μmol m⁻² s⁻¹]
    % 捕获外部变量：Type, Vcmax, Jact, Gamma_star, MM_consts, TPU, Rd, effcon, g_m
        if strcmpi(Type, 'C3')
            if any(g_m(:) < Inf)
                % 叶肉导度修正（g_m 有限时 Cc ≠ Ci）
                Vc_loc  = sel_root(1./g_m, -(MM_consts+Ci_bar+(Rd+Vcmax)./g_m),  Vcmax.*(Ci_bar-Gamma_star+Rd./g_m),     -1);
                Ve_loc  = sel_root(1./g_m, -(Ci_bar+2.*Gamma_star+(Rd+Jact.*effcon)./g_m), Jact.*effcon.*(Ci_bar-Gamma_star+Rd./g_m), -1);
                cpe_loc = Ve_loc ./ max(Jact, 1e-12);
            else
                % 标准 Farquhar C3（g_m=Inf，Cc=Ci）
                % Vc：Rubisco 限制  [μmol m⁻² s⁻¹]
                Vc_loc  = Vcmax .* (Ci_bar - Gamma_star) ./ max(MM_consts + Ci_bar, 1e-12);
                % cpe：CO₂-电子转换效率 = effcon × (Ci-Γ*)/(Ci+2Γ*)
                %      effcon=1/4（Farquhar 1980，与 CLM5 Jmax25 参数体系一致）
                cpe_loc = ((Ci_bar - Gamma_star) ./ max(Ci_bar + 2.*Gamma_star, 1e-12)) .* effcon;
                % Ve：光限制速率（Jact 已含 Jmax 共限制）[μmol m⁻² s⁻¹]
                Ve_loc  = Jact .* cpe_loc;
            end
            % Vs：TPU 限制速率 [μmol m⁻² s⁻¹]
            Vs_loc = 3 .* TPU;
            % 第一步平滑共限制：min_smooth(Vc, Ve)，曲率 0.98（Collatz 1991 / CLM5）
            V_loc   = sel_root(0.98, -(Vc_loc+Ve_loc), Vc_loc.*Ve_loc, sign(Gamma_star-Ci_bar));
            % 第二步平滑共限制：min_smooth(V, Vs)，曲率 0.95（Collatz 1991 / CLM5）
            Ag_loc  = sel_root(0.95, -(V_loc+Vs_loc), V_loc.*Vs_loc, -1);
        else
            % C4 branch: keep the original SCOPE-style Collatz implementation
            % to avoid mixing CLM5 C4 empirical limits with the SCOPE
            % fluorescence/electron-transport pathway used below.
            Vc_loc  = Vcmax;
            Vs_loc  = kpepcase .* Ci_bar;
            cpe_loc = effcon .* ones(size(Ci_bar));
            Ve_loc  = Jact .* cpe_loc;
            V_loc   = sel_root(0.80, -(Vc_loc + Ve_loc), Vc_loc .* Ve_loc, -1);
            Ag_loc  = sel_root(0.95, -(V_loc + Vs_loc), V_loc .* Vs_loc, -1);
        end
        A_loc = Ag_loc - Rd;    % 净光合 = 总光合 - 暗呼吸
    end
end

% =========================================================================
% 局部函数：medlyn_gs — Medlyn (2011) 最优气孔导度
% Local function: medlyn_gs — Medlyn (2011) optimal stomatal conductance
% =========================================================================
function Topt_C = peaked_arrhenius_topt_c(Ha, Hd, DeltaS, Rgas)
%PEAKED_ARRHENIUS_TOPT_C  Return the peak temperature implied by f_Arr*f_Peak.
%
% Inputs:
%   Ha      [J mol-1]       Activation energy
%   Hd      [J mol-1]       Deactivation energy
%   DeltaS  [J mol-1 K-1]   Entropy term, scalar or array
%   Rgas    [J mol-1 K-1]   Universal gas constant
%
% Output:
%   Topt_C  [degC]          Diagnostic optimum temperature
%
% [DIAG 2026-06-19] Temporary online validation support. Delete after the
% Medlyn/SOX/thermal-acclimation annual diagnostic is archived.
    denom = DeltaS - Rgas .* log(Hd ./ Ha - 1);
    Topt_C = Hd ./ denom - 273.15;
end

% =========================================================================
% Local function: medlyn_gs - Medlyn (2011) optimal stomatal conductance
% =========================================================================
function gs = medlyn_gs(Cs_ppm, A_umol, VPD, g0, g1, gs_model_type)
% PPM-form Medlyn implementation used consistently throughout this file:
% gs = g0 + 1.6 * (1 + g1 / sqrt(VPD)) * A_umol / Cs_ppm
% 参数：g0 [mol m⁻² s⁻¹]，g1 [kPa^0.5]，VPD [kPa]，
%       A_umol [μmol m⁻² s⁻¹]，Cs_ppm [μmol mol⁻¹]
% 由于 A/Cs 的量纲为 mol m⁻² s⁻¹，本式与下游 Ci_ppm 迭代口径保持一致。
if gs_model_type ~= 1
    error('vegetation.biochemical:UnsupportedGsModel', 'Only gs_model_type==1 (Medlyn) is implemented.');
end
Aeff = max(A_umol, 0);   % 夜间 A<0 时 gs 取最小值 g0（不驱动负扩散）
gs   = g0 + 1.6 .* (1 + g1 ./ sqrt(max(VPD, 0.05))) .* Aeff ./ max(Cs_ppm, 1e-9);
gs   = max(g0, gs);
end

function g_heat = heat_induced_gs(Tc, VPD, A_umol, wt, p)
%HEAT_INDUCED_GS  Additive heat-driven stomatal opening after Medlyn (optional path C).
%
%   g_heat = g_max * phi_T * phi_w * phi_D * phi_A  [mol m-2 s-1]
%   When p.enable is false (production default), returns zeros(size(A_umol)).
%
%   Gates (linear ramps; docs/plant/stomata_heat_induced_gs_design.md):
%     phi_T: ramp(Tc; T_on, T_full) [degC]
%     phi_w: ramp(wt; wt_off, wt_on)  [-] hydraulic sfactor
%     phi_D: 1 - ramp(D; D_lo, D_hi)  [kPa] suppresses extreme-VPD opening
%     phi_A: ramp(An; A_lo, A_hi)     [umol m-2 s-1]
%
%   Never raises Medlyn g0; drought (low wt) hard-blocks g_heat.
% Editor: Jizixing Zhu; Date: 2026-07-20.
    if ~p.enable
        g_heat = zeros(size(A_umol));
        return
    end
    An_pos = max(A_umol, 0);
    phi_T = linear_ramp(Tc, p.T_on, p.T_full);
    phi_w = linear_ramp(wt, p.wt_off, p.wt_on);
    phi_D = 1 - linear_ramp(VPD, p.D_lo, p.D_hi);
    phi_A = linear_ramp(An_pos, p.A_lo, p.A_hi);
    g_heat = p.g_max .* phi_T .* phi_w .* phi_D .* phi_A;
    g_heat = max(g_heat, 0);
end

function gs = apply_heat_gs_cap(gs0, g_heat, p)
% Apply additive heat term; apply gs_cap only when heat path is enabled.
% When disabled, gs = gs0 unchanged (production numerical identity).
% Editor: Jizixing Zhu; Date: 2026-07-20.
    if ~p.enable
        gs = gs0;
        return
    end
    gs = gs0 + g_heat;
    if isfinite(p.gs_cap) && p.gs_cap > 0
        gs = min(gs, p.gs_cap);
    end
    gs = max(gs, 0);
end

function y = linear_ramp(x, x0, x1)
% Linear ramp from 0 at x0 to 1 at x1 (elementwise).
% Editor: Jizixing Zhu; Date: 2026-07-20.
    if ~isfinite(x0) || ~isfinite(x1) || x1 <= x0
        y = double(x >= x1);
        return
    end
    y = (x - x0) ./ (x1 - x0);
    y = max(0, min(1, y));
end

function flag = local_heat_gs_flag(s, name, default_value)
    flag = local_heat_gs_num(s, name, default_value) ~= 0;
end

function v = local_heat_gs_num(s, name, default_value)
    if isfield(s, name) && ~isempty(s.(name)) && isfinite(double(s.(name)))
        v = double(s.(name));
    else
        v = default_value;
    end
end

% =========================================================================
% 局部函数：sel_root — 平滑共限制二次方程求根
% Local function: sel_root — smooth co-limitation quadratic root
% =========================================================================
function x = sel_root(a, b, c, dsign)
% 求解 a·x² + b·x + c = 0，dsign=-1 取小根，dsign=+1 取大根
% a=0 时退化为线性：x = -c/b
if isscalar(a) && a == 0
    x = -c ./ b;
else
    dsign(dsign == 0) = -1;
    disc = max(b.^2 - 4 .* a .* c, 0);     % 钳零防止数值误差导致虚根
    x    = (-b + dsign .* sqrt(disc)) ./ (2 .* a);
end
end

% =========================================================================
% 局部函数：Fluorescencemodel — van der Tol (2014) 荧光模型
% Local function: Fluorescencemodel — van der Tol (2014) fluorescence
% =========================================================================
function [eta, qE, qQ, fs, fo, fm, fo0, fm0, Kn] = Fluorescencemodel(ps, x, Kp, Kf, Kd, Knparams)
% 输入：ps [-] PSII开放比例；x=ps_rel [-] 相对光化学下降量
% Kno, alpha, beta：NPQ 经验参数（van der Tol 2014）
Kno   = Knparams(1);
alpha = Knparams(2);
beta  = Knparams(3);

% NPQ 速率常数 Kn（Hill 方程形式）
x_alpha = exp(log(max(x, 1e-12)) .* alpha);
Kn      = Kno * (1 + beta) .* x_alpha ./ (beta + x_alpha);

% 各状态荧光产额
fo0 = Kf ./ (Kf + Kp + Kd);        % 暗适应最小荧光（NPQ=0，开放中心）
fo  = Kf ./ (Kf + Kp + Kd + Kn);   % 光适应最小荧光（含NPQ，开放中心）
fm  = Kf ./ (Kf + Kd + Kn);        % 光适应最大荧光（含NPQ，关闭中心）
fm0 = Kf ./ (Kf + Kd);             % 暗适应最大荧光（NPQ=0，关闭中心）

fs   = fm .* (1 - ps);              % 稳态荧光产额（混合开/关中心）
eta  = fs ./ fo0;                   % 相对荧光（归一化到暗适应fo0）

qQ   = 1 - (fs - fo)  ./ max(fm  - fo,  1e-12);    % 光化学猝灭
qE   = 1 - (fm - fo)  ./ max(fm0 - fo0, 1e-12);    % 能量猝灭（NPQ指示）
end
