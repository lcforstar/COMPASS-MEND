function [Snit, SVnit, Cvnit_N2o, mm_factor, M_total, eT_matrix, Snit_N2o, moisture_factor] = nitrification(...
    TT, N, NL, Theta_LL, Theta_f, Theta_s, RHO_b_Kd, Vmax_nit, Km_nit, ...
    Snit, R, Tr, ePH, Cvnit_N2o, PH1, SVnit, KT, DeltZ, f_aerobic, ph_response_enabled, matric_head_cm, M_min_nit)
%NITRIFICATION  Nitrification (NH4 → NO3) — DAISY 风格 Michaelis-Menten 显式源汇.
%
%   设计思路：MM 是非线性反应, 不走 SVnit 半隐式 (会因 Picard mm_factor(C_old) vs C_new
%   失配在施肥脉冲后失真), **直接作为源汇 Snit 进 G1/G2** (跟 Shys/Smin/Svol 一样).
%
%   公式来源：DAISY (Hansen 1991), nitrification_soil.C:69
%             rate = k_10 × w_factor × T_factor × M/(K_m+M)   [μg N/(cm³·d)]
%
%   核心方程（土体基, μg N/(cm³·s)）:
%     M_total   = (θ_L + ρ·Kd) · C_NH4_liq                     [μg/cm³ soil]  DAISY 底物
%     mm_factor = M_total / (K_m + M_total)                    [-]            DAISY 饱和
%     Snit      = V_max × eT × f_aerobic × ePH × mm_factor × f_dry [μg N/(cm³·s)]
%   f_dry is the DAISY dry-side matric-head response. Wet-side oxygen
%   limitation remains solely with f_aerobic; no second wetness penalty.
%
%   单位约定 (跟 DAISY 原文):
%     V_max     [μg N/(cm³·s)]  Excel 'Vmax_nit' [μg/(cm³·d)] / 86400
%                               这是 *饱和时的最大反应速率*, 不是一阶系数. Hansen k_10.
%     K_m       [μg/cm³ soil]   Excel 'Km_nit', 总底物半饱和  (M_total 单位)
%
%   ADE 耦合 (assembleSoluteParameters.m):
%     NH4: G1 += -Snit  (硝化是 NH4 的汇, 显式扣除)
%     NO3: G2 += 0.98 × Snit; the remaining 0.02 × Snit is N2O-N.
%   SVnit 输出保留 (=0), 仅用于兼容老接口的 F1 = -SVnit 不再起作用.
%
%   INPUTS
%   ------
%   TT          [°C, NN×1]    节点温度
%   N           [μg/cm³ 液相, NN×1] 液相 NH4 浓度 C_NH4_liq
%   NL          [-]            层数
%   Theta_LL    [m³/m³, NL×2]  液相含水率
%   Theta_f, Theta_s [m³/m³, NL×1] Retained legacy inputs, not used by
%               the dry response: its independent driver is matric_head_cm.
%   RHO_b_Kd    [-, NL×1]      ρ_bulk × K_d (NH4 adsorption buffer)
%   Vmax_nit    [μg/(cm³ soil·s), NL×1]  最大硝化绝对速率 (Excel 'Vmax_nit' [μg/(cm³·d)] /86400)
%                                        DAISY Hansen 1991 k_10 风格. 注: 与 denit 'Km' [1/s] 单位不同
%   Km_nit      [μg/cm³ soil, NL×1] 半饱和总底物 M_total (Excel 'Km_nit', DAISY 风格)
%                                        注: 跟 denit 'KC' [mg/L 液相] 单位不同, 因 nit 底物是 (θ+ρKd)·C 总池子
%   Snit        [NL×2]         前一 Picard iter Snit (留接口, 实际不读)
%   R           [J/mol/K]      气体常数
%   Tr          [°C]           参考温度 (20°C)
%   ePH         [-, NL×2]      pH 因子 (覆盖输出)
%   Cvnit_N2o   [KT×1]         Column N2O-N rate history (inout), not a time integral
%   PH1         [-, NN×1]      节点 pH
%   SVnit       [NL×2]         前一 Picard iter SVnit (留接口, 覆盖输出)
%   KT          [-]            当前时步
%   DeltZ       [cm, NL×1]     层厚
%   f_aerobic   [-, NL×2]      好氧分配因子 (从 solveO2Transport)
%   ph_response_enabled [optional logical scalar] Default true retains the
%               reference pH response; false sets only its rate factor to 1.
%               PH1 remains the actual soil pH used by other chemistry.
%   matric_head_cm [optional real finite numeric NN-vector, cm]
%               Converged carrier-tile matric head, bottom-to-surface.
%               Positive/zero head is saturated. Omitted or [] retains the
%               legacy factor of one for direct/reference callers. Active
%               production O2 and neutral replay paths supply actual heads.
%   M_min_nit [optional nonnegative finite real numeric NL-vector]
%               Bottom-up protected total NH4 [ug N/cm3 soil]. Empty or zero
%               retains the original arithmetic; positive thresholds expose
%               only max(M_total-M_min_nit,0) to the kinetic equation. A zero
%               scalar is accepted for direct-call compatibility. The caller
%               must additionally limit assembled shared-node reaction loads.
%
%   OUTPUTS
%   -------
%   Snit        [NL×2]  Total NH4-N reaction rate [μg N/(cm³·s)],
%               before the caller's conservative substrate/O2 limiting.
%   SVnit       [NL×2]  Zero [1/s]; no duplicate semi-implicit NH4 loss.
%   Cvnit_N2o   [KT×1]  Column N2O-N rate [μg N/(cm²·s)] at index KT.
%   mm_factor   [NL×2]  MM 饱和因子 [0,1]            诊断
%   M_total     [NL×2]  总底物 [μg/cm³ soil]         诊断 (M = (θ+ρKd)·C, DAISY 风格底物可视化)
%   eT_matrix   [NL×2]  Dimensionless endpoint temperature response.
%   Snit_N2o    [NL×2]  0.02 × Snit [μg N/(cm³·s)].
%   moisture_factor [NL×2] Dimensionless dry-side factor, node-to-endpoint
%               mapping [head(1:NL), head(2:NL+1)] without reversing order.
%   Invalid nonempty matric heads raise nitrification:matricHeadShape or
%   nitrificationMoistureFactor:head; no missing production head is inferred.
%
%   References:
%     Hydrus-1D Technical Manual §3.5 — Michaelis-Menten on liquid + sorbed coupling
%     Hansen S et al. 1991 — DAISY 硝化半饱和动力学
%     Harris R F (1982) — Arrhenius EA=81171 J/mol
%
%   See also: soil.gas.solveO2Transport, soil.nitrogen.chemistry.denitrification
% Change: 2026-10-23 Wenjie Liu -- make the biological pH multiplier optional;
% preserve the reference equation exactly when the option is absent or true.
% Change: 2026-10-28 Wenjie Liu -- adapt the DAISY 7.1.0 dry-side water-potential
% response; retain explicit O2 wet-side control, legacy direct-call parity,
% and the shared NH4/NO3/N2O reaction stoichiometry.
% Editor: Wenjie Liu; Date: 2026-11-01.
% Change: 2026-11-01 Wenjie Liu -- add optional protected NH4 substrate while
% preserving the original zero-threshold MM arithmetic and total diagnostic.

if nargin<20,ph_response_enabled=true;end
if nargin<21,matric_head_cm=[];end
if nargin<22,M_min_nit=[];end
assert(islogical(ph_response_enabled)&&isscalar(ph_response_enabled), ...
    'nitrification:phResponse','ph_response_enabled must be a logical scalar.');
thresholdEnabled=false;
if ~isempty(M_min_nit)
    assert(isnumeric(M_min_nit)&&isreal(M_min_nit)&&isvector(M_min_nit) ...
        &&(numel(M_min_nit)==NL||(isscalar(M_min_nit)&&M_min_nit==0)) ...
        &&all(isfinite(M_min_nit(:)))&&all(M_min_nit(:)>=0), ...
        'nitrification:protectedSubstrate', ...
        'M_min_nit must be a finite nonnegative NL-vector or zero scalar.');
    thresholdEnabled=any(M_min_nit(:)>0);
end

% Activation energy for nitrification [J/mol]   Harris R F (1982)
EA = 81171;

% --- 0. Force column vectors ---------------------------------------------
TT       = TT(:);
N        = N(:);
RHO_b_Kd = RHO_b_Kd(:);
Vmax_nit = Vmax_nit(:);
Km_nit   = Km_nit(:);
PH1      = PH1(:);
DeltZ    = DeltZ(:);

% A node head is mapped onto the two endpoints of each bottom-up element.
moisture_factor = ones(NL,2);
if ~isempty(matric_head_cm)
    assert(isvector(matric_head_cm)&&numel(matric_head_cm)==NL+1, ...
        'nitrification:matricHeadShape','matric_head_cm must contain NL+1 node heads.');
    moisture_nodes = soil.nitrogen.chemistry.nitrificationMoistureFactor(matric_head_cm);
    moisture_factor = [moisture_nodes(1:NL),moisture_nodes(2:NL+1)];
end

% --- 1. Temperature factor: Arrhenius ------------------------------------
TD_nodes  = (TT - Tr) ./ R ./ (TT + 273.15) ./ (Tr + 273.15);
eT_nodes  = exp(EA * TD_nodes);
eT_matrix = [eT_nodes(1:NL), eT_nodes(2:NL+1)];   % [NL x 2]
N_matrix  = [N(1:NL), N(2:NL+1)];                  % [NL x 2]  C_NH4_liq

% --- 2. pH factor (Stange 2007 风格 + COMPASS 历史一致) ------------------
if ph_response_enabled
    ePH_col = 0.56 + atan(pi() * 0.45 * (-5 + min(PH1, 14))) / pi();
else
    ePH_col = ones(NL,1);
end
ePH     = [ePH_col, ePH_col];   % [NL x 2]

% --- 3. DAISY-style Michaelis-Menten 显式源汇 (走 G1/G2, 不走 SVnit 半隐式) ----------
% 物理:
%   M_total   = (θ+ρKd) · C_NH4_liq            DAISY 底物 [μg/cm³ soil]
%   mm_factor = M_total / (K_m + M_total)      DAISY 饱和因子 [0,1]
%   Snit      = V_max × eT × f_aero × ePH × mm_factor × f_dry [μg/(cm³·s)]
%
% 为什么不走 SVnit 半隐式 (跟 denit 不同):
%   denit 用 SVden×C 半隐式时, mm_factor(C_old) 跟 C_new 失配的影响小 (NO3 变化慢);
%   nit 在施肥脉冲后 NH4 突变 5-10×, mm_factor(C_old) 严重低估 C_new 时的 MM 速率
%   → 直接显式 G1=-Snit 物理清晰, 步长足够小 (1h outer + Picard sub) 不会越扣到负.
RHO_b_Kd_ext = [RHO_b_Kd, RHO_b_Kd];                          % [NL x 2]
Vmax_ext     = [Vmax_nit, Vmax_nit];                          % [NL x 2]  μg/(cm³·s)
Km_ext       = [Km_nit,   Km_nit];                            % [NL x 2]  μg/cm³ soil

M_total   = (Theta_LL + RHO_b_Kd_ext) .* N_matrix;            % [NL x 2]  μg/cm³ soil
if thresholdEnabled
    M_min_nit=M_min_nit(:);
    M_available=max(M_total-[M_min_nit,M_min_nit],0);
    mm_factor=M_available./(Km_ext+M_available+eps);
else
    mm_factor = M_total ./ (Km_ext + M_total + eps);          % Original zero-threshold path.
end

% Snit 定义: NH4 总消耗速率 (100% × MM rate), 含 NO3 + N2O 两条路径
% 物理: NH4 → 0.98 NO3 + 0.02 N2O (Bessou 2010)
%   G1 (NH4 ADE) = ... - Snit        (NH4 总消耗)
%   G2 (NO3 ADE) = + 0.98 × Snit ... (98% 入 NO3 池)
%   Cvnit_N2o    = 0.02 × Snit       (2% 逸出气相, 不入水溶池)
Snit      = Vmax_ext .* eT_matrix .* f_aerobic .* ePH .* mm_factor;   % [μg/(cm³·s)] NH4 总消耗
if ~isempty(matric_head_cm)
    Snit = Snit .* moisture_factor;
end
Snit_N2o  = 0.02 * Snit;                                              % [μg/(cm³·s)] N2O 部分

% SVnit 闲置 (=0): NH4 ADE 的 F1=-SVnit 不再扣硝化, Snit 完全走 G1 显式源汇.
SVnit = zeros(NL, 2);

% --- 4. N2O 深度积分诊断 ---------------------------------------------------
Snit_N2o_avg  = (Snit_N2o(:,1) + Snit_N2o(:,2)) * 0.5;
Cvnit_N2o(KT) = DeltZ' * Snit_N2o_avg;                                % [μg/(cm²·s)]

end
