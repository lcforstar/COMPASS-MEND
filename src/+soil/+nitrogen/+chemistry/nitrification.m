function [Snit, SVnit, Cvnit_N2o, mm_factor, M_total, eT_matrix, Snit_N2o] = nitrification(...
    TT, N, NL, Theta_LL, Theta_f, Theta_s, RHO_b_Kd, Vmax_nit, Km_nit, ...
    Snit, R, Tr, ePH, Cvnit_N2o, PH1, SVnit, KT, DeltZ, f_aerobic)
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
%     Snit      = 0.98 × V_max × eT × f_aerobic × ePH × mm_factor   [μg N/(cm³·s)]
%
%   单位约定 (跟 DAISY 原文):
%     V_max     [μg N/(cm³·s)]  Excel 'Vmax_nit' [μg/(cm³·d)] / 86400
%                               这是 *饱和时的最大反应速率*, 不是一阶系数. Hansen k_10.
%     K_m       [μg/cm³ soil]   Excel 'Km_nit', 总底物半饱和  (M_total 单位)
%
%   ADE 耦合 (assembleSoluteParameters.m):
%     NH4: G1 += -Snit  (硝化是 NH4 的汇, 显式扣除)
%     NO3: G2 += +Snit  (硝化产物进 NO3 池子)
%   SVnit 输出保留 (=0), 仅用于兼容老接口的 F1 = -SVnit 不再起作用.
%
%   INPUTS
%   ------
%   TT          [°C, NN×1]    节点温度
%   N           [μg/cm³ 液相, NN×1] 液相 NH4 浓度 C_NH4_liq
%   NL          [-]            层数
%   Theta_LL    [m³/m³, NL×2]  液相含水率
%   Theta_f, Theta_s            (保留兼容; MM 公式中不再使用)
%   RHO_b_Kd    [g/cm³, NL×1]  ρ_bulk × K_d (NH4 吸附常数)
%   Vmax_nit    [μg/(cm³ soil·s), NL×1]  最大硝化绝对速率 (Excel 'Vmax_nit' [μg/(cm³·d)] /86400)
%                                        DAISY Hansen 1991 k_10 风格. 注: 与 denit 'Km' [1/s] 单位不同
%   Km_nit      [μg/cm³ soil, NL×1] 半饱和总底物 M_total (Excel 'Km_nit', DAISY 风格)
%                                        注: 跟 denit 'KC' [mg/L 液相] 单位不同, 因 nit 底物是 (θ+ρKd)·C 总池子
%   Snit        [NL×2]         前一 Picard iter Snit (留接口, 实际不读)
%   R           [J/mol/K]      气体常数
%   Tr          [°C]           参考温度 (20°C)
%   ePH         [-, NL×2]      pH 因子 (覆盖输出)
%   Cvnit_N2o   [KT×1]         N2O 累计 (inout)
%   PH1         [-, NN×1]      节点 pH
%   SVnit       [NL×2]         前一 Picard iter SVnit (留接口, 覆盖输出)
%   KT          [-]            当前时步
%   DeltZ       [cm, NL×1]     层厚
%   f_aerobic   [-, NL×2]      好氧分配因子 (从 solveO2Transport)
%
%   OUTPUTS
%   -------
%   Snit        [NL×2]  瞬时硝化速率 [μg N/(cm³·s)]  诊断输出，ADE 用 SVnit×C_new 自洽重算
%   SVnit       [NL×2]  一阶等效系数 [1/s]            ADE F1 = -SVnit
%   Cvnit_N2o   [KT×1]  N2O 累计 [μg N/(cm²·s)]
%   mm_factor   [NL×2]  MM 饱和因子 [0,1]            诊断
%   M_total     [NL×2]  总底物 [μg/cm³ soil]         诊断 (M = (θ+ρKd)·C, DAISY 风格底物可视化)
%
%   References:
%     Hydrus-1D Technical Manual §3.5 — Michaelis-Menten on liquid + sorbed coupling
%     Hansen S et al. 1991 — DAISY 硝化半饱和动力学
%     Harris R F (1982) — Arrhenius EA=81171 J/mol
%
%   See also: soil.gas.solveO2Transport, soil.nitrogen.chemistry.denitrification

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

% --- 1. Temperature factor: Arrhenius ------------------------------------
TD_nodes  = (TT - Tr) ./ R ./ (TT + 273.15) ./ (Tr + 273.15);
eT_nodes  = exp(EA * TD_nodes);
eT_matrix = [eT_nodes(1:NL), eT_nodes(2:NL+1)];   % [NL x 2]
N_matrix  = [N(1:NL), N(2:NL+1)];                  % [NL x 2]  C_NH4_liq

% --- 2. pH factor (Stange 2007 风格 + COMPASS 历史一致) ------------------
ePH_col = 0.56 + atan(pi() * 0.45 * (-5 + min(PH1, 14))) / pi();
ePH     = [ePH_col, ePH_col];   % [NL x 2]

% --- 3. DAISY-style Michaelis-Menten 显式源汇 (走 G1/G2, 不走 SVnit 半隐式) ----------
% 物理:
%   M_total   = (θ+ρKd) · C_NH4_liq            DAISY 底物 [μg/cm³ soil]
%   mm_factor = M_total / (K_m + M_total)      DAISY 饱和因子 [0,1]
%   Snit      = 0.98 × V_max × eT × f_aero × ePH × mm_factor   [μg/(cm³·s)]
%
% 为什么不走 SVnit 半隐式 (跟 denit 不同):
%   denit 用 SVden×C 半隐式时, mm_factor(C_old) 跟 C_new 失配的影响小 (NO3 变化慢);
%   nit 在施肥脉冲后 NH4 突变 5-10×, mm_factor(C_old) 严重低估 C_new 时的 MM 速率
%   → 直接显式 G1=-Snit 物理清晰, 步长足够小 (1h outer + Picard sub) 不会越扣到负.
RHO_b_Kd_ext = [RHO_b_Kd, RHO_b_Kd];                          % [NL x 2]
Vmax_ext     = [Vmax_nit, Vmax_nit];                          % [NL x 2]  μg/(cm³·s)
Km_ext       = [Km_nit,   Km_nit];                            % [NL x 2]  μg/cm³ soil

M_total   = (Theta_LL + RHO_b_Kd_ext) .* N_matrix;            % [NL x 2]  μg/cm³ soil
mm_factor = M_total ./ (Km_ext + M_total + eps);              % [NL x 2]  [0,1]

% Snit 定义: NH4 总消耗速率 (100% × MM rate), 含 NO3 + N2O 两条路径
% 物理: NH4 → 0.98 NO3 + 0.02 N2O (Bessou 2010)
%   G1 (NH4 ADE) = ... - Snit        (NH4 总消耗)
%   G2 (NO3 ADE) = + 0.98 × Snit ... (98% 入 NO3 池)
%   Cvnit_N2o    = 0.02 × Snit       (2% 逸出气相, 不入水溶池)
Snit      = Vmax_ext .* eT_matrix .* f_aerobic .* ePH .* mm_factor;   % [μg/(cm³·s)] NH4 总消耗
Snit_N2o  = 0.02 * Snit;                                              % [μg/(cm³·s)] N2O 部分

% SVnit 闲置 (=0): NH4 ADE 的 F1=-SVnit 不再扣硝化, Snit 完全走 G1 显式源汇.
SVnit = zeros(NL, 2);

% --- 4. N2O 深度积分诊断 ---------------------------------------------------
Snit_N2o_avg  = (Snit_N2o(:,1) + Snit_N2o(:,2)) * 0.5;
Cvnit_N2o(KT) = DeltZ' * Snit_N2o_avg;                                % [μg/(cm²·s)]

end
