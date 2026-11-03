function SoilVariables = calculateHydraulicConductivity(SoilVariables, VanGenuchten, ~, L_f, cfg, Delt_t, TransportCoefficient)
% calculateHydraulicConductivity - 计算土壤水力传导率与各相含水量
%                                   Calculate soil hydraulic conductivity and water content
%
% 函数定义 / Function Definition:
%   SoilVariables = calculateHydraulicConductivity(SoilVariables, VanGenuchten, KIT, L_f, cfg)
%
% 输入 / Inputs:
%   SoilVariables - 土壤变量结构体 / Soil variables struct
%                   必需字段: hhU, TT, h, hh_frez, Theta_L, Theta_U, POR, XCAP, Imped, TT_CRIT, Ks
%   VanGenuchten  - van Genuchten 参数结构体 / VG parameter struct (Theta_r, Theta_s, Alpha, n, m, Ks)
%   KIT           - 当前迭代步数 (0=首次猜测) / Current iteration count (0=first guess)
%   L_f           - 冻结潜热 (J/kg) / Latent heat of freezing (J/kg)
%   cfg           - 模型配置对象 / Model settings (legacy ModelSettings)
%
% 输出 / Outputs:
%   SoilVariables - 更新字段 / Updated fields:
%     .Theta_LL, .Theta_UU, .Theta_II, .DTheta_LLh, .DTheta_UUh
%     .Se, .KL_h, .KfL_h, .KfL_T, .Ratio_ice
%     .Gamma_hh                                   ∈[0,1]，液态水对水汽扩散的有效性系数
%
% 使用示例 / Usage Example:
%   SoilVariables = soil.water.conductivity.calculateHydraulicConductivity( ...
%       SoilVariables, VanGenuchten, KIT, L_f, cfg);
%
% 引用 / References:
%   van Genuchten (1980) Soil Sci. Soc. Am. J., 44, 892-898
%   Mualem (1976) Water Resour. Res., 12, 513-522
%   Harlan (1973) Water Resour. Res., 9, 1314-1323
%   Lundin (1990) Cold Reg. Sci. Technol., 18, 47-64
%
% 修改者 / Modified by:   朱寄子星 (Zhu Jizixing)
% 修改日期 / Modified date: 2026-04-10

phys = equations.physicalConstants();
NL   = cfg.Grid.NL;
nD   = cfg.Grid.nD;
if nargin < 6 || isempty(Delt_t)
    Delt_t = 0;
end
if nargin < 7
    TransportCoefficient = [];
end

dallAmicoEnabled = cfg.lFreezing == 1 && cfg.lThmrlefc == 1 && Delt_t > 0 && ...
    getLogicalOption(cfg, 'lDallAmicoPhaseChange', false) && ...
    hasThermalTimeInputs(TransportCoefficient);

if dallAmicoEnabled
    SoilVariables.DallAmicoThermalTime = NaN(NL, nD);
    SoilVariables.DallAmicoThetaLStar = NaN(NL, nD);
    SoilVariables.DallAmicoPhaseSource = NaN(NL, nD);
    SoilVariables.DallAmicoLimiterCode = zeros(NL, nD);
end

for i = 1:NL
    for j = 1:nD
        MN = i + j - 1;

        % 构造节点标量结构体（仿源码逐节点切片）
        % Build per-node scalar SV, mirroring source sliceVector/sliceMatrix
        SV.hh      = SoilVariables.hhU(MN);
        SV.hh_prev = SoilVariables.hU(MN);   % COR 修正后的上一时步水头
        SV.hh_frez = SoilVariables.hh_frez(MN);   % 冻结压降水头，时步末（lSFCC=1 路径使用）
        SV.h_frez  = SoilVariables.h_frez(MN);    % 冻结压降水头，时步初（导数差分使用）
        SV.TT      = SoilVariables.TT(MN);
        SV.POR     = SoilVariables.POR(i);
        SV.XCAP    = SoilVariables.XCAP(i);
        SV.Ks      = SoilVariables.Ks(i);
        SV.Imped   = SoilVariables.Imped(i);
        SV.Ss      = SoilVariables.Ss;
        SV.h_trans = SoilVariables.h_trans;
        SV.TT_CRIT = SoilVariables.TT_CRIT(MN);
        SV.Theta_L  = SoilVariables.Theta_L(i, j);
        SV.Theta_U  = SoilVariables.Theta_U(i, j);
        SV.Theta_II = SoilVariables.Theta_II(i, j);

        % --- fixHeat：极端水头截断 / Clamp head to valid range ---
        SV.hh = soil.water.conductivity.fixHeat(SV.hh, cfg);

        % --- Gamma_hh：液态水有效性系数（水汽扩散）/ Liquid water availability factor for vapor diffusion ---
        % 基于基质水头干湿程度：|hh|<9899cm → 1（湿润）；|hh|≥1e7cm → 0（绝干）；中间对数插值
        % Based on matric head wetness: moist→1, oven-dry→0, log interpolation in between
        Gamma_hh = soil.water.conductivity.calculateGamma_hh(SV.hh);

        % --- Theta_m：VG 含水量有效上限（受 Gamma_hh 干湿调制）/ VG water content upper bound (Gamma_hh-weighted) ---
        Theta_m = soil.water.conductivity.calculateTheta_m(Gamma_hh, VanGenuchten, i, SV.POR);

        % --- Theta_UU：总含水量 / Total water content (VG) ---
        Theta_UU = soil.water.conductivity.calculateTheta_UU(Theta_m, Gamma_hh, SV, VanGenuchten, i, cfg);

        % --- Theta_II：冰相体积分率（lSFCC=0 正弦曲线初算，后由质量守恒覆盖）/ Volumetric ice content ---
        Theta_II = soil.water.conductivity.calculateTheta_II(SV.TT, SV.XCAP, SV.hh, SV.Theta_II, cfg);

        % --- Theta_LL：液态含水量 + 内部有效饱和度 / Liquid content + internal Se for K ---
        [Theta_LL, Se_eff] = soil.water.conductivity.calculateTheta_LL( ...
            Theta_UU, Theta_II, Theta_m, Gamma_hh, SV, VanGenuchten, i, cfg);

        % 液相守恒反算冰相体积（仅冻融启用时）/ Ice volume from mass balance (freezing only)
        if cfg.lFreezing == 0
            Theta_UU = Theta_LL;
            Theta_II = 0;
        end

        if cfg.lFreezing == 1
            Theta_II = max(0, (Theta_UU - Theta_LL) * phys.RHOL / phys.RHOI);
        end

        if dallAmicoEnabled
            rhoC_native = getThermalCapacity(TransportCoefficient, i, j);
            lambda_native = TransportCoefficient.Lambda_eff(i, j);
            phase = soil.water.conductivity.calculateDallAmicoPhaseChange( ...
                Theta_LL, Theta_II, SV.TT, VanGenuchten.Theta_s(i), VanGenuchten.Theta_r(i), ...
                VanGenuchten.Alpha(i), VanGenuchten.n(i), VanGenuchten.m(i), ...
                rhoC_native, lambda_native, cfg.Grid.DeltZ(i), Delt_t);
            Theta_LL = phase.theta_l_next;
            Theta_II = phase.theta_i_next;
            Theta_UU = Theta_LL + Theta_II * phys.RHOI / phys.RHOL;
            SoilVariables.DallAmicoThermalTime(i, j) = phase.thermal_time_s;
            SoilVariables.DallAmicoThetaLStar(i, j) = phase.theta_l_star;
            SoilVariables.DallAmicoPhaseSource(i, j) = phase.source_liquid_s;
            SoilVariables.DallAmicoLimiterCode(i, j) = limiterCode(phase.limiter);
        end

        % --- DTheta_UUh / DTheta_LLh：含水量对 h 的导数 / Derivatives ---
        DTheta_UUh = soil.water.conductivity.calculateDTheta_UUh( ...
            Theta_UU, Theta_m, Theta_LL, Gamma_hh, SV, VanGenuchten, i, cfg);
        DTheta_LLh = soil.water.conductivity.calculateDTheta_LLh( ...
            DTheta_UUh, Theta_m, Theta_UU, Theta_LL, Gamma_hh, SV, VanGenuchten, i, cfg);
        % --- Se_eff：有效饱和度 / Effective saturation ---
        % lFreezing=0：Se_eff 由 calculateTheta_LL 的 4 区模型内部已算好，直接复用
        % lFreezing=1 / lSFCC=1：calculateTheta_LL 返回默认值，需另行计算
        % lSFCC only selects the freezing formulation; when lFreezing=0,
        % keep the no-ice hydraulic path even if cfg.lSFCC is left at 1.
        if cfg.lFreezing == 1
            Se_eff = soil.water.conductivity.calculateSe(Theta_LL, SV.POR, SV, cfg);
        end

        % --- Ratio_ice：冰水比 / Ice-to-water ratio ---
        if Theta_UU ~= 0
            Ratio_ice = phys.RHOI * Theta_II / (phys.RHOL * Theta_UU);
        else
            Ratio_ice = 0;
        end

        % KIT==0 still needs physical conductivity for plant hydraulics.
            % --- KL_h, KfL_h, KfL_T：导水率与冻融修正 ---
        % Conductivity is needed during initialization too; KIT only denotes
        % first-guess iteration status and should not zero the physical K field.
            % 粘度参考值（20°C）/ Reference viscosity at 20°C
        MU_WN = phys.MU_W0 * exp(phys.MU1 / (8.31441 * (20 + 133.3)));
        if SV.TT < -20
            MU_W = 3.71e-2;
        elseif SV.TT > 150
            MU_W = 1.81e-3;
        else
            MU_W = phys.MU_W0 * exp(phys.MU1 / (8.31441 * (SV.TT + 133.3)));
        end
        % [PATCH 2026-05-08] 给 Mualem 液态导水率公式喂 VG 标准 Se = (θ-θr)/(θs-θr)
        % 旧版用 Se_eff = θ/POR (POR=1-bd/2.65)，分母不一致导致饱和段 K 偏小 9 倍
        % 干空气方程 Sa = θ_g/POR 仍用 SoilVariables.Se 写回值（不变）
        theta_r_i = VanGenuchten.Theta_r(i);
        theta_s_i = VanGenuchten.Theta_s(i);
        Se_for_K = (Theta_LL - theta_r_i) / max(theta_s_i - theta_r_i, 1e-12);
        Se_for_K = max(0, min(1, Se_for_K));
        KL_h  = soil.water.conductivity.calculateKL_h(MU_W, MU_WN, Se_for_K, SV.Ks, VanGenuchten.m(i));

        KfL_h = KL_h * 10^(-SV.Imped * Ratio_ice);
        if cfg.lFreezing == 1
            KfL_T = helpers.heaviside1(SV.TT_CRIT - (SV.TT + phys.T_freeze)) * L_f * 1e4 / (phys.g * phys.T_freeze);
        else
            KfL_T = 0;
        end
        % --- 写回 SoilVariables / Write back per-node results ---
        SoilVariables.Theta_LL(i, j)   = Theta_LL;
        SoilVariables.Theta_UU(i, j)   = Theta_UU;
        SoilVariables.Theta_II(i, j)   = Theta_II;
        SoilVariables.DTheta_LLh(i, j) = DTheta_LLh;
        SoilVariables.DTheta_UUh(i, j) = DTheta_UUh;
        SoilVariables.Se(i, j)         = Se_eff;
        SoilVariables.KL_h(i, j)       = KL_h;
        SoilVariables.KfL_h(i, j)      = KfL_h;
        SoilVariables.KfL_T(i, j)      = KfL_T;
        SoilVariables.Ratio_ice(i, j)  = Ratio_ice;
        SoilVariables.Gamma_hh(MN)     = Gamma_hh;
    end
end

end

function tf = getLogicalOption(cfg, name, defaultValue)
tf = defaultValue;
if isstruct(cfg)
    if isfield(cfg, name)
        tf = logical(cfg.(name));
    end
elseif isobject(cfg) && isprop(cfg, name)
    tf = logical(cfg.(name));
end
end

function tf = hasThermalTimeInputs(TransportCoefficient)
tf = false;
if isempty(TransportCoefficient) || ~isfield(TransportCoefficient, 'Lambda_eff')
    return;
end
if ~isfield(TransportCoefficient, 'EHCAP') && ~isfield(TransportCoefficient, 'c_unsat')
    return;
end
lambda = TransportCoefficient.Lambda_eff;
if isempty(lambda) || ~any(isfinite(lambda(:)) & lambda(:) > 0)
    return;
end
if isfield(TransportCoefficient, 'EHCAP')
    rhoC = TransportCoefficient.EHCAP;
else
    rhoC = [];
end
if isempty(rhoC) || ~any(isfinite(rhoC(:)) & rhoC(:) > 0)
    rhoC = TransportCoefficient.c_unsat;
end
tf = ~isempty(rhoC) && any(isfinite(rhoC(:)) & rhoC(:) > 0);
end

function rhoC_native = getThermalCapacity(TransportCoefficient, i, j)
rhoC_native = NaN;
if isfield(TransportCoefficient, 'EHCAP') && ...
        size(TransportCoefficient.EHCAP, 1) >= i && size(TransportCoefficient.EHCAP, 2) >= j
    rhoC_native = TransportCoefficient.EHCAP(i, j);
end
if ~isfinite(rhoC_native) || rhoC_native <= 0
    rhoC_native = TransportCoefficient.c_unsat(i, j);
end
end

function code = limiterCode(limiter)
switch limiter
    case 'freezing-limited'
        code = 1;
    case 'melting-limited'
        code = 2;
    otherwise
        code = 0;
end
end
