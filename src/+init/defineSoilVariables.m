function SoilVariables = defineSoilVariables(state, SoilProperties)
%{
        Create a SoilVariables structure and define some variables.
        SoilVariables is input for initializing Temperature, Matric potential
        and soil air pressure.
%}

% Create a SoilVariables structure
SoilVariables = struct();

% 从 state 复制初始状态字段到 SoilVariables
soil_fields = { ...
    'P_g',      'P_gg',     ...  % 气压 (NN×1)
    'h',        'hh',       ...  % 基质水头 (NN×1)
    'h_frez',   'hh_frez',  ...  % 冻结压降水头 (NN×1)
    'T',        'TT',       ...  % 温度 (NN×1)
    'T_CRIT',   'TT_CRIT',  ...  % 冻结临界温度 (NN×1)
    'Theta_L',  'Theta_LL', ...  % 液态含水量 (NL×2)
    'Theta_U',  'Theta_UU', ...  % 总含水量 (NL×2)
    'Theta_I',  'Theta_II', ...  % 冰相体积 (NL×2)
    'Theta_V',  'Theta_g',  ...  % 气相体积 (NL×2)
    'Se',                   ...  % 饱和度 (NL×2)
    'KL_h',     'KfL_h',    ...  % 导水率 (NL×2)
    'KfL_T',                ...  % 热驱动系数 (NL×2)
    'DTheta_LLh','DTheta_UUh','Ratio_ice', ...  % 导数与冰水比 (NL×2)
    'C_NH4',    'C_NH4_new',  'C_NH4_OLD',  ...  % 铵态氮浓度 (NN×1)
    'C_NO3',    'C_NO3_new',  'C_NO3_OLD',  ...  % 硝态氮浓度 (NN×1)
    'C_Urea',   'C_Urea_new', 'C_Urea_OLD', ...  % 尿素浓度 (NN×1)
    'C_NH4_ice', 'C_NO3_ice',               ...  % 冰相包合 N (NL×2)  -- P0a.5 (2026-05-13)
    'C_NH4_granule'                         ...  % NH4 颗粒肥缓冲池 (NN×1) -- 2026-05-17 dissolution pool
    };
for field = soil_fields
    SoilVariables.(field{1}) = state.(field{1});
end

% Calculate some of SoilVariables
SoilVariables.XK = SoilProperties.XK; % 0.11 This is for silt loam; For sand XK=0.025
SoilVariables.XWILT = SoilProperties.Wilt_point_mc;
SoilVariables.XCAP = SoilProperties.Field_mc;

SoilVariables.VPERSOC = SoilProperties.MSOC_volume; % fraction of soil organic matter

% Add variables to SoilVariables, see issue 141
SoilVariables.POR = SoilProperties.Porosity;
SoilVariables.VPERS = SoilProperties.Sand_fraction .* (1 - SoilVariables.POR);
SoilVariables.VPERSL = SoilProperties. Silt_fraction.* (1 - SoilVariables.POR);
SoilVariables.VPERC = SoilProperties.Clay_fraction .* (1 - SoilVariables.POR);
SoilVariables.Ks = SoilProperties.Ks;

% Index of wetting history of soil which would be assumed as dry at the first with the value of 1
NL = size(SoilProperties.Porosity, 1);
SoilVariables.IH = ones(NL, 1);

% 冻融阻抗因子 (NL+1)，Lundin 1990，固定=7（长武苹果塬面率定值）
SoilVariables.Imped = repmat(7, NL+1, 1);

% 4 区变饱和流模型参数（标量，全层统一）/ 4-zone model scalar parameters
SoilVariables.Ss      = 1.0e-5;  % 比储水率 / Specific storage (1/cm)
SoilVariables.h_trans = 0.1;     % 过渡区宽度 / Transition head (cm)

% 滞后历史参考含水量 (NL×2)，初始化为零 / Hysteresis reference water content
nD = 2;
SoilVariables.XWRE = zeros(NL, nD);

% [PATCH 2026-05-08] HYDRUS Lenhard hysteresis 状态字段
% 由 hysterIni 在主流程入口处填充（此处仅占位避免下游访问空字段）
SoilVariables.RHSW   = zeros(NL, 7);
SoilVariables.RASW   = zeros(NL, 7);
SoilVariables.IPSW   = ones(NL, 1);
SoilVariables.JJH    = zeros(NL, 1);
SoilVariables.MPSW   = 3 * ones(NL, 1);
SoilVariables.PHSW   = zeros(NL, 1);
SoilVariables.SARW   = zeros(NL, 1);
SoilVariables.RSWAW  = ones(NL, 1);
SoilVariables.AlphaW = zeros(NL, 1);
SoilVariables.SARWI  = zeros(NL, 1);

% [PATCH 2026-05-08 v15] HYDRUS dual-porosity 状态字段
% 由 dualPorIni 在主流程入口处填充
SoilVariables.Theta_im     = zeros(NL, nD);
SoilVariables.Theta_im_old = zeros(NL, nD);
SoilVariables.Theta_m_old  = zeros(NL, nD);

% 冻结分区因子 (NN×1)，初始化为1（无冻融）/ Freezing partition factor
SoilVariables.Gamma_hh = ones(NL + 1, 1);

% "OLD" state snapshots for adaptive time-stepping rollback (NN×1)
SoilVariables.hOLD      = state.h(:);
SoilVariables.hOLD_frez = state.hOLD_frez(:);
SoilVariables.TOLD      = state.T(:);
SoilVariables.P_gOLD    = state.P_g(:);

end
