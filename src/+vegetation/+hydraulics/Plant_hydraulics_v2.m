function [Tp_sun, Tp_sha, psil_sun, psil_sha, psir, psis, ...
          W_root_NEW, W_stem_NEW, W_leaf_NEW, ...
          wt_sun, wt_sha, gsr, info] = ...
    Plant_hydraulics_v2(plantHydro, Khh, hh, Tp_sun, Tp_sha, KT, H, ...
                        W_root_OLD, W_stem_OLD, W_leaf_OLD, Delt_t)
%PLANT_HYDRAULICS_V2 瞬态植物水力学版（#7） / Prognostic plant hydraulics
%
% 替代原 Plant_hydraulics.m 的"稳态 SPAC"求解，实现 ClimaLand 公式 32–37
% 的瞬态水量平衡。保留与 SCOPE 兼容的输出接口（psil_sun/sha 等），新增
% 三层储水量 W_root/W_stem/W_leaf 作为 prognostic 变量。
%
% 与 v1 区别 / Differences from v1:
%   - 新增输入：W_root_OLD, W_stem_OLD, W_leaf_OLD (上一时步储水), Delt_t
%   - 新增输出：W_root_NEW, W_stem_NEW, W_leaf_NEW
%   - 内部用 prognosticPlantWater.m 解 IMEX 3×3 系统
%   - psil 来自瞬态 W → ψ 反算（保留 Manzoni VC 计算 wt 兼容性）
%
% 输入 / Inputs:
%   plantHydro     - 植物水力初始化结构体（v1 字段 + 新增）
%       v1 兼容:
%         .gpmax, .p50, .aa     (Manzoni VC)
%         .RAI, .LAI, .h_root, .h_leaf
%       v2 新增 / v2 new:
%         .K_sat        木质部饱和导水率 [m/s]，默认 2e-8
%         .psi_63       63% 损失水势 [m]，默认 -275（也可复用 .p50 转 m）
%         .c_weibull    Weibull 形状参数，默认 4（也可复用 .aa）
%         .a_retention  [3×1, m^-1] 持水曲线斜率，默认 [4.9e-4; 4.9e-4; 4.9e-4]
%         .W_capacity   [3×1, m] 储水容量，默认 [0.003; 0.005; 0.001]
%         .dz_rs        根→茎距离 [m]，默认 H/2
%         .dz_sl        茎→叶距离 [m]，默认 H/2
%         .P_root_dist  根分布密度 [NN×1, 归一化]
%         .K_soil_root  土-根接口导度 [NN×1, m/s]
%   Khh            - 土壤水力传导率 (NL×2, cm/s) — 同 v1
%   hh             - 土壤压力水头 (NN×1, cmH₂O) — 同 v1
%   Tp_sun, Tp_sha - 阳/阴叶蒸腾速率 (Nmsrmn×1, cm/s) — 同 v1
%   KT             - 当前时步索引 — 同 v1
%   H              - 树高 (m) — 同 v1
%   W_root_OLD, W_stem_OLD, W_leaf_OLD - 上一时步储水 [m]
%   Delt_t         - 时步 [s]
%
% 输出 / Outputs:
%   Tp_sun, Tp_sha             - 水力截断后的实际蒸腾 (cm/s)
%   psil_sun, psil_sha         - 阳/阴叶水势 (MPa)，KT 位置覆盖
%   psir                       - 根水势 (MPa)
%   psis                       - 根区加权土水势 (MPa)
%   W_root_NEW, W_stem_NEW, W_leaf_NEW - 新储水 [m]
%   wt_sun, wt_sha             - 水力胁迫因子 ∈ [0,1]（Manzoni VC，兼容 v1）
%   gsr                        - 根吸水率 [m/s]
%   info                       - 诊断结构体 / Diagnostics
%       .water_balance_residual
%       .stress_factor
%       .K_root, .K_stem, .K_leaf
%       .is_hydraulic_cutoff
%
% 物理参考 / Reference:
%   - Deck et al. (2026) ClimaLand JAMES, §2.3.1 Eq. 32–37
%   - Manzoni et al. (2013) Adv. Water Resour., 51, 292–304
%
% Author: Zhu Jizixing
% Date:   2026-05-01

    NN = length(hh);

    % ============================================================
    % Step 1: 解析 v2 参数（向后兼容，缺失时给默认）
    % ============================================================
    K_sat = getfield_default(plantHydro, 'K_sat', 2e-8);
    psi_63_m = getfield_default(plantHydro, 'psi_63', -275);
    c_weibull = getfield_default(plantHydro, 'c_weibull', 4);

    a_retention = getfield_default(plantHydro, 'a_retention', ...
                                   [4.9e-4; 4.9e-4; 4.9e-4]);
    W_capacity  = getfield_default(plantHydro, 'W_capacity', ...
                                   [0.003; 0.005; 0.001]);
    dz_rs = getfield_default(plantHydro, 'dz_rs', max(H/2, 0.5));
    dz_sl = getfield_default(plantHydro, 'dz_sl', max(H/2, 0.5));

    % --- 根分布与土-根接口导度 ---
    P_root_dist = getfield_default(plantHydro, 'P_root_dist', ...
                                   ones(NN,1)/NN);   % 默认均匀
    if isfield(plantHydro, 'K_soil_root')
        K_soil_root = plantHydro.K_soil_root(:);
    else
        % 用 Khh 节点平均估算（cm/s → m/s）
        K_soil_root = zeros(NN, 1);
        if size(Khh, 1) >= 1
            K_node = mean(Khh, 2);   % 单元中点
            % 投到节点（端节点延拓）
            K_soil_root(1) = K_node(1) * 0.01;     % cm/s → m/s
            for i = 2:NN-1
                K_soil_root(i) = mean(K_node(max(1,i-1):min(end,i))) * 0.01;
            end
            K_soil_root(NN) = K_node(end) * 0.01;
        else
            K_soil_root(:) = 1e-7;
        end
    end

    plant_params = struct( ...
        'K_sat',       K_sat, ...
        'psi_63',      psi_63_m, ...
        'c_weibull',   c_weibull, ...
        'a_retention', a_retention, ...
        'W_capacity',  W_capacity, ...
        'dz_rs',       dz_rs, ...
        'dz_sl',       dz_sl);

    % ============================================================
    % Step 2: 单位转换 hh (cm) → psi_soil (m)
    % ============================================================
    psi_soil_profile = hh(:) * 0.01;   % cm → m

    % ============================================================
    % Step 3: 潜在蒸腾（cm/s → m/s，阳+阴叶之和）
    % ============================================================
    Tp_sun_curr = 0;
    Tp_sha_curr = 0;
    if KT >= 1 && KT <= numel(Tp_sun)
        Tp_sun_curr = max(Tp_sun(KT), 0) * 0.01;   % cm/s → m/s
        Tp_sha_curr = max(Tp_sha(KT), 0) * 0.01;
    end
    T_potential = Tp_sun_curr + Tp_sha_curr;

    % ============================================================
    % Step 4: 调用 IMEX 求解器
    % ============================================================
    W_old_vec = [W_root_OLD; W_stem_OLD; W_leaf_OLD];

    [W_new_vec, psi_new_vec, q_fluxes, hyd_info] = ...
        vegetation.hydraulics.prognosticPlantWater( ...
            W_old_vec, psi_soil_profile, P_root_dist, K_soil_root, ...
            T_potential, plant_params, Delt_t);

    W_root_NEW = W_new_vec(1);
    W_stem_NEW = W_new_vec(2);
    W_leaf_NEW = W_new_vec(3);

    % ============================================================
    % Step 5: 输出水势（m → MPa）
    %   1 m H2O ≈ 0.00981 MPa
    % ============================================================
    m_to_MPa = 9.81e-3;
    psir = psi_new_vec(1) * m_to_MPa;
    psis = sum(P_root_dist .* psi_soil_profile) / max(sum(P_root_dist), 1e-12) * m_to_MPa;
    psil_value_MPa = psi_new_vec(3) * m_to_MPa;

    % --- KT 位置覆盖 psil_sun/sha ---
    %   假设阳/阴叶共享 W_leaf；高级版可分开两个 W_leaf
    if KT >= 1 && KT <= numel(psil_sun)
        psil_sun(KT) = psil_value_MPa;
    end
    if KT >= 1 && KT <= numel(psil_sha)
        psil_sha(KT) = psil_value_MPa;
    end

    % ============================================================
    % Step 6: 实际蒸腾 = 水力截断后的总通量按比例分配回阳/阴叶
    % ============================================================
    if T_potential > 1e-15
        ratio = q_fluxes.T_actual / T_potential;
    else
        ratio = 0;
    end
    ratio = max(0, min(1, ratio));
    if KT >= 1 && KT <= numel(Tp_sun)
        Tp_sun(KT) = Tp_sun(KT) * ratio;
        Tp_sha(KT) = Tp_sha(KT) * ratio;
    end

    % ============================================================
    % Step 7: 水力胁迫因子（保留 Manzoni VC 兼容 v1 接口）
    %   wt = VC(psil) / gpmax = 1 / (1 + (psil/p50)^aa)
    % ============================================================
    if isfield(plantHydro, 'gpmax') && isfield(plantHydro, 'p50') && isfield(plantHydro, 'aa')
        wt_sun_val = vegetation.hydraulics.VulnerabilityCurve( ...
                         psil_value_MPa, plantHydro.gpmax, plantHydro.p50, plantHydro.aa) ...
                     / plantHydro.gpmax;
        wt_sha_val = wt_sun_val;   % 共享 W_leaf
    else
        % 缺 Manzoni 参数：直接用 stress_factor（K/K_sat）
        wt_sun_val = hyd_info.stress_factor;
        wt_sha_val = hyd_info.stress_factor;
    end
    wt_sun = max(0, min(1, wt_sun_val));
    wt_sha = max(0, min(1, wt_sha_val));

    % ============================================================
    % Step 8: 根吸水率
    % ============================================================
    gsr = q_fluxes.q_soil_to_root;

    % ============================================================
    % Step 9: 诊断信息
    % ============================================================
    info = struct( ...
        'water_balance_residual', hyd_info.water_balance_residual, ...
        'stress_factor',          hyd_info.stress_factor, ...
        'K_root',                 hyd_info.K_root, ...
        'K_stem',                 hyd_info.K_stem, ...
        'K_leaf',                 hyd_info.K_leaf, ...
        'is_hydraulic_cutoff',    hyd_info.is_hydraulic_cutoff, ...
        'q_soil_to_root',         q_fluxes.q_soil_to_root, ...
        'q_root_to_stem',         q_fluxes.q_root_to_stem, ...
        'q_stem_to_leaf',         q_fluxes.q_stem_to_leaf, ...
        'T_actual',               q_fluxes.T_actual, ...
        'T_potential',            T_potential);

end

% ----------------------------------------------------------------
% 局部辅助 / Local helper
% ----------------------------------------------------------------
function val = getfield_default(s, fname, default_val)
%GETFIELD_DEFAULT 字段不存在时返回默认值 / Get field with default
    if isstruct(s) && isfield(s, fname) && ~isempty(s.(fname))
        val = s.(fname);
    else
        val = default_val;
    end
end
