function [W_new, psi_new, q_fluxes, info] = prognosticPlantWater( ...
    W_old, psi_soil_profile, P_root_dist, K_soil_root_profile, ...
    T_potential, plant_params, Delt_t)
% prognosticPlantWater  瞬态植物水量平衡（IMEX 半隐式 3×3） / Prognostic plant water balance
%
% 求解方程组 / System (ClimaLand 2026 Eq. 32–37)：
%   d W_root / dt = q_soil_to_root - q_root_to_stem
%   d W_stem / dt = q_root_to_stem - q_stem_to_leaf
%   d W_leaf / dt = q_stem_to_leaf - T_actual
%
% 通量 (Darcy + 重力) / Fluxes:
%   q_{i->j} = -K_eff * ((psi_j - psi_i)/dz_ij + 1)
%   K_eff   = harmonic mean of Weibull conductivities
%
% 持水曲线 (线性) / Retention:
%   psi(W) = (W/W_cap - 1) / a
%   dpsi/dW = 1 / (a * W_cap)
%
% 时间离散 / Time discretization (IMEX semi-implicit):
%   解 3×3 隐式系统 A * dW = R
%   - 系数 K_eff 用 W^n（即上一步水势）评估，避免非线性迭代
%   - dpsi/dW 用线性持水曲线的解析解
%
% 输入 / Inputs:
%   W_old              - 上一时步储水 / Previous step storage [3×1, m]
%                        顺序: [W_root; W_stem; W_leaf]
%   psi_soil_profile   - 土壤剖面水势 / Soil potential profile [NN×1, m]
%   P_root_dist        - 根分布密度 / Root density distribution [NN×1, 归一化]
%   K_soil_root_profile- 土壤-根界面导度 / Soil-root interface K [NN×1, m/s]
%                        来自土壤导水率，用 SoilVariables.KfL_h 节点平均
%   T_potential        - 潜在蒸腾 / Potential transpiration [m/s] (>0)
%   plant_params       - 植物参数结构体 / Plant parameter struct，必备字段：
%     .K_sat       木质部饱和导水率 / Xylem K_sat [m/s]
%     .psi_63      63% 损失水势 / 63% loss potential [m, 负值]
%     .c_weibull   Weibull 形状参数 / Weibull shape [-]
%     .a_retention 持水曲线斜率 [3×1, m^-1]，[a_root; a_stem; a_leaf]
%     .W_capacity  储水容量 [3×1, m]，[W_r_max; W_s_max; W_l_max]
%     .dz_rs       根→茎距离 / Root-to-stem distance [m]
%     .dz_sl       茎→叶距离 / Stem-to-leaf distance [m]
%   Delt_t             - 时步 / Time step [s]
%
% 输出 / Outputs:
%   W_new              - 更新后储水 / Updated storage [3×1, m]
%   psi_new            - 新水势 / New potentials [3×1, m]，[psi_root; psi_stem; psi_leaf]
%   q_fluxes           - 通量结构体 / Flux struct
%     .q_soil_to_root  土→根净通量 [m/s]
%     .q_root_to_stem  根→茎通量 [m/s]
%     .q_stem_to_leaf  茎→叶通量 [m/s]
%     .T_actual        实际蒸腾 [m/s]
%   info               - 诊断结构体 / Diagnostic struct
%     .K_root, K_stem, K_leaf  各层 Weibull 导度 [m/s]
%     .water_balance_residual  水量守恒残差 [m]
%     .stress_factor           水力胁迫因子 ∈ [0,1]
%
% 物理参考 / Reference:
%   - Deck et al. (2026) ClimaLand JAMES, §2.3.1 Eq. 32–37
%   - Lawrence et al. (2019) CLM5 Tech Note, Ch. 11
%   - Christoffersen et al. (2016) GMD, 9, 4227–4255
%
% Author: Zhu Jizixing
% Date:   2026-05-01

    % ============================================================
    % Step 0: 输入校验 / Input validation
    % ============================================================
    if numel(W_old) ~= 3
        error('prognosticPlantWater:DimMismatch', ...
              'W_old 必须为 3×1 [W_root; W_stem; W_leaf]');
    end
    if numel(plant_params.W_capacity) ~= 3 || numel(plant_params.a_retention) ~= 3
        error('prognosticPlantWater:DimMismatch', ...
              'W_capacity 和 a_retention 必须为 3×1 数组');
    end
    if Delt_t <= 0
        error('prognosticPlantWater:InvalidDt', 'Delt_t 必须 > 0');
    end

    W_old = W_old(:);
    Wcap  = plant_params.W_capacity(:);
    a_ret = plant_params.a_retention(:);

    % ============================================================
    % Step 1: 当前水势 (来自 W_old) / Current potentials
    % ============================================================
    psi_old = vegetation.hydraulics.linearRetention(W_old, Wcap, a_ret);

    % ============================================================
    % Step 2: 当前导度 (Weibull) / Current Weibull conductivities
    % ============================================================
    K_root = vegetation.hydraulics.weibullConductance( ...
        psi_old(1), plant_params.K_sat, plant_params.psi_63, plant_params.c_weibull);
    K_stem = vegetation.hydraulics.weibullConductance( ...
        psi_old(2), plant_params.K_sat, plant_params.psi_63, plant_params.c_weibull);
    K_leaf = vegetation.hydraulics.weibullConductance( ...
        psi_old(3), plant_params.K_sat, plant_params.psi_63, plant_params.c_weibull);

    % --- 节点间有效导度 (调和平均) / Inter-node effective conductivities ---
    K_rs = vegetation.hydraulics.harmonicMean(K_root, K_stem);
    K_sl = vegetation.hydraulics.harmonicMean(K_stem, K_leaf);

    dz_rs = plant_params.dz_rs;
    dz_sl = plant_params.dz_sl;
    if dz_rs <= 0, dz_rs = 1.0; end
    if dz_sl <= 0, dz_sl = 1.0; end

    % ============================================================
    % Step 3: 土→根净通量 / Soil-to-root net flux
    %   q_sr = sum( P(z) * K_sr(z) * (psi_soil(z) - psi_root) )
    %   注意：根分布积分已通过 P_root_dist 实现
    % ============================================================
    P_root_dist = P_root_dist(:);
    psi_soil_profile = psi_soil_profile(:);
    K_soil_root_profile = K_soil_root_profile(:);

    if length(P_root_dist) ~= length(psi_soil_profile)
        error('prognosticPlantWater:RootSoilMismatch', ...
              'P_root_dist 与 psi_soil_profile 长度必须相等');
    end
    % 归一化保险（防止 P 不和为 1 引起通量偏差）
    P_root_dist_norm = P_root_dist / max(sum(P_root_dist), 1e-12);

    K_sr_eff = sum(P_root_dist_norm .* K_soil_root_profile);    % 根加权土壤导度
    psi_soil_eff = sum(P_root_dist_norm .* psi_soil_profile);   % 根加权土壤水势

    q_soil_to_root = K_sr_eff * (psi_soil_eff - psi_old(1));    % [m/s]

    % ============================================================
    % Step 4: 木质部通量 (Darcy + 重力) / Xylem fluxes
    %   注意：q > 0 表示从下往上流（与重力反向）
    % ============================================================
    q_root_to_stem = -K_rs * ((psi_old(2) - psi_old(1)) / dz_rs + 1.0);
    q_stem_to_leaf = -K_sl * ((psi_old(3) - psi_old(2)) / dz_sl + 1.0);

    % ============================================================
    % Step 5: 实际蒸腾 (水力截断) / Actual transpiration with hydraulic cutoff
    %
    % [FIX 2026-05-28] The leaf has its OWN capacitance: it transpires from
    % stored leaf water first, drawing ψ_leaf down, which THEN drives stem→leaf
    % supply. Limiting T_actual by the instantaneous q_stem_to_leaf (which is
    % even negative at full saturation due to the gravity +1 term) wrongly pins
    % T_actual=0 and flattens ψ_leaf — the exact "flat ψ_leaf" failure mode.
    % Correct cutoff: transpire at potential rate until the leaf water content
    % would fall below the wilting limit (ψ_leaf ≤ ψ_min); only then throttle.
    if isfield(plant_params,'psi_leaf_min') && ~isempty(plant_params.psi_leaf_min)
        psi_leaf_min = plant_params.psi_leaf_min;   % [m]
    else
        psi_leaf_min = -3.0/9.81e-3;                % [m] ≈ -3 MPa wilting limit
    end
    % Max water the leaf can release this step without dropping below ψ_min:
    %   W at ψ_min = Wcap·(1 + a·ψ_min);  available = (W_old − W_min)/Δt
    W_leaf_min = Wcap(3) * (1 + a_ret(3) * psi_leaf_min);
    avail_leaf = max(0, (W_old(3) - W_leaf_min)) / Delt_t;          % [m/s]
    % Supply from stem is also available (positive part), plus leaf storage:
    q_supply_max = max(q_stem_to_leaf, 0) + avail_leaf;             % [m/s]
    T_actual = min(T_potential, q_supply_max);

    % ============================================================
    % Step 6: 隐式雅可比 3×3 / Implicit Jacobian 3×3
    %   d psi / d W = 1 / (a * W_cap)
    %   方程：dW/dt - flux_jacobian * dW = R_explicit
    % ============================================================
    dpsi_dW = 1.0 ./ (a_ret .* Wcap);   % [3×1, m^-1]

    % --- 隐式矩阵 A (3×3) ---
    %   A(1,1): root 方程对 W_root 的贡献
    A = zeros(3, 3);

    % 行 1: 根方程 dW_r/dt = q_sr - q_rs
    A(1,1) = 1/Delt_t + (K_sr_eff + K_rs/dz_rs) * dpsi_dW(1);
    A(1,2) = -K_rs/dz_rs * dpsi_dW(2);
    A(1,3) = 0;

    % 行 2: 茎方程 dW_s/dt = q_rs - q_sl
    A(2,1) = -K_rs/dz_rs * dpsi_dW(1);
    A(2,2) = 1/Delt_t + (K_rs/dz_rs + K_sl/dz_sl) * dpsi_dW(2);
    A(2,3) = -K_sl/dz_sl * dpsi_dW(3);

    % 行 3: 叶方程 dW_l/dt = q_sl - T_actual
    %   T_actual 对 psi_leaf 的依赖：水力截断时 dT/dpsi ≈ K_leaf
    %   未截断时 dT/dpsi ≈ 0（T_potential 与 psi 无关）
    is_cutoff = (T_actual < T_potential - 1e-12);
    if is_cutoff
        dT_dpsi_leaf = K_leaf;
    else
        dT_dpsi_leaf = 0;
    end
    A(3,1) = 0;
    A(3,2) = -K_sl/dz_sl * dpsi_dW(2);
    A(3,3) = 1/Delt_t + (K_sl/dz_sl + dT_dpsi_leaf) * dpsi_dW(3);

    % --- 显式残差 R / Explicit RHS ---
    R = [q_soil_to_root - q_root_to_stem;     % 根方程
         q_root_to_stem - q_stem_to_leaf;     % 茎方程
         q_stem_to_leaf - T_actual];          % 叶方程

    % ============================================================
    % Step 7: 求解 dW，加防奇异保护 / Solve dW with singularity guard
    % ============================================================
    if rcond(A) < 1e-12
        % 矩阵接近奇异：退化为显式 Forward Euler（不稳定但可继续）
        dW = R * Delt_t;
        warning('prognosticPlantWater:NearSingular', ...
                '隐式矩阵接近奇异，回退到显式步');
    else
        dW = A \ R;
    end

    % ============================================================
    % Step 8: 物理钳位 / Physical clamping
    % ============================================================
    W_new = W_old + dW;
    W_new = max(W_new, 0);                 % 不能为负
    W_new = min(W_new, Wcap .* 1.001);     % 不能超容量（0.1% 容差）

    % ============================================================
    % Step 9: 更新水势 / Update potentials
    % ============================================================
    psi_new = vegetation.hydraulics.linearRetention(W_new, Wcap, a_ret);

    % ============================================================
    % Step 10: 通量结构体 / Output flux struct
    %   重新算一次（用 psi_new）以确保物理一致
    % ============================================================
    q_fluxes = struct( ...
        'q_soil_to_root', q_soil_to_root, ...
        'q_root_to_stem', q_root_to_stem, ...
        'q_stem_to_leaf', q_stem_to_leaf, ...
        'T_actual',       T_actual);

    % ============================================================
    % Step 11: 守恒诊断 / Conservation diagnostic
    %   理论：sum(W_new - W_old) = (q_sr - T_actual) * dt
    %   误差表示 IMEX 离散偏差
    % ============================================================
    cum_balance = (q_soil_to_root - T_actual) * Delt_t;
    actual_dW   = sum(W_new - W_old);
    water_balance_residual = actual_dW - cum_balance;

    % ============================================================
    % Step 12: 水力胁迫因子 / Hydraulic stress factor
    %   stress = K_leaf / K_sat ∈ [0, 1]
    %   1 = 健康，0 = 完全胁迫
    % ============================================================
    stress_factor = K_leaf / plant_params.K_sat;
    stress_factor = max(0, min(1, stress_factor));

    info = struct( ...
        'K_root',                  K_root, ...
        'K_stem',                  K_stem, ...
        'K_leaf',                  K_leaf, ...
        'K_rs',                    K_rs, ...
        'K_sl',                    K_sl, ...
        'q_supply_max',            q_supply_max, ...
        'water_balance_residual',  water_balance_residual, ...
        'stress_factor',           stress_factor, ...
        'is_hydraulic_cutoff',     is_cutoff);
end
