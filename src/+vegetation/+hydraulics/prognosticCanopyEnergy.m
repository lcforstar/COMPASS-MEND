function [T_new, dT, info] = prognosticCanopyEnergy( ...
    T_old, T_iter, Rn, H, LE, S_re_e, ...
    rho_c_h, A_canopy, ...
    rho_a, cp, Lv, r_ah, r_ae, ...
    qsat_curr, T_air, ...
    epsilon_canopy, sigmaSB, R_v, ...
    Delt_t)
% prognosticCanopyEnergy  瞬态冠层能量平衡（IMEX 半隐式） / Prognostic canopy energy
%
% 求解方程 / Equation (ClimaLand 2026 Eq. 39):
%   (rho * c * h) * A_canopy * dT/dt = R_n - LE - H + S_re_e
%
% 时间离散 / Time discretization (IMEX semi-implicit):
%   把 (H, LE, R_n) 对 T 的强非线性部分隐式化，结果是单步标量代数方程：
%
%       A_imp * (T^{n+1} - T^n) = R_explicit
%
%   A_imp = (rho_c_h * A_canopy)/dt
%         + rho_a * cp / r_ah                         (显热雅可比)
%         + (rho_a * Lv / r_ae) * (Lv * qsat) / (R_v * T^2)   (Clausius-Clapeyron)
%         + 8 * eps * sigma * T^3                      (长波辐射)
%
%   R_explicit = R_n - H - LE + S_re_e
%
% 数值优势 / Numerical advantage:
%   - 单步线性代数，零成本
%   - IMEX 无条件稳定（任意 dt）
%   - 替换原 Picard "T += Wc * EBer / (...)"，结构相同但物理正确
%
% 输入 / Inputs:
%   T_old        - 上一时步温度 / Previous step temperature [K]
%   T_iter       - 当前 Picard 迭代温度 / Current iteration temperature [K]
%   Rn, H, LE    - 净辐射、显热、潜热 / Net rad, sensible heat, latent heat [W/m^2]
%   S_re_e       - 根吸水带走的能量 / Root-extracted energy flux [W/m^2]
%   rho_c_h      - 冠层单位面积热容 / Canopy areal heat capacity [J/m^2/K]
%   A_canopy     - 冠层总面积指数 = LAI + SAI / Total area index [-]
%   rho_a, cp    - 空气密度、定压比热 / Air density, cp [kg/m^3, J/kg/K]
%   Lv           - 汽化潜热 / Latent heat of vaporization [J/kg]
%   r_ah, r_ae   - 显热/水汽阻力 / Sensible/vapor resistance [s/m]
%   qsat_curr    - 当前 T_iter 下饱和比湿 / qsat(T_iter) [kg/kg]
%   T_air        - 大气温度 / Air temperature [K]
%   epsilon_canopy, sigmaSB - 发射率、Stefan-Boltzmann / Emissivity & SB const
%   R_v          - 水汽气体常数 / Water vapor gas constant (461.5 J/kg/K)
%   Delt_t       - 时步 / Time step [s]
%
% 输出 / Outputs:
%   T_new        - 更新后温度 / Updated temperature [K]
%   dT           - 本步增量 / Increment this step [K]
%   info         - 诊断结构体 / Diagnostic struct
%     .A_implicit       隐式系数
%     .R_explicit       显式残差 [W/m^2]
%     .EB_residual      能量平衡残差（用于守恒诊断）
%     .tau_canopy       热惯性时标 [s]
%
% 物理参考 / Reference:
%   - Deck et al. (2026) ClimaLand JAMES, §2.3.2 Eq. 39
%   - Bonan G. (2019) Climate Change and Terrestrial Ecosystem Modeling, Ch. 17
%   - Lawrence et al. (2019) CLM5 Tech Note, §5.2
%
% Author: Zhu Jizixing
% Date:   2026-05-01

    % --- 输入合理性检查 / Input sanity ---
    if rho_c_h <= 0
        error('prognosticCanopyEnergy:InvalidHeatCapacity', ...
              'rho_c_h 必须 > 0');
    end
    if A_canopy < 0
        A_canopy = 0;
    end
    if Delt_t <= 0
        error('prognosticCanopyEnergy:InvalidDt', 'Delt_t 必须 > 0');
    end

    % --- 防止稀疏冠层（A_canopy → 0）时 A_imp 退化为 0 ---
    A_canopy_eff = max(A_canopy, 1e-3);

    % --- 关键阻力的安全下界 / Safe floors on resistances ---
    r_ah_safe = max(r_ah, 1.0);
    r_ae_safe = max(r_ae, 1.0);

    % --- IMEX 隐式系数 A_imp [W/m^2/K] ---
    %   各项物理对应 ClimaLand Eq. 39 的雅可比展开
    term_storage = (rho_c_h * A_canopy_eff) / Delt_t;
    term_H       = rho_a * cp / r_ah_safe;
    term_LE      = (rho_a * Lv / r_ae_safe) * (Lv * qsat_curr / (R_v * T_iter^2));
    term_LWnet   = 8.0 * epsilon_canopy * sigmaSB * T_iter^3;

    A_implicit = term_storage + term_H + term_LE + term_LWnet;

    % --- 显式残差 R_exp [W/m^2] ---
    R_explicit = Rn - H - LE + S_re_e;

    % --- IMEX 单步更新 / Single-step IMEX update ---
    dT_raw = R_explicit / max(A_implicit, 1e-6);

    % --- 物理钳位：单步增量限制 ±20 K（防止数值发散）---
    dT_max = 20.0;
    dT = max(min(dT_raw, dT_max), -dT_max);

    % --- 新温度 ---
    T_new = T_old + dT;

    % --- 物理范围钳位：合理冠层温度 [220 K, 330 K] = [-53°C, 57°C] ---
    T_new = max(min(T_new, 330.0), 220.0);

    % --- 守恒诊断 / Conservation diagnostic ---
    %   能量平衡残差（瞬态情况下应等于储能变化率）
    EB_residual = R_explicit - (rho_c_h * A_canopy_eff) * dT / Delt_t;

    % --- 热惯性时标 / Thermal inertia timescale ---
    G_h_eff = term_H + term_LE + term_LWnet;
    if G_h_eff > 0
        tau_canopy = (rho_c_h * A_canopy_eff) / G_h_eff;
    else
        tau_canopy = Inf;
    end

    info = struct( ...
        'A_implicit',  A_implicit, ...
        'R_explicit',  R_explicit, ...
        'EB_residual', EB_residual, ...
        'tau_canopy',  tau_canopy);
end
