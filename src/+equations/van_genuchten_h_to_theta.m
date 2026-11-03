function theta = van_genuchten_h_to_theta(h, theta_s, theta_r, alpha, n, T, enable_thermal)
%%VAN_GENUCHTEN_H_TO_THETAvarargin)
    %VAN_GENUCHTEN_H_TO_THETA van Genuchten: 基质势 → 含水量（可选温度修正）
    %
    % 语法：
    %   theta = van_genuchten_h_to_theta(h, theta_s, theta_r, alpha, n)
    %   theta = van_genuchten_h_to_theta(h, theta_s, theta_r, alpha, n, T)
    %   theta = van_genuchten_h_to_theta(h, theta_s, theta_r, alpha, n, T, enable_thermal)
    %
    % 输入参数：
    %   h               - 基质势 (cm)，可为标量或向量
    %   theta_s         - 饱和含水量 (cm³/cm³)
    %   theta_r         - 残余含水量 (cm³/cm³)
    %   alpha           - van Genuchten α参数 (1/cm)
    %   n               - van Genuchten n参数
    %   T               - 温度 (°C)，可选，默认为20°C
    %   enable_thermal  - 是否启用温度修正 (logical)，可选，默认为false
    %
    % 输出参数：
    %   theta           - 含水量 (cm³/cm³)
    %
    % 功能说明：
    %   使用van Genuchten方程计算基质势对应的含水量
    %   可选的温度修正考虑表面张力随温度的变化
    %
    % van Genuchten方程（反向）：
    %   Se = [1 + (α|h|)^n]^(-m)              [有效饱和度]
    %   θ = θ_r + (θ_s - θ_r) * Se            [含水量]
    %   其中 m = 1 - 1/n
    %
    % 温度修正：
    %   表面张力随温度变化，影响基质势
    %   修正因子 cor = exp(-cor_h * (T - T_ref))
    %   h_eff = h * cor
    %   然后用 h_eff 代替 h 计算θ
    %
    %   Parameter: Cw is the shared retention/flow temperature slope [K-1].
    %        T_ref = 20 (°C)
    %
    % 边界处理：
    %   h ≥ -1e-6   →  θ = θ_s (饱和)
    %   h ≤ -1e7    →  θ = θ_r (极干旱)
    %   -1e7 < h < -1e-6  →  使用van Genuchten方程
    %
    % 参考：
    %   van Genuchten, M.T. (1980). A Closed-form Equation for Predicting the
    %   Hydraulic Conductivity of Unsaturated Soils. Soil Science Society of
    %   America Journal, 44: 892-898.
    %
    % 示例：
    %   theta = equations.van_genuchten_h_to_theta(-100, 0.45, 0.07, 0.01, 1.5);
    %   theta = equations.van_genuchten_h_to_theta(-100, 0.45, 0.07, 0.01, 1.5, 25, true);
    %
    %   Editor: Jizixing Zhu; Date: 2026-07-13.
    %   Change: 2026-07-13 Jizixing Zhu -- restore the documented optional
    %           temperature and thermal-switch defaults.

    % 默认参数
    if nargin < 6 || isempty(T)
        T = equations.getPhysicalConstant('Tr');
    end
    if nargin < 7 || isempty(enable_thermal)
        enable_thermal = false;
    end

    % 确保输入维度一致
    h = h(:);
    theta = zeros(size(h));

    % 计算 m = 1 - 1/n
    m = 1 - 1/n;

    % =========================================================================
    % 温度修正（如果启用）
    % =========================================================================
    if enable_thermal
        % 表面张力修正因子
        % 根据物理模型：表面张力随温度升高而降低
        % σ(T) / σ(T_ref) ≈ exp(-cor_h * (T - T_ref))
        cor_h = equations.getPhysicalConstant('Cw');  % [K-1] Shared retention/flow slope

        % 确保T为列向量，维度与h一致
       
        T_vec = T(:);
        

        % 计算温度修正因子
        cor = exp(-cor_h * (T_vec - equations.getPhysicalConstant('Tr')));

        % 应用修正
        h_eff = h .* cor;
    else
        h_eff = h;
    end

    % =========================================================================
    % van Genuchten计算
    % =========================================================================

    % 处理饱和情况：h ≥ -1e-6
    sat_idx = h_eff >= -1e-6;
    theta(sat_idx) = theta_s;

    % 处理极端干旱情况：h ≤ -1e7
    dry_idx = h_eff <= -1e7;
    theta(dry_idx) = theta_r;

    % 处理中间情况：-1e7 < h < -1e-6
    normal_idx = ~sat_idx & ~dry_idx;

    if any(normal_idx)
        h_normal = h_eff(normal_idx);

        % van Genuchten方程：Se = [1 + (α|h|)^n]^(-m)
        abs_h = abs(h_normal);

        % 防止数值溢出
        alpha_h = alpha * abs_h;
        alpha_h = min(alpha_h, 1e6);  % 限制最大值

        term = 1 + alpha_h.^n;
        Se = term.^(-m);

        % θ = θ_r + (θ_s - θ_r) * Se
        theta(normal_idx) = theta_r + (theta_s - theta_r) * Se;
    end

    % =========================================================================
    % 最终检查与约束
    % =========================================================================

    % 确保含水量在物理范围内
    theta = max(min(theta, theta_s), theta_r);

end
