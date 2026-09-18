function h = van_genuchten_theta_to_h(theta, theta_s, theta_r, alpha, n)
    %VAN_GENUCHTEN_THETA_TO_H van Genuchten: 含水量 → 基质势
    %
    % 语法：
    %   h = van_genuchten_theta_to_h(theta, theta_s, theta_r, alpha, n)
    %
    % 输入参数：
    %   theta   - 含水量 (cm³/cm³)，可为标量或向量
    %   theta_s - 饱和含水量 (cm³/cm³)
    %   theta_r - 残余含水量 (cm³/cm³)
    %   alpha   - van Genuchten α参数 (1/cm)
    %   n       - van Genuchten n参数
    %
    % 输出参数：
    %   h       - 基质势 (cm)
    %
    % 功能说明：
    %   使用van Genuchten方程计算含水量对应的基质势
    %
    % van Genuchten方程（正向）：
    %   Se = (θ - θ_r) / (θ_s - θ_r)              [有效饱和度]
    %   h = -1/α * (Se^(-1/m) - 1)^(1/n)          [基质势]
    %   其中 m = 1 - 1/n
    %
    % 边界处理：
    %   θ ≥ θ_s  →  h = -1e-6 (接近0，饱和)
    %   θ ≤ θ_r  →  h = -1e7 (很干旱)
    %   θ_r < θ < θ_s  →  使用van Genuchten方程
    %
    % 参考：
    %   van Genuchten, M.T. (1980). A Closed-form Equation for Predicting the
    %   Hydraulic Conductivity of Unsaturated Soils. Soil Science Society of
    %   America Journal, 44: 892-898.
    %
    % 示例：
    %   h = equations.van_genuchten_theta_to_h(0.30, 0.45, 0.07, 0.01, 1.5);

    % 计算 m = 1 - 1/n
    m = 1 - 1/n;

    % 计算有效饱和度
    Se = (theta - theta_r) / (theta_s - theta_r);
    Se = max(min(Se, 1.0), 1e-12);

    % van Genuchten方程：h = -1/α * (Se^(-1/m) - 1)^(1/n)
    if Se >= 1.0
        h = -1e-6;      % 饱和
    elseif Se <= 1e-12
        h = -1e7;       % 极干旱
    else
        term = Se^(-1/m) - 1;
        h = -(1/alpha) * (term^(1/n));
    end

end
