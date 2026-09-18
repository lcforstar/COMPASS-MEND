function Theta_m = calculateTheta_m(Gamma_hh, VanGenuchten, i, POR_i)
% calculateTheta_m - 计算毛管含水量上限（第 i 层，SFCC=0 路径）
%                    Compute capillary water content upper bound for layer i
%
% 函数定义 / Function Definition:
%   Theta_m = calculateTheta_m(Gamma_hh, VanGenuchten, i, POR_i)
%
% 输入 / Inputs:
%   Gamma_hh    - 液态水有效性系数（标量）/ Liquid water availability factor (scalar)
%   VanGenuchten - VG 参数结构体（各字段为 NL×1）/ VG parameter struct (fields NL×1)
%   i           - 土层索引 / Layer index
%   POR_i       - 第 i 层孔隙度 / Porosity of layer i
%
% 输出 / Outputs:
%   Theta_m     - 毛管含水量上限（标量）/ Capillary water content upper bound (scalar)
%
% 修改者 / Modified by:   朱寄子星 (Zhu Jizixing)
% 修改日期 / Modified date: 2026-04-10

    theta_s = VanGenuchten.Theta_s(i);
    theta_r = VanGenuchten.Theta_r(i);
    alpha   = VanGenuchten.Alpha(i);
    n       = VanGenuchten.n(i);
    m       = VanGenuchten.m(i);

    % hh=-1 处的上界含水量 / Upper bound water content at hh = -1
    Theta_m = Gamma_hh * theta_r + (theta_s - Gamma_hh * theta_r) * (1 + abs(alpha * (-1))^n)^(-m);

    % 孔隙度和饱和含水量截断 / Clamp to [theta_s, POR]
    if Theta_m >= POR_i
        Theta_m = POR_i;
    elseif Theta_m <= theta_s
        Theta_m = theta_s;
    end
end
