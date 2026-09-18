function dtheta = calculateDTheta(subRoutine, heat_term, theta_s, theta_r, theta_m, gamma_hh, theta_ll, theta_l, theta_uu, theta_u, hh, h, hh_frez, h_frez, phi_s, lamda, alpha, n, m)
% calculateDTheta - 按子路径计算含水量对水头的导数
%                   Compute derivative of water content w.r.t. matric head by subroutine
%
% 函数定义 / Function Definition:
%   dtheta = calculateDTheta(subRoutine, heat_term, theta_s, theta_r, theta_m, gamma_hh, ...
%                            theta_ll, theta_l, theta_uu, theta_u, hh, h, hh_frez, h_frez, ...
%                            phi_s, lamda, alpha, n, m)
%
% 子路径索引 / Subroutine index:
%   0 - lFreezing=1 解析路径（含 hd 截断）/ lFreezing=1 analytical (with hd clamp)
%   1 - VG 解析导数，使用 theta_m / VG analytical derivative using theta_m
%   2 - theta_ll 割线；极小分母时用解析切线 / Secant with analytical fallback
%   3 - theta_uu 割线；极小分母时用解析切线 / Secant with analytical fallback
%   4 - Clapp-Hornberger 解析导数 / Clapp-Hornberger analytical derivative
%   5 - VG 解析导数，使用 theta_s（见 issue 181 item 7）/ VG analytical using theta_s
%
% 修改者 / Modified by:   朱寄子星 (Zhu Jizixing)
% 修改日期 / Modified date: 2026-04-11

SoilConstants = io.getSoilConstants();

switch subRoutine
    case 0
        % lFreezing=1 解析路径 / lFreezing=1 analytical path
        % heat_term = hh + hh_frez (for DTheta_UUh) or hh (for DTheta_LLh)
        A = (-theta_r) / (abs(heat_term) * log(abs(SoilConstants.hd / SoilConstants.hm))) ...
            * (1 - (1 + abs(alpha * heat_term)^n)^(-m));
        B = alpha * n * m * (theta_m - gamma_hh * theta_r);
        C = (1 + abs(alpha * heat_term)^n)^(-m - 1);
        D = abs(alpha * heat_term)^(n - 1);
        dtheta = A - B * C * D;
        if (hh + hh_frez) <= SoilConstants.hd
            dtheta = 0;
        end

    case 1
        % heat_term = hh or  hh + hh_frez
        A = (theta_m - theta_r) * alpha * n;
        B =  abs(alpha * heat_term)^(n - 1) * (-m);
        C = (1 + abs(alpha * heat_term)^n)^(-m - 1);
        dtheta = A * B * C;

    case 2
        % theta_ll secant when the increment is resolved; analytical tangent
        % when the denominator is effectively zero.
        A = theta_ll - theta_l;
        B = hh + hh_frez - h - h_frez;
        if abs(B) > derivativeHeadEps()
            dtheta = A / B;
        else
            dtheta = vgDerivative(theta_m, theta_r, alpha, n, m, heat_term);
        end

    case 3
        % theta_uu secant when the increment is resolved; analytical tangent
        % when the denominator is effectively zero.
        A = theta_uu - theta_u;
        B = hh - h;
        if abs(B) > derivativeHeadEps()
            dtheta = A / B;
        else
            dtheta = vgDerivative(theta_m, theta_r, alpha, n, m, heat_term);
        end

    case 4
        % Clapp-Hornberger 解析导数 / Clapp-Hornberger analytical derivative
        dtheta = theta_s / phi_s * (hh / phi_s)^(-1 * lamda - 1);

    case 5
        % VG analytical derivative using theta_s, following the legacy STEMMUS negative-sign convention.
        A = (theta_s - theta_r) * alpha * n;
        B =  abs(alpha * heat_term)^(n - 1) * (-m);
        C = (1 + abs(alpha * heat_term)^n)^(-m - 1);
        dtheta = A * B * C;
end
end

function eps_h = derivativeHeadEps()
% Below this head increment, use the analytical tangent instead of a secant.
eps_h = 1e-6;
end

function dtheta = vgDerivative(theta_m, theta_r, alpha, n, m, heat_term)
A = (theta_m - theta_r) * alpha * n;
B =  abs(alpha * heat_term)^(n - 1) * (-m);
C = (1 + abs(alpha * heat_term)^n)^(-m - 1);
dtheta = A * B * C;
end
