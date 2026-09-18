function K = weibullConductance(psi, K_sat, psi_63, c)
% weibullConductance  Weibull 木质部导度模型 / Weibull xylem conductivity model
%
% 公式 / Formula (ClimaLand 2026 Eq. 34):
%   K(psi) = K_sat * exp[-(|psi| / |psi_63|)^c]    if psi < 0 (吸力)
%   K(psi) = K_sat                                  if psi >= 0 (饱和)
%
% 输入 / Inputs:
%   psi    - 木质部水势 / Xylem water potential [m] (负值代表吸力)
%   K_sat  - 饱和导水率  / Saturated conductivity [m s^-1]
%   psi_63 - 63% 损失水势 / Potential at 63% conductivity loss [m] (负值)
%   c      - Weibull 形状参数 / Weibull shape parameter [-] (典型 2-6)
%
% 输出 / Outputs:
%   K      - 当前导水率 / Current conductivity [m s^-1]
%
% 物理参考 / Reference:
%   - Deck et al. (2026) ClimaLand JAMES 18, e2025MS005118, Eq. 34
%   - Sperry & Love (2015) New Phytol., 207, 14–27
%   - Liu et al. (2020) Glob. Change Biol., 26, 3949–3962
%
% Author: Zhu Jizixing
% Date:   2026-05-01

    % --- 输入合理性检查 / Input sanity checks ---
    if K_sat <= 0
        error('weibullConductance:InvalidKsat', 'K_sat 必须 > 0 (m/s)');
    end
    if psi_63 >= 0
        error('weibullConductance:InvalidPsi63', 'psi_63 应 < 0 (m)');
    end
    if c <= 0
        error('weibullConductance:InvalidC', 'c 必须 > 0');
    end

    % --- 向量化计算 / Vectorized evaluation ---
    K = K_sat * ones(size(psi));
    mask = psi < 0;
    if any(mask(:))
        ratio = abs(psi(mask)) / abs(psi_63);
        K(mask) = K_sat * exp(-ratio.^c);
    end

    % --- 数值下界（防止 0 除）/ Numerical floor ---
    K = max(K, 1e-30);
end
