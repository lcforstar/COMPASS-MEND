function psi = linearRetention(W, W_capacity, a)
% linearRetention  植物组织线性持水曲线 / Linear plant tissue retention curve
%
% 公式 / Formula (ClimaLand 2026 Eq. 35):
%   psi(W) = (1/a) * (W/W_cap - 1)
%
% 性质 / Properties:
%   - psi(W = W_cap) = 0      (零吸力即饱和)
%   - psi(W = 0)     = -1/a   (最大吸力)
%   - dpsi/dW        = 1 / (a * W_cap)
%
% 输入 / Inputs:
%   W          - 当前体积水含量 / Current water content [m, 等效水柱]
%   W_capacity - 储水容量 / Storage capacity [m]
%   a          - 持水曲线斜率 / Retention slope [m^-1]
%
% 输出 / Outputs:
%   psi        - 水势 / Water potential [m] (负值表示吸力)
%
% 物理参考 / Reference:
%   - Deck et al. (2026) ClimaLand JAMES, Eq. 35
%   - Christoffersen et al. (2016) GMD, 9, 4227–4255
%
% Author: Zhu Jizixing
% Date:   2026-05-01

    % --- 输入检查 / Input checks ---
    if any(W_capacity(:) <= 0)
        error('linearRetention:InvalidCapacity', 'W_capacity 必须 > 0');
    end
    if any(a(:) <= 0)
        error('linearRetention:InvalidA', 'a 必须 > 0');
    end

    % --- 线性持水曲线 / Linear retention curve ---
    psi = (W ./ W_capacity - 1) ./ a;
end
