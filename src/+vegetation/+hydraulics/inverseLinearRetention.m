function W = inverseLinearRetention(psi, W_capacity, a)
% inverseLinearRetention  线性持水曲线反函数 / Inverse linear retention curve
%
% 由 psi(W) = (W/W_cap - 1)/a 反解 W
%   W = W_cap * (psi * a + 1)
%
% 用于：从初始水势计算初始储水量 / For computing initial W from initial psi.
%
% 输入 / Inputs:
%   psi        - 水势 / Water potential [m]
%   W_capacity - 储水容量 / Storage capacity [m]
%   a          - 持水曲线斜率 / Retention slope [m^-1]
%
% 输出 / Outputs:
%   W          - 体积水含量 / Water content [m]
%
% Author: Zhu Jizixing
% Date:   2026-05-01

    if any(W_capacity(:) <= 0)
        error('inverseLinearRetention:InvalidCapacity', 'W_capacity 必须 > 0');
    end
    if any(a(:) <= 0)
        error('inverseLinearRetention:InvalidA', 'a 必须 > 0');
    end

    W = W_capacity .* (psi .* a + 1);

    % --- 物理钳位 / Physical clamping ---
    W = max(W, 0);                      % 不能为负
    W = min(W, W_capacity .* 1.001);    % 不能超容量（留 0.1% 容差）
end
