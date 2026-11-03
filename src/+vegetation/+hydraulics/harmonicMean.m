function K_eff = harmonicMean(K1, K2)
% harmonicMean  导水率调和平均 / Harmonic mean of conductivities
%
% 公式 / Formula:
%   K_eff = 2 * K1 * K2 / (K1 + K2)
%
% 物理意义 / Physical meaning:
%   两个串联节点的有效导水率，跟随更小的那个 (low K dominates).
%   ClimaLand water_flux (PlantHydraulics.jl line 274) 使用此式。
%
% 输入 / Inputs: K1, K2 - 两节点导水率 (m/s)
% 输出 / Outputs: K_eff
%
% Author: Zhu Jizixing
% Date:   2026-05-01

    K1_safe = max(K1, 1e-30);
    K2_safe = max(K2, 1e-30);
    K_eff = 2 .* K1_safe .* K2_safe ./ (K1_safe + K2_safe);
end
