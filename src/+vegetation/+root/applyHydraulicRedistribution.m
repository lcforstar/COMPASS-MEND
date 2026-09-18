function [Theta_LL_new, hr_diag] = applyHydraulicRedistribution( ...
    Theta_LL, hh, gsr_node, psir, DeltZ_node, Delt_t, VanGenuchten, NL)
% applyHydraulicRedistribution - CLM5 风格水力再分配 (Step 2)
%
% Richards 解完 → 调用本模块更新 θ，让根系做层间真实搬运。
% 总量自动守恒（HR 只在层间搬运，整体净流 = SPAC actual_total 已被 Richards 吸收）。
%
% 物理:
%   HR_flux_j = gsr_j × (ψ_soil_j_new − ψ_root)     [m/s]
%     ψ_soil_j_new 来自 Richards 解后的新 hh
%     正 → 土→根 (干层吸水), 负 → 根→土 (湿层 HR 回灌)
%   Δθ_j     = −HR_flux_j × Δt × 100 / Δz_j         [-]
%     × 100：flux [m/s] × 100 cm/m / Δz [cm] = [1/s]，再 × Δt 得 [-]
%   边界保护: θ ∈ [θr+ε, θs−ε]
%
% Inputs:
%   Theta_LL    - NLx2  Richards 后液相含水量 [-]
%   hh          - NNx1  Richards 后水头 [cmH2O]
%   gsr_node    - NNx1  土-根节点导度 [m/s/MPa]
%   psir        - 标量  根水势 [MPa]
%   DeltZ_node  - NNx1  节点厚度 [cm]
%   Delt_t      - 标量  时间步 [s]
%   VanGenuchten- 结构  含 Theta_r, Theta_s
%   NL          - 标量  层数
%
% Outputs:
%   Theta_LL_new - NLx2  HR 后含水量 [-]
%   hr_diag      - 结构  诊断

hr_diag = struct( ...
    'total_HR_pos_m_s', 0, 'total_HR_neg_m_s', 0, ...
    'net_residual_m_s', 0, 'n_clipped', 0, ...
    'mass_imbalance_cm', 0);

% 防御：输入不合法直接返回原值
if isempty(gsr_node) || ~isfinite(psir) || all(gsr_node <= 0)
    Theta_LL_new = Theta_LL;
    return;
end

% ── 1) Richards 后新水势 (cmH2O → MPa) ──
hhMPA_new = hh(:) * 9.80665e-5;     % NNx1 [MPa]

% ── 2) 节点级 HR 通量 (含正负) ──
hr_flux = gsr_node(:) .* (hhMPA_new - psir);    % NNx1 [m/s]

% ── 3) 节点 → 层 θ 变化 ──
DeltZ_node = DeltZ_node(:);
delta_theta_node = -hr_flux .* (Delt_t * 100) ./ DeltZ_node;   % NNx1 [-]
% × 100 是 m/s → cm/s 转换

% 节点端点 → 层端点（NL×2）
delta_theta_layer = zeros(NL, 2);
for ML = 1:NL
    delta_theta_layer(ML, 1) = delta_theta_node(ML);
    delta_theta_layer(ML, 2) = delta_theta_node(ML + 1);
end

% ── 4) 应用 + 边界保护 ──
Theta_r = VanGenuchten.Theta_r(:);
Theta_s = VanGenuchten.Theta_s(:);
Theta_r_layer = repmat(Theta_r, 1, 2);
Theta_s_layer = repmat(Theta_s, 1, 2);
eps_theta = 1e-6;

Theta_LL_proposed = Theta_LL + delta_theta_layer;
Theta_LL_new      = max(Theta_r_layer + eps_theta, min(Theta_s_layer - eps_theta, Theta_LL_proposed));

% ── 5) 诊断 ──
clip_loss_cm = sum(sum( (Theta_LL_proposed - Theta_LL_new) .* (repmat(DeltZ_node(1:NL), 1, 2)) )) / 2;
hr_diag.mass_imbalance_cm  = clip_loss_cm;        % 正=系统失水（限幅切了）
hr_diag.total_HR_pos_m_s   = sum(max(hr_flux, 0));
hr_diag.total_HR_neg_m_s   = sum(min(hr_flux, 0));
hr_diag.net_residual_m_s   = sum(hr_flux);
hr_diag.n_clipped          = nnz(abs(Theta_LL_proposed - Theta_LL_new) > 1e-9);
end
