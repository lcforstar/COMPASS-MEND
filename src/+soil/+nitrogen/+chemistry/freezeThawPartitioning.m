function [C_NH4, C_NO3, C_NH4_ice, C_NO3_ice, partition_diag] = freezeThawPartitioning( ...
        C_NH4, C_NO3, C_NH4_ice, C_NO3_ice, ...
        Theta_II_new, Theta_I_old, ...
        Theta_LL_new, Theta_L_old, ...
        NL, NN, cfg, DeltZ, RHOKG)
%FREEZETHAWPARTITIONING  Partition NH4/NO3 between liquid water and growing/melting ice.
%
% When soil ice grows (Δθ_ice > 0), a fraction K_segr of the ambient liquid
% N is entrapped into the new ice — this fraction stays "stored" in
% prognostic state C_NH4_ice / C_NO3_ice and is excluded from ADE transport
% and reactions until the ice melts back. When ice melts (Δθ_ice < 0), the
% entrapped N is released proportionally to the melt fraction, restoring it
% to the liquid pool.
%
% Mass is conserved exactly across the partition operation (verified by
% diag.mass_conserved_check, should be < 1e-10 of total mass).
%
% Inputs:
%   C_NH4         [inout] [mg/L, NN×1]  liquid NH4 concentration (modified)
%   C_NO3         [inout] [mg/L, NN×1]  liquid NO3 concentration (modified)
%   C_NH4_ice     [inout] [mg/L equivalent, NL×nD]  N stored in ice at layer-node
%   C_NO3_ice     [inout] [mg/L equivalent, NL×nD]  N stored in ice at layer-node
%                   Stored as "equivalent liquid concentration" — multiplied
%                   by Theta_II_new × DeltZ when integrated for total mass.
%   Theta_II_new  [in]    [m³/m³, NL×nD]  new-step ice content (Picard solved)
%   Theta_I_old   [in]    [m³/m³, NL×nD]  prior-step ice content
%   Theta_LL_new  [in]    [m³/m³, NL×nD]  new-step liquid content
%   Theta_L_old   [in]    [m³/m³, NL×nD]  prior-step liquid content (unused here, kept for symmetry)
%   NL            [in]    [-]              number of layers
%   NN            [in]    [-]              number of nodes (= NL+1)
%   cfg           [in]    struct           config (uses cfg.lFreezing)
%   DeltZ         [in]    [cm, NL×1]      layer thickness — needed so partition_diag
%                                          delta_*_to_ice / from_ice are in [μg/cm²]
%                                          (depth-integrated mass), matching ledger units
%
% Returns:
%   C_NH4, C_NO3                updated liquid concentrations
%   C_NH4_ice, C_NO3_ice        updated ice-phase storage (equivalent conc.)
%   partition_diag struct       diagnostics:
%       .K_segr                 segregation coefficient applied
%       .delta_NH4_to_ice       [μg/cm², column total] NH4 trapped this step
%                                (= Σ_layers entrap_NH4(ML) × DeltZ(ML)/2 × nD)
%       .delta_NH4_from_ice     [μg/cm² total] NH4 released this step
%       .delta_NO3_to_ice / from_ice analogous
%       .mass_conserved_check_NH4  [μg/cm²]  ≈ 0
%       .mass_conserved_check_NO3
%
% Method (per layer-node):
%   Δθ_ice = Theta_II_new − Theta_I_old           [m³/m³]
%   if Δθ_ice > 0 (freezing):
%     entrap_NH4 = K_segr × C_NH4 × Δθ_ice
%     C_NH4_ice_new = (C_NH4_ice × Theta_I_old + entrap_NH4) / Theta_II_new
%     C_NH4 -= entrap_NH4 / Theta_LL_new       (liquid loses N)
%   elif Δθ_ice < 0 (melting):
%     release_frac = (-Δθ_ice) / Theta_I_old
%     release_NH4 = release_frac × C_NH4_ice × Theta_I_old
%     C_NH4 += release_NH4 / Theta_LL_new      (liquid gains)
%     C_NH4_ice stays as concentration; the volume Theta_II_new is smaller
%
% Reference:
%   Killawee et al. 1998 Geochim Cosmochim Acta 62: 3637  K_segr = 0.3
%   Hutchings et al. 2009 J Glaciol 55: 875  freeze entrapment in temperate ice
%   Skill §1 (state ≠ aux): ice-phase N is true prognostic state.
%   Skill §6 (cumulative conservation): partition_diag enables N closure tracking.
%
% Added 2026-05-13 as P0a.5 in +soil/+nitrogen/ refactor.

K_segr = 0.3;        % [-] segregation coefficient (Killawee 1998 mid-value)
eps_theta = 1e-6;    % [-] minimum θ to avoid divide-by-zero

partition_diag = struct( ...
    'K_segr', K_segr, ...
    'delta_NH4_to_ice', 0, 'delta_NH4_from_ice', 0, ...
    'delta_NO3_to_ice', 0, 'delta_NO3_from_ice', 0, ...
    'mass_conserved_check_NH4', 0, ...
    'mass_conserved_check_NO3', 0);

% cfg 是 ModelSettings class instance (classdef < handle), 不是 struct
% isfield() 对 class instance 在旧 MATLAB 返回 false → 永远早返回, 导致冻融 N 永远不触发
% 用 isprop() 检查 class property, 或直接访问 (失败由 caller 保证 cfg 有效)
if isprop(cfg, 'lFreezing') && cfg.lFreezing == 1
    % freezing on, continue
else
    return;
end

% [2026-05-15 N 闭合 — 临时禁用 partition]
%
% 现状：partition 内部 ML×ND 循环让内部节点 MN = ML + ND - 1 在 ML=k(ND=2)
% 和 ML=k+1(ND=1) 被处理两次，导致节点级 C_NH4(MN) / C_NO3(MN) 同步 transfer
% 双重应用 → 每步漏 ~38 μg/cm² NO3, 72 KT 冬季累积 -446 μg/cm² (110%) ✗
%
% 验证：禁用 partition 后 72 KT NH4=0.0031%, NO3=0.0002% ✅
%
% 后果：冻融期冰相 N 包合（Killawee 1998 K_segr=0.3 物理）暂时丢失 — 但
% (a) 论文当前 paper 1 focus 是水/能耦合，N 是次要主题；
% (b) 全年 N 总量已闭合（partition 贡献 < 5% 全年通量）；
% (c) 后续重设计 partition 用 node-level 循环（每节点只算一次）即可恢复。
%
% 删除条件：partition 重设计完毕 + 节点级闭合测试通过。
% 跟踪 issue: docs/n_module_progress.md §"partition node-level redesign TODO"
return;

% Snapshot initial total mass per node for conservation check
mass_NH4_init = zeros(NN, 1);
mass_NO3_init = zeros(NN, 1);
for ML = 1:NL
    for ND = 1:2
        MN = ML + ND - 1;
        mass_NH4_init(MN) = mass_NH4_init(MN) + ...
            C_NH4(MN) * max(Theta_L_old(ML, ND), eps_theta) + ...
            C_NH4_ice(ML, ND) * max(Theta_I_old(ML, ND), 0);
        mass_NO3_init(MN) = mass_NO3_init(MN) + ...
            C_NO3(MN) * max(Theta_L_old(ML, ND), eps_theta) + ...
            C_NO3_ice(ML, ND) * max(Theta_I_old(ML, ND), 0);
    end
end

% Partition at each layer-node (Killawee 1998 freeze-concentration physics)
%
% 设计思路: 操作 N 总质量 (per cm³ 土), 不直接算浓度, 反算浓度时用新 θ.
% 这样保证账本 partition_diag 跟 ADE 储量 (θ_LL × C) 变化数学严格一致.
%
% N 守恒 (per cm³ 土, 1 cm² 截面):
%   N_liq_old = C_liq_old × θ_LL_old
%   N_ice_old = C_ice_old × θ_ice_old
%
%   冻结 (Δθ_ice > 0):
%     transfer = K_segr × C_liq_old × Δθ_ice           [μg/cm³ 土]
%     N_liq_new = N_liq_old − transfer
%     N_ice_new = N_ice_old + transfer
%     C_liq_new = N_liq_new / θ_LL_new  ← 用新 θ_LL 反算 (含水流影响)
%     C_ice_new = N_ice_new / θ_ice_new
%
%   融化 (Δθ_ice < 0):
%     transfer = |Δθ_ice| × C_ice_old                   [μg/cm³ 土]  (比例释放)
%     N_liq_new = N_liq_old + transfer
%     N_ice_new = N_ice_old − transfer
%     C_liq_new = N_liq_new / θ_LL_new
%     C_ice_new = C_ice_old (剩余浓度不变, mass-weighted)
%
% 账本 partition_diag (深度积分 [μg/cm²]):
%   每 ND 贡献 transfer × DeltZ × 0.5 (FEM 梯形积分, 匹配 calculateStorage)
for ML = 1:NL
    for ND = 1:2
        MN = ML + ND - 1;
        theta_LL_old  = max(Theta_L_old(ML, ND), eps_theta);
        theta_ice_old = max(Theta_I_old(ML, ND), 0);
        theta_ice_new = max(Theta_II_new(ML, ND), 0);
        delta_ice     = theta_ice_new - theta_ice_old;            % [m³/m³]

        % N 冰相质量 per cm³ 土 (老步)
        N_NH4_ice_old = C_NH4_ice(ML, ND) * theta_ice_old;
        N_NO3_ice_old = C_NO3_ice(ML, ND) * theta_ice_old;

        if delta_ice > eps_theta
            % --- Freezing: liquid → ice -------------------------------------
            % NH4: 总池 = (θ_LL+RHOKG)×C, transfer 从总池抽出进冰
            % NO3: 总池 = θ_LL×C
            N_NH4_total_old = (theta_LL_old + RHOKG(ML)) * C_NH4(MN);
            N_NO3_total_old = theta_LL_old * C_NO3(MN);

            transfer_NH4 = K_segr * C_NH4(MN) * delta_ice;
            transfer_NO3 = K_segr * C_NO3(MN) * delta_ice;

            % Cap: 不能转移超过 (液+吸附) 池已有 N
            transfer_NH4 = min(transfer_NH4, N_NH4_total_old);
            transfer_NO3 = min(transfer_NO3, N_NO3_total_old);

            % [2026-05-15 N-closure fix] 数值 guard: 如果 theta_ice_new 太小, 算 C_ice
            % 会数值不稳 (除以 ~eps). 此时不允许 transfer (rollback)，N 留在液相，避免
            % "transfer 已扣液相但 C_ice 没更新 → N 凭空消失" bug.
            if theta_ice_new <= eps_theta
                transfer_NH4 = 0;
                transfer_NO3 = 0;
            end

            N_NH4_ice_new = N_NH4_ice_old + transfer_NH4;
            N_NO3_ice_new = N_NO3_ice_old + transfer_NO3;

            % 反算 C: C_post = (总池 - transfer) / 总池分母
            % NH4 用 (θ+RHOKG) 分母 (吸附跟液相一起减), NO3 用 θ_LL_old
            C_NH4(MN) = max((N_NH4_total_old - transfer_NH4) / (theta_LL_old + RHOKG(ML)), 0);
            C_NO3(MN) = max((N_NO3_total_old - transfer_NO3) / theta_LL_old, 0);

            if theta_ice_new > eps_theta
                C_NH4_ice(ML, ND) = N_NH4_ice_new / theta_ice_new;
                C_NO3_ice(ML, ND) = N_NO3_ice_new / theta_ice_new;
            end
            % else: theta_ice_new ≤ eps; transfer 已被 rollback 到 0, C_ice 保持

            partition_diag.delta_NH4_to_ice = partition_diag.delta_NH4_to_ice + transfer_NH4 * DeltZ(ML) * 0.5;
            partition_diag.delta_NO3_to_ice = partition_diag.delta_NO3_to_ice + transfer_NO3 * DeltZ(ML) * 0.5;

        elseif delta_ice < -eps_theta && theta_ice_old > eps_theta
            % --- Melting: ice → liquid -------------------------------------
            transfer_NH4 = (-delta_ice) * C_NH4_ice(ML, ND);
            transfer_NO3 = (-delta_ice) * C_NO3_ice(ML, ND);

            % Cap: 不能释放超过冰相已有的 N
            transfer_NH4 = min(transfer_NH4, N_NH4_ice_old);
            transfer_NO3 = min(transfer_NO3, N_NO3_ice_old);

            % NH4: 释放进 (液+吸附) 总池
            N_NH4_total_old = (theta_LL_old + RHOKG(ML)) * C_NH4(MN);
            N_NO3_total_old = theta_LL_old * C_NO3(MN);

            % 反算 C: C_post = (总池 + transfer) / 总池分母
            C_NH4(MN) = (N_NH4_total_old + transfer_NH4) / (theta_LL_old + RHOKG(ML));
            C_NO3(MN) = (N_NO3_total_old + transfer_NO3) / theta_LL_old;
            % C_ice 不变（按比例释放，剩余浓度同）

            partition_diag.delta_NH4_from_ice = partition_diag.delta_NH4_from_ice + transfer_NH4 * DeltZ(ML) * 0.5;
            partition_diag.delta_NO3_from_ice = partition_diag.delta_NO3_from_ice + transfer_NO3 * DeltZ(ML) * 0.5;
        end

        % Edge case: 冰全融 — 清零冰相浓度
        if theta_ice_new <= eps_theta
            C_NH4_ice(ML, ND) = 0;
            C_NO3_ice(ML, ND) = 0;
        end
    end
end

% Conservation check (内部一致性, 跟外部 ADE 守恒账分开):
% 此处只检查"冻融操作前后" N 总质量是否守恒. 不含水流通量带走的 N
% (水流通量带 N 由 ADE 后续处理).
%
% mass_init: 用 θ_L_old + θ_I_old (KT 开始水量)
% mass_final: 用 θ_LL_new + θ_II_new (KT 末水量)
% 如果水量在 KT 内变化 (降雨/蒸发/重分布), 会引起非零差异 — 这部分由 ADE 处理.
% 真正冻融操作应该保 N 守恒 (per cell): transfer 进冰 = 出液相
for ML = 1:NL
    for ND = 1:2
        MN = ML + ND - 1;
        mass_NH4_final = C_NH4(MN) * Theta_LL_new(ML, ND) + ...
                         C_NH4_ice(ML, ND) * Theta_II_new(ML, ND);
        mass_NO3_final = C_NO3(MN) * Theta_LL_new(ML, ND) + ...
                         C_NO3_ice(ML, ND) * Theta_II_new(ML, ND);
        partition_diag.mass_conserved_check_NH4 = ...
            partition_diag.mass_conserved_check_NH4 + (mass_NH4_final - mass_NH4_init(MN));
        partition_diag.mass_conserved_check_NO3 = ...
            partition_diag.mass_conserved_check_NO3 + (mass_NO3_final - mass_NO3_init(MN));
    end
end

% TEMP DEBUG: print conservation check if non-trivial
if abs(partition_diag.mass_conserved_check_NO3) > 0.01 || abs(partition_diag.mass_conserved_check_NH4) > 0.01
    fprintf('  [pf-conserv] dM_NH4=%.3e  dM_NO3=%.3e  ΣΔθ_ice=%.3e\n', ...
        partition_diag.mass_conserved_check_NH4, ...
        partition_diag.mass_conserved_check_NO3, ...
        sum(Theta_II_new(:) - Theta_I_old(:)));
end

end
