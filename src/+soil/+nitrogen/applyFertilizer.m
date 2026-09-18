function [SoilVariables, fert_diag] = applyFertilizer(SoilVariables, BoundaryCondition, current_date, SUMDELTZ, NL, NN, DeltZ, NitrogenVariables)
%APPLYFERTILIZER  Check + apply fertilizer N injection at current_date; report mass added.
%
% Triggered at the start of each outer step. Looks up BoundaryCondition.Fertilizer.controller
% to see if any prescribed fertilizer event matches current_date. If yes, splits the event
% N amount across NH4 / NO3 / urea pools per the configured ratios and injects into
% the soil-layer concentration arrays at the specified depth band.
%
% Inputs:
%   SoilVariables    [inout]  struct - C_NH4_new / C_NO3_new / C_Urea_new updated (NN×1, μg/cm³)
%   BoundaryCondition [in]    struct - .Fertilizer.controller / DEPTHSTART / DEPTHEND
%   current_date     [in]     datetime scalar - current outer-step timestamp
%   SUMDELTZ         [in]     [cm, NL×1] - cumulative layer thickness from surface
%   NL, NN           [in]     [-]      - layer / node counts
%
% Returns:
%   SoilVariables   updated concentrations
%   fert_diag       struct with cumulative-conservation accounting:
%       .applied         [bool]            true if this step had fertilization
%       .nh4_added       [μg/cm²]          total NH4 mass added this step (column integral)
%       .no3_added       [μg/cm²]          total NO3 mass added this step
%       .urea_added      [μg/cm²]          total urea mass added this step (separate pool)
%
% Note (Skill §6 cumulative conservation):
%   Caller should accumulate fert_diag.nh4_added / no3_added into the cum_nh4_in /
%   cum_no3_in ledger in solveNitrogen.m, so the N mass balance check sees fertilizer
%   as an input source (previously it was a hidden source and broke closure).

fert_diag = struct('applied', false, 'nh4_added', 0, 'no3_added', 0, 'urea_added', 0);

[should_apply, F_conc, fR] = BoundaryCondition.Fertilizer.controller(current_date);
if ~should_apply, return; end

fert_diag.applied = true;

[~, iStart] = min(abs(SUMDELTZ - BoundaryCondition.Fertilizer.DEPTHSTART));
[~, iEnd]   = min(abs(SUMDELTZ - BoundaryCondition.Fertilizer.DEPTHEND));

V_control = zeros(NN, 1);
V_control(NL-iEnd+1 : NL-iStart+1) = 1;   % 底到顶节点编号
F_concnode = F_conc * V_control;            % NN×1, µg/cm³

% Snapshot pre-fert N storage to compute exact mass-added (matches calculateStorage definition)
% This ensures fert_diag.*_added has the SAME units as cv_nh4 / cv_no3, so conservation balance works.
%
% [2026-05-17 dissolution pool] NH4 storage now includes granule pool:
%     cv_nh4_total = cv_nh4_dissolved + cv_nh4_granule
% so that fert_diag.nh4_added reports the TOTAL NH4 mass injected into the system
% (granule + dissolved), regardless of which sub-pool currently holds it.
% The granule pool subsequently drains to dissolved pool by nh4FertDissolution.
Theta_LL_pre = SoilVariables.Theta_LL;
RHOKG = NitrogenVariables.RHOKG;
has_granule = isfield(SoilVariables, 'C_NH4_granule');
cv_nh4_pre = 0;  cv_no3_pre = 0;  cv_nh4_granule_pre = 0;
for ML_p = 1:NL
    cv_nh4_pre = cv_nh4_pre + DeltZ(ML_p) * ( ...
        (Theta_LL_pre(ML_p,1) + RHOKG(ML_p)) * SoilVariables.C_NH4_new(ML_p) + ...
        (Theta_LL_pre(ML_p,2) + RHOKG(ML_p)) * SoilVariables.C_NH4_new(ML_p+1)) / 2;
    cv_no3_pre = cv_no3_pre + DeltZ(ML_p) * ( ...
        Theta_LL_pre(ML_p,1) * SoilVariables.C_NO3_new(ML_p) + ...
        Theta_LL_pre(ML_p,2) * SoilVariables.C_NO3_new(ML_p+1)) / 2;
    if has_granule
        % Granule pool: μg/cm³ × layer thickness = μg/cm² (no Theta_LL scaling — it's a SURFACE-pool concentration)
        cv_nh4_granule_pre = cv_nh4_granule_pre + DeltZ(ML_p) * ( ...
            SoilVariables.C_NH4_granule(ML_p) + SoilVariables.C_NH4_granule(ML_p+1)) / 2;
    end
end

% [2026-05-17 dissolution pool] NH4 颗粒肥不直接进 NH4 溶解池, 先进 granule 缓冲池
% 然后通过 nh4FertDissolution.m 按一阶速率缓溶到 C_NH4 (mirrors urea→NH4 hydrolysis path).
% 物理理由: 真实 NH4 颗粒肥 ((NH4)2SO4, NH4Cl) 颗粒缓溶 → 表层 NH4 浓度局部峰值 ↑
% → Henry 驱动 NH3 挥发 ↑. 瞬时均匀稀释假设让 Svol 通量低估 10-30×.
% 参考: Probert 1998 APSIM-SoilN; Sommer & Hutchings 2001 EJA.
%
% NO3 颗粒不走此机制 (NO3 不吸附 + 无挥发驱动 → 瞬时混合物理正确).
if fR.NH4_ratio  > 0
    if isfield(SoilVariables, 'C_NH4_granule')
        SoilVariables.C_NH4_granule = SoilVariables.C_NH4_granule + F_concnode * fR.NH4_ratio;
    else
        % Backward compat: 若未启用 dissolution pool, 退回老瞬时混合
        SoilVariables.C_NH4_new  = SoilVariables.C_NH4_new  + F_concnode * fR.NH4_ratio;
    end
end
if fR.NO3_ratio  > 0
    SoilVariables.C_NO3_new  = SoilVariables.C_NO3_new  + F_concnode * fR.NO3_ratio;
end
if fR.urea_ratio > 0
    SoilVariables.C_Urea_new = SoilVariables.C_Urea_new + F_concnode * fR.urea_ratio;
end

% Recompute post-fert storage; (post - pre) = mass added, in EXACTLY the same units as cv_nh4/cv_no3
cv_nh4_post = 0;  cv_no3_post = 0;  cv_nh4_granule_post = 0;
for ML_p = 1:NL
    cv_nh4_post = cv_nh4_post + DeltZ(ML_p) * ( ...
        (Theta_LL_pre(ML_p,1) + RHOKG(ML_p)) * SoilVariables.C_NH4_new(ML_p) + ...
        (Theta_LL_pre(ML_p,2) + RHOKG(ML_p)) * SoilVariables.C_NH4_new(ML_p+1)) / 2;
    cv_no3_post = cv_no3_post + DeltZ(ML_p) * ( ...
        Theta_LL_pre(ML_p,1) * SoilVariables.C_NO3_new(ML_p) + ...
        Theta_LL_pre(ML_p,2) * SoilVariables.C_NO3_new(ML_p+1)) / 2;
    if has_granule
        cv_nh4_granule_post = cv_nh4_granule_post + DeltZ(ML_p) * ( ...
            SoilVariables.C_NH4_granule(ML_p) + SoilVariables.C_NH4_granule(ML_p+1)) / 2;
    end
end
% [2026-05-17 dissolution pool] fert_diag.nh4_added reports TOTAL NH4 injected
% (dissolved+granule). Conservation ledger then treats granule→dissolved as
% internal redistribution (Sdis source for dissolved, Sdis sink for granule),
% with mass balance: d(cv_nh4 + cv_nh4_granule)/dt = ... (external sources only).
fert_diag.nh4_added  = (cv_nh4_post - cv_nh4_pre) + (cv_nh4_granule_post - cv_nh4_granule_pre);
fert_diag.no3_added  = cv_no3_post - cv_no3_pre;
fert_diag.nh4_to_granule = cv_nh4_granule_post - cv_nh4_granule_pre;   % diagnostic split
% Urea pool tracked separately; not added to NH4/NO3 ledger directly.
% When urea hydrolyzes it adds to NH4 — that source should be tracked via DeltaUrea.
fert_diag.urea_added = sum(F_concnode * fR.urea_ratio);  % diagnostic only
end
