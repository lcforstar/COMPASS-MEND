function [Theta_Solute_LL, Theta_Solute_L] = effectiveLiquidWaterForSolute( ...
        SoilVariables, cfg)
%EFFECTIVELIQUIDWATERFORSOLUTE  Effective liquid water content available for N transport.
%
% Solute (NH4 / NO3) ADE uses an "effective" liquid θ that:
%   (a) caps at porosity (so over-saturated nodes don't inflate transport),
%   (b) handles saturated surface ponding (top node hh > 0 expands beyond θ_s).
%
% **Important convention** (COMPASS state semantics, see also defineSoilVariables.m):
%   Theta_LL = liquid-phase water content       (already excludes ice)
%   Theta_II = ice-phase water content (volume) (separate state)
%   Theta_UU = total water = Theta_LL + Theta_II × (ρ_I/ρ_L)  (composite, NOT used here)
%
% Since Theta_LL is ALREADY the liquid phase by definition, we do NOT subtract
% Theta_II again here — that would double-count the freeze-out and inflate
% solute concentration spuriously. The 2026-05-13 first attempt at "ice-subtraction"
% was reverted after clarification of UU/LL/II semantics.
%
% The proper freeze-thaw N effect is captured by:
%   - Theta_LL itself shrinking when water freezes (Richards solver handles this)
%   - Partition of N into ice phase: soil.nitrogen.chemistry.freezeThawPartitioning
%     (P0a.5 — ice-phase N as separate prognostic state C_NH4_ice / C_NO3_ice)
%
% Inputs:
%   SoilVariables [in]  struct  fields: Theta_LL, Theta_L, POR, Theta_s, hh, h
%   cfg           [in]  struct  fields: Grid.NL, Grid.NN, Grid.nD, Grid.DeltZ
%
% Returns:
%   Theta_Solute_LL [m³/m³, NL×nD]  effective liquid θ at new step for solute transport
%   Theta_Solute_L  [m³/m³, NL×nD]  effective liquid θ at prior step
%
% Reference:
%   HYDRUS-1D SOLUTE.f90 — ponding-aware top branch
%   STEMMUS state convention: UU = total water, LL = liquid only, II = ice only
%
% Added 2026-05-13. Originally tried ice-subtraction (since reverted) after
% recognizing Theta_LL semantics; kept as a function-extracted module per
% Skill §4.1 (no more inline state-clamp scattered in solveNitrogen.m).

NL    = cfg.Grid.NL;
NN    = cfg.Grid.NN;
nD    = cfg.Grid.nD;
DeltZ = cfg.Grid.DeltZ;

Theta_Solute_LL = SoilVariables.Theta_LL;
Theta_Solute_L  = SoilVariables.Theta_L;

% Idempotency note (2026-05-15): solveNitrogen now stashes Theta_LL_for_N_prev
% as the EFFECTIVE Theta_Solute_LL (post-cap/post-ponding) so cvno3OLD aligns
% with cvno3(KT-1). When this function receives such an already-effective input
% (interior nodes ≤ POR by construction; top node may exceed POR via ponding),
% we must NOT re-clip — that would destroy the ponding adjustment.
% Detection: if Theta_L(top, ND) > Theta_s(ML), it was already ponding-corrected
% upstream — leave both top and interior alone (they're consistent already).
already_effective_prev = ...
    (SoilVariables.Theta_L(NL, 2) > SoilVariables.Theta_s(NL)) || ...
    (SoilVariables.Theta_L(NL, 1) > SoilVariables.Theta_s(NL));

for ML = 1:NL
    for ND = 1:nD
        MN = ML + ND - 1;
        if MN < NN
            % Interior nodes: cap effective liquid water at POR (over-saturation guard)
            % Theta_LL is already liquid-only (ice = Theta_II is separate).
            Theta_Solute_LL(ML, ND) = min(SoilVariables.Theta_LL(ML, ND), SoilVariables.POR(ML));
            if already_effective_prev
                Theta_Solute_L(ML, ND)  = SoilVariables.Theta_L(ML, ND);   % already effective
            else
                Theta_Solute_L(ML, ND)  = min(SoilVariables.Theta_L(ML, ND),  SoilVariables.POR(ML));
            end
        elseif MN == NN
            % Top node: ponding-aware expansion when surface is saturated
            if SoilVariables.hh(NN) > 0
                Theta_Solute_LL(ML, ND) = SoilVariables.Theta_s(ML) + SoilVariables.hh(NN) / DeltZ(ML);
            else
                Theta_Solute_LL(ML, ND) = SoilVariables.Theta_LL(ML, ND);
            end
            if already_effective_prev
                % Theta_L already carries the prev-step ponding correction; keep as-is.
                Theta_Solute_L(ML, ND) = SoilVariables.Theta_L(ML, ND);
            elseif SoilVariables.h(NN) > 0
                Theta_Solute_L(ML, ND) = SoilVariables.Theta_s(ML) + SoilVariables.h(NN) / DeltZ(ML);
            else
                Theta_Solute_L(ML, ND) = SoilVariables.Theta_L(ML, ND);
            end
        end
    end
end

end
