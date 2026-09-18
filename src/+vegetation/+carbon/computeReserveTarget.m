function [target_stem, target_root] = computeReserveTarget(C_stem_wood, croot_stem, nsc_frac)
%COMPUTERESERVETARGET  Physical target for the two woody NSC reserve pools (wood-based).
%
%   Each woody reserve pool's target = nsc_frac x its OWN wood carbon. The reserve is physically a
%   fraction (nsc_frac) of the non-structural carbohydrate held in perennial woody tissue, so the
%   aboveground (stem) reserve tracks stem wood and the belowground (coarse-root) reserve tracks
%   coarse-root wood. COMPASS has no explicit coarse-root organ pool, so coarse-root wood is
%   reconstructed as croot_stem x stem wood. The partition between the two pools therefore EMERGES
%   from the actual wood biomass ratio (root:stem = croot_stem) and evolves as the plant grows -- it
%   is NOT a frozen parameter. The target is the WINTER-MAXIMUM ("full") reserve level; the running
%   pool cycles below it (spring drawdown -> autumn refill), an emergent behaviour of the
%   target-driven refill in computeCarbonAllocation, not imposed here.
%
%   Returns:
%     target_stem [g C/m2]  aboveground (stem/branch wood) reserve target
%     target_root [g C/m2]  belowground (coarse-root wood) reserve target
%
%   Inputs:
%     C_stem_wood [g C/m2]  perennial aboveground wood C = C_livestem + C_deadstem (model state)
%     croot_stem  [-]       coarse-root : stem wood ratio (PFT table; apple 0.225 Zanotelli 2013
%                           woodBG/woodAG). Reconstructs coarse-root wood (no explicit pool in COMPASS).
%     nsc_frac    [-]       woody NSC carbon fraction (PFT table; winter-max NSC%DW x 0.44/0.47)
%
%   Reference:
%     Massoud et al. 2019 GMD 12:4133 (FATES target carbon storage = cushion x leaf; target IDEA)
%     fates/parteh/PRTAllometricCarbonMod.F90:598-613 (smooth refill toward target)
%     Martinez-Vilalta et al. 2016 Ecol Monogr 86:495 (woody NSC by functional type)
%     Zanotelli et al. 2013 Biogeosciences 10:3089 Table 2 (croot_stem = woodBG/woodAG)
%
%   Variable class (per compass-coding-standard §1): C_stem_wood = state slice; croot_stem/nsc_frac =
%   parameters; target_* = derived auxiliary (consumed by the reserve refill).
%
%   Author: Zhu Jizixing 2026-06-10

    C_coarseroot_wood = croot_stem * max(C_stem_wood, 0);   % [g C/m2] reconstructed coarse-root wood
    target_stem = nsc_frac * max(C_stem_wood, 0);           % [g C/m2] aboveground reserve target
    target_root = nsc_frac * C_coarseroot_wood;             % [g C/m2] belowground reserve target
end
