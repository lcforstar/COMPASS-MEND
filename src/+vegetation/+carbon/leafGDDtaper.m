function f_taper = leafGDDtaper(gdd_budburst, gdd_lfmat, bfact)
%LEAFGDDTAPER  Leaf-allocation taper factor declining 1->0 with thermal time.
%
%   Multiplier applied to the (otherwise constant) leaf allocation fraction so
%   leaf carbon allocation ramps DOWN to ~0 once the canopy reaches thermal
%   maturity (gdd_budburst = gdd_lfmat). Fixes the monotonic C_leaf overshoot
%   of a constant f_alloc_leaf: without it leaf C grows all season; with it leaf
%   C plateaus by canopy maturity and the caller redirects the freed fraction to
%   fruit/wood. Used by BOTH the carbon allocation and the N-demand pre-pass so
%   the C-dominant organ shares stay consistent.
%
%   Returns:
%     f_taper  [-]  leaf-allocation multiplier in [0,1]; 1 at bud-break, 0 at gdd_lfmat
%
%   Inputs:
%     gdd_budburst [°C·d]  growing-degree-days since bud-break (COMPASS base 5°C)
%     gdd_lfmat    [°C·d]  GDD from bud-break to canopy maturity (taper->0 endpoint)
%     bfact        [-]     J-Norman exponential shape factor (<0 => convex decline)
%
%   Method:
%     Dombrowski et al. (2022) GMD 15:5167 Eq. 4 (J. Norman leaf-allocation decline);
%     CLM5-FruitTree NutrientCompetitionFlexibleCNMod.F90:1663-1666 (aleaf), normalized
%     to 1 at g=0 and 0 at g=gdd_lfmat:
%       f_taper = (exp(-b) - exp(-b * g/gdd_lfmat)) / (exp(-b) - 1)
%     g clamped to [0, gdd_lfmat]; bfact~0 falls back to a linear ramp.
%
%   Author: Zhu Jizixing 2026-06-03 (F2 / stage-13a apple leaf allocation)

    if ~isfinite(gdd_budburst) || gdd_lfmat <= 0
        f_taper = 1.0;            % pre-budburst / undefined GDD: no taper (caller gates use)
        return;
    end

    g = max(0, min(gdd_budburst, gdd_lfmat));     % [°C·d] clamp to [0, gdd_lfmat]

    % Dombrowski et al. 2022 GMD 15:5167 Eq.4 / CLM5-FruitTree
    % NutrientCompetitionFlexibleCNMod.F90:1663-1666 (normalized J-Norman decline).
    eb = exp(-bfact);
    if abs(eb - 1.0) < 1e-9                        % bfact ~ 0 -> linear fallback
        f_taper = 1.0 - g / gdd_lfmat;
    else
        f_taper = (eb - exp(-bfact * g / gdd_lfmat)) / (eb - 1.0);
    end
    f_taper = max(0.0, min(1.0, f_taper));         % [-] guard to [0,1]
end
