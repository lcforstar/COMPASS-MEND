function [N_root_new, dN_root_litter] = fineRootTurnover( ...
    N_root_fine, root_turnover_rate, dt, phen_gate)
%FINEROOTTURNOVER  Daily fine-root N turnover with phenology gating.
%
%   Returns:
%     N_root_new      [g N/m2 ground]
%     dN_root_litter  [g N/m2 ground]   sent to SOM (cum_N_root_litter)
%
%   Inputs:
%     N_root_fine        [g N/m2 ground]
%     root_turnover_rate [1/day]  Apple 0.00274 ~ 1/365 (Wells & Eissenstat 2002)
%     dt                 [s]
%     phen_gate          [-, 0..1]  [F5 2026-05-31] phenology gate on turnover rate;
%                        optional (missing/empty -> 1.0, year-round, backward compatible).
%                        CLM5 ties fine-root litterfall to leaf bglfr (=0 in dormancy):
%                        dormant period drops fewer roots, growing season full rate.
%
%   Reference:
%     Wells CE, Eissenstat DM (2002) Plant Soil 21:324 - Apple fine root ~year turnover
%     CLM5 CNPhenologyMod.F90:2623 - frootc_to_litter = bglfr*frootc (same bglfr as leaf)
%
%   Variable class: forcing transformer (caller owns state)
%
%   Author: Zhu Jizixing 2026-05-21 (M2-B.4 S2.3); F5 phen gate 2026-05-31

    if nargin < 4 || isempty(phen_gate)
        phen_gate = 1.0;                 % backward compatible: year-round
    end

    if N_root_fine <= 0 || root_turnover_rate <= 0
        N_root_new     = N_root_fine;
        dN_root_litter = 0;
        return;
    end

    dN_root_litter = N_root_fine * root_turnover_rate * (dt / 86400) * phen_gate;
    N_root_new     = N_root_fine - dN_root_litter;
end
