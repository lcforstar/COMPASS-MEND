function [N_fruit_new, dN_harvest, harvested] = appleFruitHarvest( ...
    N_fruit, harvest_signal, harvest_threshold, harvest_fraction, already_harvested)
%APPLEFRUITHARVEST  One-shot Apple fruit harvest event (harvest_signal > threshold).
%
%   Returns:
%     N_fruit_new       [g N/m²]    N_fruit pool after harvest
%     dN_harvest        [g N/m²]    N removed in this KT (0 if not triggered)
%     harvested         [-, {0,1}]  1 if harvest triggered this KT (set already_harvested flag)
%
%   Inputs:
%     N_fruit           [g N/m²]    current fruit N pool
%     harvest_signal    [-]         maturity signal from caller (0/1). Phase D (P5 fix):
%                                   = phen_state.harvested_flag, set when independent fruit
%                                   clock DS_fruit > harvest_DS (0.95). NOT frPHU anymore —
%                                   frPHU caps at 0.648 at Changwu so it could never trigger.
%     harvest_threshold [-]         signal threshold (caller passes 0.5: signal is 0/1)
%     harvest_fraction  [-]         fraction of N_fruit removed at harvest (Apple: 0.90, 留 10% as seed)
%     already_harvested [-, {0,1}]  state flag: 1 if harvest already happened this year
%
%   Method:
%     - First time harvest_signal > threshold && N_fruit > 0 && not already_harvested → trigger
%     - dN_harvest = N_fruit × harvest_fraction (single event)
%     - Subsequent calls (already_harvested=1) → return unchanged
%     - State flag must be reset by caller at year boundary
%
%   Reference:
%     DS_fruit maturity clock: seasonDecidPhenology.m (DS_fruit > harvest_DS=0.95);
%       harvest_DS picks EPIC harvest index (Williams 1989 Trans.ASAE 32:497) ~0.9-0.95 range.
%     Mao 2018 — Apple N budget, commercial harvest removes ~90% of fruit N
%
%   Variable class (per compass-coding-standard §1): forcing transformer (state stays in caller)
%
%   Author: Zhu Jizixing 2026-05-21 (M2-B.4 S2.1; P5 fix 2026-05-26 frPHU→DS_fruit signal)

    if already_harvested || N_fruit <= 0 || harvest_signal < harvest_threshold
        N_fruit_new = N_fruit;
        dN_harvest  = 0;
        harvested   = 0;
        return;
    end

    dN_harvest  = N_fruit * harvest_fraction;
    N_fruit_new = N_fruit - dN_harvest;
    harvested   = 1;
end
