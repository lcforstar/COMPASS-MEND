function [N_leaf_new, N_storage_new, dN_litter] = leafSenescenceRetransl( ...
    N_leaf, N_storage, sene_signal, LAI, ndays_off, retransl_fraction, ...
    sene_threshold, dt, offset_counter)
%LEAFSENESCENCERETRANSL  CLM5-style leaf senescence + N retranslocation.
%
%   Verbatim port of CLM5 CNPhenologyMod.F90 offset behavior:
%     - For first (ndays_off - 1) days of offset ramp: linear drain at rate
%       1/ndays_off per day
%     - On the LAST day (offset_counter == dt, i.e. counter just hit 0):
%       drain all remaining leafC (= leafc / dt × 1 step in CLM5)
%
%   Returns:
%     N_leaf_new        [g N/m² ground]   after this step's drain
%     N_storage_new     [g N/m² ground]   += d_N_leaf × retransl_fraction
%     dN_litter         [g N/m² ground]   = d_N_leaf × (1 - retransl_fraction)
%
%   Inputs:
%     N_leaf, N_storage [g N/m² ground]
%     sene_signal       [-]      binary 0/1, =1 during offset_flag=true window
%     LAI               [m²/m²]  [UNUSED 2026-06-05] kept for positional-call
%                                compat; drain no longer gated on LAI (see E-fix
%                                below — CLM5 gates on offset_flag, not LAI).
%     ndays_off         [days]   CLM5 ndays_off (=15 for Apple), determines rate
%     retransl_fraction [-]      Apple 0.7 (Aerts 1996)
%     sene_threshold    [-]      binary threshold, typically 0.5
%     dt                [s]      time step (typically 3600)
%     offset_counter    [days]   CLM5 phen_state.offset_counter; if = 1
%                                (or <= 1) this is the last day → full drain.
%                                Optional: if missing/empty, use linear ramp only.
%
%   Method (CLM5 CNPhenologyMod.F90:2448-2453 verbatim):
%     if offset_counter == dt (== 1 day in our discrete units):
%        leafc_to_litter = (1/dt) × leafc   ← drain ALL remaining
%     else:
%        d_leafc = leafc × (1/ndays_off) × (dt/86400)   ← linear ramp
%
%     N follows leafC proportionally; nrelocfrac to storage, rest to litter.
%
%   Reference:
%     CLM5 CNPhenologyMod.F90:2448-2453 (offset_counter == dt → flush residual)
%     CLM5 CNPhenologyMod.F90:808-820  (offset_counter linear decrement)
%     Aerts R (1996) Oikos 67:469      (retransl_fraction meta-analysis)
%
%   Author: Zhu Jizixing 2026-05-21; Phase D 2026-05-26 CLM5 strict port

    % [E-fix 2026-06-05] Gate on offset_flag (sene_signal) + leaf pool only, NOT
    % displayed LAI. CLM5 CNPhenologyMod.F90:2449-2474 drains leafc_to_litter purely
    % on offset_counter/offset_flag with no LAI guard. The old LAI<0.1 check stranded
    % leaf C in the leafc_sla engine, where LAI=phen*C_leaf*SLA collapses to ~0 as phen
    % ramps down during the very offset window we need to drain (leaf C frozen all winter).
    if sene_signal < sene_threshold || N_leaf <= 0
        N_leaf_new    = N_leaf;
        N_storage_new = N_storage;
        dN_litter     = 0;
        return;
    end

    % [leaf-fall-linear 2026-06-08] TRUE linear drain: divide the CURRENT leaf pool by the
    %   REMAINING offset days (offset_counter), not by the fixed ndays_off. This keeps the
    %   absolute amount shed per day CONSTANT (= N_leaf_at_offset_start / ndays_off) and drives
    %   the pool to exactly 0 on the last offset day, matching the linear phen ramp
    %   (seasonDecidPhenology.m: phen -= 1/ndays_off) that drives displayed LAI -> leaf C/N and
    %   LAI reach 0 together, no overwinter residual. Observed Changwu apple leaf-fall is linear
    %   (DOY243-329, R^2 0.74 linear vs 0.53 exponential, LAI xlsx 247obs).
    %   CLM5 CNPhenologyMod.F90:808-820 (offset ramp) + :2450-2453 (last-day flush).
    %   OLD form  d_N = N_leaf * (1/ndays_off)  was a fixed fraction of the SHRINKING pool
    %   (geometric decay): after ndays_off ~ (1-1/n)^n ~ 37% stranded, then frozen all winter
    %   once offset_flag cleared (v61 ndays_off=45 residual bug). Using remaining-days as the
    %   divisor turns the same per-step call into a constant absolute flux, no need to store
    %   the offset-start pool.
    if nargin >= 9 && ~isempty(offset_counter) && offset_counter >= 1
        days_left = offset_counter;                  % [days] remaining offset days incl. today
    else
        days_left = max(ndays_off, 1);               % fallback: no counter -> assume full window
    end
    d_N_leaf = N_leaf / days_left * (dt / 86400);    % [g N/m2/step] constant daily abs. flux

    % CLM5 line 2450-2453: last offset day (counter <= 1) -> flush ALL remaining N so nothing
    % strands over winter. With the constant-flux form this is already ~N_leaf; clamp for safety.
    if nargin >= 9 && ~isempty(offset_counter) && offset_counter <= 1
        d_N_leaf = N_leaf;
    end

    d_N_leaf = min(d_N_leaf, N_leaf);

    dN_retrans    = d_N_leaf * retransl_fraction;
    dN_litter     = d_N_leaf * (1 - retransl_fraction);
    N_leaf_new    = N_leaf - d_N_leaf;
    N_storage_new = N_storage + dN_retrans;
end
