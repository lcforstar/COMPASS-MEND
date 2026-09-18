function [state, diag] = seasonDecidPhenology(state, pftParams, forcing, doy, ~)
%SEASONDECIDPHENOLOGY  Seasonal-deciduous phenology (Sykes 1996 + CLM5 daylength).
%
%   Advances the per-PFT phenology state by one day for a SEASON_DECID PFT
%   (e.g., Apple PFT 23, Que_rob, Bet_pen). Combines:
%     1. LPJ-GUESS Sykes 1996 chilling-GDD coupled budburst (growth.cpp:91-99)
%     2. LPJ-GUESS aphen / APHEN_MAX leaf-on day cap (growth.cpp:148-152)
%     3. LPJ-GUESS cold-month gdd5/chilldays reset (driver.cpp:583-586)
%     4. CLM5 daylength-triggered senescence (CNPhenologyMod.F90:957)
%
%   Returns:
%     state    [inout]  PhenologyState struct, all fields advanced by 1 day:
%                       chilldays, gdd5, phen, aphen, onset_flag, offset_flag,
%                       dormant_flag, budburst_doy, senescence_doy updated.
%                       Source: see PhenologyState.m for variable definitions.
%     diag     struct   per-step diagnostics (only-write, never read by physics):
%                       .gdd0_threshold  [°C·d] today's Sykes 1996 threshold
%                       .ramp_progress   [-] (gdd5-gdd0)/phengdd5ramp clamp [0,1]
%                       .reset_triggered [bool] cold-month reset fired today
%
%   Inputs:
%     state         [inout]  PhenologyState (per-PFT, persists across days)
%     pftParams     [in]     struct, required fields (from PFT 表):
%                              .k_chilla        [°C·d]  Sykes 1996 a
%                              .k_chillb        [°C·d]  Sykes 1996 b
%                              .k_chillk        [d^-1]  Sykes 1996 k
%                              .phengdd5ramp    [°C·d]  ramp width to phen=1
%                              .crit_dayl_h     [h]     senescence dayl threshold
%     forcing       [in]     struct, current-day weather:
%                              .T_avg           [°C]    daily mean air temp
%                              .dayl_h          [h]     day length (sunrise→sunset)
%                              .mtemp           [°C]    monthly running-mean temp
%     doy           [in]     [-]   day-of-year 1-365 (for diagnostic logging)
%     year_mtemp_max [in]    [°C]  this-year max monthly temp seen so far
%                                  (used to latch ifsensechill flag, LPJ-GUESS
%                                  driver.cpp:540-545)
%
%   Method (LPJ-GUESS Sykes 1996, growth.cpp:82-101, verbatim ported):
%
%     g(T) accumulation:
%       gdd5 += max(0, T_avg - 5)              [°C·d]
%       if T_avg < 5:  chilldays += 1          [d]
%
%     ifsensechill latch (LPJ-GUESS driver.cpp:540-545):
%       if mtemp == year_mtemp_max:  ifsensechill = true
%
%     Cold-month reset (LPJ-GUESS driver.cpp:583-586):
%       if mtemp_last >= 5 AND mtemp < 5 AND ifsensechill:
%         gdd5 = 0;  chilldays = 0;  ifsensechill = false
%         dormant_flag = true;  offset_flag = false
%
%     Sykes 1996 chilling-GDD threshold (growth.cpp:92-93):
%       gdd0 = k_chilla + k_chillb × exp(-k_chillk × chilldays)
%
%     Budburst & phen ramp (growth.cpp:95-99):
%       if gdd5 > gdd0 AND aphen < APHEN_MAX (=210):
%         phen = min(1, (gdd5 - gdd0) / phengdd5ramp)
%         if dormant_flag:  onset_flag = true; dormant_flag = false; budburst_doy = doy
%       else:
%         phen = 0
%
%     CLM5 daylength senescence (CNPhenologyMod.F90:957):
%       if dayl_h < crit_dayl_h AND NOT dormant_flag AND aphen > 30 (in growing season):
%         offset_flag = true
%         senescence_doy = doy  (only set on first trigger this year)
%
%     aphen accumulation (growth.cpp:152):
%       aphen += phen  (caller resets at year boundary via APHEN_MAX hit + reset)
%
%   Hard constants:
%     APHEN_MAX        = 210   [days]   Max leaf-on days (LPJ-GUESS growth.cpp:50)
%     T_CHILL_BASE     = 5     [°C]     Chill accumulator threshold
%     T_GDD_BASE       = 5     [°C]     GDD accumulator base (LPJ-GUESS uses 5)
%     APHEN_GROW_MIN   = 30    [days]   Min aphen before daylength senescence
%                                       can trigger (prevents false trigger
%                                       early in season when dayl shortening
%                                       coincides with first leaf-out)
%
%   References:
%     LPJ-GUESS v4.1.1 (D:/OneDrive/LPJ-GUESS-v4.1.1/guess_4.1/):
%       modules/growth.cpp:50,82-101,148-152
%       modules/driver.cpp:540-551, 583-586
%       (Sykes MT, Prentice IC, Cramer W (1996) J Biogeogr 23:203-233 cited at
%        growth.cpp:89; PDF not locally available — source-code citation)
%     CLM5 CTSM release-clm5.0.37:
%       src/biogeochem/CNPhenologyMod.F90:67 (crit_dayl global), 957 (offset trigger)
%     CLM5 PFT NetCDF: data/references/clm5/clm5_params.c171117.nc
%       (crit_dayl=39300 s = 10.9 h, scalar global for all PFTs)
%
%   Created 2026-05-26 (Phase D phenology unification, Step 4).

% --- Hard constants (LPJ-GUESS / CLM5 verbatim) -----------------------------
APHEN_MAX             = 210;   % [days] LPJ-GUESS growth.cpp:50
T_CHILL_BASE          = 5;     % [°C]   chill accumulator threshold
T_GDD_BASE            = 5;     % [°C]   gdd5 base (LPJ-GUESS verbatim)
APHEN_GROW_MIN        = 30;    % [days] min leaf-on before senescence dayl trigger
% Hemisphere-aware COLDEST/WARMEST days (LPJ-GUESS framework/guess.h:178-194):
%   N hemisphere: COLDEST = Jan 15 (DOY 14), WARMEST = Jul 15 (DOY 195)
%   S hemisphere: COLDEST = Jul 15 (DOY 195), WARMEST = Jan 15 (DOY 14)
% Dispatch by sign of latitude (LPJ-GUESS growth.cpp:148 dispatch pattern).
COLDEST_DAY_NHEMISPHERE = 14;
COLDEST_DAY_SHEMISPHERE = 195;
WARMEST_DAY_NHEMISPHERE = COLDEST_DAY_SHEMISPHERE;   % = 195
WARMEST_DAY_SHEMISPHERE = COLDEST_DAY_NHEMISPHERE;   % = 14

% Resolve hemisphere from forcing.latitude_deg (default to N if absent)
if isfield(forcing, 'latitude_deg') && forcing.latitude_deg < 0
    coldest_day = COLDEST_DAY_SHEMISPHERE;
    warmest_day = WARMEST_DAY_SHEMISPHERE;
else
    coldest_day = COLDEST_DAY_NHEMISPHERE;
    warmest_day = WARMEST_DAY_NHEMISPHERE;
end

% --- Defensive guard: only for SEASON_DECID PFT ------------------------------
% Sykes 1996 params must be numeric (not NaN). NaN means this PFT is not
% season_decid and the dispatcher routed it here by mistake.
if any(isnan([pftParams.k_chilla, pftParams.k_chillb, ...
              pftParams.k_chillk, pftParams.phengdd5ramp]))
    error('seasonDecidPhenology:nonSeasonDecidPFT', ...
        ['Sykes 1996 params contain NaN — this PFT is not season_decid. ', ...
         'Dispatcher should route non-season_decid PFTs elsewhere.']);
end

% --- 1. Latch ifsensechill on WARMEST_DAY (LPJ-GUESS driver.cpp:540-545) ----
% [P2 fix 2026-05-26] Use fixed hemisphere-specific DOY (NH: Jul 15 = DOY 195;
% SH: Jan 15 = DOY 14) instead of a running max — matches LPJ-GUESS
% guess.h:190-194 + driver.cpp:540-545 exactly.
if doy == warmest_day
    state.ifsensechill = true;
end

% --- 2a. aphen reset at COLDEST_DAY (LPJ-GUESS growth.cpp:148-152) ---------
% [P2 fix 2026-05-26] aphen reset is a FIXED hemisphere-specific DOY event
% (NH: Jan 15 = DOY 14; SH: Jul 15 = DOY 195). LPJ-GUESS growth.cpp:148-149
% dispatches by latitude sign. Previously bundled into cold-month reset, which
% fired on mtemp dropping below 5°C — often DOY 320+ and inconsistent.
if doy == coldest_day
    state.aphen = 0;
end

% --- 2b. Cold-month reset (LPJ-GUESS driver.cpp:583-586) -------------------
% Reset gdd5/chilldays/onset_gdd when monthly mean temp crosses below 5°C,
% gated by ifsensechill (must have passed WARMEST_DAY=195 this year).
% Does NOT reset aphen (that's in 2a above, COLDEST_DAY=14 only).
reset_triggered = false;
% [leaf-flush 2026-06-08] leaf_flush_pending: one-shot signal to the turnover module
% that the cold-month reset forced dormancy WHILE leaf-fall was still ramping down, so
% the displayed leaf C/N pool has not finished draining. Without it, clearing offset_flag
% below stops leafSenescenceRetransl (gated on offset_flag) mid-ramp and strands the
% residual leaf C/N over winter (the v62 ndays_off=45 vs cold-reset timing conflict:
% leaf-fall needs ndays_off days but the reset fires ~day 32 at Changwu, freezing ~24 gC
% of leaf C). A cold snap defoliating a deciduous tree in one step is physically standard
% (CLM5 CNPhenologyMod.F90:2450-2453 last-day flush drains all remaining leafc at once).
state.leaf_flush_pending = false;
if ~isnan(state.mtemp_last) && ~isnan(forcing.mtemp) ...
        && state.mtemp_last >= T_CHILL_BASE ...
        && forcing.mtemp < T_CHILL_BASE ...
        && state.ifsensechill
    if state.offset_flag                 % reset interrupts an active leaf-fall ramp
        state.leaf_flush_pending = true; % -> turnover flushes residual leaf C/N this step
    end
    state.gdd5          = 0;
    state.chilldays     = 0;
    state.onset_gdd     = 0;
    state.onset_gddflag = 0;
    state.ifsensechill  = false;
    state.dormant_flag  = true;
    state.offset_flag   = false;
    state.chill_day      = 0;        % [13b] reset Cesaraccio chill accumulators (apple)
    state.anti_chill_day = 0;
    reset_triggered     = true;
end
state.mtemp_last = forcing.mtemp;

% --- 3a. ws_flag: winter→summer solstice flag (CLM5 line 802-807) -----------
% ws_flag = 1 if dayl rising (spring half-year, winter→summer solstice)
% ws_flag = 0 if dayl falling (autumn half-year, summer→winter solstice)
% First call (prev_dayl NaN): conservative ws_flag = 0 (don't trigger gddflag)
if ~isnan(state.prev_dayl)
    state.ws_flag = (forcing.dayl_h >= state.prev_dayl);
else
    state.ws_flag = 0;
end

% --- 3b. onset_gddflag transitions (CLM5 line 885-898) -----------------------
% Turn ON at winter solstice (dormant + ws_flag flip to 1 + gdd_flag was 0)
if state.dormant_flag && state.onset_gddflag == 0 && state.ws_flag == 1
    state.onset_gddflag = 1;
    state.onset_gdd     = 0;       % start fresh GDD accumulator
end

% Turn OFF at summer solstice if threshold not reached (winter→spring failed)
if state.dormant_flag && state.onset_gddflag == 1 && state.ws_flag == 0
    state.onset_gddflag = 0;
    state.onset_gdd     = 0;
end

% --- 3c. Daily accumulators -------------------------------------------------
% gdd5: LPJ-GUESS-style continuous (used for legacy diagnostics and ramp computation)
% onset_gdd: CLM5-style spring-half-year-only, gated by onset_gddflag
% chilldays: LPJ-GUESS-style winter cold accumulator
state.gdd5 = state.gdd5 + max(0, forcing.T_avg - T_GDD_BASE);

% CLM5 line 905-908: only accumulate onset_gdd when ws_flag==1 AND gddflag==1
% (i.e., during spring half-year after winter solstice latch)
% Use T_air (forcing.T_avg) above freezing as gate (CLM5 uses soil temp; here T_air proxy).
if state.onset_gddflag == 1 && forcing.T_avg > T_GDD_BASE
    state.onset_gdd = state.onset_gdd + (forcing.T_avg - T_GDD_BASE);
end

% chilldays: count days with T_avg < 5°C (no reset here, only at cold-month event)
if forcing.T_avg < T_CHILL_BASE
    state.chilldays = state.chilldays + 1;
end

% [13b] Cesaraccio chill/anti-chill accumulation for apple bud-burst (CLM5-FruitTree
% CNPhenologyMod.F90:1824-1862). chill_day (<=0) accumulates toward crequ; once met,
% anti_chill_day (>=0) accumulates the post-chilling warmth that triggers bud-burst.
% Other PFTs skip (use_cesaraccio_chill absent) and keep the Sykes GDD threshold.
if isfield(pftParams, 'use_cesaraccio_chill') && pftParams.use_cesaraccio_chill == 1 ...
        && state.dormant_flag
    Tmin_d = forcing.T_avg;  Tmax_d = forcing.T_avg;
    if isfield(forcing, 'Tmin'), Tmin_d = forcing.Tmin; end
    if isfield(forcing, 'Tmax'), Tmax_d = forcing.Tmax; end
    [Cd, Ca] = vegetation.integrated.phenology.cesaraccioChillDay( ...
        Tmin_d, Tmax_d, forcing.T_avg, pftParams.baset);
    if state.chill_day > pftParams.crequ
        state.chill_day = state.chill_day + Cd;
    else
        state.anti_chill_day = state.anti_chill_day + Ca;
    end
end

% --- 4. Sykes 1996 chilling-GDD threshold (growth.cpp:92-93) ----------------
% gdd0 = k_chilla + k_chillb × exp(-k_chillk × chilldays)
% In Sykes 1996, the GDD compared against gdd0 is the spring GDD accumulator.
% We now use CLM5-style onset_gdd (which is correctly gated by ws_flag) instead
% of the all-year gdd5. This is the P1 fix per Phase D audit 2026-05-26.
gdd0_threshold = pftParams.k_chilla ...
               + pftParams.k_chillb * exp(-pftParams.k_chillk * state.chilldays);

% --- 5. Budburst trigger (Sykes 1996 + CLM5 onset_counter init) -------------
% [P3 fix 2026-05-26] CLM5 CNPhenologyMod.F90:911-915 — when threshold passes,
% set onset_flag=1 and initialize onset_counter to ndays_on; phen ramps up
% linearly over ndays_on days (not instant jump).
% [13b] Apple bud-burst via Cesaraccio chilling (post-chilling warmth anti_chill_day
% overcomes the chill deficit -chill_day); other PFTs via the Sykes GDD threshold
% (CLM5-FruitTree CNPhenologyMod.F90:1946).
if isfield(pftParams, 'use_cesaraccio_chill') && pftParams.use_cesaraccio_chill == 1
    budburst_ready = (state.anti_chill_day > 0) && (state.anti_chill_day >= -state.chill_day);
else
    budburst_ready = (state.onset_gdd > gdd0_threshold);
end
if state.dormant_flag && state.onset_gddflag == 1 && budburst_ready ...
        && state.aphen < APHEN_MAX
    state.dormant_flag  = false;
    state.onset_flag    = true;
    state.onset_counter = pftParams.ndays_on;       % start countdown (days)
    state.budburst_doy  = doy;
    state.onset_gddflag = 0;                          % CLM5 line 914 deactivate
    state.onset_gdd     = 0;                          % CLM5 line 915 reset
end

% --- 6. CLM5 daylength senescence trigger (CNPhenologyMod.F90:957) ----------
% [P4 fix 2026-05-26] Senescence trigger initializes offset_counter to ndays_off;
% phen ramps down linearly over ndays_off days.
% Guard: only trigger when ws_flag=0 (autumn half-year, matches CLM5 line 957),
% dormant=false, aphen sufficient, not already in offset.
if state.ws_flag == 0 ...
        && forcing.dayl_h < pftParams.crit_dayl_h ...
        && ~state.dormant_flag ...
        && state.aphen > APHEN_GROW_MIN ...
        && ~state.offset_flag
    state.offset_flag    = true;
    state.offset_counter = pftParams.ndays_off;       % start countdown
    state.senescence_doy = doy;
end

% --- 7. Phen ramp logic (CLM5 onset/offset_counter, lines 833-845, 808-820) -
% [P3+P4 fix 2026-05-26] Replace instantaneous 0→1 / 1→0 jumps with linear
% ramps over ndays_on / ndays_off, matching CLM5 storage→xfer→leaf timing.
% Daily increment: 1/ndays_on per day for onset, -1/ndays_off per day for offset.
if state.onset_flag
    % Linear ramp up to 1.0
    state.phen = min(1, state.phen + 1/pftParams.ndays_on);
    state.onset_counter = state.onset_counter - 1;
    if state.onset_counter <= 0
        % Onset complete (CLM5 line 840-845)
        state.onset_flag    = false;
        state.onset_counter = 0;
        state.phen          = 1;     % clamp to 1.0
    end
elseif state.offset_flag
    % Linear ramp down to 0.0
    state.phen = max(0, state.phen - 1/pftParams.ndays_off);
    state.offset_counter = state.offset_counter - 1;
    if state.offset_counter <= 0
        % Offset complete (CLM5 line 815-820): plant fully dormant
        state.offset_flag    = false;
        state.offset_counter = 0;
        state.phen           = 0;
        state.dormant_flag   = true;  % CLM5 line 822: re-enter dormant
    end
elseif ~state.dormant_flag
    % Steady-state full leaf (between onset complete and offset trigger)
    state.phen = 1;
else
    % Dormant, no leaves
    state.phen = 0;
end

% Diagnostic: ramp_progress for output (track ramp completion %)
if state.onset_flag
    ramp_progress = 1 - state.onset_counter / pftParams.ndays_on;
elseif state.offset_flag
    ramp_progress = state.offset_counter / pftParams.ndays_off;
else
    ramp_progress = NaN;
end

% --- 7. Accumulate aphen (growth.cpp:152) -----------------------------------
state.aphen = state.aphen + state.phen;

% APHEN_MAX hit: force phen=0 — a fallback ONLY for PFTs without a daylength
% senescence trigger. PFTs with crit_dayl_h drop leaves via offset_flag (step 6,
% CLM5 CNPhenologyMod.F90:957), so APHEN_MAX would prematurely kill leaves before
% the daylength offset fires. For apple (budburst DOY79 + 210 = DOY289), the
% APHEN_MAX cap zeroed LAI ~DOY289, BEFORE daylength offset (~DOY295-300) and
% cutting off the post-harvest leaf-on photosynthesis window that CLM5-FruitTree
% uses to refill C_storage (storage_flag=1 -> fcur=0, NutrientCompetitionFlexibleCNMod.F90:552-553),
% starving the overwinter reserve and killing year-2. Skip the cap when a
% daylength offset mechanism exists (LPJ-GUESS APHEN_MAX growth.cpp:50 is a
% fallback for PFTs lacking daylength control).
has_dayl_offset = isfield(pftParams, 'crit_dayl_h') && ~isnan(pftParams.crit_dayl_h) ...
                  && pftParams.crit_dayl_h > 0;
if state.aphen >= APHEN_MAX && ~has_dayl_offset
    state.phen = 0;
end

% --- 8. Fruit sub-state machine (Apple has_fruit PFTs, P5 fix 2026-05-26) ---
% Independent fruit phenology clock: budburst → bloom → fruit ripening → harvest
% Reference:
%   - DAISY phenology_TSum.C: DS ∈ [0, 2] for crop phenology, anthesis at DS=1
%   - LPJ-GUESS cropphenology.cpp:274-309: DS sub-phase computation
%   - Sitch 2003 Eq 27: Croot,new = Croot,old × (1 - froot) — separate organ turnover
%   - Apple Wang 2009 LMWE phenology dates
% Logic only fires if GDD_bloom_after_budburst is finite (i.e., apple-class PFT).
if isfield(pftParams, 'GDD_bloom_after_budburst') ...
        && ~isnan(pftParams.GDD_bloom_after_budburst)
    % Cache gdd5 at budburst (one-time on bud break day)
    if ~isnan(state.budburst_doy) && doy == state.budburst_doy && isnan(state.gdd5_at_budburst)
        state.gdd5_at_budburst = state.gdd5;
    end

    % Detect bloom: GDD since budburst > GDD_bloom_after_budburst
    if ~isnan(state.gdd5_at_budburst) && isnan(state.bloom_doy)
        GDD_since_budburst = state.gdd5 - state.gdd5_at_budburst;
        if GDD_since_budburst > pftParams.GDD_bloom_after_budburst
            state.bloom_doy     = doy;
            state.gdd5_at_bloom = state.gdd5;
        end
    end

    % DS_fruit ramp: post-bloom progression (0 at bloom → 1 at fruit maturity)
    if ~isnan(state.gdd5_at_bloom)
        GDD_since_bloom = state.gdd5 - state.gdd5_at_bloom;
        state.DS_fruit  = min(1, max(0, GDD_since_bloom / pftParams.TSum_fruit));
    else
        state.DS_fruit = 0;
    end

    % Harvest event: DS_fruit > harvest_DS, fire once
    if state.DS_fruit > pftParams.harvest_DS && ~state.harvested_flag
        state.harvested_flag = true;
        state.harvest_doy    = doy;
    end

    % Reset fruit state at COLDEST_DAY (annual cycle reset, parallel to aphen)
    if doy == coldest_day
        state.gdd5_at_budburst = NaN;
        state.bloom_doy        = NaN;
        state.gdd5_at_bloom    = NaN;
        state.DS_fruit         = 0;
        state.harvested_flag   = false;
        state.harvest_doy      = NaN;
    end
end

% --- 9. Cache today's dayl for next-day ws_flag computation (CLM5 line 802) ---
state.prev_dayl = forcing.dayl_h;

% --- Diagnostics (only-write) ----------------------------------------------
diag.gdd0_threshold  = gdd0_threshold;
diag.ws_flag         = state.ws_flag;
diag.onset_gddflag   = state.onset_gddflag;
diag.onset_gdd       = state.onset_gdd;
diag.DS_fruit        = state.DS_fruit;
diag.bloom_doy       = state.bloom_doy;
diag.harvest_doy     = state.harvest_doy;
diag.ramp_progress   = ramp_progress;
diag.reset_triggered = reset_triggered;
diag.phen            = state.phen;
diag.gdd5            = state.gdd5;
diag.chilldays       = state.chilldays;
diag.aphen           = state.aphen;

end
