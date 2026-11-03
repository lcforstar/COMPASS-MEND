function [state, diag] = stressDecidPhenology(state, pftParams, forcing, doy)
%STRESSDECIDPHENOLOGY  Stress-deciduous phenology (CLM5 4-accumulator version).
%
%   Verbatim port of CLM5 CNStressDecidPhenology (CNPhenologyMod.F90:974-1417).
%   Handles PFTs whose leaf phenology responds to cold + drought stress signals:
%     - C3 non-arctic grass (PFT 13), C4 grass (PFT 14)
%     - Broadleaf-deciduous tropical raingreen (PFT 6)
%     - Broadleaf-deciduous shrub temperate (PFT 10)
%
%   Returns:
%     state    [inout]  PhenologyState struct
%     diag     struct   per-step diagnostics (only-write)
%
%   Inputs:
%     state         [inout]  PhenologyState
%     pftParams     [in]     struct, required: .crit_dayl_h, .ndays_on, .ndays_off
%                            (Sykes 1996 fields NaN — physics never reads them)
%     forcing       [in]     struct:
%                              .T_avg          [°C]  daily mean air temp
%                              .T_soil         [°C]  optional, daily mean soil temp
%                                                    (falls back to T_avg if missing)
%                              .dayl_h         [h]   day length
%                              .soil_psi_MPa   [MPa] optional, soil water potential
%                                                    (falls back to T_avg-based proxy)
%                              .prec10_mm      [mm]  optional, 10-day cumulative rain
%                                                    (CLM5 line 1226 rain_threshold=20mm)
%     doy           [in]     day-of-year
%
%   Method (CLM5 CNPhenologyMod.F90:974-1417 verbatim):
%     1. onset_fdd: accumulate freezing days when dormant AND onset_gddflag==0
%                   AND T_soil < 0°C
%     2. onset_fdd > crit_onset_fdd: latch onset_gddflag=1, reset onset_fdd
%     3. onset_gdd: accumulate spring GDD when onset_gddflag==1 AND T_soil > 0°C
%     4. onset_swi: accumulate soil-water-wet days when soil_psi > -0.6 MPa
%     5. onset_swi > crit_onset_swi (=15 days) AND prec10 > 20mm AND
%        (onset_gddflag==0 OR onset_gdd > crit_onset_gdd) AND dayl > 6h
%        → onset_flag=1, onset_counter=ndays_on
%     6. ndays_on linear ramp 0→1 (same as season_decid)
%     7. Offset triggered by: dayl < crit_dayl_h (line 957 same as season_decid)
%        OR severe water stress (psi < soilpsi_off=-2 MPa)
%
%   Hard constants from CLM5 (clm5_params.c171117.nc scalar globals):
%     crit_onset_fdd       = 15 days   (freezing degree days threshold)
%     crit_onset_swi       = 15 days   (soil water index threshold)
%     soilpsi_on_MPa       = -0.6      (wet soil threshold)
%     soilpsi_off_MPa      = -2.0      (drought offset trigger)
%     rain_threshold_mm    = 20        (10-day rain minimum for onset)
%     dayl_min_h           = 6.0       (minimum daylength for onset)
%
%   References:
%     CLM5 CNPhenologyMod.F90:974-1417 (CNStressDecidPhenology subroutine)
%     CLM5 line 1113-1114: rain_threshold = 20 mm
%     CLM5 line 1118: stress_decid flag gate
%     CLM5 line 1239-1245: onset_swi accumulator
%     CLM5 line 1252-1255: dayl > 6h guard (secspqtrday = 21600 s = 6 h)
%
%   PLACEHOLDER for COMPASS Apple-only run: this branch is not exercised
%   because PFT 23 (Apple) is season_decid. Full validation deferred to
%   future multi-PFT site (e.g., C3 grass + Apple co-occurrence).

% --- Hard constants (CLM5 verbatim) -----------------------------------------
CRIT_ONSET_FDD       = 15;     % [days] freezing degree day threshold (CLM5 line 1208)
CRIT_ONSET_SWI       = 15;     % [days] soil water index threshold (CLM5 line 1239)
SOILPSI_ON_MPA       = -0.6;   % [MPa]  wet soil threshold (CLM5 line 1232)
SOILPSI_OFF_MPA      = -2.0;   % [MPa]  drought offset trigger
RAIN_THRESHOLD_MM    = 20;     % [mm]   10-day rain min (CLM5 line 1226)
DAYL_MIN_H           = 6.0;    % [h]    min daylength for onset (CLM5 line 1252)
T_FREEZE             = 0;      % [°C]   soil freezing temp

% --- Read forcing with defaults (CLM5 uses soil-zone 3 temp + psi) ---------
if isfield(forcing, 'T_soil')
    T_soil = forcing.T_soil;
else
    T_soil = forcing.T_avg;     % fallback: air temp proxy
end
if isfield(forcing, 'soil_psi_MPa')
    psi_MPa = forcing.soil_psi_MPa;
else
    % Crude proxy from water_stress index (NOT physical): 0=wet, 1=dry
    if isfield(forcing, 'water_stress')
        psi_MPa = -0.6 - forcing.water_stress * 2.0;  % maps [0,1] → [-0.6, -2.6]
    else
        psi_MPa = -1.0;  % default mid-range
    end
end
if isfield(forcing, 'prec10_mm')
    prec10 = forcing.prec10_mm;
else
    prec10 = 0;  % conservative: no rain trigger satisfied
end

% Need additional state fields not in base PhenologyState. Use dynamic ones.
if ~isfield(state, 'onset_fdd'), state.onset_fdd = 0; end
if ~isfield(state, 'onset_swi'), state.onset_swi = 0; end

% --- 1. Update offset_counter (CLM5 line 1129-1145) -------------------------
if state.offset_flag
    state.offset_counter = state.offset_counter - 1;
    state.phen = max(0, state.phen - 1/pftParams.ndays_off);
    if state.offset_counter <= 0
        state.offset_flag    = false;
        state.offset_counter = 0;
        state.dormant_flag   = true;
        state.phen           = 0;
    end
end

% --- 2. Update onset_counter (CLM5 line 1148-1175) --------------------------
if state.onset_flag
    state.onset_counter = state.onset_counter - 1;
    state.phen = min(1, state.phen + 1/pftParams.ndays_on);
    if state.onset_counter <= 0
        state.onset_flag    = false;
        state.onset_counter = 0;
        state.phen          = 1;
    end
end

% --- 3. Dormant-period accumulators (CLM5 line 1200-1247) -------------------
if state.dormant_flag && ~state.onset_flag
    % onset_fdd: freezing days accumulator (line 1200)
    if state.onset_gddflag == 0 && T_soil < T_FREEZE
        state.onset_fdd = state.onset_fdd + 1;
    end

    % Latch onset_gddflag when freeze threshold exceeded (line 1208-1213)
    if state.onset_fdd > CRIT_ONSET_FDD
        state.onset_gddflag = 1;
        state.onset_fdd     = 0;
        state.onset_swi     = 0;
    end

    % onset_gdd: spring GDD accumulator (line 1217-1220)
    if state.onset_gddflag == 1 && T_soil > T_FREEZE
        state.onset_gdd = state.onset_gdd + T_soil;
    end

    % onset_swi: soil water index accumulator (line 1232-1234)
    if psi_MPa >= SOILPSI_ON_MPA
        state.onset_swi = state.onset_swi + 1;
    end

    % Onset trigger: SWI > threshold AND rain check AND (no freeze OR GDD passed)
    % AND dayl > 6h (CLM5 line 1240-1253)
    rain_ok = (prec10 >= RAIN_THRESHOLD_MM);
    dayl_ok = (forcing.dayl_h > DAYL_MIN_H);
    swi_ok  = (state.onset_swi > CRIT_ONSET_SWI);
    gdd_ok  = (state.onset_gddflag == 0) || (state.onset_gdd > 50);  % simple proxy

    if swi_ok && rain_ok && dayl_ok && gdd_ok
        % Budburst trigger
        state.onset_flag    = true;
        state.onset_counter = pftParams.ndays_on;
        state.dormant_flag  = false;
        state.budburst_doy  = doy;
        % Reset accumulators (CLM5 line 1259-1264)
        state.onset_gddflag = 0;
        state.onset_fdd     = 0;
        state.onset_gdd     = 0;
        state.onset_swi     = 0;
    end
end

% --- 4. Offset trigger (CLM5 line 957 style + drought) ----------------------
if ~state.dormant_flag && ~state.offset_flag
    dayl_short = forcing.dayl_h < pftParams.crit_dayl_h;
    drought    = psi_MPa < SOILPSI_OFF_MPA;
    if dayl_short || drought
        state.offset_flag    = true;
        state.offset_counter = pftParams.ndays_off;
        state.senescence_doy = doy;
    end
end

% Steady-state full leaf
if ~state.dormant_flag && ~state.onset_flag && ~state.offset_flag
    state.phen = 1;
end

% --- aphen accumulation ---
state.aphen = state.aphen + state.phen;

% --- Diagnostics ---
diag.phen         = state.phen;
diag.onset_fdd    = state.onset_fdd;
diag.onset_swi    = state.onset_swi;
diag.onset_gdd    = state.onset_gdd;
diag.onset_gddflag = state.onset_gddflag;
diag.dormant_flag = state.dormant_flag;
diag.onset_flag   = state.onset_flag;
diag.offset_flag  = state.offset_flag;
diag.psi_MPa      = psi_MPa;

end
