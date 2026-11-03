function [state, diag] = cropPhenology(state, pftParams, forcing, doy)
%CROPPHENOLOGY  Annual crop phenology placeholder (DS-based, DAISY/LPJ-GUESS style).
%
%   For one-season-cycle crops (PFT 15-22: corn/wheat/soybean/cotton/rice/...).
%
%   PLACEHOLDER: this is a first-pass implementation. COMPASS current
%   focus is Apple PFT 23 (season_decid), so this branch is not exercised
%   in the Apple-only Changwu site. Full implementation deferred to future
%   crop-PFT site validation.
%
%   Returns:
%     state    [inout]  PhenologyState struct
%     diag     struct   diagnostics
%
%   Inputs:
%     state         [inout]  PhenologyState (per-PFT)
%     pftParams     [in]     PFT params (Sykes 1996 fields should be NaN)
%     forcing       [in]     struct: .T_avg [°C], .dayl_h [h], .mtemp [°C]
%     doy           [in]     day-of-year
%
%   Method (LPJ-GUESS cropphenology.cpp:274-309 DS scheme, simplified):
%     DS = ΣTa / TSum_crop, range [0, 2]
%     DS<1 vegetative, DS=1 anthesis, DS=2 maturity → harvest
%     phen ramps from 0 to 1 during DS ∈ [0, 0.3], stays 1 to DS=1.5,
%     ramps down to 0 by DS=2.0.
%
%     Sowing trigger: sowing_day from config_vegetation.txt (fixed DOY)
%     or LPJ-GUESS 20-year mean detection (cropsowing.cpp:28-71).
%
%   References:
%     LPJ-GUESS modules/cropphenology.cpp:274-309 (DS computation)
%     CLM5 CNPhenologyMod.F90:1420-2031 (CropPhenology natural-veg branch)

if isfield(pftParams, 'winter_wheat_cn') && pftParams.winter_wheat_cn == 1
    [state, diag] = vegetation.integrated.phenology.winterWheatPhenology( ...
        state, pftParams, forcing, doy);
    return;
end

% --- Hard constants (placeholder defaults for crop PFTs not yet implemented) ---
TSum_crop  = 1500;   % [°C·d] generic temperate cereal default
T_BASE     = 5;      % [°C]
SOWING_DOY = 100;    % default; should come from config

% --- Accumulate GDD (DAISY DS style: ΣTa / TSum) ---
if doy >= SOWING_DOY && state.dormant_flag
    state.dormant_flag = false;
    state.onset_flag   = true;
    state.budburst_doy = doy;
end

if ~state.dormant_flag
    state.gdd5 = state.gdd5 + max(0, forcing.T_avg - T_BASE);
    DS = state.gdd5 / TSum_crop;   % [0, 2]
else
    DS = 0;
end

% --- Phen mapped from DS ---
if DS < 0.3
    state.phen = DS / 0.3;   % ramp to 1 by DS=0.3
elseif DS < 1.5
    state.phen = 1;
elseif DS < 2.0
    state.phen = max(0, 1 - (DS - 1.5) / 0.5);
else
    state.phen = 0;
    if ~state.offset_flag
        state.offset_flag    = true;
        state.senescence_doy = doy;   % maturity / harvest
    end
end

% --- aphen ---
state.aphen = state.aphen + state.phen;

% --- Diagnostics ---
diag.phen   = state.phen;
diag.DS     = DS;
diag.gdd5   = state.gdd5;

end
