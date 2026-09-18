function [state, diag] = phenologyDispatch(state, pftParams, forcing, doy, year_mtemp_max)
%PHENOLOGYDISPATCH  CLM5-style filter+flag-gate dispatcher for per-PFT phenology.
%
%   Routes a single PFT to one of 4 phenology branches based on its
%   binary flags. Exactly one flag must be 1, others 0 (validated at PFT-table
%   load time; here we trust the input).
%
%   Returns:
%     state    [inout]  PhenologyState struct (per-PFT, persists across days)
%     diag     struct   per-branch diagnostics from the active branch
%
%   Inputs:
%     state         [inout]  PhenologyState
%     pftParams     [in]     struct with .evergreen / .season_decid /
%                            .stress_decid / .crop_phenology binary flags,
%                            plus branch-specific params (Sykes 1996, etc.)
%     forcing       [in]     struct: .T_avg [°C], .dayl_h [h], .mtemp [°C],
%                              optionally .water_stress [-]
%     doy           [in]     day-of-year (1-365)
%     year_mtemp_max [in]    [°C] max monthly temp this year (used by season_decid)
%
%   Architecture: CLM5 CNPhenologyMod.F90:241-342 (CNPhenology master driver)
%     calls 4 sibling branches in sequence; each branch gates by its PFT flag.
%     COMPASS site-level: dispatcher is per-PFT, branches once per call.
%
%   References:
%     CLM5 CNPhenologyMod.F90:241-342 (master CNPhenology)
%     CLM5 pftconMod.F90:200-202 (binary flag declarations)
%     LPJ-GUESS framework/guess.h:61 (phenologytype enum)

% --- Sanity check: exactly one flag = 1 ---
flag_sum = pftParams.evergreen + pftParams.season_decid ...
         + pftParams.stress_decid + pftParams.crop_phenology;
if flag_sum ~= 1
    error('phenologyDispatch:flagSum', ...
        ['Phenology flags must sum to 1, got %d. ', ...
         'Check PFT table: evergreen + season_decid + stress_decid + crop_phenology = 1.'], ...
        flag_sum);
end

% --- Default water_stress = 0 if not provided ---
if ~isfield(forcing, 'water_stress')
    forcing.water_stress = 0;
end

% --- Dispatch to the active branch (CLM5 flag gate pattern) ---
if pftParams.evergreen == 1
    [state, diag] = vegetation.integrated.phenology.evergreenPhenology( ...
        state, pftParams, forcing, doy);
    diag.branch = 'evergreen';
elseif pftParams.season_decid == 1
    [state, diag] = vegetation.integrated.phenology.seasonDecidPhenology( ...
        state, pftParams, forcing, doy, year_mtemp_max);
    diag.branch = 'season_decid';
elseif pftParams.stress_decid == 1
    [state, diag] = vegetation.integrated.phenology.stressDecidPhenology( ...
        state, pftParams, forcing, doy);
    diag.branch = 'stress_decid';
elseif pftParams.crop_phenology == 1
    [state, diag] = vegetation.integrated.phenology.cropPhenology( ...
        state, pftParams, forcing, doy);
    diag.branch = 'crop_phenology';
end

end
