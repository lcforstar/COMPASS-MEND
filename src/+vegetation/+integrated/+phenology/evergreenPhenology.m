function [state, diag] = evergreenPhenology(state, pftParams, forcing, doy)
%EVERGREENPHENOLOGY  Evergreen phenology — constant phen=1, no budburst event.
%
%   For PFTs with continuous leaf turnover (no winter dormancy):
%   - Needleleaf evergreen (PFT 1, 2)
%   - Broadleaf evergreen (PFT 4, 5)
%   - Broadleaf evergreen shrub (PFT 9)
%
%   Returns:
%     state    [inout]  PhenologyState (phen forced to 1, no event flags)
%     diag     struct   diagnostics
%
%   Method (CLM5 CNEvergreenPhenology simplified, CNPhenologyMod.F90:502-666):
%     phen = 1.0 constant (caller drives LAI via leaf_long parameter,
%     not via phen).
%
%   For evergreen PFTs, LAI dynamics come from continuous leaf turnover
%   (bglfr = 1 / (leaf_long × 365)), not from phenology onset/offset events.
%   The phen variable here is just a marker (=1, leaves always present);
%   COMPASS LAI module should multiply phen × LAImax — for evergreens,
%   this means LAI ≡ LAImax (always full canopy, modulated only by stress).
%
%   References:
%     CLM5 CNPhenologyMod.F90:502-666 (CNEvergreenPhenology)
%     CLM5 CNPhenologyMod.F90:603 (evergreen flag gate, constant bglfr)

state.phen          = 1.0;
state.dormant_flag  = false;
state.onset_flag    = false;
state.offset_flag   = false;

% First-call initialization (budburst_doy NaN → set to 1, day-of-year start)
if isnan(state.budburst_doy)
    state.budburst_doy = 1;
end

% aphen tracks leaf-on days; evergreen always leaf-on, accumulates linearly
state.aphen = state.aphen + 1;

diag.phen          = state.phen;
diag.evergreen     = true;

end
