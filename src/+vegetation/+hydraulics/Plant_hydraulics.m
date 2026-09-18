function [Tp_sun, Tp_sha, psil_sun, psil_sha, psir, psis, wt_sun, wt_sha, gsr] = ...
    Plant_hydraulics(plantHydro, Khh, hh, Tp_sun, Tp_sha, KT, H)
%PLANT_HYDRAULICS Main plant hydraulics wrapper.
%   Calls calSPAC and computes hydraulic stress factors for sunlit and
%   shaded leaves.
%
%   Inputs:
%     plantHydro - plant hydraulic initialization/state structure.
%     Khh        - soil hydraulic conductivity at endpoints (NL x 2, cm/s).
%     hh         - soil pressure head (NN x 1, cmH2O).
%     Tp_sun     - sunlit transpiration demand time series (cm/s).
%     Tp_sha     - shaded transpiration demand time series (cm/s).
%     KT         - current time-step index.
%     H          - canopy height (m).
%
%   Outputs:
%     Tp_sun/Tp_sha     - hydraulically limited transpiration (cm/s).
%     psil_sun/psil_sha - leaf water potentials (MPa).
%     psir/psis         - root and weighted soil water potentials (MPa).
%     wt_sun/wt_sha     - hydraulic stress factors in [0,1].
%     gsr               - nodal soil-root conductance (m/s/MPa).

    [Tp_sun, Tp_sha, psil_sun, psil_sha, psir, psis, gsr] = ...
        vegetation.hydraulics.calSPAC(plantHydro, Khh, hh, Tp_sun, Tp_sha, KT, H);

    wt_sun = vegetation.hydraulics.VulnerabilityCurve( ...
                 psil_sun(KT), plantHydro.gpmax, plantHydro.p50, plantHydro.aa) ...
             / plantHydro.gpmax;

    wt_sha = vegetation.hydraulics.VulnerabilityCurve( ...
                 psil_sha(KT), plantHydro.gpmax, plantHydro.p50, plantHydro.aa) ...
             / plantHydro.gpmax;

    wt_sun = max(0, min(1, wt_sun));
    wt_sha = max(0, min(1, wt_sha));
end
