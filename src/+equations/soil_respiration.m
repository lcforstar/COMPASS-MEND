function Rh = soil_respiration(T_air, P30_cm)
%SOIL_RESPIRATION  Heterotrophic soil respiration, Raich, Potter & Bhagawati (2002).
%
%   Rh = soil_respiration(T_air, P30_cm)
%
%   Returns:
%     Rh        [kg C/m2/s]  heterotrophic (microbial) soil respiration rate
%
%   Inputs:
%     T_air     [degC]  air / near-surface temperature, this step
%     P30_cm    [cm]    running 30-day precipitation total (the moisture term)
%
%   Method:
%     Raich, Potter & Bhagawati (2002) Glob. Change Biol. 8:800-812,
%     doi:10.1046/j.1365-2486.2002.00511.x  (the temperature-&-precipitation "T&P" model):
%       Rs = F * exp(Q*T) * P/(K+P)        [g C/m2/d]   total soil respiration
%       Rh = f_h * Rs                                    heterotrophic fraction
%     F   [g C/m2/d] base rate -- global default 1.25 (their Table; site-calibrate so annual
%         Rh ~ litterfall input, ~190 g C/m2 here); Q [1/degC] T sensitivity 0.05452;
%         K [cm] precip half-saturation 4.259; f_h [-] Rh/Rs = 0.5 (Bond-Lamberth et al. 2004
%         Glob. Change Biol. 10:1756). Directly computable from forcing T + precip -- NO SOM
%         pool cascade, NO N cycle (so it also works for C-only runs).
%
%   The 1-arg form soil_respiration(Ts) now ERRORS [2026-06-14] (it used to return 0 silently,
%   which gave standalone ebal / profile / old tests a fake 0 soil respiration). The ecosystem Rh
%   for NEE is computed in COMPASS.m, where the running 30-day precip window is available, and
%   fluxes.Resp is set there; ebal*.m set fluxes.Resp = 0 as an explicit placeholder (overwritten
%   in the main loop). This REPLACES the old proxy 2.5+0.054*Ts (which gave Rh ~1177 -> NEE wrong
%   sign). DAMM (Davidson et al. 2012) moisture form is the documented upgrade path.

    if nargin < 2
        % [2026-06-14] The 1-arg legacy form used to return 0 SILENTLY, which gave standalone ebal /
        %   profile / old tests a fake 0 soil respiration. Raich needs (T_air, P30_cm); the ecosystem
        %   Rh for NEE is computed in COMPASS.m (2-arg). Error loudly so no caller silently gets 0.
        error('soil_respiration:needsP30', ...
              'soil_respiration requires (T_air, P30_cm). The 1-arg form is removed; Rh is computed in COMPASS.m via the 2-arg Raich form.');
    end

    % [SOURCE: Raich, Potter & Bhagawati 2002 GCB 8:800 globals]  Rh/Rs=0.5: Bond-Lamberth 2004 GCB 10:1756
    % [F reverted 2.10->1.25 2026-06-14] Restored the Raich 2002 published GLOBAL default F=1.25. The 2.10
    %   value (2026-06-13) was site-calibrated to the Italian CLM5-FruitTree SIMULATED Rh~444 gC/m2/yr
    %   (Dombrowski et al. 2022 GMD 15:5167 p.5179, humid Caldaro site); for the Loess Plateau apple orchard
    %   it overshot -- total soil respiration Rs = 2*Rh would be ~884 gC/m2/yr, vs the Loess-region
    %   expectation. At F=1.25 -> Rs~526 / Rh~263, consistent with the Loess apple-orchard ecosystem
    %   respiration benchmark Reco=525 (Yang et al. 2018 Huan Jing Ke Xue 39(5):2339, eddy covariance,
    %   PMID 29965535). [SOURCE UNVERIFIED 2026-06-14] chamber soil-Rs ~450-600 (Changwu, hjkx 2014) cited
    %   from search abstract only; full text not retrieved -> not used as a hard target.
    F = 1.25;  Q = 0.05452;  K = 4.259;  f_h = 0.5;
    Rs_gC_m2_d = F * exp(Q * T_air) * P30_cm / (K + P30_cm);   % [g C/m2/d] total soil respiration
    Rh = f_h * Rs_gC_m2_d / 86400 / 1000;                     % [g C/m2/d] -> [kg C/m2/s] heterotrophic
end
