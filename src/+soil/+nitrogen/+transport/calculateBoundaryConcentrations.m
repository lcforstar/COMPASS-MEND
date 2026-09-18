function [cTop_NH4_new, cTop_NO3_new, cBot_NH4_new, cBot_NO3_new] = calculateBoundaryConcentrations( ...
        NitrogenFlux, NitrogenBoundary, QL_nodes, KT, Delt_t) %#ok<INUSD>
%CALCULATEBOUNDARYCONCENTRATIONS  Compute top + bottom N boundary concentrations at step KT.
%
% Top boundary models atmospheric N deposition (wet/dry) carried in by rainfall,
% mixing with any surface ponding reservoir. Bottom boundary uses prescribed
% groundwater / drainage N concentration (user-supplied time series).
%
% Inputs:
%   NitrogenFlux       [in]   struct   current N flux time-series; we read
%                                       previous-step cTop_*/cBot_* and only
%                                       update the KT-th element.
%   NitrogenBoundary   [in]   struct   per-step inputs:
%       .cT_dep_flux  [mg/m²/h, tS×1]  NH4 deposition FLUX (Excel col 22)
%       .cT1_dep_flux [mg/m²/h, tS×1]  NO3 deposition FLUX (Excel col 23)
%       .cB_NH4_conc  [mg/L,    tS×1]  prescribed bottom NH4 CONCENTRATION
%       .cB_NO3_conc  [mg/L,    tS×1]  prescribed bottom NO3 CONCENTRATION
%       .Precip_cms   [cm/s,    tS×1]  rainfall rate (used for unit conversion)
%       .DSTMAX, DSTOR0, h_target, NBCh, NBCN  (ponding/BC switches)
%       .Infiltration [cm/s, tS×1]  actual infiltration
%       .EVAP         [cm/s, tS×1]  evaporation
%   QL_nodes           [in]   [cm/s, NN×1]  liquid flux at nodes (reserved)
%   KT                 [in]   [-]      current step index
%   Delt_t             [in]   [s]      step length
%
% Returns:
%   cTop_NH4_new / cTop_NO3_new   [mg/L, tS×1] only KT element updated
%   cBot_NH4_new / cBot_NO3_new   [mg/L, tS×1] only KT element updated
%
% Method (KT-step computation):
%   1. Convert deposition FLUX → CONCENTRATION via rainfall water flux:
%        cT_NH4_conc = NH4_dep_flux / (Precip[cm/s] × 3.6e4)   [mg/L]
%      → If Precip < 1e-7 cm/s (≈ 0.0036 mm/h, ie no rain): cT_NH4_conc = 0
%        (Dry deposition not handled here; deferred to a future surface flux BC.)
%   2. If surface ponding (DSTMAX>0 and h_target>0): mix ponding reservoir + rain:
%        cTop_NH4 = (h_target/Δt × cPond + Infil × cT_conc) / (h_target/Δt + Infil − E)
%   3. Else atmospheric water BC + flux solute BC (NBCh=3, NBCN=2): scale by net flux.
%   4. Bottom: cBot_NH4(KT) = cB_NH4_conc(KT).
%
% Reference:
%   STEMMUS-N legacy N_sub.m (Zhu Jizixing thesis Eq 3.41-3.43)
%   van Genuchten & Wagenet 1989 SSSAJ — solute boundary mixing for ADE
%
% Refactored 2026-05-13:
%   - removed redundant intermediate vectors (only KT scalar needed)
%   - explicit FLUX→CONCENTRATION unit conversion via rainfall rate
%   - field names: cT/cT1/cB/cB1 → cT_dep_flux/cT1_dep_flux/cB_NH4_conc/cB_NO3_conc

% --- 1. Copy current full time series; we'll only overwrite the KT-th element ---
cTop_NH4_new = NitrogenFlux.cTop_NH4;
cTop_NO3_new = NitrogenFlux.cTop_NO3;
cBot_NH4_new = NitrogenFlux.cBot_NH4;
cBot_NO3_new = NitrogenFlux.cBot_NO3;

% --- 2. Convert deposition FLUX [mg/(m²·h)] to in-rain CONCENTRATION [mg/L] ---
%   Water volume per m² per hour = Precip[cm/s] × 3600 s × 0.01 m/cm × 1000 L/m³
%                                = Precip × 3.6e4   L/(m²·h)
%   cT_conc = dep_flux / water_volume_per_m2_per_h
Precip_cms_now   = NitrogenBoundary.Precip_cms(KT);   % cm/s
NH4_dep_flux_now = NitrogenBoundary.cT_dep_flux(KT);  % mg/(m²·h)
NO3_dep_flux_now = NitrogenBoundary.cT1_dep_flux(KT);
if Precip_cms_now > 1e-7  % >= 0.0036 mm/h
    water_Lpm2h  = Precip_cms_now * 3.6e4;  % L/(m²·h)
    cT_NH4_conc_KT = NH4_dep_flux_now / water_Lpm2h;  % mg/L
    cT_NO3_conc_KT = NO3_dep_flux_now / water_Lpm2h;
else
    % Insufficient rainfall: dry deposition not carried by rain → no top BC contribution.
    % (Dry deposition needs a separate surface-flux source term, future work.)
    cT_NH4_conc_KT = 0;
    cT_NO3_conc_KT = 0;
end

% --- 3. Top BC: zero default ---
cTop_NH4_new(KT) = 0;
cTop_NO3_new(KT) = 0;

if NitrogenBoundary.DSTMAX > 0 && NitrogenBoundary.h_target > 0
    % --- Case A: surface ponding active → mix ponding + rain ---
    % Pond reservoir default concentration = rainfall concentration (no separate tracking)
    cPond_NH4 = cT_NH4_conc_KT;
    cPond_NO3 = cT_NO3_conc_KT;
    cTop_NH4_new(KT) = cPond_NH4;
    cTop_NO3_new(KT) = cPond_NO3;
    AVAIL0 = (NitrogenBoundary.Infiltration(KT) - NitrogenBoundary.EVAP(KT)) ...
             + NitrogenBoundary.DSTOR0 / Delt_t;
    if AVAIL0 > 0
        denom_mix = NitrogenBoundary.h_target / Delt_t ...
                    + (NitrogenBoundary.Infiltration(KT) - NitrogenBoundary.EVAP(KT));
        cTop_NH4_new(KT) = (NitrogenBoundary.h_target / Delt_t * cPond_NH4 ...
                            + NitrogenBoundary.Infiltration(KT) * cT_NH4_conc_KT) / denom_mix;
        cTop_NO3_new(KT) = (NitrogenBoundary.h_target / Delt_t * cPond_NO3 ...
                            + NitrogenBoundary.Infiltration(KT) * cT_NO3_conc_KT) / denom_mix;
    end
elseif NitrogenBoundary.NBCh == 3 && NitrogenBoundary.NBCN == 2
    % --- Case B: atmospheric water BC + flux solute BC → scale by net flux ---
    net_flux = NitrogenBoundary.Infiltration(KT) - NitrogenBoundary.EVAP(KT);
    if abs(net_flux) > 1e-6
        cTop_NH4_new(KT) = NitrogenBoundary.Infiltration(KT) / net_flux * cT_NH4_conc_KT;
        cTop_NO3_new(KT) = NitrogenBoundary.Infiltration(KT) / net_flux * cT_NO3_conc_KT;
    end
    % else: net_flux ≈ 0 (pure evaporation balance), no top BC update — keep 0
end

% --- 4. Bottom BC: prescribed Dirichlet from user (groundwater conc.) ---
cBot_NH4_new(KT) = NitrogenBoundary.cB_NH4_conc(KT);
cBot_NO3_new(KT) = NitrogenBoundary.cB_NO3_conc(KT);

end
