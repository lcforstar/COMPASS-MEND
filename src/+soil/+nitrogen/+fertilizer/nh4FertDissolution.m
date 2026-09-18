function [Sdis, DeltaGranule] = nh4FertDissolution(K_nh4_dissolve, NL, NH4_granule, WFPS, Delt_t, TT)
%NH4FERTDISSOLUTION  Granular NH4-fertilizer dissolution pool — first-order release to NH4.
%
% Physical model:
%     Surface granular NH4 fertilizer (e.g. (NH4)2SO4, NH4Cl pellets) does not
%     instantaneously dissolve into the bulk soil NH4 pool. We track an extra
%     transient state (granule pool) and let it leak into the dissolved NH4 pool
%     by first-order kinetics with Q10 temperature and WFPS moisture modulation,
%     mirroring the urea-hydrolysis convention of APSIM-SoilN (Probert 1998).
%
%     This addresses the well-known mean-field underestimation of NH3 volatilization
%     (Sommer & Hutchings 2001 EJA): when fertilizer is treated as instantaneous
%     bulk mixing, the surface NH4 concentration peak is severely diluted, and the
%     Henry-driven NH3 efflux Svol is underestimated 10-30× compared to field
%     observations.
%
% Formulation (analytically integrated, identical structure to ureaHydrolysis.m):
%     k_eff(ML,ND) = 5 × WFPS(ML,ND) × K_nh4_dissolve_T(ML,ND)        [1/s]
%     K_nh4_dissolve_T = K_nh4_dissolve × Q10^((T − T_ref)/10)         [1/s]
%     ΔG = NH4_granule(MN) × (1 − exp(−k_eff × Δt))                    [μg/cm³]
%     Sdis(ML,ND) = ΔG / Δt                                             [μg/(cm³·s)]
%
%     ADE coupling (assembleSoluteParameters.m):
%         G1_NH4 += +Sdis   (granule → dissolved NH4 pool, source term)
%
% Inputs:
%     K_nh4_dissolve [in]    [1/s, NL×1]      base dissolution rate (Excel 'K_nh4_dissolve' /86400)
%     NL             [in]    [-]              number of layers
%     NH4_granule    [in]    [μg/cm³, NN×1]  prior-step granule-pool concentration at nodes
%     WFPS           [in]    [-, NL×2]        water-filled pore space
%     Delt_t         [in]    [s]              outer step length
%     TT             [in]    [°C, NN×1]      soil temperature at nodes
%
% Returns:
%     Sdis         [μg/(cm³·s), NL×2]  dissolution rate to be summed into ADE G1 NH4 source
%     DeltaGranule [μg/cm³, NN×1]      granule concentration drop this Δt (for caller pool update)
%
% Reference:
%     Probert M E et al. 1998 Agric Systems 56(1):1-28 — APSIM-SoilN urea-hydrolysis-style
%         transient surface pool, here generalized to NH4 granular fertilizers
%     Sommer S G & Hutchings N J 2001 European Journal of Agronomy 15(1):1-15 —
%         identifies mean-field fertilizer mixing as primary cause of model
%         underestimation of NH3 volatilization
%     Lloyd J & Taylor J A 1994 Func Ecol 8:315-323 — Q10 = 2.5 temperature response
%     Sherlock R R & Goh K M 1985 Fertilizer Research 6:3-22 — surface NH3 volatilization
%         physics requiring localized peak NH4 concentration
%
% History (2026-05-17):
%     New module. Mirrors src/+soil/+nitrogen/+chemistry/ureaHydrolysis.m exactly
%     in numerical structure (analytical 1st-order with WFPS×Q10 modulation), but
%     applies to granular NH4 fertilizers rather than urea. NO3 granules are NOT
%     handled here — they bypass dissolution pool (no volatilization driver, no
%     adsorption, so instantaneous mixing is physically correct).

    Q10 = 2.5;
    T_ref = 20;

    DeltaGranule = zeros(size(NH4_granule));
    Sdis = zeros(NL, 2);

    for ML = 1:NL
        for ND = 1:2
            MN = ML + ND - 1;

            % Temperature-adjusted dissolution rate (Lloyd & Taylor 1994 Q10)
            K_dis_T = K_nh4_dissolve(ML) * Q10^((TT(MN) - T_ref) / 10);     % [1/s]

            % Effective rate with WFPS modulation (mirror ureaHydrolysis factor 5.0)
            k_eff = 5.0 * WFPS(ML, ND) * K_dis_T;                            % [1/s]

            % Analytical concentration drop over Δt
            dG = NH4_granule(MN) * (1 - exp(-k_eff * Delt_t));               % [μg/cm³]

            % Average rate over Δt (for ADE G1 NH4 source, μg/(cm³·s))
            Sdis(ML, ND) = dG / max(Delt_t, eps);

            % Concentration drop for caller's granule-pool update
            DeltaGranule(MN) = dG;
        end
    end
end
