function [NPP, R_a_plant, diag] = computeNPP(GPP_kg, Rd_leaf_kg, N_root_fine, N_livestem, N_livecroot, N_fruit, T_air, T_soil, pftParams)
%COMPUTENPP  Net primary production following Sitch et al. (2003) LPJ.
%
%   Maintenance respiration uses the Lloyd & Taylor (1994) modified Arrhenius
%   temperature response (Sitch 2003 Eq 23) rather than a fixed Q10. The
%   per-tissue formula R_organ = r × N_organ × g(T) is mathematically equivalent
%   to Sitch Eq 22 since N_organ ≡ C_organ / cn_organ. Leaf maintenance
%   respiration is taken directly from SCOPE biochemistry (Farquhar 1980) to
%   avoid double-counting with the canopy Rd already returned by ebal.
%
%   Returns:
%     NPP         [kg C/(m²·s)]  Net primary production = GPP - R_a_plant (≥ 0)
%                                 Clamped to ≥0 per Sitch Eq 24 (LPJ does not
%                                 allow negative NPP; deficit drawn from reserves).
%     R_a_plant   [kg C/(m²·s)]  Total autotrophic resp = R_m + R_g (= GPP - NPP)
%     diag        struct, fields R_m_leaf / R_m_root / R_m_stem / R_m_fruit / R_g
%                                 all in [kg C/(m²·s)] (diagnostic only).
%                   .NPP_net  [kg C/(m²·s)]  Unclamped net carbon = GPP - R_a_plant
%                                 (may be negative at night/winter). Diagnostic
%                                 only — use this (NOT the clamped NPP field) for
%                                 cumulative carbon-use efficiency CUE = ΣNPP_net/ΣGPP,
%                                 because clamping NPP to ≥0 discards winter/night
%                                 net-carbon losses and inflates the annual sum.
%
%   Inputs:
%     GPP_kg       [kg C/(m²·s)]  Canopy gross photosynthesis from ebal
%     Rd_leaf_kg   [kg C/(m²·s)]  Canopy leaf dark respiration from biochem
%                                  (Farquhar 1980, sun/shade aggregated in ebal)
%     N_root_fine  [g N/m²]       Fine-root N pool (prognostic state)
%     N_livestem   [g N/m²]       Live stem N pool [13c] (respires; herbaceous PFT = 0)
%     N_livecroot  [g N/m²]       Live coarse-root N pool [coarse-root] (respires, belowground;
%                                  uses soil-temperature g(T); non-woody PFT = 0)
%     N_fruit      [g N/m²]       Fruit N pool (respires for perennial fruit trees;
%                                  non-fruit PFT = 0, so R_m,fruit vanishes)
%     T_air        [°C]           Canopy air temperature
%     T_soil       [°C]           Root-zone mean soil temperature
%     pftParams    struct         Required fields (read from PFT 表):
%                                   .r_m10           [g C/(g N · d)] @ 10°C base
%                                                    Sitch 2003 Table 3 line
%                                                    "Plant Respiration": r
%                                                    (temperate woody/herb = 0.066,
%                                                     tropical TrBE/TrBR = 0.011)
%                                   .R_g_fraction    [-]   Growth resp fraction
%                                                    Sitch 2003 p.169 col 2 line 9:
%                                                    "25% of the remainder is
%                                                    taken as growth respiration"
%                                                    (Ryan 1991)
%
%   Method (Sitch et al. 2003 GCB 9:161-185, p.169):
%
%     1. Lloyd & Taylor (1994) Arrhenius temperature factor (Sitch Eq 23):
%
%            g(T) = exp{ 308.56 × [1/56.02 − 1/(T_C + 46.02)] }
%
%        with T_C in °C. Note: the +46.02 term acts as a low-T asymptote
%        (g → 0 as T_C → -46.02 °C), avoiding the unbounded Q10 blow-up at low
%        temperatures. T_C is clamped to [-40, 50] before evaluation.
%
%     2. Per-tissue maintenance respiration (Sitch Eq 22, math-equivalent form):
%
%            R_root    = r_m10 × N_root  × g(T_soil)   [g C/(m²·d)]
%            R_sapwood = r_m10 × N_stem  × g(T_air)    (no f_phen — wood does not
%                                                       deciduous-shed)
%            R_leaf    = Rd_leaf (from biochem, replaces Sitch Eq 22a so we
%                                  do not double-count canopy Rd)
%
%        Conversion to SI: g C/(m²·d) → kg C/(m²·s) factor = 1e-3 / 86400.
%
%     3. Growth respiration (Sitch p.169 col 2; Ryan 1991):
%
%            R_g = R_g_fraction × max(GPP - R_m, 0)        (R_g_fraction = 0.25)
%
%     4. NPP closure (Sitch Eq 24 form):
%
%            NPP       = max(GPP - R_m - R_g, 0)
%            R_a_plant = R_m + R_g
%
%   Reference:
%     Sitch S, Smith B, Prentice IC et al. (2003) Evaluation of ecosystem
%       dynamics, plant geography and terrestrial carbon cycling in the LPJ
%       dynamic global vegetation model. Global Change Biology 9:161-185.
%       Eq 22-24 + Table 3 + p.169 (PDF: data/references/Global_Change_Biology_-_2003_-_Sitch_-_*.pdf)
%     Lloyd J, Taylor JA (1994) On the temperature dependence of soil
%       respiration. Functional Ecology 8:315-323.
%       (cited in Sitch 2003 Eq 23 line "modified Arrhenius equation")
%     Ryan MG (1991) Effects of climate change on plant respiration.
%       Ecological Applications 1:157-167.
%       (cited in Sitch 2003 p.169 col 2 for both r-values and 25% R_g rule)
%     Farquhar GD, von Caemmerer S, Berry JA (1980) A biochemical model of
%       photosynthetic CO2 assimilation in leaves of C3 species. Planta
%       149:78-90. (basis of biochem Rd_leaf in ebal)
%
%   Audit history:
%     2026-05-19  Original M2-B.1b version (Q10=2 + r_m25 in s⁻¹) — values
%                 lacked published source. Switched to Sitch 2003 Eq 22-23 on
%                 2026-05-26 after literature audit (see SKILL.md §5.5).

% --- Temperature response g(T): Sitch 2003 Eq 23 / Lloyd & Taylor 1994 ---
% Clamp temperatures to physical range before evaluating Arrhenius.
% Lower bound -40°C: well above the -46.02°C singularity in Sitch Eq 23.
T_air_safe  = max(min(T_air,  50), -40);   % [°C]
T_soil_safe = max(min(T_soil, 50), -40);   % [°C]

g_T_air  = exp(308.56 * (1/56.02 - 1/(T_air_safe  + 46.02)));   % [-] Sitch Eq 23
g_T_soil = exp(308.56 * (1/56.02 - 1/(T_soil_safe + 46.02)));   % [-] Sitch Eq 23

use_clm_crop_resp = isfield(pftParams, 'clm_crop_respiration') ...
    && pftParams.clm_crop_respiration == 1;
if use_clm_crop_resp
    % CLM5 CNMRespMod.F90 / CoLM MOD_BGC_Veg_CNMResp.F90:
    % R_i = N_i * br_i * Q10^((T_i-20)/10), br_i [g C gN-1 s-1].
    q10_crop = pftParams.clm_resp_Q10;
    tc_air_crop = q10_crop^((T_air_safe - 20) / 10);
    tc_soil_crop = q10_crop^((T_soil_safe - 20) / 10);
end

% --- R_m,root (Sitch Eq 22c, leaf-out flag f absorbed in N pool dynamics) ---
% r_m10 × N has units gC/(gN·d) × gN/m² = gC/(m²·d) → ×1e-3/86400 to kg C/(m²·s)
if use_clm_crop_resp
    R_m_root_kg = N_root_fine * pftParams.clm_br_root_s * tc_soil_crop * 1e-3;
else
    R_m_root_kg = pftParams.r_m10 * N_root_fine * g_T_soil * (1e-3 / 86400);
end

% --- R_m,live stem [13c]: only live wood respires, dead wood (heartwood) does not ---
% CLM5 CNMRespMod.F90:265  livestem_mr = livestemn * br * tc  (deadstem not in eqn)
% Sitch 2003 GCB 9:161 state-var table p.163: maintenance respiration ~ "sapwood C"
% Dual-pool [13c]: use N_livestem directly (was single N_stem_wood x f_livewood proxy).
if use_clm_crop_resp
    R_m_stem_kg = N_livestem * pftParams.clm_br_s * tc_air_crop * 1e-3;
else
    R_m_stem_kg = pftParams.r_m10 * N_livestem * g_T_air * (1e-3 / 86400);
end

% --- R_m,live coarse root [coarse-root 2026-06-12]: only LIVE coarse root respires (dead does not) ---
% CLM5 CNMRespMod.F90:266  livecroot_mr = livecrootn * br_root * tc  (deadcroot NOT in eqn).
% Belowground organ -> use soil-temperature factor g_T_soil (mirrors R_m,root). N_livecroot = 0 for
% non-woody PFTs, so this term vanishes there.
if use_clm_crop_resp
    R_m_livecroot_kg = 0; % CLM crops have no live coarse-root maintenance term
else
    R_m_livecroot_kg = pftParams.r_m10 * N_livecroot * g_T_soil * (1e-3 / 86400);
end

% --- R_m,fruit [FruitTree]: perennial fruit pool respires (apple); non-fruit N_fruit=0 ---
% CLM5 FruitTree CNMRespMod.F90:267-268 (O.Dombrowski perennial addition):
%   if (perennial) grain_mr = grainn * br * tc   — same form as livestem_mr (line 265)
% Fruit is aboveground, so use canopy-air temperature factor g_T_air (mirrors R_m,stem).
% N_fruit = 0 for herbaceous / non-fruit PFTs, so this term vanishes there.
if use_clm_crop_resp
    R_m_fruit_kg = N_fruit * pftParams.clm_br_s * tc_air_crop * 1e-3;
else
    R_m_fruit_kg = pftParams.r_m10 * N_fruit * g_T_air * (1e-3 / 86400);
end

% --- R_m,leaf: from SCOPE biochem (Farquhar 1980, replaces Sitch Eq 22a) ---
% Avoids double-counting with canopy Rd already returned by ebal.
R_m_leaf_kg = Rd_leaf_kg;

% --- Total maintenance respiration ---
R_m_kg = R_m_leaf_kg + R_m_root_kg + R_m_stem_kg + R_m_livecroot_kg + R_m_fruit_kg;

% --- Growth respiration (Sitch 2003 p.169 col 2 line 9: "25% of the remainder") ---
R_g_kg = pftParams.R_g_fraction * max(GPP_kg - R_m_kg, 0);

% --- NPP closure (Sitch Eq 24 form) ---
NPP       = max(GPP_kg - R_m_kg - R_g_kg, 0);
R_a_plant = R_m_kg + R_g_kg;

% --- Unclamped net carbon for diagnostics (may be negative) ---
% NPP above is clamped ≥0 per Sitch Eq 24; cumulating that field over a year
% drops the night/winter net-carbon deficits (GPP < R_a) and overstates the
% annual total, so a CUE computed from ΣNPP/ΣGPP is biased high. NPP_net keeps
% the sign so ΣNPP_net/ΣGPP recovers the true annual CUE (= 1 − ΣR_a/ΣGPP).
NPP_net = GPP_kg - R_a_plant;   % [kg C/(m²·s)] diagnostic, unclamped

% --- Diagnostic struct (only-write, never read by physics) ---
diag.R_m_leaf  = R_m_leaf_kg;
diag.R_m_root  = R_m_root_kg;
diag.R_m_stem  = R_m_stem_kg;
diag.R_m_livecroot = R_m_livecroot_kg;   % [coarse-root 2026-06-12]
diag.R_m_fruit = R_m_fruit_kg;
diag.R_g       = R_g_kg;
diag.NPP_net  = NPP_net;
diag.g_T_air  = g_T_air;
diag.g_T_soil = g_T_soil;

end

% Editor: Jizixing Zhu; Date: 2026-08-24.
% Change: 2026-08-24 Jizixing Zhu -- add the fixed CLM5/CoLM crop
% maintenance-respiration branch without changing the Apple path.
