function PFTs = getPFTParameterSet()
%GETPFTPARAMETERSET  Plant Functional Type parameter table for COMPASS (25 PFTs).
% Editor: Wenjie Liu; Date: 2026-10-22.
% Change: 2026-10-22 Wenjie Liu -- set the Changwu winter-wheat grain
% target C:N to 20 g C/g N; other PFT parameters retain their prior values.
%
% Layout: one block per PFT (CLM5 pft-physiology style). To add a new PFT, copy
% any block and modify fields. To add a new field, add to ALL 25 PFT blocks below.
%
% Tunability legend (in front of each field):
%   [const]  CLM5 hardcoded constant (Fortran `parameter`), should not be changed
%   [tune]   Tunable per-PFT (may calibrate to site observations)
%   [site]   COMPASS site-calibrated value (Apple-orchard reverse-fit)
%   [struct] Structural / morphological, varies by species
%   [reserved] M3/future placeholder, currently unused
%
% =========================================================================
% Field groups (in order within each block):
% =========================================================================
%
% [A] Identification
%   PFT             [string]       Short tag (e.g. 'NET_temp', 'Apple') –used by getPFTValue
%   FullName        [string]       Human-readable name
%
% [B] PFT properties (M2-B.0)
%   photo_type      [string]       'C3' or 'C4' (currently all C3)
%   is_perennial    [-, {0,1}]     1=perennial (carry biomass year-to-year)
%   is_evergreen    [-, {0,1}]     1=evergreen (leaves not shed each year)
%   has_fruit       [-, {0,1}]     1=produces fruit pool (M2-B.5 harvest)
%   has_storage_pool [-, {0,1}]    1=has N_storage seasonal reserve (M2-B.4)
%
% [C] Stomatal conductance (Medlyn 2011 / Ball-Berry, CLM5 Table 8.1)
%   g1, g0 REMOVED 2026-06-02: Medlyn slope/intercept are read at runtime from
%     config_vegetation.txt (gs_g1 / gs_g0) via COMPASS.m get_veg_property, NOT
%     from this PFT table. The old per-PFT g1/g0 entries were dead values never
%     read by biochemical.m (it consumes leafbio.medlyn_g1/g0 sourced from cfg.Veg).
%   m               [-]            [tune] Ball-Berry slope (legacy)
%   alpha           [-]            [tune] Quantum yield (0.04-0.08)
%
% [D] Structure + scaling
%   Ztop            [m]            [struct] Canopy top height
%   Zbot            [m]            [struct] Canopy bottom height
%   SLA0            [m²/g C]       [struct/site] Top-canopy SLA (Apple 0.0125 = Wang 2009
%   SLAm            [m²/g C / LAI] [struct] dSLA/dLAI canopy gradient (CLM4.5 Tech Note Table 8.1 dsladlai);
%   CN_L            [g C/g N]      [tune] Leaf C:N ratio (target)
%   FLNR            [-]            [tune] Fraction of leaf N in Rubisco (CLM4 Table 8.1)
%   fN              [-]            [tune] Empirical N limitation factor (legacy, replaced by LUNA)
%
% [E] Hydraulic vulnerability (CLM5 Table 8.1)
%   psi_o           [mm]           [tune] Onset of stomatal closure
%   psi_c           [mm]           [tune] Complete stomatal closure
%
% [F] Photosynthetic capacity (legacy, replaced by LUNA when cfg.Veg.lLUNA=1)
%   Vcmax25         [µmol/m²/s]    [site] N-sufficient potential @ 25°C (Apple: 105 from Wang 2009)
%   Vcmax25_fN      [µmol/m²/s]    [tune] N-stressed potential
%   g_m             [mol/m²/s]     [tune] Mesophyll conductance at 25°C (gm25). A finite value
%                                  turns on the mesophyll CO2 drawdown (Cc = Ci − A/g_m < Ci) in
%                                  vegetation.biochemical.m; Inf disables it (standard Farquhar Cc=Ci).
%                                  Per-PFT values from Knauer et al. (2019) Glob Change Biol
%                                  25(5):1820-1838, doi:10.1111/gcb.14604, Table 3 — a JSBACH
%                                  land-surface-model literature synthesis of gm,max25 by PFT.
%                                  COMPASS PFTs mapped to Knauer functional types:
%                                    needleleaf-evergreen (NET) → ENF 0.078
%                                    needleleaf-deciduous (NDT) → DNF 0.057
%                                    broadleaf-evergreen tropical (BET-trop) → TRF 0.152
%                                    broadleaf-evergreen temperate + evergreen shrub (BET-temp,BES)
%                                                                    → EBF 0.101
%                                    broadleaf-deciduous tropical (BDT-trop) → TDF 0.058
%                                    broadleaf-deciduous temp/boreal + apple + robinia (BDT,PFT23,24)
%                                                                    → DBF 0.175
%                                    broadleaf-deciduous shrub (BDS) → DSH 0.098
%                                    C3 grass (C3G) → 0.197,  C4 grass (C4G) → 0.452
%                                    C3 crop (C3C) → 0.295,   C4 crop (C4C) → 0.739
%                                  Units [mol m⁻² s⁻¹] used directly as the per-bar form
%                                  [mol m⁻² s⁻¹ bar⁻¹] biochemical.m expects (≈ equal at ~1 bar
%                                  total pressure, mole-fraction↔partial-pressure).
%
% [G] LUNA feedback (Ali 2016 GMD 9:587; CLM5 LunaMod.F90)
%   Jmaxb1          [-]   [tune]   LUNA PAR-response (CLM5 namelist default 0.093563;
%                                  Apple [site] 0.1745 reverse-fit to Vcmax25=105 from Wang 2009)
%                                  →ONLY field actually read by physics: COMPASS.m:864 calls
%                                    vegetation.getPFTValue(pft,'Jmaxb1') to override lunaParams default.
%   % [REMOVED 2026-05-26] Rd_fraction / Jmaxb0 / N_struct_per_LMA / N_photo_ref /
%   %   tc_j0 / H_LUNA / NUE_Vcmax25 / NUE_Jmax25 were dead fields –never consumed
%   %   by physics code. LUNA uses hardcoded CLM5 constants in lunaParams.m:
%   %   Jmaxb0=0.0311, SNC=0.004, Wc2Wjb0=0.8054, relhExp=6.0999.
%   %   Rd_fraction (Rd25/Vcmax25 ratio) has no consumer in COMPASS –CLM5 uses it
%   %   in PhotosynthesisMod.F90 but COMPASS LUNA's RESP includes Rd implicitly.
%
% [H] Organ C:N + allocation (Yuan 2011, Mao 2018; Sitch 2003 LPJ Table 1)
%   CN_target_root  [g C/g N]      [tune] Fine-root C:N target
%   CN_target_stem  [g C/g N]      [tune] Stem/wood C:N target (annual PFT = 100 placeholder)
%   CN_target_fruit [g C/g N]      [tune] Fruit/seed C:N target
%   f_alloc_leaf    [-]            [tune] NPP fraction →leaf (
%   f_alloc_root    [-]            [tune] NPP fraction →fine-root
%   f_alloc_stem    [-]            [tune] NPP fraction →stem
%   f_alloc_fruit   [-]            [tune] NPP fraction →fruit (DOY 150-290 active for Apple)
%
% [I] Nacry root membrane (Nacry 2013 Plant Soil 370:1, Couturier 2007)
%   Vmax_NH4_umol_g_h [µmol N/(g-root·h)]  [tune] HATS NH4 max uptake
%   Vmax_NO3_umol_g_h [µmol N/(g-root·h)]  [tune] HATS NO3 max uptake
%   Km_NH4          [mg/L]         [tune] Half-saturation NH4 (≈0.3 for Apple)
%   Km_NO3          [mg/L]         [tune] Half-saturation NO3
%   Cmin_NH4        [mg/L]         [tune] Physiological lower bound NH4 uptake
%   Cmin_NO3        [mg/L]         [tune] Physiological lower bound NO3 uptake
%   SRL             [m-root/g-root] [struct] Specific root length (Wells & Eissenstat 2002, Apple≈0)
%   Vmax_NO3_LATS_umol_g_h [reserved] [reserved] M3 LATS transporter
%   Km_LATS_NO3     [reserved]     [reserved] M3 LATS half-saturation
%
% [J] Phenology
%   leaf_long_days  [days]         [struct] Leaf longevity (Niinemets 2010)
%   retransl_fraction [-]          [tune] Senescence N retransloc fraction (Aerts 1996: 0.5-0.7)
%
% [K] Plant maintenance + growth respiration (Sitch et al. 2003 GCB 9:161, p.169)
%   PDF: data/references/Global_Change_Biology_-_2003_-_Sitch_-_*.pdf
%
%   r_m10           [g C/(g N·d)]  [const] Maintenance respiration rate at 10°C base.
%                                  Sitch 2003 Table 3 line "Plant Respiration" –r:
%                                    Temperate / Boreal / temperate herbaceous: 0.066
%                                    Tropical (TrBE, TrBR, TrH = C4 herbaceous): 0.011
%                                  Cited basis in Sitch p.169 col 1: Fukai & Silsburg 1977;
%                                    Lechowicz 1980; Amthor 1989; Ryan 1991 (warmer-climate
%                                    plants have lower respiration at any given T).
%                                  Consumed by g(T)-Arrhenius (Eq 23) in computeNPP.m.
%
%   R_g_fraction    [-]            [const] Growth respiration as fraction of (GPP ∝R_m).
%                                  Sitch 2003 p.169 col 2 line 9: "25% of the remainder
%                                  is taken as growth respiration" (Ryan 1991). = 0.25
%                                  for ALL PFTs (Sitch makes no PFT differentiation).
%
%   [REMOVED 2026-05-26] r_m25_root / r_m25_stem / Q10_resp were COMPASS-specific
%     extensions with values not traceable to published sources (Sitch uses a single
%     r per PFT, no organ split; Sitch uses Arrhenius g(T), not Q10). Replaced by
%     Sitch 2003 Eq 22-23 strict implementation in src/+vegetation/+nitrogen/computeNPP.m.
%     See SKILL.md §5.5 (anti-hallucination rule) for the audit.
%
% [L] Pool initial values –DEPRECATED & REMOVED (2026-05-30, C-primary):
%   N_storage_init
%
% [M] Turnover + Harvest + Spring N_storage release (M2-B.4 2026-05-21; Phase D 2026-05-26)
%   harvest_fraction      [-]       [site] N_fruit × this = N
%   root_turnover_rate    [1/day]   [site]
%   spring_release_days   [days]    [site] N_storage
%
%   [REMOVED 2026-05-26 Phase D] harvest_threshold + spring_release_endphu
%     were frPHU-based site-specific hacks (
%     0.90/0.20 scaled down to 0.60/0.13). Replaced by CLM5+LPJ-GUESS phenology
%     dispatcher (src/+vegetation/+integrated/+phenology/) which drives turnover via
%     PhenologyState flags (offset_flag, aphen, budburst_doy) –site-independent,
%     climate-adaptive.
%
% Source: CLM4 Table 2.1+8.1 (Oleson 2013), CLM5 Tech Note (Lawrence 2019),
%         EPIC (Williams 1989, harvest_threshold), Aerts 1996 (retransl),
%         Phase C M2-B (Apple specific: Wang 2009, Wells & Eissenstat 2002, Mao 2018).
%
% [N] Phenology classification (2026-05-26, CLM5 architecture)
%   Source: CLM5 pftconMod.F90:200-202 (D:/OneDrive/CTSM-release-clm5.0.37/src/main/)
%           Each PFT gets 4 binary flags, EXACTLY ONE = 1, the others = 0.
%           Dispatcher (src/+vegetation/+integrated/+phenology/phenologyDispatch.m)
%           routes each PFT to the correct branch based on the active flag.
%
%   evergreen        [0/1]  Continuous leaf turnover, no budburst/senescence event
%                          (CLM5 CNEvergreenPhenology, CNPhenologyMod.F90:502-666)
%                          PFT 1,2,4,5,9 (NET temperate/boreal, BET tropical/temperate,
%                                          BES temperate)
%
%   season_decid     [0/1]  Cold-triggered budburst (Sykes 1996) + daylength senescence
%                          (CLM5 CNSeasonDecidPhenology, CNPhenologyMod.F90:669-971)
%                          (LPJ-GUESS growth.cpp:91-99 Sykes 1996 formula)
%                          PFT 3,7,8,10,11,23,24 (NDT/BDT boreal/temperate, BDS,
%                                                  Apple, Robinia)
%
%   stress_decid     [0/1]  Temperature/water-stress-triggered phenology
%                          (CLM5 CNStressDecidPhenology, CNPhenologyMod.F90:974-1417)
%                          PFT 6,12,13,14 (BDT tropical raingreen, C3/C4 grass)
%
%   crop_phenology   [0/1]  Annual sowing
%                          (CLM5 CropPhenology, CNPhenologyMod.F90:1420-2031)
%                          PFT 15-22 (corn, wheat, soybean, cotton, rice, sugarcane etc.)
%
%   [N.4] Onset/offset ramp duration (CLM5, all PFTs)
%     ndays_on   [days] Onset ramp (CLM5 CNPhenologyMod.F90:51; clm5_params.nc scalar = 30)
%     ndays_off  [days] Offset ramp (CLM5 CNPhenologyMod.F90:52; clm5_params.nc scalar = 15)
%
%   [N.5] Fruit phenology sub-state machine (Apple PFT 23 only –Wang 2009 LMWE
%         observations; other has_fruit PFTs 15-22: NaN placeholders pending
%         future calibration with their site-specific bloom/harvest data)
%     GDD_bloom_after_budburst [°C·d] Bud break →bloom thermal requirement
%                              Apple Wang 2009: DOY 80→30 = 50 days, T_base=5.59
%                              →CHU integral ≈250 °C·d
%     TSum_fruit               [°C·d] Bloom →harvest thermal requirement
%                              Apple Wang 2009: DOY 130→85 = 155 days
%                              →CHU integral ≈2132 °C·d
%     harvest_DS               [-]    Maturity threshold for harvest trigger
%                              DAISY phenology_TSum.C uses DS=2.0 as maturity;
%                              EPIC harvest index (Williams 1989) uses
%                              fphu=0.9-0.95 standard. We pick 0.95 (mid-range).
%
% Editor: Jizixing Zhu; Date: 2026-07-03.
% Change: 2026-07-03 Jizixing Zhu — add per-PFT g_m (mesophyll conductance gm25, group [F]) to all 24 blocks; values from Knauer 2019 GCB Table 3.
% Change: 2026-07-16 Jizixing Zhu — add foliage clumping index to all 24 PFTs; natural-vegetation values follow ClimaLand pfts.jl:95-463, while crop and Robinia entries are documented functional-type proxies.
% Change: 2026-07-16 Jizixing Zhu — replace the Apple broadleaf-tree proxy with the orchard measurement reported by Liu et al. (2013), Sect. 3.4.4 and Fig. 10.
% Change: 2026-08-09 Jizixing Zhu -- add PFT 25 Winter Wheat from the
%   user-supplied parameter block. Active COMPASS fields are explicit below;
%   candidate CLM5 crop-calendar fields are retained as unconsumed reserved
%   metadata and are not claimed as active phenology. [SOURCE UNVERIFIED 2026-08-09]
% Source for PFT Omega values: Chen et al. (2012), Global Biogeochemical
% Cycles 26, GB1019, Table 1, doi:10.1029/2010GB003996, as tabulated by
% ClimaLand. Crop and Robinia entries below are explicit proxies; the Apple
% entry uses a species-specific orchard measurement cited on that field.

% [P1.5-b 2026-06-16] Memoise the constant 25-PFT table — it is rebuilt from
% literals on every call (getPFTValue calls this per lookup, ~18k/run). The table
% is constant EXCEPT PFTs(23).labile_resp_days, which a parameter sweep can
% override via COMPASS_LABILE_RESP_DAYS; that one field is re-applied on every
% cache hit so the override stays live. (Constant-table cache — §3.3 module-
% constant exception.)
persistent CACHED_PFTS
if ~isempty(CACHED_PFTS)
    PFTs = CACHED_PFTS;
    PFTs(23).labile_resp_days = localLabileRespDays();
    return;
end

PFTs = repmat(struct(), 1, 25);

% ==============================================================================
% PFT 1: NET Temperate –Needleleaf evergreen tree –temperate
% ==============================================================================
PFTs(1).PFT                    = 'NET Temperate';
PFTs(1).FullName               = 'Needleleaf evergreen tree –temperate';
PFTs(1).clumping_index         = 0.74;  % [-] ClimaLand pfts.jl:95-102, temperate needleleaf evergreen
%   [HABIT] (PFT lifeform flags)
PFTs(1).photo_type             = 'C3';
PFTs(1).is_perennial           = 1;
PFTs(1).has_fruit              = 0;
PFTs(1).has_storage_pool       = 0;
%   [PHENOLOGY] (see [N.0]-[N.5] sub-stages below)
PFTs(1).evergreen            = 1;
PFTs(1).season_decid         = 0;
PFTs(1).stress_decid         = 0;
PFTs(1).crop_phenology       = 0;   % NET Temperate
%   [N.1] Sykes 1996 chilling-GDD (NaN: not season_decid, physics skips)
PFTs(1).k_chilla            = NaN;
PFTs(1).k_chillb            = NaN;
PFTs(1).k_chillk            = NaN;
PFTs(1).phengdd5ramp        = NaN;
%   [N.2] Senescence (CLM5 daylength trigger, CNPhenologyMod.F90:957)
PFTs(1).crit_dayl_h         = 10.9;    % [h] dayl threshold for offset trigger (CLM5 global default)
PFTs(1).ndays_on            = 30;      % [days] CLM5 onset ramp duration (CNPhenologyMod.F90:51, clm5_params.c171117.nc scalar)
PFTs(1).ndays_off           = 15;      % [days] CLM5 offset ramp duration (CNPhenologyMod.F90:52, clm5_params.c171117.nc scalar)
PFTs(1).fstor2tran          = 0.5;     % [-] CLM5 budburst fraction (CNPhenologyMod.F90:53, clm5_params.c171117.nc scalar = 0.5)
PFTs(1).GDD_bloom_after_budburst = NaN;     % [°C·d] only Apple PFT 23 has fruit phenology calibrated for Changwu
PFTs(1).TSum_fruit              = NaN;
PFTs(1).harvest_DS              = NaN;
%   [N.3] N retranslocation on tissue turnover (LPJ-GUESS growth.cpp:256-285)
PFTs(1).nrelocfrac          = 0.0;     % [-] N retransloc fraction to N_storage at senescence
%   [STOMATA]
PFTs(1).m                      = 4.4713;   % [CLM5 NC audit 2026-05-26] mbbopt, was 9
PFTs(1).alpha                  = 0.06;
%   [LEAF]
PFTs(1).Ztop                   = 0;
PFTs(1).Zbot                   = 0;
PFTs(1).SLA0                   = 0.01;
PFTs(1).SLAm                   = 0.00125;
PFTs(1).CN_L                   = 58.00;    % [CLM5 audit 2026-05-26] was 35
%   [PHOTOSYNTH]
PFTs(1).Vcmax25                = 62.5;
PFTs(1).g_m                    = 0.078;  % [mol m-2 s-1] gm25 (Knauer 2019 GCB25 Table 3 ENF, needleleaf evergreen)
%   [LUNA]
PFTs(1).Jmaxb1                 = 0.0936;     % [-] LUNA PAR-response, CLM5 global calibrated (namelist_defaults_clm4_5.xml:430); Ali 2016 GMD 9:587 sensitivity range 0.05-0.30
%   [ALLOC_CN]
PFTs(1).CN_target_root         = 42.00;    % [CLM5 NC audit 2026-05-26] was 60
PFTs(1).CN_target_stem         = 50.00;    % [CLM5 NC audit 2026-05-26] was 200
% [13c dual stem pool] CLM5 CN_Allocation Tbl col7-8 (CN_lw/CN_dw); clm5_params (flivewd/lwtop_ann)
PFTs(1).CN_lw                  = 50;    % [gC/gN] live stem (respires)
PFTs(1).CN_dw                  = 500;   % [gC/gN] dead stem (heartwood, ~N-free)
PFTs(1).flivewd                = 0.1;  % [-] new-stem-C frac to live wood (tree/shrub .1/grass 0/crop 1)
% croot_stem: coarse-root:stem allocation ratio -> routes coarse-root share of stem
%   allocation into belowground (root) reserve at storage (CLM5
%   NutrientCompetitionFlexibleCNMod.F90:529, natural tree/shrub value 0.30; clm5_params
%   .c171117__FruitTree.nc verified 2026-06-07). Consumed in computeCarbonAllocation.
PFTs(1).croot_stem             = 0.30;  % [-] coarse-root:stem ratio (CLM5 natural tree/shrub)
PFTs(1).froot_leaf             = 1.5;   % [-] f1 fine-root:leaf  (c_allometry alloc; CLM5 nc idx 1)
PFTs(1).croot_stem_alloc       = 0.30;  % [-] f2 coarse:stem ALLOC flux (c_allometry; vs reserve croot_stem)
PFTs(1).stem_leaf              = 2.3;   % [-] f3 stem:leaf  (c_allometry alloc; CLM5 nc idx 1)
PFTs(1).lwtop                  = 0.7;   % [1/yr] live->dead wood turnover (CLM5 lwtop_ann)
PFTs(1).CN_target_fruit        = 0;
% [dead-value guard 2026-06-28 Jizixing Zhu] evergreen woody (allometry + has_storage_pool=0):
%   f_alloc_* never read (allometry split + no spring onset-xfer); NaN so any accidental read surfaces.
PFTs(1).f_alloc_leaf           = NaN;   % [-] DEAD: allometry, not f_alloc
PFTs(1).f_alloc_root           = NaN;   % [-] DEAD
PFTs(1).f_alloc_stem           = NaN;   % [-] DEAD
PFTs(1).f_alloc_fruit          = 0;     % [-] no fruit (has_fruit=0)
%   [N_UPTAKE] (Nacry kinetics)
PFTs(1).Vmax_NH4_umol_g_h      = 8;
PFTs(1).Vmax_NO3_umol_g_h      = 4;
PFTs(1).Km_NH4                 = 0.3;
PFTs(1).Km_NO3                 = 0.3;
PFTs(1).Cmin_NH4               = 0.5;
PFTs(1).Cmin_NO3               = 0.25;
PFTs(1).SRL                    = 20;
PFTs(1).Vmax_NO3_LATS_umol_g_h = 0;
PFTs(1).Km_LATS_NO3            = 70;
%   [LEAF_LIFESPAN]
PFTs(1).leaf_long_days         = 1207.85;    % [CLM5 NC audit 2026-05-26] was 600
PFTs(1).retransl_fraction      = 0.5;
%   [RESPIRATION] Plant respiration (Sitch 2003 GCB 9:161 Eq 22-23, p.169)
PFTs(1).r_m10                = 0.066;     % [g C/(g N·d)] @ 10°C base, Sitch 2003 Table 3 "Plant Respiration" –TeNE –Sitch Table 3 default (temperate)
PFTs(1).R_g_fraction         = 0.11;      % [-] growth resp factor (CLM5 grperc; was Sitch 0.25).
PFTs(1).harvest_fraction         = 0.90;
PFTs(1).root_turnover_rate       = 0.00274;
PFTs(1).spring_release_days      = 60;

% ==============================================================================
% PFT 2: NET Boreal –Needleleaf evergreen tree –boreal
% ==============================================================================
PFTs(2).PFT                    = 'NET Boreal';
PFTs(2).FullName               = 'Needleleaf evergreen tree –boreal';
PFTs(2).clumping_index         = 0.74;  % [-] ClimaLand pfts.jl:123-130, boreal needleleaf evergreen
%   [HABIT] (PFT lifeform flags)
PFTs(2).photo_type             = 'C3';
PFTs(2).is_perennial           = 1;
PFTs(2).has_fruit              = 0;
PFTs(2).has_storage_pool       = 0;
%   [PHENOLOGY] (see [N.0]-[N.5] sub-stages below)
PFTs(2).evergreen            = 1;
PFTs(2).season_decid         = 0;
PFTs(2).stress_decid         = 0;
PFTs(2).crop_phenology       = 0;   % NET Boreal
%   [N.1] Sykes 1996 chilling-GDD (NaN: not season_decid, physics skips)
PFTs(2).k_chilla            = NaN;
PFTs(2).k_chillb            = NaN;
PFTs(2).k_chillk            = NaN;
PFTs(2).phengdd5ramp        = NaN;
%   [N.2] Senescence (CLM5 daylength trigger, CNPhenologyMod.F90:957)
PFTs(2).crit_dayl_h         = 10.9;    % [h] dayl threshold for offset trigger (CLM5 global default)
PFTs(2).ndays_on            = 30;      % [days] CLM5 onset ramp duration (CNPhenologyMod.F90:51, clm5_params.c171117.nc scalar)
PFTs(2).ndays_off           = 15;      % [days] CLM5 offset ramp duration (CNPhenologyMod.F90:52, clm5_params.c171117.nc scalar)
PFTs(2).fstor2tran          = 0.5;     % [-] CLM5 budburst fraction (CNPhenologyMod.F90:53, clm5_params.c171117.nc scalar = 0.5)
PFTs(2).GDD_bloom_after_budburst = NaN;     % [°C·d] only Apple PFT 23 has fruit phenology calibrated for Changwu
PFTs(2).TSum_fruit              = NaN;
PFTs(2).harvest_DS              = NaN;
%   [N.3] N retranslocation on tissue turnover (LPJ-GUESS growth.cpp:256-285)
PFTs(2).nrelocfrac          = 0.0;     % [-] N retransloc fraction to N_storage at senescence
%   [STOMATA]
PFTs(2).m                      = 7.1017;   % [CLM5 NC audit 2026-05-26] mbbopt, was 9
PFTs(2).alpha                  = 0.06;
%   [LEAF]
PFTs(2).Ztop                   = 0;
PFTs(2).Zbot                   = 0;
PFTs(2).SLA0                   = 0.0100;    % [CLM5 audit 2026-05-26] was 0.008
PFTs(2).SLAm                   = 0.001;
PFTs(2).CN_L                   = 58.00;    % [CLM5 audit 2026-05-26] was 40
%   [PHOTOSYNTH]
PFTs(2).Vcmax25                = 62.6;
PFTs(2).g_m                    = 0.078;  % [mol m-2 s-1] gm25 (Knauer 2019 GCB25 Table 3 ENF, needleleaf evergreen)
%   [LUNA]
PFTs(2).Jmaxb1                 = 0.0936;
%   [ALLOC_CN]
PFTs(2).CN_target_root         = 42.00;    % [CLM5 NC audit 2026-05-26] was 70
PFTs(2).CN_target_stem         = 50.00;    % [CLM5 NC audit 2026-05-26] was 200
% [13c dual stem pool] CLM5 CN_Allocation Tbl col7-8 (CN_lw/CN_dw); clm5_params (flivewd/lwtop_ann)
PFTs(2).CN_lw                  = 50;    % [gC/gN] live stem (respires)
PFTs(2).CN_dw                  = 500;   % [gC/gN] dead stem (heartwood, ~N-free)
PFTs(2).flivewd                = 0.1;  % [-] new-stem-C frac to live wood (tree/shrub .1/grass 0/crop 1)
% croot_stem: coarse-root:stem allocation ratio -> routes coarse-root share of stem
%   allocation into belowground (root) reserve at storage (CLM5
%   NutrientCompetitionFlexibleCNMod.F90:529, natural tree/shrub value 0.30; clm5_params
%   .c171117__FruitTree.nc verified 2026-06-07). Consumed in computeCarbonAllocation.
PFTs(2).croot_stem             = 0.30;  % [-] coarse-root:stem ratio (CLM5 natural tree/shrub)
PFTs(2).froot_leaf             = 1.5;   % [-] f1 fine-root:leaf  (c_allometry alloc; CLM5 nc idx 2)
PFTs(2).croot_stem_alloc       = 0.30;  % [-] f2 coarse:stem ALLOC flux (c_allometry; vs reserve croot_stem)
PFTs(2).stem_leaf              = 2.3;   % [-] f3 stem:leaf  (c_allometry alloc; CLM5 nc idx 2)
PFTs(2).lwtop                  = 0.7;   % [1/yr] live->dead wood turnover (CLM5 lwtop_ann)
PFTs(2).CN_target_fruit        = 0;
% [dead-value guard 2026-06-28 Jizixing Zhu] evergreen woody (allometry + has_storage_pool=0):
%   f_alloc_* never read (allometry split + no spring onset-xfer); NaN so any accidental read surfaces.
PFTs(2).f_alloc_leaf           = NaN;   % [-] DEAD: allometry, not f_alloc
PFTs(2).f_alloc_root           = NaN;   % [-] DEAD
PFTs(2).f_alloc_stem           = NaN;   % [-] DEAD
PFTs(2).f_alloc_fruit          = 0;     % [-] no fruit (has_fruit=0)
%   [N_UPTAKE] (Nacry kinetics)
PFTs(2).Vmax_NH4_umol_g_h      = 8;
PFTs(2).Vmax_NO3_umol_g_h      = 4;
PFTs(2).Km_NH4                 = 0.3;
PFTs(2).Km_NO3                 = 0.3;
PFTs(2).Cmin_NH4               = 0.5;
PFTs(2).Cmin_NO3               = 0.25;
PFTs(2).SRL                    = 20;
PFTs(2).Vmax_NO3_LATS_umol_g_h = 0;
PFTs(2).Km_LATS_NO3            = 70;
%   [LEAF_LIFESPAN]
PFTs(2).leaf_long_days         = 1207.85;    % [CLM5 NC audit 2026-05-26] was 800
PFTs(2).retransl_fraction      = 0.5;
%   [RESPIRATION] Plant respiration (Sitch 2003 GCB 9:161 Eq 22-23, p.169)
PFTs(2).r_m10                = 0.066;     % [g C/(g N·d)] @ 10°C base, Sitch 2003 Table 3 "Plant Respiration" –BoNE –Sitch Table 3 default (boreal)
PFTs(2).R_g_fraction         = 0.11;      % [-] growth resp factor (CLM5 grperc; was Sitch 0.25).
PFTs(2).harvest_fraction         = 0.90;
PFTs(2).root_turnover_rate       = 0.00274;
PFTs(2).spring_release_days      = 60;

% ==============================================================================
% PFT 3: NDT Boreal –Needleleaf deciduous tree –boreal
% ==============================================================================
PFTs(3).PFT                    = 'NDT Boreal';
PFTs(3).FullName               = 'Needleleaf deciduous tree –boreal';
PFTs(3).clumping_index         = 0.78;  % [-] ClimaLand pfts.jl:151-158, boreal needleleaf deciduous
%   [HABIT] (PFT lifeform flags)
PFTs(3).photo_type             = 'C3';
PFTs(3).is_perennial           = 1;
PFTs(3).has_fruit              = 0;
PFTs(3).has_storage_pool       = 1;
%   [PHENOLOGY] (see [N.0]-[N.5] sub-stages below)
PFTs(3).evergreen            = 0;
PFTs(3).season_decid         = 1;
PFTs(3).stress_decid         = 0;
PFTs(3).crop_phenology       = 0;   % NDT Boreal
%   [N.1] Sykes 1996 chilling-GDD for season_decid (LPJ-GUESS growth.cpp:91-97)
PFTs(3).k_chilla            = 0;       % [°C·d] Sykes 1996 a, LPJ-GUESS default
PFTs(3).k_chillb            = 350;     % [°C·d] Sykes 1996 b, species-specific (LPJ-GUESS europe.ins)
PFTs(3).k_chillk            = 0.05;    % [d^-1] Sykes 1996 k, LPJ-GUESS default
PFTs(3).phengdd5ramp        = 200;     % [°C·d] linear ramp to full leaf-on (LPJ-GUESS)
%   [N.2] Senescence (CLM5 daylength trigger, CNPhenologyMod.F90:957)
PFTs(3).crit_dayl_h         = 10.9;    % [h] dayl threshold for offset trigger (CLM5 global default)
PFTs(3).ndays_on            = 30;      % [days] CLM5 onset ramp duration (CNPhenologyMod.F90:51, clm5_params.c171117.nc scalar)
PFTs(3).ndays_off           = 15;      % [days] CLM5 offset ramp duration (CNPhenologyMod.F90:52, clm5_params.c171117.nc scalar)
PFTs(3).fstor2tran          = 0.5;     % [-] CLM5 budburst fraction (CNPhenologyMod.F90:53, clm5_params.c171117.nc scalar = 0.5)
PFTs(3).GDD_bloom_after_budburst = NaN;     % [°C·d] only Apple PFT 23 has fruit phenology calibrated for Changwu
PFTs(3).TSum_fruit              = NaN;
PFTs(3).harvest_DS              = NaN;
%   [N.3] N retranslocation on tissue turnover (LPJ-GUESS growth.cpp:256-285)
PFTs(3).nrelocfrac          = 0.5;     % [-] N retransloc fraction to N_storage at senescence
%   [STOMATA]
PFTs(3).m                      = 7.1017;   % [CLM5 NC audit 2026-05-26] mbbopt, was 9
PFTs(3).alpha                  = 0.06;
%   [LEAF]
PFTs(3).Ztop                   = 0;
PFTs(3).Zbot                   = 0;
PFTs(3).SLA0                   = 0.0202;    % [CLM5 audit 2026-05-26] was 0.024
PFTs(3).SLAm                   = 0.003;
PFTs(3).CN_L                   = 25.8100;    % [CLM5 audit 2026-05-26] was 25
%   [PHOTOSYNTH]
PFTs(3).Vcmax25                = 39.1;
PFTs(3).g_m                    = 0.057;  % [mol m-2 s-1] gm25 (Knauer 2019 GCB25 Table 3 DNF, needleleaf deciduous)
%   [LUNA]
PFTs(3).Jmaxb1                 = 0.0936;
%   [ALLOC_CN]
PFTs(3).CN_target_root         = 42.00;    % [CLM5 NC audit 2026-05-26] was 50
PFTs(3).CN_target_stem         = 50.00;    % [CLM5 NC audit 2026-05-26] was 250
% [13c dual stem pool] CLM5 CN_Allocation Tbl col7-8 (CN_lw/CN_dw); clm5_params (flivewd/lwtop_ann)
PFTs(3).CN_lw                  = 50;    % [gC/gN] live stem (respires)
PFTs(3).CN_dw                  = 500;   % [gC/gN] dead stem (heartwood, ~N-free)
PFTs(3).flivewd                = 0.1;  % [-] new-stem-C frac to live wood (tree/shrub .1/grass 0/crop 1)
% croot_stem: coarse-root:stem allocation ratio -> routes coarse-root share of stem
%   allocation into belowground (root) reserve at storage (CLM5
%   NutrientCompetitionFlexibleCNMod.F90:529, natural tree/shrub value 0.30; clm5_params
%   .c171117__FruitTree.nc verified 2026-06-07). Consumed in computeCarbonAllocation.
PFTs(3).croot_stem             = 0.30;  % [-] coarse-root:stem ratio (CLM5 natural tree/shrub)
PFTs(3).froot_leaf             = 1.5;   % [-] f1 fine-root:leaf  (c_allometry alloc; CLM5 nc idx 3)
PFTs(3).croot_stem_alloc       = 0.30;  % [-] f2 coarse:stem ALLOC flux (c_allometry; vs reserve croot_stem)
PFTs(3).stem_leaf              = 1.0;   % [-] f3 stem:leaf  (c_allometry alloc; CLM5 nc idx 3)
PFTs(3).lwtop                  = 0.7;   % [1/yr] live->dead wood turnover (CLM5 lwtop_ann)
PFTs(3).CN_target_fruit        = 0;
PFTs(3).f_alloc_leaf           = 0.35;
PFTs(3).f_alloc_root           = 0.25;
PFTs(3).f_alloc_stem           = 0.4;
PFTs(3).f_alloc_fruit          = 0;
%   [N_UPTAKE] (Nacry kinetics)
PFTs(3).Vmax_NH4_umol_g_h      = 8;
PFTs(3).Vmax_NO3_umol_g_h      = 4;
PFTs(3).Km_NH4                 = 0.3;
PFTs(3).Km_NO3                 = 0.3;
PFTs(3).Cmin_NH4               = 0.5;
PFTs(3).Cmin_NO3               = 0.25;
PFTs(3).SRL                    = 20;
PFTs(3).Vmax_NO3_LATS_umol_g_h = 0;
PFTs(3).Km_LATS_NO3            = 70;
%   [LEAF_LIFESPAN]
PFTs(3).leaf_long_days         = 184.93;    % [CLM5 NC audit 2026-05-26] was 200
PFTs(3).retransl_fraction      = 0.6;
%   [RESPIRATION] Plant respiration (Sitch 2003 GCB 9:161 Eq 22-23, p.169)
PFTs(3).r_m10                = 0.066;     % [g C/(g N·d)] @ 10°C base, Sitch 2003 Table 3 "Plant Respiration" –BoNS –Sitch Table 3 default (boreal needleleaf deciduous)
PFTs(3).R_g_fraction         = 0.11;      % [-] growth resp factor (CLM5 grperc; was Sitch 0.25).
PFTs(3).harvest_fraction         = 0.90;
PFTs(3).root_turnover_rate       = 0.00274;
PFTs(3).spring_release_days      = 60;

% ==============================================================================
% PFT 4: BET Tropical –Broadleaf evergreen tree –tropical
% ==============================================================================
PFTs(4).PFT                    = 'BET Tropical';
PFTs(4).FullName               = 'Broadleaf evergreen tree –tropical';
PFTs(4).clumping_index         = 0.66;  % [-] ClimaLand pfts.jl:179-186, tropical broadleaf evergreen
%   [HABIT] (PFT lifeform flags)
PFTs(4).photo_type             = 'C3';
PFTs(4).is_perennial           = 1;
PFTs(4).has_fruit              = 0;
PFTs(4).has_storage_pool       = 0;
%   [PHENOLOGY] (see [N.0]-[N.5] sub-stages below)
PFTs(4).evergreen            = 1;
PFTs(4).season_decid         = 0;
PFTs(4).stress_decid         = 0;
PFTs(4).crop_phenology       = 0;   % BET Tropical
%   [N.1] Sykes 1996 chilling-GDD (NaN: not season_decid, physics skips)
PFTs(4).k_chilla            = NaN;
PFTs(4).k_chillb            = NaN;
PFTs(4).k_chillk            = NaN;
PFTs(4).phengdd5ramp        = NaN;
%   [N.2] Senescence (CLM5 daylength trigger, CNPhenologyMod.F90:957)
PFTs(4).crit_dayl_h         = 10.9;    % [h] dayl threshold for offset trigger (CLM5 global default)
PFTs(4).ndays_on            = 30;      % [days] CLM5 onset ramp duration (CNPhenologyMod.F90:51, clm5_params.c171117.nc scalar)
PFTs(4).ndays_off           = 15;      % [days] CLM5 offset ramp duration (CNPhenologyMod.F90:52, clm5_params.c171117.nc scalar)
PFTs(4).fstor2tran          = 0.5;     % [-] CLM5 budburst fraction (CNPhenologyMod.F90:53, clm5_params.c171117.nc scalar = 0.5)
PFTs(4).GDD_bloom_after_budburst = NaN;     % [°C·d] only Apple PFT 23 has fruit phenology calibrated for Changwu
PFTs(4).TSum_fruit              = NaN;
PFTs(4).harvest_DS              = NaN;
%   [N.3] N retranslocation on tissue turnover (LPJ-GUESS growth.cpp:256-285)
PFTs(4).nrelocfrac          = 0.0;     % [-] N retransloc fraction to N_storage at senescence
%   [STOMATA]
PFTs(4).m                      = 7.2915;   % [CLM5 NC audit 2026-05-26] mbbopt, was 9
PFTs(4).alpha                  = 0.06;
%   [LEAF]
PFTs(4).Ztop                   = 0;
PFTs(4).Zbot                   = 0;
PFTs(4).SLA0                   = 0.0190;    % [CLM5 audit 2026-05-26] was 0.012
PFTs(4).SLAm                   = 0.0012;    % [CLM5 NC audit 2026-05-26] was 0.0015
PFTs(4).CN_L                   = 29.6000;    % [CLM5 audit 2026-05-26] was 30
%   [PHOTOSYNTH]
PFTs(4).Vcmax25                = 55;
PFTs(4).g_m                    = 0.152;  % [mol m-2 s-1] gm25 (Knauer 2019 GCB25 Table 3 TRF, tropical evergreen broadleaf)
%   [LUNA]
PFTs(4).Jmaxb1                 = 0.0936;
%   [ALLOC_CN]
PFTs(4).CN_target_root         = 42.00;    % [CLM5 NC audit 2026-05-26] was 40
PFTs(4).CN_target_stem         = 50.00;    % [CLM5 NC audit 2026-05-26] was 150
% [13c dual stem pool] CLM5 CN_Allocation Tbl col7-8 (CN_lw/CN_dw); clm5_params (flivewd/lwtop_ann)
PFTs(4).CN_lw                  = 50;    % [gC/gN] live stem (respires)
PFTs(4).CN_dw                  = 500;   % [gC/gN] dead stem (heartwood, ~N-free)
PFTs(4).flivewd                = 0.1;  % [-] new-stem-C frac to live wood (tree/shrub .1/grass 0/crop 1)
% croot_stem: coarse-root:stem allocation ratio -> routes coarse-root share of stem
%   allocation into belowground (root) reserve at storage (CLM5
%   NutrientCompetitionFlexibleCNMod.F90:529, natural tree/shrub value 0.30; clm5_params
%   .c171117__FruitTree.nc verified 2026-06-07). Consumed in computeCarbonAllocation.
PFTs(4).croot_stem             = 0.30;  % [-] coarse-root:stem ratio (CLM5 natural tree/shrub)
PFTs(4).froot_leaf             = 1.5;   % [-] f1 fine-root:leaf  (c_allometry alloc; CLM5 nc idx 4)
PFTs(4).croot_stem_alloc       = 0.30;  % [-] f2 coarse:stem ALLOC flux (c_allometry; vs reserve croot_stem)
PFTs(4).stem_leaf              = 2.3;   % [-] f3 stem:leaf  (c_allometry alloc; CLM5 nc idx 4)
PFTs(4).lwtop                  = 0.7;   % [1/yr] live->dead wood turnover (CLM5 lwtop_ann)
PFTs(4).CN_target_fruit        = 0;
% [dead-value guard 2026-06-28 Jizixing Zhu] evergreen woody (allometry + has_storage_pool=0):
%   f_alloc_* never read (allometry split + no spring onset-xfer); NaN so any accidental read surfaces.
PFTs(4).f_alloc_leaf           = NaN;   % [-] DEAD: allometry, not f_alloc
PFTs(4).f_alloc_root           = NaN;   % [-] DEAD
PFTs(4).f_alloc_stem           = NaN;   % [-] DEAD
PFTs(4).f_alloc_fruit          = 0;     % [-] no fruit (has_fruit=0)
%   [N_UPTAKE] (Nacry kinetics)
PFTs(4).Vmax_NH4_umol_g_h      = 8;
PFTs(4).Vmax_NO3_umol_g_h      = 4;
PFTs(4).Km_NH4                 = 0.3;
PFTs(4).Km_NO3                 = 0.3;
PFTs(4).Cmin_NH4               = 0.5;
PFTs(4).Cmin_NO3               = 0.25;
PFTs(4).SRL                    = 20;
PFTs(4).Vmax_NO3_LATS_umol_g_h = 0;
PFTs(4).Km_LATS_NO3            = 70;
%   [LEAF_LIFESPAN]
PFTs(4).leaf_long_days         = 511.91;    % [CLM5 NC audit 2026-05-26] was 500
PFTs(4).retransl_fraction      = 0.4;
%   [RESPIRATION] Plant respiration (Sitch 2003 GCB 9:161 Eq 22-23, p.169)
PFTs(4).r_m10                = 0.011;     % [g C/(g N·d)] @ 10°C base, Sitch 2003 Table 3 "Plant Respiration" –TrBE –Sitch Table 3 tropical value (broadleaf evergreen)
PFTs(4).R_g_fraction         = 0.11;      % [-] growth resp factor (CLM5 grperc; was Sitch 0.25).
PFTs(4).harvest_fraction         = 0.90;
PFTs(4).root_turnover_rate       = 0.00274;
PFTs(4).spring_release_days      = 60;

% ==============================================================================
% PFT 5: BET temperate –Broadleaf evergreen tree –temperate
% ==============================================================================
PFTs(5).PFT                    = 'BET temperate';
PFTs(5).FullName               = 'Broadleaf evergreen tree –temperate';
PFTs(5).clumping_index         = 0.66;  % [-] ClimaLand pfts.jl:207-214, temperate broadleaf evergreen
%   [HABIT] (PFT lifeform flags)
PFTs(5).photo_type             = 'C3';
PFTs(5).is_perennial           = 1;
PFTs(5).has_fruit              = 0;
PFTs(5).has_storage_pool       = 0;
%   [PHENOLOGY] (see [N.0]-[N.5] sub-stages below)
PFTs(5).evergreen            = 1;
PFTs(5).season_decid         = 0;
PFTs(5).stress_decid         = 0;
PFTs(5).crop_phenology       = 0;   % BET temperate
%   [N.1] Sykes 1996 chilling-GDD (NaN: not season_decid, physics skips)
PFTs(5).k_chilla            = NaN;
PFTs(5).k_chillb            = NaN;
PFTs(5).k_chillk            = NaN;
PFTs(5).phengdd5ramp        = NaN;
%   [N.2] Senescence (CLM5 daylength trigger, CNPhenologyMod.F90:957)
PFTs(5).crit_dayl_h         = 10.9;    % [h] dayl threshold for offset trigger (CLM5 global default)
PFTs(5).ndays_on            = 30;      % [days] CLM5 onset ramp duration (CNPhenologyMod.F90:51, clm5_params.c171117.nc scalar)
PFTs(5).ndays_off           = 15;      % [days] CLM5 offset ramp duration (CNPhenologyMod.F90:52, clm5_params.c171117.nc scalar)
PFTs(5).fstor2tran          = 0.5;     % [-] CLM5 budburst fraction (CNPhenologyMod.F90:53, clm5_params.c171117.nc scalar = 0.5)
PFTs(5).GDD_bloom_after_budburst = NaN;     % [°C·d] only Apple PFT 23 has fruit phenology calibrated for Changwu
PFTs(5).TSum_fruit              = NaN;
PFTs(5).harvest_DS              = NaN;
%   [N.3] N retranslocation on tissue turnover (LPJ-GUESS growth.cpp:256-285)
PFTs(5).nrelocfrac          = 0.0;     % [-] N retransloc fraction to N_storage at senescence
%   [STOMATA]
PFTs(5).m                      = 7.1017;   % [CLM5 NC audit 2026-05-26] mbbopt, was 9
PFTs(5).alpha                  = 0.06;
%   [LEAF]
PFTs(5).Ztop                   = 0;
PFTs(5).Zbot                   = 0;
PFTs(5).SLA0                   = 0.0190;    % [CLM5 audit 2026-05-26] was 0.012
PFTs(5).SLAm                   = 0.0015;
PFTs(5).CN_L                   = 29.6000;    % [CLM5 audit 2026-05-26] was 30
%   [PHOTOSYNTH]
PFTs(5).Vcmax25                = 61.5;
PFTs(5).g_m                    = 0.101;  % [mol m-2 s-1] gm25 (Knauer 2019 GCB25 Table 3 EBF, evergreen broadleaf)
%   [LUNA]
PFTs(5).Jmaxb1                 = 0.0936;
%   [ALLOC_CN]
PFTs(5).CN_target_root         = 42.00;    % [CLM5 NC audit 2026-05-26] was 40
PFTs(5).CN_target_stem         = 50.00;    % [CLM5 NC audit 2026-05-26] was 150
% [13c dual stem pool] CLM5 CN_Allocation Tbl col7-8 (CN_lw/CN_dw); clm5_params (flivewd/lwtop_ann)
PFTs(5).CN_lw                  = 50;    % [gC/gN] live stem (respires)
PFTs(5).CN_dw                  = 500;   % [gC/gN] dead stem (heartwood, ~N-free)
PFTs(5).flivewd                = 0.1;  % [-] new-stem-C frac to live wood (tree/shrub .1/grass 0/crop 1)
% croot_stem: coarse-root:stem allocation ratio -> routes coarse-root share of stem
%   allocation into belowground (root) reserve at storage (CLM5
%   NutrientCompetitionFlexibleCNMod.F90:529, natural tree/shrub value 0.30; clm5_params
%   .c171117__FruitTree.nc verified 2026-06-07). Consumed in computeCarbonAllocation.
PFTs(5).croot_stem             = 0.30;  % [-] coarse-root:stem ratio (CLM5 natural tree/shrub)
PFTs(5).froot_leaf             = 1.5;   % [-] f1 fine-root:leaf  (c_allometry alloc; CLM5 nc idx 5)
PFTs(5).croot_stem_alloc       = 0.30;  % [-] f2 coarse:stem ALLOC flux (c_allometry; vs reserve croot_stem)
PFTs(5).stem_leaf              = 1.5;   % [-] f3 stem:leaf  (c_allometry alloc; CLM5 nc idx 5)
PFTs(5).lwtop                  = 0.7;   % [1/yr] live->dead wood turnover (CLM5 lwtop_ann)
PFTs(5).CN_target_fruit        = 0;
% [dead-value guard 2026-06-28 Jizixing Zhu] evergreen woody (allometry + has_storage_pool=0):
%   f_alloc_* never read (allometry split + no spring onset-xfer); NaN so any accidental read surfaces.
PFTs(5).f_alloc_leaf           = NaN;   % [-] DEAD: allometry, not f_alloc
PFTs(5).f_alloc_root           = NaN;   % [-] DEAD
PFTs(5).f_alloc_stem           = NaN;   % [-] DEAD
PFTs(5).f_alloc_fruit          = 0;     % [-] no fruit (has_fruit=0)
%   [N_UPTAKE] (Nacry kinetics)
PFTs(5).Vmax_NH4_umol_g_h      = 8;
PFTs(5).Vmax_NO3_umol_g_h      = 4;
PFTs(5).Km_NH4                 = 0.3;
PFTs(5).Km_NO3                 = 0.3;
PFTs(5).Cmin_NH4               = 0.5;
PFTs(5).Cmin_NO3               = 0.25;
PFTs(5).SRL                    = 20;
PFTs(5).Vmax_NO3_LATS_umol_g_h = 0;
PFTs(5).Km_LATS_NO3            = 70;
%   [LEAF_LIFESPAN]
PFTs(5).leaf_long_days         = 511.91;    % [CLM5 NC audit 2026-05-26] was 400
PFTs(5).retransl_fraction      = 0.4;
%   [RESPIRATION] Plant respiration (Sitch 2003 GCB 9:161 Eq 22-23, p.169)
PFTs(5).r_m10                = 0.066;     % [g C/(g N·d)] @ 10°C base, Sitch 2003 Table 3 "Plant Respiration" –TeBE –Sitch Table 3 default (temperate broadleaf evergreen)
PFTs(5).R_g_fraction         = 0.11;      % [-] growth resp factor (CLM5 grperc; was Sitch 0.25).
PFTs(5).harvest_fraction         = 0.90;
PFTs(5).root_turnover_rate       = 0.00274;
PFTs(5).spring_release_days      = 60;

% ==============================================================================
% PFT 6: BDT tropical –Broadleaf deciduous tree –tropical
% ==============================================================================
PFTs(6).PFT                    = 'BDT tropical';
PFTs(6).FullName               = 'Broadleaf deciduous tree –tropical';
PFTs(6).clumping_index         = 0.70;  % [-] ClimaLand pfts.jl:235-242, tropical broadleaf deciduous
%   [HABIT] (PFT lifeform flags)
PFTs(6).photo_type             = 'C3';
PFTs(6).is_perennial           = 1;
PFTs(6).has_fruit              = 0;
PFTs(6).has_storage_pool       = 1;
%   [PHENOLOGY] (see [N.0]-[N.5] sub-stages below)
PFTs(6).evergreen            = 0;
PFTs(6).season_decid         = 0;
PFTs(6).stress_decid         = 1;
PFTs(6).crop_phenology       = 0;   % BDT tropical
%   [N.1] Sykes 1996 chilling-GDD (NaN: not season_decid, physics skips)
PFTs(6).k_chilla            = NaN;
PFTs(6).k_chillb            = NaN;
PFTs(6).k_chillk            = NaN;
PFTs(6).phengdd5ramp        = NaN;
%   [N.2] Senescence (CLM5 daylength trigger, CNPhenologyMod.F90:957)
PFTs(6).crit_dayl_h         = 10.9;    % [h] dayl threshold for offset trigger (CLM5 global default)
PFTs(6).ndays_on            = 30;      % [days] CLM5 onset ramp duration (CNPhenologyMod.F90:51, clm5_params.c171117.nc scalar)
PFTs(6).ndays_off           = 15;      % [days] CLM5 offset ramp duration (CNPhenologyMod.F90:52, clm5_params.c171117.nc scalar)
PFTs(6).fstor2tran          = 0.5;     % [-] CLM5 budburst fraction (CNPhenologyMod.F90:53, clm5_params.c171117.nc scalar = 0.5)
PFTs(6).GDD_bloom_after_budburst = NaN;     % [°C·d] only Apple PFT 23 has fruit phenology calibrated for Changwu
PFTs(6).TSum_fruit              = NaN;
PFTs(6).harvest_DS              = NaN;
%   [N.3] N retranslocation on tissue turnover (LPJ-GUESS growth.cpp:256-285)
PFTs(6).nrelocfrac          = 0.5;     % [-] N retransloc fraction to N_storage at senescence
%   [STOMATA]
PFTs(6).m                      = 7.1797;   % [CLM5 NC audit 2026-05-26] mbbopt, was 9
PFTs(6).alpha                  = 0.06;
%   [LEAF]
PFTs(6).Ztop                   = 0;
PFTs(6).Zbot                   = 0;
PFTs(6).SLA0                   = 0.0308;    % [CLM5 audit 2026-05-26] was 0.03
PFTs(6).SLAm                   = 0.0027;    % [CLM5 NC audit 2026-05-26] was 0.004
PFTs(6).CN_L                   = 23.4500;    % [CLM5 audit 2026-05-26] was 25
%   [PHOTOSYNTH]
PFTs(6).Vcmax25                = 41;
PFTs(6).g_m                    = 0.058;  % [mol m-2 s-1] gm25 (Knauer 2019 GCB25 Table 3 TDF, tropical deciduous broadleaf)
%   [LUNA]
PFTs(6).Jmaxb1                 = 0.0936;
%   [ALLOC_CN]
PFTs(6).CN_target_root         = 42.00;    % [CLM5 NC audit 2026-05-26] was 35
PFTs(6).CN_target_stem         = 50.00;    % [CLM5 NC audit 2026-05-26] was 100
% [13c dual stem pool] CLM5 CN_Allocation Tbl col7-8 (CN_lw/CN_dw); clm5_params (flivewd/lwtop_ann)
PFTs(6).CN_lw                  = 50;    % [gC/gN] live stem (respires)
PFTs(6).CN_dw                  = 500;   % [gC/gN] dead stem (heartwood, ~N-free)
PFTs(6).flivewd                = 0.1;  % [-] new-stem-C frac to live wood (tree/shrub .1/grass 0/crop 1)
% croot_stem: coarse-root:stem allocation ratio -> routes coarse-root share of stem
%   allocation into belowground (root) reserve at storage (CLM5
%   NutrientCompetitionFlexibleCNMod.F90:529, natural tree/shrub value 0.30; clm5_params
%   .c171117__FruitTree.nc verified 2026-06-07). Consumed in computeCarbonAllocation.
PFTs(6).croot_stem             = 0.30;  % [-] coarse-root:stem ratio (CLM5 natural tree/shrub)
PFTs(6).froot_leaf             = 1.5;   % [-] f1 fine-root:leaf  (c_allometry alloc; CLM5 nc idx 6)
PFTs(6).croot_stem_alloc       = 0.30;  % [-] f2 coarse:stem ALLOC flux (c_allometry; vs reserve croot_stem)
PFTs(6).stem_leaf              = 1.0;   % [-] f3 stem:leaf  (c_allometry alloc; CLM5 nc idx 6)
PFTs(6).lwtop                  = 0.7;   % [1/yr] live->dead wood turnover (CLM5 lwtop_ann)
PFTs(6).CN_target_fruit        = 0;
PFTs(6).f_alloc_leaf           = 0.4;
PFTs(6).f_alloc_root           = 0.25;
PFTs(6).f_alloc_stem           = 0.35;
PFTs(6).f_alloc_fruit          = 0;
%   [N_UPTAKE] (Nacry kinetics)
PFTs(6).Vmax_NH4_umol_g_h      = 12;
PFTs(6).Vmax_NO3_umol_g_h      = 5;
PFTs(6).Km_NH4                 = 0.3;
PFTs(6).Km_NO3                 = 0.3;
PFTs(6).Cmin_NH4               = 0.5;
PFTs(6).Cmin_NO3               = 0.25;
PFTs(6).SRL                    = 40;
PFTs(6).Vmax_NO3_LATS_umol_g_h = 0;
PFTs(6).Km_LATS_NO3            = 70;
%   [LEAF_LIFESPAN]
PFTs(6).leaf_long_days         = 176.42;    % [CLM5 NC audit 2026-05-26] was 200
PFTs(6).retransl_fraction      = 0.7;
%   [RESPIRATION] Plant respiration (Sitch 2003 GCB 9:161 Eq 22-23, p.169)
PFTs(6).r_m10                = 0.011;     % [g C/(g N·d)] @ 10°C base, Sitch 2003 Table 3 "Plant Respiration" –TrBR –Sitch Table 3 tropical value (broadleaf raingreen)
PFTs(6).R_g_fraction         = 0.11;      % [-] growth resp factor (CLM5 grperc; was Sitch 0.25).
PFTs(6).harvest_fraction         = 0.90;
PFTs(6).root_turnover_rate       = 0.00274;
PFTs(6).spring_release_days      = 60;

% ==============================================================================
% PFT 7: BDT temperate –Broadleaf deciduous tree –temperate
% ==============================================================================
PFTs(7).PFT                    = 'BDT temperate';
PFTs(7).FullName               = 'Broadleaf deciduous tree –temperate';
PFTs(7).clumping_index         = 0.70;  % [-] ClimaLand pfts.jl:263-270, temperate broadleaf deciduous
%   [HABIT] (PFT lifeform flags)
PFTs(7).photo_type             = 'C3';
PFTs(7).is_perennial           = 1;
PFTs(7).has_fruit              = 0;
PFTs(7).has_storage_pool       = 1;
%   [PHENOLOGY] (see [N.0]-[N.5] sub-stages below)
PFTs(7).evergreen            = 0;
PFTs(7).season_decid         = 1;
PFTs(7).stress_decid         = 0;
PFTs(7).crop_phenology       = 0;   % BDT temperate
%   [N.1] Sykes 1996 chilling-GDD for season_decid (LPJ-GUESS growth.cpp:91-97)
PFTs(7).k_chilla            = 0;       % [°C·d] Sykes 1996 a, LPJ-GUESS default
PFTs(7).k_chillb            = 100;     % [°C·d] Sykes 1996 b, species-specific (LPJ-GUESS europe.ins)
PFTs(7).k_chillk            = 0.05;    % [d^-1] Sykes 1996 k, LPJ-GUESS default
PFTs(7).phengdd5ramp        = 200;     % [°C·d] linear ramp to full leaf-on (LPJ-GUESS)
%   [N.2] Senescence (CLM5 daylength trigger, CNPhenologyMod.F90:957)
PFTs(7).crit_dayl_h         = 10.9;    % [h] dayl threshold for offset trigger (CLM5 global default)
PFTs(7).ndays_on            = 30;      % [days] CLM5 onset ramp duration (CNPhenologyMod.F90:51, clm5_params.c171117.nc scalar)
PFTs(7).ndays_off           = 15;      % [days] CLM5 offset ramp duration (CNPhenologyMod.F90:52, clm5_params.c171117.nc scalar)
PFTs(7).fstor2tran          = 0.5;     % [-] CLM5 budburst fraction (CNPhenologyMod.F90:53, clm5_params.c171117.nc scalar = 0.5)
PFTs(7).GDD_bloom_after_budburst = NaN;     % [°C·d] only Apple PFT 23 has fruit phenology calibrated for Changwu
PFTs(7).TSum_fruit              = NaN;
PFTs(7).harvest_DS              = NaN;
%   [N.3] N retranslocation on tissue turnover (LPJ-GUESS growth.cpp:256-285)
PFTs(7).nrelocfrac          = 0.5;     % [-] N retransloc fraction to N_storage at senescence
%   [STOMATA]
PFTs(7).m                      = 7.1017;   % [CLM5 NC audit 2026-05-26] mbbopt, was 9
PFTs(7).alpha                  = 0.06;
%   [LEAF]
PFTs(7).Ztop                   = 0;
PFTs(7).Zbot                   = 0;
PFTs(7).SLA0                   = 0.0308;    % [CLM5 audit 2026-05-26] was 0.03
PFTs(7).SLAm                   = 0.0027;    % [CLM5 NC audit 2026-05-26] was 0.004
PFTs(7).CN_L                   = 23.4500;    % [CLM5 audit 2026-05-26] was 25
%   [PHOTOSYNTH]
PFTs(7).Vcmax25                = 57.7;
PFTs(7).g_m                    = 0.175;  % [mol m-2 s-1] gm25 (Knauer 2019 GCB25 Table 3 DBF, deciduous broadleaf)
%   [LUNA]
PFTs(7).Jmaxb1                 = 0.0936;
%   [ALLOC_CN]
PFTs(7).CN_target_root         = 42.00;    % [CLM5 NC audit 2026-05-26] was 35
PFTs(7).CN_target_stem         = 50.00;    % [CLM5 NC audit 2026-05-26] was 100
% [13c dual stem pool] CLM5 CN_Allocation Tbl col7-8 (CN_lw/CN_dw); clm5_params (flivewd/lwtop_ann)
PFTs(7).CN_lw                  = 50;    % [gC/gN] live stem (respires)
PFTs(7).CN_dw                  = 500;   % [gC/gN] dead stem (heartwood, ~N-free)
PFTs(7).flivewd                = 0.1;  % [-] new-stem-C frac to live wood (tree/shrub .1/grass 0/crop 1)
% croot_stem: coarse-root:stem allocation ratio -> routes coarse-root share of stem
%   allocation into belowground (root) reserve at storage (CLM5
%   NutrientCompetitionFlexibleCNMod.F90:529, natural tree/shrub value 0.30; clm5_params
%   .c171117__FruitTree.nc verified 2026-06-07). Consumed in computeCarbonAllocation.
PFTs(7).croot_stem             = 0.30;  % [-] coarse-root:stem ratio (CLM5 natural tree/shrub)
PFTs(7).froot_leaf             = 1.5;   % [-] f1 fine-root:leaf  (c_allometry alloc; CLM5 nc idx 7)
PFTs(7).croot_stem_alloc       = 0.30;  % [-] f2 coarse:stem ALLOC flux (c_allometry; vs reserve croot_stem)
PFTs(7).stem_leaf              = 2.3;   % [-] f3 stem:leaf  (c_allometry alloc; CLM5 nc idx 7)
PFTs(7).lwtop                  = 0.7;   % [1/yr] live->dead wood turnover (CLM5 lwtop_ann)
PFTs(7).CN_target_fruit        = 0;
PFTs(7).f_alloc_leaf           = 0.4;
PFTs(7).f_alloc_root           = 0.25;
PFTs(7).f_alloc_stem           = 0.35;
PFTs(7).f_alloc_fruit          = 0;
%   [N_UPTAKE] (Nacry kinetics)
PFTs(7).Vmax_NH4_umol_g_h      = 12;
PFTs(7).Vmax_NO3_umol_g_h      = 5;
PFTs(7).Km_NH4                 = 0.3;
PFTs(7).Km_NO3                 = 0.3;
PFTs(7).Cmin_NH4               = 0.5;
PFTs(7).Cmin_NO3               = 0.25;
PFTs(7).SRL                    = 40;
PFTs(7).Vmax_NO3_LATS_umol_g_h = 0;
PFTs(7).Km_LATS_NO3            = 70;
%   [LEAF_LIFESPAN]
PFTs(7).leaf_long_days         = 176.42;    % [CLM5 NC audit 2026-05-26] was 200
PFTs(7).retransl_fraction      = 0.7;
%   [RESPIRATION] Plant respiration (Sitch 2003 GCB 9:161 Eq 22-23, p.169)
PFTs(7).r_m10                = 0.066;     % [g C/(g N·d)] @ 10°C base, Sitch 2003 Table 3 "Plant Respiration" –TeBS –Sitch Table 3 default (temperate broadleaf summergreen)
PFTs(7).R_g_fraction         = 0.11;      % [-] growth resp factor (CLM5 grperc; was Sitch 0.25).
PFTs(7).harvest_fraction         = 0.90;
PFTs(7).root_turnover_rate       = 0.00274;
PFTs(7).spring_release_days      = 60;

% ==============================================================================
% PFT 8: BDT boreal –Broadleaf deciduous tree –boreal
% ==============================================================================
PFTs(8).PFT                    = 'BDT boreal';
PFTs(8).FullName               = 'Broadleaf deciduous tree –boreal';
PFTs(8).clumping_index         = 0.70;  % [-] ClimaLand pfts.jl:291-298, boreal broadleaf deciduous
%   [HABIT] (PFT lifeform flags)
PFTs(8).photo_type             = 'C3';
PFTs(8).is_perennial           = 1;
PFTs(8).has_fruit              = 0;
PFTs(8).has_storage_pool       = 1;
%   [PHENOLOGY] (see [N.0]-[N.5] sub-stages below)
PFTs(8).evergreen            = 0;
PFTs(8).season_decid         = 1;
PFTs(8).stress_decid         = 0;
PFTs(8).crop_phenology       = 0;   % BDT boreal
%   [N.1] Sykes 1996 chilling-GDD for season_decid (LPJ-GUESS growth.cpp:91-97)
PFTs(8).k_chilla            = 0;       % [°C·d] Sykes 1996 a, LPJ-GUESS default
PFTs(8).k_chillb            = 350;     % [°C·d] Sykes 1996 b, species-specific (LPJ-GUESS europe.ins)
PFTs(8).k_chillk            = 0.05;    % [d^-1] Sykes 1996 k, LPJ-GUESS default
PFTs(8).phengdd5ramp        = 200;     % [°C·d] linear ramp to full leaf-on (LPJ-GUESS)
%   [N.2] Senescence (CLM5 daylength trigger, CNPhenologyMod.F90:957)
PFTs(8).crit_dayl_h         = 10.9;    % [h] dayl threshold for offset trigger (CLM5 global default)
PFTs(8).ndays_on            = 30;      % [days] CLM5 onset ramp duration (CNPhenologyMod.F90:51, clm5_params.c171117.nc scalar)
PFTs(8).ndays_off           = 15;      % [days] CLM5 offset ramp duration (CNPhenologyMod.F90:52, clm5_params.c171117.nc scalar)
PFTs(8).fstor2tran          = 0.5;     % [-] CLM5 budburst fraction (CNPhenologyMod.F90:53, clm5_params.c171117.nc scalar = 0.5)
PFTs(8).GDD_bloom_after_budburst = NaN;     % [°C·d] only Apple PFT 23 has fruit phenology calibrated for Changwu
PFTs(8).TSum_fruit              = NaN;
PFTs(8).harvest_DS              = NaN;
%   [N.3] N retranslocation on tissue turnover (LPJ-GUESS growth.cpp:256-285)
PFTs(8).nrelocfrac          = 0.5;     % [-] N retransloc fraction to N_storage at senescence
%   [STOMATA]
PFTs(8).m                      = 7.1017;   % [CLM5 NC audit 2026-05-26] mbbopt, was 9
PFTs(8).alpha                  = 0.06;
%   [LEAF]
PFTs(8).Ztop                   = 0;
PFTs(8).Zbot                   = 0;
PFTs(8).SLA0                   = 0.0308;    % [CLM5 audit 2026-05-26] was 0.03
PFTs(8).SLAm                   = 0.0027;    % [CLM5 NC audit 2026-05-26] was 0.004
PFTs(8).CN_L                   = 23.4500;    % [CLM5 audit 2026-05-26] was 25
%   [PHOTOSYNTH]
PFTs(8).Vcmax25                = 57.7;
PFTs(8).g_m                    = 0.175;  % [mol m-2 s-1] gm25 (Knauer 2019 GCB25 Table 3 DBF, deciduous broadleaf)
%   [LUNA]
PFTs(8).Jmaxb1                 = 0.0936;
%   [ALLOC_CN]
PFTs(8).CN_target_root         = 42.00;    % [CLM5 NC audit 2026-05-26] was 40
PFTs(8).CN_target_stem         = 50.00;    % [CLM5 NC audit 2026-05-26] was 150
% [13c dual stem pool] CLM5 CN_Allocation Tbl col7-8 (CN_lw/CN_dw); clm5_params (flivewd/lwtop_ann)
PFTs(8).CN_lw                  = 50;    % [gC/gN] live stem (respires)
PFTs(8).CN_dw                  = 500;   % [gC/gN] dead stem (heartwood, ~N-free)
PFTs(8).flivewd                = 0.1;  % [-] new-stem-C frac to live wood (tree/shrub .1/grass 0/crop 1)
% croot_stem: coarse-root:stem allocation ratio -> routes coarse-root share of stem
%   allocation into belowground (root) reserve at storage (CLM5
%   NutrientCompetitionFlexibleCNMod.F90:529, natural tree/shrub value 0.30; clm5_params
%   .c171117__FruitTree.nc verified 2026-06-07). Consumed in computeCarbonAllocation.
PFTs(8).croot_stem             = 0.30;  % [-] coarse-root:stem ratio (CLM5 natural tree/shrub)
PFTs(8).froot_leaf             = 1.5;   % [-] f1 fine-root:leaf  (c_allometry alloc; CLM5 nc idx 8)
PFTs(8).croot_stem_alloc       = 0.30;  % [-] f2 coarse:stem ALLOC flux (c_allometry; vs reserve croot_stem)
PFTs(8).stem_leaf              = 2.3;   % [-] f3 stem:leaf  (c_allometry alloc; CLM5 nc idx 8)
PFTs(8).lwtop                  = 0.7;   % [1/yr] live->dead wood turnover (CLM5 lwtop_ann)
PFTs(8).CN_target_fruit        = 0;
PFTs(8).f_alloc_leaf           = 0.4;
PFTs(8).f_alloc_root           = 0.25;
PFTs(8).f_alloc_stem           = 0.35;
PFTs(8).f_alloc_fruit          = 0;
%   [N_UPTAKE] (Nacry kinetics)
PFTs(8).Vmax_NH4_umol_g_h      = 12;
PFTs(8).Vmax_NO3_umol_g_h      = 5;
PFTs(8).Km_NH4                 = 0.3;
PFTs(8).Km_NO3                 = 0.3;
PFTs(8).Cmin_NH4               = 0.5;
PFTs(8).Cmin_NO3               = 0.25;
PFTs(8).SRL                    = 40;
PFTs(8).Vmax_NO3_LATS_umol_g_h = 0;
PFTs(8).Km_LATS_NO3            = 70;
%   [LEAF_LIFESPAN]
PFTs(8).leaf_long_days         = 176.42;    % [CLM5 NC audit 2026-05-26] was 200
PFTs(8).retransl_fraction      = 0.65;
%   [RESPIRATION] Plant respiration (Sitch 2003 GCB 9:161 Eq 22-23, p.169)
PFTs(8).r_m10                = 0.066;     % [g C/(g N·d)] @ 10°C base, Sitch 2003 Table 3 "Plant Respiration" –BoBS –Sitch Table 3 default (boreal broadleaf summergreen)
PFTs(8).R_g_fraction         = 0.11;      % [-] growth resp factor (CLM5 grperc; was Sitch 0.25).
PFTs(8).harvest_fraction         = 0.90;
PFTs(8).root_turnover_rate       = 0.00274;
PFTs(8).spring_release_days      = 60;

% ==============================================================================
% PFT 9: BES temperate –Broadleaf evergreen shrub –temperate
% ==============================================================================
PFTs(9).PFT                    = 'BES temperate';
PFTs(9).FullName               = 'Broadleaf evergreen shrub –temperate';
PFTs(9).clumping_index         = 0.75;  % [-] ClimaLand pfts.jl:319-326, temperate broadleaf evergreen shrub
%   [HABIT] (PFT lifeform flags)
PFTs(9).photo_type             = 'C3';
PFTs(9).is_perennial           = 1;
PFTs(9).has_fruit              = 0;
PFTs(9).has_storage_pool       = 0;
%   [PHENOLOGY] (see [N.0]-[N.5] sub-stages below)
PFTs(9).evergreen            = 1;
PFTs(9).season_decid         = 0;
PFTs(9).stress_decid         = 0;
PFTs(9).crop_phenology       = 0;   % BES temperate
%   [N.1] Sykes 1996 chilling-GDD (NaN: not season_decid, physics skips)
PFTs(9).k_chilla            = NaN;
PFTs(9).k_chillb            = NaN;
PFTs(9).k_chillk            = NaN;
PFTs(9).phengdd5ramp        = NaN;
%   [N.2] Senescence (CLM5 daylength trigger, CNPhenologyMod.F90:957)
PFTs(9).crit_dayl_h         = 10.9;    % [h] dayl threshold for offset trigger (CLM5 global default)
PFTs(9).ndays_on            = 30;      % [days] CLM5 onset ramp duration (CNPhenologyMod.F90:51, clm5_params.c171117.nc scalar)
PFTs(9).ndays_off           = 15;      % [days] CLM5 offset ramp duration (CNPhenologyMod.F90:52, clm5_params.c171117.nc scalar)
PFTs(9).fstor2tran          = 0.5;     % [-] CLM5 budburst fraction (CNPhenologyMod.F90:53, clm5_params.c171117.nc scalar = 0.5)
PFTs(9).GDD_bloom_after_budburst = NaN;     % [°C·d] only Apple PFT 23 has fruit phenology calibrated for Changwu
PFTs(9).TSum_fruit              = NaN;
PFTs(9).harvest_DS              = NaN;
%   [N.3] N retranslocation on tissue turnover (LPJ-GUESS growth.cpp:256-285)
PFTs(9).nrelocfrac          = 0.0;     % [-] N retransloc fraction to N_storage at senescence
%   [STOMATA]
PFTs(9).m                      = 7.1017;   % [CLM5 NC audit 2026-05-26] mbbopt, was 9
PFTs(9).alpha                  = 0.06;
%   [LEAF]
PFTs(9).Ztop                   = 0.5;
PFTs(9).Zbot                   = 0.1;
PFTs(9).SLA0                   = 0.0180;    % [CLM5 audit 2026-05-26] was 0.012
PFTs(9).SLAm                   = 0;
PFTs(9).CN_L                   = 36.4200;    % [CLM5 audit 2026-05-26] was 30
%   [PHOTOSYNTH]
PFTs(9).Vcmax25                = 61.7;
PFTs(9).g_m                    = 0.101;  % [mol m-2 s-1] gm25 (Knauer 2019 GCB25 Table 3 EBF, evergreen broadleaf; nearest to broadleaf evergreen shrub)
%   [LUNA]
PFTs(9).Jmaxb1                 = 0.0936;
%   [ALLOC_CN]
PFTs(9).CN_target_root         = 42.00;    % [CLM5 NC audit 2026-05-26] was 40
PFTs(9).CN_target_stem         = 50.00;    % [CLM5 NC audit 2026-05-26] was 120
% [13c dual stem pool] CLM5 CN_Allocation Tbl col7-8 (CN_lw/CN_dw); clm5_params (flivewd/lwtop_ann)
PFTs(9).CN_lw                  = 50;    % [gC/gN] live stem (respires)
PFTs(9).CN_dw                  = 500;   % [gC/gN] dead stem (heartwood, ~N-free)
PFTs(9).flivewd                = 0.5;  % [-] f4 live-wood frac; BES temperate shrub (CLM5 nc flivewd, vs .1 trees)
% croot_stem: coarse-root:stem allocation ratio -> routes coarse-root share of stem
%   allocation into belowground (root) reserve at storage (CLM5
%   NutrientCompetitionFlexibleCNMod.F90:529, natural tree/shrub value 0.30; clm5_params
%   .c171117__FruitTree.nc verified 2026-06-07). Consumed in computeCarbonAllocation.
PFTs(9).croot_stem             = 0.30;  % [-] coarse-root:stem ratio (CLM5 natural tree/shrub)
PFTs(9).froot_leaf             = 1.5;   % [-] f1 fine-root:leaf  (c_allometry alloc; CLM5 nc idx 9)
PFTs(9).croot_stem_alloc       = 0.30;  % [-] f2 coarse:stem ALLOC flux (c_allometry; vs reserve croot_stem)
PFTs(9).stem_leaf              = 1.4;   % [-] f3 stem:leaf  (c_allometry alloc; CLM5 nc idx 9, BES shrub)
PFTs(9).lwtop                  = 0.7;   % [1/yr] live->dead wood turnover (CLM5 lwtop_ann)
PFTs(9).CN_target_fruit        = 0;
% [dead-value guard 2026-06-28 Jizixing Zhu] evergreen woody (allometry + has_storage_pool=0):
%   f_alloc_* never read (allometry split + no spring onset-xfer); NaN so any accidental read surfaces.
PFTs(9).f_alloc_leaf           = NaN;   % [-] DEAD: allometry, not f_alloc
PFTs(9).f_alloc_root           = NaN;   % [-] DEAD
PFTs(9).f_alloc_stem           = NaN;   % [-] DEAD
PFTs(9).f_alloc_fruit          = 0;     % [-] no fruit (has_fruit=0)
%   [N_UPTAKE] (Nacry kinetics)
PFTs(9).Vmax_NH4_umol_g_h      = 10;
PFTs(9).Vmax_NO3_umol_g_h      = 5;
PFTs(9).Km_NH4                 = 0.3;
PFTs(9).Km_NO3                 = 0.3;
PFTs(9).Cmin_NH4               = 0.5;
PFTs(9).Cmin_NO3               = 0.25;
PFTs(9).SRL                    = 35;
PFTs(9).Vmax_NO3_LATS_umol_g_h = 0;
PFTs(9).Km_LATS_NO3            = 70;
%   [LEAF_LIFESPAN]
PFTs(9).leaf_long_days         = 483.02;    % [CLM5 NC audit 2026-05-26] was 500
PFTs(9).retransl_fraction      = 0.4;
%   [RESPIRATION] Plant respiration (Sitch 2003 GCB 9:161 Eq 22-23, p.169)
PFTs(9).r_m10                = 0.066;     % [g C/(g N·d)] @ 10°C base, Sitch 2003 Table 3 "Plant Respiration" –TeBE class (shrub, Sitch lumps shrubs with temperate woody)
PFTs(9).R_g_fraction         = 0.11;      % [-] growth resp factor (CLM5 grperc; was Sitch 0.25).
PFTs(9).harvest_fraction         = 0.90;
PFTs(9).root_turnover_rate       = 0.00274;
PFTs(9).spring_release_days      = 60;

% ==============================================================================
% PFT 10: BDS temperate –Broadleaf deciduous shrub –temperate
% ==============================================================================
PFTs(10).PFT                    = 'BDS temperate';
PFTs(10).FullName               = 'Broadleaf deciduous shrub –temperate';
PFTs(10).clumping_index         = 0.75;  % [-] ClimaLand pfts.jl:347-354, temperate broadleaf deciduous shrub
%   [HABIT] (PFT lifeform flags)
PFTs(10).photo_type             = 'C3';
PFTs(10).is_perennial           = 1;
PFTs(10).has_fruit              = 0;
PFTs(10).has_storage_pool       = 1;
%   [PHENOLOGY] (see [N.0]-[N.5] sub-stages below)
PFTs(10).evergreen            = 0;
PFTs(10).season_decid         = 0;
PFTs(10).stress_decid         = 1;       % [CLM5 NC audit 2026-05-26] was season_decid; CLM5 PFT 10 is stress_decid (clm5_params.c171117.nc verified)
PFTs(10).crop_phenology       = 0;   % BDS temperate
%   [N.1] Sykes 1996 chilling-GDD (NaN: not season_decid, physics skips)
PFTs(10).k_chilla            = NaN;
PFTs(10).k_chillb            = NaN;
PFTs(10).k_chillk            = NaN;
PFTs(10).phengdd5ramp        = NaN;
%   [N.2] Senescence (CLM5 daylength trigger, CNPhenologyMod.F90:957)
PFTs(10).crit_dayl_h         = 10.9;    % [h] dayl threshold for offset trigger (CLM5 global default)
PFTs(10).ndays_on            = 30;      % [days] CLM5 onset ramp duration (CNPhenologyMod.F90:51, clm5_params.c171117.nc scalar)
PFTs(10).ndays_off           = 15;      % [days] CLM5 offset ramp duration (CNPhenologyMod.F90:52, clm5_params.c171117.nc scalar)
PFTs(10).fstor2tran          = 0.5;     % [-] CLM5 budburst fraction (CNPhenologyMod.F90:53, clm5_params.c171117.nc scalar = 0.5)
PFTs(10).GDD_bloom_after_budburst = NaN;     % [°C·d] only Apple PFT 23 has fruit phenology calibrated for Changwu
PFTs(10).TSum_fruit              = NaN;
PFTs(10).harvest_DS              = NaN;
%   [N.3] N retranslocation on tissue turnover (LPJ-GUESS growth.cpp:256-285)
PFTs(10).nrelocfrac          = 0.5;     % [-] N retransloc fraction to N_storage at senescence
%   [STOMATA]
PFTs(10).m                      = 7.1017;   % [CLM5 NC audit 2026-05-26] mbbopt, was 9
PFTs(10).alpha                  = 0.06;
%   [LEAF]
PFTs(10).Ztop                   = 0.5;
PFTs(10).Zbot                   = 0.1;
PFTs(10).SLA0                   = 0.0307;    % [CLM5 audit 2026-05-26] was 0.03
PFTs(10).SLAm                   = 0;
PFTs(10).CN_L                   = 23.2600;    % [CLM5 audit 2026-05-26] was 25
%   [PHOTOSYNTH]
PFTs(10).Vcmax25                = 54;
PFTs(10).g_m                    = 0.098;  % [mol m-2 s-1] gm25 (Knauer 2019 GCB25 Table 3 DSH, deciduous shrub)
%   [LUNA]
PFTs(10).Jmaxb1                 = 0.0936;
%   [ALLOC_CN]
PFTs(10).CN_target_root         = 42.00;    % [CLM5 NC audit 2026-05-26] was 40
PFTs(10).CN_target_stem         = 50.00;    % [CLM5 NC audit 2026-05-26] was 120
% [13c dual stem pool] CLM5 CN_Allocation Tbl col7-8 (CN_lw/CN_dw); clm5_params (flivewd/lwtop_ann)
PFTs(10).CN_lw                  = 50;    % [gC/gN] live stem (respires)
PFTs(10).CN_dw                  = 500;   % [gC/gN] dead stem (heartwood, ~N-free)
PFTs(10).flivewd                = 0.5;  % [-] f4 live-wood frac; BDS temperate shrub (CLM5 nc flivewd, vs .1 trees)
% croot_stem: coarse-root:stem allocation ratio -> routes coarse-root share of stem
%   allocation into belowground (root) reserve at storage (CLM5
%   NutrientCompetitionFlexibleCNMod.F90:529, natural tree/shrub value 0.30; clm5_params
%   .c171117__FruitTree.nc verified 2026-06-07). Consumed in computeCarbonAllocation.
PFTs(10).croot_stem             = 0.30;  % [-] coarse-root:stem ratio (CLM5 natural tree/shrub)
PFTs(10).froot_leaf             = 1.5;   % [-] f1 fine-root:leaf  (c_allometry alloc; CLM5 nc idx 10)
PFTs(10).croot_stem_alloc       = 0.30;  % [-] f2 coarse:stem ALLOC flux (c_allometry; vs reserve croot_stem)
PFTs(10).stem_leaf              = 0.24;  % [-] f3 stem:leaf  (c_allometry alloc; CLM5 nc idx 10, BDS shrub)
PFTs(10).lwtop                  = 0.7;   % [1/yr] live->dead wood turnover (CLM5 lwtop_ann)
PFTs(10).CN_target_fruit        = 0;
PFTs(10).f_alloc_leaf           = 0.4;
PFTs(10).f_alloc_root           = 0.3;
PFTs(10).f_alloc_stem           = 0.3;
PFTs(10).f_alloc_fruit          = 0;
%   [N_UPTAKE] (Nacry kinetics)
PFTs(10).Vmax_NH4_umol_g_h      = 10;
PFTs(10).Vmax_NO3_umol_g_h      = 5;
PFTs(10).Km_NH4                 = 0.3;
PFTs(10).Km_NO3                 = 0.3;
PFTs(10).Cmin_NH4               = 0.5;
PFTs(10).Cmin_NO3               = 0.25;
PFTs(10).SRL                    = 35;
PFTs(10).Vmax_NO3_LATS_umol_g_h = 0;
PFTs(10).Km_LATS_NO3            = 70;
%   [LEAF_LIFESPAN]
PFTs(10).leaf_long_days         = 142.35;    % [CLM5 NC audit 2026-05-26] was 200
PFTs(10).retransl_fraction      = 0.65;
%   [RESPIRATION] Plant respiration (Sitch 2003 GCB 9:161 Eq 22-23, p.169)
PFTs(10).r_m10                = 0.066;     % [g C/(g N·d)] @ 10°C base, Sitch 2003 Table 3 "Plant Respiration" –TeBS class (shrub, Sitch lumps shrubs with temperate woody)
PFTs(10).R_g_fraction         = 0.11;      % [-] growth resp factor (CLM5 grperc; was Sitch 0.25).
PFTs(10).harvest_fraction         = 0.90;
PFTs(10).root_turnover_rate       = 0.00274;
PFTs(10).spring_release_days      = 60;

% ==============================================================================
% PFT 11: BDS boreal –Broadleaf deciduous shrub –boreal
% ==============================================================================
PFTs(11).PFT                    = 'BDS boreal';
PFTs(11).FullName               = 'Broadleaf deciduous shrub –boreal';
PFTs(11).clumping_index         = 0.75;  % [-] ClimaLand pfts.jl:375-382, boreal broadleaf deciduous shrub
%   [HABIT] (PFT lifeform flags)
PFTs(11).photo_type             = 'C3';
PFTs(11).is_perennial           = 1;
PFTs(11).has_fruit              = 0;
PFTs(11).has_storage_pool       = 1;
%   [PHENOLOGY] (see [N.0]-[N.5] sub-stages below)
PFTs(11).evergreen            = 0;
PFTs(11).season_decid         = 1;
PFTs(11).stress_decid         = 0;
PFTs(11).crop_phenology       = 0;   % BDS boreal
%   [N.1] Sykes 1996 chilling-GDD for season_decid (LPJ-GUESS growth.cpp:91-97)
PFTs(11).k_chilla            = 0;       % [°C·d] Sykes 1996 a, LPJ-GUESS default
PFTs(11).k_chillb            = 350;     % [°C·d] Sykes 1996 b, species-specific (LPJ-GUESS europe.ins)
PFTs(11).k_chillk            = 0.05;    % [d^-1] Sykes 1996 k, LPJ-GUESS default
PFTs(11).phengdd5ramp        = 200;     % [°C·d] linear ramp to full leaf-on (LPJ-GUESS)
%   [N.2] Senescence (CLM5 daylength trigger, CNPhenologyMod.F90:957)
PFTs(11).crit_dayl_h         = 10.9;    % [h] dayl threshold for offset trigger (CLM5 global default)
PFTs(11).ndays_on            = 30;      % [days] CLM5 onset ramp duration (CNPhenologyMod.F90:51, clm5_params.c171117.nc scalar)
PFTs(11).ndays_off           = 15;      % [days] CLM5 offset ramp duration (CNPhenologyMod.F90:52, clm5_params.c171117.nc scalar)
PFTs(11).fstor2tran          = 0.5;     % [-] CLM5 budburst fraction (CNPhenologyMod.F90:53, clm5_params.c171117.nc scalar = 0.5)
PFTs(11).GDD_bloom_after_budburst = NaN;     % [°C·d] only Apple PFT 23 has fruit phenology calibrated for Changwu
PFTs(11).TSum_fruit              = NaN;
PFTs(11).harvest_DS              = NaN;
%   [N.3] N retranslocation on tissue turnover (LPJ-GUESS growth.cpp:256-285)
PFTs(11).nrelocfrac          = 0.5;     % [-] N retransloc fraction to N_storage at senescence
%   [STOMATA]
PFTs(11).m                      = 7.1017;   % [CLM5 NC audit 2026-05-26] mbbopt, was 9
PFTs(11).alpha                  = 0.06;
%   [LEAF]
PFTs(11).Ztop                   = 0.5;
PFTs(11).Zbot                   = 0.1;
PFTs(11).SLA0                   = 0.0307;    % [CLM5 audit 2026-05-26] was 0.03
PFTs(11).SLAm                   = 0;
PFTs(11).CN_L                   = 23.2600;    % [CLM5 audit 2026-05-26] was 25
%   [PHOTOSYNTH]
PFTs(11).Vcmax25                = 54;
PFTs(11).g_m                    = 0.098;  % [mol m-2 s-1] gm25 (Knauer 2019 GCB25 Table 3 DSH, deciduous shrub)
%   [LUNA]
PFTs(11).Jmaxb1                 = 0.0936;
%   [ALLOC_CN]
PFTs(11).CN_target_root         = 42.00;    % [CLM5 NC audit 2026-05-26] was 40
PFTs(11).CN_target_stem         = 50.00;    % [CLM5 NC audit 2026-05-26] was 120
% [13c dual stem pool] CLM5 CN_Allocation Tbl col7-8 (CN_lw/CN_dw); clm5_params (flivewd/lwtop_ann)
PFTs(11).CN_lw                  = 50;    % [gC/gN] live stem (respires)
PFTs(11).CN_dw                  = 500;   % [gC/gN] dead stem (heartwood, ~N-free)
PFTs(11).flivewd                = 0.1;  % [-] new-stem-C frac to live wood (tree/shrub .1/grass 0/crop 1)
% croot_stem: coarse-root:stem allocation ratio -> routes coarse-root share of stem
%   allocation into belowground (root) reserve at storage (CLM5
%   NutrientCompetitionFlexibleCNMod.F90:529, natural tree/shrub value 0.30; clm5_params
%   .c171117__FruitTree.nc verified 2026-06-07). Consumed in computeCarbonAllocation.
PFTs(11).croot_stem             = 0.30;  % [-] coarse-root:stem ratio (CLM5 natural tree/shrub)
PFTs(11).froot_leaf             = 1.5;   % [-] f1 fine-root:leaf  (c_allometry alloc; CLM5 nc idx 11)
PFTs(11).croot_stem_alloc       = 0.30;  % [-] f2 coarse:stem ALLOC flux (c_allometry; vs reserve croot_stem)
PFTs(11).stem_leaf              = 0.24;  % [-] f3 stem:leaf  (c_allometry alloc; CLM5 nc idx 11, BDS shrub)
PFTs(11).lwtop                  = 0.7;   % [1/yr] live->dead wood turnover (CLM5 lwtop_ann)
PFTs(11).CN_target_fruit        = 0;
PFTs(11).f_alloc_leaf           = 0.4;
PFTs(11).f_alloc_root           = 0.3;
PFTs(11).f_alloc_stem           = 0.3;
PFTs(11).f_alloc_fruit          = 0;
%   [N_UPTAKE] (Nacry kinetics)
PFTs(11).Vmax_NH4_umol_g_h      = 10;
PFTs(11).Vmax_NO3_umol_g_h      = 5;
PFTs(11).Km_NH4                 = 0.3;
PFTs(11).Km_NO3                 = 0.3;
PFTs(11).Cmin_NH4               = 0.5;
PFTs(11).Cmin_NO3               = 0.25;
PFTs(11).SRL                    = 35;
PFTs(11).Vmax_NO3_LATS_umol_g_h = 0;
PFTs(11).Km_LATS_NO3            = 70;
%   [LEAF_LIFESPAN]
PFTs(11).leaf_long_days         = 142.35;    % [CLM5 NC audit 2026-05-26] was 200
PFTs(11).retransl_fraction      = 0.65;
%   [RESPIRATION] Plant respiration (Sitch 2003 GCB 9:161 Eq 22-23, p.169)
PFTs(11).r_m10                = 0.066;     % [g C/(g N·d)] @ 10°C base, Sitch 2003 Table 3 "Plant Respiration" –BoBS class (shrub, Sitch lumps shrubs with boreal woody)
PFTs(11).R_g_fraction         = 0.11;      % [-] growth resp factor (CLM5 grperc; was Sitch 0.25).
PFTs(11).harvest_fraction         = 0.90;
PFTs(11).root_turnover_rate       = 0.00274;
PFTs(11).spring_release_days      = 60;

% ==============================================================================
% PFT 12: C3 arctic grass –C3 arctic grass
% ==============================================================================
PFTs(12).PFT                    = 'C3 arctic grass';
PFTs(12).FullName               = 'C3 arctic grass';
PFTs(12).clumping_index         = 0.76;  % [-] ClimaLand pfts.jl:403-410, arctic C3 grass
%   [HABIT] (PFT lifeform flags)
PFTs(12).photo_type             = 'C3';
PFTs(12).is_perennial           = 1;
PFTs(12).has_fruit              = 0;
PFTs(12).has_storage_pool       = 0;
%   [PHENOLOGY] (see [N.0]-[N.5] sub-stages below)
PFTs(12).evergreen            = 0;
PFTs(12).season_decid         = 1;       % [CLM5 NC audit 2026-05-26] was stress_decid; CLM5 PFT 12 is season_decid (clm5_params.c171117.nc verified)
PFTs(12).stress_decid         = 0;
PFTs(12).crop_phenology       = 0;   % C3 arctic grass
%   [N.1] Sykes 1996 chilling-GDD for season_decid (LPJ-GUESS growth.cpp:91-97)
PFTs(12).k_chilla            = 0;       % [°C·d] Sykes 1996 a, LPJ-GUESS default
PFTs(12).k_chillb            = 350;     % [°C·d] Sykes 1996 b, arctic grass –high chilling (
PFTs(12).k_chillk            = 0.05;    % [d^-1] Sykes 1996 k, LPJ-GUESS default
PFTs(12).phengdd5ramp        = 200;     % [°C·d] linear ramp to full leaf-on (LPJ-GUESS)
%   [N.2] Senescence (CLM5 daylength trigger, CNPhenologyMod.F90:957)
PFTs(12).crit_dayl_h         = 10.9;    % [h] dayl threshold for offset trigger (CLM5 global default)
PFTs(12).ndays_on            = 30;      % [days] CLM5 onset ramp duration (CNPhenologyMod.F90:51, clm5_params.c171117.nc scalar)
PFTs(12).ndays_off           = 15;      % [days] CLM5 offset ramp duration (CNPhenologyMod.F90:52, clm5_params.c171117.nc scalar)
PFTs(12).fstor2tran          = 0.5;     % [-] CLM5 budburst fraction (CNPhenologyMod.F90:53, clm5_params.c171117.nc scalar = 0.5)
PFTs(12).GDD_bloom_after_budburst = NaN;     % [°C·d] only Apple PFT 23 has fruit phenology calibrated for Changwu
PFTs(12).TSum_fruit              = NaN;
PFTs(12).harvest_DS              = NaN;
%   [N.3] N retranslocation on tissue turnover (LPJ-GUESS growth.cpp:256-285)
PFTs(12).nrelocfrac          = 0.5;     % [-] N retransloc fraction to N_storage at senescence
%   [STOMATA]
PFTs(12).m                      = 7.1017;   % [CLM5 NC audit 2026-05-26] mbbopt, was 9
PFTs(12).alpha                  = 0.06;
%   [LEAF]
PFTs(12).Ztop                   = 0.5;
PFTs(12).Zbot                   = 0.01;
PFTs(12).SLA0                   = 0.0402;    % [CLM5 audit 2026-05-26] was 0.03
PFTs(12).SLAm                   = 0;
PFTs(12).CN_L                   = 28.0300;    % [CLM5 audit 2026-05-26] was 25
%   [PHOTOSYNTH]
PFTs(12).Vcmax25                = 78.2;
PFTs(12).g_m                    = 0.197;  % [mol m-2 s-1] gm25 (Knauer 2019 GCB25 Table 3 C3G, C3 grass/herb)
%   [LUNA]
PFTs(12).Jmaxb1                 = 0.0936;
%   [ALLOC_CN]
PFTs(12).CN_target_root         = 42.00;    % [CLM5 NC audit 2026-05-26] was 40
PFTs(12).CN_target_stem         = 0;
% [13c dual stem pool] CLM5 CN_Allocation Tbl col7-8 (CN_lw/CN_dw); clm5_params (flivewd/lwtop_ann)
PFTs(12).CN_lw                  = 50;    % [gC/gN] live stem (respires)
PFTs(12).CN_dw                  = 500;   % [gC/gN] dead stem (heartwood, ~N-free)
PFTs(12).flivewd                = 0.0;  % [-] new-stem-C frac to live wood (tree/shrub .1/grass 0/crop 1)
PFTs(12).lwtop                  = 0.7;   % [1/yr] live->dead wood turnover (CLM5 lwtop_ann)
PFTs(12).CN_target_fruit        = 0;
% [CLM5 grass allometry 2026-06-28 Jizixing Zhu] grass routes C+N through the c_allometry branch
%   (froot_leaf=1.5/stem_leaf=0 -> S=2.5 -> leaf 0.4/root 0.6; clm5_params.c171117 c3/c4 grass,
%   root-dominant), so f_alloc_* is never read for this PFT. Marked NaN (not a 0.6/0.4 split) so any
%   accidental f_alloc read on grass surfaces immediately. See buildPlantCNParams herbAllometry guard.
% Change: 2026-06-28 Jizixing Zhu — f_alloc_leaf/root/stem 0.6/0.4/0 -> NaN (dead-value guard).
PFTs(12).f_alloc_leaf           = NaN;   % [-] DEAD: grass uses allometry, not f_alloc
PFTs(12).f_alloc_root           = NaN;   % [-] DEAD
PFTs(12).f_alloc_stem           = NaN;   % [-] DEAD
PFTs(12).f_alloc_fruit          = 0;     % [-] grass has no fruit (has_fruit=0)
PFTs(12).froot_leaf             = 1.5;   % [-] f1 fine-root:leaf (CLM5 grass)
PFTs(12).stem_leaf              = 0;     % [-] f3 stem:leaf (grass: no wood)
PFTs(12).croot_stem             = 0;     % [-] coarse-root:stem (grass: none)
PFTs(12).croot_stem_alloc       = 0;     % [-] f2 coarse:stem alloc flux (grass: none)
%   [N_UPTAKE] (Nacry kinetics)
PFTs(12).Vmax_NH4_umol_g_h      = 15;
PFTs(12).Vmax_NO3_umol_g_h      = 10;
PFTs(12).Km_NH4                 = 0.3;
PFTs(12).Km_NO3                 = 0.3;
PFTs(12).Cmin_NH4               = 0.5;
PFTs(12).Cmin_NO3               = 0.25;
PFTs(12).SRL                    = 60;
PFTs(12).Vmax_NO3_LATS_umol_g_h = 0;
PFTs(12).Km_LATS_NO3            = 70;
%   [LEAF_LIFESPAN]
PFTs(12).leaf_long_days         = 117.10;    % [CLM5 NC audit 2026-05-26] was 200
PFTs(12).retransl_fraction      = 0.5;
%   [RESPIRATION] Plant respiration (Sitch 2003 GCB 9:161 Eq 22-23, p.169)
PFTs(12).r_m10                = 0.066;     % [g C/(g N·d)] @ 10°C base, Sitch 2003 Table 3 "Plant Respiration" –TeH –Sitch Table 3 default (temperate herbaceous, C3 arctic)
PFTs(12).R_g_fraction         = 0.11;      % [-] growth resp factor (CLM5 grperc; was Sitch 0.25).
PFTs(12).harvest_fraction         = 0.90;
PFTs(12).root_turnover_rate       = 0.00274;
PFTs(12).spring_release_days      = 60;

% ==============================================================================
% PFT 13: C3 grass –C3 grass
% ==============================================================================
PFTs(13).PFT                    = 'C3 grass';
PFTs(13).FullName               = 'C3 grass';
PFTs(13).clumping_index         = 0.76;  % [-] ClimaLand pfts.jl:431-438, non-arctic C3 grass
%   [HABIT] (PFT lifeform flags)
PFTs(13).photo_type             = 'C3';
PFTs(13).is_perennial           = 1;
PFTs(13).has_fruit              = 0;
PFTs(13).has_storage_pool       = 0;
%   [PHENOLOGY] (see [N.0]-[N.5] sub-stages below)
PFTs(13).evergreen            = 0;
PFTs(13).season_decid         = 0;
PFTs(13).stress_decid         = 1;
PFTs(13).crop_phenology       = 0;   % C3 grass
%   [N.1] Sykes 1996 chilling-GDD (NaN: not season_decid, physics skips)
PFTs(13).k_chilla            = NaN;
PFTs(13).k_chillb            = NaN;
PFTs(13).k_chillk            = NaN;
PFTs(13).phengdd5ramp        = NaN;
%   [N.2] Senescence (CLM5 daylength trigger, CNPhenologyMod.F90:957)
PFTs(13).crit_dayl_h         = 10.9;    % [h] dayl threshold for offset trigger (CLM5 global default)
PFTs(13).ndays_on            = 30;      % [days] CLM5 onset ramp duration (CNPhenologyMod.F90:51, clm5_params.c171117.nc scalar)
PFTs(13).ndays_off           = 15;      % [days] CLM5 offset ramp duration (CNPhenologyMod.F90:52, clm5_params.c171117.nc scalar)
PFTs(13).fstor2tran          = 0.5;     % [-] CLM5 budburst fraction (CNPhenologyMod.F90:53, clm5_params.c171117.nc scalar = 0.5)
PFTs(13).GDD_bloom_after_budburst = NaN;     % [°C·d] only Apple PFT 23 has fruit phenology calibrated for Changwu
PFTs(13).TSum_fruit              = NaN;
PFTs(13).harvest_DS              = NaN;
%   [N.3] N retranslocation on tissue turnover (LPJ-GUESS growth.cpp:256-285)
PFTs(13).nrelocfrac          = 0.5;     % [-] N retransloc fraction to N_storage at senescence
%   [STOMATA]
PFTs(13).m                      = 9.7575;   % [CLM5 NC audit 2026-05-26] mbbopt, was 9
PFTs(13).alpha                  = 0.06;
%   [LEAF]
PFTs(13).Ztop                   = 0.5;
PFTs(13).Zbot                   = 0.01;
PFTs(13).SLA0                   = 0.0402;    % [CLM5 audit 2026-05-26] was 0.03
PFTs(13).SLAm                   = 0;
PFTs(13).CN_L                   = 28.0300;    % [CLM5 audit 2026-05-26] was 25
%   [PHOTOSYNTH]
PFTs(13).Vcmax25                = 78.2;
PFTs(13).g_m                    = 0.197;  % [mol m-2 s-1] gm25 (Knauer 2019 GCB25 Table 3 C3G, C3 grass/herb)
%   [LUNA]
PFTs(13).Jmaxb1                 = 0.0936;
%   [ALLOC_CN]
PFTs(13).CN_target_root         = 42.00;    % [CLM5 NC audit 2026-05-26] was 40
PFTs(13).CN_target_stem         = 0;
% [13c dual stem pool] CLM5 CN_Allocation Tbl col7-8 (CN_lw/CN_dw); clm5_params (flivewd/lwtop_ann)
PFTs(13).CN_lw                  = 50;    % [gC/gN] live stem (respires)
PFTs(13).CN_dw                  = 500;   % [gC/gN] dead stem (heartwood, ~N-free)
PFTs(13).flivewd                = 0.0;  % [-] new-stem-C frac to live wood (tree/shrub .1/grass 0/crop 1)
PFTs(13).lwtop                  = 0.7;   % [1/yr] live->dead wood turnover (CLM5 lwtop_ann)
PFTs(13).CN_target_fruit        = 0;
% [CLM5 grass allometry 2026-06-28 Jizixing Zhu] grass routes C+N through the c_allometry branch
%   (froot_leaf=1.5/stem_leaf=0 -> leaf 0.4/root 0.6; clm5_params.c171117 c3 grass). f_alloc_* never
%   read here; marked NaN so any accidental f_alloc read on grass surfaces immediately.
% Change: 2026-06-28 Jizixing Zhu — f_alloc_leaf/root/stem 0.6/0.4/0 -> NaN (dead-value guard).
PFTs(13).f_alloc_leaf           = NaN;   % [-] DEAD: grass uses allometry, not f_alloc
PFTs(13).f_alloc_root           = NaN;   % [-] DEAD
PFTs(13).f_alloc_stem           = NaN;   % [-] DEAD
PFTs(13).f_alloc_fruit          = 0;     % [-] grass has no fruit (has_fruit=0)
PFTs(13).froot_leaf             = 1.5;   % [-] f1 fine-root:leaf (CLM5 grass)
PFTs(13).stem_leaf              = 0;     % [-] f3 stem:leaf (no wood)
PFTs(13).croot_stem             = 0;     % [-] coarse-root:stem (none)
PFTs(13).croot_stem_alloc       = 0;     % [-] f2 coarse:stem alloc (none)
%   [N_UPTAKE] (Nacry kinetics)
PFTs(13).Vmax_NH4_umol_g_h      = 15;
PFTs(13).Vmax_NO3_umol_g_h      = 10;
PFTs(13).Km_NH4                 = 0.3;
PFTs(13).Km_NO3                 = 0.3;
PFTs(13).Cmin_NH4               = 0.5;
PFTs(13).Cmin_NO3               = 0.25;
PFTs(13).SRL                    = 60;
PFTs(13).Vmax_NO3_LATS_umol_g_h = 0;
PFTs(13).Km_LATS_NO3            = 70;
%   [LEAF_LIFESPAN]
PFTs(13).leaf_long_days         = 117.10;    % [CLM5 NC audit 2026-05-26] was 200
PFTs(13).retransl_fraction      = 0.5;
%   [RESPIRATION] Plant respiration (Sitch 2003 GCB 9:161 Eq 22-23, p.169)
PFTs(13).r_m10                = 0.066;     % [g C/(g N·d)] @ 10°C base, Sitch 2003 Table 3 "Plant Respiration" –TeH –Sitch Table 3 default (temperate herbaceous, C3)
PFTs(13).R_g_fraction         = 0.11;      % [-] growth resp factor (CLM5 grperc; was Sitch 0.25).
PFTs(13).harvest_fraction         = 0.90;
PFTs(13).root_turnover_rate       = 0.00274;
PFTs(13).spring_release_days      = 60;

% ==============================================================================
% PFT 14: C4 grass –C4 grass
% ==============================================================================
PFTs(14).PFT                    = 'C4 grass';
PFTs(14).FullName               = 'C4 grass';
PFTs(14).clumping_index         = 0.75;  % [-] ClimaLand pfts.jl:459-466, C4 grass
%   [HABIT] (PFT lifeform flags)
PFTs(14).photo_type             = 'C4';
PFTs(14).is_perennial           = 1;
PFTs(14).has_fruit              = 0;
PFTs(14).has_storage_pool       = 0;
%   [PHENOLOGY] (see [N.0]-[N.5] sub-stages below)
PFTs(14).evergreen            = 0;
PFTs(14).season_decid         = 0;
PFTs(14).stress_decid         = 1;
PFTs(14).crop_phenology       = 0;   % C4 grass
%   [N.1] Sykes 1996 chilling-GDD (NaN: not season_decid, physics skips)
PFTs(14).k_chilla            = NaN;
PFTs(14).k_chillb            = NaN;
PFTs(14).k_chillk            = NaN;
PFTs(14).phengdd5ramp        = NaN;
%   [N.2] Senescence (CLM5 daylength trigger, CNPhenologyMod.F90:957)
PFTs(14).crit_dayl_h         = 10.9;    % [h] dayl threshold for offset trigger (CLM5 global default)
PFTs(14).ndays_on            = 30;      % [days] CLM5 onset ramp duration (CNPhenologyMod.F90:51, clm5_params.c171117.nc scalar)
PFTs(14).ndays_off           = 15;      % [days] CLM5 offset ramp duration (CNPhenologyMod.F90:52, clm5_params.c171117.nc scalar)
PFTs(14).fstor2tran          = 0.5;     % [-] CLM5 budburst fraction (CNPhenologyMod.F90:53, clm5_params.c171117.nc scalar = 0.5)
PFTs(14).GDD_bloom_after_budburst = NaN;     % [°C·d] only Apple PFT 23 has fruit phenology calibrated for Changwu
PFTs(14).TSum_fruit              = NaN;
PFTs(14).harvest_DS              = NaN;
%   [N.3] N retranslocation on tissue turnover (LPJ-GUESS growth.cpp:256-285)
PFTs(14).nrelocfrac          = 0.5;     % [-] N retransloc fraction to N_storage at senescence
%   [STOMATA]
PFTs(14).m                      = 9.7575;   % [CLM5 NC audit 2026-05-26] mbbopt, was 4
PFTs(14).alpha                  = 0.05;
%   [LEAF]
PFTs(14).Ztop                   = 0.5;
PFTs(14).Zbot                   = 0.01;
PFTs(14).SLA0                   = 0.0385;    % [CLM5 audit 2026-05-26] was 0.03
PFTs(14).SLAm                   = 0;
PFTs(14).CN_L                   = 35.3600;    % [CLM5 audit 2026-05-26] was 25
%   [PHOTOSYNTH]
PFTs(14).Vcmax25                = 51.6;
PFTs(14).g_m                    = 0.452;  % [mol m-2 s-1] gm25 (Knauer 2019 GCB25 Table 3 C4G, C4 grass)
%   [LUNA]
PFTs(14).Jmaxb1                 = 0.0936;
%   [ALLOC_CN]
PFTs(14).CN_target_root         = 42.00;    % [CLM5 NC audit 2026-05-26] was 40
PFTs(14).CN_target_stem         = 0;
% [13c dual stem pool] CLM5 CN_Allocation Tbl col7-8 (CN_lw/CN_dw); clm5_params (flivewd/lwtop_ann)
PFTs(14).CN_lw                  = 50;    % [gC/gN] live stem (respires)
PFTs(14).CN_dw                  = 500;   % [gC/gN] dead stem (heartwood, ~N-free)
PFTs(14).flivewd                = 0.0;  % [-] new-stem-C frac to live wood (tree/shrub .1/grass 0/crop 1)
PFTs(14).lwtop                  = 0.7;   % [1/yr] live->dead wood turnover (CLM5 lwtop_ann)
PFTs(14).CN_target_fruit        = 0;
% [CLM5 grass allometry 2026-06-28 Jizixing Zhu] grass routes C+N through the c_allometry branch
%   (froot_leaf=1.5/stem_leaf=0 -> leaf 0.4/root 0.6; clm5_params.c171117 c4 grass). f_alloc_* never
%   read here; marked NaN so any accidental f_alloc read on grass surfaces immediately.
% Change: 2026-06-28 Jizixing Zhu — f_alloc_leaf/root/stem 0.6/0.4/0 -> NaN (dead-value guard).
PFTs(14).f_alloc_leaf           = NaN;   % [-] DEAD: grass uses allometry, not f_alloc
PFTs(14).f_alloc_root           = NaN;   % [-] DEAD
PFTs(14).f_alloc_stem           = NaN;   % [-] DEAD
PFTs(14).f_alloc_fruit          = 0;     % [-] grass has no fruit (has_fruit=0)
PFTs(14).froot_leaf             = 1.5;   % [-] f1 fine-root:leaf (CLM5 grass)
PFTs(14).stem_leaf              = 0;     % [-] f3 stem:leaf (no wood)
PFTs(14).croot_stem             = 0;     % [-] coarse-root:stem (none)
PFTs(14).croot_stem_alloc       = 0;     % [-] f2 coarse:stem alloc (none)
%   [N_UPTAKE] (Nacry kinetics)
PFTs(14).Vmax_NH4_umol_g_h      = 15;
PFTs(14).Vmax_NO3_umol_g_h      = 10;
PFTs(14).Km_NH4                 = 0.3;
PFTs(14).Km_NO3                 = 0.3;
PFTs(14).Cmin_NH4               = 0.5;
PFTs(14).Cmin_NO3               = 0.25;
PFTs(14).SRL                    = 60;
PFTs(14).Vmax_NO3_LATS_umol_g_h = 0;
PFTs(14).Km_LATS_NO3            = 70;
%   [LEAF_LIFESPAN]
PFTs(14).leaf_long_days         = 51.1000;    % [CLM5 NC audit 2026-05-26] was 200
PFTs(14).retransl_fraction      = 0.5;
%   [RESPIRATION] Plant respiration (Sitch 2003 GCB 9:161 Eq 22-23, p.169)
PFTs(14).r_m10                = 0.011;     % [g C/(g N·d)] @ 10°C base, Sitch 2003 Table 3 "Plant Respiration" –TrH –Sitch Table 3 tropical herbaceous (C4 in Sitch Eq 21 footnote line 375)
PFTs(14).R_g_fraction         = 0.11;      % [-] growth resp factor (CLM5 grperc; was Sitch 0.25).
PFTs(14).harvest_fraction         = 0.90;
PFTs(14).root_turnover_rate       = 0.00274;
PFTs(14).spring_release_days      = 60;

% ==============================================================================
% PFT 15: Temperate Corn –Temperate Corn (Maize)
% ==============================================================================
PFTs(15).PFT                    = 'Temperate Corn';
PFTs(15).FullName               = 'Temperate Corn (Maize)';
PFTs(15).clumping_index         = 0.75;  % [-] Proxy: ClimaLand C4-grass value; no crop-specific entry in its PFT table
%   [HABIT] (PFT lifeform flags)
PFTs(15).photo_type             = 'C4';
PFTs(15).is_perennial           = 0;
PFTs(15).has_fruit              = 1;
PFTs(15).has_storage_pool       = 0;
%   [PHENOLOGY] (see [N.0]-[N.5] sub-stages below)
PFTs(15).evergreen            = 0;
PFTs(15).season_decid         = 0;
PFTs(15).stress_decid         = 0;
PFTs(15).crop_phenology       = 1;   % Corn
%   [N.1] Sykes 1996 chilling-GDD (NaN: not season_decid, physics skips)
PFTs(15).k_chilla            = NaN;
PFTs(15).k_chillb            = NaN;
PFTs(15).k_chillk            = NaN;
PFTs(15).phengdd5ramp        = NaN;
%   [N.2] Senescence (CLM5 daylength trigger, CNPhenologyMod.F90:957)
PFTs(15).crit_dayl_h         = 10.9;    % [h] dayl threshold for offset trigger (CLM5 global default)
PFTs(15).ndays_on            = 30;      % [days] CLM5 onset ramp duration (CNPhenologyMod.F90:51, clm5_params.c171117.nc scalar)
PFTs(15).ndays_off           = 15;      % [days] CLM5 offset ramp duration (CNPhenologyMod.F90:52, clm5_params.c171117.nc scalar)
PFTs(15).fstor2tran          = 0.5;     % [-] CLM5 budburst fraction (CNPhenologyMod.F90:53, clm5_params.c171117.nc scalar = 0.5)
PFTs(15).GDD_bloom_after_budburst = NaN;     % [°C·d] only Apple PFT 23 has fruit phenology calibrated for Changwu
PFTs(15).TSum_fruit              = NaN;
PFTs(15).harvest_DS              = NaN;
%   [N.3] N retranslocation on tissue turnover (LPJ-GUESS growth.cpp:256-285)
PFTs(15).nrelocfrac          = 0.5;     % [-] N retransloc fraction to N_storage at senescence
%   [STOMATA]
PFTs(15).m                      = 4;
PFTs(15).alpha                  = 0.05;
%   [LEAF]
PFTs(15).Ztop                   = 0.5;
PFTs(15).Zbot                   = 0.01;
PFTs(15).SLA0                   = 0.05;
PFTs(15).SLAm                   = 0;
PFTs(15).CN_L                   = 25;
%   [PHOTOSYNTH]
PFTs(15).Vcmax25                = 100.7;
PFTs(15).g_m                    = 0.739;  % [mol m-2 s-1] gm25 (Knauer 2019 GCB25 Table 3 C4C, C4 crop; temperate corn/maize is C4)
%   [LUNA]
PFTs(15).Jmaxb1                 = 0.0936;
%   [ALLOC_CN]
PFTs(15).CN_target_root         = 42.00;    % [CLM5 NC audit 2026-05-26] was 30
PFTs(15).CN_target_stem         = 50.00;    % [CLM5 NC audit 2026-05-26] was 0
% [13c dual stem pool] CLM5 CN_Allocation Tbl col7-8 (CN_lw/CN_dw); clm5_params (flivewd/lwtop_ann)
PFTs(15).CN_lw                  = 50;    % [gC/gN] live stem (respires)
PFTs(15).CN_dw                  = 500;   % [gC/gN] dead stem (heartwood, ~N-free)
PFTs(15).flivewd                = 1.0;  % [-] new-stem-C frac to live wood (tree/shrub .1/grass 0/crop 1)
PFTs(15).lwtop                  = 0.7;   % [1/yr] live->dead wood turnover (CLM5 lwtop_ann)
PFTs(15).CN_target_fruit        = 15;
PFTs(15).f_alloc_leaf           = 0.45;
PFTs(15).f_alloc_root           = 0.3;
PFTs(15).f_alloc_stem           = 0.15;
PFTs(15).f_alloc_fruit          = 0.1;
%   [N_UPTAKE] (Nacry kinetics)
PFTs(15).Vmax_NH4_umol_g_h      = 20;
PFTs(15).Vmax_NO3_umol_g_h      = 15;
PFTs(15).Km_NH4                 = 0.3;
PFTs(15).Km_NO3                 = 0.3;
PFTs(15).Cmin_NH4               = 0.5;
PFTs(15).Cmin_NO3               = 0.25;
PFTs(15).SRL                    = 50;
PFTs(15).Vmax_NO3_LATS_umol_g_h = 0;
PFTs(15).Km_LATS_NO3            = 70;
%   [LEAF_LIFESPAN]
PFTs(15).leaf_long_days         = 365.00;    % [CLM5 NC audit 2026-05-26] was 90
PFTs(15).retransl_fraction      = 0.4;
%   [RESPIRATION] Plant respiration (Sitch 2003 GCB 9:161 Eq 22-23, p.169)
PFTs(15).r_m10                = 0.011;     % [g C/(g N·d)] @ 10°C base, Sitch 2003 Table 3 "Plant Respiration" –TrH-class –C4 crop (Sitch treats all C4 as TrH=0.011)
PFTs(15).R_g_fraction         = 0.11;      % [-] growth resp factor (CLM5 grperc; was Sitch 0.25).
PFTs(15).harvest_fraction         = 0.90;
PFTs(15).root_turnover_rate       = 0.00274;
PFTs(15).spring_release_days      = 60;

% ==============================================================================
% PFT 16: Spring Wheat –Spring Wheat
% ==============================================================================
PFTs(16).PFT                    = 'Spring Wheat';
PFTs(16).FullName               = 'Spring Wheat';
PFTs(16).clumping_index         = 0.76;  % [-] Proxy: ClimaLand C3-grass value; no crop-specific entry in its PFT table
%   [HABIT] (PFT lifeform flags)
PFTs(16).photo_type             = 'C3';
PFTs(16).is_perennial           = 0;
PFTs(16).has_fruit              = 1;
PFTs(16).has_storage_pool       = 0;
%   [PHENOLOGY] (see [N.0]-[N.5] sub-stages below)
PFTs(16).evergreen            = 0;
PFTs(16).season_decid         = 0;
PFTs(16).stress_decid         = 0;
PFTs(16).crop_phenology       = 1;   % Wheat
%   [N.1] Sykes 1996 chilling-GDD (NaN: not season_decid, physics skips)
PFTs(16).k_chilla            = NaN;
PFTs(16).k_chillb            = NaN;
PFTs(16).k_chillk            = NaN;
PFTs(16).phengdd5ramp        = NaN;
%   [N.2] Senescence (CLM5 daylength trigger, CNPhenologyMod.F90:957)
PFTs(16).crit_dayl_h         = 10.9;    % [h] dayl threshold for offset trigger (CLM5 global default)
PFTs(16).ndays_on            = 30;      % [days] CLM5 onset ramp duration (CNPhenologyMod.F90:51, clm5_params.c171117.nc scalar)
PFTs(16).ndays_off           = 15;      % [days] CLM5 offset ramp duration (CNPhenologyMod.F90:52, clm5_params.c171117.nc scalar)
PFTs(16).fstor2tran          = 0.5;     % [-] CLM5 budburst fraction (CNPhenologyMod.F90:53, clm5_params.c171117.nc scalar = 0.5)
PFTs(16).GDD_bloom_after_budburst = NaN;     % [°C·d] only Apple PFT 23 has fruit phenology calibrated for Changwu
PFTs(16).TSum_fruit              = NaN;
PFTs(16).harvest_DS              = NaN;
%   [N.3] N retranslocation on tissue turnover (LPJ-GUESS growth.cpp:256-285)
PFTs(16).nrelocfrac          = 0.5;     % [-] N retransloc fraction to N_storage at senescence
%   [STOMATA]
PFTs(16).m                      = 9;
PFTs(16).alpha                  = 0.06;
%   [LEAF]
PFTs(16).Ztop                   = 0.5;
PFTs(16).Zbot                   = 0.01;
PFTs(16).SLA0                   = 0.0350;    % [CLM5 audit 2026-05-26] was 0.07
PFTs(16).SLAm                   = 0;
PFTs(16).CN_L                   = 20.00;    % [CLM5 audit 2026-05-26] was 25
%   [PHOTOSYNTH]
PFTs(16).Vcmax25                = 100.7;
PFTs(16).g_m                    = 0.295;  % [mol m-2 s-1] gm25 (Knauer 2019 GCB25 Table 3 C3C, C3 crop; spring wheat is C3)
%   [LUNA]
PFTs(16).Jmaxb1                 = 0.0936;
%   [ALLOC_CN]
PFTs(16).CN_target_root         = 42.00;    % [CLM5 NC audit 2026-05-26] was 30
PFTs(16).CN_target_stem         = 50.00;    % [CLM5 NC audit 2026-05-26] was 0
% [13c dual stem pool] CLM5 CN_Allocation Tbl col7-8 (CN_lw/CN_dw); clm5_params (flivewd/lwtop_ann)
PFTs(16).CN_lw                  = 50;    % [gC/gN] live stem (respires)
PFTs(16).CN_dw                  = 500;   % [gC/gN] dead stem (heartwood, ~N-free)
PFTs(16).flivewd                = 1.0;  % [-] new-stem-C frac to live wood (tree/shrub .1/grass 0/crop 1)
PFTs(16).lwtop                  = 0.7;   % [1/yr] live->dead wood turnover (CLM5 lwtop_ann)
PFTs(16).CN_target_fruit        = 30;
PFTs(16).f_alloc_leaf           = 0.45;
PFTs(16).f_alloc_root           = 0.3;
PFTs(16).f_alloc_stem           = 0.15;
PFTs(16).f_alloc_fruit          = 0.1;
%   [N_UPTAKE] (Nacry kinetics)
PFTs(16).Vmax_NH4_umol_g_h      = 20;
PFTs(16).Vmax_NO3_umol_g_h      = 15;
PFTs(16).Km_NH4                 = 0.3;
PFTs(16).Km_NO3                 = 0.3;
PFTs(16).Cmin_NH4               = 0.5;
PFTs(16).Cmin_NO3               = 0.25;
PFTs(16).SRL                    = 50;
PFTs(16).Vmax_NO3_LATS_umol_g_h = 0;
PFTs(16).Km_LATS_NO3            = 70;
%   [LEAF_LIFESPAN]
PFTs(16).leaf_long_days         = 365.00;    % [CLM5 NC audit 2026-05-26] was 120
PFTs(16).retransl_fraction      = 0.4;
%   [RESPIRATION] Plant respiration (Sitch 2003 GCB 9:161 Eq 22-23, p.169)
PFTs(16).r_m10                = 0.066;     % [g C/(g N·d)] @ 10°C base, Sitch 2003 Table 3 "Plant Respiration" –TeH-class –C3 temperate crop (Sitch TeH=0.066)
PFTs(16).R_g_fraction         = 0.11;      % [-] growth resp factor (CLM5 grperc; was Sitch 0.25).
PFTs(16).harvest_fraction         = 0.90;
PFTs(16).root_turnover_rate       = 0.00274;
PFTs(16).spring_release_days      = 60;

% ==============================================================================
% PFT 17: Temperate Soybean –Temperate Soybean
% ==============================================================================
PFTs(17).PFT                    = 'Temperate Soybean';
PFTs(17).FullName               = 'Temperate Soybean';
PFTs(17).clumping_index         = 0.76;  % [-] Proxy: ClimaLand C3-grass value; no crop-specific entry in its PFT table
%   [HABIT] (PFT lifeform flags)
PFTs(17).photo_type             = 'C3';
PFTs(17).is_perennial           = 0;
PFTs(17).has_fruit              = 1;
PFTs(17).has_storage_pool       = 0;
%   [PHENOLOGY] (see [N.0]-[N.5] sub-stages below)
PFTs(17).evergreen            = 0;
PFTs(17).season_decid         = 0;
PFTs(17).stress_decid         = 0;
PFTs(17).crop_phenology       = 1;   % Soybean
%   [N.1] Sykes 1996 chilling-GDD (NaN: not season_decid, physics skips)
PFTs(17).k_chilla            = NaN;
PFTs(17).k_chillb            = NaN;
PFTs(17).k_chillk            = NaN;
PFTs(17).phengdd5ramp        = NaN;
%   [N.2] Senescence (CLM5 daylength trigger, CNPhenologyMod.F90:957)
PFTs(17).crit_dayl_h         = 10.9;    % [h] dayl threshold for offset trigger (CLM5 global default)
PFTs(17).ndays_on            = 30;      % [days] CLM5 onset ramp duration (CNPhenologyMod.F90:51, clm5_params.c171117.nc scalar)
PFTs(17).ndays_off           = 15;      % [days] CLM5 offset ramp duration (CNPhenologyMod.F90:52, clm5_params.c171117.nc scalar)
PFTs(17).fstor2tran          = 0.5;     % [-] CLM5 budburst fraction (CNPhenologyMod.F90:53, clm5_params.c171117.nc scalar = 0.5)
PFTs(17).GDD_bloom_after_budburst = NaN;     % [°C·d] only Apple PFT 23 has fruit phenology calibrated for Changwu
PFTs(17).TSum_fruit              = NaN;
PFTs(17).harvest_DS              = NaN;
%   [N.3] N retranslocation on tissue turnover (LPJ-GUESS growth.cpp:256-285)
PFTs(17).nrelocfrac          = 0.5;     % [-] N retransloc fraction to N_storage at senescence
%   [STOMATA]
PFTs(17).m                      = 9;
PFTs(17).alpha                  = 0.06;
%   [LEAF]
PFTs(17).Ztop                   = 0.5;
PFTs(17).Zbot                   = 0.01;
PFTs(17).SLA0                   = 0.0350;    % [CLM5 audit 2026-05-26] was 0.07
PFTs(17).SLAm                   = 0;
PFTs(17).CN_L                   = 20.00;    % [CLM5 audit 2026-05-26] was 25
%   [PHOTOSYNTH]
PFTs(17).Vcmax25                = 100.7;
PFTs(17).g_m                    = 0.295;  % [mol m-2 s-1] gm25 (Knauer 2019 GCB25 Table 3 C3C, C3 crop; soybean is C3)
%   [LUNA]
PFTs(17).Jmaxb1                 = 0.0936;
%   [ALLOC_CN]
PFTs(17).CN_target_root         = 42.00;    % [CLM5 NC audit 2026-05-26] was 30
PFTs(17).CN_target_stem         = 50.00;    % [CLM5 NC audit 2026-05-26] was 0
% [13c dual stem pool] CLM5 CN_Allocation Tbl col7-8 (CN_lw/CN_dw); clm5_params (flivewd/lwtop_ann)
PFTs(17).CN_lw                  = 50;    % [gC/gN] live stem (respires)
PFTs(17).CN_dw                  = 500;   % [gC/gN] dead stem (heartwood, ~N-free)
PFTs(17).flivewd                = 1.0;  % [-] new-stem-C frac to live wood (tree/shrub .1/grass 0/crop 1)
PFTs(17).lwtop                  = 0.7;   % [1/yr] live->dead wood turnover (CLM5 lwtop_ann)
PFTs(17).CN_target_fruit        = 15;
PFTs(17).f_alloc_leaf           = 0.45;
PFTs(17).f_alloc_root           = 0.3;
PFTs(17).f_alloc_stem           = 0.15;
PFTs(17).f_alloc_fruit          = 0.1;
%   [N_UPTAKE] (Nacry kinetics)
PFTs(17).Vmax_NH4_umol_g_h      = 20;
PFTs(17).Vmax_NO3_umol_g_h      = 15;
PFTs(17).Km_NH4                 = 0.3;
PFTs(17).Km_NO3                 = 0.3;
PFTs(17).Cmin_NH4               = 0.5;
PFTs(17).Cmin_NO3               = 0.25;
PFTs(17).SRL                    = 50;
PFTs(17).Vmax_NO3_LATS_umol_g_h = 0;
PFTs(17).Km_LATS_NO3            = 70;
%   [LEAF_LIFESPAN]
PFTs(17).leaf_long_days         = 365.00;    % [CLM5 NC audit 2026-05-26] was 120
PFTs(17).retransl_fraction      = 0.4;
%   [RESPIRATION] Plant respiration (Sitch 2003 GCB 9:161 Eq 22-23, p.169)
PFTs(17).r_m10                = 0.066;     % [g C/(g N·d)] @ 10°C base, Sitch 2003 Table 3 "Plant Respiration" –TeH-class –C3 temperate crop
PFTs(17).R_g_fraction         = 0.11;      % [-] growth resp factor (CLM5 grperc; was Sitch 0.25).
PFTs(17).harvest_fraction         = 0.90;
PFTs(17).root_turnover_rate       = 0.00274;
PFTs(17).spring_release_days      = 60;

% ==============================================================================
% PFT 18: Cotton –Cotton
% ==============================================================================
PFTs(18).PFT                    = 'Cotton';
PFTs(18).FullName               = 'Cotton';
PFTs(18).clumping_index         = 0.76;  % [-] Proxy: ClimaLand C3-grass value; no crop-specific entry in its PFT table
%   [HABIT] (PFT lifeform flags)
PFTs(18).photo_type             = 'C3';
PFTs(18).is_perennial           = 0;
PFTs(18).has_fruit              = 1;
PFTs(18).has_storage_pool       = 0;
%   [PHENOLOGY] (see [N.0]-[N.5] sub-stages below)
PFTs(18).evergreen            = 0;
PFTs(18).season_decid         = 0;
PFTs(18).stress_decid         = 0;
PFTs(18).crop_phenology       = 1;   % Cotton
%   [N.1] Sykes 1996 chilling-GDD (NaN: not season_decid, physics skips)
PFTs(18).k_chilla            = NaN;
PFTs(18).k_chillb            = NaN;
PFTs(18).k_chillk            = NaN;
PFTs(18).phengdd5ramp        = NaN;
%   [N.2] Senescence (CLM5 daylength trigger, CNPhenologyMod.F90:957)
PFTs(18).crit_dayl_h         = 10.9;    % [h] dayl threshold for offset trigger (CLM5 global default)
PFTs(18).ndays_on            = 30;      % [days] CLM5 onset ramp duration (CNPhenologyMod.F90:51, clm5_params.c171117.nc scalar)
PFTs(18).ndays_off           = 15;      % [days] CLM5 offset ramp duration (CNPhenologyMod.F90:52, clm5_params.c171117.nc scalar)
PFTs(18).fstor2tran          = 0.5;     % [-] CLM5 budburst fraction (CNPhenologyMod.F90:53, clm5_params.c171117.nc scalar = 0.5)
PFTs(18).GDD_bloom_after_budburst = NaN;     % [°C·d] only Apple PFT 23 has fruit phenology calibrated for Changwu
PFTs(18).TSum_fruit              = NaN;
PFTs(18).harvest_DS              = NaN;
%   [N.3] N retranslocation on tissue turnover (LPJ-GUESS growth.cpp:256-285)
PFTs(18).nrelocfrac          = 0.5;     % [-] N retransloc fraction to N_storage at senescence
%   [STOMATA]
PFTs(18).m                      = 9;
PFTs(18).alpha                  = 0.06;
%   [LEAF]
PFTs(18).Ztop                   = 0.5;
PFTs(18).Zbot                   = 0.01;
PFTs(18).SLA0                   = 0.0350;    % [CLM5 audit 2026-05-26] was 0.03
PFTs(18).SLAm                   = 0;
PFTs(18).CN_L                   = 20.00;    % [CLM5 audit 2026-05-26] was 25
%   [PHOTOSYNTH]
PFTs(18).Vcmax25                = 100.7;
PFTs(18).g_m                    = 0.295;  % [mol m-2 s-1] gm25 (Knauer 2019 GCB25 Table 3 C3C, C3 crop; cotton is C3)
%   [LUNA]
PFTs(18).Jmaxb1                 = 0.0936;
%   [ALLOC_CN]
PFTs(18).CN_target_root         = 42.00;    % [CLM5 NC audit 2026-05-26] was 30
PFTs(18).CN_target_stem         = 50.00;    % [CLM5 NC audit 2026-05-26] was 0
% [13c dual stem pool] CLM5 CN_Allocation Tbl col7-8 (CN_lw/CN_dw); clm5_params (flivewd/lwtop_ann)
PFTs(18).CN_lw                  = 50;    % [gC/gN] live stem (respires)
PFTs(18).CN_dw                  = 500;   % [gC/gN] dead stem (heartwood, ~N-free)
PFTs(18).flivewd                = 1.0;  % [-] new-stem-C frac to live wood (tree/shrub .1/grass 0/crop 1)
PFTs(18).lwtop                  = 0.7;   % [1/yr] live->dead wood turnover (CLM5 lwtop_ann)
PFTs(18).CN_target_fruit        = 30;
PFTs(18).f_alloc_leaf           = 0.45;
PFTs(18).f_alloc_root           = 0.3;
PFTs(18).f_alloc_stem           = 0.15;
PFTs(18).f_alloc_fruit          = 0.1;
%   [N_UPTAKE] (Nacry kinetics)
PFTs(18).Vmax_NH4_umol_g_h      = 20;
PFTs(18).Vmax_NO3_umol_g_h      = 15;
PFTs(18).Km_NH4                 = 0.3;
PFTs(18).Km_NO3                 = 0.3;
PFTs(18).Cmin_NH4               = 0.5;
PFTs(18).Cmin_NO3               = 0.25;
PFTs(18).SRL                    = 50;
PFTs(18).Vmax_NO3_LATS_umol_g_h = 0;
PFTs(18).Km_LATS_NO3            = 70;
%   [LEAF_LIFESPAN]
PFTs(18).leaf_long_days         = 365.00;    % [CLM5 NC audit 2026-05-26] was 150
PFTs(18).retransl_fraction      = 0.4;
%   [RESPIRATION] Plant respiration (Sitch 2003 GCB 9:161 Eq 22-23, p.169)
PFTs(18).r_m10                = 0.066;     % [g C/(g N·d)] @ 10°C base, Sitch 2003 Table 3 "Plant Respiration" –TeH-class –C3 temperate crop
PFTs(18).R_g_fraction         = 0.11;      % [-] growth resp factor (CLM5 grperc; was Sitch 0.25).
PFTs(18).harvest_fraction         = 0.90;
PFTs(18).root_turnover_rate       = 0.00274;
PFTs(18).spring_release_days      = 60;

% ==============================================================================
% PFT 19: Rice –Rice
% ==============================================================================
PFTs(19).PFT                    = 'Rice';
PFTs(19).FullName               = 'Rice';
PFTs(19).clumping_index         = 0.76;  % [-] Proxy: ClimaLand C3-grass value; no crop-specific entry in its PFT table
%   [HABIT] (PFT lifeform flags)
PFTs(19).photo_type             = 'C3';
PFTs(19).is_perennial           = 0;
PFTs(19).has_fruit              = 1;
PFTs(19).has_storage_pool       = 0;
%   [PHENOLOGY] (see [N.0]-[N.5] sub-stages below)
PFTs(19).evergreen            = 0;
PFTs(19).season_decid         = 0;
PFTs(19).stress_decid         = 0;
PFTs(19).crop_phenology       = 1;   % Rice
%   [N.1] Sykes 1996 chilling-GDD (NaN: not season_decid, physics skips)
PFTs(19).k_chilla            = NaN;
PFTs(19).k_chillb            = NaN;
PFTs(19).k_chillk            = NaN;
PFTs(19).phengdd5ramp        = NaN;
%   [N.2] Senescence (CLM5 daylength trigger, CNPhenologyMod.F90:957)
PFTs(19).crit_dayl_h         = 10.9;    % [h] dayl threshold for offset trigger (CLM5 global default)
PFTs(19).ndays_on            = 30;      % [days] CLM5 onset ramp duration (CNPhenologyMod.F90:51, clm5_params.c171117.nc scalar)
PFTs(19).ndays_off           = 15;      % [days] CLM5 offset ramp duration (CNPhenologyMod.F90:52, clm5_params.c171117.nc scalar)
PFTs(19).fstor2tran          = 0.5;     % [-] CLM5 budburst fraction (CNPhenologyMod.F90:53, clm5_params.c171117.nc scalar = 0.5)
PFTs(19).GDD_bloom_after_budburst = NaN;     % [°C·d] only Apple PFT 23 has fruit phenology calibrated for Changwu
PFTs(19).TSum_fruit              = NaN;
PFTs(19).harvest_DS              = NaN;
%   [N.3] N retranslocation on tissue turnover (LPJ-GUESS growth.cpp:256-285)
PFTs(19).nrelocfrac          = 0.5;     % [-] N retransloc fraction to N_storage at senescence
%   [STOMATA]
PFTs(19).m                      = 9;
PFTs(19).alpha                  = 0.06;
%   [LEAF]
PFTs(19).Ztop                   = 0.5;
PFTs(19).Zbot                   = 0.01;
PFTs(19).SLA0                   = 0.0350;    % [CLM5 audit 2026-05-26] was 0.07
PFTs(19).SLAm                   = 0;
PFTs(19).CN_L                   = 20.00;    % [CLM5 audit 2026-05-26] was 25
%   [PHOTOSYNTH]
PFTs(19).Vcmax25                = 100.7;
PFTs(19).g_m                    = 0.295;  % [mol m-2 s-1] gm25 (Knauer 2019 GCB25 Table 3 C3C, C3 crop; rice is C3)
%   [LUNA]
PFTs(19).Jmaxb1                 = 0.0936;
%   [ALLOC_CN]
PFTs(19).CN_target_root         = 42.00;    % [CLM5 NC audit 2026-05-26] was 30
PFTs(19).CN_target_stem         = 50.00;    % [CLM5 NC audit 2026-05-26] was 0
% [13c dual stem pool] CLM5 CN_Allocation Tbl col7-8 (CN_lw/CN_dw); clm5_params (flivewd/lwtop_ann)
PFTs(19).CN_lw                  = 50;    % [gC/gN] live stem (respires)
PFTs(19).CN_dw                  = 500;   % [gC/gN] dead stem (heartwood, ~N-free)
PFTs(19).flivewd                = 1.0;  % [-] new-stem-C frac to live wood (tree/shrub .1/grass 0/crop 1)
PFTs(19).lwtop                  = 0.7;   % [1/yr] live->dead wood turnover (CLM5 lwtop_ann)
PFTs(19).CN_target_fruit        = 30;
PFTs(19).f_alloc_leaf           = 0.45;
PFTs(19).f_alloc_root           = 0.3;
PFTs(19).f_alloc_stem           = 0.15;
PFTs(19).f_alloc_fruit          = 0.1;
%   [N_UPTAKE] (Nacry kinetics)
PFTs(19).Vmax_NH4_umol_g_h      = 20;
PFTs(19).Vmax_NO3_umol_g_h      = 15;
PFTs(19).Km_NH4                 = 0.3;
PFTs(19).Km_NO3                 = 0.3;
PFTs(19).Cmin_NH4               = 0.5;
PFTs(19).Cmin_NO3               = 0.25;
PFTs(19).SRL                    = 50;
PFTs(19).Vmax_NO3_LATS_umol_g_h = 0;
PFTs(19).Km_LATS_NO3            = 70;
%   [LEAF_LIFESPAN]
PFTs(19).leaf_long_days         = 365.00;    % [CLM5 NC audit 2026-05-26] was 120
PFTs(19).retransl_fraction      = 0.4;
%   [RESPIRATION] Plant respiration (Sitch 2003 GCB 9:161 Eq 22-23, p.169)
PFTs(19).r_m10                = 0.066;     % [g C/(g N·d)] @ 10°C base, Sitch 2003 Table 3 "Plant Respiration" –TeH-class –C3 crop
PFTs(19).R_g_fraction         = 0.11;      % [-] growth resp factor (CLM5 grperc; was Sitch 0.25).
PFTs(19).harvest_fraction         = 0.90;
PFTs(19).root_turnover_rate       = 0.00274;
PFTs(19).spring_release_days      = 60;

% ==============================================================================
% PFT 20: Sugarcane –Sugarcane
% ==============================================================================
PFTs(20).PFT                    = 'Sugarcane';
PFTs(20).FullName               = 'Sugarcane';
PFTs(20).clumping_index         = 0.75;  % [-] Proxy: ClimaLand C4-grass value; no crop-specific entry in its PFT table
%   [HABIT] (PFT lifeform flags)
PFTs(20).photo_type             = 'C4';
PFTs(20).is_perennial           = 0;
PFTs(20).has_fruit              = 1;
PFTs(20).has_storage_pool       = 0;
%   [PHENOLOGY] (see [N.0]-[N.5] sub-stages below)
PFTs(20).evergreen            = 0;
PFTs(20).season_decid         = 0;
PFTs(20).stress_decid         = 0;
PFTs(20).crop_phenology       = 1;   % Sugarcane
%   [N.1] Sykes 1996 chilling-GDD (NaN: not season_decid, physics skips)
PFTs(20).k_chilla            = NaN;
PFTs(20).k_chillb            = NaN;
PFTs(20).k_chillk            = NaN;
PFTs(20).phengdd5ramp        = NaN;
%   [N.2] Senescence (CLM5 daylength trigger, CNPhenologyMod.F90:957)
PFTs(20).crit_dayl_h         = 10.9;    % [h] dayl threshold for offset trigger (CLM5 global default)
PFTs(20).ndays_on            = 30;      % [days] CLM5 onset ramp duration (CNPhenologyMod.F90:51, clm5_params.c171117.nc scalar)
PFTs(20).ndays_off           = 15;      % [days] CLM5 offset ramp duration (CNPhenologyMod.F90:52, clm5_params.c171117.nc scalar)
PFTs(20).fstor2tran          = 0.5;     % [-] CLM5 budburst fraction (CNPhenologyMod.F90:53, clm5_params.c171117.nc scalar = 0.5)
PFTs(20).GDD_bloom_after_budburst = NaN;     % [°C·d] only Apple PFT 23 has fruit phenology calibrated for Changwu
PFTs(20).TSum_fruit              = NaN;
PFTs(20).harvest_DS              = NaN;
%   [N.3] N retranslocation on tissue turnover (LPJ-GUESS growth.cpp:256-285)
PFTs(20).nrelocfrac          = 0.5;     % [-] N retransloc fraction to N_storage at senescence
%   [STOMATA]
PFTs(20).m                      = 4;
PFTs(20).alpha                  = 0.05;
%   [LEAF]
PFTs(20).Ztop                   = 0.5;
PFTs(20).Zbot                   = 0.01;
PFTs(20).SLA0                   = 0.05;
PFTs(20).SLAm                   = 0;
PFTs(20).CN_L                   = 25;
%   [PHOTOSYNTH]
PFTs(20).Vcmax25                = 100.7;
PFTs(20).g_m                    = 0.739;  % [mol m-2 s-1] gm25 (Knauer 2019 GCB25 Table 3 C4C, C4 crop; sugarcane is C4)
%   [LUNA]
PFTs(20).Jmaxb1                 = 0.0936;
%   [ALLOC_CN]
PFTs(20).CN_target_root         = 42.00;    % [CLM5 NC audit 2026-05-26] was 30
PFTs(20).CN_target_stem         = 50.00;    % [CLM5 NC audit 2026-05-26] was 0
% [13c dual stem pool] CLM5 CN_Allocation Tbl col7-8 (CN_lw/CN_dw); clm5_params (flivewd/lwtop_ann)
PFTs(20).CN_lw                  = 50;    % [gC/gN] live stem (respires)
PFTs(20).CN_dw                  = 500;   % [gC/gN] dead stem (heartwood, ~N-free)
PFTs(20).flivewd                = 1.0;  % [-] new-stem-C frac to live wood (tree/shrub .1/grass 0/crop 1)
PFTs(20).lwtop                  = 0.7;   % [1/yr] live->dead wood turnover (CLM5 lwtop_ann)
PFTs(20).CN_target_fruit        = 40;
PFTs(20).f_alloc_leaf           = 0.45;
PFTs(20).f_alloc_root           = 0.3;
PFTs(20).f_alloc_stem           = 0.15;
PFTs(20).f_alloc_fruit          = 0.1;
%   [N_UPTAKE] (Nacry kinetics)
PFTs(20).Vmax_NH4_umol_g_h      = 20;
PFTs(20).Vmax_NO3_umol_g_h      = 15;
PFTs(20).Km_NH4                 = 0.3;
PFTs(20).Km_NO3                 = 0.3;
PFTs(20).Cmin_NH4               = 0.5;
PFTs(20).Cmin_NO3               = 0.25;
PFTs(20).SRL                    = 50;
PFTs(20).Vmax_NO3_LATS_umol_g_h = 0;
PFTs(20).Km_LATS_NO3            = 70;
%   [LEAF_LIFESPAN]
PFTs(20).leaf_long_days         = 365.00;    % [CLM5 NC audit 2026-05-26] was 300
PFTs(20).retransl_fraction      = 0.4;
%   [RESPIRATION] Plant respiration (Sitch 2003 GCB 9:161 Eq 22-23, p.169)
PFTs(20).r_m10                = 0.011;     % [g C/(g N·d)] @ 10°C base, Sitch 2003 Table 3 "Plant Respiration" –TrH-class –C4 tropical crop (Sitch TrH=0.011)
PFTs(20).R_g_fraction         = 0.11;      % [-] growth resp factor (CLM5 grperc; was Sitch 0.25).
PFTs(20).harvest_fraction         = 0.90;
PFTs(20).root_turnover_rate       = 0.00274;
PFTs(20).spring_release_days      = 60;

% ==============================================================================
% PFT 21: Tropical Corn –Tropical Corn
% ==============================================================================
PFTs(21).PFT                    = 'Tropical Corn';
PFTs(21).FullName               = 'Tropical Corn';
PFTs(21).clumping_index         = 0.75;  % [-] Proxy: ClimaLand C4-grass value; no crop-specific entry in its PFT table
%   [HABIT] (PFT lifeform flags)
PFTs(21).photo_type             = 'C4';
PFTs(21).is_perennial           = 0;
PFTs(21).has_fruit              = 1;
PFTs(21).has_storage_pool       = 0;
%   [PHENOLOGY] (see [N.0]-[N.5] sub-stages below)
PFTs(21).evergreen            = 0;
PFTs(21).season_decid         = 0;
PFTs(21).stress_decid         = 0;
PFTs(21).crop_phenology       = 1;   % Trop Corn
%   [N.1] Sykes 1996 chilling-GDD (NaN: not season_decid, physics skips)
PFTs(21).k_chilla            = NaN;
PFTs(21).k_chillb            = NaN;
PFTs(21).k_chillk            = NaN;
PFTs(21).phengdd5ramp        = NaN;
%   [N.2] Senescence (CLM5 daylength trigger, CNPhenologyMod.F90:957)
PFTs(21).crit_dayl_h         = 10.9;    % [h] dayl threshold for offset trigger (CLM5 global default)
PFTs(21).ndays_on            = 30;      % [days] CLM5 onset ramp duration (CNPhenologyMod.F90:51, clm5_params.c171117.nc scalar)
PFTs(21).ndays_off           = 15;      % [days] CLM5 offset ramp duration (CNPhenologyMod.F90:52, clm5_params.c171117.nc scalar)
PFTs(21).fstor2tran          = 0.5;     % [-] CLM5 budburst fraction (CNPhenologyMod.F90:53, clm5_params.c171117.nc scalar = 0.5)
PFTs(21).GDD_bloom_after_budburst = NaN;     % [°C·d] only Apple PFT 23 has fruit phenology calibrated for Changwu
PFTs(21).TSum_fruit              = NaN;
PFTs(21).harvest_DS              = NaN;
%   [N.3] N retranslocation on tissue turnover (LPJ-GUESS growth.cpp:256-285)
PFTs(21).nrelocfrac          = 0.5;     % [-] N retransloc fraction to N_storage at senescence
%   [STOMATA]
PFTs(21).m                      = 4;
PFTs(21).alpha                  = 0.05;
%   [LEAF]
PFTs(21).Ztop                   = 0.5;
PFTs(21).Zbot                   = 0.01;
PFTs(21).SLA0                   = 0.05;
PFTs(21).SLAm                   = 0;
PFTs(21).CN_L                   = 25;
%   [PHOTOSYNTH]
PFTs(21).Vcmax25                = 100.7;
PFTs(21).g_m                    = 0.739;  % [mol m-2 s-1] gm25 (Knauer 2019 GCB25 Table 3 C4C, C4 crop; tropical corn/maize is C4)
%   [LUNA]
PFTs(21).Jmaxb1                 = 0.0936;
%   [ALLOC_CN]
PFTs(21).CN_target_root         = 42.00;    % [CLM5 NC audit 2026-05-26] was 30
PFTs(21).CN_target_stem         = 50.00;    % [CLM5 NC audit 2026-05-26] was 0
% [13c dual stem pool] CLM5 CN_Allocation Tbl col7-8 (CN_lw/CN_dw); clm5_params (flivewd/lwtop_ann)
PFTs(21).CN_lw                  = 50;    % [gC/gN] live stem (respires)
PFTs(21).CN_dw                  = 500;   % [gC/gN] dead stem (heartwood, ~N-free)
PFTs(21).flivewd                = 1.0;  % [-] new-stem-C frac to live wood (tree/shrub .1/grass 0/crop 1)
PFTs(21).lwtop                  = 0.7;   % [1/yr] live->dead wood turnover (CLM5 lwtop_ann)
PFTs(21).CN_target_fruit        = 15;
PFTs(21).f_alloc_leaf           = 0.45;
PFTs(21).f_alloc_root           = 0.3;
PFTs(21).f_alloc_stem           = 0.15;
PFTs(21).f_alloc_fruit          = 0.1;
%   [N_UPTAKE] (Nacry kinetics)
PFTs(21).Vmax_NH4_umol_g_h      = 20;
PFTs(21).Vmax_NO3_umol_g_h      = 15;
PFTs(21).Km_NH4                 = 0.3;
PFTs(21).Km_NO3                 = 0.3;
PFTs(21).Cmin_NH4               = 0.5;
PFTs(21).Cmin_NO3               = 0.25;
PFTs(21).SRL                    = 50;
PFTs(21).Vmax_NO3_LATS_umol_g_h = 0;
PFTs(21).Km_LATS_NO3            = 70;
%   [LEAF_LIFESPAN]
PFTs(21).leaf_long_days         = 365.00;    % [CLM5 NC audit 2026-05-26] was 90
PFTs(21).retransl_fraction      = 0.4;
%   [RESPIRATION] Plant respiration (Sitch 2003 GCB 9:161 Eq 22-23, p.169)
PFTs(21).r_m10                = 0.011;     % [g C/(g N·d)] @ 10°C base, Sitch 2003 Table 3 "Plant Respiration" –TrH-class –C4 tropical crop
PFTs(21).R_g_fraction         = 0.11;      % [-] growth resp factor (CLM5 grperc; was Sitch 0.25).
PFTs(21).harvest_fraction         = 0.90;
PFTs(21).root_turnover_rate       = 0.00274;
PFTs(21).spring_release_days      = 60;

% ==============================================================================
% PFT 22: Tropical Soybean –Tropical Soybean
% ==============================================================================
PFTs(22).PFT                    = 'Tropical Soybean';
PFTs(22).FullName               = 'Tropical Soybean';
PFTs(22).clumping_index         = 0.76;  % [-] Proxy: ClimaLand C3-grass value; no crop-specific entry in its PFT table
%   [HABIT] (PFT lifeform flags)
PFTs(22).photo_type             = 'C3';
PFTs(22).is_perennial           = 0;
PFTs(22).has_fruit              = 1;
PFTs(22).has_storage_pool       = 0;
%   [PHENOLOGY] (see [N.0]-[N.5] sub-stages below)
PFTs(22).evergreen            = 0;
PFTs(22).season_decid         = 0;
PFTs(22).stress_decid         = 0;
PFTs(22).crop_phenology       = 1;   % Trop Soybean
%   [N.1] Sykes 1996 chilling-GDD (NaN: not season_decid, physics skips)
PFTs(22).k_chilla            = NaN;
PFTs(22).k_chillb            = NaN;
PFTs(22).k_chillk            = NaN;
PFTs(22).phengdd5ramp        = NaN;
%   [N.2] Senescence (CLM5 daylength trigger, CNPhenologyMod.F90:957)
PFTs(22).crit_dayl_h         = 10.9;    % [h] dayl threshold for offset trigger (CLM5 global default)
PFTs(22).ndays_on            = 30;      % [days] CLM5 onset ramp duration (CNPhenologyMod.F90:51, clm5_params.c171117.nc scalar)
PFTs(22).ndays_off           = 15;      % [days] CLM5 offset ramp duration (CNPhenologyMod.F90:52, clm5_params.c171117.nc scalar)
PFTs(22).fstor2tran          = 0.5;     % [-] CLM5 budburst fraction (CNPhenologyMod.F90:53, clm5_params.c171117.nc scalar = 0.5)
PFTs(22).GDD_bloom_after_budburst = NaN;     % [°C·d] only Apple PFT 23 has fruit phenology calibrated for Changwu
PFTs(22).TSum_fruit              = NaN;
PFTs(22).harvest_DS              = NaN;
%   [N.3] N retranslocation on tissue turnover (LPJ-GUESS growth.cpp:256-285)
PFTs(22).nrelocfrac          = 0.5;     % [-] N retransloc fraction to N_storage at senescence
%   [STOMATA]
PFTs(22).m                      = 9;
PFTs(22).alpha                  = 0.06;
%   [LEAF]
PFTs(22).Ztop                   = 0.5;
PFTs(22).Zbot                   = 0.01;
PFTs(22).SLA0                   = 0.0350;    % [CLM5 audit 2026-05-26] was 0.07
PFTs(22).SLAm                   = 0;
PFTs(22).CN_L                   = 20.00;    % [CLM5 audit 2026-05-26] was 25
%   [PHOTOSYNTH]
PFTs(22).Vcmax25                = 100.7;
PFTs(22).g_m                    = 0.295;  % [mol m-2 s-1] gm25 (Knauer 2019 GCB25 Table 3 C3C, C3 crop; soybean is C3)
%   [LUNA]
PFTs(22).Jmaxb1                 = 0.0936;
%   [ALLOC_CN]
PFTs(22).CN_target_root         = 42.00;    % [CLM5 NC audit 2026-05-26] was 30
PFTs(22).CN_target_stem         = 50.00;    % [CLM5 NC audit 2026-05-26] was 0
% [13c dual stem pool] CLM5 CN_Allocation Tbl col7-8 (CN_lw/CN_dw); clm5_params (flivewd/lwtop_ann)
PFTs(22).CN_lw                  = 50;    % [gC/gN] live stem (respires)
PFTs(22).CN_dw                  = 500;   % [gC/gN] dead stem (heartwood, ~N-free)
PFTs(22).flivewd                = 1.0;  % [-] new-stem-C frac to live wood (tree/shrub .1/grass 0/crop 1)
PFTs(22).lwtop                  = 0.7;   % [1/yr] live->dead wood turnover (CLM5 lwtop_ann)
PFTs(22).CN_target_fruit        = 15;
PFTs(22).f_alloc_leaf           = 0.45;
PFTs(22).f_alloc_root           = 0.3;
PFTs(22).f_alloc_stem           = 0.15;
PFTs(22).f_alloc_fruit          = 0.1;
%   [N_UPTAKE] (Nacry kinetics)
PFTs(22).Vmax_NH4_umol_g_h      = 20;
PFTs(22).Vmax_NO3_umol_g_h      = 15;
PFTs(22).Km_NH4                 = 0.3;
PFTs(22).Km_NO3                 = 0.3;
PFTs(22).Cmin_NH4               = 0.5;
PFTs(22).Cmin_NO3               = 0.25;
PFTs(22).SRL                    = 50;
PFTs(22).Vmax_NO3_LATS_umol_g_h = 0;
PFTs(22).Km_LATS_NO3            = 70;
%   [LEAF_LIFESPAN]
PFTs(22).leaf_long_days         = 365.00;    % [CLM5 NC audit 2026-05-26] was 120
PFTs(22).retransl_fraction      = 0.4;
%   [RESPIRATION] Plant respiration (Sitch 2003 GCB 9:161 Eq 22-23, p.169)
PFTs(22).r_m10                = 0.066;     % [g C/(g N·d)] @ 10°C base, Sitch 2003 Table 3 "Plant Respiration" –TrH-C3-class –tropical C3 crop, Sitch makes no tropical-C3 herb dist., use temperate default
PFTs(22).R_g_fraction         = 0.11;      % [-] growth resp factor (CLM5 grperc; was Sitch 0.25).
PFTs(22).harvest_fraction         = 0.90;
PFTs(22).root_turnover_rate       = 0.00274;
PFTs(22).spring_release_days      = 60;

% ==============================================================================
% PFT 23: Apple –Apple Tree (Malus domestica)
% ==============================================================================
PFTs(23).PFT                    = 'Apple';
PFTs(23).FullName               = 'Apple Tree (Malus domestica)';
PFTs(23).clumping_index         = 0.80;  % [-] Apple orchard measurement: Liu et al. (2013), Sect. 3.4.4 and Fig. 10, doi:10.1016/j.scienta.2013.10.009
%   [HABIT] (PFT lifeform flags)
PFTs(23).photo_type             = 'C3';
PFTs(23).is_perennial           = 1;
PFTs(23).has_fruit              = 1;
PFTs(23).has_storage_pool       = 1;
% [dual-storage 2026-06-06] INITIAL split of the overwinter NSC reserve between coarse-root
% and stem(sapwood) pools. ONLY used to partition the spin-up initial C/N_storage at KT=1
% (no history then). At RUNTIME the two reserves evolve INDEPENDENTLY by source organ
% (root growth's (1-fcur) -> root reserve; stem/leaf/fruit -> stem reserve; spring growth
% drains each), so the ratio drifts freely and is NOT pinned to this value.
% [regrounded 2026-06-10] Reserve lives in perennial woody tissue; total pool = NSC% x woody biomass.
% Aboveground wood biomass (Zanotelli et al. 2013 Biogeosciences 10:3089 Table 2: woodAG 0.84 kgC/m2)
% far exceeds coarse root (woodBG 0.19 kgC/m2), while starch concentration is COMPARABLE in root vs
% stem/branch (apple root starch 15-18% DW, Stutte et al. 1994 JASHS 119:492 Tables 4-5; peach winter
% root 11.6-14.5% DW, Lopez et al. 2007 Tree Physiol 27:1619 Fig.3; branch/shoot starch in the same
% range, Stutte 1994 p.495). So total reserve is STEM-dominant: f_storage_root = woodBG/(woodAG+woodBG)
% = 0.19/1.03 = 0.18. This REVERSES the prior root-dominant 0.55 (which wrongly read root-starch
% concentration as total pool). Spring demand agrees: stem reserve funds the large leaf flush, root
% reserve the smaller fine-root flush.
PFTs(23).f_storage_root         = 0.18;   % [-] INITIAL fraction of reserve in coarse root (stem-dominant)

% [target-reserve 2026-06-10] nsc_frac: woody NSC CARBON fraction = winter-MAX NSC%(dry weight) x
% 0.44/0.47 (starch C-frac / wood C-frac = x0.936). Sets each woody reserve pool's TARGET = nsc_frac x
% its own wood C (target-driven refill, replacing the drift-prone fixed-fraction fcur loading; FATES
% target idea Massoud et al. 2019 GMD 12:4133, made wood-based per-pool). The 14% is the WINTER MAXIMUM
% ("full" reserve); actual NSC cycles below it seasonally (spring draw -> autumn refill), emergent from
% the target dynamics, not imposed. Apple 0.131 = 14%DW grounded (Lopez et al. 2007 Tree Physiol
% 27:1619 Fig.3 peach root 11.6-14.5%; Stutte et al. 1994 JASHS 119:492 Tables 4-5 apple root 15-18%).
% Wild forest woody NSC is much LOWER: decid broadleaf tree ~8%DW, decid needleleaf/shrub ~7%DW
% (Martinez-Vilalta et al. 2016 Ecol Monogr 86:495 global 121-study DB: stem 70 mg/g=7%, belowground
% 106 mg/g=10.6%, conifers lowest; Furze et al. 2018 New Phytol 221:1466 forest woody NSC ~3.6% ->
% models overpredict, so the managed-orchard 14% sits at the TOP of the range). Used only where
% has_storage_pool=1 (the 8 deciduous/stress-decid woody PFTs); evergreen/grass/crop default 0.07 unused.
% Consolidated here (per-pool reserve-target params); one value per PFT (apple root~stem, MV2016 root>stem
% but a single woody value is a defensible 1st-order choice).
PFTs(23).nsc_frac = 0.131;   % Apple — managed fruit tree, 14%DW (grounded)
PFTs(7).nsc_frac  = 0.075;   % BDT temperate — decid broadleaf tree, 8%DW (MV2016)
PFTs(8).nsc_frac  = 0.075;   % BDT boreal    — decid broadleaf tree, 8%DW (MV2016)
PFTs(6).nsc_frac  = 0.075;   % BDT tropical  — stress-decid broadleaf, 8%DW (MV2016)
PFTs(24).nsc_frac = 0.075;   % Black Locust (Robinia) — decid broadleaf tree, 8%DW (MV2016)
PFTs(3).nsc_frac  = 0.066;   % NDT Boreal (larch) — decid needleleaf, 7%DW (MV2016; larch-specific not found)
PFTs(10).nsc_frac = 0.066;   % BDS temperate — decid broadleaf shrub, 7%DW (MV2016)
PFTs(11).nsc_frac = 0.066;   % BDS boreal    — decid broadleaf shrub, 7%DW (MV2016)
% ============================================================================================
%   [PHENOLOGY] PHENOLOGY — ordered by phenological calendar: dormancy-release -> budburst -> spring
%       remob -> bloom -> fruit fill -> harvest -> leaf-fall -> dormancy. Carbon allocation
%       phases that this calendar drives are in [N.6]. Reorganized 2026-06-09 (was scattered).
% ============================================================================================
%   [N.0] Habit classification (exactly one flag = 1; dispatched in phenologyDispatch.m)
PFTs(23).evergreen            = 0;     % [0/1] evergreen leaf habit
PFTs(23).season_decid         = 1;     % [0/1] seasonal-deciduous (apple) <-- ACTIVE
PFTs(23).stress_decid         = 0;     % [0/1] drought-deciduous
PFTs(23).crop_phenology       = 0;     % [0/1] annual-crop sow/harvest cycle
%
%   [N.1] DORMANCY RELEASE & BUDBURST — two chilling models; apple uses Cesaraccio (Model B).
%   --- Model A: Sykes 1996 chilling-GDD (LPJ-GUESS growth.cpp:91-97) — fallback, NOT used by apple ---
PFTs(23).k_chilla            = 0;       % [°C·d] Sykes 1996 a, LPJ-GUESS default
PFTs(23).k_chillb            = 1150;    % [°C·d] Sykes 1996 b, Apple –Sykes 1996 original Fag_syl value (LPJ-GUESS europe.ins:581 notes halved to 600 for European LPJ-GUESS calibration; apple needs unhalved due to high chilling 600-1200 chill hours requirement)
PFTs(23).k_chillk            = 0.05;    % [d^-1] Sykes 1996 k, LPJ-GUESS default
%   --- Model B: Cesaraccio et al. 2004 chilling (apple-specific, nc idx35) — ACTIVE for apple ---
%   Budburst needs two gates: (1) accumulate chill_day (cold days, T below baset) down to crequ;
%   (2) then accumulate anti_chill_day (warmth) until it offsets the chill deficit. Blocks
%   warm-winter false budburst. Apple chilling requirement is high (600-1200 chill hours), so it
%   cannot use the halved European Sykes value.
PFTs(23).use_cesaraccio_chill = 1;      % [0/1] 1=apple uses Cesaraccio chill model (other PFTs Sykes)
PFTs(23).crequ               = -68;     % crequ = "chilling requirement" [°C·d] chill_day accumulates down to this (Cesaraccio 2004 symbol)
PFTs(23).baset               = 4;       % baset = "base temperature" [°C] Cesaraccio chilling base temp (cold days below this count)
PFTs(23).phengdd5ramp        = 200;     % phengdd5ramp = "phenology GDD(base 5°C) ramp" [°C·d] post-budburst GDD over which phen ramps 0->1 (full leaf-on). LPJ-GUESS symbol.
PFTs(23).ndays_on            = 30;      % ndays_on = "number of days, onset" [days] budburst leaf-display ramp duration (CLM5 CNPhenologyMod.F90:51)
%
%   [N.2] SPRING NSC REMOBILISATION — at budburst, fraction of each organ reserve sent to xfer.
% [field-validated 2026-06-09] Apple 0.5 confirmed: Cheng et al. 2021 Front. Plant Sci. 12:582661
%   Fig.7 + Sect.3, root N declined 43.2-49.7% of total tree N over spring budburst (roots = main
%   remobilisation source), bracketing 0.5. = CLM5 CNPhenologyMod.F90:53 scalar 0.5. The root-
%   reserve decline was a refill shortfall (fixed via croot_stem in [ALLOC_CN]), NOT an over-large drawdown.
PFTs(23).fstor2tran          = 0.5;     % fstor2tran = "fraction storage to transfer" [-] spring reserve->xfer fraction at budburst (apple, Cheng 2021; = CLM5 0.5)
PFTs(23).ndays_stor          = 50;      % ndays_stor = "number of days, storage growth" [days] perennial storage-growth onset duration (CLM5 FruitTree citrus idx35; base trees use ndays_on)
%
%   [N.3] BLOOM -> FRUIT FILL -> HARVEST (GDD calendar; two GDD origins: 'budburst' or 'bloom')
%   Stages: budburst --GDD_bloom--> bloom --grnfill--> fruit-fill --(DS=GDD_bloom/TSum_fruit)-->
%           DS>=harvest_DS -> harvest (removes harvest_fraction of fruit, see [L]).
PFTs(23).GDD_bloom_after_budburst = 369;     % [°C·d from budburst] budburst(DOY~80)->bloom(DOY~130). 25-yr Changwu Tmax/Tmin (2000-2024), mean 368.7
PFTs(23).grnfill             = 400;     % grnfill = "grain-fill start" [°C·d from budburst] budburst -> fruit-fill start (carbon begins routing to fruit; fruitTreeAllocation phase switch). CLM5 AgroIBIS symbol.
PFTs(23).grnrp               = 1100;    % grnrp = "ripening start (huiripe offset)" [°C·d from budburst] P2 fruit-fill -> P3 ripening boundary; root reaches arootf2=0.08 here (FruitTree.nc citrus; CNPhenologyMod.F90:1952). [5-phase 2026-06-13]
PFTs(23).hybgdd              = 2880;    % hybgdd = "hybrid GDD" (cultivar total-season GDD) [°C·d from budburst] budburst -> physiological maturity (fruitTreeAllocation fruit-phase end). CLM5 AgroIBIS symbol.
PFTs(23).TSum_fruit          = 2208;    % TSum_fruit = "temperature sum, fruit" [°C·d from bloom] bloom(DOY~130)->commercial harvest(DOY~285); drives DS_fruit = GDD_since_bloom/TSum_fruit. 25-yr Changwu mean 2208.3
PFTs(23).harvest_DS          = 0.95;    % harvest_DS = "harvest development stage" [-] fruit DS that triggers harvest (DAISY/EPIC maturity; Williams 1989 EPIC harvest index). harvest_fraction is in [TURNOVER_HARVEST].
%
%   [N.4] LEAF SENESCENCE / FALL (CLM5 daylength offset trigger, CNPhenologyMod.F90:957)
PFTs(23).crit_dayl_h         = 10.9;    % [h] daylength below which leaf-fall (offset) triggers (CLM5 global default; apple obs onset ~DOY250, this gives DOY~294 — known late bias)
PFTs(23).ndays_off           = 45;      % [days] leaf-fall ramp. CLM5 scalar default = 15 (CNPhenologyMod.F90:52),
                                         % but Changwu apple leaf-fall spans ~75 d (DOY250-325, LAI xlsx 247obs).
                                         % Capped at 45 here: senescence onset stays DOY~293 (crit_dayl_h), so 45 d
                                         % falls out by ~DOY338 without crossing year-end; full 75 d needs earlier
                                         % onset (crit_dayl_h) which is deferred. [site-calibrated 2026-06-07]
PFTs(23).nrelocfrac          = 0.7;     % nrelocfrac = "N relocation fraction" [-] N retranslocated to N_storage at leaf senescence (the "leaf yellowing" recovery; LPJ-GUESS growth.cpp:256-285)
PFTs(23).prune_fr            = 0.85;    % prune_fr = "prune fraction" [-] fraction of annual deadstem (wood) increment pruned to litter each autumn (CLM5 FruitTree citrus idx35; mulch_pruning=1 -> litter)
%
%   [N.5] FRUITTREE 4-PHASE GDD CARBON ALLOCATION (drives leaf taper + organ shares vs GDD).
%   Leaf alloc fraction = f_alloc_leaf x leafGDDtaper(g), g = GDD since budburst (base 5). Taper
%   declines 1->0 over [0, ft_GDDlfmat] so C_leaf plateaus at canopy maturity. Phase root shares
%   step arooti->arootf->arootf2 across the season (see fruitTreeAllocation.m). Source: CLM5-
%   FruitTree (Dombrowski et al. 2022 GMD 15:5167 Eq.4 + Table C1; NutrientCompetitionFlexibleCN
%   Mod.F90:1663-1665; clm5_params.c171117__FruitTree.nc idx35 pftname 'citrus' = apple).
PFTs(23).use_fruittree_alloc = 1;       % [0/1] route apple leaf alloc through leafGDDtaper (other PFTs = 0)
PFTs(23).ft_GDDlfmat         = 1350;    % ft_GDDlfmat = "FruitTree GDD to leaf maturity" [°C·d] budburst -> canopy maturity; leaf-alloc taper -> 0 here (nc idx35 'lfmat')
PFTs(23).ft_bfact            = -0.5;    % ft_bfact = "FruitTree b-factor" [-] J-Norman exponential leaf-taper shape factor (nc idx35 'bfact')
PFTs(23).fleafi              = 0.85;    % fleafi = "fraction leaf, initial" [-] initial leaf fraction of (1-aroot) at budburst (CLM5 AgroIBIS symbol)
PFTs(23).arooti              = 0.7;     % arooti = "allocation root, initial" [-] budburst root allocation fraction, phase-1 start (CLM5 AgroIBIS symbol)
PFTs(23).arootf              = 0.2;     % arootf = "allocation root, final" [-] phase-1-end root fraction (CLM5 AgroIBIS symbol)
PFTs(23).arootf2             = 0.08;    % arootf2 = "allocation root, final-2" [-] fruit-phase root fraction (CLM5 AgroIBIS symbol)
PFTs(23).astemf              = 0.22;    % astemf = "allocation stem, final" [-] mature stem fraction floor (CLM5 AgroIBIS symbol)
PFTs(23).aleafstor           = 0.3;     % aleafstor = "allocation leaf to storage" [-] post-harvest leaf->storage refill fraction
%   [STOMATA]
PFTs(23).m                      = 9;
PFTs(23).alpha                  = 0.06;
%   [LEAF]
PFTs(23).Ztop                   = 0;
PFTs(23).Zbot                   = 0;
% SLA0 recalibration 2026-06-02: 0.0125 (Wang 2009 site) -> 0.028 (Dabbaghi et al.
%   2022 GMD 15:5167 Table C1 slatop, Zanotelli 2013 growing-season LAI+leaf biomass).
%   Reason: LUNA lnc = 1/(SLA0 x CN_L) (lnc_opt=false, Climate10dBuffer.m:234). With
%   SLA0=0.0125, CN_L=19.7 -> lnc=4.06 gN/m2 leaf, ~2x Greer 2018 measured 2.0-2.3,
%   inflating Vcmax/Rd -> CUE 0.43 (obs 0.71). SLA0=0.028 -> lnc=1.81, in measured range.
PFTs(23).SLA0                   = 0.028;     % [m²/g C] Apple top-canopy SLA (Dabbaghi 2022; was 0.0125 Wang 2009)
PFTs(23).SLAm                   = 0.00133;   % [m²/g C / LAI, site] dSLA/dLAI canopy gradient (CLM5 BDT_temperate; CLM4.5 Tech Note Table 8.1) –
% Dabbaghi et al. 2022 GMD 15:5167 Table C1 leafCN=19.7 (Zanotelli 2010, 6 measurements, apple)
PFTs(23).CN_L                   = 19.7;   % [gC/gN] leaf C:N
                                          % Source: CLM5-FruitTree (Dabbaghi 2022 GMD 15:5167) Table C1 leafCN=19.7, Zanotelli 2010
                                          % (data/references/allocation_allometry/Dabbaghi_2022_CLM5_FruitTree_GMD.pdf p.20).
%   [PHOTOSYNTH]
PFTs(23).Vcmax25                = 105;
PFTs(23).g_m                    = 0.175;  % [mol m-2 s-1] gm25 (Knauer 2019 GCB25 Table 3 DBF, deciduous broadleaf; apple is a deciduous broadleaf tree)
%   [LUNA]
% Jmaxb1 reverted to 0.244 (2026-06-02): the 0.18 trial barely moved Vcmax (95.6->89.8,
%   LUNA compensates by re-allocating N to carboxylation), so tuning Jmaxb1 to lower Vcmax
%   does not work. Instead, Vcmax is lowered upstream via SLA0 0.0125->0.028 (lnc 4.06->1.81,
%   see SLA0 note above) — LUNA Vcmax scales with lnc. Jmaxb1 kept at its Wang-2009-based
%   reverse-fit value (Vcmax25=105 target).
PFTs(23).Jmaxb1                 = 0.244;      % [-] Apple-orchard site-calibrated (reverse-fit Vcmax25=105, Wang 2009)
%   [ALLOC_CN]
% Dabbaghi et al. 2022 GMD 15:5167 Table C1 frootCN=32 (Zanotelli 2010, 6 measurements, apple)
PFTs(23).CN_target_root         = 32;    % [gC/gN] fine-root C:N
PFTs(23).CN_target_stem         = 60;    % [gC/gN] (single-pool, kept for compat; dual-pool uses CN_lw/CN_dw below)
% [13c dual stem pool, recalibrated 2026-06-02]
% Live wood: Dabbaghi et al. 2022 GMD 15:5167 Table C1 livewdCN=60 (Zanotelli 2010, 6 measurements, apple).
%   (was 50 from CLM5 generic livestem; Dabbaghi apple-measured value is 60.)
% Dead wood: no apple measurement; heartwood ~N-free, CLM5 deadwdcn high-C:N convention.
PFTs(23).CN_lw                  = 60;    % [gC/gN] live stem (phloem+ray parenchyma, respires)
PFTs(23).CN_dw                  = 500;   % [gC/gN] dead stem (heartwood, ~N-free, no respiration)
% Dabbaghi et al. 2022 GMD 15:5167 Table C1 flivewd=0.15 (Same as BDT in CLM5)
%   (was 0.22; Dabbaghi apple value 0.15 = CLM5 broadleaf-deciduous-tree default.)
% lwtop: CLM5 clm5_params.c171117.nc lwtop_ann (no apple-specific value; CLM5 default kept).
PFTs(23).flivewd                = 0.15;  % [-] new-stem-C fraction to live wood (apple, Dabbaghi 2022)
% croot_stem: coarse-root:stem (belowground-woody : aboveground-woody) allocation ratio.
%   Routes the coarse-root share of the stem allocation stream into the belowground (root)
%   reserve at storage time, matching CLM5 cpool_to_livecrootc_storage = nlc*croot_stem*f3*f4*
%   (1-fcur) (NutrientCompetitionFlexibleCNMod.F90:530,532). Consumed in computeCarbonAllocation.
% [site-calibrated 2026-06-09] Apple value 0.225 from measured Changwu-comparable orchard
%   allocation: Zanotelli et al. 2013 Biogeosciences 10:3089 Table 2 + Sect.3.3, woodBG (coarse
%   root + belowground trunk) = 16.3% vs woodAG (aboveground woody) = 72.3% of NPP ->
%   woodBG/woodAG = 16.3/72.3 = 0.225. Replaces the CLM5 citrus proxy 0.08 (clm5_params
%   .c171117__FruitTree.nc idx35), which is a generic evergreen-citrus value ~2.8x too small for
%   apple: with 0.08 the belowground reserve refilled only from the thin fine-root flux and
%   decayed 112->61->38 gC over 2 yr (root N spring drawdown 43-50%, Cheng et al. 2021 Front.
%   Plant Sci. 12:582661 Fig.7, far exceeds the refill). Natural tree/shrub PFTs keep 0.30.
PFTs(23).croot_stem             = 0.225; % [-] coarse-root:stem STANDING ratio (apple, Zanotelli 2013 woodBG/woodAG=16.3/72.3); reserve targeting
% [source-verified 2026-06-13] Values read DIRECTLY from the CLM5-FruitTree fork parameter file
%   clm5_params.c171117__FruitTree.nc, citrus PFT (idx 35) -- the apple/fruit-tree PFT Dombrowski
%   et al. 2022 GMD 15:5167 use. Allocation flux f1/f2 (CNAllocation): cpool_to_frootc = nlc*f1,
%   cpool_to_(live)crootc = nlc*f2*f3*f4 (NutrientCompetitionCLM45defaultMod.F90:351,358).
PFTs(23).croot_stem_alloc       = 0.08;  % [-] f2 coarse-root:stem (FruitTree.nc citrus; was 0.83 = misread Zanotelli, 10x too high)
PFTs(23).froot_leaf             = 2.0;   % [-] f1 fine-root:leaf (FruitTree.nc citrus; was 0.09, 22x too low -> fine root starved)
PFTs(23).lwtop                  = 0.7;   % [1/yr] live->dead wood turnover (FruitTree.nc lwtop_ann=0.7)
PFTs(23).r_mort                 = 0.02;  % [1/yr] background gap mortality (FruitTree.nc r_mort=0.02; CNGapMortalityMod.F90:196-200 removes deadstem/deadcroot -> litter/CWD; COMPASS was MISSING this)
PFTs(23).CN_target_fruit        = 139;   % [gC/gN]
                                         % Source: CLM5-FruitTree (Dabbaghi/Dombrowski et al. 2022 GMD 15:5167) Table C1 grainCN=139,
                                         % Zanotelli 2010
% NPP allocation fractions calibrated to mature apple orchard observations.
% [NOTE 2026-06-09] Apple has use_fruittree_alloc=1, so it goes through fruitTreeAllocation.m
%   (GDD-phase a_ft.leaf/root/stem/repr), NOT this static f_alloc_* set. These f_alloc_* are
%   the fallback for use_fruittree_alloc~=1 (other PFTs / legacy). In fruitTreeAllocation,
%   a_ft.root = fine root and a_ft.stem = aboveground stem only (coarse root is added separately
%   to the belowground reserve via croot_stem in computeCarbonAllocation, see [ALLOC_CN]).
% Zanotelli et al. 2013 Biogeosciences 10:3089 (Fuji/M9, total NPP 960 gC/m2/yr):
%   fruit 52% (p.3098), leaf 11% (p.3097), woody AG+BG ~19% (p.3097), fine root standing ~1%.
% Cross-check Wu et al. 2012 PLoS ONE 7(6):e38883 (Beijing 18-yr): fruit 46/wood 37/leaf 14/root 3%.
% Prior values (leaf 0.40/fruit 0.15) inverted leaf<->fruit -> C_leaf 2x high + yield 3x low.
% fruit fraction is the in-window instantaneous share (DOY gate in computeOrganDemand
% thins the annual mean); see Zanotelli abstract "largely allocated to fruit: 49%".
PFTs(23).f_alloc_leaf           = 0.12;   % Zanotelli 2013 leaf 11% (fallback only; apple uses fruitTreeAllocation)
PFTs(23).f_alloc_root           = 0.08;   % Zanotelli 2013 fine-root standing ~1% + margin (fallback only)
PFTs(23).f_alloc_stem           = 0.30;   % Zanotelli 2013 woody AG+BG ~19% (fallback only; in fruitTreeAllocation a_ft.stem is AG stem only)
PFTs(23).f_alloc_fruit          = 0.50;   % Zanotelli 2013 fruit 52% in-window (fallback only)
%   [N_UPTAKE] (Nacry kinetics)
PFTs(23).Vmax_NH4_umol_g_h      = 12;
PFTs(23).Vmax_NO3_umol_g_h      = 5;
PFTs(23).Km_NH4                 = 0.3;
PFTs(23).Km_NO3                 = 0.3;
PFTs(23).Cmin_NH4               = 0.5;
PFTs(23).Cmin_NO3               = 0.25;
PFTs(23).SRL                    = 40;
PFTs(23).Vmax_NO3_LATS_umol_g_h = 0;
PFTs(23).Km_LATS_NO3            = 70;
%   [LEAF_LIFESPAN]
PFTs(23).leaf_long_days         = 200;
PFTs(23).retransl_fraction      = 0.7;
%   [RESPIRATION] Plant respiration (Sitch 2003 GCB 9:161 Eq 22-23, p.169)
PFTs(23).r_m10                = 0.066;     % [g C/(g N·d)] @ 10°C base, Sitch 2003 Table 3 "Plant Respiration" –TeBS-class –temperate broadleaf deciduous tree (apple, Wang 2009)
PFTs(23).R_g_fraction         = 0.11;      % [-] growth resp factor (CLM5 grperc; was Sitch 0.25).
% [resp scheme promoted 2026-06-15] tc_mm_nodebt selected as the winning maintenance-respiration
%   substrate scheme after the controlled comparison runs; the former env-switch (COMPASS_RESP_SCHEME
%   none|bruhn|tc_mm|tc_mm_nodebt|tc_mm_fastrec|jones) is RETIRED -- apple now always uses scheme 4.
%   Substrate Michaelis-Menten turn-down (Thornley & Cannell 2000 Ann Bot 85:55 Eq 17) on TOTAL NSC
%   (sugar + starch), with NO xsmrpool debt: when reserves run low the plant throttles maintenance
%   respiration rather than carrying a deficit. Consumed in applyNSCCoupling Step 1b/2 (scheme==4).
%   The other schemes' dead branches in applyNSCCoupling (xsmrpool debt / bruhn / jones) are pending a
%   separate cleanup (verify plant_c_residual after removal).
PFTs(23).resp_scheme = 4;   % tc_mm_nodebt (T&C 2000 Eq 17 turn-down on total NSC, no debt pool)
PFTs(23).resp_Km     = 2;   % [g C/m2] M-M half-saturation for the substrate turn-down
PFTs(23).resp_jones_kf = 0.227;   % [1/hr] Jones et al. 2024 fast-pool turnover (50th pct)
PFTs(23).resp_jones_Yg = 0.75;    % [-] McCree 1982 Crop Sci 22 growth yield (Jones et al. 2024 default)
% [labile sugar floor 2026-06-14] Soluble-sugar working-pool floor (computeLabileSugarTarget Level 2).
%   labile_resp_days = days of potential maintenance respiration held as soluble sugar; calibration knob
%   (sensitivity set [0.25 0.5 1 2], docs/nsc_labile_sugar_floor_design.md). Env override for sweeps.
PFTs(23).labile_resp_days = localLabileRespDays();   % [d] env-overridable (COMPASS_LABILE_RESP_DAYS); see cache note at top
PFTs(23).C_labile_abs_min = 0.01;   % [g C/m2] small nonzero floor
%   [TURNOVER_HARVEST] Turnover + Harvest + Spring release (M2-B.4 2026-05-21)
%   --------------------------------------------------------------------------
%   --------------------------------------------------------------------------
%   frPHU
%
%
%     1.
%     2. spring_release_endphu = 0.20 × frPHU_yr_end
%     3. harvest_threshold     = 0.90 × frPHU_yr_end
%     4. DLAI
%
%   --------------------------------------------------------------------------
%
%     ·
%     ·
%     · ∝
%
%   PHU=3921
%   PHU
%
%   "
%
%     EPIC
%     -----------------------------------------------------------------
%     spring_release 0.20   → 0.13  (= 0.20 × 0.65)     → ~DOY 107 (4
%     DLAI            0.45  → 0.592 (LAI
%     harvest        0.90  → 0.60  (= 0.90 × 0.65)     → ~DOY 285 (10
%
%
PFTs(23).harvest_fraction         = 0.90;     % [-, site] N_fruit removed at harvest (Wang/Mao commercial ~80-95%)
PFTs(23).root_turnover_rate       = 0.00274;  % [1/day, site] Annual ≈1/365 (Wells & Eissenstat 2002 Apple fine root)
PFTs(23).spring_release_days      = 60;       % [days, site] N_storage →N_labile decay constant

% ==============================================================================
% PFT 24: Robinia pseudoacacia –Black Locust (Robinia pseudoacacia)
% ==============================================================================
%   f_livewood/CN/
PFTs(24).PFT                    = 'Robinia pseudoacacia';
PFTs(24).FullName               = 'Black Locust (Robinia pseudoacacia)';
PFTs(24).clumping_index         = 0.70;  % [-] Proxy: ClimaLand temperate broadleaf-deciduous tree; not a species measurement
%   [HABIT] (PFT lifeform flags)
PFTs(24).photo_type             = 'C3';
PFTs(24).is_perennial           = 1;
PFTs(24).has_fruit              = 0;
PFTs(24).has_storage_pool       = 1;
%   [PHENOLOGY] (see [N.0]-[N.5] sub-stages below)
PFTs(24).evergreen            = 0;
PFTs(24).season_decid         = 1;
PFTs(24).stress_decid         = 0;
PFTs(24).crop_phenology       = 0;   % Robinia
%   [N.1] Sykes 1996 chilling-GDD for season_decid (LPJ-GUESS growth.cpp:91-97)
PFTs(24).k_chilla            = 0;       % [°C·d] Sykes 1996 a, LPJ-GUESS default
PFTs(24).k_chillb            = 100;     % [°C·d] Sykes 1996 b, species-specific (LPJ-GUESS europe.ins)
PFTs(24).k_chillk            = 0.05;    % [d^-1] Sykes 1996 k, LPJ-GUESS default
PFTs(24).phengdd5ramp        = 200;     % [°C·d] linear ramp to full leaf-on (LPJ-GUESS)
%   [N.2] Senescence (CLM5 daylength trigger, CNPhenologyMod.F90:957)
PFTs(24).crit_dayl_h         = 10.9;    % [h] dayl threshold for offset trigger (CLM5 global default)
PFTs(24).ndays_on            = 30;      % [days] CLM5 onset ramp duration (CNPhenologyMod.F90:51, clm5_params.c171117.nc scalar)
PFTs(24).ndays_off           = 15;      % [days] CLM5 offset ramp duration (CNPhenologyMod.F90:52, clm5_params.c171117.nc scalar)
PFTs(24).fstor2tran          = 0.5;     % [-] CLM5 budburst fraction (CNPhenologyMod.F90:53, clm5_params.c171117.nc scalar = 0.5)
PFTs(24).GDD_bloom_after_budburst = NaN;     % [°C·d] only Apple PFT 23 has fruit phenology calibrated for Changwu
PFTs(24).TSum_fruit              = NaN;
PFTs(24).harvest_DS              = NaN;
%   [N.3] N retranslocation on tissue turnover (LPJ-GUESS growth.cpp:256-285)
PFTs(24).nrelocfrac          = 0.5;     % [-] N retransloc fraction to N_storage at senescence
%   [STOMATA]
PFTs(24).m                      = 9;
PFTs(24).alpha                  = 0.06;
%   [LEAF]
PFTs(24).Ztop                   = 0;
PFTs(24).Zbot                   = 0;
PFTs(24).SLA0                   = 0.03;
PFTs(24).SLAm                   = 0.004;
PFTs(24).CN_L                   = 25;
%   [PHOTOSYNTH]
PFTs(24).Vcmax25                = 109;
PFTs(24).g_m                    = 0.175;  % [mol m-2 s-1] gm25 (Knauer 2019 GCB25 Table 3 DBF, deciduous broadleaf; robinia is a deciduous broadleaf tree)
%   [LUNA]
PFTs(24).Jmaxb1                 = 0.0936;
%   [ALLOC_CN]
PFTs(24).CN_target_root         = 35;
PFTs(24).CN_target_stem         = 100;
% [13c dual stem pool] CLM5 CN_Allocation Tbl col7-8 (CN_lw/CN_dw); clm5_params (flivewd/lwtop_ann)
PFTs(24).CN_lw                  = 50;    % [gC/gN] live stem (respires)
PFTs(24).CN_dw                  = 500;   % [gC/gN] dead stem (heartwood, ~N-free)
PFTs(24).flivewd                = 0.1;  % [-] new-stem-C frac to live wood (tree/shrub .1/grass 0/crop 1)
% croot_stem: coarse-root:stem allocation ratio -> routes coarse-root share of stem
%   allocation into belowground (root) reserve at storage (CLM5
%   NutrientCompetitionFlexibleCNMod.F90:529, natural tree/shrub value 0.30; clm5_params
%   .c171117__FruitTree.nc verified 2026-06-07). Consumed in computeCarbonAllocation.
PFTs(24).croot_stem             = 0.30;  % [-] coarse-root:stem ratio (CLM5 natural tree/shrub)
PFTs(24).froot_leaf             = 1.5;   % [-] f1 fine-root:leaf  (c_allometry alloc; CLM5 nc idx 7, BDT temperate)
PFTs(24).croot_stem_alloc       = 0.30;  % [-] f2 coarse:stem ALLOC flux (c_allometry; vs reserve croot_stem)
PFTs(24).stem_leaf              = 2.3;   % [-] f3 stem:leaf  (c_allometry alloc; CLM5 nc idx 7, BDT temperate)
PFTs(24).lwtop                  = 0.7;   % [1/yr] live->dead wood turnover (CLM5 lwtop_ann)
PFTs(24).CN_target_fruit        = 0;
% [fallback only 2026-06-28] BDT routes normal NPP/N + spring onset through the
%   c_allometry branch (froot_leaf=1.5/stem_leaf=2.3), so f_alloc_* below is NOT the
%   main allocation path. Kept as a finite fallback/legacy bridge (not NaN, to avoid
%   poisoning any fruittree-off test/legacy path); do NOT read it as the BDT allocation.
PFTs(24).f_alloc_leaf           = 0.4;
PFTs(24).f_alloc_root           = 0.25;
PFTs(24).f_alloc_stem           = 0.35;
PFTs(24).f_alloc_fruit          = 0;
%   [N_UPTAKE] (Nacry kinetics)
PFTs(24).Vmax_NH4_umol_g_h      = 12;
PFTs(24).Vmax_NO3_umol_g_h      = 5;
PFTs(24).Km_NH4                 = 0.3;
PFTs(24).Km_NO3                 = 0.3;
PFTs(24).Cmin_NH4               = 0.5;
PFTs(24).Cmin_NO3               = 0.25;
PFTs(24).SRL                    = 40;
PFTs(24).Vmax_NO3_LATS_umol_g_h = 0;
PFTs(24).Km_LATS_NO3            = 70;
%   [LEAF_LIFESPAN]
PFTs(24).leaf_long_days         = 200;
PFTs(24).retransl_fraction      = 0.7;
%   [RESPIRATION] Plant respiration (Sitch 2003 GCB 9:161 Eq 22-23, p.169)
PFTs(24).r_m10                = 0.066;     % [g C/(g N·d)] @ 10°C base, Sitch 2003 Table 3 "Plant Respiration" –TeBS-class –temperate broadleaf deciduous tree (Robinia)
PFTs(24).R_g_fraction         = 0.11;      % [-] growth resp factor (CLM5 grperc; was Sitch 0.25).
PFTs(24).harvest_fraction         = 0.90;
PFTs(24).root_turnover_rate       = 0.00274;
PFTs(24).spring_release_days      = 60;

% ==============================================================================
% PFT 25: Winter Wheat
% ==============================================================================
% Provenance: CoLM MOD_BGC_Veg_CNPhenology.F90 and
% MOD_BGC_Veg_NutrientCompetition.F90, which retain the CLM5 crop equations.
% Winter-wheat equations are dispatched independently; Apple/FruitTree remains
% unchanged. Parameter values marked proxy still require Changwu calibration.
PFTs(25).PFT                    = 'Winter Wheat';
PFTs(25).FullName               = 'Winter Wheat';
PFTs(25).clumping_index         = 0.76;  % [-] C3 grass/crop proxy
%   [HABIT]
PFTs(25).photo_type             = 'C3';
PFTs(25).is_perennial           = 0;
PFTs(25).has_fruit              = 1;
PFTs(25).has_storage_pool       = 0;
%   [PHENOLOGY CLASS]
PFTs(25).evergreen              = 0;
PFTs(25).season_decid           = 0;
PFTs(25).stress_decid           = 0;
PFTs(25).crop_phenology         = 1;
PFTs(25).winter_wheat_cn        = 1;
PFTs(25).k_chilla               = NaN;
PFTs(25).k_chillb               = NaN;
PFTs(25).k_chillk               = NaN;
PFTs(25).phengdd5ramp           = NaN;
PFTs(25).crit_dayl_h            = 10.9;
PFTs(25).ndays_on               = 30;
PFTs(25).ndays_off              = 15;
PFTs(25).fstor2tran             = 0.5;
PFTs(25).GDD_bloom_after_budburst = NaN;
PFTs(25).TSum_fruit               = NaN;
PFTs(25).harvest_DS               = NaN;
PFTs(25).nrelocfrac               = 0.5;
PFTs(25).baset                    = 0.0;  % [degC] wheat GDD base, CLM5 crop convention
%   [STOMATA / LEAF]
PFTs(25).m                      = 9.0;    % [-] supplied CLM5 mbbopt
PFTs(25).alpha                  = 0.06;
PFTs(25).Ztop                   = 0.5;
PFTs(25).Zbot                   = 0.01;
PFTs(25).SLA0                   = 0.035;  % [m2 leaf g-1 C] supplied CLM5 slatop
PFTs(25).SLAm                   = 0.0;    % [m2 leaf g-1 C per LAI] supplied CLM5 dsladlai
PFTs(25).CN_L                   = 20.0;   % [g C g-1 N] supplied CLM5 leafcn
%   [PHOTOSYNTH / LUNA]
PFTs(25).Vcmax25                = 100.7;  % [umol CO2 m-2 s-1] Spring Wheat proxy
PFTs(25).g_m                    = 0.295;  % [mol m-2 s-1] C3-crop proxy
PFTs(25).Jmaxb1                 = 0.0936;
%   [ORGAN C:N]
PFTs(25).CN_target_root         = 42.0;   % [g C g-1 N] supplied CLM5 frootcn
PFTs(25).CN_target_stem         = 50.0;   % [g C g-1 N] supplied CLM5 livewdcn
PFTs(25).CN_lw                  = 50.0;
PFTs(25).CN_dw                  = 500.0;
PFTs(25).flivewd                = 1.0;
PFTs(25).lwtop                  = 0.7;
PFTs(25).CN_target_fruit        = 20.0;   % [g C g-1 N] prescribed Changwu grain C:N target
% Change: 2026-09-30 Wenjie Liu -- explicit winter-wheat grain-fill donor
% C:N limits, CoLM-SYSU/CoLM commit 626347a9156ab97db4ca846f47cd1b7f2728f705,
% MOD_Const_PFT.F90 winter-wheat index 21 (fleafcn/fstemcn/ffrootcn).
% These are reference defaults, not Changwu-fitted measurements. The actual-N
% one-way adapter cannot export N from a root already below C/40.
PFTs(25).CN_leaf_grainfill      = 65.0;
PFTs(25).CN_stem_grainfill      = 100.0;
PFTs(25).CN_root_grainfill      = 40.0;
%   [LEGACY FALLBACK ALLOCATION]
% The active winter-wheat path uses the CoLM HUI-dependent coefficients below.
PFTs(25).f_alloc_leaf           = 0.45;
PFTs(25).f_alloc_root           = 0.30;
PFTs(25).f_alloc_stem           = 0.15;
PFTs(25).f_alloc_fruit          = 0.10;
%   [N UPTAKE]
% C3-crop proxy; these kinetic values are not supplied by the CLM5 crop table.
PFTs(25).Vmax_NH4_umol_g_h      = 20;
PFTs(25).Vmax_NO3_umol_g_h      = 15;
PFTs(25).Km_NH4                 = 0.3;
PFTs(25).Km_NO3                 = 0.3;
PFTs(25).Cmin_NH4               = 0.5;
PFTs(25).Cmin_NO3               = 0.25;
PFTs(25).SRL                    = 50;
PFTs(25).Vmax_NO3_LATS_umol_g_h = 0;
PFTs(25).Km_LATS_NO3            = 70;
%   [TURNOVER / RESPIRATION]
PFTs(25).leaf_long_days         = 365;
PFTs(25).retransl_fraction      = 0.4;
PFTs(25).r_m10                  = 0.066;   % legacy fallback; active crop path uses CLM br/Q10 below
PFTs(25).R_g_fraction           = 0.25;    % [-] LPJ/Sitch fraction of substrate after maintenance; not CoLM grperc
PFTs(25).clm_crop_respiration   = 1;       % [0/1] use CLM5/CoLM crop maintenance respiration
PFTs(25).clm_br_s               = 2.525e-6;% [g C gN-1 s-1] Ryan 1991 conversion used by CLM5/CoLM
PFTs(25).clm_br_root_s          = 2.525e-6;% [g C gN-1 s-1] CLM5/CoLM belowground base rate
PFTs(25).clm_resp_Q10           = 1.5;     % [-] CLM5/CoLM maintenance-respiration temperature factor
PFTs(25).resp_scheme            = 0;       % potential Rm unchanged; winter_wheat_cn enforces physical NSC funding
% Change: 2026-09-30 Wenjie Liu -- clarify that this is the potential
% respiration scheme, not permission to create a signed crop carbon debt.
PFTs(25).xsmr_recovery_days     = 30.0;    % [d] CLM5/CoLM deficit recovery timescale
PFTs(25).harvest_fraction       = 0.90;
PFTs(25).root_turnover_rate     = 0.00274;
PFTs(25).spring_release_days    = 60;
%   [ACTIVE LOCAL CoLM WINTER-CROP PARAMETERS]
% Change: 2026-09-30 Wenjie Liu -- use the verified local CoLM
% MOD_Const_PFT.F90 winter_wheat index 21 reference defaults (CROP enabled),
% identical to commit 626347a9156ab97db4ca846f47cd1b7f2728f705.
% Source declarations: arooti L1376, fleafi L1406, bfact L1421,
% grnfill L1330, mxmat L1315. These are model reference values, not
% Changwu measurements; observed_lai_zero still owns physical harvest.
PFTs(25).arooti                 = 0.10;
PFTs(25).arootf                 = 0.00;
PFTs(25).astemf                 = 0.05;
PFTs(25).allconss               = 1.00;
PFTs(25).bfact                  = 0.10;
PFTs(25).fleafi                 = 0.90;
PFTs(25).lfemerg                = 0.03;
PFTs(25).grnfill                = 0.67;
PFTs(25).mxmat                  = 270;
PFTs(25).min_NH_planting_date   = 901;
PFTs(25).max_NH_planting_date   = 1130;
PFTs(25).min_SH_planting_date   = 301;
PFTs(25).max_SH_planting_date   = 530;
PFTs(25).min_planting_temp_K    = 283.15;
PFTs(25).planting_temp_K        = 1000.0;
PFTs(25).initial_seed_c         = 3.0;   % [g C m-2], CoLM CropPhenology
PFTs(25).croot_stem_alloc       = 0.0;
PFTs(25).use_fruittree_alloc    = 0;

% Change: 2026-08-24 Jizixing Zhu -- map Winter Wheat to the fixed CLM5/CoLM
% crop respiration, growth-respiration, and cpool-xsmrpool parameters.

CACHED_PFTS = PFTs;   % [P1.5-b] cache the built constant table for subsequent calls
end

function v = localLabileRespDays()
% Env-overridable labile-respiration-days knob (COMPASS_LABILE_RESP_DAYS sweep);
% default 1.0 d. Single source shared by the table build and the cache-hit path.
    lrd_env = str2double(getenv('COMPASS_LABILE_RESP_DAYS'));
    if ~isnan(lrd_env) && lrd_env > 0
        v = lrd_env;
    else
        v = 1.0;
    end
end
