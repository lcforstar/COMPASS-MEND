function vod = compute_VOD(state, params)
%COMPUTE_VOD  Component VOD at L-band, prefers prognostic W storage when available.
%
%   Returns:
%     vod.total     [-]  Total L-band vegetation optical depth (nadir)
%     vod.leaf      [-]  Leaf component
%     vod.stem      [-]  Stem/wood component
%     vod.int       [-]  Intercepted-water component
%     vod.b_leaf_eff [-] Effective b_leaf after ψ_leaf and freeze modulation
%     vod.source    str  'prognostic_W' or 'empirical_LAI_AGB' (which VWC path used)
%
%   Inputs:
%     state.LAI            [m²/m²]   Leaf area index
%     state.h2ocan         [mm]      Canopy intercepted water
%     state.psi_leaf_MPa   [MPa]     Leaf water potential
%     state.Tc_K           [K]       Canopy temperature
%
%     [Optional, preferred path] If prognostic v2 plant water is available,
%     pass these and they override the empirical VWC formulas:
%     state.W_stem_m       [m]       Prognostic stem water storage (per m² ground)
%     state.W_leaf_m       [m]       Prognostic leaf water storage
%
%     params.calib.b_leaf      [-]        Reference leaf b coefficient at ψ_ref
%     params.calib.b_stem      [-]        Stem b coefficient
%     params.calib.b_int       [-]        Intercepted-water b
%     params.calib.alpha_psi   [MPa⁻¹]    Slope of b_leaf vs ψ_leaf
%     params.calib.psi_ref     [MPa]      Reference ψ_leaf
%     params.calib.f_frozen    [-]        L-band ice multiplier
%     params.calib.stem_dynamic_fraction [-] Fraction of stem tissue water
%                                      that follows capacitance RWC.
%     params.site.Cw           [g/cm²]    [empirical path only] leaf water thickness
%     params.site.AGB_static   [kg/m²]    [empirical path only] static AGB
%
%   Method:
%     Two VWC paths, chosen automatically:
%
%     Path A — prognostic (recommended, used when W_stem_m / W_leaf_m present):
%       VWC_leaf [kg/m²] = W_leaf_m × 1000   (1 m water column = 1000 kg/m²)
%       VWC_stem [kg/m²] = W_stem_m × 1000
%       VWC_int  [kg/m²] = h2ocan_mm
%       Source: storage_from_psi inverse retention W ← v1 ψ (updateOrganStorageFromPsi)

%
%     Path B — empirical (legacy fallback, used when W fields absent):
%       VWC_leaf = LAI · Cw · 10     (Fluspect convention)
%       VWC_stem = AGB_static        (constant)
%       VWC_int  = h2ocan
%
%     Stem-water closure:
%       Woody L-band VOD includes structural/free tissue water and a smaller
%       exchangeable capacitance part. Therefore the RWC perturbation is
%       applied only to a configurable dynamic fraction of stem water:
%         M_stem = M_stem,max * [(1 - f_dyn) + f_dyn * RWC_stem]
%       rather than scaling the full woody water pool by RWC_stem.
%
%     b_leaf modulation (both paths):
%       b_leaf(ψ, T) = b_leaf_0 · g(ψ_leaf) · h(T_canopy)
%       g(ψ) = max(0.3, 1 + α·(ψ - ψ_ref))      Momen et al. 2017 JGR Eq 6
%       h(T) = f_frozen if T < 273.15 K else 1   Holtzman 2023 WRR Table 1
%
%   Design basis (project PPT slide 17, microwave VOD operator):
%     VOD_total = VOD_leaf + VOD_stem + VOD_int
%     VOD_organ = β_organ · M_organ · G_organ(θ, p)
%     where M_organ [kg/m²] is the organ water content taken DIRECTLY from
%     COMPASS prognostic plant hydraulics (not LAI proxy / static VWC), β is
%     the microwave mass extinction coefficient, G the angle/polarisation
%     factor. This prognostic-water driving is the core novelty.
%
%   References (verified against local readable PDFs only):
%     Holtzman et al. 2023 WRR — CliMA Land + plant hydraulics + VOD;
%       Eq 3 VOD = (1+a·ψ)·(b+c·LAI), first term ≡ VWC under linear P-V curve;
%       line 83: total canopy water = aboveground biomass × VWC.
%       (data/references/vod_literature/Holtzman_2023_WRR_CliMA_microwave.pdf)
%     [SOURCE UNVERIFIED] Wigneron 2007 L-MEB, Konings 2019, Momen 2017 are
%       relevant but their local PDFs are corrupt (Cloudflare HTML, not PDF);
%       β temperature dependence below follows Kirdyashev via compute_b_kirdyashev.

% ─── iBiomass switch: decide source of VWC (CMEM-style) ───────────────────
%   0 = static (AGB_static constant + LAI·Cw·10)           — backward compatible
%   1 = strict COMPASS CN pools (errors if absent)         — research/debug
%   2 = auto: use CN pools if available, else fall back    — recommended (default)
%
% VWC = wet vegetation mass (kg water / m²) ≠ dry biomass.
% Conversion: VWC = dry_biomass × (1 - m_d) / m_d
%   where m_d is the dry-matter fraction. For stem (m_d ≈ 0.5),
%   wet_mass / dry_mass = 2, so VWC = dry_biomass × 1.0 (water = dry).
%   For leaf (m_d ≈ 0.3), wet/dry = 3.33, so VWC = dry × 2.33.
% AGB_static = 5.0 kg/m² is already wet biomass (Konings 2019 convention),
% so the static path uses it directly. CN-pool path converts dry → wet
% using m_d_leaf and m_d_stem.
if ~isfield(params.calib, 'iBiomass'), params.calib.iBiomass = 2; end
if ~isfield(params.calib, 'b_int'), params.calib.b_int = 0.10; end

has_C_pools = isfield(state, 'plant_c_stem_kgm2') && isfield(state, 'plant_c_leaf_kgm2') ...
              && ~isempty(state.plant_c_stem_kgm2) && ~isempty(state.plant_c_leaf_kgm2);

% (m_d is consumed inside compute_b_kirdyashev for the Kirdyashev branch;
%  the WDR ratios below are species-specific tissue water-to-dry ratios.)
if ~isfield(params.calib, 'WDR_leaf'), params.calib.WDR_leaf = 2.5; end
if ~isfield(params.calib, 'WDR_stem'), params.calib.WDR_stem = 0.8; end
if ~isfield(params.calib, 'use_dual_stem_wdr'), params.calib.use_dual_stem_wdr = false; end
WDR_leaf = params.calib.WDR_leaf;     % [-] saturated leaf water/dry-mass ratio
WDR_stem = params.calib.WDR_stem;     % [-] saturated stem water/dry-mass ratio
stem_water_max = [];     % [kg/m2] stem water at RWC = 1

switch params.calib.iBiomass
    case 0
        agb_stem = params.site.AGB_static;
        agb_leaf_proxy = state.LAI * params.site.Cw * 10;
        stem_water_max = agb_stem;
        vod.biomass_source = 'static_AGB';
    case 1
        if ~has_C_pools
            error('iBiomass=1 requires plant_c_stem_kgm2 / plant_c_leaf_kgm2 from COMPASS CN');
        end
        % VWC = dry biomass x WDR (water-to-dry-mass ratio). WDR_stem
        % and WDR_leaf are configurable site parameters; use local tissue
        % measurements or constrained calibration when available.
        % Current closure: WDR_stem and WDR_leaf are configurable site
        % parameters. The old hard-coded stem WDR=1.5 is no longer used.
        agb_stem      = state.plant_c_stem_kgm2;
        stem_water_max = agb_stem * WDR_stem;
        agb_leaf_proxy = state.plant_c_leaf_kgm2 * WDR_leaf;
        vod.biomass_source = 'compass_cn_strict';
        vod.WDR_stem = WDR_stem; vod.WDR_leaf = WDR_leaf;
    case 2
        if has_C_pools
            % Use dual live/dead C pools for dry-mass accounting. Different
            % live/dead WDR values are used only when explicitly enabled;
            % otherwise dual stem pools are summed with uniform WDR.
            has_dual = isfield(state, 'plant_c_livestem_kgm2') && isfield(state, 'plant_c_deadstem_kgm2') ...
                       && ~isempty(state.plant_c_livestem_kgm2) && ~isempty(state.plant_c_deadstem_kgm2);
            if has_dual && params.calib.use_dual_stem_wdr
                if ~isfield(params.calib, 'WDR_live'), params.calib.WDR_live = WDR_stem; end
                if ~isfield(params.calib, 'WDR_dead'), params.calib.WDR_dead = WDR_stem; end
                WDR_live = params.calib.WDR_live;   % [-] live wood water/dry-mass ratio
                WDR_dead = params.calib.WDR_dead;   % [-] dead wood water/dry-mass ratio
                agb_stem = state.plant_c_livestem_kgm2 + state.plant_c_deadstem_kgm2;
                stem_water_max = state.plant_c_livestem_kgm2 * WDR_live ...
                               + state.plant_c_deadstem_kgm2 * WDR_dead;
                vod.biomass_source = 'compass_cn_auto_dualstem';
                vod.WDR_live = WDR_live; vod.WDR_dead = WDR_dead;
            else
                agb_stem = state.plant_c_stem_kgm2;
                stem_water_max = agb_stem * WDR_stem;
                vod.biomass_source = 'compass_cn_auto';
                vod.WDR_stem = WDR_stem;
            end
            agb_leaf_proxy = state.plant_c_leaf_kgm2 * WDR_leaf;
            vod.WDR_leaf = WDR_leaf;
            vod.use_dual_stem_wdr = has_dual && params.calib.use_dual_stem_wdr;
        else
            agb_stem = params.site.AGB_static;
            agb_leaf_proxy = state.LAI * params.site.Cw * 10;
            stem_water_max = agb_stem;
            vod.biomass_source = 'fallback_static_AGB';
        end
    otherwise
        error('compute_VOD: unsupported iBiomass = %d', params.calib.iBiomass);
end
vod.agb_stem = agb_stem;
vod.agb_leaf_proxy = agb_leaf_proxy;
vod.stem_water_max = stem_water_max;

% ─── M_organ: organ water content [kg/m²] ────────────────────────────────
% PPT slide 17: VOD_organ = β · M_organ · G. The full tissue water (structural
% + free) is what microwaves respond to, so M is the saturated tissue-water
% capacity split into a structural baseline and an exchangeable RWC perturbation:
%   M_stem = (AGB_stem_dry · WDR_max) · rwc_stem        [kg/m²]
%            └── max tissue water ──┘   └ COMPASS prognostic (0.83–0.95) ┘
% AGB_stem_dry·WDR_max sets the magnitude (a few kg/m², so β stays physical);
% rwc_stem (w_stem_v2/w_stem_max) perturbs only the configured f_dyn fraction.
% This avoids assigning the full seasonal capacitance signal to structural wood.
% Current code below uses f_dyn; older whole-pool RWC wording is superseded.
% This keeps the bulk structural water (≈90%) yet lets the WHOLE pool breathe
% with the true RWC — no empirical dormancy factor needed.
% Leaf carbon dry-mass [kg/m²] = C_leaf / carbon_fraction (0.466, CLM5). agb_leaf_proxy
% above already folded one WDR into it; recover dry mass so M_leaf uses the same
% (dry x WDR_max x rwc) form as M_stem and stays symmetric. WDR_leaf is
% configurable and defaults to 2.5.
WDR_leaf_max = WDR_leaf;   % [-] saturated leaf water/dry-mass ratio (carried from iBiomass block)
leaf_dry = agb_leaf_proxy / max(WDR_leaf, 1e-9);   % [kg/m²] back out dry leaf mass
use_rwc = isfield(state, 'rwc_stem') && ~isempty(state.rwc_stem) && ~isnan(state.rwc_stem);
use_rwc_leaf = isfield(state, 'rwc_leaf') && ~isempty(state.rwc_leaf) && ~isnan(state.rwc_leaf);
use_prognostic_W = isfield(state, 'W_stem_m') && isfield(state, 'W_leaf_m') ...
                   && ~isempty(state.W_stem_m) && ~isempty(state.W_leaf_m);
if ~isfield(params.calib, 'stem_dynamic_fraction')
    params.calib.stem_dynamic_fraction = 0.5;
end
stem_dynamic_fraction = min(max(params.calib.stem_dynamic_fraction, 0.0), 1.0);

if use_rwc
    stem_rwc_factor = (1 - stem_dynamic_fraction) + stem_dynamic_fraction * state.rwc_stem;
    VWC_stem = stem_water_max * stem_rwc_factor;             % [kg/m2] M_stem, structural + dynamic
    % M_leaf = leaf_dry x WDR_leaf_max x rwc_leaf, symmetric with M_stem: carbon-driven
    % saturated water modulated by capacitance RWC. Falls back to the (uncscaled) wet
    % proxy if rwc_leaf is unavailable (older runs without rwc_leaf_v2).
    if use_rwc_leaf
        VWC_leaf = leaf_dry * WDR_leaf_max * state.rwc_leaf; % [kg/m²] M_leaf (RWC-scaled)
        vod.rwc_leaf     = state.rwc_leaf;
        vod.VWC_leaf_max = leaf_dry * WDR_leaf_max;
    else
        VWC_leaf = agb_leaf_proxy;                           % [kg/m²] wet proxy (no rwc_leaf)
    end
    VWC_int  = state.h2ocan;                                 % [kg/m²] W_interception
    vod.source = 'rwc_scaled_biomass';
    vod.VWC_stem_max = stem_water_max;
    vod.rwc_stem     = state.rwc_stem;
    vod.stem_dynamic_fraction = stem_dynamic_fraction;
    vod.stem_rwc_factor = stem_rwc_factor;
elseif use_prognostic_W
    VWC_leaf = state.W_leaf_m * 1000;       % [kg/m²] M_leaf
    VWC_stem = state.W_stem_m * 1000;       % [kg/m²] M_stem (prognostic, seasonal)
    VWC_int  = state.h2ocan;                % [kg/m²] W_interception
    vod.source = 'prognostic_W';
else
    % ─── Fallback: empirical biomass proxy (no prognostic water available) ─
    VWC_leaf = agb_leaf_proxy;                       % [kg/m²]
    VWC_stem = stem_water_max;                        % [kg/m²]
    VWC_int  = state.h2ocan;                          % [kg/m²]
    vod.source = 'empirical_AGB_only';
end

% ─── iCanopy switch: choose source of mass extinction coefficients b ─────
%   0 = empirical (params.calib.b_leaf/b_stem, fmincon-tuned)   — default
%   1 = Kirdyashev physical (CMEM vegkird, uses Tc + m_d + LIDF)
%   (b_int stays empirical since intercepted water is free liquid water,
%    not bound in tissue — empirical b_int ≈ 0.5 from Saleh 2007.)
if ~isfield(params.calib, 'iCanopy'), params.calib.iCanopy = 0; end

% g(ψ): leaf-water-potential modulation of b_leaf (instantaneous, physical).
g_psi = max(0.3, 1 + params.calib.alpha_psi * ...
                     (state.psi_leaf_MPa - params.calib.psi_ref));
% h(T): freeze factor. DISABLED (=1.0) — at Changwu the observed VODCA-L is
% HIGHEST in winter (Nov-Dec ≈ 0.32-0.34), so a freeze cut (f_frozen=0.3) put
% a winter VOD collapse that does not exist in the data and flipped the
% seasonal phase to summer-high (R=-0.6 vs VODCA-L). Apple stems do not freeze
% to the dielectric-collapse degree the factor assumed. Kept switchable via
% params.calib.f_frozen for sites where winter freeze truly matters.
if isfield(params.calib,'f_frozen') && params.calib.f_frozen < 1 && state.Tc_K < 273.15
    h_T = params.calib.f_frozen;
else
    h_T = 1.0;
end

switch params.calib.iCanopy
    case 0
        % Empirical: use calibrated b_leaf, b_stem (current default)
        b_leaf_base = params.calib.b_leaf;
        b_stem_base = params.calib.b_stem;
        vod.canopy_source = 'empirical_b';
    case 1
        % Kirdyashev physical b (CMEM vegkird), per-organ a_geo:
        %   b = a_geo · k · Im(ε_veg(T_C, m_d)) / cosθ / ρ_w
        % Use:
        %   Tc_C       from state.Tc_K
        %   m_d_leaf   from site.m_d (PROSPECT Cdm/(Cdm+Cw); typical leaf 0.3)
        %   m_d_stem   = 0.5 (woody tissue, Wegmueller 1990)
        %   a_geo_leaf = 2/3 (spherical/horizontal broadleaf, CMEM veg_sub.F90
        %                     header "isotropic leaves: a_h=2/3, a_v=2/3")
        %   a_geo_stem = 1.0 (vertical-leaning branches, CMEM table
        %                     "vertical stalks: a_h=0, a_v=sin²θ → close to 1.0)
        %                     For apple orchard with upright branches, this is
        %                     the physically appropriate per-organ structure.
        Tc_C = state.Tc_K - 273.15;
        if isfield(params.site, 'm_d') && ~isempty(params.site.m_d)
            m_d_leaf = params.site.m_d;
        else
            m_d_leaf = 0.3;
        end
        m_d_stem  = 0.5;
        a_geo_leaf = 2/3;
        a_geo_stem = 1.0;
        b_leaf_base = surface.radiation.microwave.compute_b_kirdyashev( ...
                        Tc_C, m_d_leaf, params.calib.Fre, params.calib.Inc, a_geo_leaf);
        b_stem_base = surface.radiation.microwave.compute_b_kirdyashev( ...
                        Tc_C, m_d_stem, params.calib.Fre, params.calib.Inc, a_geo_stem);
        % Structural enhancement factors (Konings 2019). Kirdyashev effective-
        % medium under-predicts b for both leaf and stem (ignores fiber
        % directionality, cell-wall & vascular scattering). Empirical b from
        % Wigneron 2007 / Konings 2019 ≈ 0.06-0.10 vs Kirdyashev ≈ 0.04.
        % Per-organ scale factors let calibration restore the right magnitude
        % while keeping the physical b(T) seasonal/diurnal shape.
        if ~isfield(params.calib, 'b_stem_scale'), params.calib.b_stem_scale = 1.0; end
        if ~isfield(params.calib, 'b_leaf_scale'), params.calib.b_leaf_scale = 1.0; end
        b_stem_base = b_stem_base * params.calib.b_stem_scale;
        b_leaf_base = b_leaf_base * params.calib.b_leaf_scale;
        vod.canopy_source = 'kirdyashev_physical';
        vod.a_geo_leaf = a_geo_leaf;
        vod.a_geo_stem = a_geo_stem;
        vod.b_stem_scale = params.calib.b_stem_scale;
        vod.b_leaf_scale = params.calib.b_leaf_scale;
    otherwise
        error('compute_VOD: unsupported iCanopy = %d', params.calib.iCanopy);
end

b_leaf_eff = b_leaf_base * g_psi * h_T;
% The seasonal VOD cycle comes from prognostic M_stem (rwc_stem, summer-low/
% winter-high) — no dormancy factor and no freeze factor (h_T=1) needed; the
% real water content carries the phenological signal and the observed phase.
b_stem_eff = b_stem_base * h_T;

% ─── Nadir VOD components ─────────────────────────────────────────────────
vod.leaf  = b_leaf_eff      * VWC_leaf;
vod.stem  = b_stem_eff      * VWC_stem;
vod.int   = params.calib.b_int  * VWC_int;
vod.total = vod.leaf + vod.stem + vod.int;

% Diagnostics
vod.b_leaf_base = b_leaf_base;
vod.b_stem_base = b_stem_base;

% ─── Polarization-split VOD (Wigneron 2007 L-MEB Eq 3 / CMEM vegwign) ─────
% Anisotropy factors tt_H, tt_V encode canopy structure via the relation
%     τ_p(θ) = τ_N · (cos²θ + tt_p · sin²θ) / cosθ
% Theoretical limits (CMEM veg_sub.F90 header lines 19–26):
%   horizontal leaves:  tt_H = 1,    tt_V = (cosθ/cosθ)²·a_geo_v ≈ 0
%   isotropic (sphere): tt_H = 1,    tt_V = 1     ← classical Wigneron, no anisotropy
%   vertical stalks:    tt_H = 0,    tt_V = high
% For C3 broadleaf crops Wigneron 2007 recommends tt_H = 1.0, tt_V = 2.0.
% Apple orchard (broadleaf, nearly horizontal) → use tt_H ≈ 1, tt_V ≈ 1.5 prior.
if ~isfield(params.calib, 'tt_H'), params.calib.tt_H = 1.0; end
if ~isfield(params.calib, 'tt_V'), params.calib.tt_V = 1.0; end

vod.tt_H    = params.calib.tt_H;
vod.tt_V    = params.calib.tt_V;
vod.tau_N   = vod.total;                                  % nadir, divided by cosθ later

% Diagnostics
vod.b_leaf_eff = b_leaf_eff;
vod.g_psi      = g_psi;
vod.h_T        = h_T;
vod.VWC_leaf   = VWC_leaf;
vod.VWC_stem   = VWC_stem;
vod.VWC_int    = VWC_int;

end
