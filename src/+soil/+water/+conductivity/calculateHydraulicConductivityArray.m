function SoilVariables = calculateHydraulicConductivityArray(SoilVariables, VanGenuchten, ~, L_f, cfg, Delt_t, TransportCoefficient)
%CALCULATEHYDRAULICCONDUCTIVITYARRAY  Compute coupled soil hydraulic properties.
%
% Updates liquid water, total water, ice content, hydraulic derivatives,
% effective saturation, liquid conductivity, freeze-impeded conductivity,
% temperature-gradient conductivity, ice ratio, and Gamma_hh for the
% fully coupled van Genuchten + freezing + strict SFCC + thermal-response
% path.
%
%
% Validity: implements only the fully coupled main path. The caller
% (UpdateSoilWaterContent > useArrayHydraulicKernel) guarantees this is invoked
% ONLY when lSWCC=1 (van Genuchten), lFreezing=1, lSFCC=1, and lThmrlefc=1.
% Other switch combinations use the general hydraulic routine. The
% entry assert below fails loud if a misroute ever violates that contract.
%
% Inputs:
%   SoilVariables        [inout] soil state and hydraulic auxiliary fields
%   VanGenuchten         [in]    VG parameters: theta_s, theta_r, alpha, n, m
%   L_f                  [in]    latent heat of freezing [J/kg]
%   cfg                  [in]    model settings; requires lSWCC/lFreezing/lSFCC/lThmrlefc = 1
%   Delt_t               [in]    outer time-step length [s]
%   TransportCoefficient [in]    thermal capacity and conductivity fields
%
% Outputs:
%   Writes
%   back Theta_LL, Theta_UU, Theta_II, DTheta_LLh, DTheta_UUh, Se, KL_h, KfL_h,
%   KfL_T, Ratio_ice, Gamma_hh.
%
% Method:
%   1. Apply thermal head correction and build layer-endpoint heads.
%   2. Compute Gamma_hh, Theta_m, Theta_UU, Theta_LL, and Theta_II from the
%      VG/SFCC water-retention relations.
%   3. Apply Dall'Amico thermal-time phase-change relaxation when enabled.
%   4. Compute dTheta/dh, ice ratio, viscosity-corrected KL_h, KfL_h,
%      KfL_T, and Se from the post-phase-change hydraulic state.
%
% References:
%   van Genuchten (1980) Soil Sci. Soc. Am. J., 44, 892-898.
%   Mualem (1976) Water Resour. Res., 12, 513-522.
%   Lundin (1990) Cold Reg. Sci. Technol., 18, 47-64.
%
% Author: Zhu Jizixing, 2026-06-15.

    phys = equations.physicalConstants();
    NL   = cfg.Grid.NL;
    nD   = cfg.Grid.nD;
    if nargin < 6 || isempty(Delt_t)
        Delt_t = 0;
    end
    if nargin < 7
        TransportCoefficient = [];
    end

    % This formulation is defined only for the fully coupled path.
    assert(cfg.lSWCC == 1 && cfg.lFreezing == 1 && cfg.lSFCC == 1 && cfg.lThmrlefc == 1, ...
        'calculateHydraulicConductivityArray:unsupportedConfig', ...
        ['Coupled hydraulic-property update requires lSWCC/lFreezing/lSFCC/lThmrlefc = 1; ', ...
         'use the general hydraulic routine for any other config.']);

    dallAmicoEnabled = cfg.lFreezing == 1 && cfg.lThmrlefc == 1 && Delt_t > 0 && ...
        getLogicalOption(cfg, 'lDallAmicoPhaseChange', false) && ...
        hasThermalTimeInputs(TransportCoefficient);

    SoilConstants = io.getSoilConstants();
    hd = SoilConstants.hd;   % [cm] oven-dry reference head (-1e7)
    hm = SoilConstants.hm;   % [cm] moist-zone upper bound (-9899)

    % =====================================================================
    % Build layer-endpoint hydraulic states.
    %   Node-indexed fields map to (i,j) via MN = i + j - 1 (column j slices
    %   nodes j : j+NL-1). Layer-indexed VG parameters and porosity are applied
    %   to both endpoints of each soil layer.
    % =====================================================================
    hh_a      = zeros(NL, nD);   % [cmH2O] clamped head (fixHeat)
    hh_frez_a = zeros(NL, nD);   % [cmH2O] freezing depression head, end of step
    for j = 1:nD
        nodes           = j:(j + NL - 1);
        hh_a(:, j)      = fixHeatArray(SoilVariables.hhU(nodes), cfg);
        hh_frez_a(:, j) = SoilVariables.hh_frez(nodes);
    end

    one_row   = ones(1, nD);
    theta_s_a = VanGenuchten.Theta_s(:) * one_row;   % [m3/m3]
    theta_r_a = VanGenuchten.Theta_r(:) * one_row;   % [m3/m3]
    alpha_a   = VanGenuchten.Alpha(:)   * one_row;   % [1/cm]
    n_a       = VanGenuchten.n(:)       * one_row;   % [-]
    m_a       = VanGenuchten.m(:)       * one_row;   % [-]
    POR_a     = SoilVariables.POR(:)    * one_row;   % [m3/m3]

    % ---- Gamma_hh: liquid-water availability for vapor diffusion -------------
    % calculateGamma_hh: 3 segments; |hh|>=|hd| dry (0); |hm|<=|hh|<|hd| log
    % transition; |hh|<|hm| moist (1).
    absh    = abs(hh_a);
    Gamma_a = ones(NL, nD);
    Gamma_a(absh >= abs(hd)) = 0;
    mask_tr = (absh >= abs(hm)) & (absh < abs(hd));
    Gamma_a(mask_tr) = log(abs(hd) ./ absh(mask_tr)) ./ log(abs(hd) / abs(hm));

    % ---- Theta_m: VG capillary upper bound at hh=-1, clamped to [theta_s,POR] -
    % calculateTheta_m: Theta_m = Gamma*thr + (ths - Gamma*thr)*(1+|alpha*(-1)|^n)^(-m)
    Theta_m_a = Gamma_a .* theta_r_a + (theta_s_a - Gamma_a .* theta_r_a) .* ...
        (1 + abs(alpha_a .* (-1)).^n_a).^(-m_a);
    hi_m = Theta_m_a >= POR_a;                    % apply pore-volume upper bound first
    lo_m = ~hi_m & (Theta_m_a <= theta_s_a);      % then enforce saturated-water lower bound
    Theta_m_a(hi_m) = POR_a(hi_m);
    Theta_m_a(lo_m) = theta_s_a(lo_m);

    % ---- Theta_UU: total water content (calculateTheta_UU, lSFCC=1 path) ------
    %   hh>=-1 -> theta_s ; hh<=hd -> theta_r ; else (lThmrlefc=1):
    %   Gamma==0 -> 0 ; else calculateTheta(case 0) with heat_term = hh.
    sat_uu  = hh_a >= -1;
    dry_uu  = ~sat_uu & (hh_a <= hd);
    zero_uu = ~sat_uu & ~dry_uu & (Gamma_a == 0);
    vg_uu   = ~sat_uu & ~dry_uu & (Gamma_a ~= 0);
    tuu_vg  = Gamma_a .* theta_r_a + (Theta_m_a - Gamma_a .* theta_r_a) ./ ...
        (1 + abs(alpha_a .* hh_a).^n_a).^m_a;
    Theta_UU_a = zeros(NL, nD);
    Theta_UU_a(sat_uu)  = theta_s_a(sat_uu);
    Theta_UU_a(dry_uu)  = theta_r_a(dry_uu);
    Theta_UU_a(zero_uu) = 0;
    Theta_UU_a(vg_uu)   = tuu_vg(vg_uu);

    % ---- Theta_LL: liquid water content (calculateTheta_LL, lSFCC=1 path) -----
    %   hh>=-1 & hh_frez>=0 -> theta_s ; (hh+hh_frez)<=hd -> theta_r ;
    %   else (lThmrlefc=1): calculateTheta(case 0) with heat_term = hh+hh_frez.
    heatTerm = hh_a + hh_frez_a;
    sat_ll = (hh_a >= -1) & (hh_frez_a >= 0);
    dry_ll = ~sat_ll & (heatTerm <= hd);
    vg_ll  = ~sat_ll & ~dry_ll;
    tll_vg = Gamma_a .* theta_r_a + (Theta_m_a - Gamma_a .* theta_r_a) ./ ...
        (1 + abs(alpha_a .* heatTerm).^n_a).^m_a;
    Theta_LL_a = zeros(NL, nD);
    Theta_LL_a(sat_ll) = theta_s_a(sat_ll);
    Theta_LL_a(dry_ll) = theta_r_a(dry_ll);
    Theta_LL_a(vg_ll)  = tll_vg(vg_ll);

    % ---- Ice from liquid-phase mass balance (lFreezing=1) --------------------
    Theta_II_a = max(0, (Theta_UU_a - Theta_LL_a) * phys.RHOL / phys.RHOI);

    % =====================================================================
    % Endpoint temperatures for phase change, viscosity, and KfL_T.
    % MN = i + j - 1.
    % =====================================================================
    TT_node      = SoilVariables.TT(:);
    TT_CRIT_node = SoilVariables.TT_CRIT(:);
    TT_e      = zeros(NL, nD);
    TT_CRIT_e = zeros(NL, nD);
    for j = 1:nD
        TT_e(:, j)      = TT_node(j:(j + NL - 1));
        TT_CRIT_e(:, j) = TT_CRIT_node(j:(j + NL - 1));
    end

    % =====================================================================
    % Dall'Amico thermal-time phase-change relaxation.
    % The water-retention derivatives, conductivity, and Se below use the
    % post-phase-change liquid and ice contents. EHCAP/c_unsat/Lambda_eff are
    % layer-endpoint fields on this coupled path; rhoC uses EHCAP where valid
    % and c_unsat where EHCAP is unavailable.
    % =====================================================================
    if dallAmicoEnabled
        if isfield(TransportCoefficient, 'EHCAP') && ...
                size(TransportCoefficient.EHCAP, 1) >= NL && size(TransportCoefficient.EHCAP, 2) >= nD
            rhoC_a = TransportCoefficient.EHCAP(1:NL, 1:nD);
        else
            rhoC_a = NaN(NL, nD);
        end
        bad = ~isfinite(rhoC_a) | (rhoC_a <= 0);
        if any(bad(:))
            cu = TransportCoefficient.c_unsat(1:NL, 1:nD);
            rhoC_a(bad) = cu(bad);
        end
        lambda_a = TransportCoefficient.Lambda_eff(1:NL, 1:nD);
        dzCol    = cfg.Grid.DeltZ(:);
        dz_a     = dzCol(1:NL) * one_row;

        phaseA = soil.water.conductivity.calculateDallAmicoPhaseChangeArray( ...
            Theta_LL_a, Theta_II_a, TT_e, theta_s_a, theta_r_a, alpha_a, n_a, m_a, ...
            rhoC_a, lambda_a, dz_a, Delt_t);
        Theta_LL_a = phaseA.theta_l_next;
        Theta_II_a = phaseA.theta_i_next;
        Theta_UU_a = Theta_LL_a + Theta_II_a * phys.RHOI / phys.RHOL;
        SoilVariables.DallAmicoThermalTime  = phaseA.thermal_time_s;
        SoilVariables.DallAmicoThetaLStar   = phaseA.theta_l_star;
        SoilVariables.DallAmicoPhaseSource  = phaseA.source_liquid_s;
        SoilVariables.DallAmicoLimiterCode  = phaseA.limiter_code;
    end

    % =====================================================================
    % DTheta_UUh / DTheta_LLh for the lThmrlefc analytical branch.
    % The derivative follows calculateDTheta case 0: dtheta = A - B*C*D, with
    % the oven-dry head clamp heat_term <= hd. The derivative depends on the
    % retention curve and does not read Theta_UU/Theta_LL directly.
    % =====================================================================
    log_hd_hm = log(abs(hd / hm));                                  % case-0 constant
    B_d = alpha_a .* n_a .* m_a .* (Theta_m_a - Gamma_a .* theta_r_a);   % shared B term

    % DTheta_UUh: heat_term = hh+hh_frez. Zero in the saturated corner
    % (hh>=-1 & hh_frez>=0) and wherever heat_term<=hd (case-0 hd clamp).
    aUU   = abs(alpha_a .* heatTerm);
    powUU = 1 + aUU.^n_a;
    A_uu  = (-theta_r_a) ./ (abs(heatTerm) .* log_hd_hm) .* (1 - powUU.^(-m_a));
    DTheta_UUh_a = A_uu - B_d .* powUU.^(-m_a - 1) .* aUU.^(n_a - 1);
    DTheta_UUh_a((heatTerm <= hd) | ((hh_a >= -1) & (hh_frez_a >= 0))) = 0;

    % DTheta_LLh: heat_term = hh (formula), but the case-0 clamp still tests
    % hh+hh_frez. Zero where hh>=-1, where hh+hh_frez<=hd, or where Gamma==0
    % (calculateDTheta_LLh has no Gamma~=0 fallthrough — Gamma==0 -> 0).
    aLL   = abs(alpha_a .* hh_a);
    powLL = 1 + aLL.^n_a;
    A_ll  = (-theta_r_a) ./ (abs(hh_a) .* log_hd_hm) .* (1 - powLL.^(-m_a));
    DTheta_LLh_a = A_ll - B_d .* powLL.^(-m_a - 1) .* aLL.^(n_a - 1);
    DTheta_LLh_a((hh_a >= -1) | (heatTerm <= hd) | (Gamma_a == 0)) = 0;

    % =====================================================================
    % Ice ratio, viscosity-corrected KL_h, freeze-impeded KfL_h, and KfL_T.
    % =====================================================================
    % Ratio_ice = RHOI*Theta_II / (RHOL*Theta_UU); 0 where Theta_UU == 0.
    Ratio_ice_a = zeros(NL, nD);
    nz = Theta_UU_a ~= 0;
    Ratio_ice_a(nz) = phys.RHOI * Theta_II_a(nz) ./ (phys.RHOL * Theta_UU_a(nz));

    % Temperature-dependent water viscosity; reference at 20C (calculateKL_h CKT).
    % (TT_e / TT_CRIT_e were built once above, before the Dall'Amico stage.)
    % NOTE (Jizixing Zhu, 2026-07-03): the gas-constant literal 8.31441 here differs from the
    % new-chain kernel soil.water.computeEndpointHydraulics>calculateKLArray, which uses
    % phys.R = 8.314468. MU_WN/MU_W cancels at the 20C reference, so the gap induces only a benign
    % ~ppm KfL_h difference (0 at 20C; -3.4e-6 rel at 2C, growing with |T-20|). Left UNALIGNED by
    % decision (2026-07-03): phys.R is canonical, but the legacy literal is not unified to it here
    % to avoid re-baselining COMPASS.m. Verified in tests/test_runStepEbal_parity.m.
    MU_WN  = phys.MU_W0 * exp(phys.MU1 / (8.31441 * (20 + 133.3)));    % [g/(cm s)]
    MU_W_a = phys.MU_W0 .* exp(phys.MU1 ./ (8.31441 * (TT_e + 133.3)));
    MU_W_a(TT_e < -20) = 3.71e-2;
    MU_W_a(TT_e > 150) = 1.81e-3;

    % Mualem-van Genuchten K uses Se=(theta-thr)/(ths-thr).
    % SoilVariables.Se keeps the theta/POR convention for gas-phase diagnostics.
    Se_K = (Theta_LL_a - theta_r_a) ./ max(theta_s_a - theta_r_a, 1e-12);
    Se_K = max(0, min(1, Se_K));
    % Per-layer fields may be stored length NN (per node); take the first NL
    % layers used by the hydraulic-conductivity relation.
    KsCol    = SoilVariables.Ks(:);
    ImpedCol = SoilVariables.Imped(:);
    Ks_a = KsCol(1:NL) * one_row;
    KL_h_a = (MU_WN ./ MU_W_a) .* Ks_a .* sqrt(Se_K) .* (1 - (1 - Se_K.^(1 ./ m_a)).^m_a).^2;
    KL_h_a(KL_h_a <= 1e-20) = 1e-20;             % calculateKL_h numerical floor
    KL_h_a(isnan(KL_h_a)) = 0;                   % calculateKL_h NaN guard
    KL_h_a(Se_K == 0) = 0;                        % Se==0 early return (overrides floor)

    Imped_a = ImpedCol(1:NL) * one_row;
    KfL_h_a = KL_h_a .* 10.^(-Imped_a .* Ratio_ice_a);

    % KfL_T (lFreezing=1): Heaviside(TT_CRIT - (TT+T_freeze)) * L_f*1e4/(g*T_freeze)
    KfL_T_a = helpers.heaviside1(TT_CRIT_e - (TT_e + phys.T_freeze)) * L_f * 1e4 / (phys.g * phys.T_freeze);

    % ---- Se: effective saturation (calculateSe, lSFCC=1 / lThmrlefc=1 path) ----
    % seCalc=Theta_LL/POR; force 1 in the saturated corner, 0 in the dry corner,
    % then clamp [0,1]. Reads the post-phase-change Theta_LL.
    Se_a   = Theta_LL_a ./ POR_a;
    sat_se = (hh_a >= -1) & (hh_frez_a >= 0);
    dry_se = ~sat_se & (heatTerm <= hd);
    Se_a(sat_se) = 1;
    Se_a(dry_se) = 0;
    Se_a = max(0, min(1, Se_a));

    % Gamma_hh per node: nodes 1..NL use the left endpoint column; the final
    % boundary node uses the last layer's right endpoint.
    gamma_node = [Gamma_a(1:NL, 1); Gamma_a(NL, 2:nD).'];   % nodes 1..NL+nD-1

    % =====================================================================
    % Write back hydraulic state and auxiliary fields.
    % =====================================================================
    SoilVariables.Theta_LL   = Theta_LL_a;
    SoilVariables.Theta_UU   = Theta_UU_a;
    SoilVariables.Theta_II   = Theta_II_a;
    SoilVariables.DTheta_UUh = DTheta_UUh_a;
    SoilVariables.DTheta_LLh = DTheta_LLh_a;
    SoilVariables.KL_h       = KL_h_a;
    SoilVariables.KfL_h      = KfL_h_a;
    SoilVariables.KfL_T      = KfL_T_a;
    SoilVariables.Ratio_ice  = Ratio_ice_a;
    SoilVariables.Se         = Se_a;
    SoilVariables.Gamma_hh(1:NL + nD - 1) = gamma_node;
end

% =========================================================================
% Local helper for thermally corrected matric head.
% =========================================================================
function hh = fixHeatArray(hh, cfg)
% Apply the same lower-head clamp convention as fixHeat for non-strict SFCC.
% The strict SFCC path leaves the thermally corrected head unchanged.
    hh = hh(:);
    if cfg.lSWCC == 1
        if ~(cfg.lFreezing == 1 && cfg.lSFCC == 1)
            hh(hh < -1e7) = -1e7;
        end
    end
end

function tf = getLogicalOption(cfg, name, defaultValue)
tf = defaultValue;
if isstruct(cfg)
    if isfield(cfg, name)
        tf = logical(cfg.(name));
    end
elseif isobject(cfg) && isprop(cfg, name)
    tf = logical(cfg.(name));
end
end

function tf = hasThermalTimeInputs(TransportCoefficient)
tf = false;
if isempty(TransportCoefficient) || ~isfield(TransportCoefficient, 'Lambda_eff')
    return;
end
if ~isfield(TransportCoefficient, 'EHCAP') && ~isfield(TransportCoefficient, 'c_unsat')
    return;
end
lambda = TransportCoefficient.Lambda_eff;
if isempty(lambda) || ~any(isfinite(lambda(:)) & lambda(:) > 0)
    return;
end
if isfield(TransportCoefficient, 'EHCAP')
    rhoC = TransportCoefficient.EHCAP;
else
    rhoC = [];
end
if isempty(rhoC) || ~any(isfinite(rhoC(:)) & rhoC(:) > 0)
    rhoC = TransportCoefficient.c_unsat;
end
tf = ~isempty(rhoC) && any(isfinite(rhoC(:)) & rhoC(:) > 0);
end
