function phase = calculateDallAmicoPhaseChangeArray(theta_l, theta_i, T_C, porosity, theta_r, ...
    alpha_cm, n_vg, m_vg, rhoC_native, lambda_native, dz_cm, dt_s)
%CALCULATEDALLAMICOPHASECHANGEARRAY  Dall'Amico thermal-time phase-change source.
%
% Applies the soil freezing/thawing relaxation to each layer endpoint. The
% limiter is returned as a numeric phase-state code:
%   0 = no bound active, 1 = freezing limited, 2 = melting limited.
%
% Inputs (each NL x nD unless noted; COMPASS native units):
%   theta_l        [-]        liquid water volume fraction
%   theta_i        [-]        ice volume fraction
%   T_C            [degC]     soil (endpoint) temperature
%   porosity       [-]        effective pore volume (theta_s / POR)
%   theta_r        [-]        residual water content
%   alpha_cm       [1/cm]     van Genuchten alpha
%   n_vg, m_vg     [-]        van Genuchten shape parameters
%   rhoC_native    [J/(cm3 K)] volumetric heat capacity (= MJ/(m3 K))
%   lambda_native  [W/(cm K)] effective thermal conductivity
%   dz_cm          [cm]       layer thickness
%   dt_s           [s]        update interval
%
% Outputs (struct, each NL x nD):
%   phase.theta_l_next     [-]    theta_l + bounded dtheta_l
%   phase.theta_i_next     [-]    theta_i + bounded dtheta_i
%   phase.thermal_time_s   [s]    3*rhoC*dz^2/lambda
%   phase.theta_l_star     [-]    liquid equilibrium from freezing-point depression
%   phase.source_liquid_s  [1/s]  (theta_l - theta_l_star)/thermal_time
%   phase.limiter_code     [-]    0 none / 1 freezing-limited / 2 melting-limited
%
% Method:
%   thermal_time = 3 * rhoC * dz^2 / lambda
%   theta_total  = min(rho_i/rho_l * theta_i + theta_l, porosity)
%   psi_T        = L_f/g * log(T/Tf_depressed) for T < Tf_depressed, else 0
%   theta_l_star = retention_inverse(psi_w0 + psi_T)
%   source       = (theta_l - theta_l_star) / thermal_time
% then a mass/pore-bounded finite step over dt_s (boundedPhaseStep).
%
% Reference note:
%   The formulation follows the existing COMPASS Dall'Amico/ClimaLand-style
%   thermal-time phase-change structure. No new literature citation is claimed
%   here beyond the source formulation already used by the model.
%
% Author: Zhu Jizixing, 2026-06-15.

    if nargin < 12 || isempty(dt_s)
        dt_s = 0;
    end
    if nargin < 8 || isempty(m_vg)
        m_vg = 1 - 1 ./ n_vg;
    end

    phys = equations.physicalConstants();

    % --- Input sanitisation ---
    theta_l  = max(0, theta_l);
    theta_i  = max(0, theta_i);
    porosity = max(theta_r + eps, porosity);
    theta_r  = max(0, min(theta_r, porosity - eps));
    alpha_cm = max(alpha_cm, realmin);
    n_vg     = max(n_vg, 1 + eps);
    m_vg     = max(m_vg, eps);
    rhoC_native   = max(rhoC_native, realmin);
    lambda_native = max(lambda_native, realmin);
    dz_cm    = max(dz_cm, realmin);
    dt_s     = max(dt_s, 0);

    rho_l = phys.RHOL;
    rho_i = phys.RHOI;
    rho_i_over_l = rho_i / rho_l;
    rho_l_over_i = rho_l / rho_i;

    rhoC_SI   = rhoC_native * 1e6;
    lambda_SI = lambda_native * 100;
    dz_m      = dz_cm / 100;
    thermal_time_s = 3 * rhoC_SI .* dz_m.^2 ./ lambda_SI;

    % --- Equilibrium liquid content from freezing-point depression ---
    theta_total = min(rho_i_over_l * theta_i + theta_l, porosity);
    Se_total = (theta_total - theta_r) ./ (porosity - theta_r);
    Se_total = min(max(Se_total, eps), 1);                 % clamp(eps,1)

    % vgHeadFromSe: h = -((Se^(-1/m)-1)^(1/n))/alpha; 0 where Se>=1.
    psi_w0_cm = -((Se_total.^(-1 ./ m_vg) - 1).^(1 ./ n_vg)) ./ alpha_cm;
    psi_w0_cm(Se_total >= 1) = 0;
    psi_w0_m = psi_w0_cm / 100;

    T_K = max(T_C + phys.T_freeze, eps);
    Tf_depressed_K = phys.T_freeze * exp(phys.g_SI * psi_w0_m / phys.lambda_fus);
    psi_T_m = phys.lambda_fus / phys.g_SI .* log(T_K ./ Tf_depressed_K);
    psi_T_m(T_K >= Tf_depressed_K) = 0;                    % only depress below Tf

    psi_star_cm = (psi_w0_m + psi_T_m) * 100;
    % vgSeFromHead: Se = (1+|alpha*h|^n)^(-m), clamp[0,1].
    Se_star = (1 + abs(alpha_cm .* psi_star_cm).^n_vg).^(-m_vg);
    Se_star = min(max(Se_star, 0), 1);
    theta_l_star = theta_r + Se_star .* (porosity - theta_r);
    theta_l_star = min(max(theta_l_star, theta_r), porosity);   % clamp(theta_r,porosity)

    source_liquid_s = (theta_l - theta_l_star) ./ thermal_time_s;

    % --- Bounded finite step over dt_s ---
    sz = size(source_liquid_s);
    dtheta_l = -source_liquid_s * dt_s;
    dtheta_i = rho_l_over_i * source_liquid_s * dt_s;
    limiter_code = zeros(sz);

    inactive = (dt_s <= 0) | (source_liquid_s == 0);
    freezing = ~inactive & (source_liquid_s > 0);
    melting  = ~inactive & (source_liquid_s < 0);

    dtheta_l(inactive) = 0;
    dtheta_i(inactive) = 0;

    % Freezing (source>0): liquid -> ice, bounded by available liquid and pore.
    expansion = rho_l_over_i - 1;                          % > 0 (water denser than ice)
    liquid_loss = source_liquid_s * dt_s;
    max_loss_by_liquid = max(0, theta_l - theta_r);
    max_loss_by_pore   = max(0, porosity - theta_l - theta_i) / expansion;
    bounded_loss = min(min(liquid_loss, max_loss_by_liquid), max_loss_by_pore);
    fr_limited = freezing & (bounded_loss < liquid_loss);
    dtheta_l(freezing) = -bounded_loss(freezing);
    dtheta_i(freezing) = rho_l_over_i * bounded_loss(freezing);
    limiter_code(fr_limited) = 1;

    % Melting (source<0): ice -> liquid, bounded by available ice.
    ice_loss = -rho_l_over_i * source_liquid_s * dt_s;
    bounded_ice_loss = min(ice_loss, theta_i);
    mt_limited = melting & (bounded_ice_loss < ice_loss);
    dtheta_i(melting) = -bounded_ice_loss(melting);
    dtheta_l(melting) = rho_i_over_l * bounded_ice_loss(melting);
    limiter_code(mt_limited) = 2;

    phase = struct();
    phase.theta_l_next    = theta_l + dtheta_l;
    phase.theta_i_next    = theta_i + dtheta_i;
    phase.thermal_time_s  = thermal_time_s;
    phase.theta_l_star    = theta_l_star;
    phase.source_liquid_s = source_liquid_s;
    phase.limiter_code    = limiter_code;
end
