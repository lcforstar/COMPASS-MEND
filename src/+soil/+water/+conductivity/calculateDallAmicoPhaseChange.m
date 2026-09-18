function phase = calculateDallAmicoPhaseChange(theta_l, theta_i, T_C, porosity, theta_r, ...
    alpha_cm, n_vg, m_vg, rhoC_native, lambda_native, dz_cm, dt_s)
% calculateDallAmicoPhaseChange - Dall'Amico/ClimaLand-style phase change source.
%
% Inputs use COMPASS native units:
%   theta_l       liquid water volume fraction [-]
%   theta_i       ice volume fraction [-]
%   T_C           soil temperature [degC]
%   porosity      effective pore volume, normally theta_s or POR [-]
%   theta_r       residual water content [-]
%   alpha_cm      van Genuchten alpha [1/cm]
%   n_vg, m_vg    van Genuchten shape parameters [-]
%   rhoC_native   volumetric heat capacity EHCAP/c_unsat [J/(cm3 K)]
%                 This is numerically equal to [MJ/(m3 K)].
%   lambda_native effective thermal conductivity Lambda_eff [W/(cm K)]
%   dz_cm         layer thickness [cm]
%   dt_s          optional update interval [s]
%
% Output:
%   phase.thermal_time_s      3*rhoC*dz^2/lambda [s]
%   phase.theta_l_star        liquid equilibrium from freezing-point depression [-]
%   phase.source_liquid_s     (theta_l - theta_l_star)/thermal_time [1/s]
%   phase.dtheta_l, dtheta_i  bounded finite-step increments over dt_s [-]
%   phase.rhoC_SI             rhoC_native converted to [J/(m3 K)]
%   phase.lambda_SI           lambda_native converted to [W/(m K)]
%
% This is intentionally independent from the current STEMMUS freezing path.
% It follows the ClimaLand/Dall'Amico structure:
%   thermal_time = 3 * rhoC * dz^2 / lambda
%   theta_total  = min(rho_i/rho_l * theta_i + theta_l, porosity)
%   psi_T        = L_f/g * log(T/Tf_depressed) for T < Tf_depressed
%   theta_l_star = retention_inverse(psi_w0 + psi_T)
%   source       = (theta_l - theta_l_star) / thermal_time

if nargin < 12 || isempty(dt_s)
    dt_s = 0;
end
if nargin < 8 || isempty(m_vg)
    m_vg = 1 - 1 / n_vg;
end

phys = equations.physicalConstants();

theta_l = max(0, theta_l);
theta_i = max(0, theta_i);
porosity = max(theta_r + eps, porosity);
theta_r = max(0, min(theta_r, porosity - eps));
alpha_cm = max(alpha_cm, realmin);
n_vg = max(n_vg, 1 + eps);
m_vg = max(m_vg, eps);
rhoC_native = max(rhoC_native, realmin);
lambda_native = max(lambda_native, realmin);
dz_cm = max(dz_cm, realmin);
dt_s = max(dt_s, 0);

rho_l = phys.RHOL;
rho_i = phys.RHOI;
rho_i_over_l = rho_i / rho_l;
rho_l_over_i = rho_l / rho_i;

rhoC_SI = rhoC_native * 1e6;
lambda_SI = lambda_native * 100;
dz_m = dz_cm / 100;
thermal_time_s = 3 * rhoC_SI * dz_m^2 / lambda_SI;

theta_total = min(rho_i_over_l * theta_i + theta_l, porosity);
Se_total = (theta_total - theta_r) / (porosity - theta_r);
Se_total = clamp(Se_total, eps, 1);

psi_w0_cm = vgHeadFromSe(Se_total, alpha_cm, n_vg, m_vg);
psi_w0_m = psi_w0_cm / 100;

T_K = max(T_C + phys.T_freeze, eps);
Tf_depressed_K = phys.T_freeze * exp(phys.g_SI * psi_w0_m / phys.lambda_fus);
if T_K < Tf_depressed_K
    psi_T_m = phys.lambda_fus / phys.g_SI * log(T_K / Tf_depressed_K);
else
    psi_T_m = 0;
end

psi_star_cm = (psi_w0_m + psi_T_m) * 100;
Se_star = vgSeFromHead(psi_star_cm, alpha_cm, n_vg, m_vg);
theta_l_star = theta_r + Se_star * (porosity - theta_r);
theta_l_star = clamp(theta_l_star, theta_r, porosity);

source_liquid_s = (theta_l - theta_l_star) / thermal_time_s;
[dtheta_l, dtheta_i, limiter] = boundedPhaseStep(source_liquid_s, theta_l, theta_i, ...
    theta_r, porosity, rho_l_over_i, rho_i_over_l, dt_s);

phase = struct();
phase.theta_total = theta_total;
phase.Se_total = Se_total;
phase.psi_w0_cm = psi_w0_cm;
phase.psi_w0_m = psi_w0_m;
phase.T_K = T_K;
phase.Tf_depressed_K = Tf_depressed_K;
phase.psi_T_m = psi_T_m;
phase.psi_star_cm = psi_star_cm;
phase.theta_l_star = theta_l_star;
phase.source_liquid_s = source_liquid_s;
phase.thermal_time_s = thermal_time_s;
phase.dtheta_l = dtheta_l;
phase.dtheta_i = dtheta_i;
phase.theta_l_next = theta_l + dtheta_l;
phase.theta_i_next = theta_i + dtheta_i;
phase.limiter = limiter;
phase.rhoC_SI = rhoC_SI;
phase.lambda_SI = lambda_SI;
phase.dz_m = dz_m;
phase.units = 'native inputs: EHCAP J/cm3/K (=MJ/m3/K), Lambda_eff W/cm/K, dz cm';
end

function h_cm = vgHeadFromSe(Se, alpha_cm, n_vg, m_vg)
if Se >= 1
    h_cm = 0;
else
    h_cm = -((Se^(-1 / m_vg) - 1)^(1 / n_vg)) / alpha_cm;
end
end

function Se = vgSeFromHead(h_cm, alpha_cm, n_vg, m_vg)
Se = (1 + abs(alpha_cm * h_cm)^n_vg)^(-m_vg);
Se = clamp(Se, 0, 1);
end

function [dtheta_l, dtheta_i, limiter] = boundedPhaseStep(source_liquid_s, theta_l, theta_i, ...
    theta_r, porosity, rho_l_over_i, rho_i_over_l, dt_s)
limiter = 'none';
dtheta_l = -source_liquid_s * dt_s;
dtheta_i = rho_l_over_i * source_liquid_s * dt_s;

if dt_s <= 0 || source_liquid_s == 0
    dtheta_l = 0;
    dtheta_i = 0;
    return;
end

if source_liquid_s > 0
    liquid_loss = -dtheta_l;
    max_loss_by_liquid = max(0, theta_l - theta_r);
    expansion = rho_l_over_i - 1;
    if expansion > 0
        max_loss_by_pore = max(0, porosity - theta_l - theta_i) / expansion;
    else
        max_loss_by_pore = Inf;
    end
    bounded_loss = min([liquid_loss, max_loss_by_liquid, max_loss_by_pore]);
    if bounded_loss < liquid_loss
        limiter = 'freezing-limited';
    end
    dtheta_l = -bounded_loss;
    dtheta_i = rho_l_over_i * bounded_loss;
else
    ice_loss = -dtheta_i;
    bounded_ice_loss = min(ice_loss, theta_i);
    if bounded_ice_loss < ice_loss
        limiter = 'melting-limited';
    end
    dtheta_i = -bounded_ice_loss;
    dtheta_l = rho_i_over_l * bounded_ice_loss;
end
end

function value = clamp(value, lo, hi)
value = min(max(value, lo), hi);
end
