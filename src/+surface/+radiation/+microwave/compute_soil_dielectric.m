function aux = compute_soil_dielectric(state, params)
%COMPUTE_SOIL_DIELECTRIC  Effective soil dielectric and rough-surface reflectivity for L-band.
%
%   Returns:
%     aux.smc_eff   [m³/m³]   EM-weighted effective soil moisture, iteratively
%                              consistent with eps_soil (Wilheit 1978)
%     aux.delta_p   [cm]      Wilheit penetration depth, λ / (2π √ε″)
%     aux.eps_soil  [-]       Complex dielectric constant of soil
%     aux.r_V       [-]       Rough-surface reflectivity, V polarization
%     aux.r_H       [-]       Rough-surface reflectivity, H polarization
%     aux.eS_V      [-]       Rough-surface emissivity = 1 - r_V
%     aux.eS_H      [-]       Rough-surface emissivity = 1 - r_H
%
%   Inputs:
%     state.theta_prof  [m³/m³, NN×1]   Layer liquid water content (top→bottom)
%     state.T_prof_K    [K, NN×1]       Layer soil temperature
%     state.depth_cm    [cm, NN×1]      Node depth, top = 0 cm
%     state.sand_prof   [-,  NN×1]      Sand fraction per node (decimal, 0–1)
%     state.clay_prof   [-,  NN×1]      Clay fraction per node (decimal, 0–1)
%     state.bd_prof     [g/cm³, NN×1]   Bulk density per node
%       NOTE: caller must convert from SoilProperties (1=bottom) to top-down
%             ordering matching depth_cm(1)=0 surface. See load_compass_state.
%     params.calib.sig, .cl             RMS roughness height, correlation length (cm)
%     params.calib.Fre, .Inc            Frequency (GHz), incidence angle (deg)
%     params.calib.iDC_warm, .iDC_cold  Dielectric model indices for T>0 / T<0
%     params.calib.iPR                  Roughness model index (1=Wigneron 2011)
%
%   Method (four optimisations folded in here):
%     ① Per-layer dielectric — every node uses its own (θ, T, sand, clay, bd)
%        so texture stratification + freeze-state heterogeneity propagate
%        into the dielectric profile.
%     ② Per-layer freeze-aware model switch — Mironov (iDC_warm) for T ≥ 0 °C
%        nodes, Zhang (iDC_cold) for T < 0 °C; the latter explicitly tracks
%        ice fraction vi via unfrozen-water curve (Zhang & Sun 2011 IGARSS).
%     ③ Wilheit-weighted effective ε — δ_p = λ / (2π √Im(ε)). Iterate twice
%        to make ε_eff and δ_p mutually consistent (Wilheit 1978 RSE).
%     ④ Rough-surface emissivity via PassiveRoughness (Wigneron 2011 default).
%
%   References:
%     Wilheit, T.T. 1978, RSE 7(1), 67–86.        — penetration depth
%     Mironov et al. 2009, IEEE-TGRS 47(7).        — dielectric, T > 0 °C
%     Zhang & Sun 2011, IGARSS.                   — dielectric with ice (iDC=4)
%     Wigneron et al. 2011, IEEE-TGRS 49(3).      — rough surface (iPR=1)
%
%   Notes:
%     This function is part of the standalone +surface/+radiation/+microwave subpackage
%     and does NOT modify any COMPASS state. It reads precomputed COMPASS
%     output (state struct) and returns an aux struct that the caller
%     (microwave_forward) routes into the τ-ω synthesis. No globals.

import surface.radiation.microwave.compute_dielectric_constant
import surface.radiation.microwave.compute_fresnel_reflectivity
import surface.radiation.microwave.compute_rough_surface_emissivity

% ─────────────────────────────────────────────────────────────────────────
% Profile vectors (top-down: index 1 = 0 cm surface, end = bottom)
% ─────────────────────────────────────────────────────────────────────────
depth_cm   = state.depth_cm(:);
theta_prof = state.theta_prof(:);
T_prof_K   = state.T_prof_K(:);
sand_prof  = state.sand_prof(:);
clay_prof  = state.clay_prof(:);
bd_prof    = state.bd_prof(:);

dz_cm = diff([0; depth_cm]);                                 % [cm] node interval
dz_cm = max(dz_cm, 1e-3);                                    % guard zero

% ─────────────────────────────────────────────────────────────────────────
% Optimisation ②: freeze-aware dielectric model selection (per layer)
%   Mironov is unfrozen-only; Zhang explicitly partitions unfrozen liquid
%   vs ice via the freezing curve.  Each node chooses its model from its
%   own temperature, so the warm/cold interface near 0 °C is handled
%   physically rather than via a top-node binary switch.
% ─────────────────────────────────────────────────────────────────────────
NN = numel(depth_cm);
iDC_prof  = zeros(NN, 1);
for i = 1:NN
    if T_prof_K(i) < 273.15
        iDC_prof(i) = params.calib.iDC_cold;                 % e.g. 4 Zhang
    else
        iDC_prof(i) = params.calib.iDC_warm;                 % e.g. 1 Mironov
    end
end

% ─────────────────────────────────────────────────────────────────────────
% Per-layer complex dielectric ε(z)
%   Each node uses its own (θ, T, sand, clay, bd) so the dielectric profile
%   reflects all texture stratification and freeze-state heterogeneity.
% ─────────────────────────────────────────────────────────────────────────
eps_prof = complex(zeros(NN, 1));
for i = 1:NN
    theta_i = max(min(theta_prof(i), 0.55), 0.02);
    T_i_C   = T_prof_K(i) - 273.15;
    eps_prof(i) = compute_dielectric_constant(params.calib.Fre, ...
                                     sand_prof(i), clay_prof(i), bd_prof(i), ...
                                     T_i_C, theta_i, iDC_prof(i));
end

% ─────────────────────────────────────────────────────────────────────────
% Optimisation ③: Wilheit 1978 EM-weighted effective dielectric
%   δ_p(z=0) = λ / (2π √Im(ε_top)).  Iterate two passes:
%     1. seed δ_p from top-node ε
%     2. weight nodes by  w_i ∝ exp(-2 z_i / δ_p) · dz_i
%     3. re-derive δ_p from ε_eff
%   This produces ε_eff that respects both the texture profile and the
%   sensing depth (~5 cm in moist soil, ~10 cm dry).  Two passes converge
%   within < 0.2 cm.
% ─────────────────────────────────────────────────────────────────────────
c_mps     = 2.998e8;                                         % [m/s]
lambda_cm = (c_mps / (params.calib.Fre * 1e9)) * 100;        % [cm] wavelength

eps_top  = eps_prof(1);
eps_imag = max(imag(eps_top), 1e-3);
delta_p_cm = lambda_cm / (2 * pi * sqrt(eps_imag));

for it = 1:2
    w  = exp(-2 * depth_cm / delta_p_cm) .* dz_cm;
    w  = w / max(sum(w), 1e-9);
    eps_soil = sum(w .* eps_prof);                            % complex-weighted
    eps_imag = max(imag(eps_soil), 1e-3);
    delta_p_cm = lambda_cm / (2 * pi * sqrt(eps_imag));
end

% Diagnostic: EM-weighted effective soil moisture (for cross-check)
smc_eff = sum(w .* theta_prof);

% ─────────────────────────────────────────────────────────────────────────
% Fresnel reflectivity at the flat planar boundary (air → soil)
% ─────────────────────────────────────────────────────────────────────────
eps_air = complex(1.0, 0.0);
[~, ~, gV_flat, gH_flat] = compute_fresnel_reflectivity(params.calib.Inc, eps_air, eps_soil);

% ─────────────────────────────────────────────────────────────────────────
% Rough-surface emissivity (Wigneron 2011 by default; iPR selectable)
% Returns rough-surface emissivities (eS_V, eS_H); reflectivity = 1 - emiss.
% ─────────────────────────────────────────────────────────────────────────
[eS_V, eS_H] = compute_rough_surface_emissivity(params.calib.Fre, params.calib.Inc, ...
                                params.calib.sig, params.calib.cl, ...
                                gV_flat, gH_flat, params.calib.iPR);

% ─────────────────────────────────────────────────────────────────────────
% Pack auxiliary output
% ─────────────────────────────────────────────────────────────────────────
aux.smc_eff   = smc_eff;          % [m³/m³] EM-weighted effective SMC
aux.delta_p   = delta_p_cm;       % [cm] Wilheit penetration depth
aux.eps_soil  = eps_soil;         % [-] EM-weighted complex dielectric
aux.eS_V      = eS_V;
aux.eS_H      = eS_H;
aux.r_V       = 1 - eS_V;
aux.r_H       = 1 - eS_H;
aux.iDC_top   = iDC_prof(1);      % which DC model the surface node used
aux.iDC_bot   = iDC_prof(end);
aux.n_frozen  = sum(iDC_prof == params.calib.iDC_cold);  % count of frozen nodes

end
