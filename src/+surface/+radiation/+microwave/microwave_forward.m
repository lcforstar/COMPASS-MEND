function out = microwave_forward(state, params)
%MICROWAVE_FORWARD  Single time-step L-band τ-ω forward simulator.
%
%   Returns:
%     out.VOD       [-]    Total vegetation optical depth (nadir)
%     out.Tb_V      [K]    Brightness temperature, V polarization
%     out.Tb_H      [K]    Brightness temperature, H polarization
%     out.gamma_V   [-]    Canopy transmissivity V (= exp(-VOD/cos θ))
%     out.gamma_H   [-]    Canopy transmissivity H
%     out.eS_V/eS_H [-]    Rough-surface emissivity
%     out.eps_soil  [-]    Effective complex dielectric (single-layer mode) or
%                          top-node dielectric (Wilheit mode, diagnostic only)
%     out.smc_eff   [m³/m³]  EM-weighted effective surface moisture
%     out.delta_p   [cm]   Wilheit penetration depth
%     out.vod_comp  struct VOD component breakdown (leaf/stem/int)
%     out.Tg_eff    [K]    Effective ground temperature (top T or Wilheit-weighted)
%
%   Inputs:
%     state.{theta_prof, T_prof_K, depth_cm, sand_prof, clay_prof, bd_prof}
%        — full soil profile (top-down) for this time step
%     state.{LAI, h2ocan, psi_leaf_MPa, Tc_K, Ts_top_K}
%        — canopy state for this time step
%     params.calib.{Fre, Inc, omega, sig, cl,
%                   iDC_warm, iDC_cold, iPR,
%                   iRough,                                 (default 1)
%                   iLayer,                                 (default 1)
%                   b_leaf, b_stem, b_int,
%                   alpha_psi, psi_ref, f_frozen}
%     params.site.{Cw, AGB_static}
%
%   Mode switches:
%     iRough  1 = Wigneron 2011 empirical (fast, default)
%             2 = I2EM Fung 2002 (Improved IEM, slow ~50 ms/call)
%             Both invoked AFTER smooth-surface reflectivity is obtained.
%
%     iLayer  1 = single-layer EM-weighted effective ε (default, fast)
%             2 = Wilheit 1978 multilayer recursive Fresnel (slow but
%                 honours dielectric/temperature stratification rigorously)
%
%   Pipeline (τ-ω two-source model, Wigneron 2007 L-MEB Eq 1):
%       Tb_p = (1-ω_p)(1-γ_p)(1 + γ_p·r_p)·Tc + (1-r_p)·γ_p·Tg
%
%     1. Soil branch:
%        — single-layer:  per-layer ε → Wilheit weighted ε_eff → Fresnel
%        — multilayer:    per-layer ε → Wilheit_multilayer recursion (r,Tg)
%        Then iRough applies rough-surface correction on top.
%     2. Vegetation branch: component VOD with ψ-modulated b_leaf
%     3. τ-ω synthesis: combine using state.Tc_K and Tg_eff
%
%   References:
%     Wigneron et al. 2007, Remote Sens. Environ. 107, 639–655 — L-MEB
%     Mo et al. 1982, JGR-Oceans 87, 11229–11237 — τ-ω formulation
%     Wilheit 1978, IEEE TGRS GE-16(2), 138–143 — stratified-ε model
%     Fung, Liu, Chen, Tsang 2002, JEMWA 16, 689 — I2EM

import surface.radiation.microwave.compute_soil_dielectric
import surface.radiation.microwave.compute_VOD
import surface.radiation.microwave.compute_atmospheric_attenuation
import surface.radiation.microwave.compute_dielectric_constant
import surface.radiation.microwave.compute_rough_surface_emissivity
import surface.radiation.microwave.compute_i2em_emissivity
import surface.radiation.microwave.compute_wilheit_multilayer

% Defaults for new switches (backward compatible)
if ~isfield(params.calib, 'iRough'), params.calib.iRough = 1; end
if ~isfield(params.calib, 'iLayer'), params.calib.iLayer = 1; end
if ~isfield(params.calib, 'iAtm'),   params.calib.iAtm   = 0; end

% ─── 1. Soil branch ───────────────────────────────────────────────────────
switch params.calib.iLayer
    case 1
        % Single-layer EM-weighted ε (existing path).
        % compute_soil_dielectric internally applies PassiveRoughness, but
        % we need the smooth-Fresnel reflectivities here so the iRough branch
        % below can choose Wigneron vs I2EM uniformly. Re-derive from eps_soil.
        aux_soil = compute_soil_dielectric(state, params);
        % Tg: prefer SCOPE-computed skin temperature (Ts_top_K = Tsave =
        % LAI-weighted sun/shade soil surface T) over the top soil node
        % T_prof(1). Tsave is the physically correct microwave-emitting
        % surface temperature; T_prof(1) is the layer-mean of 0-1 cm bulk
        % soil. Fall back to T_prof(1) only if Ts_top_K is missing.
        if isfield(state, 'Ts_top_K') && ~isnan(state.Ts_top_K)
            Tg_eff = state.Ts_top_K;
        else
            Tg_eff = state.T_prof_K(1);
        end
        eps_air  = complex(1, 0);
        [~, ~, gV_smooth, gH_smooth] = ...
            surface.radiation.microwave.compute_fresnel_reflectivity(params.calib.Inc, ...
                                                       eps_air, aux_soil.eps_soil);

    case 2
        % Wilheit multilayer recursive Fresnel
        depth_cm   = state.depth_cm(:);
        theta_prof = state.theta_prof(:);
        T_prof_K   = state.T_prof_K(:);
        sand_prof  = state.sand_prof(:);
        clay_prof  = state.clay_prof(:);
        bd_prof    = state.bd_prof(:);
        NN = numel(depth_cm);

        % Per-node complex dielectric (same logic as compute_soil_dielectric)
        eps_prof = complex(zeros(NN, 1));
        for i = 1:NN
            if T_prof_K(i) < 273.15
                iDC_i = params.calib.iDC_cold;
            else
                iDC_i = params.calib.iDC_warm;
            end
            theta_i = max(min(theta_prof(i), 0.55), 0.02);
            T_i_C   = T_prof_K(i) - 273.15;
            eps_prof(i) = compute_dielectric_constant(params.calib.Fre, ...
                                             sand_prof(i), clay_prof(i), bd_prof(i), ...
                                             T_i_C, theta_i, iDC_i);
        end

        % Build Wilheit n/dz/T profile: prepend air layer (top), append
        % semi-infinite bottom replica.
        n_prof_w   = [complex(1, 0); sqrt(eps_prof); sqrt(eps_prof(end))];
        dz_node_cm = diff([0; depth_cm]);
        dz_prof_w  = [1e9; max(dz_node_cm, 1e-3); 1e9];
        T_prof_w   = [T_prof_K(1); T_prof_K; T_prof_K(end)];

        [gV_smooth, TbV_layer] = compute_wilheit_multilayer(params.calib.Fre, params.calib.Inc, ...
                                                     n_prof_w, dz_prof_w, T_prof_w, 'V');
        [gH_smooth, TbH_layer] = compute_wilheit_multilayer(params.calib.Fre, params.calib.Inc, ...
                                                     n_prof_w, dz_prof_w, T_prof_w, 'H');

        % Effective Tg from Wilheit weighted emission (avg V,H)
        Tg_eff = 0.5 * (TbV_layer + TbH_layer) / max(1 - 0.5*(gV_smooth+gH_smooth), 1e-6);
        % (TbV_layer already includes (1-R) factor; restore plain Tg)

        % Diagnostics
        c_mps = 2.998e8;
        lambda_cm = (c_mps / (params.calib.Fre * 1e9)) * 100;
        delta_p_cm = lambda_cm / (2*pi * sqrt(max(imag(eps_prof(1)), 1e-3)));
        w_em = exp(-2 * depth_cm / delta_p_cm) .* max(dz_node_cm, 1e-3);
        w_em = w_em / max(sum(w_em), 1e-9);
        aux_soil.eps_soil  = eps_prof(1);
        aux_soil.smc_eff   = sum(w_em .* theta_prof);
        aux_soil.delta_p   = delta_p_cm;
        aux_soil.iDC_top   = (T_prof_K(1) < 273.15) * params.calib.iDC_cold ...
                           + (T_prof_K(1) >= 273.15) * params.calib.iDC_warm;
        aux_soil.n_frozen  = sum(T_prof_K < 273.15);

    otherwise
        error('microwave_forward: unsupported iLayer = %d', params.calib.iLayer);
end

% ─── 1b. Apply rough-surface emissivity model ─────────────────────────────
switch params.calib.iRough
    case 1
        % Wigneron 2011 empirical (existing)
        [eS_V, eS_H] = compute_rough_surface_emissivity(params.calib.Fre, params.calib.Inc, ...
                                        params.calib.sig, params.calib.cl, ...
                                        gV_smooth, gH_smooth, params.calib.iPR);

    case 2
        % I2EM Fung 2002 — uses eps_soil directly, not gV/gH
        % I2EM expects sig, L in meters; we have them in cm (config convention)
        sig_m = params.calib.sig * 0.01;
        L_m   = params.calib.cl  * 0.01;
        sp    = 1;                                   % exponential spectrum (default)
        [eS_V, eS_H] = compute_i2em_emissivity(params.calib.Fre, sig_m, L_m, ...
                                              params.calib.Inc, aux_soil.eps_soil, sp);
        % Clip to physical range
        eS_V = max(0.1, min(eS_V, 1.0));
        eS_H = max(0.1, min(eS_H, 1.0));

    otherwise
        error('microwave_forward: unsupported iRough = %d', params.calib.iRough);
end
r_V = 1 - eS_V;
r_H = 1 - eS_H;

% Tg used in τ-ω synthesis.
%   iLayer=1: Tg_eff was set above to Ts_top_K (Tsave from SCOPE)
%   iLayer=2: Tg_eff is the Wilheit weighted multilayer-emission T
% Both already represent the correct microwave emitting temperature.
Tg = Tg_eff;

% ─── 2. Vegetation branch ─────────────────────────────────────────────────
vod = compute_VOD(state, params);

% ─── 3. τ-ω synthesis (Wigneron 2007 Eq 1, polarization-split) ────────────
% Wigneron L-MEB polarization-corrected VOD (CMEM vegwign formula):
%   τ_p = τ_N · (cos²θ + tt_p · sin²θ) / cosθ
% Falls back to symmetric VOD/cosθ when tt_H = tt_V = 1 (Jackson default).
cos_theta = cosd(params.calib.Inc);
sin2 = sind(params.calib.Inc)^2;
cos2 = cos_theta^2;
tau_H = vod.tau_N * (cos2 + vod.tt_H * sin2) / cos_theta;
tau_V = vod.tau_N * (cos2 + vod.tt_V * sin2) / cos_theta;
gamma_V = exp(-tau_V);
gamma_H = exp(-tau_H);

omega = params.calib.omega;

Tc = state.Tc_K;
Tb_V_toc = (1-omega) * (1-gamma_V) * (1 + gamma_V*r_V) * Tc ...
         + (1 - r_V) * gamma_V * Tg;
Tb_H_toc = (1-omega) * (1-gamma_H) * (1 + gamma_H*r_H) * Tc ...
         + (1 - r_H) * gamma_H * Tg;

% ─── 4. Atmospheric correction (Pellarin 2003 L-band, CMEM ATMPELLARIN) ──
% TOC (top of canopy) → TOA (top of atmosphere):
%   Tb_TOA = τ_atm · Tb_TOC + Tb_up + τ_atm · r_TOC · γ² · (Tb_dn + τ_atm·T_cosmic)
% where r_TOC = (1-eS)·γ² is the effective surface reflectivity seen from TOC.
% iAtm = 0 (default) skips this step → Tb_TOA = Tb_TOC (back-compat).
state_atm.Tair_K = state.Tair_K;
atm = compute_atmospheric_attenuation(params, state_atm);

if params.calib.iAtm == 0
    Tb_V = Tb_V_toc;
    Tb_H = Tb_H_toc;
else
    r_toc_V = r_V * gamma_V^2;
    r_toc_H = r_H * gamma_H^2;
    Tb_V = atm.gossat * Tb_V_toc + atm.tb_up + ...
           atm.gossat * r_toc_V * (atm.tb_dn + atm.gossat * atm.T_cosmic);
    Tb_H = atm.gossat * Tb_H_toc + atm.tb_up + ...
           atm.gossat * r_toc_H * (atm.tb_dn + atm.gossat * atm.T_cosmic);
end

% ─── Pack output ──────────────────────────────────────────────────────────
out.VOD      = vod.total;
out.Tb_V     = Tb_V;
out.Tb_H     = Tb_H;
out.gamma_V  = gamma_V;
out.gamma_H  = gamma_H;
out.eS_V     = eS_V;
out.eS_H     = eS_H;
out.r_V      = r_V;
out.r_H      = r_H;
out.eps_soil = aux_soil.eps_soil;
out.smc_eff  = aux_soil.smc_eff;
out.delta_p  = aux_soil.delta_p;
out.vod_comp = vod;
out.iDC_top  = aux_soil.iDC_top;
out.n_frozen = aux_soil.n_frozen;
out.Tg_eff   = Tg;
out.iRough   = params.calib.iRough;
out.iLayer   = params.calib.iLayer;
out.iAtm     = params.calib.iAtm;
out.tau_H    = tau_H;
out.tau_V    = tau_V;
out.Tb_V_toc = Tb_V_toc;
out.Tb_H_toc = Tb_H_toc;
out.tau_atm  = atm.tau_atm;
out.tb_up    = atm.tb_up;
out.tb_dn    = atm.tb_dn;

end
