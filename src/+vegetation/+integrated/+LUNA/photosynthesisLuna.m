function [ci_out, Kc, Kj, A] = photosynthesisLuna( ...
    forc_pbot, tleafd, relh, CO2a, O2a, rb, Vcmax, JmeanL, p)
%PHOTOSYNTHESISLUNA  Internal Farquhar + Ball-Berry solver (LUNA NUE optimization).
%
%   Solves Farquhar (A vs Ci) + Ball-Berry (gs vs A) + diffusion
%   simultaneously by iterating ci. Two-phase approach:
%     Phase 1: Rubisco-limited (Wc ≤ Wj) — iterate ci using Wc
%     Phase 2: Light-limited (Wj > Wc)   — iterate ci using Wj
%
%   This is LUNA's INTERNAL photosynthesis used during NUE optimization
%   (Nitrogen_investments calls). NOT the production biochem.m.
%
%   Inputs:
%     forc_pbot  [Pa]     air pressure
%     tleafd     [°C]     daytime leaf temperature
%     relh       [-]      relative humidity (0-1)
%     CO2a       [Pa]     atmospheric CO2 partial pressure
%     O2a        [Pa]     atmospheric O2 partial pressure
%     rb         [s/m]    boundary layer resistance
%     Vcmax      [μmol/m²/s]  max carboxylation rate
%     JmeanL     [μmol/m²/s]  average electron transport rate
%     p          struct   LUNA parameters (Kc25, Ko25, Cp25, mp, minrelh)
%
%   Returns:
%     ci_out  [Pa]   intercellular CO2 partial pressure
%     Kc      [-]    Vcmax→Wc conversion factor
%     Kj      [-]    J→Wj conversion factor
%     A       [μmol/m²/s] gross photosynthesis rate
%
%   Reference:
%     CLM5 LunaMod.F90 line 1037-1147 (Photosynthesis_luna) — 1:1 port
%     Farquhar GD et al. (1980) Planta 149:78 — A-ci model
%     Ball JT, Berry JA (1987) — stomatal conductance
%
%   Author: Zhu Jizixing 2026-05-20

    tfrz = 273.15;
    rgas = 8.314;

    theta_cj = 0.95;
    rsmax0 = 2.0e4;
    bp = 2000.0;
    tleaf  = tleafd;
    tleafk = tleaf + tfrz;
    aquad  = 1.0;
    relhc  = max(p.minrelh, relh);
    bbb    = 1.0 / bp;
    mbb    = p.mp;

    CO2c = CO2a;
    O2c  = O2a;
    ci   = 0.7 * CO2c;
    ciold = ci - 0.02;

    cf = forc_pbot / (rgas * tleafk) * 1.0e6;
    gb_mol = cf / rb;

    % Bernacchi 2001 T responses
    k_c = p.Kc25 * exp((79430.0 / (rgas * (25.0 + tfrz))) * (1.0 - (tfrz + 25.0) / (tfrz + tleaf)));
    k_o = p.Ko25 * exp((36380.0 / (rgas * (25.0 + tfrz))) * (1.0 - (tfrz + 25.0) / (tfrz + tleaf)));
    c_p = p.Cp25 * exp((37830.0 / (rgas * (25.0 + tfrz))) * (1.0 - (tfrz + 25.0) / (tfrz + tleaf)));
    awc = k_c * (1.0 + O2c / k_o);

    % --- Phase 1: Rubisco-limited ci iteration ---
    i = 1;
    while abs(ci - ciold) > 0.01 && i < 100
        i = i + 1;
        ciold = ci;
        Kc = max(ci - c_p, 0.0) / (ci + awc);
        Wc = Kc * Vcmax;
        gs_mol = bbb + mbb * Wc / CO2c * forc_pbot * relhc;
        phi = forc_pbot * (1.37 * gs_mol + 1.6 * gb_mol) / (gb_mol * gs_mol);
        bquad = awc - CO2c + phi * Vcmax;
        cquad = -(c_p * phi * Vcmax + awc * CO2c);
        [r1, r2] = vegetation.integrated.LUNA.quadraticRoots(aquad, bquad, cquad);
        ci = max(r1, r2);
        if ci < 0.0
            ci = c_p + 0.5 * ciold;
        end
    end

    Kj = max(ci - c_p, 0.0) / (4.0 * ci + 8.0 * c_p);
    Kc = max(ci - c_p, 0.0) / (ci + awc);
    Wc = Kc * Vcmax;
    Wj = Kj * JmeanL;
    ciold = ci - 0.02;

    % --- Phase 2: Light-limited ci iteration (only if Wj < Wc) ---
    if Wj < Wc
        i = 1;
        while abs(ci - ciold) > 0.01 && i < 100
            i = i + 1;
            ciold = ci;
            gs_mol = bbb + mbb * Wj / CO2c * forc_pbot * relhc;
            phi = forc_pbot * (1.37 * gs_mol + 1.6 * gb_mol) / (gb_mol * gs_mol);
            bquad = 2.0 * c_p - CO2c + phi * JmeanL / 4.0;
            cquad = -(c_p * phi * JmeanL / 4.0 + 2.0 * c_p * CO2c);
            [r1, r2] = vegetation.integrated.LUNA.quadraticRoots(aquad, bquad, cquad);
            ci = max(r1, r2);
            if ci < 0.0
                ci = c_p + 0.5 * ciold;
            end
            Kj = max(ci - c_p, 0.0) / (4.0 * ci + 8.0 * c_p);
            Wj = Kj * JmeanL;
        end
        Kj = max(ci - c_p, 0.0) / (4.0 * ci + 8.0 * c_p);
        Kc = max(ci - c_p, 0.0) / (ci + awc);
        Wc = Kc * Vcmax;
        Wj = Kj * JmeanL;
    end

    % Hyperbolic colimitation (Wc, Wj)
    A = (1.0 - theta_cj) * max(Wc, Wj) + theta_cj * min(Wc, Wj);

    % rs (computed for completeness, not returned)
    % rs = cf / gs_mol; rs = min(rsmax0, rs);

    ci_out = ci;
end
