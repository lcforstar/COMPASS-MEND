function [Kc, Kj, ci, Vcmax, Jmax, JmeanL, JmaxL, Net, Ncb, Nresp, PSN, RESP] = ...
    nitrogenInvestments( ...
        KcKjFlag, FNCa, Nlc, forc_pbot10, relh10, CO2a10, O2a10, ...
        PARi10, PARimx10, rb10, hourpd, tair10, tleafd10, tleafn10, ...
        Kj2Kc, JmaxCoef, Fc, Fj, NUEc, NUEj, NUEcref, NUEjref, NUEr, ...
        Kc_in, Kj_in, ci_in, p)
%NITROGENINVESTMENTS  Compute Vcmax/Jmax/N investments from N_lc + climate (10d means).
%
%   Given current N_lc (light-capture N pool) and 10-day mean climate, compute:
%     - Vcmax (current-T) via Wc:Wj balance + NUE
%     - Jmax  (current-T) via Jmaxb0 baseline + PAR-driven coefficient
%     - PSN   (daily integrated photosynthesis)
%     - RESP  (daily respiration, includes night Vcmax)
%     - N investments: Net, Ncb, Nresp
%
%   Inputs:
%     KcKjFlag   [-]      0 → call photosynthesisLuna to update Kc, Kj, ci
%                          1 → use prescribed Kc, Kj, ci
%     FNCa       [g N/m² leaf]  total functional N (after structural N removed)
%     Nlc        [g N/m² leaf]  light-capture N (this is the optimization variable)
%     forc_pbot10 [Pa]    10-day mean air pressure
%     relh10     [-]      10-day mean relative humidity
%     CO2a10     [Pa]     10-day mean atmospheric CO2 partial pressure
%     O2a10      [Pa]     10-day mean atmospheric O2 partial pressure
%     PARi10     [μmol/m²/s]  10-day mean PAR in canopy layer
%     PARimx10   [μmol/m²/s]  10-day mean 24hr-max PAR
%     rb10       [s/m]    10-day mean boundary layer resistance
%     hourpd     [hr]     daylight hours
%     tair10     [°C]     10-day mean 2m temperature
%     tleafd10   [°C]     10-day mean daytime leaf temperature
%     tleafn10   [°C]     10-day mean nighttime leaf temperature
%     Kj2Kc      [-]      current NUE ratio
%     JmaxCoef   [-]      PAR & RH coefficient (computed in nitrogenAllocation)
%     Fc, Fj     [-]      T adjustment factors × Fc25/Fj25 already applied
%     NUEc, NUEj [-]      current carboxylation & electron transport NUE
%     NUEcref, NUEjref [-]  reference NUE @ 25°C, 380ppm
%     NUEr       [μmol C/(g N·day)]  respiration NUE
%     Kc_in, Kj_in, ci_in  [in/inout]  used if KcKjFlag=1, else freshly computed
%     p          struct   LUNA constants
%
%   Returns:
%     Kc, Kj, ci      [-/Pa]  Vcmax→Wc & J→Wj factors and ci
%     Vcmax, Jmax     [μmol/m²/s]  current-T capacities (NOT 25°C)
%     JmeanL, JmaxL   [μmol/m²/s]  electron transport at mean & max PAR
%     Net, Ncb        [g N/m² leaf]  electron-transport & carboxylation N pools
%     Nresp           [g N/m² leaf]  respiration N pool
%     PSN             [g C/(m² leaf · day)]  daily photosynthesis
%     RESP            [g C/(m² leaf · day)]  daily respiration
%
%   Reference:
%     CLM5 LunaMod.F90 line 948-1028 (Nitrogen_investments) — 1:1 port
%     Ali et al. 2016 GMD 9:587 — Eq.9-11 NUE-driven N investment
%
%   Author: Zhu Jizixing 2026-05-20

    leaf_mr_vcm = 0.015;   % Scalar: leaf maintenance respiration / Vcmax (CLM5 constant)
    theta_cj = 0.95;       % co-limitation interpolation coefficient

    Kc = Kc_in; Kj = Kj_in; ci = ci_in;

    % --- Light absorption (Evans 1989 form, Cb=1.78) ---
    theta = 0.292 / (1.0 + 0.076 / (Nlc * p.Cb));
    ELTRNabsorb = theta * PARi10;

    % --- Jmax baseline + PAR-driven (LunaMod.F90 line 1006-1007) ---
    Jmaxb0act = p.Jmaxb0 * FNCa * Fj;
    Jmax = Jmaxb0act + JmaxCoef * ELTRNabsorb;
    JmaxL = theta * PARimx10 / sqrt(1.0 + (theta * PARimx10 / Jmax)^2);

    % --- Vcmax from Wc:Wj balance (NUE change) ---
    NUEchg = (NUEc / NUEcref) * (NUEjref / NUEj);
    Wc2Wj = p.Wc2Wjb0 * sqrt(NUEchg);
    Wc2Wj = min(1.0, Wc2Wj);
    Vcmax = Wc2Wj * JmaxL * Kj2Kc;

    JmeanL = theta * PARi10 / sqrt(1.0 + (ELTRNabsorb / Jmax)^2);

    % --- Photosynthesis: either solve Farquhar+BB or use prescribed Kc/Kj ---
    if KcKjFlag == 0
        [ci, Kc, Kj, A] = vegetation.integrated.LUNA.photosynthesisLuna( ...
            forc_pbot10, tleafd10, relh10, CO2a10, O2a10, rb10, Vcmax, JmeanL, p);
    else
        Wc = Kc * Vcmax;
        Wj = Kj * JmeanL;
        A = (1.0 - theta_cj) * max(Wc, Wj) + theta_cj * min(Wc, Wj);
    end

    % --- Daily PSN (gross) and RESP (with night Vcmax T-adjustment) ---
    PSN = p.Cv * A * hourpd;
    Vcmaxnight = vegetation.integrated.LUNA.tempresp.vcmxTKattge(tair10, tleafn10) ...
               / vegetation.integrated.LUNA.tempresp.vcmxTKattge(tair10, tleafd10) ...
               * Vcmax;
    RESP = p.Cv * leaf_mr_vcm * (Vcmax * hourpd + Vcmaxnight * (24.0 - hourpd));

    % --- N investments (back-solved from Vcmax/Jmax/RESP) ---
    Net   = Jmax  / Fj;
    Ncb   = Vcmax / Fc;
    Nresp = RESP  / NUEr;
end
