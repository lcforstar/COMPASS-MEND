function [PNstoreopt, PNlcopt, PNetopt, PNrespopt, PNcbopt] = nitrogenAllocation( ...
        FNCa, forc_pbot10, relh10, CO2a10, O2a10, PARi10, PARimx10, rb10, ...
        hourpd, tair10, tleafd10, tleafn10, ...
        PNlcold, PNetold, PNrespold, PNcbold, p)
%NITROGENALLOCATION  Hill-climb optimization of 5-pool N partition (Ali 2016 Eq.6).
%
%   Solves max_{N_lc} [PSN(N_lc) - C_cost(N_lc)] by hill-climbing:
%     - try N_lc + 0.02·FNCa (increase) or N_lc - 0.02·FNCa (decrease)
%     - accept the direction that yields ΔPSN > ΔC_cost
%     - continue until PNlc stable (PNlcoldi == PNlc) or 100 iterations
%
%   Inputs:
%     FNCa      [g N/m² leaf]   total functional N (= leaf N - structural N)
%     forc_pbot10, relh10, CO2a10, O2a10, PARi10, PARimx10, rb10, hourpd,
%     tair10, tleafd10, tleafn10  -- 10-day mean climate (see Nitrogen_investments)
%     PNlcold, PNetold, PNrespold, PNcbold  [-]  previous-step pool proportions
%     p         struct  LUNA constants
%
%   Returns:
%     PNstoreopt  [-]  storage N proportion at end of optimization
%     PNlcopt     [-]  light-capture N proportion (optimization variable)
%     PNetopt     [-]  electron-transport N proportion (= Net/FNCa)
%     PNrespopt   [-]  respiration N proportion (= Nresp/FNCa)
%     PNcbopt     [-]  carboxylation N proportion (= Ncb/FNCa)
%
%   Reference:
%     CLM5 LunaMod.F90 line 757-942 (NitrogenAllocation) — 1:1 port
%     Ali et al. 2016 GMD 9:587 — Eq.6 optimization objective
%
%   Author: Zhu Jizixing 2026-05-20

    % --- Setup ---
    [NUEjref, NUEcref, Kj2Kcref] = vegetation.integrated.LUNA.nue.nueRef(p); %#ok<ASGLU>
    theta_cj = 0.95;  %#ok<NASGU>

    Nlc   = PNlcold   * FNCa;
    Net   = PNetold   * FNCa;
    Nresp = PNrespold * FNCa;
    Ncb   = PNcbold   * FNCa;

    if Nlc > FNCa * 0.5
        Nlc = 0.5 * FNCa;
    end
    chg_per_step = 0.02 * FNCa;
    PNlc = PNlcold;
    PNlcoldi = PNlcold - 0.001;

    PARi10c   = max(p.PARLowLim, PARi10);
    PARimx10c = max(p.PARLowLim, PARimx10);
    increase_flag = 0;
    jj = 1;

    tleafd10c = min(max(tleafd10, p.Trange1), p.Trange2);
    tleafn10c = min(max(tleafn10, p.Trange1), p.Trange2);
    ci = 0.7 * CO2a10;
    Kc = 0; Kj = 0;   % initial placeholder, recomputed in first call

    JmaxCoef = p.Jmaxb1 * ((hourpd / 12.0)^2) * (1.0 - exp(-p.relhExp * ...
               max(relh10 - p.minrelh, 0.0) / (1.0 - p.minrelh)));

    % --- Main hill-climb loop (max 100 iterations) ---
    while (PNlcoldi ~= PNlc) && (jj < 100)
        Fc = vegetation.integrated.LUNA.tempresp.vcmxTKattge(tair10, tleafd10c) * p.Fc25;
        Fj = vegetation.integrated.LUNA.tempresp.jmxTKattge(tair10, tleafd10c) * p.Fj25;
        NUEr = p.Cv * p.NUEr25 * ( ...
            vegetation.integrated.LUNA.tempresp.respTBernacchi(tleafd10c) * hourpd + ...
            vegetation.integrated.LUNA.tempresp.respTBernacchi(tleafn10c) * (24.0 - hourpd));

        % Reference call: compute Kc, Kj, ci at current Nlc
        KcKjFlag = 0;
        [NUEj, NUEc, Kj2Kc] = vegetation.integrated.LUNA.nue.nueCurrent( ...
            O2a10, ci, tair10, tleafd10c, p);

        [Kc, Kj, ci, ~, ~, ~, ~, Net, Ncb, Nresp, PSN, ~] = ...
            vegetation.integrated.LUNA.nitrogenInvestments( ...
                KcKjFlag, FNCa, Nlc, forc_pbot10, relh10, CO2a10, O2a10, ...
                PARi10c, PARimx10c, rb10, hourpd, tair10, tleafd10c, tleafn10c, ...
                Kj2Kc, JmaxCoef, Fc, Fj, NUEc, NUEj, NUEcref, NUEjref, NUEr, ...
                Kc, Kj, ci, p);

        Npsntarget = Nlc + Ncb + Net;
        PNlcoldi = Nlc / FNCa;
        Nstore = FNCa - Npsntarget - Nresp;

        % --- Try INCREASE of Nlc ---
        if Nstore > 0.0 && (increase_flag == 1 || jj == 1)
            Nlc2 = Nlc + chg_per_step;
            if Nlc2 / FNCa > 0.95
                Nlc2 = 0.95 * FNCa;
            end
            KcKjFlag = 1;
            [~, ~, ~, ~, ~, ~, ~, Net2, Ncb2, Nresp2, PSN2, ~] = ...
                vegetation.integrated.LUNA.nitrogenInvestments( ...
                    KcKjFlag, FNCa, Nlc2, forc_pbot10, relh10, CO2a10, O2a10, ...
                    PARi10c, PARimx10c, rb10, hourpd, tair10, tleafd10c, tleafn10c, ...
                    Kj2Kc, JmaxCoef, Fc, Fj, NUEc, NUEj, NUEcref, NUEjref, NUEr, ...
                    Kc, Kj, ci, p);

            Npsntarget2 = Nlc2 + Ncb2 + Net2;
            Carboncost2 = (Npsntarget2 - Npsntarget) * p.NMCp25 * p.Cv * ( ...
                vegetation.integrated.LUNA.tempresp.respTBernacchi(tleafd10c) * hourpd + ...
                vegetation.integrated.LUNA.tempresp.respTBernacchi(tleafn10c) * (24.0 - hourpd));
            Carbongain2 = PSN2 - PSN;
            if (Carbongain2 > Carboncost2) && ((Npsntarget2 + Nresp2) < 0.95 * FNCa)
                Nlc = Nlc2;
                Net = Net2;
                Ncb = Ncb2;
                Nstore = FNCa - Npsntarget2 - Nresp2;
                if jj == 1
                    increase_flag = 1;
                end
            end
        end

        % --- Try DECREASE of Nlc ---
        if increase_flag == 0
            if Nstore < 0.0
                Nlc1 = Nlc * 0.8;     % bigger step if negative
            else
                Nlc1 = Nlc - chg_per_step;
            end
            if Nlc1 < 0.05
                Nlc1 = 0.05;
            end
            KcKjFlag = 1;
            [~, ~, ~, ~, ~, ~, ~, Net1, Ncb1, Nresp1, PSN1, ~] = ...
                vegetation.integrated.LUNA.nitrogenInvestments( ...
                    KcKjFlag, FNCa, Nlc1, forc_pbot10, relh10, CO2a10, O2a10, ...
                    PARi10c, PARimx10c, rb10, hourpd, tair10, tleafd10c, tleafn10c, ...
                    Kj2Kc, JmaxCoef, Fc, Fj, NUEc, NUEj, NUEcref, NUEjref, NUEr, ...
                    Kc, Kj, ci, p);

            Npsntarget1 = Nlc1 + Ncb1 + Net1;
            Carboncost1 = (Npsntarget - Npsntarget1) * p.NMCp25 * p.Cv * ( ...
                vegetation.integrated.LUNA.tempresp.respTBernacchi(tleafd10c) * hourpd + ...
                vegetation.integrated.LUNA.tempresp.respTBernacchi(tleafn10c) * (24.0 - hourpd));
            Carbongain1 = PSN - PSN1;
            if ((Carbongain1 < Carboncost1) && Nlc1 > 0.05) || ...
                    (Npsntarget + Nresp) > 0.95 * FNCa
                Nlc = Nlc1;
                Net = Net1;
                Ncb = Ncb1;
                Nstore = FNCa - Npsntarget1 - Nresp1;
            end
        end

        PNlc = Nlc / FNCa;
        jj = jj + 1;
    end

    PNlcopt    = Nlc   / FNCa;
    PNstoreopt = Nstore / FNCa;
    PNcbopt    = Ncb   / FNCa;
    PNetopt    = Net   / FNCa;
    PNrespopt  = Nresp / FNCa;
end
