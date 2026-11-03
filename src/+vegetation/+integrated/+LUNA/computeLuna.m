function lunaOut = computeLuna(lunaIn, lunaState, p)
%COMPUTELUNA  Main entry: CLM5 LUNA daily Vcmax25/Jmax25 update.
%
%   Compute 5-pool N allocation (light-capture / electron-transport /
%   carboxylation / respiration / storage) by hill-climbing optimization,
%   then back out Vcmax25_opt = N_cb × Fc25 and Jmax25_opt = N_et × Fj25.
%   Apply daily rate limit (≤10%/d × Q10 thermal correction).
%
%   Inputs:
%     lunaIn struct  (10-day mean climate + canopy state)
%       .lnc          [g N/m² leaf]    top leaf N concentration (= N_leaf/LAI)
%       .tlai         [m²/m² ground]   total LAI
%       .nrad         [-]              # canopy layers for LUNA
%       .tlai_z       [nrad×1]         LAI per layer
%       .slatop       [m²/g C]         specific leaf area at top
%       .par240d_z    [nrad×1, W/m²]   10-day mean PAR per layer (visible band)
%       .par240x_z    [nrad×1, W/m²]   10-day max PAR per layer
%       .rabsorb      [-]              leaf absorption (1 - rho_leaf - tau_leaf)
%       .forc_pbot10  [Pa]             10-day mean air pressure
%       .relh10       [-]              10-day mean RH (0-1)
%       .CO2a10       [Pa]             10-day mean atmospheric CO2 partial pressure
%       .O2a10        [Pa]             10-day mean atmospheric O2 partial pressure
%       .rb10         [s/m]            10-day mean boundary layer resistance
%       .hourpd       [hr]             current daylight hours
%       .tair10       [°C]             10-day mean 2m T
%       .tleafd10     [°C]             10-day mean daytime leaf T
%       .tleafn10     [°C]             10-day mean nighttime leaf T
%       .tleaf10      [°C]             10-day mean leaf T (for EnzTurnoverTFactor)
%       .gpp_day      [kg C/m²/day]    yesterday GPP (≥0 → optimize; ≤0 → decay)
%       .c3psn        [-]              1 for C3 (LUNA only optimizes C3)
%
%     lunaState struct  (state carried across days)
%       .Vcmax25_z   [nrad×1]   previous Vcmax25 [μmol/m²/s]
%       .Jmax25_z    [nrad×1]   previous Jmax25 [μmol/m²/s]
%       .PNlc_z      [nrad×1]   previous N_lc proportion [-]
%       .enzs_z      [nrad×1]   enzyme abundance multiplier [0,1]
%       .first_day   [logical]  true → seed with Vcmax25=50, Jmax25=85
%
%     p struct  LUNA constants (from lunaParams.m)
%
%   Returns:
%     lunaOut struct
%       .Vcmax25_z   [nrad×1, μmol/m²/s]  updated 25°C Vcmax (current after daily limit)
%       .Jmax25_z    [nrad×1, μmol/m²/s]  updated 25°C Jmax
%       .PNlc_z      [nrad×1, -]          new light-capture N proportion
%       .enzs_z      [nrad×1, -]          enzyme abundance
%       .FNCa_z      [nrad×1, g N/m² leaf]  functional N per layer (diagnostic)
%       .PNet_z, .PNcb_z, .PNresp_z, .PNstore_z   [nrad×1, -]  4 other pool props (diag)
%
%   Reference:
%     CLM5 LunaMod.F90 line 280-449 (Update_Photosynthesis_Capacity) — 1:1 port
%     Ali et al. 2016 GMD 9:587 — LUNA v1.0
%
%   Author: Zhu Jizixing 2026-05-20
% Change: 2026-09-30 Wenjie Liu -- actual-leaf-N mode bounds all functional
% investments after capacity rate limits; compatibility mode is unchanged.
% Change: 2026-09-30 Wenjie Liu -- roll back this round's leaf-N LUNA
%   extension to the frozen COMPASS implementation; historical records above
%   are retained. Winter-wheat C/N allocation is handled separately.

    nrad = lunaIn.nrad;

    % First-day initialization (matches CLM5 cold-start: PNlc=0, see LunaMod.F90:364)
    if lunaState.first_day
        lunaOut.Vcmax25_z = ones(nrad,1) * p.Vcmax25_fallback;
        lunaOut.Jmax25_z  = ones(nrad,1) * p.Jmax25_fallback;
        lunaOut.PNlc_z    = zeros(nrad,1);   % cold start (Bug #9 audit 2026-05-20)
        lunaOut.enzs_z    = ones(nrad,1) * 1.0;
        lunaOut.FNCa_z    = zeros(nrad,1);
        lunaOut.PNet_z    = zeros(nrad,1);
        lunaOut.PNcb_z    = zeros(nrad,1);
        lunaOut.PNresp_z  = zeros(nrad,1);
        lunaOut.PNstore_z = zeros(nrad,1);
        return;
    end

    % Default carry-over
    lunaOut.Vcmax25_z = lunaState.Vcmax25_z;
    lunaOut.Jmax25_z  = lunaState.Jmax25_z;
    lunaOut.PNlc_z    = lunaState.PNlc_z;
    lunaOut.enzs_z    = lunaState.enzs_z;
    lunaOut.FNCa_z    = zeros(nrad,1);
    lunaOut.PNet_z    = zeros(nrad,1);
    lunaOut.PNcb_z    = zeros(nrad,1);
    lunaOut.PNresp_z  = zeros(nrad,1);
    lunaOut.PNstore_z = zeros(nrad,1);

    % --- Enzyme turnover rate, daily-limited capacity change ---
    EnzTurnoverTFactor = p.Q10Enz ^ (0.1 * (min(40.0, lunaIn.tleaf10) - 25.0));
    max_daily_pchg = EnzTurnoverTFactor * p.Enzyme_turnover_daily;

    % --- Guard: only optimize when LAI>0, leaf N>0, C3, and GPP>0 ---
    if lunaIn.tlai <= 0.0 || lunaIn.lnc <= 0.0
        return;   % no leaves, keep state
    end
    if round(lunaIn.c3psn) ~= 1
        return;   % skip C4
    end

    if lunaIn.gpp_day <= 0.0
        % --- Decay branch (drought or dormancy) ---
        max_daily_decay = min(0.5, 0.1 * max_daily_pchg);
        for z = 1:nrad
            if lunaOut.enzs_z(z) > 0.5
                lunaOut.enzs_z(z)    = lunaOut.enzs_z(z)    * (1.0 - max_daily_decay);
                lunaOut.Jmax25_z(z)  = lunaOut.Jmax25_z(z)  * (1.0 - max_daily_decay);
                lunaOut.Vcmax25_z(z) = lunaOut.Vcmax25_z(z) * (1.0 - max_daily_decay);
            end
        end
        return;
    end

    % --- FNCa per layer (canopy-position correction) ---
    RadTop = lunaIn.par240d_z(1) / lunaIn.rabsorb;
    PARTop = RadTop * 4.6;     % W/m² visible → μmol photons/m²/s
    SNCa = 1.0 / lunaIn.slatop * p.SNC;

    for z = 1:nrad
        if lunaIn.tlai_z(z) > 0.0
            qabs   = lunaIn.par240d_z(z) / lunaIn.rabsorb;
            PARi10 = qabs * 4.6;
        else
            PARi10 = 0.01;
        end
        relRad   = PARi10 / max(PARTop, 1e-6);
        relCLNCa = 0.1802 * log(max(relRad, 1e-6)) + 1.0;
        relCLNCa = max(0.2, min(1.0, relCLNCa));

        if 0.9 * lunaIn.lnc > SNCa
            FNCa = relCLNCa * (lunaIn.lnc - SNCa);
        else
            FNCa = relCLNCa * 0.1 * lunaIn.lnc;
        end

        % Boundary: cap unrealistically high N
        if FNCa > 15.0
            FNCa = 15.0;
        end
        lunaOut.FNCa_z(z) = FNCa;
    end

    % --- Hill-climb optimization per layer ---
    for z = 1:nrad
        FNCa = lunaOut.FNCa_z(z);
        if FNCa <= 0.01
            continue;   % skip empty layer
        end

        radmax2mean = lunaIn.par240x_z(z) / max(lunaIn.par240d_z(z), 1e-6);
        if lunaIn.tlai_z(z) > 0.0
            qabs   = lunaIn.par240d_z(z) / lunaIn.rabsorb;
            PARi10 = qabs * 4.6;
        else
            PARi10 = 0.01;
        end
        PARimx10 = PARi10 * radmax2mean;

        PNlcold   = lunaOut.PNlc_z(z);
        PNetold   = 0.0;
        PNrespold = 0.0;
        PNcbold   = 0.0;

        [PNstoreopt, PNlcopt, PNetopt, PNrespopt, PNcbopt] = ...
            vegetation.integrated.LUNA.nitrogenAllocation( ...
                FNCa, lunaIn.forc_pbot10, lunaIn.relh10, ...
                lunaIn.CO2a10, lunaIn.O2a10, PARi10, PARimx10, lunaIn.rb10, ...
                lunaIn.hourpd, lunaIn.tair10, lunaIn.tleafd10, lunaIn.tleafn10, ...
                PNlcold, PNetold, PNrespold, PNcbold, p);

        vcmx25_opt = PNcbopt * FNCa * p.Fc25;
        jmx25_opt  = PNetopt * FNCa * p.Fj25;

        % --- Daily change rate limit (≤10%/d × Q10) ---
        chg = vcmx25_opt - lunaOut.Vcmax25_z(z);
        chg_constrn = min(abs(chg), lunaOut.Vcmax25_z(z) * max_daily_pchg);
        lunaOut.Vcmax25_z(z) = lunaOut.Vcmax25_z(z) + sign(chg) * chg_constrn;

        chg = jmx25_opt - lunaOut.Jmax25_z(z);
        chg_constrn = min(abs(chg), lunaOut.Jmax25_z(z) * max_daily_pchg);
        lunaOut.Jmax25_z(z) = lunaOut.Jmax25_z(z) + sign(chg) * chg_constrn;

        lunaOut.PNlc_z(z)    = PNlcopt;
        lunaOut.PNet_z(z)    = PNetopt;
        lunaOut.PNcb_z(z)    = PNcbopt;
        lunaOut.PNresp_z(z)  = PNrespopt;
        lunaOut.PNstore_z(z) = PNstoreopt;

        if lunaOut.enzs_z(z) < 1.0
            lunaOut.enzs_z(z) = lunaOut.enzs_z(z) * (1.0 + max_daily_pchg);
        end

        % --- Fallbacks (LunaMod.F90 line 392-419) ---
        if isnan(lunaOut.Vcmax25_z(z)) || lunaOut.Vcmax25_z(z) > p.Vcmax25_max ...
                || lunaOut.Vcmax25_z(z) < 0
            lunaOut.Vcmax25_z(z) = p.Vcmax25_fallback;
        end
        if isnan(lunaOut.Jmax25_z(z)) || lunaOut.Jmax25_z(z) > p.Jmax25_max ...
                || lunaOut.Jmax25_z(z) < 0
            lunaOut.Jmax25_z(z) = p.Jmax25_fallback;
        end
    end
end
