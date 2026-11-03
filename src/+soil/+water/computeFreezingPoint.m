function [TT_CRIT, h_frez, partials] = computeFreezingPoint(hh, TT, T0, g, L_f, hd, Tmin)
%COMPUTEFREEZINGPOINT  Compute Clapeyron-consistent soil freezing temperature and head.
%
%   [TT_CRIT, h_frez] = soil.water.computeFreezingPoint(hh, TT, T0, g, L_f, hd, Tmin)
%   returns the critical freezing temperature and the additional matric head
%   used by the soil-freezing characteristic curve (SFCC). The calculation is
%   pure and vectorized; it does not update model state.
%
%   Inputs:
%     hh    [cmH2O, NN x 1] Trial matric head; negative when unsaturated
%     TT    [degC, NN x 1]  Trial soil temperature relative to T0
%     T0    [K]             Pure-water freezing temperature
%     g     [cm/s2]         Gravitational acceleration in the model CGS convention
%     L_f   [J/kg]          Latent heat of fusion; must be positive
%     hd    [cmH2O]         Oven-dry head; hh <= hd does not retain a freezing head
%     Tmin  [degC]          Lower numerical bound for freezing-point depression
%
%   Returns:
%     TT_CRIT [K, NN x 1]     Critical soil freezing temperature
%     h_frez  [cmH2O, NN x 1] Freezing head added to hh in the SFCC
%     partials [struct]        Piecewise derivatives of the returned h_frez:
%       .dh_frez_dh [cm/cm, NN x 1] Head derivative
%       .dh_frez_dT [cm/K, NN x 1]  Temperature derivative
%
%   Mathematical closure:
%     s_hT = L_f*1e4/(g*T0)                   [cmH2O/K]
%     Tcrit(h) = T0 + h/s_hT                 [K]
%     h_eq(T) = s_hT*((TT+T0)-T0)            [cmH2O]
%     h_frez = s_hT*((TT+T0)-Tcrit)          [cmH2O], if TT+T0 < Tcrit
%   On the unclamped frozen branch, h_frez = h_eq(T)-hh, hence
%   hh+h_frez = h_eq(T), dh_frez/dhh = -1, and dh_frez/dTT = s_hT.
%   The T0+min(Tmin,0), oven-dry, tiny-head, and -1e6 cm limits are retained
%   as numerical guards and have zero derivative on their clamped branches.
%
%   Reference:
%     Yu et al. 2018 JGR Atmospheres 123:7391-7414, Appendix Eqs. A2-A4,
%     p. 7409, DOI:10.1029/2018JD028502. Equation A4 has no factor 1/2.
%     The numerical guards originate from src/HT_frez.m:17-40.
%   Verification: tests/test_computeFreezingPoint_oracle.m checks the
%   analytical closure, derivatives, clamps, and dry/floor overlap.
%   Editor: Jizixing Zhu; Date: 2026-07-07.
%   Change: 2026-07-11 Jizixing Zhu -- expose piecewise derivatives of the
%           legacy freezing-head map for the conservative storage Jacobian.
%   Change: 2026-07-12 Jizixing Zhu -- remove the legacy factor 1/2, restore
%           exact Clapeyron closure, and make the oven-dry mask dominate floors.
    hh = hh(:);
    TT = TT(:);
    NN = numel(hh);

    % Yu et al. 2018 JGR Atmospheres 123:7391-7414, Appendix Eq. A4, p. 7409.
    % One kelvin of cooling corresponds to freezeHeadPerKelvin of suction head.
    freezeHeadPerKelvin = L_f * 1e4 / (g * T0);  % [cmH2O/K]
    kelvinPerHead = 1 / freezeHeadPerKelvin;      % [K/cmH2O]
    soilTemperatureK = TT + T0;                   % [K]

    % Step 1: map matric head to its equilibrium freezing temperature.
    rawTcrit = T0 + kelvinPerHead .* hh;           % [K]
    nearSaturated = hh >= -1e-6;
    rawTcrit(nearSaturated) = T0;

    % Retain the legacy numerical temperature bounds without changing the
    % interior Clapeyron slope. A bound is constant, so its derivative is zero.
    TcritFloor = T0 + min(Tmin, 0);                % [K]
    interiorTcrit = ~nearSaturated & rawTcrit > TcritFloor & rawTcrit < T0;
    TT_CRIT = min(max(rawTcrit, TcritFloor), T0);
    dTcritDh = zeros(NN, 1);
    dTcritDh(interiorTcrit) = kelvinPerHead;

    % Yu et al. 2018 JGR Atmospheres 123:7391-7414, Appendix Eq. A3, p. 7409.
    % Step 2: below Tcrit, add exactly the head needed to reach h_eq(T).
    temperatureDeficit = soilTemperatureK - TT_CRIT;  % [K]
    freezeActive = hh > hd & temperatureDeficit < 0;
    h_frez = zeros(NN, 1);
    h_frez(freezeActive) = freezeHeadPerKelvin .* temperatureDeficit(freezeActive);

    % Step 3: apply output guards. The oven-dry condition is already part of
    % freezeActive, so an oven-dry node cannot be overwritten by the -1e6 floor.
    floored = freezeActive & h_frez <= -1e6;
    tiny = freezeActive & h_frez >= -1e-6;
    h_frez(floored) = -1e6;
    h_frez(tiny) = 0;

    % Step 4: differentiate only the smooth, non-clamped branch returned above.
    differentiable = freezeActive & ~floored & ~tiny;
    dhFrezDh = zeros(NN, 1);
    dhFrezDT = zeros(NN, 1);
    dhFrezDh(differentiable) = ...
        -freezeHeadPerKelvin .* dTcritDh(differentiable);
    dhFrezDT(differentiable) = freezeHeadPerKelvin;

    partials = struct('dh_frez_dh', dhFrezDh, 'dh_frez_dT', dhFrezDT);
end
