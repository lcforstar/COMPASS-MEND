function f = vcmxTKattge(tgrow, tleaf)
%VCMXTKATTGE  Vcmax temperature response (Kattge & Knorr 2007, with acclimation).
%
%   Inputs:
%     tgrow  [°C]  daytime + nighttime growth temperature (10-day mean)
%     tleaf  [°C]  current leaf temperature
%
%   Returns:
%     f      [-]   Vcmax / Vcmax25 multiplicative factor
%
%   Reference:
%     Kattge J & Knorr W (2007) PCE 30:1176 — temperature acclimation
%     CLM5 LunaMod.F90 line 1234-1249 (VcmxTKattge) — 1:1 port
%
%   Author: Zhu Jizixing 2026-05-20

    tfrz = 273.15;
    % CLM5 uses rgas in J/(K·kmole) = 8314, so rgas*1e-3 in F90 = 8.314 J/(K·mol).
    % MATLAB rgas already in J/(K·mol) — no extra ×1e-3 needed.
    rgas_J = 8.314;    % [J/(K·mol)]  matches F90 `rgas*1.e-3`

    % Acclimation: limit tgrow to [11, 35]°C
    tgrow_c = min(max(tgrow, 11.0), 35.0);
    TlimVcmx = 668.39 - 1.07 * tgrow_c;

    Vcmxf1 = 1.0 + exp((TlimVcmx * (25.0 + tfrz) - 200000.0) / (rgas_J * (25.0 + tfrz)));
    Vcmxf2 = exp((72000.0 / (rgas_J * (25.0 + tfrz))) * (1.0 - (tfrz + 25.0) / (tfrz + tleaf)));
    Vcmxf3 = 1.0 + exp((TlimVcmx * (tleaf + tfrz) - 200000.0) / (rgas_J * (tleaf + tfrz)));
    f = Vcmxf1 * Vcmxf2 / Vcmxf3;
end
