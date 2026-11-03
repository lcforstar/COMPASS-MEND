function f = jmxTKattge(tgrow, tleaf)
%JMXTKATTGE  Jmax temperature response (Kattge & Knorr 2007, with acclimation).
%
%   Inputs:
%     tgrow  [°C]  daytime + nighttime growth temperature (10-day mean)
%     tleaf  [°C]  current leaf temperature
%
%   Returns:
%     f      [-]   Jmax / Jmax25 multiplicative factor
%
%   Reference:
%     Kattge J & Knorr W (2007) PCE 30:1176
%     CLM5 LunaMod.F90 line 1254-1269 (JmxTKattge) — 1:1 port
%
%   Author: Zhu Jizixing 2026-05-20

    tfrz = 273.15;
    rgas_J = 8.314;   % [J/(K·mol)]  matches F90 `rgas*1.e-3` where rgas in J/K/kmole

    tgrow_c = min(max(tgrow, 11.0), 35.0);
    TlimJmx = 659.7 - 0.75 * tgrow_c;

    Jmxf1 = 1.0 + exp((TlimJmx * (25.0 + tfrz) - 200000.0) / (rgas_J * (25.0 + tfrz)));
    Jmxf2 = exp((50000.0 / (rgas_J * (25.0 + tfrz))) * (1.0 - (tfrz + 25.0) / (tleaf + tfrz)));
    Jmxf3 = 1.0 + exp((TlimJmx * (tleaf + tfrz) - 200000.0) / (rgas_J * (tleaf + tfrz)));
    f = Jmxf1 * Jmxf2 / Jmxf3;
end
