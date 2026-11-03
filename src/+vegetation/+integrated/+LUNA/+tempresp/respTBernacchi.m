function f = respTBernacchi(tleaf)
%RESPTBERNACCHI  Dark respiration temperature response (Bernacchi 2001).
%
%   Inputs:
%     tleaf  [°C]  leaf temperature
%
%   Returns:
%     f      [-]   Rd / Rd25 multiplicative factor
%
%   Reference:
%     Bernacchi CJ et al. (2001) PCE 24:253-259
%     CLM5 LunaMod.F90 line 1313-1318 (RespTBernacchi) — 1:1 port
%
%   Author: Zhu Jizixing 2026-05-20

    tfrz = 273.15;
    % F90 `rgas*1.e-6` where rgas in J/K/kmole = 8314 → 8.314e-3 kJ/(K·mol)
    % Bernacchi (2001) Ea = 46.39 kJ/mol (matches eq form 46.39 / (R·T))
    rgas_kJ = 8.314e-3;   % [kJ/(K·mol)]

    f = exp(18.72 - 46.39 / (rgas_kJ * (tleaf + tfrz)));
end
