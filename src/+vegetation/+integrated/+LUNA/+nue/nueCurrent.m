function [NUEj, NUEc, Kj2Kc] = nueCurrent(O2a, ci, tgrow, tleaf, p)
%NUECURRENT  NUE for J and Vcmax at current conditions (CO2, T_leaf, O2).
%
%   Inputs:
%     O2a    [Pa]    air O2 partial pressure
%     ci     [Pa]    leaf inter-cellular CO2 partial pressure
%     tgrow  [°C]    10-day mean growth temperature
%     tleaf  [°C]    current leaf temperature
%     p      struct  LUNA parameters (Fc25, Fj25, Kc25, Ko25, Cp25)
%
%   Returns:
%     NUEj   [μmol CO2/(g N · s)]  electron transport NUE at current conditions
%     NUEc   [μmol CO2/(g N · s)]  carboxylation NUE at current conditions
%     Kj2Kc  [-]                    Kj/Kc ratio at current conditions
%
%   Reference:
%     CLM5 LunaMod.F90 line 1198-1229 (NUE) — 1:1 port
%
%   Author: Zhu Jizixing 2026-05-20

    tfrz = 273.15;
    rgas_J = 8.314;   % [J/(K·mol)]  matches F90 rgas*1.e-3

    Fc = vegetation.integrated.LUNA.tempresp.vcmxTKattge(tgrow, tleaf) * p.Fc25;
    Fj = vegetation.integrated.LUNA.tempresp.jmxTKattge(tgrow, tleaf) * p.Fj25;

    k_c = p.Kc25 * exp((79430.0 / (rgas_J * (25.0 + tfrz))) * (1.0 - (tfrz + 25.0) / (tfrz + tleaf)));
    k_o = p.Ko25 * exp((36380.0 / (rgas_J * (25.0 + tfrz))) * (1.0 - (tfrz + 25.0) / (tfrz + tleaf)));
    c_p = p.Cp25 * exp((37830.0 / (rgas_J * (25.0 + tfrz))) * (1.0 - (tfrz + 25.0) / (tfrz + tleaf)));
    awc = k_c * (1.0 + O2a / k_o);

    Kj = max(ci - c_p, 0.0) / (4.0 * ci + 8.0 * c_p);
    Kc = max(ci - c_p, 0.0) / (ci + awc);

    NUEj  = Kj * Fj;
    NUEc  = Kc * Fc;
    Kj2Kc = Kj / Kc;
end
