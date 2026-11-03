function [NUEjref, NUEcref, Kj2Kcref] = nueRef(p)
%NUEREF  Reference NUE for J and Vcmax at 25°C, 380 ppm CO2 (LUNA reference state).
%
%   Inputs:
%     p  struct  LUNA parameters (from lunaParams.m): Fc25, Fj25, Kc25, Ko25,
%                Cp25, O2ref, CO2ref, forc_pbot_ref
%
%   Returns:
%     NUEjref  [μmol CO2 / (g N · s) at reference conditions]  for electron transport
%     NUEcref  [μmol CO2 / (g N · s) at reference conditions]  for carboxylation
%     Kj2Kcref [-]                                              Kj/Kc at reference
%
%   Reference:
%     CLM5 LunaMod.F90 line 1154-1192 (NUEref) — 1:1 port
%     Ali et al. 2016 GMD 9:587 — Eq.7-8 reference NUE
%
%   Author: Zhu Jizixing 2026-05-20

    tfrz = 273.15;
    rgas_J = 8.314;   % [J/(K·mol)]  matches F90 rgas*1.e-3 where rgas in J/K/kmole

    tgrow = 25.0;
    tleaf = 25.0;
    Fc = vegetation.integrated.LUNA.tempresp.vcmxTKattge(tgrow, tleaf) * p.Fc25;
    Fj = vegetation.integrated.LUNA.tempresp.jmxTKattge(tgrow, tleaf) * p.Fj25;

    CO2c = p.CO2ref * p.forc_pbot_ref * 1e-6;   % [Pa]
    O2c  = p.O2ref  * p.forc_pbot_ref * 1e-6;   % [Pa]

    k_c = p.Kc25 * exp((79430.0 / (rgas_J * (25.0 + tfrz))) * (1.0 - (tfrz + 25.0) / (tfrz + tleaf)));
    k_o = p.Ko25 * exp((36380.0 / (rgas_J * (25.0 + tfrz))) * (1.0 - (tfrz + 25.0) / (tfrz + tleaf)));
    c_p = p.Cp25 * exp((37830.0 / (rgas_J * (25.0 + tfrz))) * (1.0 - (tfrz + 25.0) / (tfrz + tleaf)));
    awc = k_c * (1.0 + O2c / k_o);

    ci = 0.7 * CO2c;
    Kj = max(ci - c_p, 0.0) / (4.0 * ci + 8.0 * c_p);
    Kc = max(ci - c_p, 0.0) / (ci + awc);

    NUEjref  = Kj * Fj;
    NUEcref  = Kc * Fc;
    Kj2Kcref = Kj / Kc;
end
