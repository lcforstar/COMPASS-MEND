function I_membrane = nacryMembraneFlux(C_root, Vmax_umol_g_h, Km, C_min, SRL)
%NACRYMEMBRANEFLUX  Root-membrane N transport flux (Michaelis-Menten on AMT/NRT).
%
%   Returns the per-unit-root-length flux through the root membrane based on
%   physiological-level Michaelis-Menten kinetics of NH4 (AMT) or NO3 (NRT2)
%   transporters per Nacry et al. 2013. Internally converts Vmax from
%   physiological units to flux units using SRL.
%
%   Returns:
%     I_membrane  [kg N / (m-root · s)]  Root-membrane N flux. Always ≥ 0.
%                                          Couples DAISY single-root flux at the
%                                          root surface (C_root) — see
%                                          solveRootSurface for the coupling.
%
%   Inputs:
%     C_root          [kg/m³]              Root-surface N concentration.
%     Vmax_umol_g_h   [μmol N/(g-root·h)]  Maximum physiological uptake rate.
%                                            Apple defaults: NH4 12, NO3 5 (Nacry 2013).
%     Km              [kg/m³]              Half-saturation concentration.
%                                            Apple default: 0.3 mg/L = 3e-4 kg/m³.
%     C_min           [kg/m³]              Membrane physical lower bound:
%                                            below this, no further uptake possible.
%                                            Apple default: NH4 0.5 mg/L, NO3 0.25 mg/L.
%     SRL             [m-root/g-root]      Specific root length, used in unit conversion.
%                                            Apple default: 40 (Wells & Eissenstat 2002).
%
%   Method:
%     1) Unit conversion (μmol/(g-root·h) → kg/(m-root·s)):
%          Vmax_kg_m_s = Vmax_umol_g_h × 14e-9 / (SRL × 3600)
%     2) Michaelis-Menten with threshold:
%          ΔC = max(C_root - C_min, 0)
%          I_membrane = Vmax_kg_m_s × ΔC / (Km + ΔC)
%
%   Reference:
%     Nacry P, Bouguyon E, Gojon A (2013) Nitrogen acquisition by roots:
%     physiological and developmental mechanisms. Plant and Soil 370:1-29.
%
%   Phase C M1 (2026-05-18): new module. Coupled with daisySingleRootFlux
%   via solveRootSurface (bisection finds C_root such that DAISY == Nacry).

% ---- Unit conversion ----
% 14e-9 = molar mass N (14 g/mol) × 1e-9 (kg/μg)... wait:
%   1 μmol N = 14 μg N = 14e-9 kg N    ✓
%   1 μmol/(g·h) × 14e-9 kg/μmol = 14e-9 kg/(g·h)
%   × 1/SRL (g-root/m-root) = 14e-9 / SRL  kg/(m-root·h)
%   / 3600 (h/s) = 14e-9 / (SRL × 3600)  kg/(m-root·s)
Vmax_kg_m_s = Vmax_umol_g_h * 14e-9 / (SRL * 3600);

% ---- Michaelis-Menten with C_min threshold ----
dC = C_root - C_min;
if dC <= 0
    I_membrane = 0;
    return;
end

I_membrane = Vmax_kg_m_s * dC / (Km + dC);

end
