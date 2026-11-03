function F_conc = fertilizerFlux(F_kg_ha, DEPTH)
%FERTILIZERFLUX  Convert fertilizer mass-loading (kg-N/ha) to per-volume concentration (μg/cm³).
%
% Inputs:
%   F_kg_ha [in]  [kg/ha]   single-event fertilizer N amount
%                            (legacy name 'F_tree' was misleading — was kg/ha all along;
%                             corrected 2026-05-13 by user)
%   DEPTH    [in]  [cm]      injection band thickness (= DEPTHEND − DEPTHSTART)
%
% Returns:
%   F_conc   [μg/cm³]        uniform concentration increment within the injection band.
%                              Per-node C += F_conc × ratio  (NH4 / NO3 / urea split by ratios).
%
% Unit math:
%   1 kg = 1e9 μg
%   1 ha = 1e8 cm² (area)
%   Mass per area = (F_kg_ha × 1e9) / 1e8 = F_kg_ha × 10  [μg/cm²]
%   Concentration  = mass/area / depth     = F_kg_ha × 10 / DEPTH  [μg/cm³]
%
% Example: 50 kg-N/ha across 5 cm band → F_conc = 500/5 = 100 μg/cm³
%
% Editor: Jizixing Zhu; Date: 2026-07-02.
% Change: 2026-07-02 Jizixing Zhu — align function name to file (fertilizer_flux -> fertilizerFlux).
    F_conc = (F_kg_ha * 1e9) / (1e8 * DEPTH);   % μg/cm³
end
