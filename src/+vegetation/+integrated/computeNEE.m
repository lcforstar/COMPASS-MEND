function [NEE, diag] = computeNEE(GPP, R_h, R_a_plant)
%COMPUTENEE  Net Ecosystem Exchange = R_h + R_a_plant − GPP.
%
%   Returns:
%     NEE     [kg C/(m²·s)]  Net ecosystem exchange (positive = ecosystem source)
%     diag    struct         Component breakdown for diagnostic:
%                              .GPP_signed       [kg C/m²/s] −GPP (uptake, source-positive frame)
%                              .R_a_signed       [kg C/m²/s] +R_a_plant
%                              .R_h_signed       [kg C/m²/s] +R_h
%                              .NEE_total        [kg C/m²/s] same as NEE
%
%   Inputs:
%     GPP         [kg C/(m²·s)]  Canopy gross primary production (from SCOPE biochem)
%     R_h         [kg C/(m²·s)]  Heterotrophic (soil) respiration. In COMPASS this
%                                comes from fluxes.Resp (CO₂ from soil microbial
%                                activity). Sign convention: positive = source.
%     R_a_plant   [kg C/(m²·s)]  Autotrophic plant respiration (R_m_total + R_g)
%                                from vegetation.nitrogen.computeNPP.
%                                Sign convention: positive = source.
%
%   Method (Chapin et al. 2002 "Principles of Terrestrial Ecosystem Ecology"):
%
%     NEE = ecosystem CO₂ emitted to atmosphere
%         = total respiration − total assimilation
%         = (R_a_plant + R_h) − GPP
%
%     Sign convention:
%       NEE > 0: ecosystem is a CO₂ source (net emission)
%       NEE < 0: ecosystem is a CO₂ sink (net uptake)
%
%   Reference:
%     Chapin FS et al. (2002) Principles of Terrestrial Ecosystem Ecology
%       Ch. 6 §6.4 — NEE = NPP − R_h = (GPP − R_a) − R_h
%       so −NEE = GPP − R_a − R_h, i.e. NEE = R_h + R_a − GPP (this implementation)
%     CLM5 Tech Note (Lawrence 2019) §15.1 (NEE accounting)
%
%   Phase C M2-B.1b S3 (2026-05-19): moved out of ebal.m where NEE was
%   incorrectly = (Resp − Actot) × 12 / 1e9 (missing R_a_plant). With M2-B.1b
%   LPJ NPP, full R_a_plant is now available and NEE can be closed.
%
%   Author: Zhu Jizixing 2026-05-19

    NEE = R_h + R_a_plant - GPP;

    diag.GPP_signed = -GPP;
    diag.R_a_signed = R_a_plant;
    diag.R_h_signed = R_h;
    diag.NEE_total  = NEE;
end
