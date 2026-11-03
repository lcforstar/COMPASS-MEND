function Cab_nlc = nLeafToCab(N_lc)
%NLEAFTOCAB  Convert LUNA light-capture leaf N to leaf chlorophyll content Cab.
%
%   Computes leaf chlorophyll content (Cab, the PROSPECT/Fluspect input) from
%   LUNA's DYNAMIC light-capture N pool (N_lc). This carries the LUNA optimal-
%   allocation signal — under N stress, the partition of N toward light-capture
%   vs Rubisco shifts and Cab responds accordingly (physically self-consistent).
%   (A total-leaf-N route was tried but removed 2026-06-02: lnc is a fixed
%    constant 1/(SLA×CN_L) under lnc_opt=false, so it produced a near-constant
%    Cab ≈ prescribed value — no dynamics, useless.)
%
%   Inputs:
%     N_lc   [g N/m² leaf]   LUNA light-capture N pool = PNlc × FNCa.
%                            Interpreted as total thylakoid N (light-harvesting
%                            + electron-transport complexes), per Evans 1989.
%
%   Returns:
%     Cab_nlc  [μg/cm²]   chlorophyll content (from N_lc, 50:1 thylakoid N:Chl).
%
%   Method (unit chain):
%     mol N      = N_g / M_N                  M_N  = 14.007 g/mol (N atomic mass)
%     mol Chl    = mol N × (Chl:N molar ratio = 1/50)
%     Chl_g_m2   = mol Chl × M_Chl            M_Chl = 893.5 g/mol (Chl a, C55H72MgN4O5)
%     Cab[μg/cm²]= Chl_g_m2 × 1e6(g→μg) / 1e4(m²→cm²) = Chl_g_m2 × 100
%
%   Reference:
%     Evans JR (1989) Oecologia 78:9-19, p.12 col.1 & Table 2
%       "50 mol thylakoid N mol⁻¹ Chl can be assumed" (p.12, also Park & Port
%       1963 in Table 2 note h) → Chl:N = 1/50 mol/mol.
%     M_Chl = 893.5 g/mol: Chl a molecular formula C55H72MgN4O5 (standard
%       chemistry constant; Evans 1989 works in mol and does not give mass).
%
%   Note:
%     N_lc is treated as total thylakoid N (the LUNA light-capture pool covers
%     light-harvesting + ETC complexes), which is the quantity the 50:1 ratio
%     refers to — NOT the 4 mol N/mol Chl of the bare pigment (Evans 1989 p.12).
%
%   Author: Zhu Jizixing 2026-06-02 (M2-B.3 LUNA → Cab coupling)

    % --- Physical constants ---
    M_N   = 14.007;     % [g/mol] nitrogen atomic mass
    M_Chl = 893.5;      % [g/mol] chlorophyll a (C55H72MgN4O5)

    % --- Evans 1989 Chl:N molar ratio (thylakoid) ---
    chlN = 1.0 / 50.0;   % [mol Chl/mol N]  Evans 1989 p.12 (50:1 thylakoid)

    % --- N_lc (LUNA dynamic light-capture N) → Cab ---
    molN    = max(N_lc, 0.0) / M_N;     % [mol N/m² leaf]
    Chl     = molN * chlN * M_Chl;      % [g Chl/m² leaf]
    Cab_nlc = Chl * 100.0;              % [μg/cm²]
end
