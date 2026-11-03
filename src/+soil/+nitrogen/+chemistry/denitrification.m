function [Sden, SVden, Cvdenit_N2o, Cvdenit_N2] = denitrification(...
    TT, Vmax_denit, Theta_LL, No3, NL, eT, eTheta, Sden, R, Tr, ...     %#ok<INUSD>
    Cvdenit_N2o, Cvdenit_N2, SVden, KT, DeltZ, Km_denit, f_anaerobic)   %#ok<INUSD>
%DENITRIFICATION  Denitrification rate (NO3 → N2/N2O) — DAISY-style MM, soil-basis.
%
%   Formulation (A2 — symmetric with nitrification.m):
%     M_total   = θ · C_NO3_liq                        [μg N/(cm³ 土)]   substrate (NO3 不吸附, ρKd≡0)
%     mm_factor = M_total / (M_total + Km_denit)       [-]               saturation in [0,1]
%     Sden      = V_max · eT · f_anaerobic · mm_factor [μg N/(cm³ 土·s)] explicit source/sink
%
%   Unit convention (mirrors nitrification.m, ADE consistency):
%     V_max_denit  [μg N/(cm³ 土·s)]   Excel 'Vmax_denit' [μg/(cm³·d)] /86400
%                                       Soil-basis maximum rate (θ 已含). 与 Vmax_nit 同约定.
%     Km_denit     [μg N/(cm³ 土)]    Excel 'Km_denit', total substrate half-saturation.
%                                       Soil-basis (= θ · Km_liquid). 与 Km_nit 同约定.
%
%   Why M_total = θ·C (not (θ+ρKd)·C as in nitrification):
%     NO3⁻ does not adsorb to soil minerals (negative ion, repelled by clay). So
%     M_total reduces to liquid pool only: M = θ·C. NH4⁺ does adsorb, hence
%     nitrification uses (θ + ρKd)·C as substrate. Structure is identical;
%     only the substrate definition reflects ion-specific chemistry.
%
%   ADE coupling (assembleSoluteParameters.m):
%     G2 = +0.98·Snit − Sden − SinkS(:,:,2)    (NO3 ADE source/sink, explicit)
%     SVden = 0 (legacy interface, F2 = -SVden no longer active)
%
%   Pre-A2 history (FIXED 2026-05-17):
%     Old code had Sden = V·θ·eT·eθ·mm·C_NO3 which expands to V·θ·eT·eθ·C²/(C+Km)
%     — i.e. a C² 2nd-order reaction, not Michaelis-Menten. Additionally V_max was
%     [1/s] (first-order coefficient, liquid-basis), inconsistent with ADE soil-basis
%     source-term convention. A2 fix: soil-basis substrate + standard MM = symmetric
%     with nitrification, dimensionally consistent with ADE.
%
%   INPUTS
%   ------
%   TT          [°C, NN×1]    Node temperature
%   Vmax_denit  [μg N/(cm³ 土·s), NL×1]  Max denitrification rate (Excel 'Vmax_denit' /86400)
%   Theta_LL    [m³/m³, NL×2] Liquid water content at layer endpoints
%   No3         [μg N/cm³ liq, NN×1]  Node-level NO3 liquid concentration
%   NL          [-]           Number of layers
%   eT          [inout, NN×1] Temperature modulator (computed below, kept for interface)
%   eTheta      [inout, NL×2] Moisture modulator (now = f_anaerobic, O2-driven)
%   Sden        [inout, NL×2] Output: denitrification rate [μg/(cm³ 土·s)]
%   R           [J/(mol·K)]   Gas constant 8.314
%   Tr          [°C]          Reference temperature = 20°C
%   Cvdenit_N2o [inout, KT×1] Cumulative N2O from denit [μg N/(cm²·s)]
%   Cvdenit_N2  [inout, KT×1] Cumulative N2  from denit [μg N/(cm²·s)]
%   SVden       [out, NL×2]   Legacy 1st-order coefficient (= 0, F2 idle)
%   KT          [-]           Outer time index
%   DeltZ       [cm, NL×1]    Layer thickness
%   Km_denit    [μg N/cm³ 土, NL×1]  Half-saturation (Excel 'Km_denit', soil-basis)
%   f_anaerobic [-, NL×2]     O2-derived anaerobic factor [0,1] from solveO2Transport
%
%   OUTPUTS
%   -------
%   Sden        [NL×2]   Denitrification rate, soil-basis [μg N/(cm³ 土·s)]
%   SVden       [NL×2]   = 0 (legacy)
%   Cvdenit_N2o [KT×1]   N2O cumulative [μg N/(cm²·s)]
%   Cvdenit_N2  [KT×1]   N2  cumulative [μg N/(cm²·s)]
%
%   References:
%     Hansen S et al. 1991 — DAISY MM nitrification template applied symmetrically
%     Harris R F 1982 — Energetics of Nitrogen Transformations (Ea = 68576.3 J/mol)
%     Bessou C et al. 2010 EJSS — N2/N2O = 0.89/0.11 partitioning
%     Šimůnek J et al. 2008 HYDRUS-1D Technical Manual §3 — ADE soil-basis source convention
%
%   See also: soil.nitrogen.chemistry.nitrification, soil.gas.solveO2Transport,
%             soil.nitrogen.transport.assembleSoluteParameters

% --- 0. Activation energy & column vectors ------------------------------
EA = 68576.3;   % [J/mol]  Harris 1982

TT          = TT(:);
No3         = No3(:);
Vmax_denit  = Vmax_denit(:);
Km_denit    = Km_denit(:);
DeltZ       = DeltZ(:);

% --- 1. Arrhenius temperature factor ------------------------------------
TD_nodes  = (TT - Tr) ./ R ./ (TT + 273.15) ./ (Tr + 273.15);
eT_nodes  = exp(EA * TD_nodes);
eT_matrix = [eT_nodes(1:NL), eT_nodes(2:NL+1)];   % [NL×2]
No3_matrix = [No3(1:NL), No3(2:NL+1)];             % [NL×2]  C_NO3 liquid

% --- 2. Anaerobic moisture factor (replaces empirical WFPS) -------------
% f_anaerobic from solveO2Transport is the physical driver (= 1 - C_O2/(C_O2+Km_O2)).
eTheta = f_anaerobic;   % [NL×2]

% --- 3. DAISY-style MM on soil-basis substrate (A2 fix) -----------------
% NO3 不吸附 → M_total = θ·C [μg/cm³ 土]
%   mm = M_total / (M_total + Km_denit)   Km in soil-basis units
%   Sden = V_max · eT · f_anaerobic · mm
Vmax_ext = [Vmax_denit, Vmax_denit];                                 % [NL×2]  μg/(cm³ 土·s)
Km_ext   = [Km_denit,   Km_denit];                                   % [NL×2]  μg/cm³ 土

M_total   = Theta_LL .* No3_matrix;                                  % [NL×2]  μg/cm³ 土
mm_factor = M_total ./ (M_total + Km_ext + eps);                     % [NL×2]  [0,1]

Sden  = Vmax_ext .* eT_matrix .* eTheta .* mm_factor;                % [μg/(cm³ 土·s)] explicit
SVden = zeros(NL, 2);                                                % legacy interface idle

% --- 4. N2/N2O partitioning (Bessou et al. 2010 EJSS) -------------------
Sden_N2o = 0.11 * Sden;
Sden_N2  = 0.89 * Sden;

% --- 5. Depth-integrated cumulative fluxes ------------------------------
Sden_N2o_avg    = (Sden_N2o(:,1) + Sden_N2o(:,2)) * 0.5;
Sden_N2_avg     = (Sden_N2(:,1)  + Sden_N2(:,2))  * 0.5;
Cvdenit_N2o(KT) = DeltZ' * Sden_N2o_avg;
Cvdenit_N2(KT)  = DeltZ' * Sden_N2_avg;

end
