function [MatricCoeff, SoilVariables, TransportCoefficient] = calculateMatricCoefficientsArray(SoilVariables, TransportCoefficient, cfg)
%CALCULATEMATRICCOEFFICIENTSARRAY  Vectorized matric coefficients (Richards eq).
%   Loop-free equivalent of soil.water.calculateMatricCoefficients: same inputs/
%   outputs, same arithmetic and statement order. Selected via cfg.MatricKernel =
%   'array'; 'scalar' runs the reference loop.
%
%   The scalar double loop indexes node fields by MN = i + j - 1 (endpoint j=1 ->
%   node i, j=2 -> node i+1). The array version slices each NN-length node field
%   into an NL x nD matrix with column j = X(j : j+NL-1), so element (i,j) again
%   reads node i+j-1. Element/layer fields are already NL x nD (or NL x 1).
%
%   Inputs / outputs (same as the scalar):
%     SoilVariables        [inout] updates DTheta_LLh/UUh, SAVEDTheta_LLh/UUh
%     TransportCoefficient [inout] updates DTheta_LLT
%     cfg                  [in]    legacy ModelSettings
%     MatricCoeff          [out]   struct: Chh, ChT, Khh, KhT, Kha, Vvh, VvT, Chg
%                                  [1/cm or cm/s per the Richards stiffness terms]
%
%   Physics:
%     Freeze correction + thermal coupling with the legacy hPARM sign convention
%     (unsaturated h < 0, so h*(-CORh) gives positive dTheta/dT). DTheta_LLT is
%     formed from the un-limited DTheta_LLh, then all three DTheta terms are
%     clamped; the coefficients use the clamped values. Rows below
%     cfg.Grid.indxBotmLayer are left untouched (groundwater coupling), i.e. the
%     loop range i = indxBotm:NL.
%
%   See also soil.water.calculateMatricCoefficients (scalar reference path).
%   Editor: Jizixing Zhu; Date: 2026-07-13.
%   Change: 2026-07-13 Jizixing Zhu -- rename the gas pressure mobility input
%   to Kg without changing the vectorized arithmetic.

    logMatricKernelOnce();   % log the active matric kernel once

    RHOL    = equations.getPhysicalConstant('RHOL');
    Gamma0  = equations.getPhysicalConstant('Gamma0');
    Gamma_w = equations.getPhysicalConstant('Gamma_w');

    NL       = cfg.Grid.NL;
    nD       = cfg.Grid.nD;
    invRHOL  = 1 / RHOL;
    indxBotm = cfg.Grid.indxBotmLayer;   % 1 when no groundwater coupling

    DTheta_h_limit = 0.02;  % [1/cm] water-capacity derivative cap (freeze-front guard)
    DTheta_T_limit = 0.5;   % [1/K]

    % --- node fields (NN x 1) -> NL x nD via column j = X(j:j+NL-1) ---
    hh_a      = sliceNode(SoilVariables.hh(:),      NL, nD);
    h_a       = sliceNode(SoilVariables.h(:),       NL, nD);
    TT_a      = sliceNode(SoilVariables.TT(:),      NL, nD);
    T_a       = sliceNode(SoilVariables.T(:),       NL, nD);
    RHOV_a    = sliceNode(SoilVariables.RHOV(:),    NL, nD);
    DRHOVh_a  = sliceNode(SoilVariables.DRHOVh(:),  NL, nD);
    DRHOVT_a  = sliceNode(SoilVariables.DRHOVT(:),  NL, nD);
    COR_a     = sliceNode(SoilVariables.COR(:),     NL, nD);
    CORh_a    = sliceNode(SoilVariables.CORh(:),    NL, nD);
    h_frez_a  = sliceNode(SoilVariables.h_frez(:),  NL, nD);
    hh_frez_a = sliceNode(SoilVariables.hh_frez(:), NL, nD);

    % --- element/layer fields (already NL x nD) ---
    Theta_UU = SoilVariables.Theta_UU;
    Theta_U  = SoilVariables.Theta_U;
    Theta_LL = SoilVariables.Theta_LL;
    Theta_L  = SoilVariables.Theta_L;
    Theta_V  = SoilVariables.Theta_V;
    KfLh     = SoilVariables.KfL_h;

    DUUh = SoilVariables.DTheta_UUh;   % prior values (updated below where masked)
    DLLh = SoilVariables.DTheta_LLh;

    if cfg.lhThmrl
        CORhh_a = -CORh_a;
        DhUU_a = COR_a .* (hh_a + hh_frez_a - h_a - h_frez_a + ...
                 (hh_a + hh_frez_a) .* CORhh_a .* (TT_a - T_a));
        DhU_a  = COR_a .* (hh_a - h_a + hh_a .* CORhh_a .* (TT_a - T_a));
        % minDh = max(1e-6, 1e-10*max(1,|hh|,|h|,|hh_frez|,|h_frez|))  per element
        M_a = max(max(max(max(abs(hh_a), abs(h_a)), abs(hh_frez_a)), abs(h_frez_a)), 1);
        minDh_a = max(1e-6, 1e-10 .* M_a);

        cond_d = (abs(DhU_a) > minDh_a) & (abs(DhUU_a) > minDh_a);
        m_uu = cond_d & (abs(Theta_UU - Theta_U) > 1e-6) & (DUUh ~= 0);
        m_ll = cond_d & (abs(Theta_LL - Theta_L) > 1e-6) & (DLLh ~= 0);
        cand_uu = (Theta_UU - Theta_U) .* COR_a ./ DhUU_a;
        cand_ll = (Theta_LL - Theta_L) .* COR_a ./ DhU_a;
        DUUh(m_uu) = cand_uu(m_uu);
        DLLh(m_ll) = cand_ll(m_ll);

        DLLT = DLLh .* (hh_a .* (-CORh_a));   % uses un-limited DLLh (clamped below)
    else
        delta_T = TT_a - T_a;
        small   = abs(delta_T) < 1e-6;
        DLLT       = (Theta_LL - Theta_L) ./ delta_T;
        DLLT_small = DLLh .* (hh_a ./ Gamma0) .* (-0.1425 - 4.76e-4 .* TT_a);
        DLLT(small) = DLLT_small(small);
        % DUUh, DLLh unchanged in this branch
    end

    % inline limitDerivative (clamp + zero non-finite), same order as scalar 94-96
    DUUh = limitDerivArray(DUUh, DTheta_h_limit);
    DLLh = limitDerivArray(DLLh, DTheta_h_limit);
    DLLT = limitDerivArray(DLLT, DTheta_T_limit);

    % --- coefficients (use clamped DLLh/DLLT) ---
    RHOV_r = RHOV_a .* invRHOL;
    D_V    = TransportCoefficient.D_V;
    D_Vg   = TransportCoefficient.D_Vg(:);   % NL x 1, broadcast over nD

    Chh = (1 - RHOV_r) .* DLLh + Theta_V .* DRHOVh_a .* invRHOL;
    Khh = (D_V + D_Vg) .* DRHOVh_a .* invRHOL + KfLh;
    Chg = KfLh;

    ChT = zeros(NL, nD); KhT = zeros(NL, nD);
    Kha = zeros(NL, nD);
    Vvh = zeros(NL, nD); VvT = zeros(NL, nD);

    if cfg.lThmrlefc
        ChT = (1 - RHOV_r) .* DLLT + Theta_V .* DRHOVT_a .* invRHOL;
        KhT = (D_V .* TransportCoefficient.Eta + D_Vg) .* DRHOVT_a .* invRHOL + ...
              TransportCoefficient.KL_T + TransportCoefficient.D_Ta;
    end
    if cfg.lKLa_Switch
        Kha = RHOV_a .* TransportCoefficient.Kg .* invRHOL + KfLh ./ Gamma_w;
    end
    if cfg.lDVa_Switch
        V_A = TransportCoefficient.V_A(:);   % NL x 1, broadcast over nD
        Vvh = -V_A .* DRHOVh_a .* invRHOL;
        VvT = -V_A .* DRHOVT_a .* invRHOL;
    end

    % --- write back only active rows (indxBotm:NL); rows below keep prior/zero ---
    r = indxBotm:NL;
    Chh_o = zeros(NL, nD); ChT_o = zeros(NL, nD);
    Khh_o = zeros(NL, nD); KhT_o = zeros(NL, nD);
    Kha_o = zeros(NL, nD); Vvh_o = zeros(NL, nD);
    VvT_o = zeros(NL, nD); Chg_o = zeros(NL, nD);
    Chh_o(r, :) = Chh(r, :);  ChT_o(r, :) = ChT(r, :);
    Khh_o(r, :) = Khh(r, :);  KhT_o(r, :) = KhT(r, :);
    Kha_o(r, :) = Kha(r, :);  Vvh_o(r, :) = Vvh(r, :);
    VvT_o(r, :) = VvT(r, :);  Chg_o(r, :) = Chg(r, :);

    SoilVariables.DTheta_UUh(r, :)        = DUUh(r, :);
    SoilVariables.DTheta_LLh(r, :)        = DLLh(r, :);
    TransportCoefficient.DTheta_LLT(r, :) = DLLT(r, :);

    SoilVariables.SAVEDTheta_LLh = SoilVariables.DTheta_LLh;
    SoilVariables.SAVEDTheta_UUh = SoilVariables.DTheta_UUh;

    MatricCoeff.Chh = Chh_o;  MatricCoeff.ChT = ChT_o;
    MatricCoeff.Khh = Khh_o;  MatricCoeff.KhT = KhT_o;
    MatricCoeff.Kha = Kha_o;  MatricCoeff.Vvh = Vvh_o;
    MatricCoeff.VvT = VvT_o;  MatricCoeff.Chg = Chg_o;
end

% ============================================================================
function Xa = sliceNode(X, NL, nD)
% Slice an NN x 1 node field into NL x nD with column j = X(j : j+NL-1).
    Xa = zeros(NL, nD);
    for j = 1:nD
        Xa(:, j) = X(j : j + NL - 1);
    end
end

function v = limitDerivArray(v, lim)
% Vectorized limitDerivative: zero non-finite, then clamp to [-lim, lim].
    v(~isfinite(v)) = 0;
    v = min(max(v, -lim), lim);
end

function logMatricKernelOnce()
% One-time log of the active matric kernel (per process).
    persistent done
    if isempty(done)
        fprintf('[MatricKernel] array (vectorized matric-coefficient kernel) active\n');
        done = true;
    end
end
