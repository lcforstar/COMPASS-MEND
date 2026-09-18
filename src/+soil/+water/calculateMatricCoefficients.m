function [MatricCoeff, SoilVariables, TransportCoefficient] = calculateMatricCoefficients(SoilVariables, TransportCoefficient, cfg)
% calculateMatricCoefficients  Compute matric coefficients for the Richards equation.
%                              Includes variably saturated flow, freeze
%                              correction, and thermal coupling.
%
% Inputs / outputs:
%   SoilVariables        - updated in place for DTheta_LLh, DTheta_UUh,
%                          SAVEDTheta_LLh, and SAVEDTheta_UUh
%   TransportCoefficient - updated in place for DTheta_LLT
%   cfg                  - model settings (legacy ModelSettings)
%   MatricCoeff          - coefficient struct with fields
%                          Chh, ChT, Khh, KhT, Kha, Vvh, VvT, Chg
%
% Notes:
%   - Freeze correction follows the legacy calculateMatricCoefficients logic.
%   - The thermal correction follows the legacy hPARM sign convention:
%     unsaturated h is negative, so h*(-CORh) gives positive dTheta/dT.
% Editor: Jizixing Zhu; Date: 2026-07-13.
% Change: 2026-07-13 Jizixing Zhu -- rename the gas pressure mobility input
% to Kg without changing the water coefficients.

RHOL    = equations.getPhysicalConstant('RHOL');
Gamma0  = equations.getPhysicalConstant('Gamma0');
Gamma_w = equations.getPhysicalConstant('Gamma_w');

NL         = cfg.Grid.NL;
nD         = cfg.Grid.nD;
DeltZ    = cfg.Grid.DeltZ;

% Extract NNx1 vectors from SoilVariables.
hh     = SoilVariables.hh(:);
h      = SoilVariables.h(:);
TT     = SoilVariables.TT(:);
T      = SoilVariables.T(:);
RHOV   = SoilVariables.RHOV(:);
DRHOVh = SoilVariables.DRHOVh(:);
DRHOVT = SoilVariables.DRHOVT(:);
COR    = SoilVariables.COR(:);
CORh   = SoilVariables.CORh(:);
h_frez = SoilVariables.h_frez(:);
hh_frez= SoilVariables.hh_frez(:);

invRHOL = 1 / RHOL;

% Guard against freeze-front secant explosions. Normal VG water-capacity
% derivatives in this setup are O(1e-4 to 1e-3) 1/cm; these caps are broad
% enough to leave regular values untouched while preventing singular matrix
% coefficients near saturation/freezing interfaces.
DTheta_h_limit = 0.02;  % 1/cm
DTheta_T_limit = 0.5;   % 1/K

% Preallocate coefficient matrices.
Chh = zeros(NL, nD); ChT = zeros(NL, nD);
Khh = zeros(NL, nD); KhT = zeros(NL, nD);
Kha = zeros(NL, nD); Vvh = zeros(NL, nD);
VvT = zeros(NL, nD); Chg = zeros(NL, nD);
indxBotm = cfg.Grid.indxBotmLayer;  % 1 when no groundwater coupling; set by AMainloop when active
for i = indxBotm:NL
    for j = 1:nD
        MN     = i + j - 1;
        h_curr = hh(MN);

        if cfg.lhThmrl
            CORhh = -1 * CORh(MN);
            DhUU = COR(MN) * (hh(MN) + hh_frez(MN) - h(MN) - h_frez(MN) + ...
                (hh(MN) + hh_frez(MN)) * CORhh * (TT(MN) - T(MN)));
            DhU  = COR(MN) * (hh(MN) - h(MN) + hh(MN) * CORhh * (TT(MN) - T(MN)));
            minDh = max(1e-6, 1e-10 * max([1, abs(hh(MN)), abs(h(MN)), ...
                abs(hh_frez(MN)), abs(h_frez(MN))]));

            if abs(DhU) > minDh && abs(DhUU) > minDh && ...
                    abs(SoilVariables.Theta_UU(i,j) - SoilVariables.Theta_U(i,j)) > 1e-6 && ...
                    SoilVariables.DTheta_UUh(i,j) ~= 0
                SoilVariables.DTheta_UUh(i,j) = ...
                    (SoilVariables.Theta_UU(i,j) - SoilVariables.Theta_U(i,j)) * COR(MN) / DhUU;
            end

            if abs(DhU) > minDh && abs(DhUU) > minDh && ...
                    abs(SoilVariables.Theta_LL(i,j) - SoilVariables.Theta_L(i,j)) > 1e-6 && ...
                    SoilVariables.DTheta_LLh(i,j) ~= 0
                SoilVariables.DTheta_LLh(i,j) = ...
                    (SoilVariables.Theta_LL(i,j) - SoilVariables.Theta_L(i,j)) * COR(MN) / DhU;
            end

            TransportCoefficient.DTheta_LLT(i,j) = SoilVariables.DTheta_LLh(i,j) * (hh(MN) * (-CORh(MN)));
        else
            delta_T = TT(MN) - T(MN);
            if abs(delta_T) < 1e-6
                TransportCoefficient.DTheta_LLT(i,j) = SoilVariables.DTheta_LLh(i,j) * ...
                    (h_curr / Gamma0) * (-0.1425 - 4.76e-4 * TT(MN));
            else
                TransportCoefficient.DTheta_LLT(i,j) = ...
                    (SoilVariables.Theta_LL(i,j) - SoilVariables.Theta_L(i,j)) / delta_T;
            end
        end

        SoilVariables.DTheta_UUh(i,j) = limitDerivative(SoilVariables.DTheta_UUh(i,j), DTheta_h_limit);
        SoilVariables.DTheta_LLh(i,j) = limitDerivative(SoilVariables.DTheta_LLh(i,j), DTheta_h_limit);
        TransportCoefficient.DTheta_LLT(i,j) = limitDerivative(TransportCoefficient.DTheta_LLT(i,j), DTheta_T_limit);

        % Matrix coefficients: Chh, Khh, Chg, ChT, KhT, Kha, Vvh, VvT.
        RHOV_r  = RHOV(MN) * invRHOL;
        Dh      = SoilVariables.DTheta_LLh(i, j);
        DT      = TransportCoefficient.DTheta_LLT(i, j);
        D_V_ij  = TransportCoefficient.D_V(i, j);
        D_Vg_i  = TransportCoefficient.D_Vg(i);
        KfLh_ij = SoilVariables.KfL_h(i, j);

        Chh(i, j) = (1 - RHOV_r) * Dh + SoilVariables.Theta_V(i,j) * DRHOVh(MN) * invRHOL;
        Khh_vapor = (D_V_ij + D_Vg_i) * DRHOVh(MN) * invRHOL;
        Khh(i, j) = Khh_vapor + KfLh_ij;
        Chg(i, j) = KfLh_ij;

        if cfg.lThmrlefc
            ChT(i, j) = (1 - RHOV_r) * DT + SoilVariables.Theta_V(i,j) * DRHOVT(MN) * invRHOL;
            KhT(i, j) = (D_V_ij * TransportCoefficient.Eta(i,j) + D_Vg_i) * DRHOVT(MN) * invRHOL + ...
                TransportCoefficient.KL_T(i,j) + TransportCoefficient.D_Ta(i,j);
        end

        if cfg.lKLa_Switch
            Kha(i, j) = RHOV(MN) * TransportCoefficient.Kg(i,j) * invRHOL + KfLh_ij / Gamma_w;
        end

        if cfg.lDVa_Switch
            Vvh(i, j) = -TransportCoefficient.V_A(i) * DRHOVh(MN) * invRHOL;
            VvT(i, j) = -TransportCoefficient.V_A(i) * DRHOVT(MN) * invRHOL;
        end

        % [DISABLED 2026-05-08] SUPG 删除（用户确认 v3 删 SUPG 跑得动）
        % Khh(i, j) = max(Khh(i, j), abs(Vvh(i, j)) * DeltZ(i) / 2);
        % KhT(i, j) = max(KhT(i, j), abs(VvT(i, j)) * DeltZ(i) / 2);
    end
end

% Save derivative snapshots for later reuse.
SoilVariables.SAVEDTheta_LLh = SoilVariables.DTheta_LLh;
SoilVariables.SAVEDTheta_UUh = SoilVariables.DTheta_UUh;

% Assemble output struct.
MatricCoeff.Chh = Chh;  MatricCoeff.ChT = ChT;
MatricCoeff.Khh = Khh;  MatricCoeff.KhT = KhT;
MatricCoeff.Kha = Kha;  MatricCoeff.Vvh = Vvh;
MatricCoeff.VvT = VvT;  MatricCoeff.Chg = Chg;
end

function value = limitDerivative(value, limitValue)
if ~isfinite(value)
    value = 0;
elseif value > limitValue
    value = limitValue;
elseif value < -limitValue
    value = -limitValue;
end
end
