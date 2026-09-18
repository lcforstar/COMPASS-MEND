function TransportCoefficient = calculateThermalLiquidDispersion(SoilVariables, SoilProperties, VanGenuchten, cfg, Delt_t, TransportCoefficient)
%CALCULATETHERMALLIQUIDDISPERSION Thermal/wetting liquid transport coefficients.
%   Transport coefficients for absorbed liquid flow due to temperature
%   gradient and heat of wetting.
%
%   Groenevelt, P. H., and B. D. Kay (1974), Soil Sci. Soc. Am. Proc.,
%   38, 400-404.
%   Zeng et al. (2011), J. Geophys. Res. Atmos., 116(D20),
%   doi:10.1029/2011JD015835.
% Editor: Jizixing Zhu; Date: 2026-07-14.
% Change: 2026-07-14 Jizixing Zhu -- precompute dW_dh (Prunty) and
% dW_dTheta (Milly) for the CTh wetting-heat chain rule (§4.3.1d).

NL         = cfg.Grid.NL;
nD         = cfg.Grid.nD;
W_Chg      = cfg.lW_Chg;
SSUR       = SoilProperties.SSUR;
RHO_bulk   = SoilProperties.Bulk_density;
Theta_s    = VanGenuchten.Theta_s;

RHOL  = equations.getPhysicalConstant('RHOL');
MU_W0 = equations.getPhysicalConstant('MU_W0');
MU1   = equations.getPhysicalConstant('MU1');
b     = equations.getPhysicalConstant('b');
W0    = equations.getPhysicalConstant('W0');

W         = zeros(NL, nD);
WW        = zeros(NL, nD);
MU_W      = zeros(NL, nD);
D_Ta      = zeros(NL, nD);
dW_dh     = zeros(NL, nD);   % ∂W/∂h  [J g^-1 cm^-1] for Prunty branch
dW_dTheta = zeros(NL, nD);   % ∂W/∂θ_l [J g^-1] for Milly branch

for ML = 1:NL
    for ND = 1:nD
        MN = ML + ND - 1;
        if W_Chg == 0
            WARG = SoilVariables.Theta_LL(ML,ND) * 1e7 / SSUR;
            if WARG < 80
                W(ML,ND)  = W0 * exp(-WARG);
                WW(ML,ND) = W0 * exp(-WARG);
            end
        else
            W(ML,ND)  = -0.2932 * SoilVariables.h(MN)  / 1000;
            WW(ML,ND) = -0.2932 * SoilVariables.hh(MN) / 1000;
            if W(ML,ND) >= 2932
                W(ML,ND) = 2932;
            end
            if WW(ML,ND) >= 2932
                WW(ML,ND) = 2932;
            end
        end

        % Wetting-heat derivatives for CTh chain-rule consistency (§4.3.1d).
        % ∂W/∂h = dW_dh + dW_dTheta * ∂θ_l/∂h; the missing term
        % -ρ_l θ_l ∂W/∂h is assembled in calculateEnergyParameters.
        if W_Chg == 0
            % Milly: W = W0 * exp(-θ_l * 1e7 / SSUR), depends on h through θ_l(h).
            if WARG < 80
                dW_dTheta(ML,ND) = -WW(ML,ND) * 1e7 / SSUR;   % [J g^-1]
            end
        else
            % Prunty: W = -0.2932 * hh / 1000, clamped at 2932.
            if WW(ML,ND) < 2932
                dW_dh(ML,ND) = -0.2932 / 1000;                % [J g^-1 cm^-1]
            end
        end

        f0   = SoilVariables.Theta_g(ML,ND)^(7/3) / Theta_s(ML)^2;
        H_W  = RHOL * WW(ML,ND) * (SoilVariables.Theta_LL(ML,ND) - SoilVariables.Theta_L(ML,ND)) / ((SSUR/RHO_bulk(ML)) * Delt_t);

        TT_MN = SoilVariables.TT(MN);
        if TT_MN < -20
            MU_W(ML,ND) = 3.71e-2;
        elseif TT_MN > 150
            MU_W(ML,ND) = 1.81e-3;
        else
            MU_W(ML,ND) = MU_W0 * exp(MU1 / (8.31441 * (TT_MN + 133.3)));
        end

        L_WT = f0 * 1e7 * 1.5550e-13 * SoilVariables.POR(ML) * H_W / (b * MU_W(ML,ND));
        if cfg.lKLT_Switch == 1
            D_Ta(ML,ND) = L_WT / (RHOL * (TT_MN + 273.15));
            D_Ta(ML,ND) = max(D_Ta(ML,ND), 0);
        end
    end
end

TransportCoefficient.W         = W;
TransportCoefficient.WW        = WW;
TransportCoefficient.MU_W      = MU_W;
TransportCoefficient.D_Ta      = D_Ta;
TransportCoefficient.dW_dh     = dW_dh;
TransportCoefficient.dW_dTheta = dW_dTheta;
end
