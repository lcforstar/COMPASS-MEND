function [Smin] = mineralization(TT, K_mineralization, Theta_LL, Theta_f, NL, Tr, R, eT, eTheta, MSOC_volume, Smin, PH1, ePH)
%MINERALIZATION  N mineralization rate (organic N → NH4) — single SOM pool, DAISY-style.
%
% Formulation:
%   Smin = K_min × MSOC_volume × e_θ × e_T × e_pH      [μg/(cm³ 土·s)]
%
% Single-pool first-order mineralization scaled by SOM volume fraction (DAISY-inspired,
% see Hansen et al. 1991). The "substrate amount" is taken as the soil organic matter
% volumetric content MSOC_volume; the constants tying MSOC to N mass and time
% (organic N fraction in SOM, SOM bulk density, etc.) are absorbed into K_mineralization
% during calibration. K_mineralization Excel values were tuned with this formulation.
%
% Each modulator (eT, eθ, ePH) clipped to [0, 1]. Activation energy Ea = 57142.1 J/mol
% (Stenger et al. 1995, Plant Soil 171:323).
%
% Inputs:
%   TT                [in]    [°C, NN×1]      soil temperature at nodes
%   K_mineralization  [in]    [1/s, NL×1]     reference mineralization rate (Excel /86400)
%   Theta_LL          [in]    [m³/m³, NL×2]   liquid water content at layer nodes
%   Theta_f           [in]    [m³/m³, NL×1]   field capacity per layer
%   NL                [in]    [-]              number of layers
%   Tr                [in]    [°C]             reference temperature (= 20 °C)
%   R                 [in]    [J/(mol·K)]     gas constant = 8.314
%   eT                [inout] [-, NN×1]       temperature modulator (overwritten)
%   eTheta            [inout] [-, NL×2]       moisture modulator (overwritten)
%   MSOC_volume       [in]    [m³/m³, NL×1]   SOM volumetric fraction (substrate basis)
%   Smin              [inout] [μg/(cm³·s), NL×2]  mineralization rate matrix (overwritten)
%   PH1               [in]    [-, NL×1]       soil pH per layer
%   ePH               [inout] [-, NL×1]       pH modulator (overwritten)
%
% Returns:
%   Smin              [μg/(cm³·s), NL×2]  mineralization rate at each layer-node
%
% Modulators:
%   Arrhenius:  TD = (T − Tr) / (R × (T+273.15) × (Tr+273.15))
%               e_T = exp(Ea × TD)                      [-, clipped to 1]
%   Walker:     e_θ = min(1, (θ / θ_f)^0.7)
%   pH sigmoid: e_pH = 1 / (1 + exp(−2.5 × (pH − 5.0)))  (Eckersten 1996)
%
% Reference:
%   Hansen et al. 1991 DAISY model — SOM pool basis for mineralization
%   Stenger, Priesack, Beese 1995 Plant and Soil 171:323 (Ea calibration)
%   Eckersten, Jansson, Johnsson 1996 SOILN — pH sigmoid
%
% History:
%   2026-05-14: rewrote substrate basis from (θ + ρ_b·Kd) → MSOC_volume (organic
%     substrate) per HYDRUS/DAISY convention. NH4 adsorption (ρ_b·Kd) belongs to the
%     NH4 transport equation (calculateRetardationFactors), not to mineralization rate.

    EA = 57142.1;
    TD = zeros(NL+1, 1);

    for ML = 1:NL
        J = ML;
        for ND = 1:2
            MN = ML + ND - 1;

            % Arrhenius temperature factor
            TD(MN) = (TT(MN) - Tr) / R / (TT(MN) + 273.15) ./ (Tr + 273.15);
            eT(MN) = exp(EA * TD(MN));

            % Walker moisture factor (capped at field capacity)
            eTheta(ML, ND) = min(1, (Theta_LL(ML, ND) / Theta_f(ML)).^0.7);

            % Eckersten pH sigmoid (Eckersten et al. 1996)
            ePH(ML) = 1 / (1 + exp(-2.5 * (PH1(J) - 5.0)));

            if eT(MN) > 1,        eT(MN)        = 1; end
            if eTheta(ML, ND) > 1, eTheta(ML, ND) = 1; end
            if ePH(ML) > 1,        ePH(ML)        = 1; end

            % Single SOM pool first-order mineralization
            Smin(ML, ND) = K_mineralization(J) * MSOC_volume(J) ...
                * eTheta(ML, ND) * eT(MN) * ePH(ML);
        end
    end
end
