function [Retard, RetardBAR] = calculateRetardationFactors(Theta_LL, RHO_b_Kd, Retard, RetardBAR, NL)
%CALCULATERETARDATIONFACTORS  Compute NH4 retardation factors from linear adsorption.
%
% R_d = 1 + (ρ_b × Kd) / θ  — the retardation factor for solutes undergoing
% linear equilibrium adsorption. Used by the NH4 transport equation to
% modify the effective storage capacity per layer.
%
% Inputs:
%   Theta_LL    [in]  [m³/m³, NL×2]  liquid water content at layer nodes (top, bottom)
%   RHO_b_Kd      [in]  [g/cm³, NL×1]  bulk density × distribution coefficient (ρ_b·Kd)
%   Retard      [inout] [-, NL×2]   retardation factor matrix (overwritten)
%   RetardBAR   [inout] [-, NL×1]   layer-mean retardation (overwritten)
%   NL          [in]  [-]            number of soil layers
%
% Returns:
%   Retard      [-, NL×2]  R_d(ML, ND) = ρ_b·Kd / θ + 1 at each layer-node
%   RetardBAR   [-, NL×1]  layer-mean = (R_d(top) + R_d(bottom)) / 2
%
% Reference:
%   van Genuchten & Wagenet 1989 SSSAJ Eq 5
%   STEMMUS-N legacy CalculateRetardationFactors.m (renamed camelCase, refactored doc only)

    for ML = 1:NL   % Loop through each soil layer
        J = ML;  % Link J with ML, where J represents the soil type
        for ND = 1:2   % Handle the two nodes (upper and lower boundary) for each soil element
            % Calculate the retardation factor: ratio of adsorption coefficient to water content + 1
            Retard(ML, ND) = RHO_b_Kd(J) / Theta_LL(ML, ND) + 1;
        end
    end
    
    for ML = 1:NL
        % Calculate the average retardation factor for each soil layer
        RetardBAR(ML) = (Retard(ML, 1) + Retard(ML, 2)) / 2;
    end
end
