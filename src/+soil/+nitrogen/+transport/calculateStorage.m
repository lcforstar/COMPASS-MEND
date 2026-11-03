function [cvnh4_KT, cvno3_KT] = calculateStorage(Theta_Solute_LL, C_NH4, C_NO3, RHOKG, DeltZ, NL, varargin)
%CALCULATESTORAGE  Compute total NH4 + NO3 column storage at the current step.
%
% Returns the depth-integrated mass of N across the soil column, split into:
%   - liquid phase  (θ_LL × C_liq)
%   - adsorbed phase (RHOKG × C_liq, NH4 only, linear isotherm)
%   - ice phase (θ_II × C_ice, optional, only if Theta_II/C_NH4_ice/C_NO3_ice provided)
%
% Used by the N mass balance check in solveNitrogen.
%
% Inputs:
%   Theta_Solute_LL  [in] [m³/m³, NL×nD]  effective liquid water content at solute nodes
%   C_NH4            [in] [mg/L, NN×1]    NH4 liquid-phase concentration
%   C_NO3            [in] [mg/L, NN×1]    NO3 liquid-phase concentration
%   RHOKG            [in] [g/cm³, NL×1]   bulk density × adsorption coeff KG
%   DeltZ            [in] [cm, NL×1]      layer thickness
%   NL               [in] [-]             number of layers
%   varargin: optional ice-phase storage (must come together as triplet)
%     Theta_II    [m³/m³, NL×nD]  ice content
%     C_NH4_ice   [mg/L, NL×nD]   N stored in ice
%     C_NO3_ice   [mg/L, NL×nD]
%
% Returns:
%   cvnh4_KT  [μg/cm²]  total NH4 column storage (liquid + adsorbed + ice if provided)
%   cvno3_KT  [μg/cm²]  total NO3 column storage (liquid + ice if provided)
%
% Notes:
%   - Trapezoidal rule per FE layer (top ND=1, bot ND=2).
%   - Ice-phase N is prognostic state (Skill §1); when freeze-thaw partitioning
%     is active, storage MUST include ice or NH4/NO3 closure breaks (Skill §6).
%   - Backward compatible: callers passing only liquid args get liquid+adsorbed
%     storage, same as before 2026-05-14.
%
% Reference:
%   STEMMUS_N.m (legacy) lines 890-902; van Genuchten & Wagenet 1989 SSSAJ.
%   Killawee et al. 1998 GCA 62: 3637 — ice-phase N entrapment.

has_ice = (length(varargin) >= 3);
has_granule = (length(varargin) >= 4);   % [2026-05-17 dissolution pool] optional NH4 granule pool
if has_ice
    Theta_II  = varargin{1};
    C_NH4_ice = varargin{2};
    C_NO3_ice = varargin{3};
end
if has_granule
    C_NH4_granule = varargin{4};
end

cvnh4_KT = 0;
cvno3_KT = 0;
for ML = 1:NL
    term1_nh4 = (Theta_Solute_LL(ML, 1) + RHOKG(ML)) * C_NH4(ML);
    term2_nh4 = (Theta_Solute_LL(ML, 2) + RHOKG(ML)) * C_NH4(ML + 1);
    cvnh4_KT = cvnh4_KT + DeltZ(ML) * (term1_nh4 + term2_nh4) / 2;

    term1_no3 = Theta_Solute_LL(ML, 1) * C_NO3(ML);
    term2_no3 = Theta_Solute_LL(ML, 2) * C_NO3(ML + 1);
    cvno3_KT = cvno3_KT + DeltZ(ML) * (term1_no3 + term2_no3) / 2;

    if has_ice
        % Ice-phase N: trapezoidal over ND=1 / ND=2 (both layer-indexed, not node)
        ice_NH4 = Theta_II(ML, 1) * C_NH4_ice(ML, 1) + Theta_II(ML, 2) * C_NH4_ice(ML, 2);
        ice_NO3 = Theta_II(ML, 1) * C_NO3_ice(ML, 1) + Theta_II(ML, 2) * C_NO3_ice(ML, 2);
        cvnh4_KT = cvnh4_KT + DeltZ(ML) * ice_NH4 / 2;
        cvno3_KT = cvno3_KT + DeltZ(ML) * ice_NO3 / 2;
    end

    if has_granule
        % [2026-05-17 dissolution pool] Granule pool: NH4 储量内置 granule
        % granule = μg/cm³ (no θ_LL factor — 表层 surface 颗粒, 不是液相溶解)
        gran_nh4 = C_NH4_granule(ML) + C_NH4_granule(ML + 1);
        cvnh4_KT = cvnh4_KT + DeltZ(ML) * gran_nh4 / 2;
    end
end

end
