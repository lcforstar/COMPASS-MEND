function C_O2_top = calcAtmO2(P_atm_Pa, T_surface_C, Oa_mmolmol, ea_Pa)
%CALCATMO2  Atmospheric O2 molar density for the soil-surface boundary.
%
%   C_O2_top = soil.gas.calcAtmO2(P_atm_Pa, T_surface_C, Oa_mmolmol, ea_Pa)
%   converts atmospheric O2 forcing into a gas-phase molar density for the top
%   Dirichlet boundary used by solveO2Transport. This is only a boundary value;
%   the soil O2 profile remains the prognostic state solved by solveO2Transport.
%
%   Inputs:
%     P_atm_Pa    [Pa]        atmospheric pressure; pass STEMMUS P_g/10
%     T_surface_C [degC]      boundary air/soil-gas temperature
%     Oa_mmolmol  [mmol/mol]  dry-air O2 mole fraction, default 209.5
%     ea_Pa       [Pa]        water-vapour pressure; optional
%
%   Returns:
%     C_O2_top    [mol/cm3]   atmospheric O2 molar density in gas phase
%
%   Physical closure:
%     pO2 = xO2 * (P_atm - ea) for moist air, then cO2 = pO2/(R*T). The final
%     1e-6 factor converts mol/m3 to mol/cm3. If ea_Pa is omitted, P_atm is used
%     directly, reproducing the legacy dry-air approximation.
%
%   Editor: Jizixing Zhu; Date: 2026-06-27.

    phys = equations.physicalConstants();

    if nargin < 1
        P_atm_Pa    = phys.Pstd;
        T_surface_C = 20;
    elseif nargin < 2
        T_surface_C = P_atm_Pa;
        P_atm_Pa    = phys.Pstd;
    end
    if nargin < 3 || isempty(Oa_mmolmol) || ~isfinite(Oa_mmolmol)
        Oa_mmolmol = phys.Oa_mmolmol_default;
    end
    if nargin < 4 || isempty(ea_Pa) || ~isfinite(ea_Pa)
        ea_Pa = 0;
    end
    if isempty(P_atm_Pa)
        P_atm_Pa = phys.Pstd;
    end
    badPressure = ~isfinite(P_atm_Pa) | P_atm_Pa <= 0;
    if any(badPressure(:))
        P_atm_Pa(badPressure) = phys.Pstd;
    end

    T_K = T_surface_C + phys.C_to_K;        % [K]
    xO2 = Oa_mmolmol * 1e-3;                % [-]
    Pdry = max(P_atm_Pa - ea_Pa, 1);        % [Pa]
    pO2 = xO2 .* Pdry;                        % [Pa]
    C_O2_top = pO2 ./ (phys.R .* T_K) * 1e-6; % [mol/cm3]
end
