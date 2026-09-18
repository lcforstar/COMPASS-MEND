function gas = buildAtmosphericGasDrivers(hourly, derived)
%BUILDATMOSPHERICGASDRIVERS  Derive atmospheric gas partial pressures.
%
%   gas = init.buildAtmosphericGasDrivers(hourly, derived) builds the
%   forcing-owned atmospheric gas driver layer from complete hourly drivers and
%   thermodynamic derivatives. These are atmospheric boundary/forcing quantities,
%   not prognostic soil-gas states.
%
%   Inputs:
%     hourly.Pg   [Pa]        atmospheric pressure
%     hourly.Ca   [ppm]       atmospheric CO2 mole fraction
%     hourly.Oa   [mmol/mol]  atmospheric O2 mole fraction
%     derived.ea  [Pa]        atmospheric water-vapour pressure
%     hourly.Ta   [degC]      air temperature, used for molar density
%
%   Returns:
%     gas.Pdry        [Pa]       dry-air pressure, Pg - ea
%     gas.CO2a_Pa     [Pa]       atmospheric CO2 partial pressure
%     gas.O2a_Pa      [Pa]       atmospheric O2 partial pressure
%     gas.O2a_molcm3  [mol/cm3]  atmospheric O2 molar density
%
%   Physical closure:
%     Atmospheric mole fractions are converted to partial pressures against
%     dry-air pressure, so water vapour displaces CO2/O2 consistently:
%       Pdry = Pg - ea
%       pCO2 = Ca*1e-6 * Pdry
%       pO2  = Oa*1e-3 * Pdry
%       cO2  = pO2 / (R*T) converted from mol/m3 to mol/cm3.
%
%   Editor: Jizixing Zhu; Date: 2026-06-27.

    phys = equations.physicalConstants();

    Pg = hourly.Pg(:);       % [Pa]
    ea = derived.ea(:);      % [Pa]
    Ta = hourly.Ta(:);       % [degC]
    Ca = hourly.Ca(:);       % [ppm]
    Oa = hourly.Oa(:);       % [mmol/mol]

    gas = struct();
    gas.Pdry       = max(Pg - ea, 1);                         % [Pa]
    gas.CO2a_Pa    = Ca * 1e-6 .* gas.Pdry;                   % [Pa]
    gas.O2a_Pa     = Oa * 1e-3 .* gas.Pdry;                   % [Pa]
    gas.O2a_molcm3 = gas.O2a_Pa ./ (phys.R .* (Ta + phys.C_to_K)) * 1e-6;
end
