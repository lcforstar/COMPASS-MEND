function TransportCoefficient = calculateThermalConductivityCapacity(SoilVariables, SoilProperties, ThermalConductivity, VanGenuchten, cfg, L, TransportCoefficient)
%CALCULATETHERMALCONDUCTIVITYCAPACITY  Select soil conductivity and heat capacity.
%
%   Inputs provide liquid, gas, and actual ice volume fractions [-], gas
%   densities [g cm^-3], soil bulk density [g cm^-3], method parameters,
%   vaporization heat L [J g^-1], and the inout coefficient container.
%
%   c_unsat and EHCAP are volumetric heat capacities [J cm^-3 K^-1];
%   Lambda_eff and ETCON are effective conductivities in the legacy thermal
%   kernel units; ZETA is dimensionless.
%
% lThermCond methods:
%   1 = de Vries (1963)           → Lambda_eff = ETCON
%   2 = Peters-Lidard / Johansen  → Lambda_eff = EfTCON / 100
%   3 = Tian (2016) de Vries form → Lambda_eff = TETCON
%   4 = Farouki capacity + PL     → Lambda_eff = EfTCON / 100
%
%   Reference: Zeng (2013), Coupled Dynamics in Soil, Eq. 4.43, p. 75
%   (PDF p. 89). The three-phase extension includes rho_i*c_i*theta_i.
%   Editor: Jizixing Zhu; Date: 2026-07-13.
%   Change: 2026-07-13 Jizixing Zhu -- retain ice sensible heat capacity in
%           the simplified thermal-property branch.

NL           = cfg.Grid.NL;
nD           = cfg.Grid.nD;
ThmrlCondCap = cfg.lThmrlCondCap;
ThermCond    = cfg.lThermCond;

Theta_LL = SoilVariables.Theta_LL;
Theta_II = SoilVariables.Theta_II;
Theta_g  = SoilVariables.Theta_g;
RHODA    = SoilVariables.RHODA;
RHOV     = SoilVariables.RHOV;

RHO_bulk = SoilProperties.Bulk_density;

phys    = equations.physicalConstants();
Lambda1 = equations.getPhysicalConstant('Lambda1');
Lambda2 = equations.getPhysicalConstant('Lambda2');
Lambda3 = equations.getPhysicalConstant('Lambda3');
c_a     = phys.c_a;
c_V     = equations.getPhysicalConstant('c_V');
c_L     = phys.c_L;
c_i     = equations.getPhysicalConstant('Cp_ice_SI') * 1e-3;
RHOI    = equations.getPhysicalConstant('RHOI');

c_unsat    = zeros(NL, nD);
Lambda_eff = zeros(NL, nD);
ETCON      = zeros(NL, nD);
EHCAP      = zeros(NL, nD);
ZETA       = zeros(NL, nD);

if ThmrlCondCap == 1
    [ETCON, EHCAP, TETCON, EfTCON, ZETA] = soil.thermal.conductivity.calculateSoilThermalProperites( ...
        SoilVariables, ThermalConductivity, VanGenuchten, cfg, L);

    for ML = 1:NL
        for ND = 1:nD
            switch ThermCond
                case {2, 4}
                    Lambda_eff(ML, ND) = EfTCON(ML, ND) / 100;
                case 3
                    Lambda_eff(ML, ND) = TETCON(ML, ND);
                otherwise  % 1 and default
                    Lambda_eff(ML, ND) = ETCON(ML, ND);
            end
            if Lambda_eff(ML, ND) <= 0
                Lambda_eff(ML, ND) = 0.0008;
            elseif Lambda_eff(ML, ND) > 0.02
                Lambda_eff(ML, ND) = 0.02;
            end
            c_unsat(ML, ND) = EHCAP(ML, ND);
        end
    end
else
    for ML = 1:NL
        for ND = 1:nD
            MN = ML + ND - 1;
            Lambda_eff(ML, ND) = Lambda1 + Lambda2 * Theta_LL(ML, ND) + Lambda3 * Theta_LL(ML, ND)^0.5;
            % Zeng (2013), Coupled Dynamics in Soil, Eq. 4.43, p. 75.
            c_unsat(ML, ND)    = 837 * RHO_bulk(ML) / 1000 + Theta_LL(ML, ND) * c_L + ...
                                  Theta_II(ML, ND) * RHOI * c_i + ...
                                  Theta_g(ML, ND) * (RHODA(MN) * c_a + RHOV(MN) * c_V);
        end
    end
end

TransportCoefficient.c_unsat    = c_unsat;
TransportCoefficient.Lambda_eff = Lambda_eff;
TransportCoefficient.ZETA       = ZETA;
TransportCoefficient.ETCON      = ETCON;
TransportCoefficient.EHCAP      = EHCAP;
end
