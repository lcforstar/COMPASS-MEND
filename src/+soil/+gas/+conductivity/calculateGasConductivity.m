function TransportCoefficient = calculateGasConductivity(SoilVariables, SoilProperties, VanGenuchten, TransportCoefficient, cfg)
%CALCULATEGASCONDUCTIVITY  Compute gas-filled porosity and intrinsic gas permeability.
%
%   Inputs:
%     SoilVariables.Theta_g    [-, NL x nD] Gas-filled pore volume.
%     SoilVariables.POR        [-, NL x 1]  Total geometric porosity.
%     SoilVariables.Imped      [-, NN x 1]  Frozen-soil impedance exponent.
%     SoilVariables.Ratio_ice  [-, NL x nD] Ice-to-water mass ratio.
%     SoilProperties.Ks        [cm s^-1, NL x 1] Saturated water conductivity.
%     VanGenuchten.Theta_s     [-, NL x 1]  Saturated VG water content.
%     VanGenuchten.Theta_r     [-, NL x 1]  Residual VG water content.
%     VanGenuchten.m           [-, NL x 1]  VG-Mualem shape parameter, 1-1/n.
%     TransportCoefficient.MU_W [g cm^-1 s^-1, NL x nD] Water viscosity.
%
%   Updated outputs:
%     TransportCoefficient.Sa  [-, NL x nD] Actual gas-filled fraction,
%         Theta_g/POR, retained for gas dispersion and pore-volume geometry.
%     TransportCoefficient.k_g [cm^2, NL x nD] Intrinsic gas permeability.
%
%   Method:
%     The permeability closure uses a separate effective mobile-gas
%     saturation. The offset POR-Theta_s is the non-water-filled pore volume
%     that remains when the VG water curve reaches Theta_s; it is excluded
%     from the connected-gas interval. In unfrozen soil this gives exactly
%     SaEff = 1-Se, with Se=(Theta_l-Theta_r)/(Theta_s-Theta_r). Using Theta_g
%     in the numerator also prevents pore-filling ice from being counted as
%     mobile gas. The existing frozen-soil impedance is then applied as a
%     separate empirical reduction.
%
%   Reference:
%     Lenhard and Garcia-Rincon (2023), Historical Development of
%     Constitutive Relations for Addressing Subsurface LNAPL Contamination,
%     Sect. 2.4, Eqs. (2.9) and (2.16), pp. 21-39. Equation (2.16)
%     reproduces the Parker et al. (1987) VG-Mualem air relative
%     permeability, DOI: 10.1029/WR023i004p00618.
%   Editor: Jizixing Zhu; Date: 2026-07-13.
%   Change: 2026-07-13 Jizixing Zhu -- use layerwise VG effective mobile-gas
%   saturation for Parker-Mualem permeability while preserving Theta_g/POR
%   as the actual gas-filled fraction used by downstream geometry closures.

NL = cfg.Grid.NL;
nD = cfg.Grid.nD;

Ks        = SoilProperties.Ks;
m         = VanGenuchten.m;
Theta_s   = VanGenuchten.Theta_s;
Theta_r   = VanGenuchten.Theta_r;
MU_W      = TransportCoefficient.MU_W;
Theta_g   = SoilVariables.Theta_g;
POR       = SoilVariables.POR;
Imped     = SoilVariables.Imped;
Ratio_ice = SoilVariables.Ratio_ice;

g    = equations.getPhysicalConstant('g');
RHOL = equations.getPhysicalConstant('RHOL');

Sa  = zeros(NL, nD);
k_g = zeros(NL, nD);

if cfg.lDEBUG_NAN_TRACE
    fprintf('[calculateGasConductivity] Input: Theta_g NaN=%d, MU_W NaN=%d\n', ...
        sum(isnan(Theta_g(:))), sum(isnan(MU_W(:))));
end

for ML = 1:NL
    thetaRange = max(Theta_s(ML) - Theta_r(ML), 1e-12);
    trappedGas = max(POR(ML) - Theta_s(ML), 0);
    mLayer = m(ML);
    for ND = 1:nD
        MN = ML + ND - 1;
        Sa(ML, ND) = Theta_g(ML, ND) / POR(ML);

        % Lenhard and Garcia-Rincon (2023), Sect. 2.4, Eqs. (2.9), (2.16).
        % SaEff is the VG effective saturation of connected gas. In unfrozen
        % soil it is algebraically identical to 1-Se; Theta_g also excludes
        % ice volume in the frozen extension used by COMPASS.
        SaEff = (Theta_g(ML, ND) - trappedGas) / thetaRange;
        SaEff = min(1, max(0, SaEff));

        if cfg.lDEBUG_NAN_TRACE
            if isnan(Sa(ML,ND)) || isnan(MU_W(ML,ND))
                fprintf('  [NaN@ML=%d,ND=%d] Sa=%.4f, MU_W=%.2e\n', ML, ND, Sa(ML,ND), MU_W(ML,ND));
            end
        end

        if cfg.lSWCC == 1
            % Parker et al. (1987), WRR 23:618-624; the closed form is
            % reproduced explicitly by Lenhard and Garcia-Rincon (2023),
            % Sect. 2.4, Eq. (2.16): krg=sqrt(SaEff)*
            % [1-(1-SaEff)^(1/m)]^(2m), with m=1-1/n.
            krg = sqrt(SaEff) * ...
                (1 - (1 - SaEff)^(1 / mLayer))^(2 * mLayer);
            k_g(ML, ND) = Ks(ML) * MU_W(ML, ND) * krg / (g * RHOL);
        else
            k_g(ML, ND) = 0;
        end

        if cfg.lSoilairefc == 1
            k_g(ML, ND) = k_g(ML, ND) * 10^(-Imped(MN) * Ratio_ice(ML, ND));
        else
            k_g(ML, ND) = 0;
        end

        if cfg.lDEBUG_NAN_TRACE && isnan(k_g(ML,ND))
            fprintf('  [CRITICAL@ML=%d,ND=%d] k_g=NaN. Sa=%.4f, SaEff=%.4f, MU_W=%.2e, m=%.4f\n', ...
                ML, ND, Sa(ML,ND), SaEff, MU_W(ML,ND), mLayer);
        end
    end
end

if cfg.lDEBUG_NAN_TRACE
    fprintf('[calculateGasConductivity] Output: Sa NaN=%d, k_g NaN=%d\n', ...
        sum(isnan(Sa(:))), sum(isnan(k_g(:))));
end

TransportCoefficient.Sa  = Sa;
TransportCoefficient.k_g = k_g;
end
