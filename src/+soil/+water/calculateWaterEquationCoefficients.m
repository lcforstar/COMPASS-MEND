function [MatricCoeff, SoilVariables, TransportCoefficient] = ...
        calculateWaterEquationCoefficients(SoilVariables, TransportCoefficient, cfg)
%CALCULATEWATEREQUATIONCOEFFICIENTS  Assemble conservative water-equation coefficients.
%
%   [MatricCoeff, SoilVariables, TransportCoefficient] =
%   soil.water.calculateWaterEquationCoefficients(SoilVariables,
%   TransportCoefficient, cfg) builds the water-equation storage and transport
%   coefficients without legacy time-layer secants, derivative sign changes,
%   or absolute-value repairs.
%
%   Inputs:
%     SoilVariables [in] struct containing endpoint arrays [NL x nD]:
%                     Theta_V [-], KfL_h [cm/s],
%                     DThetaLiquidDh [1/cm], DThetaLiquidDT [1/K],
%                     DThetaWaterEquivDh [1/cm], DThetaWaterEquivDT [1/K],
%                     DThetaGasDh [1/cm], DThetaGasDT [1/K]; and node vectors
%                     RHOV [kg/m3], DRHOVh [kg/m3/cm], DRHOVT [kg/m3/K].
%     TransportCoefficient [inout] transport auxiliary struct. DTheta_LLT is
%                     replaced by the named physical liquid-temperature partial.
%     cfg           [in] solver settings with Grid and transport switches.
%
%   Returns:
%     MatricCoeff [out] endpoint coefficient struct: Chh [1/cm], ChT [1/K],
%                     Khh/KhT/Kha [cm/s], Vvh/VvT [1/s], Chg [cm/s].
%                     Kha belongs to the conservative physical flux;
%                     Vvh/VvT are correction-Jacobian coefficients only.
%     SoilVariables [out] unchanged compatibility container.
%     TransportCoefficient [out] with DTheta_LLT equal to DThetaLiquidDT.
%
%   Method:
%     Mw = theta_u + theta_g*rho_v/rho_l. Therefore
%       Mw,x = theta_u,x + (rho_v/rho_l)*theta_g,x
%              + (theta_g/rho_l)*rho_v,x, x in {h,T}.
%     This identity is valid for liquid-vapor water and for liquid-ice-vapor
%     water because theta_u is liquid-water-equivalent total condensed water.
%
%   References:
%     Zeng 2013, Coupled Dynamics in Soil, Eqs. 4.10 and 4.57-4.60.
%     Yu et al. 2018 JGR Atmospheres 123:7391, Appendix Eqs. A1-A5.
%   Editor: Jizixing Zhu; Date: 2026-07-11.
%   Change: 2026-07-12 Jizixing Zhu -- rename the function so its water-
%           equation scope is explicit at every call site.
%   Change: 2026-07-13 Jizixing Zhu -- read gas pressure mobility from Kg;
%           coefficient values and water fluxes are unchanged.
%   Change: 2026-07-13 Jizixing Zhu -- apply gas dispersivity over the
%           gas-filled pore fraction required by the conserved vapor flux.
%   Change: 2026-08-12 Jizixing Zhu -- separate pressure-driven vapor and
%           liquid switches and classify Vvh/VvT as Jacobian-only
%           derivatives of the conservative rho_v*Kg*grad(P_g) flux.

    RHOL = equations.getPhysicalConstant('RHOL');
    GammaW = equations.getPhysicalConstant('Gamma_w');
    invRHOL = 1 / RHOL;
    NL = cfg.Grid.NL;
    nD = cfg.Grid.nD;
    indxBotm = cfg.Grid.indxBotmLayer;

    requireEndpointFields(SoilVariables, NL, nD);
    RHOV = SoilVariables.RHOV(:);
    DRHOVh = SoilVariables.DRHOVh(:);
    DRHOVT = SoilVariables.DRHOVT(:);

    Chh = zeros(NL, nD); ChT = zeros(NL, nD);
    Khh = zeros(NL, nD); KhT = zeros(NL, nD);
    Kha = zeros(NL, nD); Vvh = zeros(NL, nD);
    VvT = zeros(NL, nD); Chg = zeros(NL, nD);

    for i = indxBotm:NL
        for j = 1:nD
            MN = i + j - 1;
            thetaGas = SoilVariables.Theta_V(i, j);
            dUdh = SoilVariables.DThetaWaterEquivDh(i, j);
            dUdT = SoilVariables.DThetaWaterEquivDT(i, j);
            dGdh = SoilVariables.DThetaGasDh(i, j);
            dGdT = SoilVariables.DThetaGasDT(i, j);
            rhovRatio = RHOV(MN) * invRHOL;

            % Zeng 2013 Eqs. 4.57-4.60, generalized by explicit ice volume.
            Chh(i, j) = dUdh + rhovRatio * dGdh ...
                + thetaGas * DRHOVh(MN) * invRHOL;
            ChT(i, j) = dUdT + rhovRatio * dGdT ...
                + thetaGas * DRHOVT(MN) * invRHOL;

            Dv = TransportCoefficient.D_V(i, j);
            % Zeng (2013), Coupled Dynamics in Soil, Eqs. 4.8-4.10,
            % pp. 65-66: dispersive vapor flux uses theta_g*D_Vg.
            Dvg = thetaGas * TransportCoefficient.D_Vg(i);
            KfLh = SoilVariables.KfL_h(i, j);
            Khh(i, j) = (Dv + Dvg) * DRHOVh(MN) * invRHOL + KfLh;
            Chg(i, j) = KfLh;

            if cfg.lThmrlefc
                KhT(i, j) = (Dv * TransportCoefficient.Eta(i, j) + Dvg) ...
                    * DRHOVT(MN) * invRHOL + TransportCoefficient.KL_T(i, j) ...
                    + TransportCoefficient.D_Ta(i, j);
            end
            if cfg.lDVa_Switch
                Kha(i, j) = Kha(i, j) + RHOV(MN) ...
                    * TransportCoefficient.Kg(i, j) * invRHOL;
                Vvh(i, j) = -TransportCoefficient.V_A(i) * DRHOVh(MN) * invRHOL;
                VvT(i, j) = -TransportCoefficient.V_A(i) * DRHOVT(MN) * invRHOL;
            end
            if cfg.lKLa_Switch
                Kha(i, j) = Kha(i, j) + KfLh / GammaW;
            end
        end
    end

    assert(all(isfinite(Chh(:))) && all(Chh(:) >= -1e-12), ...
        'calculateWaterEquationCoefficients:invalidCapacity', ...
        'Total-water pressure-head capacity must be finite and non-negative.');
    Chh(Chh < 0) = 0; % roundoff only; the assertion rejects physical negatives

    TransportCoefficient.DTheta_LLT = SoilVariables.DThetaLiquidDT;
    MatricCoeff = struct('Chh', Chh, 'ChT', ChT, 'Khh', Khh, 'KhT', KhT, ...
        'Kha', Kha, 'Vvh', Vvh, 'VvT', VvT, 'Chg', Chg);
end

function requireEndpointFields(SV, NL, nD)
% Validate the new physical interface before matrix assembly.
    names = {'Theta_V','KfL_h','DThetaLiquidDh','DThetaLiquidDT', ...
        'DThetaWaterEquivDh','DThetaWaterEquivDT','DThetaGasDh','DThetaGasDT'};
    for k = 1:numel(names)
        name = names{k};
        assert(isfield(SV, name), ...
            'calculateWaterEquationCoefficients:badField', ...
            '%s is required by the conservative storage interface.', name);
        value = SV.(name);
        assert(isequal(size(value), [NL, nD]) && all(isfinite(value(:))), ...
            'calculateWaterEquationCoefficients:badField', ...
            '%s must be a finite NL x nD endpoint array.', name);
    end
end
