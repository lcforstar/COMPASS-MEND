function EnergyVariables = calculateEnergyParameters( ...
    SoilVariables, TransportCoefficient, MatricCoeff, AirVariables, cfg, ...
    Srt, L_f, airHeatImpedance, sharedFlux)
%CALCULATEENERGYPARAMETERS  Build heat storage and transport coefficients.
%
%   EnergyVariables = soil.thermal.calculateEnergyParameters(...) evaluates
%   the endpoint coefficients in
%
%     CTh*h_t + CTT*T_t + CTa*P_g,t
%       = d/dz(KTh*h_z + KTT*T_z + KTa*P_g,z + CTg)
%         + VTh*h_z + VTT*T_z + VTa*P_g,z.
%
%   The storage Jacobian differentiates one multiphase energy state. Liquid,
%   ice, and gas-volume derivatives come from the same strict SFCC closure
%   used by the water equation. The fixed-phase sensible capacity c_unsat is
%   supplemented by phase redistribution, vapor and dry-air density changes,
%   vaporization energy, fusion energy, and wetting energy. Thus the water,
%   dry-air, and heat equations use one phase state and one set of derivatives.
%
%   Inputs and outputs retain the STEMMUS cm-g-s convention: h and z [cm],
%   T [degC], density [g cm^-3], heat capacity [J g^-1 K^-1], L_f [J kg^-1],
%   and time [s]. The returned C, K, V, and independent-flux arrays are NL-by-2
%   endpoint fields consumed by the process-specific Galerkin assembler.
%
%   References:
%     Zeng (2013), Coupled Dynamics in Soil, Eqs. 4.14-4.15 and 4.61-4.65,
%       book pp. 68-69 and 81-83 (PDF pp. 82-83 and 95-97).
%     docs/compass_soil_water_heat_physics.md, Sect. 4.3.
% Editor: Jizixing Zhu; Date: 2026-04-25.
% Change: 2026-07-13 Jizixing Zhu -- rename the gas pressure mobility field
% to Kg without changing the heat-equation coefficients.
% Change: 2026-07-13 Jizixing Zhu -- differentiate the complete liquid-ice-gas
% energy storage instead of applying the legacy freeze-only capacity override.
% Change: 2026-07-13 Jizixing Zhu -- expose actual endpoint internal energy
% and assemble root-carried liquid energy as a volumetric sink.
% Change: 2026-07-13 Jizixing Zhu -- weight gas mechanical dispersion by
% the gas-filled pore fraction in vapor and dry-air energy transport.
% Change: 2026-07-14 Jizixing Zhu -- add missing wetting-heat chain-rule
% term -ρ_l θ_l ∂W/∂h to CTh (Prunty dW_dh + Milly dW_dTheta).
% Change: 2026-07-14 Jizixing Zhu -- delegate U, CTh, CTT, CTa to the
% unified internal-energy kernel evaluateInternalEnergyEndpoints.
% Change: 2026-08-14 Jizixing Zhu -- consume
% the shared liquid, gas, vapor, dry-air, and energy flux family instead of
% reconstructing process-specific element velocities.

    NL    = cfg.Grid.NL;
    nD    = cfg.Grid.nD;
    DeltZ = cfg.Grid.DeltZ(:);
    if nargin < 8 || isempty(airHeatImpedance)
        airHeatImpedance = 0;
    end
    if nargin < 9 || isempty(sharedFlux)
        sharedFlux = soil.coupling.evaluateCoupledElementFlux( ...
            SoilVariables, TransportCoefficient, cfg);
    end

    % Physical constants in the repository's cm-g-s thermal convention.
    phys = equations.physicalConstants();
    c_a  = phys.c_a;
    c_L  = phys.c_L;
    RHOL = phys.RHOL;
    Hc   = phys.Hc;

    % Soil state and phase-consistent endpoint derivatives.
    hh             = SoilVariables.hh(:);
    hh_frez        = SoilVariables.hh_frez(:);
    TT             = SoilVariables.TT(:);
    KfL_h          = SoilVariables.KfL_h;
    Theta_LL       = SoilVariables.Theta_LL;
    Theta_g        = SoilVariables.Theta_g;
    Theta_II       = SoilVariables.Theta_II;
    DThetaLiquidDh = SoilVariables.DThetaLiquidDh;
    DThetaLiquidDT = SoilVariables.DThetaLiquidDT;
    DThetaIceDh    = SoilVariables.DThetaIceDh;
    DThetaIceDT    = SoilVariables.DThetaIceDT;
    DThetaGasDh    = SoilVariables.DThetaGasDh;
    DThetaGasDT    = SoilVariables.DThetaGasDT;
    RHODA          = SoilVariables.RHODA(:);
    Xah            = SoilVariables.Xah(:);
    XaT            = SoilVariables.XaT(:);
    Xaa            = SoilVariables.Xaa(:);
    RHOV           = SoilVariables.RHOV(:);
    DRHOVh         = SoilVariables.DRHOVh(:);
    DRHOVT         = SoilVariables.DRHOVT(:);
    L_vap          = SoilVariables.L(:);         % Vaporization energy [J g^-1]
    DRHODAz        = SoilVariables.DRHODAz(:);

    % Transport coefficients.
    D_Ta       = TransportCoefficient.D_Ta;
    D_V        = TransportCoefficient.D_V;
    D_Vg       = TransportCoefficient.D_Vg(:);
    V_A        = TransportCoefficient.V_A(:);
    Kg         = TransportCoefficient.Kg;
    Eta        = TransportCoefficient.Eta;
    KL_T       = TransportCoefficient.KL_T;
    WW         = TransportCoefficient.WW;
    dW_dh      = TransportCoefficient.dW_dh;       % ∂W/∂h  [J g^-1 cm^-1], Prunty
    dW_dTheta  = TransportCoefficient.dW_dTheta;   % ∂W/∂θ_l [J g^-1], Milly
    c_unsat    = TransportCoefficient.c_unsat;
    Lambda_eff = TransportCoefficient.Lambda_eff;

    % Water-equation coefficients reused by transported energy.
    Khh = MatricCoeff.Khh;
    KhT = MatricCoeff.KhT;
    Kha = MatricCoeff.Kha;
    Vvh = MatricCoeff.Vvh;
    VvT = MatricCoeff.VvT;

    % Element means used by the STEMMUS flux expansion.
    DTDBAR   = (D_Ta(:,1)   + D_Ta(:,2))   * 0.5;
    DEhBAR   = (D_V(:,1)    + D_V(:,2))    * 0.5;
    DETBAR   = (D_V(:,1).*Eta(:,1) + D_V(:,2).*Eta(:,2)) * 0.5;
    DRHOVhDz = (DRHOVh(2:end) + DRHOVh(1:end-1)) * 0.5;
    DRHOVTDz = (DRHOVT(2:end) + DRHOVT(1:end-1)) * 0.5;
    RHOVBAR  = (RHOV(2:end)  + RHOV(1:end-1))  * 0.5;
    EtaBAR   = (Eta(:,1)     + Eta(:,2))     * 0.5;

    % One shared element owner supplies every phase/component flux. Coefficient
    % linearizations below retain the current endpoint FEM, but no consumer
    % may reconstruct or clip a different q_l, q_g, q_v, or q_a.
    KLhBAR = mean(KfL_h, 2);
    KLTBAR = mean(KL_T, 2);
    DhDZ = sharedFlux.gradients.liquid_pressure_head_per_cm;
    DDhDZ = sharedFlux.gradients.head_per_cm;
    DTDZ = sharedFlux.gradients.temperature_K_per_cm;
    DPgDZ = sharedFlux.gradients.gas_pressure_per_cm;
    QL = sharedFlux.phase.liquid_cm_s;
    QL_h = sharedFlux.decomposition.liquid_head_gravity_cm_s;
    QL_T = sharedFlux.decomposition.liquid_temperature_cm_s;
    QL_a = sharedFlux.decomposition.liquid_gas_pressure_cm_s;
    V_A = sharedFlux.phase.gas_cm_s;
    QV = sharedFlux.component.water_vapor_mass_g_cm2_s;
    QV_h = sharedFlux.decomposition.vapor_head_g_cm2_s;
    QV_T = sharedFlux.decomposition.vapor_temperature_g_cm2_s;
    QVA = sharedFlux.decomposition.vapor_gas_advection_g_cm2_s;
    Qa = sharedFlux.component.total_dry_air_mass_g_cm2_s;
    if cfg.lSoilairefc
        Kaa = AirVariables.Kaa;
        Vaa = AirVariables.Vaa;
    else
        Kaa = zeros(NL, nD);
        Vaa = zeros(NL, nD);
    end

    % Interpolate element fluxes to coefficient endpoints.
    D_Vgg = equations.interpolateElementToNode(D_Vg, DeltZ);   % NN-by-1
    V_AA  = equations.interpolateElementToNode(V_A,  DeltZ);
    QLL   = equations.interpolateElementToNode(QL,   DeltZ);
    QVV   = equations.interpolateElementToNode(QV,   DeltZ);
    Qaa   = equations.interpolateElementToNode(Qa,   DeltZ);

    % Endpoint storage and transport coefficients.
    % CTh, CTT, CTa come from the unified kernel evaluateInternalEnergyEndpoints.
    CTT_liq = zeros(NL, nD); CTT_vap = zeros(NL, nD); CTT_air = zeros(NL, nD);
    KTh = zeros(NL, nD);  KTT = zeros(NL, nD);  KTa = zeros(NL, nD);
    VTh = zeros(NL, nD);  VTT = zeros(NL, nD);  VTa = zeros(NL, nD);
    CTg = zeros(NL, nD);
    Kcva = zeros(NL, nD);  Kcah = zeros(NL, nD);
    KcaT = zeros(NL, nD);  Kcaa = zeros(NL, nD);
    Ccah = zeros(NL, nD);  CcaT = zeros(NL, nD);  Ccaa = zeros(NL, nD);

    invRHOL = 1 / RHOL;

    for ML = 1:NL
        for ND = 1:nD
            MN = ML + ND - 1;

            tt_val    = TT(MN);
            l_val     = L_vap(MN);
            rhov_val  = RHOV(MN);
            rhoda_val = RHODA(MN);

            % Dry-air energy transport is active only with the air equation.
            if cfg.lSoilairefc
                % Zeng (2013), Coupled Dynamics in Soil, Eq. 4.13, p. 67.
                gasDispersion = Theta_g(ML,ND) * D_Vgg(MN);
                diff_dv = D_V(ML,ND) + gasDispersion;
                Kcah(ML,ND) = c_a * tt_val * (diff_dv * Xah(MN) + Hc * rhoda_val * KfL_h(ML,ND));
                KcaT(ML,ND) = c_a * tt_val * (diff_dv * XaT(MN) + Hc * rhoda_val * (KL_T(ML,ND) + D_Ta(ML,ND)));
                Kcaa(ML,ND) = c_a * tt_val * Kaa(ML,ND);

                if cfg.lDVa_Switch == 1
                    Kcva(ML,ND) = l_val * rhov_val * Kg(ML,ND);
                else
                    Kcva(ML,ND) = 0;
                end

                adv_common   = c_a * tt_val * (-V_AA(MN) - Hc * QLL(MN) * invRHOL);
                Ccah(ML,ND)  = adv_common * Xah(MN);
                CcaT(ML,ND)  = adv_common * XaT(MN);
                Ccaa(ML,ND)  = c_a * tt_val * Vaa(ML,ND);
            end

            % Storage is differentiated from one phase-resolved energy state.
            % The reference state is liquid water at 0 degC. Hence ice carries
            % -L_f and vapor carries its temperature-dependent vaporization
            % Per-phase energy densities for transport coefficients.
            % CTh, CTT, CTa and the exact internal energy U now come
            % from the unified kernel evaluateInternalEnergyEndpoints
            % (§5.9.9) — no duplicate of the energy formula remains
            % in this file.
            vaporEnergyPerMass = c_L * tt_val + l_val;          % [J/g]
            airEnergyPerMass   = c_a * tt_val;                  % [J/g]
            liquidEnergyPerVolume = RHOL * (c_L * tt_val ...
                - WW(ML,ND));                                   % [J/cm3]

            % Diagnostic decomposition of CTT (telemetry only; not used
            % by the solver).
            CTT_liq(ML,ND) = liquidEnergyPerVolume ...
                * DThetaLiquidDT(ML,ND);
            CTT_vap(ML,ND) = Theta_g(ML,ND) ...
                * vaporEnergyPerMass * DRHOVT(MN);
            CTT_air(ML,ND) = Theta_g(ML,ND) ...
                * airEnergyPerMass * XaT(MN);


            % Divergence-form transport coefficients. D_Vgg is evaluated at
            % the same endpoint as every other coefficient in this row.
            % Zeng (2013), Coupled Dynamics in Soil, Eqs. 4.9 and 4.15,
            % pp. 65 and 68-69.
            gasDispersion = Theta_g(ML,ND) * D_Vgg(MN);
            KTh(ML,ND) = l_val  * (D_V(ML,ND) + gasDispersion) * DRHOVh(MN) ...
                + c_L * tt_val * RHOL * Khh(ML,ND) ...
                + Kcah(ML,ND);

            KTT(ML,ND) = Lambda_eff(ML,ND) ...
                + c_L * tt_val * RHOL * KhT(ML,ND) ...
                + KcaT(ML,ND) ...
                + l_val * (D_V(ML,ND)*Eta(ML,ND) + gasDispersion) * DRHOVT(MN);
            KTT(ML,ND) = max(KTT(ML,ND), Lambda_eff(ML,ND));

            KTa(ML,ND) = Kcva(ML,ND) + Kcaa(ML,ND) + c_L * tt_val * RHOL * Kha(ML,ND);

            % First-order terms from the coefficient-frozen STEMMUS
            % product-rule expansion (Zeng 2013, Eq. 4.65).
            if cfg.lDVa_Switch == 1
                VTh(ML,ND) = c_L*tt_val*RHOL*Vvh(ML,ND) + Ccah(ML,ND) ...
                    - l_val * V_AA(MN) * DRHOVh(MN);
                VTT(ML,ND) = c_L*tt_val*RHOL*VvT(ML,ND) + CcaT(ML,ND) ...
                    - l_val * V_AA(MN) * DRHOVT(MN) ...
                    - (c_L*(QLL(MN)+QVV(MN)) + c_a*Qaa(MN) - 2.369*QVV(MN));
            else
                VTh(ML,ND) = c_L*tt_val*RHOL*Vvh(ML,ND) + Ccah(ML,ND);
                VTT(ML,ND) = c_L*tt_val*RHOL*VvT(ML,ND) + CcaT(ML,ND) ...
                    - (c_L*(QLL(MN)+QVV(MN)) + c_a*Qaa(MN) - 2.369*QVV(MN));
            end
            % Bound the element Peclet number without adding diffusivity;
            % this retains the accepted physical KTT.
            vtt_limit = 2 * max(abs(KTT(ML,ND)), 1e-12) / max(DeltZ(ML), 1e-12);
            VTT(ML,ND) = min(max(VTT(ML,ND), -vtt_limit), vtt_limit);
            VTa(ML,ND) = Ccaa(ML,ND);

            % Gravity drainage is a flux and therefore belongs in CTg.
            % Root uptake is a volumetric energy sink; it is stored separately
            % and integrated with N_i over the element by the FEM assembler.
            CTg(ML,ND) = (c_L*RHOL + c_a*Hc*rhoda_val) ...
                * KfL_h(ML,ND) * tt_val;

        end
    end

    % Apply optional frozen-air impedance, then expose all terms.
    if cfg.lSoilairefc && airHeatImpedance > 0
        [KTa, VTa, Kcva, Kcaa, Ccaa, airPressureTransportRatio, ...
            airPressureTransportScale, airHeatIceScale] = ...
            applyAirPressureHeatImpedance(KTa, VTa, Kcva, Kcaa, Ccaa, ...
            KTT, VTT, SoilVariables, DeltZ, airHeatImpedance);
    else
        airPressureTransportRatio = NaN;
        airPressureTransportScale = 1;
        airHeatIceScale = 1;
    end

    % ---- unified internal-energy kernel (§5.9.9) ----
    % Single source of truth for U, CTh, CTT, CTa.  All consumers
    % (coefficient refresh, time residual, column transform, T inversion)
    % call this same kernel.  No copy of the energy formula remains.
    estate = struct('T', [TT(1:NL), TT(2:NL+1)], ...
        'Theta_LL', Theta_LL, 'Theta_II', Theta_II, 'Theta_g', Theta_g, ...
        'RHOV', RHOV, 'RHODA', RHODA, 'L_vap', L_vap, ...
        'WW', WW, 'c_unsat', c_unsat);
    ederiv = struct('DThetaLiquidDh', DThetaLiquidDh, ...
        'DThetaLiquidDT', DThetaLiquidDT, ...
        'DThetaIceDh', DThetaIceDh, 'DThetaIceDT', DThetaIceDT, ...
        'DThetaGasDh', DThetaGasDh, 'DThetaGasDT', DThetaGasDT, ...
        'DRHOVh', DRHOVh, 'DRHOVT', DRHOVT, ...
        'Xah', Xah, 'XaT', XaT, 'Xaa', Xaa, ...
        'dW_dh', dW_dh, 'dW_dTheta', dW_dTheta);
    ctl = struct('L_f', L_f, 'validate', false);
    [U_endpoint, CTh_k, CTT_k, CTa_k] = soil.thermal.evaluateInternalEnergyEndpoints( ...
        estate, ederiv, ctl);

    EnergyVariables.CTh = CTh_k;   EnergyVariables.CTT = CTT_k;   EnergyVariables.CTa = CTa_k;
    EnergyVariables.CTT_liq = CTT_liq; EnergyVariables.CTT_vap = CTT_vap; EnergyVariables.CTT_air = CTT_air;
    EnergyVariables.DTheta_LLT = DThetaLiquidDT; EnergyVariables.c_unsat = c_unsat; EnergyVariables.WW = WW;
    EnergyVariables.KTh = KTh;   EnergyVariables.KTT = KTT;   EnergyVariables.KTa = KTa;
    EnergyVariables.VTh = VTh;   EnergyVariables.VTT = VTT;   EnergyVariables.VTa = VTa;
    EnergyVariables.CTg = CTg;
    EnergyVariables.rootEnergySink = RHOL * c_L * Srt .* ...
        [TT(1:NL), TT(2:NL+1)];
    EnergyVariables.energyStorage = U_endpoint;
    EnergyVariables.QL  = QL;
    EnergyVariables.Qa  = Qa;
    EnergyVariables.Qtot = QL_h + QL_T + QL_a + QV_h + QV_T + QVA;
    EnergyVariables.QL_h = QL_h; EnergyVariables.QL_T = QL_T; EnergyVariables.QL_a = QL_a;
    EnergyVariables.QVH = QV_h;  EnergyVariables.QVT = QV_T;  EnergyVariables.QVa = QVA;
    EnergyVariables.KLhBAR = KLhBAR;  EnergyVariables.KLTBAR = KLTBAR;
    EnergyVariables.DTDBAR = DTDBAR;
    EnergyVariables.DhDZ   = DhDZ;    EnergyVariables.DDhDZ  = DDhDZ;
    EnergyVariables.DTDZ   = DTDZ;    EnergyVariables.DPgDZ  = DPgDZ;
    EnergyVariables.RHOVBAR = RHOVBAR; EnergyVariables.EtaBAR = EtaBAR;
    EnergyVariables.DEhBAR  = DEhBAR;  EnergyVariables.DETBAR = DETBAR;
    EnergyVariables.Kcva = Kcva; EnergyVariables.Kcah = Kcah;
    EnergyVariables.KcaT = KcaT; EnergyVariables.Kcaa = Kcaa;
    EnergyVariables.Ccah = Ccah; EnergyVariables.CcaT = CcaT; EnergyVariables.Ccaa = Ccaa;
    EnergyVariables.airPressureTransportRatio = airPressureTransportRatio;
    EnergyVariables.airPressureTransportScale = airPressureTransportScale;
    EnergyVariables.airHeatIceScale = airHeatIceScale;
    EnergyVariables.sharedEnergyFlux = sharedFlux.energy.total_J_cm2_s;
    EnergyVariables.flux_family_schema_version = sharedFlux.schema_version;
end

function [KTa, VTa, Kcva, Kcaa, Ccaa, ratio, scale, iceScaleMin] = ...
    applyAirPressureHeatImpedance(KTa, VTa, Kcva, Kcaa, Ccaa, KTT, VTT, ...
    SoilVariables, DeltZ, airHeatImpedance)
% Apply the configured ice impedance to pressure-driven heat transport.
% Editor: Jizixing Zhu; Date: 2026-07-14.

NL = numel(DeltZ);
P = SoilVariables.P_gg(:);
DPgDZ = abs((P(2:NL+1) - P(1:NL)) ./ DeltZ);

if isfield(SoilVariables, 'Ratio_ice') && ~isempty(SoilVariables.Ratio_ice)
    iceRatio = max(0, min(1, max(SoilVariables.Ratio_ice, [], 2)));
else
    iceRatio = zeros(NL, 1);
end
iceScaleElem = exp(-airHeatImpedance * iceRatio);
iceScale = [iceScaleElem, iceScaleElem];

KTa = iceScale .* KTa;
VTa = iceScale .* VTa;
Kcva = iceScale .* Kcva;
Kcaa = iceScale .* Kcaa;
Ccaa = iceScale .* Ccaa;

pressureTransport = max(abs(KTa), [], 2) .* DPgDZ ...
    + max(abs(VTa), [], 2) .* max(abs(P(1:NL)), abs(P(2:NL+1)));
thermalTransport = max(abs(KTT), [], 2) ./ max(DeltZ, 1e-12) ...
    + max(abs(VTT), [], 2);

ratio = max(pressureTransport ./ max(thermalTransport, 1e-12));
softRatio = 0.1;
scale = min(1, softRatio / max(ratio, eps));
iceScaleMin = min(iceScaleElem);

if scale < 1
    KTa = scale * KTa;
    VTa = scale * VTa;
    Kcva = scale * Kcva;
    Kcaa = scale * Kcaa;
    Ccaa = scale * Ccaa;
end
end
