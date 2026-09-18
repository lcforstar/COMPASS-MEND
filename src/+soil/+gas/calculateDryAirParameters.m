function AirVariables = calculateDryAirParameters( ...
        SoilVariables, TransportCoefficient, cfg, DSTMAX, sharedFlux)
%CALCULATEDRYAIRPARAMETERS  Build dry-air storage and transport coefficients.
%
%   AirVariables = soil.gas.calculateDryAirParameters(SoilVariables,
%       TransportCoefficient, cfg, DSTMAX) evaluates the coefficients in
%
%     Cah*h_t + CaT*T_t + Caa*P_g,t
%       = d/dz(Kah*h_z + KaT*T_z + Kaa*P_g,z + Cag)
%         + Vah*h_z + VaT*T_z + Vaa*P_g,z.
%
%   Dry air is stored both in gas-filled pores and as dissolved air:
%
%     B_a = rho_a*(theta_g + Hc*theta_l)
%         = rho_a*(porosity + (Hc-1)*theta_l - theta_i).
%
%   Consequently Cah and CaT must differentiate liquid water and ice. This
%   is essential in frozen soil: using theta_i in B_a while omitting its
%   derivative would make the assembled Jacobian inconsistent with the
%   storage function. Gas pressure uses STEMMUS internal pressure units
%   (Pa*10 = g cm^-1 s^-2); length is cm and time is s. Kg is the gas
%   pressure mobility in V_A=-Kg*dP_g/dz.
%
%   Reference: Zeng (2013), Coupled Dynamics in Soil, Eqs. 4.11-4.13 and
%   4.56, book pp. 67-68 and 79 (PDF pp. 81-82 and 93).
% Editor: Jizixing Zhu; Date: 2026-07-13.
% Change: 2026-07-13 Jizixing Zhu -- rename the gas mobility field to Kg without changing
% the dry-air coefficients.
% Change: 2026-07-13 Jizixing Zhu -- differentiate the complete three-phase
% dry-air storage, including the ice-volume contribution.
% Change: 2026-07-13 Jizixing Zhu -- weight mechanical gas dispersion by
% the gas-filled pore fraction in every dry-air diffusion coefficient.
% Change: 2026-08-14 Jizixing Zhu -- consume
% one shared liquid/gas/component flux family, remove the gas-only pressure
% gradient clip, and make the KLa switch consistent in the pressure column.

    NL    = cfg.Grid.NL;
    nD    = cfg.Grid.nD;
    DeltZ = cfg.Grid.DeltZ(:);

    RHOL    = equations.getPhysicalConstant('RHOL');
    Hc      = equations.getPhysicalConstant('Hc');
    Gamma_w = equations.getPhysicalConstant('Gamma_w');

    hh         = SoilVariables.hh(:);
    hh_frez    = SoilVariables.hh_frez(:);
    TT         = SoilVariables.TT(:);
    P_gg       = SoilVariables.P_gg(:);
    Theta_LL   = SoilVariables.Theta_LL;
    Theta_II   = SoilVariables.Theta_II;
    DTheta_LLh = SoilVariables.DThetaLiquidDh;
    DTheta_LLT = SoilVariables.DThetaLiquidDT;
    DTheta_IIh = SoilVariables.DThetaIceDh;
    DTheta_IIT = SoilVariables.DThetaIceDT;
    POR        = SoilVariables.POR(:);
    KfL_h      = SoilVariables.KfL_h;
    KLa_Switch = SoilVariables.KLa_Switch;
    Xah        = SoilVariables.Xah(:);
    XaT        = SoilVariables.XaT(:);
    Xaa        = SoilVariables.Xaa(:);
    RHODA      = SoilVariables.RHODA(:);

    D_Ta        = TransportCoefficient.D_Ta;
    KL_T        = TransportCoefficient.KL_T;
    D_V         = TransportCoefficient.D_V;
    D_Vg        = TransportCoefficient.D_Vg(:);
    V_A         = TransportCoefficient.V_A(:);
    Kg          = TransportCoefficient.Kg;

    KLhBAR = (KfL_h(:,1) + KfL_h(:,2)) * 0.5;
    KLTBAR = (KL_T(:,1)  + KL_T(:,2))  * 0.5;
    DTDBAR = (D_Ta(:,1)  + D_Ta(:,2))  * 0.5;

    if nargin < 5 || isempty(sharedFlux)
        sharedFlux = soil.coupling.evaluateCoupledElementFlux( ...
            SoilVariables, TransportCoefficient, cfg);
    end
    DDhDZ = sharedFlux.gradients.head_per_cm;
    % h is the actual liquid pressure head used by the water equation.
    % h+h_frez is a retention/phase-equilibrium coordinate, not a second
    % Darcy potential. The shared owner gives every consumer the same grad(h).
    DhDZ = sharedFlux.gradients.liquid_pressure_head_per_cm;
    DTDZ = sharedFlux.gradients.temperature_K_per_cm;
    DPgDZ = sharedFlux.gradients.gas_pressure_per_cm;
    QL = sharedFlux.phase.liquid_cm_s;
    QL_h = sharedFlux.decomposition.liquid_head_gravity_cm_s;
    QL_T = sharedFlux.decomposition.liquid_temperature_cm_s;
    QL_A = sharedFlux.decomposition.liquid_gas_pressure_cm_s;
    V_A = sharedFlux.phase.gas_cm_s;

    D_Vg_n = equations.interpolateElementToNode(D_Vg, DeltZ);
    V_A_n  = equations.interpolateElementToNode(V_A,  DeltZ);
    QL_n   = equations.interpolateElementToNode(QL,   DeltZ);

    Cah = zeros(NL, nD);  CaT = zeros(NL, nD);  Caa = zeros(NL, nD);
    Kah = zeros(NL, nD);  KaT = zeros(NL, nD);  Kaa = zeros(NL, nD);
    Vah = zeros(NL, nD);  VaT = zeros(NL, nD);  Vaa = zeros(NL, nD);
    Cag = zeros(NL, nD);

    Hc_m1    = Hc - 1;
    inv_RHOL = 1 / RHOL;
    liquidPressureSwitch = double(~isempty(KLa_Switch) ...
        && KLa_Switch == 1);

    for ML = 1:NL
        for ND = 1:nD
            MN = ML + ND - 1;

            % Three-phase dry-air storage term. Ice volume is not treated
            % as gas pore space, while dissolved air in liquid water is
            % represented by Henry's coefficient.
            theta_g_node = max(POR(ML) - Theta_LL(ML, ND) - Theta_II(ML, ND), 1e-14);
            stor = theta_g_node + Hc * Theta_LL(ML, ND);
            rhoda = RHODA(MN);

            dStorageVolumeDh = Hc_m1 * DTheta_LLh(ML, ND) ...
                - DTheta_IIh(ML, ND);
            dStorageVolumeDT = Hc_m1 * DTheta_LLT(ML, ND) ...
                - DTheta_IIT(ML, ND);

            Cah(ML,ND) = Xah(MN) * stor + rhoda * dStorageVolumeDh;
            CaT(ML,ND) = XaT(MN) * stor + rhoda * dStorageVolumeDT;
            Caa(ML,ND) = Xaa(MN) * stor;

            % In ponded conditions gas pathways may be saturated; then no
            % artificial positive floor is imposed on mobility or diffusivity.
            % Zeng (2013), Coupled Dynamics in Soil, Eqs. 4.12-4.13,
            % p. 67: dispersive dry-air flux uses theta_g*D_Vg.
            gasDispersion = theta_g_node * D_Vg_n(MN);
            if DSTMAX > 1e-6 && hh(end) > 1e-6
                safe_Kg = Kg(ML, ND);
                safe_diff = max(D_V(ML,ND) + gasDispersion, 0);
            else
                safe_Kg = max(Kg(ML, ND), 1e-15);
                safe_diff = max(D_V(ML,ND) + gasDispersion, 1e-18);
            end

            Kah(ML,ND) = Xah(MN) * safe_diff + Hc * rhoda * KfL_h(ML, ND);
            KaT(ML,ND) = XaT(MN) * safe_diff + Hc * rhoda * (KL_T(ML,ND) + D_Ta(ML,ND));
            Kaa(ML,ND) = Xaa(MN) * safe_diff + rhoda * (safe_Kg ...
                + liquidPressureSwitch * Hc * KfL_h(ML,ND) / Gamma_w);

            Cag(ML,ND) = Hc * rhoda * KfL_h(ML, ND);

            adv = -(V_A_n(MN) + Hc * QL_n(MN) * inv_RHOL);
            Vah(ML,ND) = adv * Xah(MN);
            VaT(ML,ND) = adv * XaT(MN);
            Vaa(ML,ND) = adv * Xaa(MN);

            % Disabled experiment retained for traceability: forcing a unit
            % element Peclet number changes the governing transport law.
            % Kah(ML,ND) = max(Kah(ML,ND), abs(Vah(ML,ND)) * DeltZ(ML) / 2);
            % KaT(ML,ND) = max(KaT(ML,ND), abs(VaT(ML,ND)) * DeltZ(ML) / 2);
            % Kaa(ML,ND) = max(Kaa(ML,ND), abs(Vaa(ML,ND)) * DeltZ(ML) / 2);
        end
    end

    AirVariables.Cah    = Cah;     AirVariables.CaT    = CaT;     AirVariables.Caa    = Caa;
    AirVariables.Kah    = Kah;     AirVariables.KaT    = KaT;     AirVariables.Kaa    = Kaa;
    AirVariables.Vah    = Vah;     AirVariables.VaT    = VaT;     AirVariables.Vaa    = Vaa;
    AirVariables.Cag    = Cag;
    AirVariables.dryAirStorage = ...
        [RHODA(1:NL), RHODA(2:NL+1)] .* ...
        (POR + (Hc - 1) * Theta_LL - Theta_II);
    AirVariables.QL     = QL;      AirVariables.QL_h   = QL_h;
    AirVariables.QL_T   = QL_T;    AirVariables.QL_A   = QL_A;
    AirVariables.QL_a   = QL_A;
    AirVariables.KLhBAR = KLhBAR;  AirVariables.KLTBAR = KLTBAR;  AirVariables.DTDBAR = DTDBAR;
    AirVariables.DhDZ   = DhDZ;    AirVariables.DDhDZ  = DDhDZ;
    AirVariables.retentionShiftedHeadGradient = ...
        sharedFlux.gradients.retention_shifted_head_per_cm;
    AirVariables.DTDZ   = DTDZ;    AirVariables.DPgDZ  = DPgDZ;
    AirVariables.flux_family_schema_version = sharedFlux.schema_version;
end
