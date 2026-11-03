function EnergyMatrices = assembleThermalEquationMatrices(coeff, cfg)
%ASSEMBLETHERMALEQUATIONMATRICES  Assemble the thermal-equation FEM matrix bank.
%
%   EnergyMatrices = soil.thermal.assembleThermalEquationMatrices(coeff, cfg)
%   maps thermal-equation endpoint coefficients into three storage blocks
%   C1-C3 and three spatial blocks C4-C6 for the unknowns h, T, and P_g.
%   Spatial assembly is completed here; time terms, boundary conditions, and
%   the tridiagonal solve remain in the thermal trial solver.
%
%   Inputs:
%     coeff [struct] Thermal-equation coefficients at both endpoints of
%                    every element (NL x 2):
%                    CTh/CTT/CTa storage derivatives for h/T/P_g
%                    KTh/KTT/KTa gradient-flux coefficients for h/T/P_g
%                    VTh/VTT/VTa first-order transport coefficients
%                    CTg         state-independent flux coefficient
%                    rootEnergySink volumetric root-energy sink [J cm^-3 s^-1]
%     cfg   [struct] Grid.NL [-], Grid.DeltZ [cm, NL x 1], and
%                    lSoilairefc [-] pressure-coupling switch
%
%   Returns:
%     EnergyMatrices.C1/C2/C3 [equation-specific] Lumped h/T/P_g storage
%     EnergyMatrices.C4/C5/C6 [equation-specific] h/T/P_g spatial blocks
%     EnergyMatrices.C4_a/C5_a/C6_a  Strict lower diagonals
%     EnergyMatrices.C7 [equation-specific] State-independent flux vector
%     EnergyMatrices.C9 [J cm^-2 s^-1] Consistent root-energy sink vector
%
%   Numerical origin:
%     The same two-node Galerkin element used by the water equation is
%     applied to each thermal-equation column. Endpoint-lumped storage gives
%     L_e/2. The gradient-flux integral gives (K1+K2)/(2*L_e). Exact
%     integration of linearly interpolated V gives the 1/3 and 1/6 weights.
%     C7 is the signed endpoint-mean independent flux. When the dry-air
%     equation is disabled, P_g is not an active unknown, so C3/C6/C6_a must
%     remain zero even if the coefficient builder contains placeholder data.
%     The volumetric root-energy sink uses the consistent source integral
%     L_e*[2 1;1 2]/6 and is not mixed into the independent flux C7.
%
%   References:
%     Zeng (2013), Coupled Dynamics in Soil, Sect. 4.2.3,
%       Eqs. 4.47-4.53, book pp. 76-78 (PDF pp. 90-92).
%     docs/compass_soil_water_heat_physics.md, Sects. 5.1 and 5.4.
%   Editor: Jizixing Zhu; Date: 2026-07-13.
%   Change: 2026-07-13 Jizixing Zhu -- assemble root-carried liquid energy as
%   a volumetric sink instead of treating it as an independent face flux.

    NL = cfg.Grid.NL;
    NN = NL + 1;
    dz = cfg.Grid.DeltZ(:);

    C1 = zeros(NN, 2);
    C2 = zeros(NN, 2);
    C3 = zeros(NN, 2);
    C4 = zeros(NN, 2);
    C5 = zeros(NN, 2);
    C6 = zeros(NN, 2);
    C4_a = zeros(NL, 1);
    C5_a = zeros(NL, 1);
    C6_a = zeros(NL, 1);
    C7 = zeros(NN, 1);
    C9 = zeros(NN, 1);

    for e = 1:NL
        Le = dz(e);

        % Endpoint-lumped storage for the h and T columns.
        C1(e,   1) = C1(e,   1) + coeff.CTh(e, 1) * Le / 2;
        C1(e+1, 1) = C1(e+1, 1) + coeff.CTh(e, 2) * Le / 2;
        C2(e,   1) = C2(e,   1) + coeff.CTT(e, 1) * Le / 2;
        C2(e+1, 1) = C2(e+1, 1) + coeff.CTT(e, 2) * Le / 2;

        % h block: gradient-flux contribution KTh plus first-order VTh.
        C4ARG1 = (coeff.KTh(e, 1) + coeff.KTh(e, 2)) / (2 * Le);
        C4ARG2_1 = coeff.VTh(e, 1) / 3 + coeff.VTh(e, 2) / 6;
        C4ARG2_2 = coeff.VTh(e, 1) / 6 + coeff.VTh(e, 2) / 3;
        C4(e,   1) = C4(e,   1) + C4ARG1 - C4ARG2_1;
        C4(e,   2) = C4(e,   2) - C4ARG1 - C4ARG2_2;
        C4(e+1, 1) = C4(e+1, 1) + C4ARG1 + C4ARG2_2;
        C4_a(e) = -C4ARG1 + C4ARG2_1;

        % T block: gradient-flux contribution KTT plus first-order VTT.
        C5ARG1 = (coeff.KTT(e, 1) + coeff.KTT(e, 2)) / (2 * Le);
        C5ARG2_1 = coeff.VTT(e, 1) / 3 + coeff.VTT(e, 2) / 6;
        C5ARG2_2 = coeff.VTT(e, 1) / 6 + coeff.VTT(e, 2) / 3;
        C5(e,   1) = C5(e,   1) + C5ARG1 - C5ARG2_1;
        C5(e,   2) = C5(e,   2) - C5ARG1 - C5ARG2_2;
        C5(e+1, 1) = C5(e+1, 1) + C5ARG1 + C5ARG2_2;
        C5_a(e) = -C5ARG1 + C5ARG2_1;

        if cfg.lSoilairefc == 1
            % P_g is an active column only in the coupled dry-air scope.
            C3(e,   1) = C3(e,   1) + coeff.CTa(e, 1) * Le / 2;
            C3(e+1, 1) = C3(e+1, 1) + coeff.CTa(e, 2) * Le / 2;

            C6ARG1 = (coeff.KTa(e, 1) + coeff.KTa(e, 2)) / (2 * Le);
            C6ARG2_1 = coeff.VTa(e, 1) / 3 + coeff.VTa(e, 2) / 6;
            C6ARG2_2 = coeff.VTa(e, 1) / 6 + coeff.VTa(e, 2) / 3;
            C6(e,   1) = C6(e,   1) + C6ARG1 - C6ARG2_1;
            C6(e,   2) = C6(e,   2) - C6ARG1 - C6ARG2_2;
            C6(e+1, 1) = C6(e+1, 1) + C6ARG1 + C6ARG2_2;
            C6_a(e) = -C6ARG1 + C6ARG2_1;
        end

        % Independent thermal flux: integral(B_i*CTg dz) gives opposite
        % endpoint signs and the endpoint arithmetic mean.
        C7ARG = (coeff.CTg(e, 1) + coeff.CTg(e, 2)) / 2;
        C7(e) = C7(e) - C7ARG;
        C7(e+1) = C7(e+1) + C7ARG;

        % Consistent element load for liquid energy removed by root uptake.
        C9ARG1 = (2 * coeff.rootEnergySink(e, 1) ...
            + coeff.rootEnergySink(e, 2)) * Le / 6;
        C9ARG2 = (coeff.rootEnergySink(e, 1) ...
            + 2 * coeff.rootEnergySink(e, 2)) * Le / 6;
        C9(e) = C9(e) + C9ARG1;
        C9(e+1) = C9(e+1) + C9ARG2;
    end

    EnergyMatrices.C1 = C1;
    EnergyMatrices.C2 = C2;
    EnergyMatrices.C3 = C3;
    EnergyMatrices.C4 = C4;
    EnergyMatrices.C5 = C5;
    EnergyMatrices.C6 = C6;
    EnergyMatrices.C4_a = C4_a;
    EnergyMatrices.C5_a = C5_a;
    EnergyMatrices.C6_a = C6_a;
    EnergyMatrices.C7 = C7;
    EnergyMatrices.C9 = C9;
end
