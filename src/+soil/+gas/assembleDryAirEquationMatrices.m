function AirMatrices = assembleDryAirEquationMatrices(coeff, cfg)
%ASSEMBLEDRYAIREQUATIONMATRICES  Assemble the dry-air FEM matrix bank.
%
%   AirMatrices = soil.gas.assembleDryAirEquationMatrices(coeff, cfg)
%   maps dry-air endpoint coefficients into storage blocks C1-C3 and spatial
%   blocks C4-C6 for h, T, and P_g. It performs one element loop and does not
%   apply time discretization, gas boundary conditions, or the linear solve.
%
%   Inputs:
%     coeff [struct] Dry-air coefficients at both endpoints of every element
%                    (NL x 2):
%                    Cah/CaT/Caa storage derivatives for h/T/P_g
%                    Kah/KaT/Kaa gradient-flux coefficients for h/T/P_g
%                    Vah/VaT/Vaa first-order transport coefficients
%                    Cag         state-independent gas flux coefficient
%     cfg   [struct] Grid.NL [-] and Grid.DeltZ [cm, NL x 1]
%
%   Returns:
%     AirMatrices.C1/C2/C3 [equation-specific] Lumped h/T/P_g storage
%     AirMatrices.C4/C5/C6 [equation-specific] h/T/P_g spatial blocks
%     AirMatrices.C4_a/C5_a/C6_a  Strict lower diagonals
%     AirMatrices.C7 [equation-specific] State-independent gas-flux vector
%
%   Numerical origin:
%     The dry-air conservation equation uses the same two-node Galerkin
%     element as the water equation. Endpoint-lumped storage contributes
%     L_e/2. Each gradient-flux column contributes (K1+K2)/(2*L_e).
%     Exact integration of linearly interpolated first-order transport gives
%     V1/3+V2/6 and V1/6+V2/3. The independent gas flux enters C7 with
%     opposite endpoint signs because the shape-function derivatives are
%     B1=-1/L_e and B2=1/L_e.
%
%   References:
%     Zeng (2013), Coupled Dynamics in Soil, Sect. 4.2.3,
%       Eqs. 4.47-4.53, book pp. 76-78 (PDF pp. 90-92).
%     docs/compass_soil_water_heat_physics.md, Sects. 5.1 and 5.3.
%   Editor: Jizixing Zhu; Date: 2026-07-13.

    NL = cfg.Grid.NL;
    NN = NL + 1;
    dz = cfg.Grid.DeltZ(:);

    C1 = zeros(NN, 2);
    C2 = zeros(NN, 2);
    C3 = zeros(NN, 2);
    C4 = zeros(NN, 2);
    C5 = zeros(NN, 2);
    C6 = zeros(NN, 2);
    C4_a = zeros(NN, 1);
    C5_a = zeros(NN, 1);
    C6_a = zeros(NN, 1);
    C7 = zeros(NN, 1);

    for e = 1:NL
        Le = dz(e);

        % Endpoint-lumped storage for the h, T, and P_g columns.
        C1(e,   1) = C1(e,   1) + coeff.Cah(e, 1) * Le / 2;
        C1(e+1, 1) = C1(e+1, 1) + coeff.Cah(e, 2) * Le / 2;
        C2(e,   1) = C2(e,   1) + coeff.CaT(e, 1) * Le / 2;
        C2(e+1, 1) = C2(e+1, 1) + coeff.CaT(e, 2) * Le / 2;
        C3(e,   1) = C3(e,   1) + coeff.Caa(e, 1) * Le / 2;
        C3(e+1, 1) = C3(e+1, 1) + coeff.Caa(e, 2) * Le / 2;

        % h block: gradient-flux contribution Kah plus first-order Vah.
        C4ARG1 = (coeff.Kah(e, 1) + coeff.Kah(e, 2)) / (2 * Le);
        C4ARG2_1 = coeff.Vah(e, 1) / 3 + coeff.Vah(e, 2) / 6;
        C4ARG2_2 = coeff.Vah(e, 1) / 6 + coeff.Vah(e, 2) / 3;
        C4(e,   1) = C4(e,   1) + C4ARG1 - C4ARG2_1;
        C4(e,   2) = C4(e,   2) - C4ARG1 - C4ARG2_2;
        C4(e+1, 1) = C4(e+1, 1) + C4ARG1 + C4ARG2_2;
        C4_a(e) = -C4ARG1 + C4ARG2_1;

        % T block: gradient-flux contribution KaT plus first-order VaT.
        C5ARG1 = (coeff.KaT(e, 1) + coeff.KaT(e, 2)) / (2 * Le);
        C5ARG2_1 = coeff.VaT(e, 1) / 3 + coeff.VaT(e, 2) / 6;
        C5ARG2_2 = coeff.VaT(e, 1) / 6 + coeff.VaT(e, 2) / 3;
        C5(e,   1) = C5(e,   1) + C5ARG1 - C5ARG2_1;
        C5(e,   2) = C5(e,   2) - C5ARG1 - C5ARG2_2;
        C5(e+1, 1) = C5(e+1, 1) + C5ARG1 + C5ARG2_2;
        C5_a(e) = -C5ARG1 + C5ARG2_1;

        % P_g block: gradient-flux contribution Kaa plus first-order Vaa.
        C6ARG1 = (coeff.Kaa(e, 1) + coeff.Kaa(e, 2)) / (2 * Le);
        C6ARG2_1 = coeff.Vaa(e, 1) / 3 + coeff.Vaa(e, 2) / 6;
        C6ARG2_2 = coeff.Vaa(e, 1) / 6 + coeff.Vaa(e, 2) / 3;
        C6(e,   1) = C6(e,   1) + C6ARG1 - C6ARG2_1;
        C6(e,   2) = C6(e,   2) - C6ARG1 - C6ARG2_2;
        C6(e+1, 1) = C6(e+1, 1) + C6ARG1 + C6ARG2_2;
        C6_a(e) = -C6ARG1 + C6ARG2_1;

        % Independent gas flux: integral(B_i*Cag dz) gives opposite endpoint
        % signs and the endpoint arithmetic mean.
        C7ARG = (coeff.Cag(e, 1) + coeff.Cag(e, 2)) / 2;
        C7(e) = C7(e) - C7ARG;
        C7(e+1) = C7(e+1) + C7ARG;
    end

    AirMatrices.C1 = C1;
    AirMatrices.C2 = C2;
    AirMatrices.C3 = C3;
    AirMatrices.C4 = C4;
    AirMatrices.C5 = C5;
    AirMatrices.C6 = C6;
    AirMatrices.C4_a = C4_a;
    AirMatrices.C5_a = C5_a;
    AirMatrices.C6_a = C6_a;
    AirMatrices.C7 = C7;
end
