function WaterMatrix = assembleWaterEquationMatrices(coeff, rootSink, cfg)
%ASSEMBLEWATEREQUATIONMATRICES  Assemble the water-equation FEM matrix bank.
%
%   WaterMatrix = soil.water.assembleWaterEquationMatrices(coeff, rootSink, cfg)
%   maps endpoint coefficients of the coupled water equation into the C1-C9
%   tridiagonal storage contract. This function performs spatial assembly
%   only. Time discretization, boundary conditions, and the linear solve are
%   applied later.
%
%   Inputs:
%     coeff    [struct] Water-equation coefficients at both endpoints of
%                       every element (NL x 2):
%                       Chh/ChT     storage derivatives for h/T
%                       Khh/KhT/Kha gradient-flux coefficients for h/T/P_g
%                       Vvh/VvT     h/T derivatives of the pressure-driven
%                                   vapor flux, used only in correction rows
%                       Chg         state-independent flux coefficient
%     rootSink [1/s, NL x 2] Distributed root-water sink at element endpoints
%     cfg      [struct] Grid settings: NL [-], NN [-], DeltZ [cm, NL x 1]
%
%   Returns:
%     WaterMatrix.C1/C2 [equation-specific] Lumped storage blocks for h/T
%     WaterMatrix.C4/C5 [equation-specific] Spatial blocks for h/T
%     WaterMatrix.C6    [equation-specific] Spatial block for P_g
%     WaterMatrix.C4_a/C5_a       Strict lower diagonals of C4/C5
%     WaterMatrix.C4_vapor_jacobian/C5_vapor_jacobian and *_a
%                       Jacobian-only h/T blocks for vapor advection
%     WaterMatrix.C7    [equation-specific] State-independent flux vector
%     WaterMatrix.C9    [equation-specific] Consistent distributed-sink vector
%
%   Numerical origin:
%     On element e=[z_e,z_{e+1}], N1=1-x/L_e and N2=x/L_e, with
%     B1=dN1/dz=-1/L_e and B2=dN2/dz=1/L_e. C1/C2 use endpoint-lumped
%     storage, L_e/2 at each endpoint. C4/C5/C6 come from the weak-form
%     gradient operator. The physical residual keeps the pressure-driven
%     vapor flux in conservative C6 form. Vvh/VvT use the exact weights
%     V1/3+V2/6 and V1/6+V2/3 in separate correction-Jacobian banks; they
%     are never applied to the absolute h/T state. C7 integrates B_i times
%     the endpoint-mean independent flux. C9 uses the consistent source
%     integral L_e*[2 1;1 2]/6.
%
%   References:
%     Zeng (2013), Coupled Dynamics in Soil, Sect. 4.2.3,
%       Eqs. 4.47-4.53, book pp. 76-78 (PDF pp. 90-92).
%     docs/compass_soil_water_heat_physics.md, Eqs. 5.W6-5.W12.
%   Editor: Jizixing Zhu; Date: 2026-07-13.
%   Change: 2026-08-12 Jizixing Zhu -- keep the rho_v*Kg*grad(P_g) term
%           conservative in the physical residual and expose Vvh/VvT only
%           as correction-Jacobian banks, preventing product-rule overlap.

    NL = cfg.Grid.NL;
    NN = cfg.Grid.NN;
    dz = cfg.Grid.DeltZ;

    Chh = coeff.Chh;
    ChT = coeff.ChT;
    Khh = coeff.Khh;
    KhT = coeff.KhT;
    Kha = coeff.Kha;
    Vvh = coeff.Vvh;
    VvT = coeff.VvT;
    Chg = coeff.Chg;

    C1 = zeros(NN, 2);
    C2 = zeros(NN, 2);
    C4 = zeros(NN, 2);
    C4_a = zeros(NN, 1);
    C5 = zeros(NN, 2);
    C5_a = zeros(NN, 1);
    C4_vapor_jacobian = zeros(NN, 2);
    C4_vapor_jacobian_a = zeros(NN, 1);
    C5_vapor_jacobian = zeros(NN, 2);
    C5_vapor_jacobian_a = zeros(NN, 1);
    C6 = zeros(NN, 2);
    C7 = zeros(NN, 1);
    C9 = zeros(NN, 1);

    for e = 1:NL
        Le = dz(e);

        % Endpoint-lumped storage: endpoint Lobatto quadrature gives Le/2
        % on the main diagonal and zero storage off-diagonals.
        C1(e,   1) = C1(e,   1) + Chh(e, 1) * Le / 2;
        C1(e+1, 1) = C1(e+1, 1) + Chh(e, 2) * Le / 2;
        C2(e,   1) = C2(e,   1) + ChT(e, 1) * Le / 2;
        C2(e+1, 1) = C2(e+1, 1) + ChT(e, 2) * Le / 2;

        % Physical h-gradient block. The pressure-driven vapor flux is
        % already conservative in C6, so Vvh must not enter this residual.
        C4ARG1 = (Khh(e, 1) + Khh(e, 2)) / (2 * Le);
        C4(e,   1) = C4(e,   1) + C4ARG1;
        C4(e,   2) = C4(e,   2) - C4ARG1;
        C4(e+1, 1) = C4(e+1, 1) + C4ARG1;
        C4_a(e) = -C4ARG1;

        % Partial Jacobian of rho_v*Kg*grad(P_g) with respect to h.
        % This multiplies a correction delta_h, never the absolute h state.
        C4ARG2_1 = Vvh(e, 1) / 3 + Vvh(e, 2) / 6;
        C4ARG2_2 = Vvh(e, 1) / 6 + Vvh(e, 2) / 3;
        C4_vapor_jacobian(e,   1) = C4_vapor_jacobian(e,   1) - C4ARG2_1;
        C4_vapor_jacobian(e,   2) = C4_vapor_jacobian(e,   2) - C4ARG2_2;
        C4_vapor_jacobian(e+1, 1) = C4_vapor_jacobian(e+1, 1) + C4ARG2_2;
        C4_vapor_jacobian_a(e) = C4ARG2_1;

        % Physical T-gradient block and its separate vapor-flux Jacobian.
        C5ARG1 = (KhT(e, 1) + KhT(e, 2)) / (2 * Le);
        C5(e,   1) = C5(e,   1) + C5ARG1;
        C5(e,   2) = C5(e,   2) - C5ARG1;
        C5(e+1, 1) = C5(e+1, 1) + C5ARG1;
        C5_a(e) = -C5ARG1;
        C5ARG2_1 = VvT(e, 1) / 3 + VvT(e, 2) / 6;
        C5ARG2_2 = VvT(e, 1) / 6 + VvT(e, 2) / 3;
        C5_vapor_jacobian(e,   1) = C5_vapor_jacobian(e,   1) - C5ARG2_1;
        C5_vapor_jacobian(e,   2) = C5_vapor_jacobian(e,   2) - C5ARG2_2;
        C5_vapor_jacobian(e+1, 1) = C5_vapor_jacobian(e+1, 1) + C5ARG2_2;
        C5_vapor_jacobian_a(e) = C5ARG2_1;

        % P_g block: Kha multiplies dP_g/dz inside the conservative flux;
        % no separate pressure first-order coefficient exists in this bank.
        C6ARG = (Kha(e, 1) + Kha(e, 2)) / (2 * Le);
        C6(e,   1) = C6(e,   1) + C6ARG;
        C6(e,   2) = C6(e,   2) - C6ARG;
        C6(e+1, 1) = C6(e+1, 1) + C6ARG;

        % Independent flux: integral(B_i*Chg dz) gives opposite endpoint
        % signs and the arithmetic mean of the endpoint coefficients.
        C7ARG = (Chg(e, 1) + Chg(e, 2)) / 2;
        C7(e) = C7(e) - C7ARG;
        C7(e+1) = C7(e+1) + C7ARG;

        % Distributed sink: integral(N_i*rootSink dz) gives the consistent
        % two-node source matrix Le*[2 1;1 2]/6.
        C9ARG1 = (2 * rootSink(e, 1) + rootSink(e, 2)) * Le / 6;
        C9ARG2 = (rootSink(e, 1) + 2 * rootSink(e, 2)) * Le / 6;
        C9(e) = C9(e) + C9ARG1;
        C9(e+1) = C9(e+1) + C9ARG2;
    end

    WaterMatrix.C1 = C1;
    WaterMatrix.C2 = C2;
    WaterMatrix.C4 = C4;
    WaterMatrix.C4_a = C4_a;
    WaterMatrix.C5 = C5;
    WaterMatrix.C5_a = C5_a;
    WaterMatrix.C4_vapor_jacobian = C4_vapor_jacobian;
    WaterMatrix.C4_vapor_jacobian_a = C4_vapor_jacobian_a;
    WaterMatrix.C5_vapor_jacobian = C5_vapor_jacobian;
    WaterMatrix.C5_vapor_jacobian_a = C5_vapor_jacobian_a;
    WaterMatrix.C6 = C6;
    WaterMatrix.C7 = C7;
    WaterMatrix.C9 = C9;
end
