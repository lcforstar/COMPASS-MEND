function [AirMatrices, P_gg, RHS] = solveTridiagonalMatrixEquations(RHS, AirMatrices, cfg)
% solveTridiagonalMatrixEquations Solve dry-air tridiagonal system.

    NL     = cfg.Grid.NL;
    NN     = cfg.Grid.NN;
    my_eps = 1e-20;
    C6     = AirMatrices.C6;
    C6_a   = AirMatrices.C6_a;

    val = C6(1,1);
    s = sign(val); if s == 0, s = 1; end
    C6(1,1) = max(abs(val), my_eps) * s;
    RHS(1)  = RHS(1) / C6(1,1);

    for i = 2:NN
        val_prev = C6(i-1, 1);
        s = sign(val_prev); if s == 0, s = 1; end
        denom    = max(abs(val_prev), my_eps) * s;
        C6(i,1)  = C6(i,1) - C6_a(i-1) * C6(i-1,2) / denom;
        val_curr = C6(i,1);
        s = sign(val_curr); if s == 0, s = 1; end
        C6(i,1) = max(abs(val_curr), my_eps) * s;
        RHS(i)  = (RHS(i) - C6_a(i-1) * RHS(i-1)) / C6(i,1);
    end

    for i = NL:-1:1
        RHS(i) = RHS(i) - C6(i,2) * RHS(i+1) / C6(i,1);
    end

    P_gg = RHS;
    AirMatrices.C6 = C6;
end
