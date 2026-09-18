function [RHS, AirMatrices] = calculateBoundaryConditions(BoundaryCondition, AirMatrices, RHS, KT)
%CALCULATEBOUNDARYCONDITIONS  Apply retained dry-air pressure or mass-flux BCs.
%
%   The bottom-up grid has node 1 at the lower boundary and node NN at the
%   surface. The supported mathematical rows are
%
%     P_g(1)  = P_g,bottom,       P_g(NN) = P_g,top       (Dirichlet),
%     R_1    <- R_1  + q_a,bottom^in,
%     R_NN   <- R_NN + q_a,top^in                         (Neumann).
%
%   Public pressure is converted from Pa to g cm^-1 s^-2 before this routine.
%   Public mass flux is positive into the soil. The retained top BCP value is
%   outward-positive, so its row update is `R_NN-BCP`; bottom BCPB is already
%   inward-positive and uses `R_1+BCPB`.
%
%   For a prescribed pressure, eliminating the boundary unknown gives
%
%     R_2     <- R_2     - A(2,1) P_g,bottom,
%     R_(N-1) <- R_(N-1) - A(N-1,N) P_g,top.
%
%   Atmospheric and hydrostatic pressures are frozen in the typed Request.
%   The runtime fallback below is retained only for direct legacy callers.
%   Editor: Jizixing Zhu; Date: 2026-07-15.
%   Change: 2026-07-16 Jizixing Zhu -- document the pressure/mass-flux rows
%           in equation order and make the public sign conversion explicit.

    NN    = length(RHS);
    C6    = AirMatrices.C6;
    C6_a  = AirMatrices.C6_a;

    TopPg = BoundaryCondition.TopPg;
    NBCPB = BoundaryCondition.NBCPB;
    BCPB  = BoundaryCondition.BCPB;
    NBCP  = BoundaryCondition.NBCP;
    BCP   = BoundaryCondition.BCP;

    if BoundaryCondition.GroundwaterCoupling == 1
        NBCPB = 1;
    end

    if numel(TopPg) >= KT && isfinite(TopPg(KT))
        top_pg_now = TopPg(KT);
    else
        top_pg_now = BoundaryCondition.BtmPg;
    end
    if isfield(BoundaryCondition, 'BtmPg') && isfinite(BoundaryCondition.BtmPg)
        BtmPg = BoundaryCondition.BtmPg;
    else
        BtmPg = top_pg_now;
    end
    if isfield(BoundaryCondition, 'UseHydrostaticAirBottom') && BoundaryCondition.UseHydrostaticAirBottom
        depth_cm = 0;
        if isfield(BoundaryCondition, 'Tot_Depth')
            depth_cm = BoundaryCondition.Tot_Depth;
        end
        rho_air = 1.2e-3; % g/cm3
        g_cgs = equations.getPhysicalConstant('g');
        BtmPg = top_pg_now + rho_air * g_cgs * depth_cm;
    end

    if NBCPB == 1
        % P_g(1)=P_g,bottom; eliminate the node-1 contribution from row 2.
        RHS(1)   = BtmPg;
        C6(1,1)  = 1;
        RHS(2)   = RHS(2) - C6_a(1) * RHS(1);
        C6(1,2)  = 0;
        C6_a(1)  = 0;
    elseif NBCPB == 2
        % Bottom BCPB is q_a,bottom^in.
        RHS(1) = RHS(1) + BCPB;
    end

    if NBCP == 1
        % P_g(NN)=P_g,top; eliminate the surface contribution from row NN-1.
        RHS(NN)      = BCP;
        C6(NN,1)     = 1;
        RHS(NN-1)    = RHS(NN-1) - C6(NN-1,2) * RHS(NN);
        C6(NN-1,2)   = 0;
        C6_a(NN-1)   = 0;
    elseif NBCP == 2
        % -BCP is the public q_a,top^in.
        RHS(NN) = RHS(NN) - BCP;
    else
        % Legacy atmospheric-pressure code; typed requests use NBCP=1.
        RHS(NN)  = TopPg(KT);
        C6(NN,1) = 1;
        coupling  = C6(NN-1,2);
        if ~isnan(coupling) && ~isinf(coupling)
            RHS(NN-1) = RHS(NN-1) - coupling * RHS(NN);
        else
            fprintf('Warning: C6(NN-1,2) invalid at KT=%d, skipping correction.\n', KT);
        end
        C6(NN-1,2) = 0;
        C6_a(NN-1) = 0;
    end

    RHS = min(max(RHS, -1e30), 1e30);

    AirMatrices.C6   = C6;
    AirMatrices.C6_a = C6_a;
end
