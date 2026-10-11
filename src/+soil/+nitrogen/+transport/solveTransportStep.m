function [NitrogenVariables, SoilVariables, NitrogenFlux, FrozenADE] = solveTransportStep( ...
        NitrogenVariables, SoilVariables, NitrogenFlux, NitrogenBoundary, AirVariables, ...
        Theta_Solute_L, Theta_Solute_LL, cfg, PeCr, KT, Delt_t)
%SOLVETRANSPORTSTEP  One outer-step ADE solve for NH4 and NO3 column transport.
%
% Runs the standard 7-step ADE assembly-solve pipeline:
%   1. node-interpolate water flux components
%   2. compute top/bottom solute boundary concentrations
%   3. assemble per-layer transport parameters (R, E, B, F, G) for NH4 and NO3
%   4. compute Peclet/Courant numbers (warn on instability)
%   5. assemble tridiagonal matrix coefficients (C1-C6)
%   6. build right-hand-side and apply BC
%   7. Thomas-algorithm solve + realized boundary flux
%
% Inputs:
%   NitrogenVariables [inout]  struct  static + dynamic N solver state
%   SoilVariables     [inout]  struct  C_NH4_new / C_NO3_new updated in place
%   NitrogenFlux      [inout]  struct  cvTop / cvBot / cv{S}* updated
%   NitrogenBoundary  [in]     struct  per-step boundary fields (cT/cB/NBCN/Infiltration/...)
%   AirVariables      [in]     struct  liquid flux QL + sub-components
%   Theta_Solute_L    [in]     [NL×2]  prior-step liquid content (with pore clamp)
%   Theta_Solute_LL   [in]     [NL×2]  new-step liquid content
%   cfg               [in]     struct  Grid / lUpW / lThmrlefc / lBNF switches
%   PeCr              [in]     [-]     critical Peclet number (= 2)
%   KT                [in]     [-]     outer step index
%   Delt_t            [in]     [s]     outer time step length
%
% Returns:
%   NitrogenVariables  updated CHK2, CHK3, N_SAVE
%   SoilVariables      C_NH4_new, C_NO3_new (nonnegative, mass-conservative)
%   NitrogenFlux       all transport-related flux fields updated
%   FrozenADE          read-only final Picard ADE operator and total-solution
%                      projection map for one-way source tracers.  Existing
%                      three-output callers retain identical total dynamics.
%
% Reference:
%   STEMMUS-N legacy N_sub.m (Zhu Jizixing thesis Ch 3.3)
%   van Genuchten & Wagenet 1989 SSSAJ — solute ADE FE assembly
%
% Extracted as standalone module from solveNitrogenEquations.m on 2026-05-13
% as part of +soil/+nitrogen/ refactor (Skill §4: separate composition from physics).
% Change: 2026-08-30 Codex -- replace non-conservative zero clipping with a
%   trapezoidal-storage-weighted nonnegative projection after the unchanged ADE solve.
% Change: 2026-09-21 Wenjie Liu -- expose the converged total ADE operator,
%   raw solution and conservative projection map for one-way fertilizer-N
%   source tracing without feeding tracer state back into the total solve.

NL = cfg.Grid.NL;
NN = cfg.Grid.NN;
DeltZ = cfg.Grid.DeltZ;

% Step 1: node-interpolate liquid-flux components
QL_nodes   = equations.interpolateElementToNode(AirVariables.QL,   DeltZ);
QL_nodes_h = equations.interpolateElementToNode(AirVariables.QL_h, DeltZ);
QL_nodes_T = equations.interpolateElementToNode(AirVariables.QL_T, DeltZ);
QL_nodes_A = equations.interpolateElementToNode(AirVariables.QL_A, DeltZ);

% Step 2: top/bottom boundary concentrations
[NitrogenFlux.cTop_NH4, NitrogenFlux.cTop_NO3, NitrogenFlux.cBot_NH4, NitrogenFlux.cBot_NO3] = ...
    soil.nitrogen.transport.calculateBoundaryConcentrations( ...
        NitrogenFlux, NitrogenBoundary, QL_nodes, KT, Delt_t);

% Step 3 prelim: pre-allocate FE param scratch
vv  = zeros(NL, 1);
vv1 = zeros(NL, 1);
wc1 = zeros(NL, 1);
wc2 = zeros(NL, 1);
Peclet = zeros(NL, 1);
R1 = zeros(NL, cfg.Grid.nD); R1O = R1; E1 = R1; B1 = R1; F1 = R1; G1 = R1;
R2 = R1; R2O = R1; E2 = R1; B2 = R1; F2 = R1; G2 = R1;

% Step 3: assemble solute parameters (per-layer R, E, B, F, G + flux components)
[R1, R1O, E1, B1, F1, G1, R2, R2O, E2, B2, F2, G2, NitrogenFlux.cvSnit, NitrogenFlux.cvSmin, ...
    NitrogenFlux.cvSden, NitrogenFlux.cvSinkSnh4, NitrogenFlux.cvSinkSno3, NitrogenFlux.cvSvol, ...
    NitrogenFlux.QNH4, NitrogenFlux.QNO3, NitrogenFlux.DNH4DZ, NitrogenFlux.DNO3DZ, ...
    NitrogenFlux.QNH4_h, NitrogenFlux.QNH4_T, NitrogenFlux.QNH4_A, NitrogenFlux.QNH4_convection, ...
    NitrogenFlux.QNH4_convection_h, NitrogenFlux.QNH4_convection_T, NitrogenFlux.QNH4_convection_A, ...
    NitrogenFlux.QNH4_dispersion, NitrogenFlux.QNH4_dispersion_h, NitrogenFlux.QNH4_dispersion_T, ...
    NitrogenFlux.QNH4_dispersion_A, NitrogenFlux.QNO3_h, NitrogenFlux.QNO3_T, NitrogenFlux.QNO3_A, ...
    NitrogenFlux.QNO3_convection, NitrogenFlux.QNO3_convection_h, NitrogenFlux.QNO3_convection_T, ...
    NitrogenFlux.QNO3_convection_A, NitrogenFlux.QNO3_dispersion, NitrogenFlux.QNO3_dispersion_h, ...
    NitrogenFlux.QNO3_dispersion_T, NitrogenFlux.QNO3_dispersion_A] = ...
    soil.nitrogen.transport.assembleSoluteParameters(NL, NitrogenVariables.RHOKG, NitrogenVariables.Dsh, NitrogenVariables.Dsh_h, NitrogenVariables.Dsh_T, ...
    NitrogenVariables.Dsh_A, QL_nodes, QL_nodes_h, QL_nodes_T, QL_nodes_A, Theta_Solute_L, Theta_Solute_LL, ...
    NitrogenVariables.Snit, NitrogenVariables.SVnit, NitrogenVariables.SVden, NitrogenVariables.Smin, NitrogenVariables.Svol, ...
    NitrogenVariables.SinkS, NitrogenVariables.Sden, NitrogenFlux.cvSnit, NitrogenFlux.cvSvol, NitrogenFlux.cvSmin, ...
    NitrogenFlux.cvSden, NitrogenFlux.cvSinkSnh4, NitrogenFlux.cvSinkSno3, DeltZ, R1, R1O, E1, B1, F1, G1, ...
    R2, R2O, E2, B2, F2, G2, NitrogenVariables.Shys, SoilVariables.C_NH4_new, SoilVariables.C_NO3_new, KT, ...
    NitrogenVariables.Sdis, NitrogenVariables.SorgNO3);

% Step 4: Peclet/Courant numbers (NH4 uses Retard, NO3 uses ones)
[~, wc1, ~, PECR1] = soil.nitrogen.transport.pecletCourantNumber( ...
    cfg.lUpW, E1, NL, Theta_Solute_LL, QL_nodes, DeltZ, Delt_t, vv, vv1, wc1, Peclet, NitrogenVariables.Retard);
[~, wc2, ~, PECR2] = soil.nitrogen.transport.pecletCourantNumber( ...
    cfg.lUpW, E2, NL, Theta_Solute_LL, QL_nodes, DeltZ, Delt_t, vv, vv1, wc2, Peclet, ones(NL, 2));

% Step 5: assemble tridiagonal matrix coefficients
[C1, C1O, C2, C2O, C3, C4, C3_a, C4_a, C5, C6, M, K, ~, ~, ~, F] = ...
    soil.nitrogen.transport.assembleMatrixCoefficients(R1, R1O, E1, B1, F1, G1, R2, R2O, E2, B2, F2, G2, ...
                                                       DeltZ, NL, NN, wc2, wc1, Theta_Solute_L, Theta_Solute_LL);

% Step 6: build RHS and apply boundary conditions
[RHS, RHS1, C3, C4, N_SAVE] = soil.nitrogen.transport.assembleEquation( ...
    C1, C1O, C2, C2O, C3, C4, C5, C6, NL, NN, Delt_t, SoilVariables.C_NH4, SoilVariables.C_NO3, M, K, F, NitrogenVariables.N_SAVE);
preBoundaryRhsNH4=RHS;
preBoundaryRhsNO3=RHS1;
preBoundaryMatrixNH4=C3;
preBoundaryMatrixNO3=C4;
preBoundaryLowerNH4=C3_a;
preBoundaryLowerNO3=C4_a;
[RHS, RHS1, C4, C4_a, C3, C3_a, NitrogenFlux.cTop_NO3, NitrogenFlux.cTop_NH4] = ...
    soil.nitrogen.transport.applyBoundaryConditions(RHS1, C3, C3_a, RHS, NitrogenBoundary.NBCN, NitrogenBoundary.NBCNB, KT, C4, NN, C4_a, QL_nodes, ...
    cfg.lThmrlefc, NitrogenBoundary.QMB, NitrogenFlux.cBot_NH4, NitrogenFlux.cTop_NH4, NitrogenFlux.cBot_NO3, NitrogenFlux.cTop_NO3, ...
    cfg.lBNF, NitrogenBoundary.BNF);

% Step 7: Thomas solve + non-negativity clip + realized BC flux
[rawNH4, rawNO3, ~, ~, NitrogenVariables.CHK2, NitrogenVariables.CHK3] = ...
    soil.nitrogen.transport.solveSoluteLinearSystem(C4, C3, SoilVariables.C_NH4_new, SoilVariables.C_NO3_new, NN, NL, C4_a, C3_a, RHS, RHS1);
try
    [SoilVariables.C_NH4_new, NitrogenVariables.nh4_positivity] = ...
        soil.nitrogen.transport.projectNonnegativeConservative( ...
            rawNH4,Theta_Solute_LL,NitrogenVariables.RHOKG,DeltZ);
catch err
    if strcmp(err.identifier,'soil:nitrogen:materialNegativeConcentration')
        cause=MException('soil:nitrogen:negativeNH4', ...
            'NH4 positivity failure at KT=%d and transport dt=%.17g s.',KT,Delt_t);
        err=addCause(err,cause);
    end
    rethrow(err);
end
try
    [SoilVariables.C_NO3_new, NitrogenVariables.no3_positivity] = ...
        soil.nitrogen.transport.projectNonnegativeConservative( ...
            rawNO3,Theta_Solute_LL,zeros(NL,1),DeltZ);
catch err
    if strcmp(err.identifier,'soil:nitrogen:materialNegativeConcentration')
        cause=MException('soil:nitrogen:negativeNO3', ...
            'NO3 positivity failure at KT=%d and transport dt=%.17g s.',KT,Delt_t);
        err=addCause(err,cause);
    end
    rethrow(err);
end
NitrogenVariables.N_SAVE = N_SAVE;

% Freeze only data already used by the accepted total solve.  Three-output
% production callers do not pay this diagnostic/export cost.
if nargout>=4
    FrozenADE=soil.nitrogen.transport.buildFrozenADEOperator( ...
        C3,C3_a,C4,C4_a,C1O(:,1),C2O(:,1), ...
        preBoundaryMatrixNH4,preBoundaryLowerNH4, ...
        preBoundaryMatrixNO3,preBoundaryLowerNO3, ...
        preBoundaryRhsNH4,preBoundaryRhsNO3, ...
        SoilVariables.C_NH4,SoilVariables.C_NO3,rawNH4,rawNO3, ...
        SoilVariables.C_NH4_new,SoilVariables.C_NO3_new, ...
        NitrogenVariables.nh4_positivity,NitrogenVariables.no3_positivity, ...
        R1,R1O,E1,B1,F1,R2,R2O,E2,B2,F2,QL_nodes, ...
        Theta_Solute_LL,NitrogenVariables.RHOKG,DeltZ,Delt_t, ...
        NitrogenBoundary,cfg,NitrogenVariables,KT);
end

% Step 7b: 反应通量诊断, 跟 ADE 实际扣除自洽
% [2026-05-16] nit/denit 均改显式源汇 (G1 = -Snit, G2 = -Sden), 账本直接用 NitrogenVariables.*, 不重算
Snit_post = NitrogenVariables.Snit;                                  % [μg/(cm³·s)]
Sden_post = NitrogenVariables.Sden;                                  % [μg/(cm³·s)]

Sden_post_avg = (Sden_post(:,1) + Sden_post(:,2)) * 0.5;
Snit_post_avg = (Snit_post(:,1) + Snit_post(:,2)) * 0.5;
NitrogenFlux.cvSden(KT) = DeltZ(:)' * Sden_post_avg;                 % 跟 ADE 用的 G2 = -Sden 一致
NitrogenFlux.cvSnit(KT) = DeltZ(:)' * Snit_post_avg;                 % 跟 ADE 用的 G1 = -Snit 一致
% NitrogenVariables.Snit/Sden/nit_mm_factor/nit_M_total 已在 nitrification.m / denitrification.m 填好, 不动
if all(PECR1 > PeCr) || all(PECR2 > PeCr)
    warning('soil:nitrogen:soluteOscillation', ...
        'May be oscillatory; consider decreasing Delt_t or DeltZ');
end
[NitrogenFlux.cvBot, NitrogenFlux.cvTop, NitrogenFlux.cvBot1, NitrogenFlux.cvTop1] = ...
    soil.nitrogen.transport.calculateBoundaryFlux(N_SAVE, SoilVariables.C_NH4_new, SoilVariables.C_NO3_new, NN, KT, NitrogenFlux.cvTop, NitrogenFlux.cvBot, ...
    NitrogenFlux.cvTop1, NitrogenFlux.cvBot1, QL_nodes, NitrogenFlux.cBot_NH4, NitrogenFlux.cTop_NH4, NitrogenFlux.cTop_NO3, ...
    NitrogenFlux.cBot_NO3, NitrogenBoundary.NBCN, NitrogenBoundary.NBCNB, 1e-30, cfg.lThmrlefc, NitrogenBoundary.QMB);

end
