function [SoilVariables, TransportCoefficient, CHKWATER, ...
    QMT, QMB, AVAIL0, Throughfall, Infiltration, Infiltration_act, ...
    TopFlux_potential, SAVE, runoffHort, LinearSystem] = ...
    solveSoilMoistureBalance(SoilVariables, TransportCoefficient, cfg, BoundaryCondition, ...
    Delt_t, KT, SAVE, ...
    Throughfall, EVAP, ...
    hN, DSTMAX, h_target, NBChh, hcrit, thetaEqN, ...
    precomputedCoefficients, precomputedMatrix)
% solveSoilMoistureBalance  Solve one sub-step of the soil moisture equation (Richards equation).
%
% Caller must write the latest local state into SoilVariables / TransportCoefficient before calling.
% All sub-functions pass parameters through structs; no global variables used.
%   Editor: Jizixing Zhu; Date: 2026-07-10.
%   Change: 2026-07-10 Jizixing Zhu -- forward the water time-term control to calculateTimeDerivatives.
%   Change: 2026-07-11 Jizixing Zhu -- add the conservative phase-storage kernel
%           selected only by the new transactional adapter.
%   Change: 2026-07-12 Jizixing Zhu -- require accepted thetaEqN and use the
%           conservative mixed water equation as the sole time formulation.
%   Change: 2026-07-13 Jizixing Zhu -- accept optional coefficient/common-FEM
%           matrix banks from the new runOneStep chain; legacy callers that
%           omit them retain the original internal calculation path.
%   Change: 2026-07-14 Jizixing Zhu -- optionally return the fully assembled
%           post-boundary linear system without executing its absolute-state
%           solve, for the h/U/P_g correction-coordinate path.
assert(nargin >= 15 && ~isempty(thetaEqN), ...
    'solveSoilMoistureBalance:missingThetaEqN', ...
    'Accepted-state thetaEqN is required by the mixed water equation.');
if cfg.lDEBUG_NAN_TRACE
    fprintf('[water] Before calculateMatricCoefficients: P_gg NaN=%d, hh NaN=%d\n', ...
        sum(isnan(SoilVariables.P_gg(:))), sum(isnan(SoilVariables.hh(:))));
end

% --- 1. Matric coefficients --------------------------------------------------
if nargin >= 16 && ~isempty(precomputedCoefficients)
    MatricCoeff = precomputedCoefficients;
else
    switch lower(cfg.MatricKernel)
        case 'conservative'
            [MatricCoeff, SoilVariables, TransportCoefficient] = ...
                soil.water.calculateWaterEquationCoefficients( ...
                    SoilVariables, TransportCoefficient, cfg);
        case 'array'
            [MatricCoeff, SoilVariables, TransportCoefficient] = ...
                soil.water.calculateMatricCoefficientsArray(SoilVariables, TransportCoefficient, cfg);
        otherwise
            [MatricCoeff, SoilVariables, TransportCoefficient] = ...
                soil.water.calculateMatricCoefficients(SoilVariables, TransportCoefficient, cfg);
    end
end

if cfg.lDEBUG_NAN_TRACE && KT == 1
    fprintf('[water coeff checkpoint]\\n');
    fprintf('  DTheta_LLh range: [%e, %e]\\n', ...
        min(SoilVariables.DTheta_LLh(:)), max(SoilVariables.DTheta_LLh(:)));
    fprintf('  Se range: [%e, %e]\\n', ...
        min(SoilVariables.Se(:)), max(SoilVariables.Se(:)));
    fprintf('  KL_h range: [%e, %e]\\n', ...
        min(SoilVariables.KL_h(:)), max(SoilVariables.KL_h(:)));
    fprintf('  KfL_h range: [%e, %e]\\n', ...
        min(SoilVariables.KfL_h(:)), max(SoilVariables.KfL_h(:)));
end

if cfg.lDEBUG_NAN_TRACE
    fprintf('[water] After calculateMatricCoefficients: Chh NaN=%d, Khh NaN=%d\n', ...
        sum(isnan(MatricCoeff.Chh(:))), sum(isnan(MatricCoeff.Khh(:))));
end

% --- 2. Assemble stiffness matrix --------------------------------------------
% [PATCH 2026-05-08 v16] HYDRUS dual-porosity SinkIm 作为额外汇项进入 Richards RHS
% HYDRUS WATFLOW.f line 223/255/286: P -= SinkIm·dx
% 关键守恒约束：SinkIm 已在 substep 开始时由 dualPorExchange 用步初 (ThOld, ThOldIm)
% 算定并缓存到 SoilVariables.SinkIm_current，整个 Picard 迭代里保持常数。
% 这保证 Richards RHS 扣的 SinkIm·Δt 与 immobile 域加的 SinkIm·Δt 完全相等。
% 单位：SinkIm [1/s] 与 Srt [1/s] 量纲一致，直接相加。
SinkSrt = SoilVariables.Srt;
try
    if cfg.lDualPor == 1 && isfield(SoilVariables, 'SinkIm_current') ...
            && ~isempty(SoilVariables.SinkIm_current)
        SinkSrt = SinkSrt + SoilVariables.SinkIm_current;
    end
catch
end
if nargin >= 17 && ~isempty(precomputedMatrix)
    WaterMatrix = precomputedMatrix;
else
    WaterMatrix = soil.water.assembleWaterEquationMatrices(MatricCoeff, SinkSrt, cfg);
end

if cfg.lDEBUG_NAN_TRACE && KT == 1
    fprintf('  Srt range: [%e, %e]\n', min(SoilVariables.Srt(:)), max(SoilVariables.Srt(:)));
    fprintf('  C9 range: [%e, %e]\n', min(WaterMatrix.C9(:)), max(WaterMatrix.C9(:)));
end

if cfg.lDEBUG_NAN_TRACE
    fprintf('[water] After h_MAT: C4 NaN=%d\n', sum(isnan(WaterMatrix.C4(:))));
end

% --- 3. Construct right-hand side --------------------------------------------
[WaterMatrix, RHS, SAVE] = soil.water.calculateTimeDerivatives( ...
    WaterMatrix, SoilVariables, SoilVariables.P_gg, Delt_t, cfg, SAVE, thetaEqN);

if cfg.lDEBUG_NAN_TRACE && KT == 1
    fprintf('  RHS range: [%e, %e]\n', min(RHS(:)), max(RHS(:)));
end

if cfg.lDEBUG_NAN_TRACE
    fprintf('[water] After h_EQ: RHS NaN=%d\n', sum(isnan(RHS(:))));
end

% --- 4. Apply boundary conditions --------------------------------------------
[WaterMatrix, RHS, AVAIL0, Throughfall, Infiltration, Infiltration_act, TopFlux_potential, runoffHort] = ...
    soil.water.calculateBoundaryConditions(WaterMatrix, RHS, SoilVariables, BoundaryCondition, cfg, ...
    KT, hN, NBChh, h_target, Throughfall, EVAP, Delt_t);

if cfg.lDEBUG_NAN_TRACE
    fprintf('[water] After calculateBoundaryConditions: RHS NaN=%d\n', sum(isnan(RHS(:))));
end

LinearSystem = struct('matrix', WaterMatrix, 'rhs', RHS(:), 'SAVE', SAVE);

% --- 5. Solve tridiagonal system ---------------------------------------------
% 暴雨期（Throughfall>0）放行 hh<hcrit 钉死，让雨水通量自然推 hh 上升。
isRaining = Throughfall(KT) > 1e-10;  % cm/s; 1e-10 ≈ 0.0036 mm/h
assembleOnly = isfield(cfg, 'AssembleOnly') && cfg.AssembleOnly;
if assembleOnly
    CHKWATER = zeros(cfg.Grid.NN, 1);
else
    [SoilVariables, CHKWATER] = ...
        soil.water.solveTridiagonalMatrixEquations(WaterMatrix, RHS, SoilVariables, DSTMAX, hcrit, cfg, isRaining);
end

% Optional hook for neural warm-start training (PFEM-style PINO).
% Lives in +deep_learning/+warmstart so the physics solver does not depend on
% the deep-learning path. Disabled by default; enabled by environment variable.
if ~assembleOnly
    solver.warmstart.saveWaterResidualSnapshotIfRequested( ...
        WaterMatrix, RHS, SoilVariables, CHKWATER, cfg, KT, Delt_t, 'after_solve');
end

if cfg.lDEBUG_NAN_TRACE && KT == 1
    fprintf('  hh range after solve: [%e, %e]\n', ...
        min(SoilVariables.hh(:)), max(SoilVariables.hh(:)));
    fprintf('  Theta_LL range after solve: [%e, %e]\n', ...
        min(SoilVariables.Theta_LL(:)), max(SoilVariables.Theta_LL(:)));
    fprintf('  CHKWATER range: [%e, %e]\n', min(CHKWATER(:)), max(CHKWATER(:)));

end

if cfg.lDEBUG_NAN_TRACE
    fprintf('[water] After solveTridiagonalMatrixEquations: hh NaN=%d\n', sum(isnan(SoilVariables.hh(:))));
end

% --- 6. Compute boundary fluxes ----------------------------------------------
HBoundaryFlux = soil.water.calculatesSoilWaterFluxes(SAVE, SoilVariables);
QMT(KT) = HBoundaryFlux.QMT;
QMB(KT) = HBoundaryFlux.QMB;

% ---- 诊断日志（暴雨期 SAVE/QMT 反算），cfg.lDEBUG_RAIN_BC 控制 ----
if cfg.lDEBUG_RAIN_BC
    rain_thresh = 1e-6;  % cm/s
    if Throughfall(KT) > rain_thresh
        NN = cfg.Grid.NN;
        fprintf(['[qmt-debug KT=%d] hh(NN-1)=%.4e hh(NN)=%.4e ' ...
                 'SAVE(2,1,1)=%.4e SAVE(2,2,1)=%.4e SAVE(2,3,1)=%.4e ' ...
                 'QMT=%.4e QMT*dt*10=%.4f mm  QMB=%.4e\n'], ...
                KT, SoilVariables.hh(NN-1), SoilVariables.hh(NN), ...
                SAVE(2,1,1), SAVE(2,2,1), SAVE(2,3,1), ...
                QMT(KT), -QMT(KT) * Delt_t * 10, QMB(KT));
    end
end

end
