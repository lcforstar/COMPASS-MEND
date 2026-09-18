function residual = calculateLinearWaterResidual(WaterMatrix, RHS, hh, cfg)
%CALCULATELINEARWATERRESIDUAL Residual of assembled COMPASS water equation.
% Part of the +deep_learning/+warmstart package: this lives outside the
% +soil/+water solver so the physics solver has no dependency on the neural
% warm-start path. Soil water solve is the ground truth; this function only
% reads its assembled matrices to score how good a candidate hh profile is.
%
% Call:
%   r = solver.warmstart.calculateLinearWaterResidual(WaterMatrix, RHS, hh, cfg)
%
% Evaluates the same tridiagonal linear system that
% +soil/+water/solveTridiagonalMatrixEquations solves before Thomas elimination:
%
%   C4_a(i-1)*hh(i-1) + C4(i,1)*hh(i) + C4(i,2)*hh(i+1) = RHS(i)
%
% The returned residual is LHS - RHS. It is intended for:
%   - diagnostics during Picard iteration (compare Python PINO loss to MATLAB)
%   - residual gating in the warm-start hook: if r(hh_neural) > r(hh_default)
%     then fall back to the default initial guess.
% It does not alter solver state.

C4   = WaterMatrix.C4;
C4_a = WaterMatrix.C4_a;
hh   = hh(:);
RHS  = RHS(:);

NN = length(hh);
if nargin >= 4 && ~isempty(cfg) && isfield_or_prop(cfg.Grid, 'indxBotmLayer')
    indxBotm = cfg.Grid.indxBotmLayer;
else
    indxBotm = 1;
end

residual = zeros(NN, 1);
for MN = indxBotm:NN
    lhs = C4(MN, 1) * hh(MN);
    if MN < NN
        lhs = lhs + C4(MN, 2) * hh(MN + 1);
    end
    if MN > indxBotm
        lhs = lhs + C4_a(MN - 1) * hh(MN - 1);
    end
    residual(MN) = lhs - RHS(MN);
end
end

function tf = isfield_or_prop(obj, name)
if isstruct(obj)
    tf = isfield(obj, name);
elseif isobject(obj)
    tf = isprop(obj, name);
else
    tf = false;
end
end
