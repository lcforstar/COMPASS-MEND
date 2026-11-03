function HBoundaryFlux = calculatesSoilWaterFluxes(SAVE, SoilVariables, bottomNode)
% calculatesSoilWaterFluxes  计算顶底边界水通量（利用 SAVE 快照）
%   Compute top/bottom boundary water fluxes from SAVE snapshot.
%   Editor: Jizixing Zhu; Date: 2026-07-20.
%   Change: 2026-07-20 Jizixing Zhu -- recover the bottom reaction at a MODFLOW
%           internal Richards interface when bottomNode is supplied.

hh = SoilVariables.hh(:);
NN = length(hh);
if nargin < 3 || isempty(bottomNode)
    bottomNode = 1;
end
assert(isscalar(bottomNode) && isfinite(bottomNode) ...
        && bottomNode == floor(bottomNode) ...
        && bottomNode >= 1 && bottomNode < NN, ...
    'calculatesSoilWaterFluxes:badBottomNode', ...
    'bottomNode must be an integer in 1:NN-1.');

HBoundaryFlux.QMT =  SAVE(2,1,1) - SAVE(2,2,1)*hh(NN-1) - SAVE(2,3,1)*hh(NN);
HBoundaryFlux.QMB = -SAVE(1,1,1) + SAVE(1,2,1)*hh(bottomNode) ...
    + SAVE(1,3,1)*hh(bottomNode + 1);
end
