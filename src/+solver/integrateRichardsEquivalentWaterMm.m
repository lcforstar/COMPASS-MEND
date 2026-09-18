function storageMm = integrateRichardsEquivalentWaterMm(thetaEq, dz, firstNode)
%INTEGRATERICHARDSEQUIVALENTWATERMM  Integrate active-domain water content.
%
%   Inputs:
%     thetaEq [NL x 2] [-] Total liquid/ice plus vapor-equivalent water.
%     dz [NL x 1] [cm] Bottom-up finite-element thickness.
%     firstNode [-] First active Richards element/node, bottom-up.
%
%   Returns:
%     storageMm [mm] Water storage over finite elements firstNode:NL.
%
%   Reference: docs/architecture/modflow_external_boundary_contract.md,
%   section "Water quantity".
% Editor: Jizixing Zhu; Date: 2026-07-20.
    NL = numel(dz);
    assert(isequal(size(thetaEq), [NL, 2]) ...
            && isscalar(firstNode) && isfinite(firstNode) ...
            && firstNode == floor(firstNode) ...
            && firstNode >= 1 && firstNode <= NL, ...
        'solver:integrateRichardsEquivalentWaterMm:badDomain', ...
        'thetaEq and firstNode must define elements firstNode:NL.');
    storageMm = sum(0.5 * (thetaEq(firstNode:end, 1) ...
        + thetaEq(firstNode:end, 2)) .* dz(firstNode:end)) * 10;
end
