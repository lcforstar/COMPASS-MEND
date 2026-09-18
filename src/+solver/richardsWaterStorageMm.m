function storageMm = richardsWaterStorageMm( ...
        Theta_UU, Theta_V, RHOV, dz, firstNode)
%RICHARDSWATERSTORAGEMM  Integrate water over the active Richards elements.
%
%   Inputs:
%     Theta_UU [NL x 2] [-] Condensed-water equivalent at layer endpoints.
%     Theta_V  [NL x 2] [-] Vapor-filled volume at layer endpoints.
%     RHOV     [NN x 1] [kg m^-3] Water-vapor density at nodes.
%     dz       [NL x 1] [cm] Bottom-up finite-element thickness.
%     firstNode [-] First active Richards element/node, bottom-up.
%
%   Returns:
%     storageMm [mm] Liquid/ice plus vapor-equivalent water in firstNode:NL.
%
%   The external saturated segment below firstNode is not COMPASS storage.
%   Vapor mass is converted with the model water density before trapezoidal
%   finite-element integration.
%
%   Reference: docs/architecture/modflow_external_boundary_contract.md,
%   section "Water quantity".
% Editor: Jizixing Zhu; Date: 2026-07-20.
    NL = numel(dz);
    RHOV = RHOV(:);
    assert(isequal(size(Theta_UU), [NL, 2]) ...
            && isequal(size(Theta_V), [NL, 2]) ...
            && numel(RHOV) == NL + 1 ...
            && isscalar(firstNode) && isfinite(firstNode) ...
            && firstNode == floor(firstNode) ...
            && firstNode >= 1 && firstNode <= NL, ...
        'solver:richardsWaterStorageMm:badDomain', ...
        'Water endpoints and firstNode must define elements firstNode:NL.');
    vaporEndpoints = [RHOV(1:NL), RHOV(2:NL + 1)];
    thetaEq = Theta_UU ...
        + Theta_V .* vaporEndpoints / equations.getPhysicalConstant('RHOL');
    storageMm = solver.integrateRichardsEquivalentWaterMm( ...
        thetaEq, dz, firstNode);
end
