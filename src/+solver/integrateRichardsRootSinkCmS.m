function sinkCmS = integrateRichardsRootSinkCmS(Srt, dz, firstNode)
%INTEGRATERICHARDSROOTSINKCMS  Integrate root uptake in the Richards domain.
%
%   Inputs:
%     Srt [NL x 2] [s^-1] Volumetric soil-water root sink at endpoints.
%     dz [NL x 1] [cm] Bottom-up finite-element thickness.
%     firstNode [-] First active Richards element/node, bottom-up.
%
%   Returns:
%     sinkCmS [cm s^-1] Positive water removal over firstNode:NL.
%
%   Roots below an external water-table interface are outside the COMPASS
%   unsaturated control volume and must not enter its water ledger.
%
%   Reference: docs/architecture/modflow_external_boundary_contract.md,
%   section "Water quantity".
% Editor: Jizixing Zhu; Date: 2026-07-20.
    NL = numel(dz);
    assert(isequal(size(Srt), [NL, 2]) ...
            && isscalar(firstNode) && isfinite(firstNode) ...
            && firstNode == floor(firstNode) ...
            && firstNode >= 1 && firstNode <= NL, ...
        'solver:integrateRichardsRootSinkCmS:badDomain', ...
        'Srt and firstNode must define elements firstNode:NL.');
    sinkCmS = sum(0.5 * (Srt(firstNode:end, 1) ...
        + Srt(firstNode:end, 2)) .* dz(firstNode:end));
end
