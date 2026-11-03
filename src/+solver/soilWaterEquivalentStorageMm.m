function mm = soilWaterEquivalentStorageMm(Theta_UU, Theta_V, RHOV, DeltZ)
%SOILWATEREQUIVALENTSTORAGEMM  Integrate conserved column water storage.
%
%   mm = solver.soilWaterEquivalentStorageMm(Theta_UU, Theta_V, RHOV,
%   DeltZ) integrates total condensed plus vapor water with the same endpoint
%   quadrature used by the mixed Richards time term.
%
%   Inputs:
%     Theta_UU [in] [-, NL x 2] total condensed water-equivalent content.
%     Theta_V  [in] [-, NL x 2] gas-filled pore volume.
%     RHOV     [in] [kg/m3, NN x 1] node vapor density.
%     DeltZ    [in] [cm, NL x 1] element thickness, bottom-up.
%
%   Returns:
%     mm [out] [mm] depth-integrated total water-equivalent storage.
%
%   Method:
%     W = sum_e dz_e*(theta_eq,left + theta_eq,right)/2*10.
%     Theta_UU already contains liquid, ice-equivalent, and Ss*h storage.
%
%   Editor: Jizixing Zhu; Date: 2026-07-10.
%   Change: 2026-07-12 Jizixing Zhu -- use Theta_UU directly instead of rebuilding
%           condensed storage from liquid and ice phase volumes.
%   Change: 2026-07-12 Jizixing Zhu -- keep endpoint conversion local instead
%           of exposing a separate public water utility.

    NL = size(Theta_UU, 1);
    assert(numel(DeltZ) == NL && all(isfinite(DeltZ(:))) ...
        && all(DeltZ(:) > 0), ...
        'soilWaterEquivalentStorageMm:dz', ...
        'DeltZ must contain NL finite positive values.');
    thetaEq = waterEquivalentEndpoints(Theta_UU, Theta_V, RHOV);
    mm = sum(0.5 * (thetaEq(:, 1) + thetaEq(:, 2)) .* DeltZ(:)) * 10;
end

function thetaEq = waterEquivalentEndpoints(Theta_UU, Theta_V, RHOV)
% Convert condensed and vapor storage to endpoint total-water equivalent.
% Editor: Jizixing Zhu; Date: 2026-07-12.
    NL = size(Theta_UU, 1);
    RHOV = RHOV(:);
    assert(isequal(size(Theta_UU), [NL, 2]) ...
        && isequal(size(Theta_V), [NL, 2]) && numel(RHOV) == NL + 1, ...
        'soilWaterEquivalentStorageMm:shape', ...
        'Theta_UU/Theta_V must be NL-by-2 and RHOV must contain NL+1 nodes.');
    RHOV_endpoints = [RHOV(1:NL), RHOV(2:NL + 1)];
    thetaEq = Theta_UU ...
        + Theta_V .* RHOV_endpoints / equations.getPhysicalConstant('RHOL');
end
