function [QAT, QAB] = calculateBoundaryFlux(SAVE, pressure, cfg)
%CALCULATEBOUNDARYFLUX  Recover solved top and bottom dry-air mass fluxes.
%
%   QAT follows the retained top outward-positive convention; QAB is positive
%   into the soil at the bottom. Both are g cm^-2 s^-1. Public sign conversion
%   belongs to soil.gas.boundary.buildResult.
%   This pressure-only reaction is retained for legacy direct callers. The
%   production h/U/P_g path uses calculateCoupledBoundaryFlux so h/T
%   corrections and exact dry-air storage change are included.
%
%   Editor: Jizixing Zhu; Date: 2026-07-15.

    NN = cfg.Grid.NL + 1;
    pressure = pressure(:);
    assert(numel(pressure) == NN && all(isfinite(pressure)), ...
        'calculateDryAirBoundaryFlux:badPressure', ...
        'pressure must be a finite NN-by-1 vector.');
    QAT = SAVE(2,1,3) - SAVE(2,2,3) * pressure(NN-1) ...
        - SAVE(2,3,3) * pressure(NN);
    QAB = -SAVE(1,1,3) + SAVE(1,2,3) * pressure(1) ...
        + SAVE(1,3,3) * pressure(2);
end
