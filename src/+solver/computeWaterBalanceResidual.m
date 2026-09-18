function wb = computeWaterBalanceResidual( ...
        waterFlux, storageBeforeMm, storageAfterMm, rootSinkCmS, dtStep)
%COMPUTEWATERBALANCERESIDUAL  Close the soil-box total-water balance.
%
%   wb = solver.computeWaterBalanceResidual(waterFlux, storageBeforeMm,
%   storageAfterMm, rootSinkCmS, dtStep) uses the actual Richards
%   boundary fluxes and the exact depth-integrated internal root sink.
%
%   Inputs:
%     waterFlux      [in] struct of accepted mean rates [cm/s], including
%                    QMT (top, positive outward) and QMB (bottom, negative
%                    for downward outward drainage under the kernel sign).
%     storageBeforeMm [in] [mm] accepted column storage at step start.
%     storageAfterMm  [in] [mm] accepted column storage at step end.
%     rootSinkCmS     [in] [cm/s] integrated Richards internal sink.
%     dtStep          [in] [s] accepted outer-step duration.
%
%   Returns:
%     wb [out] struct [mm] with signed actual top input, signed bottom loss,
%        root sink, soil-storage change, and residual.
%
%   Method:
%     residual = (-QMT)*dt - root_sink*dt - (-QMB)*dt - DeltaStorage.
%     The top term is the boundary-contract flux for every boundary type:
%     exact request data for prescribed Neumann rows and the solved reaction
%     for state-dependent rows. Therefore a prescribed-head experiment does
%     not incorrectly substitute rainfall availability, while a fixed-flux
%     experiment independently tests storage against its imposed input.
%
%   References:
%     soil.water.calculateCoupledBoundaryFlux (accepted h/T/P_g/storage
%       reaction for state-dependent boundaries).
%     soil.water.resolveCandidateBoundaryFlux (exact prescribed Neumann
%       values retained before conservation scoring).
%     ClimaLand energy_hydrology.jl:360-383 applies actual top and bottom
%     boundary fluxes directly to the prognostic water tendency.
%   Editor: Jizixing Zhu; Date: 2026-07-08.
%   Change: 2026-07-12 Jizixing Zhu -- replace the legacy infiltration-minus-
%           evaporation ruler with actual top/bottom flux and exact Srt sink.
%   Change: 2026-07-17 Jizixing Zhu -- remove obsolete potential/actual
%           infiltration aliases; QMT is the sole realized top-interface flux.
%   Change: 2026-07-17 Jizixing Zhu -- point the ledger at the complete coupled
%           water-boundary reaction rather than the legacy h-only recovery.
%   Change: 2026-08-09 Jizixing Zhu -- distinguish prescribed Neumann data
%           from recovered state-dependent boundary reactions.

    requireScalarFinite(dtStep, 'dtStep');
    assert(dtStep > 0, 'computeWaterBalanceResidual:dt', ...
        'dtStep must be > 0 [s].');
    assert(isfield(waterFlux, 'QMT') && isfield(waterFlux, 'QMB'), ...
        'computeWaterBalanceResidual:flux', ...
        'waterFlux must contain the realized QMT and QMB boundary fluxes.');

    wb = struct();
    wb.soil_top_input_mm = -waterFlux.QMT * dtStep * 10;             % [mm]
    wb.soil_bottom_loss_mm = -waterFlux.QMB * dtStep * 10;           % [mm]
    wb.root_sink_mm = rootSinkCmS * dtStep * 10;                     % [mm]
    wb.soil_storage_change_mm = storageAfterMm - storageBeforeMm;    % [mm]
    wb.soil_storage_before_mm = storageBeforeMm;
    wb.soil_storage_after_mm = storageAfterMm;
    wb.residual_mm = wb.soil_top_input_mm - wb.root_sink_mm ...
        - wb.soil_bottom_loss_mm - wb.soil_storage_change_mm;
end

function requireScalarFinite(value, name)
% Validate scalar numerical inputs at the accounting boundary.
    assert(isnumeric(value) && isscalar(value) && isfinite(value), ...
        'computeWaterBalanceResidual:scalar', '%s must be finite scalar.', name);
end
