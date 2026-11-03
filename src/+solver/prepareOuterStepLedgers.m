function [built, storage_before_mm] = prepareOuterStepLedgers(built, dt_s, kt)
%PREPAREOUTERSTEPLEDGERS  Initialize outer-step conservation accumulators.
%
%   Inputs:
%     built [inout] Built model before the soil outer-step scheduler.
%     dt_s  [s] Actual outer-step duration.
%     kt    [-] One-based outer-step index.
%
%   Returns:
%     built             [out] Built model with empty water, heat, dry-air,
%                       and boundary ledgers.
%     storage_before_mm [mm] Accepted soil-column water before the solve.
%
%   Rejected steps retain NaN water closure; only accepted tile accumulation
%   produces a defined balance in solver.commitOuterStepLedgers.
%   Editor: Jizixing Zhu; Date: 2026-07-16.
%   Change: 2026-07-17 Jizixing Zhu -- clear the output-only surface
%           complementarity summary before each outer transaction.
%   Change: 2026-07-17 Jizixing Zhu -- clear the combined surface-soil water
%           ledger before each outer transaction.
%   Change: 2026-07-18 Jizixing Zhu -- clear the accepted zero-storage
%           surface-plus-soil energy ledger before each transaction.
%   Change: 2026-07-20 Jizixing Zhu — use the same active Richards control
%           volume as the MODFLOW candidate acceptance gate.
    built.workspace.flux.water_hourly = solver.finalizeWaterFluxAccum( ...
        solver.initWaterFluxAccum());
    built.workspace.flux.surface_boundary_hourly = struct('has_data', false);
    built.workspace.flux.surface_soil_water_balance = ...
        solver.emptySurfaceSoilWaterBalance();
    built.workspace.flux.heat_hourly = struct('QET', 0, 'QEB', 0);
    built.workspace.flux.soil_energy_balance = ...
        soil.thermal.emptyColumnEnergyBalance();
    built.workspace.flux.surface_soil_energy_balance = ...
        struct('has_data', false);
    built.workspace.flux.soil_dry_air_mass_balance = ...
        soil.gas.emptyColumnMassBalance();
    built.workspace.flux.soil_boundary_hourly = ...
        soil.boundary.emptyFluxSummary();
    firstNode = solver.activeRichardsFirstNode(built);
    storage_before_mm = solver.retentionStorageMm( ...
        built, dt_s, kt, firstNode);
    built.workspace.flux.water_balance = solver.emptyWaterBalance();
end
