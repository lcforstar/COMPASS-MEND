function built = commitOuterStepLedgers(built, soil_advance, storage_before_mm, dt_s)
%COMMITOUTERSTEPLEDGERS  Finalize conservation ledgers after outer acceptance.
%
%   Inputs:
%     built             [inout] Built model after prognostic commit.
%     soil_advance      [in] Soil scheduler result and accepted accumulators.
%     storage_before_mm [mm] Soil water storage before the outer solve.
%     dt_s              [s] Accepted outer-step duration.
%
%   Return:
%     built [out] Built model with accepted water, heat, dry-air, boundary,
%           and water-closure ledgers.
%
%   Boundary signs follow the existing inward-positive public convention.
%   A water closure is defined only when a real accepted water solve supplied
%   positive accumulation weight.
%   Editor: Jizixing Zhu; Date: 2026-07-16.
%   Change: 2026-07-16 Jizixing Zhu -- commit accepted layer flux
%           diagnostics for write-only solver-internal output.
%   Change: 2026-07-17 Jizixing Zhu -- commit the accepted finite-ponding
%           complementarity summary for output-only diagnostics.
%   Change: 2026-07-17 Jizixing Zhu -- commit a non-double-counted combined
%           surface-reservoir plus soil-column water ledger.
%   Change: 2026-07-18 Jizixing Zhu -- commit the exact no-snow,
%           zero-storage surface-plus-soil energy ledger without changing
%           any accepted state or boundary condition.
%   Change: 2026-07-19 Jizixing Zhu -- replace the restricted prototype ledger by
%           one surface-snow-soil storage and external-energy ledger.
%   Change: 2026-07-20 Jizixing Zhu — use the MODFLOW Richards active domain
%           consistently for storage and root-sink accounting.
    water_accum = soil_advance.water_flux_accum;
    heat_accum = soil_advance.heat_flux_accum;
    energy_balance = soil_advance.energy_balance;
    dry_air_balance = soil_advance.dry_air_mass_balance;
    boundary_summary = soil_advance.boundary_flux_summary;
    built.workspace.solver_diagnostics = soil_advance.solver_diagnostics;

    built.workspace.flux.water_hourly = ...
        solver.finalizeWaterFluxAccum(water_accum);
    if isfield(soil_advance, 'surface_boundary_summary') ...
            && soil_advance.surface_boundary_summary.has_data
        built.workspace.flux.surface_boundary_hourly = ...
            soil_advance.surface_boundary_summary;
    end
    if heat_accum.weight_sum > 0
        built.workspace.flux.heat_hourly = struct( ...
            'QET', heat_accum.QET / heat_accum.weight_sum, ...
            'QEB', heat_accum.QEB / heat_accum.weight_sum);
    end
    if energy_balance.has_data
        built.workspace.flux.soil_energy_balance = energy_balance;
    end
    if dry_air_balance.has_data
        built.workspace.flux.soil_dry_air_mass_balance = dry_air_balance;
    end
    if boundary_summary.has_data
        built.workspace.flux.soil_boundary_hourly = ...
            soil.boundary.finalizeFluxSummary(boundary_summary);
    end
    built.workspace.flux.surface_soil_energy_balance = ...
        surface.energy.buildSurfaceSoilEnergyBalance( ...
            built, soil_advance, dt_s);
    if water_accum.weight_sum <= 0
        return;
    end

    derived = built.workspace.derived;
    firstNode = solver.activeRichardsFirstNode(built);
    storage_after_mm = solver.richardsWaterStorageMm( ...
        derived.Theta_UU, derived.Theta_V, derived.RHOV, ...
        built.domain.dz, firstNode);
    root_sink_cm_s = solver.integrateRichardsRootSinkCmS( ...
        built.workspace.flux.Srt, built.domain.dz, firstNode);
    built.workspace.flux.water_balance = ...
        solver.computeWaterBalanceResidual( ...
        built.workspace.flux.water_hourly, storage_before_mm, ...
        storage_after_mm, root_sink_cm_s, dt_s);
    if isfield(built.workspace.flux, 'surface_boundary_hourly') ...
            && built.workspace.flux.surface_boundary_hourly.has_data
        built.workspace.flux.surface_soil_water_balance = ...
            solver.computeSurfaceSoilWaterBalance( ...
                built.workspace.flux.surface_boundary_hourly, ...
                built.workspace.flux.water_balance);
    end
end
