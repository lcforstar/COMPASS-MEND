function [model, candidates] = advanceCanopyWaterSnow( ...
        model, dt_outer_s, snowTransactionOverride)
%ADVANCECANOPYWATERSNOW  Stage canopy water and ClimaLand bulk snow.
%
%   Canopy interception is evaporated first. Ground throughfall is partitioned
%   into rain and snow with the configured ClimaLand forcing adapter.  SCOPE's
%   current-step predicted/corrected snow surface supplies cover, vapor, and
%   snow-soil heat auxiliaries, while prognostic S/S_l/U is always advanced
%   from accepted Y.snow. Snow drainage becomes atmospheric soil-water supply,
%   bare-ground evaporation remains a soil sink, and the SCOPE closure heat
%   flux is replaced by the one grid-mean flux entering the shared soil column.
%
%   All prognostic outputs remain candidates until the complete outer step is
%   accepted. Frozen depression storage is an internal staged transfer: a
%   failed outer step commits neither its removal nor its addition to snow.
%   workspace.snow.surface_state is an auxiliary SCOPE view;
%   workspace.snow.boundary_state is the current-step candidate consumed by
%   soil heat-boundary materialization and is not accepted storage by itself.
%   SCOPE's empirical ground-heat closure remains diagnostic; the soil
%   consumes the physical affine surface-energy response assembled here.
%
%   Method (body order):
%     (1) Canopy interception evaporation from SCOPE fluxes and guarded Tss.
%     (2) Throughfall rain/snow partition; frozen pond stage amount; snow heat
%         and vapor drivers; advance bulk S/S_l/U.
%     (3) Stage the affine soil-top response and compatibility diagnostics.
%     (4) Build canopy_water / water_flux / snow candidates for commitPending.
%
%   Editor: Jizixing Zhu; Date: 2026-07-16.
%   Change: 2026-07-17 Jizixing Zhu -- replace the temperature-step snow
%     bucket with the ClimaLand S/S_l/U model and couple its liquid, vapor,
%     ground-heat, and frozen-pond fluxes to the soil transaction.
%   Change: 2026-07-17 Jizixing Zhu -- consume the explicit SCOPE bare/snow
%     latent and surface-closure partitions without a second area weighting.
%   Change: 2026-07-17 Jizixing Zhu -- add numbered step comments in the body.
%   Change: 2026-07-17 Jizixing Zhu -- use the current-step snow surface seen
%     by SCOPE for snow/soil heat and area weighting, and stage the resulting
%     candidate for the same-step soil thermal boundary.
%   Change: 2026-07-17 Jizixing Zhu -- document separate auxiliary surface,
%     thermal-boundary candidate, and accepted prognostic snow ownership.
%   Change: 2026-07-18 Jizixing Zhu -- route signed snow vapor only to the
%     snow store and bare condensation as nonnegative soil liquid supply.
%   Change: 2026-07-18 Jizixing Zhu -- stage exact bare/snow contributions
%     of the single aggregate energy flux entering the shared soil column.
%   Change: 2026-07-19 Jizixing Zhu -- replace the empirical/frozen soil heat
%     target with one CLM-style affine surface-energy response at z=0.
%   Change: 2026-07-19 Jizixing Zhu -- retain the immutable snow advance drivers so
%     post-soil interface reconciliation can replay the same transition.
%   Change: 2026-07-19 Jizixing Zhu -- accept a pre-reconciled terminal-meltout
%     transaction during the one bounded surface/soil event restage.

    if nargin < 3
        snowTransactionOverride = struct();
    end

    % --- Step 0: local handles (guarded Tss is already on scope_step.thermal) ---
    built = model.built;
    scopeStep = built.workspace.scope_step;
    interception = scopeStep.interception;
    surfaceTemperature = scopeStep.thermal.Tss;  % [degC] SCOPE-GUARD output
    maxWetFraction = readField( ...
        built.parameters.veg.canopy, 'Ps_max_diurnal', 0.7);

    % --- Step 1: evaporate canopy interception (updates h2ocan candidates) ---
    evaporation = vegetation.computeCanopyWaterEvap( ...
        interception, scopeStep.fluxes, surfaceTemperature, built.Y.T, ...
        maxWetFraction, dt_outer_s);

    % --- Step 2a: throughfall -> rain vs snow rates [cm/s] ---
    phase = surface.snow.partitionPrecipitation( ...
        interception.throughfall_cms, scopeStep.meteo.Ta, ...
        scopeStep.meteo.RH, built.parameters.snow);
    snowSurface = built.Y.snow;
    if isfield(built.workspace, 'snow') ...
            && isfield(built.workspace.snow, 'surface_state')
        snowSurface = built.workspace.snow.surface_state;
    end
    surfaceSnow = surface.snow.diagnoseState( ...
        snowSurface, built.parameters.snow);

    % --- Step 2b: stage frozen depression storage -> snow ice [cm] (commit later) ---
    frozenPond = frozenPondTransfer(built, surfaceTemperature);

    % --- Step 2c: snow-soil conductive heat and SCOPE surface-closure partitions ---
    snowGroundResponse = surface.snow.computeGroundHeatResponse( ...
        built, snowSurface, built.parameters.snow);
    snowToSoil = snowGroundResponse.reference_flux_into_soil_W_m2;
    surfaceClosure = readField(scopeStep.fluxes, ...
        'soil_top_downward_energy_flux_W_m2', ...
        readField(scopeStep.fluxes, 'Gtot', 0));
    snowSurfaceClosure = readField( ...
        scopeStep.fluxes, 'Gtot_snow_surface', surfaceClosure);
    bareSurfaceClosure = readField(scopeStep.fluxes, ...
        'bare_ground_external_energy_flux_W_m2', ...
        readField(scopeStep.fluxes, 'Gtot_bare_surface', surfaceClosure));

    % --- Step 2d: snow-pack vapor loss/gain [cm/s water equivalent] ---
    snowVapor = snowVaporRate( ...
        scopeStep, snowSurface, built.parameters.snow);

    % --- Step 2e: pack drivers and advance bulk snow (S, S_l, U) ---
    drivers = struct( ...
        'rainfall_cm_s', phase.liquid_cm_s, ...
        'snowfall_cm_s', phase.snow_cm_s, ...
        'snow_to_atmosphere_vapor_cm_s', snowVapor, ...
        'surface_closure_flux_W_m2', snowSurfaceClosure, ...
        'bare_surface_closure_flux_W_m2', bareSurfaceClosure, ...
        'snow_to_soil_heat_flux_W_m2', snowToSoil, ...
        'surface_cover_fraction', surfaceSnow.cover_fraction, ...
        'air_temperature_C', scopeStep.meteo.Ta, ...
        'frozen_pond_transfer_cm', frozenPond);
    hasOverride = isstruct(snowTransactionOverride) ...
        && isfield(snowTransactionOverride, 'state') ...
        && isfield(snowTransactionOverride, 'flux');
    if hasOverride
        snowNext = snowTransactionOverride.state;
        snowFlux = snowTransactionOverride.flux;
        assert(isstruct(snowNext) && isstruct(snowFlux) ...
                && isfield(snowFlux, 'terminal_meltout_triggered') ...
                && logical(snowFlux.terminal_meltout_triggered), ...
            'advanceCanopyWaterSnow:badSnowOverride', ...
            'Only a conservative terminal-meltout transaction may be reused.');
        % The terminal water/enthalpy transfer was already integrated exactly
        % once. Refresh only the exposed-bare contribution from the restaged
        % SCOPE pass; do not advance S/S_l/U or precipitation a second time.
        fEvent = snowGroundResponse.cover_fraction;
        bareContribution = (1 - fEvent) * bareSurfaceClosure;
        snowFlux.bare_soil_heat_contribution_W_m2 = bareContribution;
        snowFlux.bare_ground_downward_soil_heat_flux_contribution_W_m2 = ...
            bareContribution;
        snowFlux.soil_heat_flux_W_m2 = bareContribution ...
            + snowFlux.snow_soil_heat_contribution_W_m2;
        snowFlux.soil_top_downward_energy_flux_W_m2 = ...
            snowFlux.soil_heat_flux_W_m2;
    else
        [snowNext, snowFlux] = surface.snow.advanceBulkSnow( ...
            built.Y.snow, built.parameters.snow, drivers, dt_outer_s);
    end
    topHeatResponse = surface.energy.buildSoilTopHeatResponse( ...
        scopeStep, snowGroundResponse, snowFlux);
    drivers.snow_energy_redistribution_W_m2 = ...
        topHeatResponse.snow_energy_redistribution_W_m2;

    % --- Step 3: heat flux for soil = snow-adjusted; keep SCOPE G for diagnostics ---
    scopeStep.fluxes.Gtot_surface_closure = surfaceClosure;
    scopeStep.fluxes.Gtot = snowFlux.soil_heat_flux_W_m2;
    scopeStep.fluxes.Gtot_bare = ...
        snowFlux.bare_soil_heat_contribution_W_m2;
    scopeStep.fluxes.Gtot_snow = ...
        snowFlux.snow_soil_heat_contribution_W_m2;
    scopeStep.fluxes.soil_top_downward_energy_flux_W_m2 = ...
        snowFlux.soil_top_downward_energy_flux_W_m2;
    scopeStep.fluxes.bare_ground_downward_soil_heat_flux_contribution_W_m2 = ...
        snowFlux.bare_ground_downward_soil_heat_flux_contribution_W_m2;
    scopeStep.fluxes.snow_downward_soil_energy_flux_contribution_W_m2 = ...
        snowFlux.snow_downward_soil_energy_flux_contribution_W_m2;
    scopeStep.fluxes.soil_top_energy_flux_reference_W_m2 = ...
        topHeatResponse.reference_flux_into_soil_W_m2;
    scopeStep.fluxes.soil_top_energy_flux_temperature_derivative_W_m2_K = ...
        topHeatResponse.dflux_dtemperature_W_m2_K;
    model.built.workspace.scope_step = scopeStep;
    model.built.workspace.surface.frozen_transfer_cm = frozenPond;
    snowWorkspace = struct( ...
        'surface_state', snowSurface, ...
        'boundary_state', snowNext, ...
        'phase', phase, ...
        'flux', snowFlux, ...
        'advance_drivers', drivers, ...
        'ground_heat_response', snowGroundResponse);
    snowWorkspace.pre_reconciled_terminal_event = hasOverride;
    if isfield(built.workspace, 'snow') ...
            && isfield(built.workspace.snow, 'surface_correction_pass')
        snowWorkspace.surface_correction_pass = ...
            built.workspace.snow.surface_correction_pass;
    end
    model.built.workspace.snow = snowWorkspace;
    model.built.workspace.surface.top_heat_response = topHeatResponse;

    % --- Step 4a: canopy water state candidate (not committed until outer accept) ---
    canopyWaterNext = built.Y.canopy_water;
    canopyWaterNext.state_h2ocan = evaporation.state_h2ocan;
    canopyWaterNext.state_h2ocan_sun = evaporation.state_h2ocan_sun;
    canopyWaterNext.state_h2ocan_sha = evaporation.state_h2ocan_sha;
    historyRow = struct( ...
        'h2ocan', evaporation.state_h2ocan, ...
        'h2ocan_sun', evaporation.state_h2ocan_sun, ...
        'h2ocan_sha', evaporation.state_h2ocan_sha, ...
        'fwet', interception.fwet, ...
        'fwet_sun', interception.fwet_sun, ...
        'fwet_sha', interception.fwet_sha, ...
        'lE_CI', evaporation.lE_CI_W, ...
        'residual', evaporation.canopy_water_residual);

    % --- Step 4b: soil-facing water rates (throughfall/melt + ground EVAP) ---
    soilEvaporation = evaporation.EVAP;
    % If SCOPE already split bare-ground LE, do not re-weight by snow cover.
    if ~isfield(scopeStep.fluxes, 'lEstot_bare')
        soilEvaporation = soilEvaporation ...
            * snowFlux.soil_evaporation_fraction;
    end
    waterFlux = struct( ...
        'Throughfall', snowFlux.throughfall_liquid_cm_s ...
            + evaporation.DEW, ...
        'CanopyInterception', evaporation.CanopyInterception, ...
        'EVAP', soilEvaporation ...
            + snowFlux.additional_soil_evaporation_cm_s);
    surfaceSoilWater = coupling.materializeSurfaceSoilWaterRequest( ...
        waterFlux);
    topHeatResponse.soil_evaporation_demand_cm_s = waterFlux.EVAP;
    model.built.workspace.surface.top_heat_response = topHeatResponse;

    % --- Step 4c: candidate bundle for stagePostEbal / commitPending ---
    candidates = struct( ...
        'canopy_water', canopyWaterNext, ...
        'snow', snowNext, ...
        'water_flux', waterFlux, ...
        'surface_soil_water', surfaceSoilWater, ...
        'history_row', historyRow, ...
        'snow_diagnostics', snowFlux);
end

function value = readField(container, name, fallback)
% Read a nonempty struct field or return a local fallback.
% Editor: Jizixing Zhu; Date: 2026-07-17.
    value = fallback;
    if isfield(container, name) && ~isempty(container.(name))
        value = container.(name);
    end
end

function amount = frozenPondTransfer(built, surfaceTemperature)
% Stage liquid depression storage as snow ice when surface is at/below freezing.
% Amount is staged only; outer-step accept commits removal from DSTOR and add to snow.
% Editor: Jizixing Zhu; Date: 2026-07-17.
% Change: 2026-07-18 Jizixing Zhu -- an explicitly disabled snow process
%   cannot take ownership of frozen surface water.
    amount = 0;
    if isfield(built.parameters, 'snow') ...
            && isfield(built.parameters.snow, 'enabled') ...
            && ~built.parameters.snow.enabled
        return;
    end
    if surfaceTemperature <= 0 && isfield(built.Y, 'surface_water') ...
            && isfield(built.Y.surface_water, 'DSTOR0')
        amount = max(double(built.Y.surface_water.DSTOR0), 0);  % [cm]
    end
end

function rate = snowVaporRate(scopeStep, snowState, parameters)
% Convert SCOPE ground latent heat (or explicit sublimation) to water rate [cm/s].
% Prefer the direction-explicit snow vapor field; retain old aliases only for
% compatibility with saved SCOPE structs.
% Editor: Jizixing Zhu; Date: 2026-07-17.
% Change: 2026-07-18 Jizixing Zhu -- prefer signed tile vapor so deposition
%         adds mass to snow instead of being discarded by a positive-only alias.
    d = surface.snow.diagnoseState(snowState, parameters);
    if isfield(scopeStep.fluxes, ...
            'snow_upward_water_vapor_flux_kg_m2_s')
        if d.cover_fraction <= eps
            rate = 0;
            return;
        end
        rate = scopeStep.fluxes.snow_upward_water_vapor_flux_kg_m2_s ...
            / d.cover_fraction * 0.1; % signed kg m-2 s-1 -> cm s-1
        return;
    end
    if isfield(scopeStep.fluxes, 'qv_snow_net_kg_m2_s')
        if d.cover_fraction <= eps
            rate = 0;
            return;
        end
        rate = scopeStep.fluxes.qv_snow_net_kg_m2_s ...
            / d.cover_fraction * 0.1; % signed kg m-2 s-1 -> cm s-1
        return;
    end
    if isfield(scopeStep.fluxes, 'sublimation_snow_kg_m2_s')
        if d.cover_fraction <= eps
            rate = 0;
            return;
        end
        rate = scopeStep.fluxes.sublimation_snow_kg_m2_s ...
            / d.cover_fraction * 0.1; % kg m-2 s-1 -> cm s-1
        return;
    end
    latentW = readField(scopeStep.fluxes, 'lEstot', 0);
    if isfield(scopeStep.fluxes, 'lEstot_snow')
        if d.cover_fraction <= eps
            rate = 0;
            return;
        end
        latentW = scopeStep.fluxes.lEstot_snow / d.cover_fraction;
    end
    phys = equations.physicalConstants();
    latent = phys.lambda_vap + (1 - d.liquid_fraction) ...
        * (phys.lambda_sub - phys.lambda_vap);
    rate = latentW / latent * 0.1; % kg m-2 s-1 -> cm s-1
end
