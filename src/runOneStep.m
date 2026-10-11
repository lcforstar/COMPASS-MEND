function [model, report] = runOneStep(model, kt)
% Change: 2026-09-30 Wenjie Liu -- selectively restore COMPASS canopy
% behavior; retain winter-wheat phenology and CoLM C/N allocation.
%RUNONESTEP  Advance one coupled land-surface step as one transaction.
%
%   [model, report] = runOneStep(model, kt) coordinates one outer model step.
%   Process packages own the physics, solver owns soil tiling and retry, and
%   +step owns validation, pending-state commit, clock, and report assembly.
%
%   Inputs:
%     model [inout] struct  Prepared COMPASS model state and runtime.
%     kt    [in]    [-]     One-based outer-step index.
%
%   Returns:
%     model  [inout] struct  Accepted model state.  On rejection, prognostic Y
%                            and accepted time are unchanged; disposable trial
%                            auxiliaries may remain available for diagnostics.
%     report [out]   struct  Acceptance, timing, solver, and artifact summary.
%
%   Coupling and transaction contract:
%     (1) Stage step-entry boundary/management auxiliaries from accepted state.
%     (2) Build current canopy/SCOPE inputs.  Canopy interception supplies the
%         current-hour ground precipitation; a precipitation-only snow surface
%         predictor then updates SCOPE optics/temperature before RTMo/ebal.
%     (3) Advance canopy water and prognostic snow from accepted state.  If the
%         predicted and completed snow regimes disagree, rerun SCOPE/ebal at
%         most once; the snow state is always recomputed from the same accepted
%         state and therefore advances only once.
%     (4) Pack finalized canopy/surface/plant candidates.  The soil scheduler
%         materializes water, heat, and dry-air requests per adaptive substep;
%         the heat request reads the current-step snow candidate, not lagged
%         accepted snow, while h/U/P_g remain one coupled soil solve.
%     (5) Commit all prognostic candidates, ledgers, time, and artifacts only
%         after the complete outer step accepts.  A rejected step leaves Y and
%         accepted time unchanged; staged auxiliaries remain disposable audit
%         data and are never treated as accepted state.
%
%   Editor: Jizixing Zhu; Date: 2026-06-29.
%   Change: 2026-07-16 Jizixing Zhu -- reduce the file to orchestration and move
%     step infrastructure and process implementations to their owner packages.
%   Change: 2026-07-17 Jizixing Zhu -- document same-step snow/SCOPE/thermal-
%     boundary coupling and its single bounded correction transaction.
%   Change: 2026-07-17 Jizixing Zhu -- skip the harmonic bottom-temperature refresh
%     when the selected bottom heat policy does not consume that driver.
%   Change: 2026-07-18 Jizixing Zhu -- retain one accepted SCOPE-to-soil pass; the
%     discarded A1 reaction correction did not materially improve the state.
%   Change: 2026-07-19 Jizixing Zhu -- reconcile snow enthalpy with the solved
%     affine soil-top heat response and allow one transactional event restage.
%   Change: 2026-07-19 Jizixing Zhu -- replay the authoritative snow transition on
%     reconciliation, carry event states into restaging, and report bounded
%     coupling exhaustion separately from a nonlinear minimum-dt rejection.
%   Change: 2026-07-19 Jizixing Zhu -- close frozen surface energy against
%     strict-Dirichlet soil reaction flux with a bounded scalar T* loop.
%   Change: 2026-07-20 Jizixing Zhu -- use the atmospheric surface tangent
%     and upper-20-cm soil response to prescribe one matched top heat flux;
%     pair bottom free drainage with liquid internal-energy export.
%   Change: 2026-07-20 Jizixing Zhu -- advance soil N only after an accepted soil
%     endpoint and install it through the same outer-step transaction.
%   Change: 2026-08-16 Codex -- pass the absolute simulation clock to the
%     organic-fertilizer matcher while retaining elapsed seconds for stepping.
%   Change: 2026-08-21 Jizixing Zhu -- stage DOM, MEMC, and CoLM litter
%     candidates behind the accepted soil endpoint and commit them atomically.
%   Change: 2026-08-30 Codex -- pass staged manure TAN to the accepted soil-N
%     transaction without inserting TAN into the organic AOM state.
%   Change: 2026-09-13 Jizixing Zhu -- report the accepted source-MEND
%     direct DOC/DON fractions of PlantCN litter input.
%   Change: 2026-10-06 Wenjie Liu -- retain zero material-specific report
%     fields when the native soil-carbon/organic-fertilizer branch is disabled.

% ===== 1. Validate the requested outer interval from accepted model time =====
kt = step.validateOneStep(model, kt);

built = model.built;
dtOuter = built.control.time.Delt_t;       % [s] configured outer step
dtMin = built.control.time.Delt_t_min;     % [s] minimum substep floor
strategy = built.control.numerics.substep;
prevTime = step.readAcceptedTime(model.runtime);
dtOuterActual = step.computeOuterDt( ...
    dtOuter, built.control.time, prevTime);

% ===== 2. Stage step-entry solver, surface-water, and management auxiliaries =====
% These values are derived from accepted state.  Surface-water staging creates
% the atmospheric water-boundary context; it does not yet choose a final
% substep active set or commit fertilizer/prognostic storage.
model.built.workspace.step = solver.warmStartTrial( ...
    model.built.workspace.step, model.built.Y);
model.built.workspace.step.L_f = ...
    soil.thermal.freezeLatentHeat(model.built.control.switches);
model.built.workspace.surface = ...
    soil.boundary.stageSurfaceWaterBoundary(model.built, kt);
% Freeze the external MODFLOW/BMI snapshot once for this outer interval.
% It is an auxiliary boundary driver, never a COMPASS-owned groundwater
% state; every adaptive soil retry receives the same snapshot.
model.built.workspace.boundary.groundwater_driver = model.runtime.groundwater;
[model, fertilizerDiag, fertilizerApplied] = ...
    soil.nitrogen.fertilizer.stageStep(model, kt, dtOuterActual);

% ===== 3. Stage the coupled atmosphere--vegetation--surface candidate graph =====
% The coupling layer owns this ordered, non-committing graph.  It returns
% domain candidates and typed soil requests plus the one context permitted for
% a bounded snow-regime restage.
[model, pending, surfaceVegetationContext] = ...
    coupling.stageSurfaceVegetation( ...
    model, kt, dtOuterActual, fertilizerApplied);

% ===== 4. Advance coupled soil h/U/P_g with adaptive boundary-aware substeps =====
% The scheduler re-materializes one immutable water/heat/dry-air request for
% every trial.  Water may switch supply/storage/runoff active sets at fixed dt;
% atmospheric heat is matched once against the upper-20-cm soil predictor and
% then remains one fixed Neumann flux inside h/U/P_g. Rejected trials do not
% enter accepted flux/storage ledgers.
solveControl = solver.gatherSolveControls(built);
solveControl.water_boundary_flux = pending.request.surface_to_soil.water;
solveControl.water_root_sink = ...
    pending.request.vegetation_to_soil.water.element_endpoint_rate_1_s;
if helpers.getfield_safe_numeric( ...
        model.built.control.switches, 'lThmrlefc', 0) == 1 ...
        && usesSeasonalBottomTemperature(model.built)
    model.built = soil.refreshThermalBoundaryAuxiliary( ...
        model.built, kt, kt);
end

solveFunction = solver.resolveAttemptFunction(model);
outerSnapshot = solver.snapshotSubstepState(model.built);
[model, soilAdvance, storageBeforeMm] = soil.advance( ...
    model, kt, dtOuterActual, dtMin, strategy, solveControl, ...
    pending.request.surface_to_soil.water, prevTime, solveFunction);
[model, pending, soilAdvance, storageBeforeMm] = ...
    reconcileSurfaceSoilHeatTransaction( ...
    model, pending, soilAdvance, storageBeforeMm, outerSnapshot, ...
    kt, dtOuterActual, dtMin, strategy, solveControl, prevTime, ...
    solveFunction, surfaceVegetationContext, fertilizerApplied);

substep = soilAdvance.substep;
elapsed = soilAdvance.elapsed_s;
accepted = elapsed >= dtOuterActual - 1e-9; % [s] completion tolerance

% ===== 5. Commit all accepted candidates or restore the soil trial snapshot =====
% pending canopy/surface/plant values never enter Y on failure.  The snapshot
% restores solver-owned prognostic state and ledgers after a partial soil hour.
if accepted
    soilAdvance.surface_energy_storage_before = struct( ...
        'snow_energy_J_cm2', snowEnergy(outerSnapshot.Y));
    pending.candidate.surface.surface_water = ...
        model.built.workspace.step.trial.surface_water;
    soilCandidate = interfaces.captureSoilCandidate(model.built);
    organicAccepted = [];
    organicFertilizerDiag = struct('triggered',false,'event_index',0, ...
        'carbon_added_g_cm2',0,'nitrogen_added_g_cm2',0, ...
        'straw_carbon_added_g_cm2',0,'straw_nitrogen_added_g_cm2',0, ...
        'manure_carbon_added_g_cm2',0,'manure_nitrogen_added_g_cm2',0, ...
        'tan_nitrogen_added_g_cm2',0,'tan_layer_nitrogen_g_cm2',[]);
    organicApplied = [];
    domAccepted = [];
    microbeAccepted = [];
    litterDiag = struct('triggered',false,'carbon_added_g_cm2',0, ...
        'nitrogen_added_g_cm2',0,'direct_doc_carbon_g_cm2',0, ...
        'direct_don_nitrogen_g_cm2',0,'mend_dom_input_fraction',0);
    % Change: 2026-10-10 Wenjie Liu -- stage direct manure into the same
    % detached MEND/DOM owners used by plant litter and the soil C/N solve.
    if helpers.getfield_safe_numeric( ...
            model.built.control.switches,'lSoilDOM',0)==1
        domAccepted=model.built.Y.soil_dom;
    end
    if helpers.getfield_safe_numeric( ...
            model.built.control.switches,'lExplicitSoilMicrobes',0)==1
        microbeAccepted=model.built.Y.soil_microbes;
    end
    if helpers.getfield_safe_numeric( ...
            model.built.control.switches,'lSoilOrganic',0) == 1
        [organicAccepted,organicFertilizerDiag,organicApplied, ...
            microbeAccepted,domAccepted] = ...
            soil.organic.fertilizer.stageStep( ...
                model.built.Y.soil_organic, ...
                model.built.sources.organic_fertilization_schedule, ...
                model.runtime.management.organic_fertilizer.applied, ...
                model.built.control.time.start_time + seconds(prevTime), ...
                dtOuterActual, ...
                model.built.parameters.soil_organic.grid_zplus_cm, ...
                flip(model.built.parameters.soil.bulk_density(:)), ...
                microbeAccepted,domAccepted,model.built.parameters.soil_microbes);
    end
    if helpers.getfield_safe_numeric( ...
            model.built.control.switches,'lLitterInput',0)==1
        litter=soil.litter.stagePlantInput(organicAccepted,microbeAccepted, ...
            domAccepted,model.built.Y.plant_cn, ...
            pending.candidate.vegetation.plant_cn,model.built);
        organicAccepted=litter.organic;
        microbeAccepted=litter.microbes;
        domAccepted=litter.dom;
        litterDiag=litter;
    end
    if helpers.getfield_safe_numeric( ...
            model.built.control.switches, 'lNefc', 0) == 1
        soilCandidate.nitrogen = soil.nitrogen.advanceAcceptedCandidate( ...
            model.built, pending.request.vegetation_to_soil.nitrogen, ...
            soilAdvance, pending.request.surface_to_soil.water, ...
            fertilizerDiag, kt, dtOuterActual, ...
            pending.candidate.surface.surface_water,organicAccepted, ...
            domAccepted,microbeAccepted,organicFertilizerDiag);
        if isfield(soilCandidate.nitrogen,'root_actual') ...
                && isfield(pending.candidate.vegetation,'plant_cn') ...
                && ~isempty(pending.candidate.vegetation.plant_cn) ...
                && (isfield(pending.candidate.vegetation,'root_nitrogen') ...
                || isfield(pending.candidate.vegetation,'SinkS'))
            if isfield(pending.candidate.vegetation,'root_nitrogen')
                requestedU = pending.candidate.vegetation.root_nitrogen.U_N_total;
            else
                sinkReq = pending.candidate.vegetation.SinkS;
                requestedU = sum(sum(mean(sinkReq,2),3).*model.built.domain.dz(:))*0.01;
            end
            actualU = soilCandidate.nitrogen.root_actual.U_N_total;
            pending.candidate.vegetation = vegetation.reconcileRootNitrogenCandidate( ...
                pending.candidate.vegetation,requestedU,actualU,dtOuterActual);
        end
    end
    % Actual ecosystem exchange is materialized only after final root-N and
    % microbial N-feedback reconciliation, before the atomic owner commit.
    % Change: 2026-09-30 Wenjie Liu -- replace empirical Rh with AOM/MEND CO2.
    if model.built.control.switches.lPlantCN == 1 ...
            && model.built.control.switches.lNefc == 1
        pending.ledger.vegetation_to_atmosphere.carbon_result = ...
            interfaces.buildEcosystemCarbonResult(model.built, ...
            pending.candidate.vegetation.plant_cn,soilCandidate.nitrogen,dtOuterActual);
    end
    transaction = interfaces.createTransaction( ...
        pending, soilCandidate, soilAdvance, storageBeforeMm, dtOuterActual);
    model = step.commitTransaction(model, transaction, kt);
    if ~isempty(organicApplied)
        model.runtime.management.organic_fertilizer.applied = organicApplied;
    end
    model.built = solver.commitOuterStepLedgers( ...
        model.built, soilAdvance, storageBeforeMm, dtOuterActual);
else
    failedInterfaceClosure = struct();
    if isfield(soilAdvance, 'interface_temperature_closure')
        failedInterfaceClosure = ...
            soilAdvance.interface_temperature_closure;
    end
    model.built = solver.restoreSubstepState( ...
        model.built, outerSnapshot);
    if ~isempty(fieldnames(failedInterfaceClosure))
        model.built.workspace.surface. ...
            soil_interface_temperature_closure = failedInterfaceClosure;
    end
end

model.runtime.substep = substep;

% ===== 6. Publish clock, diagnostics, and artifacts with accept-only semantics =====
stepClock = step.computeClock( ...
    model.runtime, built, kt, dtOuterActual);
report = step.buildReport( ...
    kt, stepClock, dtOuterActual, accepted, soilAdvance.iters, ...
    soilAdvance.residual, solveControl, soilAdvance);

if accepted
    report.fert_applied = fertilizerDiag.applied;
    report.fert_nh4_added = fertilizerDiag.nh4_added;
    report.fert_no3_added = fertilizerDiag.no3_added;
    report.fert_urea_added = fertilizerDiag.urea_added;
    report.organic_fertilizer_applied = organicFertilizerDiag.triggered;
    report.organic_fertilizer_carbon_added_g_cm2 = ...
        organicFertilizerDiag.carbon_added_g_cm2;
    report.organic_fertilizer_nitrogen_added_g_cm2 = ...
        organicFertilizerDiag.nitrogen_added_g_cm2;
    report.organic_fertilizer_tan_added_g_cm2 = ...
        organicFertilizerDiag.tan_nitrogen_added_g_cm2;
    report.straw_return_applied = ...
        organicFertilizerDiag.straw_carbon_added_g_cm2>0;
    report.straw_return_carbon_added_g_cm2 = ...
        organicFertilizerDiag.straw_carbon_added_g_cm2;
    report.straw_return_nitrogen_added_g_cm2 = ...
        organicFertilizerDiag.straw_nitrogen_added_g_cm2;
    report.manure_applied = organicFertilizerDiag.manure_carbon_added_g_cm2>0;
    report.manure_carbon_added_g_cm2 = ...
        organicFertilizerDiag.manure_carbon_added_g_cm2;
    report.manure_nitrogen_added_g_cm2 = ...
        organicFertilizerDiag.manure_nitrogen_added_g_cm2;
    report.litter_input_applied = litterDiag.triggered;
    report.litter_carbon_added_g_cm2 = litterDiag.carbon_added_g_cm2;
    report.litter_nitrogen_added_g_cm2 = litterDiag.nitrogen_added_g_cm2;
    report.litter_direct_doc_carbon_added_g_cm2 = ...
        litterDiag.direct_doc_carbon_g_cm2;
    report.litter_direct_don_nitrogen_added_g_cm2 = ...
        litterDiag.direct_don_nitrogen_g_cm2;

    model.runtime = step.commitRuntime( ...
        model.runtime, stepClock, kt);
    diagnostics.writeStep(model, kt, report);
    artifacts = io.writePeriodicStepArtifacts(model, kt);
    report = step.finalizeReport( ...
        report, model.runtime, artifacts);
else
    diagnostics.writeStep(model, kt, report);
end
end

function value = snowEnergy(Y)
% Read accepted bulk-snow internal energy before the outer transaction.
% Editor: Jizixing Zhu; Date: 2026-07-19.
    value = 0;
    if isstruct(Y) && isfield(Y, 'snow') && isstruct(Y.snow) ...
            && isfield(Y.snow, 'energy') && isscalar(Y.snow.energy) ...
            && isfinite(Y.snow.energy)
        value = double(Y.snow.energy);
    end
end

function [model, pending, advance, storageBeforeMm] = ...
        reconcileSurfaceSoilHeatTransaction( ...
        model, pending, advance, storageBeforeMm, outerSnapshot, ...
        kt, dtOuter, dtMin, strategy, solveControl, tStart, solveFunction, ...
        surfaceVegetationContext, fertilizerApplied)
% Reconcile prognostic snow energy with the solved z=0 heat response.
%
% Ordinary snow states receive an exact accept-before-commit enthalpy
% correction. A meltout/phase/cover event is recomputed once from the same
% accepted outer-step state using the first solve's surface temperature, then
% the complete h/U/P_g hour is rerun. A second event request rejects the hour
% instead of creating an unbounded surface-soil outer iteration.
% Editor: Jizixing Zhu; Date: 2026-07-19.
    if advance.elapsed_s < dtOuter - 1e-9 ...
            || ~hasSnowGroundResponse(model.built)
        return;
    end
    [model, pending, diagnostic] = applySnowGroundReconciliation( ...
        model, pending, advance, dtOuter);
    if ~diagnostic.requires_restage
        model.built.workspace.surface.thermal_reconciliation = diagnostic;
        return;
    end

    predictorTemperature = model.built.Y.T(end);
    eventSnowSurface = pending.candidate.surface.snow;
    snowTransactionOverride = struct();
    if isfield(diagnostic, 'terminal_meltout_triggered') ...
            && diagnostic.terminal_meltout_triggered
        snowTransactionOverride = struct( ...
            'state', pending.candidate.surface.snow, ...
            'flux', pending.diagnostics.snow);
    end
    model.built = solver.restoreSubstepState(model.built, outerSnapshot);
    model.built.workspace.surface. ...
        soil_surface_temperature_predictor_C = predictorTemperature;
    % A snow-domain/cover event changes ground optics, radiation, turbulent
    % exchange, plant demand, and the soil boundary together. Rebuild that
    % complete candidate graph once from the same accepted Y state; merely
    % swapping the snow bucket would leave stale snow-covered ebal fluxes on
    % an exposed soil surface.
    [model, pending] = coupling.restageSurfaceVegetation( ...
        model, kt, dtOuter, surfaceVegetationContext, fertilizerApplied, ...
        eventSnowSurface, snowTransactionOverride);
    if isfield(model.built.workspace.surface, ...
            'soil_surface_temperature_predictor_C')
        model.built.workspace.surface = rmfield( ...
            model.built.workspace.surface, ...
            'soil_surface_temperature_predictor_C');
    end

    solveControl.water_boundary_flux = pending.request.surface_to_soil.water;
    solveControl.water_root_sink = solver.resolveCurrentSrt(model.built);
    [model, advance, storageBeforeMm] = soil.advance( ...
        model, kt, dtOuter, dtMin, strategy, solveControl, ...
        pending.request.surface_to_soil.water, tStart, solveFunction);
    if advance.elapsed_s < dtOuter - 1e-9
        return;
    end
    if ~hasSnowGroundResponse(model.built)
        [model, pending] = synchronizeAcceptedSurfaceHeatOutputs( ...
            model, pending, advance, pending.candidate.surface.snow, ...
            pending.diagnostics.snow);
        diagnostic.requires_restage = false;
        diagnostic.restage_applied = true;
        diagnostic.reason = 'terminal_meltout_restage_accepted';
        model.built.workspace.surface.thermal_reconciliation = diagnostic;
        return;
    end
    [model, pending, secondDiagnostic] = ...
        applySnowGroundReconciliation(model, pending, advance, dtOuter);
    secondDiagnostic.restage_applied = true;
    secondDiagnostic.first_restage_reason = diagnostic.reason;
    if secondDiagnostic.requires_restage
        advance.elapsed_s = 0;
        secondDiagnostic.reason = 'bounded_surface_soil_correction_exhausted';
        advance.accepted = false;
        advance.outer_termination_reason = ...
            'bounded_surface_soil_correction_exhausted';
    end
    model.built.workspace.surface.thermal_reconciliation = ...
        secondDiagnostic;
end

function tf = hasSnowGroundResponse(built)
% Test for a staged snow-to-soil response owned by the current transaction.
% Editor: Jizixing Zhu; Date: 2026-07-19.
    tf = isfield(built.workspace, 'snow') ...
        && isfield(built.workspace.snow, 'ground_heat_response') ...
        && isstruct(built.workspace.snow.ground_heat_response) ...
        && isfield(built.workspace.snow.ground_heat_response, ...
            'cover_fraction') ...
        && built.workspace.snow.ground_heat_response.cover_fraction > 0 ...
        && isfield(built.workspace.snow, 'advance_drivers') ...
        && isstruct(built.workspace.snow.advance_drivers);
end

function [model, pending, diagnostic] = ...
        applySnowGroundReconciliation(model, pending, advance, dtOuter)
% Replay snow at the accepted interface T and synchronize staged outputs.
% Editor: Jizixing Zhu; Date: 2026-07-19.
% Change: 2026-07-19 Jizixing Zhu -- use the accepted substep-time-mean soil
%         surface temperature, matching the hourly affine heat response.
% Change: 2026-07-19 Jizixing Zhu -- retain the corrected event state even when a
%         full surface/soil restage is required and publish actual heat parts.
% Change: 2026-07-20 Jizixing Zhu -- before commit, read the solved trial
%         temperature rather than the still-unchanged accepted Y state.
    response = model.built.workspace.snow.ground_heat_response;
    drivers = model.built.workspace.snow.advance_drivers;
    [boundary, soilSurfaceTemperature_C] = ...
        acceptedSurfaceHeatState(model.built, advance);
    [state, flux, diagnostic] = ...
        surface.energy.reconcileSnowGroundHeat( ...
        model.built.Y.snow, drivers, response, ...
        model.built.parameters.snow, soilSurfaceTemperature_C, dtOuter);
    [frozenFlux, frozenParts] = ...
        surface.energy.evaluateFrozenSurfaceSoilHeat( ...
            model.built, advance, soilSurfaceTemperature_C);
    reconciledFlux = frozenParts.bare_W_m2 ...
        + flux.snow_soil_heat_contribution_W_m2;
    diagnostic.pre_reconcile_surface_flux_W_m2 = frozenFlux;
    diagnostic.pre_reconcile_bare_flux_W_m2 = frozenParts.bare_W_m2;
    diagnostic.pre_reconcile_potential_bare_latent_W_m2 = ...
        frozenParts.potential_bare_latent_W_m2;
    diagnostic.pre_reconcile_actual_bare_latent_W_m2 = ...
        frozenParts.actual_bare_latent_W_m2;
    diagnostic.post_reconcile_surface_flux_W_m2 = reconciledFlux;
    diagnostic.surface_flux_change_W_m2 = reconciledFlux - frozenFlux;
    [model, pending] = synchronizeAcceptedSurfaceHeatOutputs( ...
        model, pending, advance, state, flux, boundary, ...
        soilSurfaceTemperature_C);
end

function [model, pending] = synchronizeAcceptedSurfaceHeatOutputs( ...
        model, pending, advance, state, flux, boundary, ...
        soilSurfaceTemperature_C)
% Publish one accepted total heat flux and an exact bare/snow decomposition.
% Editor: Jizixing Zhu; Date: 2026-07-19.
    if nargin < 6
        [boundary, soilSurfaceTemperature_C] = ...
            acceptedSurfaceHeatState(model.built, advance);
    end
    pending.candidate.surface.snow = state;
    model.built.workspace.snow.boundary_state = state;
    topResponse = model.built.workspace.surface.top_heat_response;
    actualTopFlux = boundary.heat.top_flux_into_domain * 1e4;
    if ~isfinite(actualTopFlux)
        actualTopFlux = boundary.heat. ...
            top_response_evaluated_flux_into_domain * 1e4;
    end
    if ~isfinite(actualTopFlux)
        actualTopFlux = surface.energy.evaluateFrozenSurfaceSoilHeat( ...
            model.built, advance, soilSurfaceTemperature_C);
    end
    snowContribution = flux.snow_soil_heat_contribution_W_m2;
    bareContribution = actualTopFlux - snowContribution;
    flux.bare_soil_heat_contribution_W_m2 = bareContribution;
    flux.bare_ground_downward_soil_heat_flux_contribution_W_m2 = ...
        bareContribution;
    flux.soil_heat_flux_W_m2 = actualTopFlux;
    flux.soil_top_downward_energy_flux_W_m2 = actualTopFlux;
    pending.diagnostics.snow = flux;
    model.built.workspace.snow.flux = flux;
    topResponse.evaluated_temperature_C = soilSurfaceTemperature_C;
    topResponse.evaluated_flux_into_soil_W_m2 = actualTopFlux;
    model.built.workspace.surface.top_heat_response = topResponse;
    model.built.workspace.scope_step.fluxes.Gtot = actualTopFlux;
    model.built.workspace.scope_step.fluxes.Gtot_bare = bareContribution;
    model.built.workspace.scope_step.fluxes. ...
        soil_top_downward_energy_flux_W_m2 = actualTopFlux;
    model.built.workspace.scope_step.fluxes.Gtot_snow = ...
        snowContribution;
    model.built.workspace.scope_step.fluxes. ...
        bare_ground_downward_soil_heat_flux_contribution_W_m2 = ...
        bareContribution;
    model.built.workspace.scope_step.fluxes. ...
        snow_downward_soil_energy_flux_contribution_W_m2 = ...
        flux.snow_downward_soil_energy_flux_contribution_W_m2;
end

function [boundary, temperature_C] = acceptedSurfaceHeatState(built, advance)
% Read the time-mean surface temperature and evaluated heat response.
% Editor: Jizixing Zhu; Date: 2026-07-19.
    boundary = soil.boundary.finalizeFluxSummary( ...
        advance.boundary_flux_summary);
    temperature_C = NaN;
    if strcmp(boundary.heat.top_applied_kind, 'fixed_temperature') ...
            && isfinite(boundary.heat.top_requested_temperature_C)
        temperature_C = boundary.heat.top_requested_temperature_C;
    elseif isfinite(boundary.heat.top_response_evaluation_temperature_C)
        temperature_C = ...
            boundary.heat.top_response_evaluation_temperature_C;
    end
    if ~isfinite(temperature_C)
        temperature_C = built.workspace.step.trial.T(end);
    end
end

function tf = usesSeasonalBottomTemperature(built)
% Return true only when the current heat request consumes harmonic Tbh.
% Editor: Jizixing Zhu; Date: 2026-07-17.
    tf = isfield(built.boundary, 'heat') ...
        && isfield(built.boundary.heat, 'bottom') ...
        && strcmp(built.boundary.heat.bottom.kind, ...
            'seasonal_bottom_temperature');
end
