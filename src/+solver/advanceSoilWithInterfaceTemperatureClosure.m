function [model, advance, storageBeforeMm] = ...
        advanceSoilWithInterfaceTemperatureClosure( ...
        model, kt, dtOuter_s, dtMin_s, strategy, solveControl, ...
        currentFlux, tStart_s, solveFunction)
%ADVANCESOILWITHINTERFACETEMPERATURECLOSURE  Close frozen surface/soil heat.
%
%   Inputs:
%     model        [inout] Model shell at the pre-soil transaction snapshot.
%     kt           [-] One-based outer step.
%     dtOuter_s    [s] Complete atmospheric interval.
%     dtMin_s      [s] Minimum adaptive soil substep.
%     strategy     struct Adaptive soil scheduler controls.
%     solveControl struct Coupled h/U/P_g convergence controls.
%     currentFlux  struct Current water/root boundary drivers.
%     tStart_s     [s] Accepted physical start time.
%     solveFunction function handle for one coupled soil attempt.
%
%   Returns:
%     model           [inout] Best converged soil trial.
%     advance         struct Soil scheduler result plus scalar closure audit.
%     storageBeforeMm [mm] Soil-water storage before the selected trial.
%
%   Method:
%     The atmospheric production policy intersects the frozen post-ebal
%     surface tangent with an upper-20-cm implicit soil heat predictor, then
%     runs h/U/P_g once with the matched fixed Neumann flux. The optional
%     legacy surface-temperature policy retains its bounded scalar replay.
%     Neither path reruns ebal inside Picard.
%
%   Reference:
%     docs/plans/soil_interface_temperature_closure_plan.md, sections 4-6.
% Editor: Jizixing Zhu; Date: 2026-07-19.
% Change: 2026-07-20 Jizixing Zhu -- for the flux-driven production policy,
%         intersect the frozen surface-energy tangent with the upper-20-cm
%         implicit soil heat response and pass one matched Neumann flux to
%         h/U/P_g; never impose the diagnostic skin temperature on soil.
% Change: 2026-07-20 Jizixing Zhu -- retry a discontinuous freeze/thaw heat match on
%         finer accepted soil tiles, then permit an explicitly bounded,
%         ledger-visible residual; otherwise reject the complete outer step.
% Change: 2026-08-15 Jizixing Zhu -- reserve part of the configured combined
%         energy tolerance for the accepted soil and snow equation residuals,
%         so interface convergence cannot consume the entire ledger budget.

    if dynamicFluxClosureApplies(model.built)
        [model, advance, storageBeforeMm] = advanceDynamicFluxClosure( ...
            model, kt, dtOuter_s, dtMin_s, strategy, solveControl, ...
            currentFlux, tStart_s, solveFunction);
        return;
    end

    controls = closureControls(model.built);
    if ~controls.enabled || ~closureApplies(model.built)
        [model.built, storageBeforeMm] = solver.prepareOuterStepLedgers( ...
            model.built, dtOuter_s, kt);
        [model, advance] = solver.advanceSoilEquationOuterStep( ...
            model, kt, dtOuter_s, dtMin_s, strategy, solveControl, ...
            currentFlux, tStart_s, solveFunction);
        advance.interface_temperature_closure = disabledDiagnostic(model.built);
        return;
    end

    initialTemperature_C = model.built.workspace.scope_step.thermal.Tss;
    responseReference_C = model.built.workspace.surface. ...
        top_heat_response.reference_temperature_C;
    controls.rejected_probe_direction = ...
        sign(responseReference_C - initialTemperature_C);
    if controls.rejected_probe_direction == 0
        controls.rejected_probe_direction = -1;
    end
    baseline = model;
    [closure, payload] = surface.energy.closeSoilInterfaceTemperature( ...
        initialTemperature_C, @evaluateTrial, controls);
    if isempty(payload)
        [model.built, storageBeforeMm] = solver.prepareOuterStepLedgers( ...
            model.built, dtOuter_s, kt);
        [model, advance] = solver.advanceSoilEquationOuterStep( ...
            model, kt, dtOuter_s, dtMin_s, strategy, solveControl, ...
            currentFlux, tStart_s, solveFunction);
    else
        model = payload.model;
        advance = payload.advance;
        storageBeforeMm = payload.storage_before_mm;
    end
    model.built.workspace.surface.soil_interface_temperature_C = ...
        closure.final_temperature_C;
    model.built.workspace.surface.soil_interface_temperature_closure = closure;
    advance.interface_temperature_closure = closure;
    if ~closure.converged
        advance.accepted = false;
        advance.elapsed_s = 0;
        advance.outer_termination_reason = ...
            'soil_interface_temperature_closure_failed';
    end

    function [evaluation, trialPayload] = evaluateTrial(temperature_C)
        trialModel = baseline;
        trialModel.built.workspace.surface.soil_interface_temperature_C = ...
            temperature_C;
        [trialModel.built, trialStorage] = ...
            solver.prepareOuterStepLedgers(trialModel.built, dtOuter_s, kt);
        [trialModel, trialAdvance] = solver.advanceSoilEquationOuterStep( ...
            trialModel, kt, dtOuter_s, dtMin_s, strategy, solveControl, ...
            currentFlux, tStart_s, solveFunction);
        accepted = trialAdvance.elapsed_s >= dtOuter_s - 1e-9;
        evaluation = struct('accepted', accepted, ...
            'soil_flux_W_m2', NaN, 'surface_flux_W_m2', NaN);
        if accepted
            boundary = soil.boundary.finalizeFluxSummary( ...
                trialAdvance.boundary_flux_summary);
            evaluation.soil_flux_W_m2 = ...
                boundary.heat.top_flux_into_domain * 1e4;
            evaluation.surface_flux_W_m2 = ...
                surface.energy.evaluateFrozenSurfaceSoilHeat( ...
                trialModel.built, trialAdvance, temperature_C);
        end
        trialPayload = struct('model', trialModel, ...
            'advance', trialAdvance, ...
            'storage_before_mm', trialStorage);
    end
end

function [model, advance, storageBeforeMm] = advanceDynamicFluxClosure( ...
        model, kt, dtOuter_s, dtMin_s, strategy, solveControl, ...
        currentFlux, tStart_s, solveFunction)
% Predict one atmospheric G, refine phase-crossing time tiles when needed,
% then use the least-residual physically valid trial as a bounded last resort.
% Editor: Jizixing Zhu; Date: 2026-07-20.
% Change: 2026-07-20 Jizixing Zhu -- add phase-aware temporal retry and bounded
%         residual acceptance without weakening the strict closure flag.
    baseline = model;
    provisional = solver.materializeSoilBoundaryTrial( ...
        baseline.built, currentFlux, kt, tStart_s, dtOuter_s);
    soilResponse = soil.thermal.estimateDynamicSurfaceAdmittance( ...
        baseline.built, kt, dtOuter_s, provisional.heat);
    surfaceResponse = baseline.built.workspace.surface.top_heat_response;
    heatExchange = coupling.solveSurfaceSoilHeatExchange( ...
        surfaceResponse, soilResponse);
    analyticMatch = heatExchange.match;
    controls = closureControls(baseline.built);
    tolerance = max(0, double(controls.tolerance_W_m2));
    maxEvaluations = 1 + controls.max_iterations;
    substepCaps_s = closureSubstepCaps( ...
        baseline, dtOuter_s, dtMin_s, strategy, controls);
    maxTotalEvaluations = maxEvaluations * numel(substepCaps_s);
    fluxTrials = nan(maxTotalEvaluations, 1);
    residuals = nan(maxTotalEvaluations, 1);
    soilFluxes = nan(maxTotalEvaluations, 1);
    surfaceFluxes = nan(maxTotalEvaluations, 1);
    temperatures = nan(maxTotalEvaluations, 1);
    accepted = false(maxTotalEvaluations, 1);
    payloads = cell(maxTotalEvaluations, 1);
    storages = nan(maxTotalEvaluations, 1);
    evaluationCaps_s = nan(maxTotalEvaluations, 1);
    best = 0;
    converged = false;
    used = 0;
    groupsUsed = 0;

    for group = 1:numel(substepCaps_s)
        groupsUsed = group;
        groupFirst = used + 1;
        fluxTrials(groupFirst) = heatExchange.flux_into_soil_W_m2;
        for localEvaluation = 1:maxEvaluations
            used = used + 1;
            trial = baseline;
            if group > 1
                trial = forceClosureSubstep( ...
                    trial, substepCaps_s(group), dtOuter_s, strategy);
            end
            [trial.built, trialStorage] = solver.prepareOuterStepLedgers( ...
                trial.built, dtOuter_s, kt);
            trial.built.workspace.surface.dynamic_soil_heat_response = ...
                soilResponse;
            trial.built.workspace.surface.dynamic_soil_heat_match = analyticMatch;
            trial.built.workspace.scope_step.fluxes. ...
                soil_top_downward_energy_flux_W_m2 = fluxTrials(used);
            [trial, trialAdvance] = solver.advanceSoilEquationOuterStep( ...
                trial, kt, dtOuter_s, dtMin_s, strategy, solveControl, ...
                currentFlux, tStart_s, solveFunction);
            payloads{used} = struct('model', trial, ...
                'advance', trialAdvance);
            storages(used) = trialStorage;
            evaluationCaps_s(used) = substepCaps_s(group);
            accepted(used) = ...
                trialAdvance.elapsed_s >= dtOuter_s - 1e-9;
            if ~accepted(used)
                break;
            end

            boundary = soil.boundary.finalizeFluxSummary( ...
                trialAdvance.boundary_flux_summary);
            soilFluxes(used) = ...
                boundary.heat.top_flux_into_domain * 1e4;
            temperatures(used) = ...
                boundary.heat.top_response_evaluation_temperature_C;
            if ~isfinite(temperatures(used))
                temperatures(used) = ...
                    trial.built.workspace.step.trial.T(end);
            end
            surfaceFluxes(used) = ...
                surface.energy.evaluateFrozenSurfaceSoilHeat( ...
                    trial.built, trialAdvance, temperatures(used));
            residuals(used) = soilFluxes(used) - surfaceFluxes(used);
            if isfinite(residuals(used)) ...
                    && (best == 0 ...
                        || abs(residuals(used)) < abs(residuals(best)))
                best = used;
            end
            interfaceTolerance = remainingInterfaceTolerance( ...
                tolerance, trialAdvance, trial.built, dtOuter_s);
            if abs(residuals(used)) <= interfaceTolerance
                converged = true;
                break;
            end
            if localEvaluation == 1
                fluxTrials(used + 1) = surfaceFluxes(used);
            elseif localEvaluation < maxEvaluations
                denominator = residuals(used) - residuals(used - 1);
                if abs(denominator) > 1e-12
                    fluxTrials(used + 1) = fluxTrials(used) ...
                        - residuals(used) ...
                        * (fluxTrials(used) - fluxTrials(used - 1)) ...
                        / denominator;
                else
                    fluxTrials(used + 1) = surfaceFluxes(used);
                end
            end
        end
        if converged
            break;
        end
    end

    if best == 0
        chosen = used;
    else
        chosen = best;
    end
    model = payloads{chosen}.model;
    advance = payloads{chosen}.advance;
    storageBeforeMm = storages(chosen);
    boundedTolerance = boundedResidualTolerance( ...
        tolerance, soilFluxes(chosen), surfaceFluxes(chosen), controls);
    boundedAccepted = ~converged && best > 0 ...
        && logical(controls.bounded_acceptance_enabled) ...
        && abs(residuals(chosen)) <= boundedTolerance ...
        && isPhysicalClosureCandidate(model, advance, dtOuter_s, controls);
    finalMatch = analyticMatch;
    finalMatch.predictor_interface_temperature_C = ...
        analyticMatch.interface_temperature_C;
    finalMatch.predictor_flux_into_soil_W_m2 = ...
        analyticMatch.matched_flux_into_soil_W_m2;
    finalMatch.interface_temperature_C = temperatures(chosen);
    finalMatch.matched_flux_into_soil_W_m2 = soilFluxes(chosen);
    finalMatch.surface_flux_at_match_W_m2 = surfaceFluxes(chosen);
    finalMatch.linear_closure_residual_W_m2 = residuals(chosen);
    model.built.workspace.surface.dynamic_soil_heat_match = finalMatch;

    diagnostic = disabledDiagnostic(model.built);
    diagnostic.enabled = true;
    diagnostic.converged = converged;
    diagnostic.iterations = max(used - 1, 0);
    diagnostic.evaluations = used;
    diagnostic.reason = 'dynamic_upper_20cm_flux_match';
    if converged && evaluationCaps_s(chosen) < substepCaps_s(1)
        diagnostic.acceptance_mode = 'accepted_refined_substep';
        diagnostic.acceptance_mode_code = 2;
    elseif converged
        diagnostic.acceptance_mode = 'accepted_strict';
        diagnostic.acceptance_mode_code = 1;
    elseif boundedAccepted
        diagnostic.acceptance_mode = 'accepted_bounded_residual';
        diagnostic.acceptance_mode_code = 3;
    else
        diagnostic.acceptance_mode = 'rejected_rollback';
        diagnostic.acceptance_mode_code = -1;
    end
    diagnostic.accepted_bounded_residual = boundedAccepted;
    diagnostic.bounded_tolerance_W_m2 = boundedTolerance;
    diagnostic.substep_retry_used = ...
        evaluationCaps_s(chosen) < substepCaps_s(1);
    diagnostic.substep_retry_levels_attempted = max(groupsUsed - 1, 0);
    diagnostic.chosen_substep_cap_s = evaluationCaps_s(chosen);
    diagnostic.initial_temperature_C = ...
        analyticMatch.interface_temperature_C;
    diagnostic.final_temperature_C = temperatures(chosen);
    diagnostic.initial_residual_W_m2 = residuals(1);
    diagnostic.final_residual_W_m2 = residuals(chosen);
    diagnostic.soil_flux_W_m2 = soilFluxes(chosen);
    diagnostic.surface_flux_W_m2 = surfaceFluxes(chosen);
    diagnostic.temperatures_C = temperatures(1:used);
    diagnostic.residuals_W_m2 = residuals(1:used);
    diagnostic.trial_accepted = accepted(1:used);
    diagnostic.flux_trials_W_m2 = fluxTrials(1:used);
    diagnostic.soil_fluxes_W_m2 = soilFluxes(1:used);
    diagnostic.surface_fluxes_W_m2 = surfaceFluxes(1:used);
    diagnostic.trial_substep_caps_s = evaluationCaps_s(1:used);
    advance.interface_temperature_closure = diagnostic;
    model.built.workspace.surface.soil_interface_temperature_closure = ...
        diagnostic;
    if ~converged && ~boundedAccepted
        advance.accepted = false;
        advance.elapsed_s = 0;
        advance.outer_termination_reason = ...
            'atmospheric_soil_heat_flux_closure_failed';
    end
end

% -----------------------------------------------------------------------------
function tolerance = remainingInterfaceTolerance( ...
        configuredTolerance, advance, built, dt_s)
% Reserve the combined-ledger budget already consumed by soil/snow equations.
%
% The surface-soil identity is the sum of the soil equation residual, snow
% equation residual, and interface flux mismatch (apart from roundoff). Using
% the full combined tolerance for the interface alone can therefore accept a
% total residual above the configured gate. The triangle-inequality remainder
% is conservative and does not rely on cancellation between independent errors.
% Editor: Jizixing Zhu; Date: 2026-08-15.
    tolerance = max(0, configuredTolerance);
    if tolerance == 0 || ~isstruct(advance) ...
            || ~isfield(advance, 'energy_balance') ...
            || ~isstruct(advance.energy_balance) ...
            || ~isfield(advance.energy_balance, 'residual_J_cm2') ...
            || ~(isscalar(dt_s) && isfinite(dt_s) && dt_s > 0)
        return;
    end
    soilResidual_W_m2 = ...
        abs(double(advance.energy_balance.residual_J_cm2)) * 1e4 / dt_s;
    snowResidual_W_m2 = 0;
    if isfield(built, 'diagnostics') && isfield(built.diagnostics, 'snow') ...
            && isfield(built.diagnostics.snow, ...
                'energy_balance_residual_J_cm2')
        value = double(built.diagnostics.snow.energy_balance_residual_J_cm2);
        if isscalar(value) && isfinite(value)
            snowResidual_W_m2 = abs(value) * 1e4 / dt_s;
        end
    end
    if isfinite(soilResidual_W_m2)
        tolerance = max(0, tolerance ...
            - soilResidual_W_m2 - snowResidual_W_m2);
    end
end

function caps_s = closureSubstepCaps( ...
        model, dtOuter_s, dtMin_s, strategy, controls)
% Build a bounded geometric retry ladder for phase-crossing outer steps.
% Editor: Jizixing Zhu; Date: 2026-07-20.
    firstCap_s = min(dtOuter_s, model.runtime.substep.Delt_t);
    caps_s = firstCap_s;
    enabled = logical(controls.substep_retry_enabled) ...
        && logical(strategy.enabled);
    factor = double(controls.substep_retry_factor);
    levels = max(0, round(double(controls.substep_retry_max_levels)));
    if ~enabled || ~(isfinite(factor) && factor > 0 && factor < 1)
        return;
    end
    for level = 1:levels
        nextCap_s = max(dtMin_s, caps_s(end) * factor);
        if nextCap_s >= caps_s(end) - 1e-9
            break;
        end
        caps_s(end + 1) = nextCap_s; %#ok<AGROW>
        if nextCap_s <= dtMin_s + 1e-9
            break;
        end
    end
end

function model = forceClosureSubstep(model, cap_s, dtOuter_s, strategy)
% Hold a finer soil tile through the hour and the post-event recovery window.
% Editor: Jizixing Zhu; Date: 2026-07-20.
    sub = model.runtime.substep;
    sub.Delt_t = min(sub.Delt_t, cap_s);
    sub.Delt_old = min(sub.Delt_old, cap_s);
    sub.success_count = 0;
    recoveryTiles = max(0, round(double(strategy.hold_after_reject)));
    sub.hold_count = max(sub.hold_count, ...
        ceil(dtOuter_s / cap_s) + recoveryTiles);
    model.runtime.substep = sub;
end

function tolerance_W_m2 = boundedResidualTolerance( ...
        strictTolerance_W_m2, soilFlux_W_m2, surfaceFlux_W_m2, controls)
% Require both an absolute cap and a flux-relative cap for fallback acceptance.
% Editor: Jizixing Zhu; Date: 2026-07-20.
    absoluteCap = max(strictTolerance_W_m2, ...
        double(controls.bounded_tolerance_W_m2));
    relativeCap = max(0, double(controls.bounded_relative_tolerance)) ...
        * max(abs([soilFlux_W_m2, surfaceFlux_W_m2]));
    tolerance_W_m2 = max(strictTolerance_W_m2, ...
        min(absoluteCap, relativeCap));
end

function tf = isPhysicalClosureCandidate(model, advance, dtOuter_s, controls)
% Last-resort acceptance still requires a complete, finite soil endpoint.
% Editor: Jizixing Zhu; Date: 2026-07-20.
    tf = isstruct(advance) && advance.elapsed_s >= dtOuter_s - 1e-9 ...
        && isfield(model.built.workspace, 'step') ...
        && isfield(model.built.workspace.step, 'trial');
    if ~tf, return; end
    trial = model.built.workspace.step.trial;
    for name = {'h','U','T','p_gas'}
        if isfield(trial, name{1}) ...
                && any(~isfinite(trial.(name{1})(:)))
            tf = false;
            return;
        end
    end
    tf = isfield(trial, 'T') ...
        && all(trial.T(:) >= controls.min_temperature_C) ...
        && all(trial.T(:) <= controls.max_temperature_C);
    if ~tf || ~isfield(model.built.workspace, 'derived'), return; end
    derived = model.built.workspace.derived;
    for name = {'Theta_LL','Theta_II'}
        if isfield(derived, name{1})
            phase = derived.(name{1});
            if any(~isfinite(phase(:))) || any(phase(:) < -1e-10) ...
                    || any(phase(:) > 1 + 1e-10)
                tf = false;
                return;
            end
        end
    end
end

function tf = dynamicFluxClosureApplies(built)
% Select the upper-soil dynamic match only for the atmospheric heat policy.
% Editor: Jizixing Zhu; Date: 2026-07-20.
    tf = helpers.getfield_safe_numeric( ...
        built.control.switches, 'lThmrlefc', 0) == 1 ...
        && isfield(built.boundary, 'heat') ...
        && strcmp(built.boundary.heat.top.kind, 'atmospheric') ...
        && isfield(built.workspace, 'surface') ...
        && isfield(built.workspace.surface, 'top_heat_response') ...
        && isstruct(built.workspace.surface.top_heat_response) ...
        && ~isempty(fieldnames(built.workspace.surface.top_heat_response));
end

function controls = closureControls(built)
% Read normalized scalar closure controls from the numerical configuration.
% Editor: Jizixing Zhu; Date: 2026-07-19.
    controls = struct('enabled', false, 'tolerance_W_m2', 5, ...
        'relative_tolerance', 0, 'max_iterations', 6, ...
        'perturbation_C', 1, 'min_temperature_C', -60, ...
        'max_temperature_C', 90, 'probe_response_scale_W_m2_K', 40, ...
        'max_probe_C', 10, 'max_step_C', 20, ...
        'substep_retry_enabled', true, 'substep_retry_factor', 0.5, ...
        'substep_retry_max_levels', 3, ...
        'bounded_acceptance_enabled', true, ...
        'bounded_tolerance_W_m2', 15, ...
        'bounded_relative_tolerance', 0.20);
    if isfield(built.control, 'numerics') ...
            && isfield(built.control.numerics, 'soil_interface_closure')
        configured = built.control.numerics.soil_interface_closure;
        names = fieldnames(configured);
        for i = 1:numel(names)
            controls.(names{i}) = configured.(names{i});
        end
    end
end

function tf = closureApplies(built)
% Require enabled soil heat and the explicit strict-temperature top policy.
% Editor: Jizixing Zhu; Date: 2026-07-19.
    tf = helpers.getfield_safe_numeric( ...
        built.control.switches, 'lThmrlefc', 0) == 1 ...
        && isfield(built.boundary, 'heat') ...
        && strcmp(built.boundary.heat.top.kind, 'surface_temperature') ...
        && isfield(built.workspace, 'scope_step') ...
        && isfield(built.workspace.scope_step, 'thermal') ...
        && isfinite(built.workspace.scope_step.thermal.Tss) ...
        && isfield(built.workspace, 'surface') ...
        && isfield(built.workspace.surface, 'top_heat_response');
end

function diagnostic = disabledDiagnostic(built)
% Return a fixed-shape audit when scalar closure is inactive or inapplicable.
% Editor: Jizixing Zhu; Date: 2026-07-19.
    temperature_C = NaN;
    if isfield(built.workspace, 'scope_step') ...
            && isfield(built.workspace.scope_step, 'thermal') ...
            && isfield(built.workspace.scope_step.thermal, 'Tss')
        temperature_C = built.workspace.scope_step.thermal.Tss;
    end
    diagnostic = struct('enabled', false, 'converged', true, ...
        'iterations', 0, 'evaluations', 0, ...
        'acceptance_mode', 'disabled', ...
        'acceptance_mode_code', 0, ...
        'accepted_bounded_residual', false, ...
        'bounded_tolerance_W_m2', NaN, ...
        'substep_retry_used', false, ...
        'substep_retry_levels_attempted', 0, ...
        'chosen_substep_cap_s', NaN, ...
        'initial_temperature_C', temperature_C, ...
        'final_temperature_C', temperature_C, ...
        'initial_residual_W_m2', NaN, 'final_residual_W_m2', NaN, ...
        'soil_flux_W_m2', NaN, 'surface_flux_W_m2', NaN, ...
        'temperatures_C', zeros(0,1), 'residuals_W_m2', zeros(0,1), ...
        'trial_accepted', false(0,1), ...
        'flux_trials_W_m2', zeros(0,1), ...
        'soil_fluxes_W_m2', zeros(0,1), ...
        'surface_fluxes_W_m2', zeros(0,1), ...
        'trial_substep_caps_s', zeros(0,1));
end
