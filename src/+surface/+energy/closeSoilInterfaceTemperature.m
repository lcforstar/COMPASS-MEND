function [closure, bestPayload] = closeSoilInterfaceTemperature( ...
        initialTemperature_C, evaluateTrial, options)
%CLOSESOILINTERFACETEMPERATURE  Match soil reaction to frozen surface heat.
%
%   Inputs:
%     initialTemperature_C [degC] Guarded SCOPE Tss initial guess.
%     evaluateTrial [in] function handle returning [evaluation,payload] for
%                   one strict-Dirichlet trial temperature. evaluation fields:
%                   accepted, soil_flux_W_m2, surface_flux_W_m2.
%     options [in] struct with tolerance_W_m2, max_iterations,
%              perturbation_C, min_temperature_C, max_temperature_C.
%
%   Returns:
%     closure struct  Scalar closure audit: temperatures, fluxes, residuals,
%                     iteration count, and convergence flag.
%     bestPayload     Opaque payload from the smallest-|F| accepted trial.
%
%   Method:
%     Evaluate the initial guess and one directed probe, then use at most
%     max_iterations guarded secant updates. No atmospheric energy-balance
%     function is called here.
%
%   Reference:
%     docs/plans/soil_interface_temperature_closure_plan.md, section 5.
% Editor: Jizixing Zhu; Date: 2026-07-19.

    assert(isscalar(initialTemperature_C) && isfinite(initialTemperature_C), ...
        'closeSoilInterfaceTemperature:badInitialTemperature', ...
        'The initial interface temperature must be finite.');
    assert(isa(evaluateTrial, 'function_handle'), ...
        'closeSoilInterfaceTemperature:badEvaluator', ...
        'evaluateTrial must be a function handle.');
    opts = normalizeOptions(options);

    maxEvaluations = opts.max_iterations + 2;
    temperatures = nan(maxEvaluations, 1);
    residuals = nan(maxEvaluations, 1);
    soilFluxes = nan(maxEvaluations, 1);
    surfaceFluxes = nan(maxEvaluations, 1);
    accepted = false(maxEvaluations, 1);
    payloads = cell(maxEvaluations, 1);
    bestIndex = 0;
    converged = false;

    for k = 1:maxEvaluations
        if k == 1
            temperature = clamp(initialTemperature_C, opts);
        elseif k == 2
            direction = opts.rejected_probe_direction;
            if accepted(1) && isfinite(residuals(1))
                direction = -sign(residuals(1));
                if direction == 0, direction = 1; end
            end
            probe_C = opts.perturbation_C;
            if accepted(1) && isfinite(residuals(1))
                probe_C = min(max(probe_C, abs(residuals(1)) ...
                    / opts.probe_response_scale_W_m2_K), ...
                    opts.max_probe_C);
            end
            temperature = clamp(temperatures(1) ...
                + direction * probe_C, opts);
        else
            temperature = secantTemperature( ...
                temperatures, residuals, accepted, k, opts);
        end
        [evaluation, payloads{k}] = evaluateTrial(temperature);
        validateEvaluation(evaluation);
        temperatures(k) = temperature;
        accepted(k) = logical(evaluation.accepted);
        if accepted(k)
            soilFluxes(k) = double(evaluation.soil_flux_W_m2);
            surfaceFluxes(k) = double(evaluation.surface_flux_W_m2);
            residuals(k) = soilFluxes(k) - surfaceFluxes(k);
            if bestIndex == 0 || abs(residuals(k)) < abs(residuals(bestIndex))
                bestIndex = k;
            end
            tolerance = max(opts.tolerance_W_m2, ...
                opts.relative_tolerance ...
                * max(1, abs(surfaceFluxes(k))));
            if abs(residuals(k)) <= tolerance
                converged = true;
                break;
            end
        end
    end

    used = find(~isnan(temperatures), 1, 'last');
    if isempty(used), used = 0; end
    if bestIndex == 0
        bestPayload = [];
        bestTemperature = initialTemperature_C;
        bestSoilFlux = NaN;
        bestSurfaceFlux = NaN;
        initialResidual = NaN;
        finalResidual = NaN;
    else
        bestPayload = payloads{bestIndex};
        bestTemperature = temperatures(bestIndex);
        bestSoilFlux = soilFluxes(bestIndex);
        bestSurfaceFlux = surfaceFluxes(bestIndex);
        firstAccepted = find(accepted, 1, 'first');
        initialResidual = residuals(firstAccepted);
        finalResidual = residuals(bestIndex);
    end
    closure = struct( ...
        'enabled', true, ...
        'converged', converged, ...
        'iterations', max(used - 2, 0), ...
        'evaluations', used, ...
        'initial_temperature_C', initialTemperature_C, ...
        'final_temperature_C', bestTemperature, ...
        'initial_residual_W_m2', initialResidual, ...
        'final_residual_W_m2', finalResidual, ...
        'soil_flux_W_m2', bestSoilFlux, ...
        'surface_flux_W_m2', bestSurfaceFlux, ...
        'temperatures_C', temperatures(1:used), ...
        'residuals_W_m2', residuals(1:used), ...
        'trial_accepted', accepted(1:used));
end

function opts = normalizeOptions(options)
% Validate scalar closure controls and supply bounded defaults.
% Editor: Jizixing Zhu; Date: 2026-07-19.
    if nargin < 1 || ~isstruct(options), options = struct(); end
    opts = struct( ...
        'tolerance_W_m2', readOption(options, 'tolerance_W_m2', 5), ...
        'relative_tolerance', readOption(options, 'relative_tolerance', 0.01), ...
        'max_iterations', round(readOption(options, 'max_iterations', 6)), ...
        'perturbation_C', readOption(options, 'perturbation_C', 1), ...
        'min_temperature_C', readOption(options, 'min_temperature_C', -60), ...
        'max_temperature_C', readOption(options, 'max_temperature_C', 90), ...
        'probe_response_scale_W_m2_K', readOption( ...
            options, 'probe_response_scale_W_m2_K', 40), ...
        'max_probe_C', readOption(options, 'max_probe_C', 10), ...
        'rejected_probe_direction', readOption( ...
            options, 'rejected_probe_direction', -1), ...
        'max_step_C', readOption(options, 'max_step_C', 20));
    assert(opts.tolerance_W_m2 >= 0 && opts.relative_tolerance >= 0 ...
            && opts.max_iterations >= 1 && opts.max_iterations <= 20 ...
            && opts.perturbation_C > 0 && opts.max_step_C > 0 ...
            && opts.probe_response_scale_W_m2_K > 0 ...
            && opts.max_probe_C >= opts.perturbation_C ...
            && abs(opts.rejected_probe_direction) == 1 ...
            && opts.min_temperature_C < opts.max_temperature_C, ...
        'closeSoilInterfaceTemperature:badOptions', ...
        'Closure options are outside their valid ranges.');
end

function value = secantTemperature(temperatures, residuals, accepted, k, opts)
% Return a guarded secant update from the last two accepted evaluations.
% Editor: Jizixing Zhu; Date: 2026-07-19.
    indices = find(accepted(1:k-1));
    if isempty(indices)
        magnitude = (k - 1) * opts.perturbation_C;
        value = temperatures(1) ...
            + opts.rejected_probe_direction * magnitude;
        value = clamp(value, opts);
        return;
    elseif isscalar(indices)
        acceptedTemperature = temperatures(indices(1));
        direction = sign(acceptedTemperature - temperatures(1));
        if direction == 0
            direction = -sign(residuals(indices(1)));
            if direction == 0, direction = 1; end
        end
        value = acceptedTemperature ...
            + direction * opts.perturbation_C;
        value = clamp(value, opts);
        return;
    end
    i1 = indices(end-1);
    i2 = indices(end);
    denominator = residuals(i2) - residuals(i1);
    if ~isfinite(denominator) || abs(denominator) < 1e-12
        step = -sign(residuals(i2)) * opts.perturbation_C;
    else
        step = -residuals(i2) ...
            * (temperatures(i2) - temperatures(i1)) / denominator;
    end
    step = min(max(step, -opts.max_step_C), opts.max_step_C);
    value = clamp(temperatures(i2) + step, opts);
    if any(abs(value - temperatures(indices)) < 1e-10)
        value = clamp(value - sign(residuals(i2)) ...
            * opts.perturbation_C, opts);
    end
end

function validateEvaluation(evaluation)
% Require one fixed-shape trial result; rejected trials may omit fluxes.
% Editor: Jizixing Zhu; Date: 2026-07-19.
    assert(isstruct(evaluation) && isfield(evaluation, 'accepted') ...
            && isscalar(evaluation.accepted), ...
        'closeSoilInterfaceTemperature:badEvaluation', ...
        'Each trial evaluation must provide scalar accepted.');
    if evaluation.accepted
        required = {'soil_flux_W_m2','surface_flux_W_m2'};
        assert(all(isfield(evaluation, required)) ...
                && all(isfinite([evaluation.soil_flux_W_m2, ...
                    evaluation.surface_flux_W_m2])), ...
            'closeSoilInterfaceTemperature:badEvaluation', ...
            'Accepted trials must provide finite soil and surface fluxes.');
    end
end

function value = readOption(options, name, fallback)
% Read one optional finite scalar control.
% Editor: Jizixing Zhu; Date: 2026-07-19.
    value = fallback;
    if isfield(options, name) && isscalar(options.(name)) ...
            && isfinite(options.(name))
        value = double(options.(name));
    end
end

function value = clamp(value, opts)
% Clamp an interface-temperature trial to the guarded physical band [degC].
% Editor: Jizixing Zhu; Date: 2026-07-19.
    value = min(max(value, opts.min_temperature_C), opts.max_temperature_C);
end
