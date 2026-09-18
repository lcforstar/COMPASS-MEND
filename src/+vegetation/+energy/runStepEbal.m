function model = runStepEbal(model, kt, dt_outer_s)
%RUNSTEPEBAL  Stage one canopy-soil energy-balance solution.
%
%   Inputs:
%     model      [inout] Model shell with finalized RTMo and interception
%                auxiliaries in workspace.scope_step.
%     kt         [-] One-based outer-step index.
%     dt_outer_s [s] Actual outer-step duration.
%
%   Return:
%     model      [out] Model with energy, radiation, root-water-uptake, and
%                hydraulic candidates staged in workspace.scope_step.
%
%   Method:
%     Once per outer step: assembleEbalInput -> solveCanopySoilEnergy ->
%     optional Anderson-off retry -> applyEbalOutput. Does not commit Y.
%     I/O names: plant/+energy/EBAL_IO_CONTRACT.md §5.
%
%   Editor: Jizixing Zhu; Date: 2026-07-16.
%   Change: 2026-07-16 Jizixing Zhu — phase-A target-name path via assemble /
%           solve / apply.
%   Change: 2026-07-16 Jizixing Zhu — document that dispatcher and kernels now
%           contain no legacy soil bridge.
%   Change: 2026-07-18 Jizixing Zhu -- remove the retired Fourier-G feedback gate.

    scope_static = model.built.p.scope_static;
    if ~readField(scope_static.options, 'calc_ebal', 1)
        return;
    end

    in = vegetation.energy.assembleEbalInput(model, kt, dt_outer_s);
    out = vegetation.energy.solveCanopySoilEnergy(in);

    first_pass = out.iteration;
    retry_used = false;
    tol = in.control.tol_eb_W_m2;
    if ~isConverged(out.iteration, tol) && logical(in.control.use_anderson)
        retry_used = true;
        in_retry = in;
        in_retry.control.use_anderson = false;
        in_retry.iteration.max_iter = in.control.max_iter;
        in_retry.iteration.tol_W_m2 = tol;
        in_retry.iteration.relax_weight = 1;
        out = vegetation.energy.solveCanopySoilEnergy(in_retry);
    end

    out.iteration.ebalRetryUsed = retry_used;
    out.iteration.ebalFirstPassCounter = ...
        readField(first_pass, 'counter', NaN);
    out.iteration.ebalFirstPassMaxResidual = ...
        readField(first_pass, 'maxEBer', NaN);
    out.iteration.ebalFirstPassSunlitResidual = ...
        readField(first_pass, 'maxEBercu', NaN);
    out.iteration.ebalFirstPassShadedResidual = ...
        readField(first_pass, 'maxEBerch', NaN);
    out.iteration.ebalFirstPassSoilResidual = ...
        readField(first_pass, 'maxEBers', NaN);

    model = vegetation.energy.applyEbalOutput(model, out);
end

function converged = isConverged(iteration, tolerance)
% Test the finite component maximum against the configured tolerance.
    converged = isstruct(iteration) && isfield(iteration, 'maxEBer') ...
        && isscalar(iteration.maxEBer) && isfinite(iteration.maxEBer) ...
        && iteration.maxEBer <= tolerance;
end

function value = readField(container, name, fallback)
% Read a nonempty struct field or return a local fallback.
    value = fallback;
    if isfield(container, name) && ~isempty(container.(name))
        value = container.(name);
    end
end
