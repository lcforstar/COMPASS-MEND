function printStepEbalDiagnostics(model, kt)
%PRINTSTEPEBALDIAGNOSTICS  Print gated canopy-energy convergence details.
%
%   Inputs:
%     model [in] Model shell containing staged scope-step energy auxiliaries.
%     kt    [-] One-based outer-step index.
%
%   This read-only diagnostic produces no registry field and mutates no model
%   state. The underlying printer self-gates on non-convergence or invalid
%   temperature ranges.
%   Editor: Jizixing Zhu; Date: 2026-07-16.
    scope_static = model.built.p.scope_static;
    if ~readField(scope_static.options, 'calc_ebal', 1)
        return;
    end
    scope_step = model.built.workspace.scope_step;
    if ~isfield(scope_step, 'iter') || ~isfield(scope_step, 'thermal') ...
            || ~isfield(scope_step, 'fluxes')
        return;
    end
    solver.printEbal45Diagnostics( ...
        kt, scope_step.iter, scope_step.canopy, scope_step.meteo, ...
        scope_step.rad, scope_step.gap, scope_step.thermal, ...
        scope_step.fluxes, scope_step.hydraulic, scope_step.profiles);
end

function value = readField(container, name, fallback)
% Read a nonempty struct field or return a local fallback.
    value = fallback;
    if isfield(container, name) && ~isempty(container.(name))
        value = container.(name);
    end
end
