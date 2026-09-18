function kt = validateOneStep(model, kt)
%VALIDATEONESTEP  Validate the model shell and one-based outer-step index.
%
%   Inputs:
%     model [in]  Model transaction shell containing built state and runtime.
%     kt    [in]  [-] One-based outer-step index.
%
%   Returns:
%     kt    [-] Validated integer outer-step index.
%
%   This is a structural contract check; it performs no physical update.
%   Editor: Jizixing Zhu; Date: 2026-07-16.
    if nargin < 1 || ~isstruct(model) || ~isfield(model, 'built')
        error('runOneStep:badModel', ...
            'runOneStep expects a model struct with model.built.');
    end
    if ~isfield(model, 'runtime') || ~isstruct(model.runtime)
        error('runOneStep:missingRuntime', ...
            'runOneStep requires init.prepareRuntime(model) before stepping.');
    end
    if ~isfield(model.built, 'Y') || ~isfield(model.built, 'workspace') || ...
            ~isfield(model.built.workspace, 'step')
        error('runOneStep:badState', ...
            'model.built must contain Y and workspace.step before stepping.');
    end
    if nargin < 2 || isempty(kt) || ~isscalar(kt) || ~isfinite(kt) ...
            || kt ~= round(kt)
        error('runOneStep:badKt', 'kt must be a finite integer scalar.');
    end
    nSteps = model.built.control.time.tS;
    if isfield(model.runtime, 'run_control') ...
            && isfield(model.runtime.run_control, 'n_steps')
        nSteps = model.runtime.run_control.n_steps;
    end
    if kt < 1 || kt > nSteps
        error('runOneStep:badKt', ...
            'kt must be in [1,%d], got %g.', nSteps, kt);
    end
end

