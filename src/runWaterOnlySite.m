function model = runWaterOnlySite(siteDir, options)
%RUNWATERONLYSITE  Run a configured point-laboratory Richards experiment.
%
%   model = runWaterOnlySite(siteDir, options) initializes a laboratory site,
%   prepares its minimal runtime, and advances every configured outer step with
%   solver.runWaterOnlyExperimentStep. Boundary values come only from the
%   site's water_boundary.csv. This loop does not call field canopy, energy,
%   plant, forcing, output, or management services.
%
%   options.stop_after_initialization defaults to false. Set it true to inspect
%   the initialized task without advancing water state.
%
%   Editor: Jizixing Zhu; Date: 2026-07-11.

    if nargin < 2 || isempty(options), options = struct(); end
    stopAfterInitialization = readOption( ...
        options, 'stop_after_initialization', false);

    model = init.initialize(siteDir);
    if ~model.built.control.simulation.is_laboratory
        error('runWaterOnlySite:notLaboratory', ...
            'site.toml must set simulation.point_context = "laboratory".');
    end
    model = init.prepareRuntime(model);
    model.main_loop = struct( ...
        'context', "point_laboratory", ...
        'step_entry', "solver.runWaterOnlyExperimentStep", ...
        'status', "ready", 'step_reports', struct([]));
    if stopAfterInitialization
        model.status = "runtime_prepared";
        return;
    end

    nSteps = model.runtime.run_control.n_steps;
    reports = cell(nSteps, 1);
    completed = 0;
    for kt = 1:nSteps
        remaining = model.built.control.time.TEND - model.runtime.step.TIME;
        dt_s = min(model.built.control.time.Delt_t, remaining);
        [model, reports{kt}] = solver.runWaterOnlyExperimentStep( ...
            model, kt, struct('dt_s', dt_s));
        if ~reports{kt}.accepted
            model.main_loop.status = "stopped_step_failed";
            break;
        end
        completed = kt;
        model.runtime.run_control.last_completed_kt = kt;
        model.runtime.run_control.next_kt = kt + 1;
        model.runtime.step.clock_time = ...
            model.built.control.time.start_time + seconds(model.runtime.step.TIME);
        model.runtime.step.current_day = day(model.runtime.step.clock_time);
        model.runtime.step.current_hour = hour(model.runtime.step.clock_time);
    end
    if completed == nSteps
        model.main_loop.status = "completed";
        model.status = "run_completed";
    else
        model.status = "run_stopped_step_failed";
    end
    if ~isempty(reports)
        usedCount = completed;
        if completed < nSteps && ~isempty(reports{completed + 1})
            usedCount = completed + 1;
        end
        if usedCount > 0
            reports = reports(1:usedCount);
            model.main_loop.step_reports = [reports{:}].';
        end
    end
    model.main_loop.last_completed_kt = completed;
    model.main_loop.next_kt = completed + 1;
end

% -------------------------------------------------------------------------
function value = readOption(options, name, fallback)
% Read one optional runner setting.
% Editor: Jizixing Zhu; Date: 2026-07-11.
    value = fallback;
    if isfield(options, name) && ~isempty(options.(name))
        value = options.(name);
    end
end
