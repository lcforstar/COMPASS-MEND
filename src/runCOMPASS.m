%% 
% Editor: Wenjie Liu; Date: 2026-10-17.
function model = runCOMPASS(siteDir, options)
%RUNCOMPASS  Canonical COMPASS entry: initialization-first facade.
%
%   model = runCOMPASS(siteDir) runs the new initialization stack for one site
%   package, prepares pre-loop runtime services, and by default stops before the
%   first time step. It is intentionally a facade, not a physics kernel:
%
%   Five-layer run flow:
%       1) model = init.initialize(siteDir)    build state
%       2) model = init.prepareRuntime(model)  pre-loop services/slots (pure allocator)
%       3) model = beginRun(model)             pre-loop one-time actions (resume/mkdir/GW)
%       4) for kt: [model,report] = runOneStep(model, kt)   per-step shell
%       5) model = finalizeRun(model)          post-loop diagnostics flush shell
%
%   This is the new canonical run entry (replaces the legacy script COMPASS.m,
%   which stays in parallel and is retired step by step).
%
%   Inputs:
%     siteDir [char]   Site package path. Default: data/sites/changwu.
%     options [struct] Optional controls:
%       prepare_runtime           [logical] default true. Allocate pre-loop
%                                 runtime services after initialization.
%       stop_after_initialization [logical] default false. When false, run the
%                                 full layer 3-5 loop; set true to stop right
%                                 after initialization/runtime preparation.
%       regional_input            [struct|char] optional. If present, call
%                                 init.initializeRegional(siteDir, regional_input)
%                                 instead of init.initialize(siteDir).
%       regional_options          [struct] optional. Passed to initializeRegional.
%
%   Returns:
%     model struct with model.built.Y as accepted initial state, model.built.
%     workspace.step.trial as the first trial state, and model.main_loop as a
%     readable run boundary record.
%
%   Physical closure:
%     No prognostic state is advanced here. The accepted state is produced by
%     initialization only, including the pre-step soil auxiliary refresh. Runtime
%     services must not advance model.built.Y. The future stepping kernel must
%     update Y only after an accepted solve.
%
%   Verification:
%     See tests/test_runCOMPASS.m. The test checks that this facade reaches
%     runtime_prepared status, preserves the Y/trial boundary, and can run a
%     short no-physics layer-4 loop once runOneStep is wired.
%
%   Editor: Jizixing Zhu; Date: 2026-06-29.
%   Change: 2026-06-29 Jizixing Zhu — renamed COMPASS_main -> runCOMPASS (new canonical entry; legacy COMPASS.m stays parallel).
%   Change: 2026-06-29 Jizixing Zhu — wire layer 3 beginRun into the ~stopAfterInit branch (skeleton); layers 4-5 pending.
%   Change: 2026-07-01 Jizixing Zhu — default stop_after_initialization=false (runCOMPASS() now runs the full loop; pass true to stop after init).
%   Change: 2026-07-20 Jizixing Zhu -- save a restart-safe accepted checkpoint when a step rejects.
%   Change: 2026-08-20 Codex -- expose per-step wall time and adaptive-attempt
%   counts in driver-generation runs so isolated slow forcing hours can be
%   diagnosed without changing solver controls.

    % Change: 2026-10-17 Wenjie Liu -- resolved CSV runs audit accepted
    % interval/cumulative C/N, DIN and O2 through the standard driver.
    if nargin < 1 || isempty(siteDir)
        rootDir = fileparts(fileparts(mfilename('fullpath')));
        siteDir = fullfile(rootDir, 'data', 'sites', 'changwu');
    end
    if nargin < 2 || isempty(options)
        options = struct();
    end

    prepareRuntime = readOption(options, 'prepare_runtime', true);
    stopAfterInit  = readOption(options, 'stop_after_initialization', false);

    if isfield(options, 'regional_input') && ~isempty(options.regional_input)
        if isfield(options, 'regional_options')
            model = init.initializeRegional(siteDir, options.regional_input, options.regional_options);
        else
            model = init.initializeRegional(siteDir, options.regional_input);
        end
    else
        model = init.initialize(siteDir);
    end

    if prepareRuntime
        model = init.prepareRuntime(model);
        if isfield(options, 'attempt_solve') && ~isempty(options.attempt_solve)
            model.runtime.debug.attempt_solve = options.attempt_solve;   % test/debug solve injection
        end
    end

    model.main_loop = buildMainLoopRecord(model, prepareRuntime, stopAfterInit);

    if ~stopAfterInit
        if ~isfield(model, 'runtime')
            error('runCOMPASS:missingRuntime', ...
                'stop_after_initialization=false requires prepare_runtime=true.');
        end
        model = beginRun(model);
        nSteps = model.runtime.run_control.n_steps;
        startKt = model.runtime.run_control.next_kt;   % resume continues here; else 1
        fprintf('[COMPASS] KT loop started: KT=%d/%d\n', startKt, nSteps);
        reports = repmat(solver.emptyStepReport(), nSteps, 1);
        driverTiming = strcmpi(strtrim(io.getEnv('DIAGNOSTIC_PROFILE')), 'driver');
        stopped = false;
        auditOn=isfield(model.built.parameters,'soil_organic') ...
            &&isfield(model.built.parameters.soil_organic,'initial_dom');
        if auditOn
            auditInitialY=model.built.Y;auditRows=zeros(nSteps,10);auditDIN=[0 0];
            if startKt>1
                % Resume rebuilds disposable workspace; derive the prior
                % residual from the accepted N ledger, as solveNitrogen does.
                f=model.built.Y.soil_nitrogen.flux;
                previousNH4=(f.cvnh4(end)-f.cv_nh4_init_mg_cm2) - ...
                    (f.cum_nh4_in_mg_cm2-f.cum_nh4_out_mg_cm2);
                previousNO3=(f.cvno3(end)-f.cv_no3_init_mg_cm2) - ...
                    (f.cum_no3_in_mg_cm2-f.cum_no3_out_mg_cm2);
                auditDIN=[previousNH4 previousNO3];
            end
            cumulative=[];
        end
        for kt = startKt:nSteps
            stepWallClock = tic;
            if auditOn,auditPreviousY=model.built.Y;end
            [model, reports(kt)] = runOneStep(model, kt);
            stepWallSeconds = toc(stepWallClock);
            if driverTiming
                fprintf(['[COMPASS] KT=%d/%d status=%s accepted=%d elapsed=%.1fs ' ...
                    'wall=%.3fs attempts=%d rejected=%d min_dt=%.3fs\n'], ...
                    kt, nSteps, char(reports(kt).status), ~reports(kt).step_failed, ...
                    reports(kt).elapsed_seconds, stepWallSeconds, ...
                    reports(kt).attempt_count, reports(kt).rejected_attempt_count, ...
                    reports(kt).min_attempt_dt_s);
            else
                fprintf('[COMPASS] KT=%d/%d status=%s accepted=%d elapsed=%.1fs\n', ...
                    kt, nSteps, char(reports(kt).status), ~reports(kt).step_failed, ...
                    reports(kt).elapsed_seconds);
            end
            if reports(kt).step_failed
                % A rejected forcing interval has no accepted physical endpoint.
                % Do not copy legacy COMPASS.m's "continue at kt+1" policy: that
                % would silently skip forcing, fertilizer, snow, and C/N events.
                % Preserve the useful legacy behaviour instead: persist exactly the
                % last accepted state, together with the failure reason, so a repair
                % can resume at this same kt without reconstructing state by hand.
                model = recordFailedStep(model, reports(kt));
                stopped = true;
                break;
            end
            if auditOn
                [auditRows(kt,:),cumulative]=diagnostics.acceptedOrganicAudit( ...
                    model.built,auditPreviousY,auditInitialY,auditDIN,kt,startKt==1,kt==startKt);
                auditDIN=auditRows(kt,6:7);
                assert(abs(sum(auditRows(startKt:kt,8)))<1e-14, ...
                    'acceptedAudit:oxygenCumulative','Cumulative O2 residual failed.');
            end
        end
        if auditOn&&~isempty(cumulative)
            last=kt;if stopped,last=kt-1;end
            rows=auditRows(startKt:last,:);
            summary=struct('status','PASS','first_kt',startKt,'last_kt',last, ...
                'initialization_source','resolved_organic_initial_pools_csv', ...
                'full_budget',cumulative,'max_hourly_CN_residual_g_cm2',max(abs(rows(:,2:3)),[],1), ...
                'max_interval_DIN_residual_ug_cm2',max(abs(rows(:,4:5)),[],1), ...
                'max_cumulative_DIN_residual_ug_cm2',max(abs(rows(:,6:7)),[],1), ...
                'max_O2_residual_mol_cm2',max(abs(rows(:,8))), ...
                'cumulative_O2_residual_mol_cm2',sum(rows(:,8)), ...
                'final_soc_0_20_mgkg',rows(end,9),'final_doc_free_0_20_mgkg',rows(end,10));
            if stopped,summary.status='STOPPED';end
            out=model.runtime.output.output_dir;
            names={'KT','C_g_cm2','N_g_cm2','NH4_interval_ug_cm2','NO3_interval_ug_cm2', ...
                'NH4_cumulative_ug_cm2','NO3_cumulative_ug_cm2','O2_mol_cm2','SOC_0_20_mgkg','DOC_free_0_20_mgkg'};
            writetable(array2table(rows,'VariableNames',names),fullfile(out,'hourly_conservation_carbon.csv'));
            fid=fopen(fullfile(out,'conservation_summary.json'),'w');assert(fid>=0);
            fprintf(fid,'%s\n',jsonencode(summary,PrettyPrint=true));fclose(fid);
            model.main_loop.conservation=summary;
        end
        model.main_loop.step_reports = reports;
        if stopped
            model.main_loop.status = "stopped_step_failed";
            model.main_loop.last_completed_kt = kt - 1;
            model.main_loop.next_kt = kt;
        else
            model.main_loop.status = "completed";
            model.main_loop.last_completed_kt = nSteps;
            model.main_loop.next_kt = nSteps + 1;
        end
        model = finalizeRun(model);
        if stopped
            model.status = "run_stopped_step_failed";
        else
            model.status = "run_completed";
        end
    end
end

% -------------------------------------------------------------------------
function model = recordFailedStep(model, report)
%RECORDFailedSTEP  Preserve the last accepted restart point after rejection.
%
%   A failed adaptive outer step restores Y to its prior accepted endpoint in
%   runOneStep.  This function records the failed forcing index only as runtime
%   metadata, then writes a normal checkpoint at `last_completed_kt`.  Thus the
%   resume file always represents a physical time level, never a partial trial.
%
%   Editor: Jizixing Zhu; Date: 2026-07-20.
    runtime = model.runtime;
    lastAccepted = runtime.run_control.last_completed_kt;
    assert(runtime.run_control.next_kt == report.kt, ...
        'runCOMPASS:failureClockMismatch', ...
        'A rejected step must leave next_kt at the failed forcing index.');
    runtime.run_control.failure = struct( ...
        'failed_kt', report.kt, ...
        'termination_reason', string(report.termination_reason), ...
        'status', string(report.status), ...
        'last_completed_kt', lastAccepted);
    model.runtime = runtime;
    checkpointFile = "";
    checkpointError = "";
    try
        checkpointFile = string(io.saveCompassCheckpoint(model, lastAccepted));
    catch checkpointException
        checkpointError = string(checkpointException.message);
        warning('runCOMPASS:failedStepCheckpoint', ...
            'Rejected KT=%d; unable to write restart checkpoint: %s', ...
            report.kt, checkpointException.message);
    end
    model.main_loop.failure = struct( ...
        'failed_kt', report.kt, ...
        'termination_reason', string(report.termination_reason), ...
        'restart_checkpoint', checkpointFile, ...
        'checkpoint_error', checkpointError);
end

% -------------------------------------------------------------------------
function value = readOption(options, name, defaultValue)
% Read an option field with a default.
    if isfield(options, name) && ~isempty(options.(name))
        value = options.(name);
    else
        value = defaultValue;
    end
end

% -------------------------------------------------------------------------
function loop = buildMainLoopRecord(model, prepareRuntime, stopAfterInit)
% Build a small record that documents the initialization/run boundary.
    built = model.built;
    loop = struct();
    loop.status = "ready_before_time_loop";
    loop.initialization_entry = "init.initialize";
    if isfield(built.domain, 'kind') && ~strcmp(string(built.domain.kind), "column")
        loop.initialization_entry = "init.initializeRegional";
    end
    loop.runtime_prepared = prepareRuntime && isfield(model, 'runtime');
    loop.stop_after_initialization = logical(stopAfterInit);
    loop.next_step_entry = "runOneStep(model, kt)";
    loop.next_kt = 1;
    loop.n_steps = built.control.time.tS;
    if isfield(model, 'runtime') && isfield(model.runtime, 'run_control')
        loop.n_steps = model.runtime.run_control.n_steps;
    end
    loop.accepted_state = "model.built.Y";
    loop.trial_state = "model.built.workspace.step.trial";
    loop.previous_state = "model.built.workspace.step.previous";
    loop.note = "Initialization and pre-step soil auxiliary refresh are complete; no time step has been advanced.";
end

% -------------------------------------------------------------------------
% (report template moved to src/+solver/emptyStepReport.m — single source shared
%  by runCOMPASS preallocation and runOneStep.buildReport.)
