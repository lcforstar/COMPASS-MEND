function report = buildReport(kt, step_clock, dt_s, accepted, iters, residual, solve_control, soil_advance)
%BUILDREPORT  Build the fixed report returned by one outer-step attempt.
%
%   Inputs:
%     kt            [-] One-based outer-step index.
%     step_clock    [in] Candidate end-of-step clock from step.computeClock.
%     dt_s          [s] Actual outer-step duration.
%     accepted      [logical] True when the full outer step completed.
%     iters         [-] Solver iteration count reported by the scheduler.
%     residual      [in] Solver residual report.
%     solve_control [in] Numerical controls used for the soil solve.
%
%   Return:
%     report [struct] Fixed-schema step report.
%
%   Editor: Jizixing Zhu; Date: 2026-07-16.
%   Change: 2026-07-17 Jizixing Zhu -- publish the final coupled trial's
%           convergence gates instead of returning only a generic rejection.
%   Change: 2026-07-17 Jizixing Zhu -- expose the compact surface active-set
%           trace without mutating the solver or diagnostic registry.
%   Change: 2026-07-19 Jizixing Zhu -- distinguish bounded surface-soil coupling
%           exhaustion from a genuine nonlinear minimum-substep rejection.
%   Change: 2026-07-19 Jizixing Zhu -- publish rejected/accepted attempt counts
%           and compact per-trial vectors from soil_advance.attempts.
%   Change: 2026-07-20 Jizixing Zhu -- replace obsolete no-physics/pending labels
%     with the actual coupled-physics completion state.
    report = solver.emptyStepReport();
    report.kt = kt;
    report.clock_time = step_clock.clock_time;
    report.dt = dt_s;
    report.accepted = accepted;
    report.iters = iters;
    report.residual = residual;
    report.solve_ctl = solve_control;
    report.step_failed = ~accepted;
    report.forcing_index = step_clock.forcing_index;
    report.day_hour = step_clock.day_hour;
    report.elapsed_seconds = step_clock.TIME;
    if accepted
        report.physics_status = "completed";
        report.status = "accepted";
    else
        report.physics_status = "rejected";
        report.status = "rejected_at_min_dt";
    end
    if nargin >= 8 && isstruct(soil_advance) ...
            && isfield(soil_advance, 'attempts')
        attempts = soil_advance.attempts;
        report.attempt_count = numel(attempts);
        if ~isempty(attempts)
            acceptedFlags = false(numel(attempts), 1);
            dtList = zeros(numel(attempts), 1);
            termList = strings(numel(attempts), 1);
            iterList = zeros(numel(attempts), 1);
            residList = nan(numel(attempts), 1);
            sInit = strings(numel(attempts), 1);
            sFinal = strings(numel(attempts), 1);
            sSwitch = zeros(numel(attempts), 1);
            for i = 1:numel(attempts)
                a = attempts(i);
                acceptedFlags(i) = logical(a.accepted);
                dtList(i) = a.dt_s;
                termList(i) = string(a.termination_reason);
                iterList(i) = a.iters;
                residList(i) = a.residual;
                sInit(i) = string(a.surface_initial_regime);
                sFinal(i) = string(a.surface_final_regime);
                sSwitch(i) = a.surface_switch_count;
            end
            report.attempt_accepted = acceptedFlags;
            report.attempt_dt_s = dtList;
            report.attempt_termination_reason = termList;
            report.attempt_iters = iterList;
            report.attempt_residual = residList;
            report.attempt_surface_initial_regime = sInit;
            report.attempt_surface_final_regime = sFinal;
            report.attempt_surface_switch_count = sSwitch;
            report.accepted_attempt_count = nnz(acceptedFlags);
            report.rejected_attempt_count = nnz(~acceptedFlags);
            report.min_attempt_dt_s = min(dtList);

            last = attempts(end);
            report.last_attempt_dt_s = last.dt_s;
            report.termination_reason = string(last.termination_reason);
            report.nonlinear_converged = logical(last.nonlinear_converged);
            report.mass_converged = logical(last.mass_converged);
            report.theta_increment = last.theta_increment;
            report.head_increment_cm = last.head_increment_cm;
            report.water_mass_residual_mm = last.water_mass_residual_mm;
            report.heat_increment_C = last.heat_increment_C;
            report.energy_balance_residual_J_cm2 = ...
                last.energy_balance_residual_J_cm2;
            report.air_relative_increment = last.air_relative_increment;
            report.air_mass_residual_g_cm2 = last.air_mass_residual_g_cm2;
            report.surface_mass_residual_cm = last.surface_mass_residual_cm;
            report.surface_converged = logical(last.surface_converged);
            report.surface_initial_regime = ...
                string(last.surface_initial_regime);
            report.surface_final_regime = ...
                string(last.surface_final_regime);
            report.surface_trace = string(last.surface_trace);
        end
        if isfield(soil_advance, 'outer_termination_reason') ...
                && ~isempty(soil_advance.outer_termination_reason)
            report.termination_reason = string( ...
                soil_advance.outer_termination_reason);
            report.status = "rejected_surface_soil_coupling";
        end
    end
end
