function report = emptyStepReport()
%EMPTYSTEPREPORT  Single-source template for the per-step runOneStep report.
%
%   report = solver.emptyStepReport() returns the canonical per-step report struct
%   with every field defined and neutrally initialised. It is the ONE place the
%   report schema lives: runCOMPASS preallocates the report vector with it, and
%   runOneStep.buildReport fills a copy of it. Both deriving from this template
%   guarantees the field set can never drift, so the struct-array assignment
%   reports(kt) = report never hits "dissimilar structures".
%
%   Returns:
%     report  struct  fields: kt [-], clock_time [datetime], dt [s], accepted [-],
%                     iters [-], residual [-], solve_ctl [struct], step_failed [-],
%                     physics_status [string], forcing_index [-], day_hour [-],
%                     elapsed_seconds [s], status [string], current_day [-],
%                     current_hour [-], checkpoint_file [string], partial_output_written [-].
%
%   Editor: Jizixing Zhu; Date: 2026-07-01.
%   Change: 2026-07-01 Jizixing Zhu — add step-end side-effect fields (current_day/current_hour/checkpoint_file/partial_output_written).
%   Change: 2026-07-02 Jizixing Zhu — rename report.iter -> report.solve_ctl (holds the per-step solve controls).
%   Change: 2026-07-17 Jizixing Zhu -- expose the last coupled trial's
%           termination reason and gate residuals on failed outer steps.
%   Change: 2026-07-17 Jizixing Zhu -- retain the rejected trial's compact
%           surface active-set trace for deterministic cycle diagnosis.
%   Change: 2026-07-19 Jizixing Zhu -- publish per-trial accept/reject counts and
%           compact attempt vectors so rejected_attempts is not a pseudo-zero.
%   Change: 2026-07-20 Jizixing Zhu -- use not_run as the neutral physics status;
%     completed/rejected are assigned by step.buildReport.
%   Change: 2026-08-16 Codex -- add organic-fertilizer accepted-step fields.
%   Change: 2026-08-21 Jizixing Zhu -- add live litter accepted-step fields.
%   Change: 2026-09-12 Jizixing Zhu -- add the manure TAN diagnostic written
%     by runOneStep so canonical report preallocation remains schema-identical.
%   Change: 2026-09-13 Jizixing Zhu -- add source-MEND direct litter DOC/DON
%     fields so annual report preallocation matches runOneStep.
    report = struct();
    report.kt              = 0;
    report.clock_time      = NaT;
    report.dt              = NaN;      % [s] accepted outer-step length
    report.accepted        = false;
    report.iters           = NaN;
    report.residual        = NaN;
    report.solve_ctl       = struct('maxit', NaN, 'maxEBer', NaN, 'Wc', NaN);
    report.step_failed     = false;
    report.physics_status  = "not_run";
    report.forcing_index   = 0;        % [-] forcing row (legacy NoTime)
    report.day_hour        = NaN;      % [-] fraction of day
    report.elapsed_seconds = NaN;      % [s] cumulative elapsed time
    report.status          = "empty";
    report.current_day     = NaN;      % [-] day-of-month of clock_time (legacy current_day)
    report.current_hour    = NaN;      % [-] hour-of-day of clock_time (legacy current_hour)
    report.checkpoint_file = "";       % [string] path if a checkpoint was written this step
    report.partial_output_written = false;   % [-] partial NetCDF snapshot written this step
    report.fert_applied    = false;    % [-] fertilizer event applied this step
    report.fert_nh4_added  = 0;        % [ug/cm2] NH4 mass injected this step
    report.fert_no3_added  = 0;        % [ug/cm2] NO3 mass injected this step
    report.fert_urea_added = 0;        % [ug/cm2] urea mass injected this step
    report.organic_fertilizer_applied = false;
    report.organic_fertilizer_carbon_added_g_cm2 = 0;
    report.organic_fertilizer_nitrogen_added_g_cm2 = 0;
    report.organic_fertilizer_tan_added_g_cm2 = 0;
    report.straw_return_applied = false;
    report.straw_return_carbon_added_g_cm2 = 0;
    report.straw_return_nitrogen_added_g_cm2 = 0;
    report.manure_applied = false;
    report.manure_carbon_added_g_cm2 = 0;
    report.manure_nitrogen_added_g_cm2 = 0;
    report.litter_input_applied = false;
    report.litter_carbon_added_g_cm2 = 0;
    report.litter_nitrogen_added_g_cm2 = 0;
    report.litter_direct_doc_carbon_added_g_cm2 = 0;
    report.litter_direct_don_nitrogen_added_g_cm2 = 0;
    report.attempt_count = 0;          % [-] outer-scheduler trials (accept+reject+surface switch)
    report.accepted_attempt_count = 0; % [-] trials with accepted==true
    report.rejected_attempt_count = 0; % [-] trials with accepted==false (shrink / surface fail)
    report.min_attempt_dt_s = NaN;     % [s] smallest attempted substep this hour
    report.attempt_dt_s = [];          % [s] per-trial substep length
    report.attempt_accepted = [];      % [-] per-trial accept flag
    report.attempt_termination_reason = strings(0, 1);
    report.attempt_iters = [];         % [-] Picard iters per trial
    report.attempt_residual = [];      % residual per trial
    report.attempt_surface_initial_regime = strings(0, 1);
    report.attempt_surface_final_regime = strings(0, 1);
    report.attempt_surface_switch_count = [];
    report.last_attempt_dt_s = NaN;    % [s] final attempted substep
    report.termination_reason = "";    % [-] accepted/nit_exhausted/etc.
    report.nonlinear_converged = false;
    report.mass_converged = false;
    report.theta_increment = NaN;
    report.head_increment_cm = NaN;
    report.water_mass_residual_mm = NaN;
    report.heat_increment_C = NaN;
    report.energy_balance_residual_J_cm2 = NaN;
    report.air_relative_increment = NaN;
    report.air_mass_residual_g_cm2 = NaN;
    report.surface_mass_residual_cm = NaN;
    report.surface_converged = false;
    report.surface_initial_regime = "";
    report.surface_final_regime = "";
    report.surface_trace = "";
end
