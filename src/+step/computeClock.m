function step_clock = computeClock(runtime, built, kt, dt_step_s)
%COMPUTECLOCK  Compute accepted-step clock values without mutating runtime.
%
%   Inputs:
%     runtime   [in] Runtime service state at the previous accepted step.
%     built     [in] Built model containing the time configuration.
%     kt        [-] One-based outer-step index.
%     dt_step_s [s] Actual outer-step duration.
%
%   Return:
%     step_clock [struct] End-of-step clock and forcing-index view.
%
%   TIME advances by the full accepted outer-step duration, never by an
%   individual Picard or retry substep. The caller commits this value only
%   after acceptance.
%   Editor: Jizixing Zhu; Date: 2026-07-16.
    time_control = built.control.time;
    dt_reference_s = time_control.Delt_t;
    elapsed_s = step.readAcceptedTime(runtime) + dt_step_s;
    step_clock = struct();
    step_clock.KT = kt;
    step_clock.TIME = elapsed_s;
    step_clock.time_step = dt_step_s;
    step_clock.elapsed_seconds = elapsed_s;
    step_clock.clock_time = time_control.start_time + seconds(elapsed_s);
    step_clock.NoTime = max(1, ceil(elapsed_s / dt_reference_s));
    step_clock.forcing_index = step_clock.NoTime;
    step_clock.day_hour = mod(elapsed_s / 86400, 1);
end

