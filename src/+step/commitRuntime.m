function runtime = commitRuntime(runtime, step_clock, kt)
%COMMITRUNTIME  Commit clock and resume pointers for an accepted outer step.
%
%   Inputs:
%     runtime    [inout] Runtime service state.
%     step_clock [in] Candidate clock returned by step.computeClock.
%     kt         [-] Accepted one-based outer-step index.
%
%   Return:
%     runtime    [out] Runtime advanced to the accepted step end.
%
%   A rejected step must not call this function; its clock, calendar, and
%   resume pointer remain at the last accepted step.
%   Editor: Jizixing Zhu; Date: 2026-07-16.
    runtime.step = step_clock;
    runtime.step.current_day = day(step_clock.clock_time);
    runtime.step.current_hour = hour(step_clock.clock_time);
    runtime.calendar.current_date = ...
        dateshift(step_clock.clock_time, 'start', 'day');
    runtime.calendar.current_year = year(step_clock.clock_time);
    runtime.run_control.last_completed_kt = kt;
    runtime.run_control.next_kt = kt + 1;
end
