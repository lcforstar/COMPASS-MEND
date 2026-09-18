function ax = buildTimeAxis(built, nSteps)
%BUILDTIMEAXIS  Per-step clock/UTC time axis over the forcing timeline.
%
%   ax = init.buildTimeAxis(built, nSteps) builds the model time axis: step kt
%   (1..nSteps) sits at forcing_start + (kt-1)*Delt_t. UTC is the internal
%   absolute time used by solar geometry; clockTime mirrors the forcing-file
%   convention and is used to index daily forcing rows.
%
%   Inputs:
%     built    struct  needs control.time (forcing_start_clock/utc, Delt_t)
%     nSteps   [-]      number of forcing steps (e.g. hourly length Nh)
%
%   Returns ax ([nSteps x 1] each):
%     utcTime / clockTime   [datetime]  absolute UTC / forcing-clock time
%     utcHour / clockHour   [h]         hour-of-day 0..24
%     clockDayIndex         [-]         1-based forcing day (= which daily row)
%
%   Solar fields (solarHour, solarDayIndex, sunrise/sunset) are added by
%   init.buildSolar from UTC + lon/lat. clockDayIndex matches the legacy
%   to_hourly day index floor((kt-1)*Delt_t/86400)+1 (forcing row 1 = day 1).
%
%   Editor: Jizixing Zhu; Date: 2026-06-26.

    tcfg = built.control.time;
    dt   = tcfg.Delt_t;
    k    = (0:nSteps-1)';

    ax.utcTime   = tcfg.forcing_start_utc   + seconds(k * dt);
    ax.clockTime = tcfg.forcing_start_clock + seconds(k * dt);

    ax.utcHour   = hour(ax.utcTime)   + minute(ax.utcTime)   / 60 + second(ax.utcTime)   / 3600;
    ax.clockHour = hour(ax.clockTime) + minute(ax.clockTime) / 60 + second(ax.clockTime) / 3600;

    ax.clockDayIndex = floor(k * dt / 86400) + 1;
end
