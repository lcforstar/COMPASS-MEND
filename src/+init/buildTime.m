function t = buildTime(rawConfig)
%BUILDTIME  Parse [time] and derive run-control time quantities.
%
%   t = init.buildTime(rawConfig) parses the start/end datetimes, step sizes and
%   forcing timestamp convention from site.toml [time]. It deliberately keeps
%   three related anchors:
%       start/end       simulation window
%       forcing_start   first row of the forcing file (may precede start)
%       UTC             internal absolute time used by solar geometry
%   The clock fields preserve the forcing-file convention: either station
%   local-standard time or UTC. Pure config translation -- no physics.
%
%   Inputs:
%     rawConfig  struct  raw.config (uses .time)
%
%   Returns t  struct:
%     clock_basis 'local_standard'|'utc', utc_offset [h]
%     start_clock/end_clock/forcing_start_clock [datetime, forcing clock]
%     start_utc/end_utc/forcing_start_utc       [datetime, UTC]
%     start_time/end_time [datetime] legacy aliases of start_clock/end_clock
%     Delt_t/Delt_t_max/Delt_t_min [s], lTimeScale (daily=1, hourly=0),
%     TEND [s], tS (total outer steps), Nmsrmn (storage length), JN (day numbers).
%
%   clock->UTC uses init/private/clockToUtc. forcing_start defaults to start.
%   Nmsrmn is ALWAYS recomputed from start/end (+10% buffer), never read from
%   config: a leftover smoke-run Nmsrmn would size arrays too short and overrun.
%   JN is the legacy day-of-year lookup table over the clock-time storage window.
%   Reference: legacy TimeSettings.update_derived_variables (TEND/tS/Nmsrmn/JN).
%
%   Editor: Jizixing Zhu; Date: 2026-06-26.

    rawTime = rawConfig.time;

    % --- time-base convention: how the forcing clock relates to UTC ---
    basis = 'utc';
    if isfield(rawTime, 'clock_basis'); basis = lower(strtrim(char(rawTime.clock_basis))); end
    utc_offset = 0;
    if isfield(rawTime, 'utc_offset'); utc_offset = rawTime.utc_offset; end   % [h] clock = UTC + offset
    t.clock_basis = basis;
    t.utc_offset  = utc_offset;

    % --- clock-time anchors (as written in the forcing/config) ---
    t.start_clock = parseIsoDatetime(rawTime.start);
    t.end_clock   = parseIsoDatetime(rawTime.end);
    if isfield(rawTime, 'forcing_start') && ~isempty(rawTime.forcing_start)
        t.forcing_start_clock = parseIsoDatetime(rawTime.forcing_start);
    else
        t.forcing_start_clock = t.start_clock;     % default: sim starts at forcing start
    end

    % --- internal UTC anchors (canonical absolute time) ---
    t.start_utc         = clockToUtc(t.start_clock,         basis, utc_offset);
    t.end_utc           = clockToUtc(t.end_clock,           basis, utc_offset);
    t.forcing_start_utc = clockToUtc(t.forcing_start_clock, basis, utc_offset);

    % --- legacy clock aliases (existing consumers read start_time/end_time) ---
    t.start_time = t.start_clock;
    t.end_time   = t.end_clock;

    t.Delt_t = rawTime.dt;                                      % [s] outer step
    if isfield(rawTime, 'dt_max'); t.Delt_t_max = rawTime.dt_max; else; t.Delt_t_max = rawTime.dt; end
    if isfield(rawTime, 'dt_min'); t.Delt_t_min = rawTime.dt_min; else; t.Delt_t_min = 10;          end

    scale = 'hourly';
    if isfield(rawTime, 'forcing_timescale')
        scale = lower(strtrim(char(rawTime.forcing_timescale)));
    end
    t.lTimeScale = double(strcmp(scale, 'daily'));             % daily=1, hourly=0

    % Duration-based quantities are time-base independent (use clock anchors).
    t.TEND = seconds(t.end_clock - t.start_clock);             % [s] total duration
    t.tS   = ceil(t.TEND / t.Delt_t);                          % total outer steps

    total_hours = ceil(hours(t.end_clock - t.start_clock));
    if total_hours <= 0
        total_hours = 8760;
    end
    t.Nmsrmn = ceil(total_hours * 1.1);                        % storage length, +10% buffer

    % Legacy day-of-year table over the storage window in forcing-clock time.
    t.JN = init.calculateDayNumbers(t.start_clock, t.start_clock + hours(t.Nmsrmn));
end

% ---------------------------------------------------------------------------
function dt = parseIsoDatetime(s)
% Parse an ISO 8601 string (site.toml quotes datetimes, e.g. 2021-01-01T00:00:00).
% Editor: Jizixing Zhu; Date: 2026-06-27.
    s = char(s);
    try
        dt = datetime(s, 'InputFormat', 'yyyy-MM-dd''T''HH:mm:ss');
    catch
        dt = datetime(s);   % fall back to MATLAB's default ISO parser
    end
end
