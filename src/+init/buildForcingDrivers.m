function forcing = buildForcingDrivers(forcing, built)
%BUILDFORCINGDRIVERS Build source-timeline drivers and align their common start.
%
%   Inputs:
%     forcing [struct] complete raw source series and forcing configuration;
%             unrelated fields (e.g. runoffDunnian) are retained.
%     built   [struct] control.time/switches, domain.site and forcing.config.
%
%   Returns forcing with timeAxis, solar and hourly on one identical suffix
%   beginning exactly at control.time.start_clock. raw is not shortened and
%   clockDayIndex keeps its original source-day meaning. n_hours records the
%   selected length. derived/gas are rebuilt from the aligned meteorology.
%
%   Daily input is first disaggregated over the complete source dates; this
%   retains neighbouring-day temperature and daily rain/radiation allocation.
%   Instantaneous drivers retain their interval-start timestamps. No change
%   is made to the solar algorithm, accepted-state clock or output labels.
%   An off-grid/out-of-source start is rejected (invalidStart); every time
%   series must have the common source length (invalidSeriesLength).
%
%   Editor: Wenjie Liu; Date: 2026-10-28.
    raw = forcing.raw;
    if built.control.time.lTimeScale == 0
        Nh = raw.n_raw;
    else
        Nh = raw.n_raw * 24;
    end

    sourceTime = init.buildTimeAxis(built, Nh);
    first = find(sourceTime.clockTime == built.control.time.start_clock, 1);
    if isempty(first)
        error('buildForcingDrivers:invalidStart', ...
            'Simulation start must match a forcing timestamp.');
    end
    sourceSolar = buildSolar(sourceTime, built);
    sourceHourly = init.buildHourlyForcing(raw, built, sourceTime, sourceSolar);

    rows = first:Nh;
    forcing.timeAxis = selectTimeRows(sourceTime, rows, Nh);
    forcing.solar = selectTimeRows(sourceSolar, rows, Nh);
    forcing.hourly = selectTimeRows(sourceHourly, rows, Nh);
    forcing.derived = init.buildDerivedForcing(forcing.hourly, built);
    forcing.gas = init.buildAtmosphericGasDrivers(forcing.hourly, forcing.derived);
end

function out = selectTimeRows(source, rows, nSource)
% Select one shared time suffix, retaining source-day indexes and vector shape.
% Editor: Wenjie Liu; Date: 2026-10-28.
    out = source;
    names = fieldnames(source);
    for k = 1:numel(names)
        name = names{k};
        if strcmp(name, 'n_hours')
            continue
        end
        values = source.(name);
        if ~isvector(values) || numel(values) ~= nSource
            error('buildForcingDrivers:invalidSeriesLength', ...
                'Field %s must contain %d time samples.', name, nSource);
        end
        out.(name) = values(rows);
    end
    if isfield(source, 'n_hours')
        out.n_hours = numel(rows);
    end
end

% ===========================================================================
function s = buildSolar(timeAxis, built)
%BUILDSOLAR  Compute per-step solar geometry from UTC time and site location.
%
%   Inputs:
%     timeAxis  struct  output from init.buildTimeAxis, especially utcTime/clockHour
%     built     struct  needs built.domain.site.latitude/longitude
%
%   Returns s with solar zenith angle, solar altitude, azimuth, Sun-Earth
%   distance factor, hour angle, apparent solar hour, day length, and sunrise/
%   sunset expressed in both apparent solar time and forcing-clock time.
%
%   Method:
%     Reuses surface.radiation.calculatesunpos_ultimate_azi_atan22 (method 0). The
%     time base is timeAxis.utcTime, so no administrative timezone approximation
%     such as round(lon/15) enters the geometry. solarHour is computed as
%     mod(hour_angle / 15 + 12, 24), where hour_angle is zero at solar noon.
%
%   Reference: replaces legacy MeteoData.compute_solar (legacy used round(lon/15)).
%
%   Editor: Jizixing Zhu; Date: 2026-06-26.
    lat = built.domain.site.latitude;     %#ok<NASGU> referenced inside evalc string
    lon = built.domain.site.longitude;    %#ok<NASGU> referenced inside evalc string
    utc = timeAxis.utcTime;
    Nh  = numel(utc);

    s.tts = zeros(Nh, 1); s.sun_altitude = zeros(Nh, 1); s.azi = zeros(Nh, 1);
    s.sun_R = zeros(Nh, 1); s.hour_angle = zeros(Nh, 1); s.day_length_h = zeros(Nh, 1);

    dl_cache = containers.Map('KeyType', 'double', 'ValueType', 'double');
    for kt = 1:Nh
        t  = utc(kt);
        yr = year(t);  mo = month(t);  dy = day(t);                       %#ok<NASGU> used in evalc
        gmtime = hour(t) + minute(t) / 60 + second(t) / 3600;            %#ok<NASGU> used in evalc
        % evalc silences the function's internal fprintf
        [~] = evalc(['[solarz, altitude, ~, azi_north, R, hour_angle, ~, ~, day_length, ~] = ' ...
            'surface.radiation.calculatesunpos_ultimate_azi_atan22(yr, mo, dy, gmtime, lat, lon, 0);']);
        s.tts(kt)          = max(0, min(90, solarz));
        s.sun_altitude(kt) = altitude;
        s.azi(kt)          = azi_north;
        s.sun_R(kt)        = R;
        s.hour_angle(kt)   = hour_angle;
        dk = timeAxis.clockDayIndex(kt);            % cache day length per clock day
        if ~isKey(dl_cache, dk); dl_cache(dk) = day_length; end
        s.day_length_h(kt) = dl_cache(dk);
    end

    % apparent (true) solar hour, and the AST<->clock offset (small, ~const/day)
    s.solarHour = trueSolarHour(s.hour_angle);                  % mod(ha/15 + 12, 24)
    offset = mod(s.solarHour - timeAxis.clockHour + 12, 24) - 12;   % AST - clock [h]
    solarNoonClock = 12 - offset;                              % clock time of solar noon

    s.sunrise_solar_h = 12 - s.day_length_h / 2;               % symmetric about solar noon (AST)
    s.sunset_solar_h  = 12 + s.day_length_h / 2;
    s.sunrise_clock_h = solarNoonClock - s.day_length_h / 2;   % same events in clock time
    s.sunset_clock_h  = solarNoonClock + s.day_length_h / 2;
    s.n_hours = Nh;
end

