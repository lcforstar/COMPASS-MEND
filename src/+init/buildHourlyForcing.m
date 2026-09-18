function h = buildHourlyForcing(raw, built, timeAxis, solar)
%BUILDHOURLYFORCING  Expand raw measurements to hourly meteorological drivers.
%
%   Inputs:
%     raw       struct  raw measurement layer from buildForcing (F1)
%     built     struct  control.time, domain.site and forcing.config
%     timeAxis  struct  forcing-clock/UTC axis for the forcing horizon
%     solar     struct  solar geometry and clock-time sunrise/sunset
%
%   Returns h with hourly Ta/Ts/RH/wind/radiation/precipitation, the resolved
%   atmospheric drivers Pg/Ca/Oa (measured where present, otherwise simulated --
%   barometric Pg, background Ca/Oa) and nitrogen deposition fields.
%   Non-temperature fields reproduce legacy MeteoData.to_hourly. Air temperature
%   uses either legacy_sine (fixed cosine peak at 15:00 clock time) or
%   sunrise_anchored (three-segment curve anchored to clock-time sunrise/sunset
%   on the clock day index).
%
%   Public (F2) so a groundwater/Dunnian forcing adapter can rebuild hourly
%   drivers after modifying daily raw.Precip (init.applyGroundwaterForcingAdjustment).
%
%   Editor: Jizixing Zhu; Date: 2026-06-26.
%   Change: 2026-06-27 Jizixing Zhu -- resolve Pg/Ca/Oa here (measured|simulated)
%     so F3 derived no longer duplicates them; Oa added to the hourly drivers.
%   Change: 2026-07-01 Jizixing Zhu -- extracted verbatim from buildForcing nested
%     buildHourly to a public init.buildHourlyForcing (rebuildable by the adapter).
%   Change: 2026-08-11 Jizixing Zhu -- derive daily shortwave DOY/year length
%     from forcing-clock dates and map solar-hour centers into clock intervals.
    tcfg       = built.control.time;
    lTimeScale = tcfg.lTimeScale;
    lat        = built.domain.site.latitude;
    fcfg       = getForcingConfig(built);
    rinWatt    = logical(fcfg.Rin_unit_is_watt);
    nraw       = raw.n_raw;
    Nh         = numel(timeAxis.utcTime);

    h = struct();
    if lTimeScale == 0
        % --- hourly input: direct mapping ---
        h.Ta = raw.Ta; h.Tdew = raw.Tdew; h.Tmax = raw.Tmax; h.Tmin = raw.Tmin;
        h.Ts = raw.Ts; h.Tb = raw.Tb; h.RH = raw.RH; h.RHmax = raw.RHmax; h.RHmin = raw.RHmin;
        h.Uu = raw.Uu; h.Uv = raw.Uv;
        if ~all(isnan(raw.U)) && ~all(raw.U == 0)
            U = raw.U;
        else
            U = sqrt(raw.Uu.^2 + raw.Uv.^2);
        end
        h.Rin = raw.Rin; h.Rli = raw.Rli; h.Rn = raw.Rn; h.n = raw.n;
        h.Precip = raw.Precip / 3600 / 10;            % mm/h -> cm/s
        h.Pg = raw.Pg; h.h_SUR = raw.h_SUR; h.Ca = raw.Ca;
        h.NH4_dep = raw.NH4_dep; h.NO3_dep = raw.NO3_dep;
        h.NH4_bot = raw.NH4_bot; h.NO3_bot = raw.NO3_bot;

    elseif lTimeScale == 1
        % --- daily -> hourly conversion ---
        [Ta_, Tmax_, Tmin_, Tdew_, Ts_, Tb_] = deal(zeros(Nh, 1));
        [RH_, RHmax_, RHmin_]                = deal(zeros(Nh, 1));
        [Uu_, Uv_, U_]                       = deal(zeros(Nh, 1));
        [Rn_, Rin_, Rli_]                    = deal(zeros(Nh, 1));
        Precip_ = zeros(Nh, 1); p_ = zeros(Nh, 1);
        n_data  = zeros(Nh, 1); h_SUR_ = zeros(Nh, 1); Ca_data = zeros(Nh, 1);

        % air-temperature method; for sunrise_anchored gather per-day CLOCK
        % sunrise/sunset (constant within a clock day -> take each day's 1st step)
        method = getDiurnalMethod(built);
        if strcmp(method, 'sunrise_anchored')
            srDay = zeros(nraw, 1);  ssDay = zeros(nraw, 1);
            for dd = 1:nraw
                k0 = (dd - 1) * 24 + 1;
                srDay(dd) = solar.sunrise_clock_h(k0);
                ssDay(dd) = solar.sunset_clock_h(k0);
            end
        end

        % shortwave: Gregorian date + apparent-solar centers mapped to clock bins
        Rin_hour = zeros(nraw, 24);
        for d_pre = 1:nraw
            k0 = (d_pre - 1) * 24 + 1;
            clockDate = timeAxis.clockTime(k0);
            actualDoy = day(clockDate, 'dayofyear');
            daysInYear = day(datetime(year(clockDate), 12, 31), 'dayofyear');
            solarNoonClock = ...
                (solar.sunrise_clock_h(k0) + solar.sunset_clock_h(k0)) / 2;
            clockHourCenters = (0.5:1:23.5)';
            solarHourCenters = mod( ...
                clockHourCenters - solarNoonClock + 12, 24);
            Rin_hour(d_pre, :) = surface.radiation.distribute_radiation_fao56( ...
                raw.Rin(d_pre), actualDoy, daysInYear, lat, rinWatt, ...
                solarHourCenters)';
        end

        % precipitation: deterministic SWAT downscale (fixed T_dur/peak/start)
        for day_idx = 1:nraw
            daily_rain_cm = raw.Precip(day_idx) / 10;
            if daily_rain_cm > 0.001
                hourly_storm_cm = equations.downscale_rain_swat(daily_rain_cm, 6, 0.35, 10);
            else
                hourly_storm_cm = zeros(24, 1);
            end
            idx_start = (day_idx - 1) * 24 + 1;
            idx_end   = min(day_idx * 24, Nh);
            valid_len = idx_end - idx_start + 1;
            Precip_(idx_start:idx_end) = hourly_storm_cm(1:valid_len) / 3600;   % cm/h -> cm/s
        end

        for kt = 1:Nh
            d  = min(timeAxis.clockDayIndex(kt), nraw);    % forcing (local-standard) day
            hr = timeAxis.clockHour(kt);                   % clock hour of day 0~24

            % air temperature: legacy fixed sine (peak 15:00) or sunrise-anchored
            Tx = raw.Tmax(d);  Tn = raw.Tmin(d);
            if strcmp(method, 'sunrise_anchored')
                % clock day index + clock hour + clock sunrise/sunset (same base)
                Ta_(kt) = anchoredAirTemp(d, hr, raw.Tmin, raw.Tmax, srDay, ssDay, nraw);
            else
                Ta_(kt) = (Tx + Tn) / 2 + (Tx - Tn) / 2 * cos(2*pi*(hr - 15) / 24);
            end
            Tmax_(kt) = Tx;  Tmin_(kt) = Tn;
            Tdew_(kt) = raw.Tdew(d);                        % tiled (small diurnal range)
            % surface temperature: sine, amplitude 0.6x air, peak lag to 16:00
            Ts_(kt)   = raw.Ts(d) + 0.6 * (Tx - Tn) / 2 * cos(2*pi*(hr - 16) / 24);
            Tb_(kt)   = raw.Tb(d);                          % deep temperature tiled

            % humidity: RHmax/RHmin sine if present, else ea-conservation from Ta
            if ~isnan(raw.RHmax(d)) && ~isnan(raw.RHmin(d)) && raw.RHmax(d) ~= 0 && raw.RHmin(d) ~= 0
                Rx = raw.RHmax(d);  Rn2 = raw.RHmin(d);
                RH_(kt)    = (Rx + Rn2) / 2 + (Rx - Rn2) / 2 * cos(2*pi*(hr - 6) / 24);
                RHmax_(kt) = Rx;  RHmin_(kt) = Rn2;
            else
                ea_day  = raw.RH(d) .* 0.611 .* exp(17.27 .* raw.Ta(d) ./ (raw.Ta(d) + 237.3));
                esat_h  = 0.611 .* exp(17.27 .* Ta_(kt) ./ (Ta_(kt) + 237.3));
                RH_(kt) = min(1, ea_day / esat_h);
                RHmax_(kt) = nan;  RHmin_(kt) = nan;
            end

            % wind: adaptive diurnal cycle (peak 15:00), Ur by daily-mean magnitude
            Uu_(kt) = raw.Uu(d);  Uv_(kt) = raw.Uv(d);
            if ~isnan(raw.U(d)) && raw.U(d) ~= 0
                Um_d = raw.U(d);
            else
                Um_d = sqrt(raw.Uu(d)^2 + raw.Uv(d)^2);
            end
            if isnan(Um_d) || Um_d < 0; Um_d = 0; end
            if     Um_d < 80,  Ur = 2.5;
            elseif Um_d < 200, Ur = 9;
            else,              Ur = 15;
            end
            Umax_d = 2 * Ur * Um_d / (1 + Ur);
            Umin_d = 2        * Um_d / (1 + Ur);
            U_(kt) = Um_d + (Umax_d - Umin_d) / 2 * cos(2*pi*(hr - 15) / 24);

            % radiation
            h_idx    = min(floor(hr) + 1, 24);
            Rin_(kt) = Rin_hour(d, h_idx);
            Rn_(kt)  = raw.Rn(d);  Rli_(kt) = raw.Rli(d);   % tiled (small diurnal range)

            % other daily values tiled
            p_(kt)     = raw.Pg(d);
            n_data(kt) = raw.n(d);
            h_SUR_(kt) = raw.h_SUR(d);
            Ca_data(kt) = raw.Ca(d);
        end

        h.Ta = Ta_; h.Tmax = Tmax_; h.Tmin = Tmin_; h.Tdew = Tdew_; h.Ts = Ts_; h.Tb = Tb_;
        h.RH = RH_; h.RHmax = RHmax_; h.RHmin = RHmin_; h.Uu = Uu_; h.Uv = Uv_;
        U = U_;
        h.Rin = Rin_; h.Rn = Rn_; h.Rli = Rli_; h.n = n_data; h.Precip = Precip_;
        h.Pg = p_; h.h_SUR = h_SUR_; h.Ca = Ca_data;

        % nitrogen: deposition split /24 per hour; bottom BC concentration tiled
        NH4d = zeros(Nh, 1); NO3d = zeros(Nh, 1); NH4b = zeros(Nh, 1); NO3b = zeros(Nh, 1);
        for d_idx = 1:nraw
            hs = (d_idx - 1) * 24 + 1;  he = d_idx * 24;
            if he <= Nh
                NH4d(hs:he) = raw.NH4_dep(d_idx) / 24;
                NO3d(hs:he) = raw.NO3_dep(d_idx) / 24;
                NH4b(hs:he) = raw.NH4_bot(d_idx);
                NO3b(hs:he) = raw.NO3_bot(d_idx);
            end
        end
        h.NH4_dep = NH4d; h.NO3_dep = NO3d; h.NH4_bot = NH4b; h.NO3_bot = NO3b;
    else
        error('buildForcing:badTimeScale', 'lTimeScale must be 0 or 1, got %d', lTimeScale);
    end

    % --- resolve atmospheric drivers to complete series (measured | simulated) ---
    % Pg: measured (csv pg_hpa) where finite, else barometric estimate from site
    %     elevation + hourly Ta. Ca: measured (csv ca_ppm) where finite, else the
    %     [environment] background. Oa: never measured -> always background O2.
    h.Pg = resolvePressure(h.Pg, h.Ta, built.domain.site.elevation);
    if isempty(h.Ca) || ~any(isfinite(h.Ca) & h.Ca > 0)
        h.Ca = fcfg.Ca * ones(Nh, 1);
    end
    h.Oa = fcfg.Oa * ones(Nh, 1);

    % wind height conversion + gustiness 1 m/s floor
    h.U = max(1.0, convertWindHeight(U, fcfg));

    h.n_hours = Nh;

    % start offset: trim hourly fields so hourly(1) aligns to the simulation start
    % (forcing may begin earlier than start). For changwu start == forcing_start
    % (offset 0). NOTE: timeAxis/solar stay on the forcing timeline; align them
    % too if a site sets start > forcing_start (offset > 0).  [TODO offset>0]
    hour_offset = round(seconds(tcfg.start_clock - tcfg.forcing_start_clock) / 3600);
    if hour_offset > 0
        hf = {'Ta','Tmax','Tmin','Tdew','Ts','Tb','RH','RHmax','RHmin', ...
              'Uu','Uv','U','Rin','Rn','Rli','n','Precip','Pg','h_SUR','Ca','Oa', ...
              'NH4_dep','NO3_dep','NH4_bot','NO3_bot'};
        for i = 1:numel(hf)
            f = hf{i};
            if isfield(h, f) && ~isempty(h.(f)) && numel(h.(f)) > hour_offset
                h.(f) = h.(f)(hour_offset+1 : end);
            end
        end
    end
end

% ===========================================================================
function m = getDiurnalMethod(built)
%GETDIURNALMETHOD  Read the hourly air-temperature disaggregation method.
%
%   Inputs:
%     built  struct  optional control.switches.diurnal_temp_method
%
%   Returns:
%     m      char    'legacy_sine' by default, or a configured method string
%
%   Editor: Jizixing Zhu; Date: 2026-06-26.
%   Change: 2026-07-01 Jizixing Zhu — moved with buildHourly into buildHourlyForcing (F2-only local).
    m = 'legacy_sine';
    if isfield(built, 'control') && isfield(built.control, 'switches') ...
            && isfield(built.control.switches, 'diurnal_temp_method')
        m = lower(strtrim(char(built.control.switches.diurnal_temp_method)));
    end
end

% ===========================================================================
function Ta = anchoredAirTemp(d, hourOfDay, Tmin, Tmax, sunriseHour, sunsetHour, nraw)
%ANCHOREDAIRTEMP  Reconstruct hourly air temperature from daily Tmin/Tmax.
%
%   Inputs:
%     d            [-]  forcing-day index
%     hourOfDay    [h]  clock hour within day
%     Tmin/Tmax    [degC] daily minimum/maximum air-temperature vectors
%     sunriseHour  [h] clock-time sunrise vector
%     sunsetHour   [h] clock-time sunset vector
%     nraw         [-] number of raw forcing days
%
%   Returns:
%     Ta           [degC] hourly air temperature for the requested day/hour
%
%   hourOfDay, sunriseHour and sunsetHour must share one time base. For station
%   daily-row Tmax/Tmin this is forcing-clock time, so the day index and hour
%   agree. Boundary days clamp missing neighbours to the current day.
%
%   Segment definition:
%     pre-dawn:       previous sunset temperature linearly cools to today's Tmin
%     daytime:        sine curve from Tmin at sunrise+1.5 h to Tmax near noon+1.5 h
%     post-sunset:    sunset temperature linearly cools to the next day's Tmin
%
%   Editor: Jizixing Zhu; Date: 2026-06-26.
%   Change: 2026-07-01 Jizixing Zhu — moved with buildHourly into buildHourlyForcing (F2-only local).
    sr  = sunriseHour;  ss = sunsetHour;
    dp  = max(d - 1, 1);         % previous day (clamp at first day)
    dn  = min(d + 1, nraw);      % next day (clamp at last day)
    sr_d = sr(d);  ss_d = ss(d);
    minTime_d = sr_d + 1.5;
    if hourOfDay < minTime_d
        % pre-dawn: previous day's sunset temperature -> today's Tmin
        ss_p      = ss(dp);
        Tset_prev = sunsetTemp(Tmin(dp), Tmax(dp), sr(dp), ss(dp));
        Ta = Tset_prev + (Tmin(d) - Tset_prev) * (24 - ss_p + hourOfDay) / (24 - ss_p + minTime_d);
    elseif hourOfDay <= ss_d
        % daytime sine: Tmin -> Tmax -> sunset temperature
        Ta = Tmin(d) + (Tmax(d) - Tmin(d)) * sin(pi * (hourOfDay - sr_d - 1.5) / (ss_d - sr_d));
    else
        % post-sunset: today's sunset temperature -> next day's Tmin
        minTime_next = sr(dn) + 1.5;
        Tset_d = sunsetTemp(Tmin(d), Tmax(d), sr_d, ss_d);
        Ta = Tset_d + (Tmin(dn) - Tset_d) * (hourOfDay - ss_d) / (24 - ss_d + minTime_next);
    end
end

% ===========================================================================
function Ts = sunsetTemp(Tmin, Tmax, sr, ss)
%SUNSETTEMP  Evaluate the daytime sine temperature curve at sunset.
%
%   Inputs:
%     Tmin/Tmax  [degC] daily minimum/maximum air temperature
%     sr/ss      [h]    sunrise and sunset time in the same clock base
%
%   Returns:
%     Ts         [degC] sunset air temperature, the start of the night segment
%
%   Editor: Jizixing Zhu; Date: 2026-06-26.
%   Change: 2026-07-01 Jizixing Zhu — moved with buildHourly into buildHourlyForcing (F2-only local).
    Ts = Tmin + (Tmax - Tmin) * sin(pi * (ss - sr - 1.5) / (ss - sr));
end

% ===========================================================================
function Pg = resolvePressure(Pg, Ta, Z)
%RESOLVEPRESSURE  Measured surface pressure, or barometric estimate when absent.
%
%   Inputs:
%     Pg  [Pa]    measured pressure series (NaN/empty/<=0 where unmeasured)
%     Ta  [degC]  hourly air temperature (drives the barometric estimate)
%     Z   [m]     site elevation
%
%   Returns Pg with valid measured values preserved and missing/invalid entries
%   replaced by the barometric estimate Pstd*exp(-g*M*Z/(R*TaK)). This pointwise
%   fill is stricter than the legacy all-or-nothing fallback and prevents mixed
%   measured/missing forcing files from carrying NaNs into Esat.
%
%   Editor: Jizixing Zhu; Date: 2026-06-27.
%   Change: 2026-07-01 Jizixing Zhu — moved with buildHourly into buildHourlyForcing (F2-only local).
    phys = equations.physicalConstants();
    Pg_est = phys.Pstd .* exp(-phys.g_barometric * phys.MW_da_kg_mol * Z ./ ...
        (phys.R_universal_SI .* (Ta + phys.C_to_K)));   % [Pa]

    if isempty(Pg)
        Pg = Pg_est;
        return
    end

    Pg = Pg(:);
    bad = ~isfinite(Pg) | Pg <= 0;
    Pg(bad) = Pg_est(bad);
end
