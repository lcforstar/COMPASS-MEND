function forcing = buildForcing(rawForcing, built)
%BUILDFORCING  Build the meteorological driver container from raw.forcing.
%
%   forcing = init.buildForcing(rawForcing, built) turns the raw named forcing
%   table (raw.forcing from io.loadSite, one row per forcing step at native
%   cadence) into the model driver container. The caller stores the result at
%   built.forcing.
%
%   Time handling is centralized (init.buildTimeAxis / buildTime): the model's
%   internal time is UTC; solar geometry uses UTC + lon/lat ONLY (no
%   administrative timezone round(lon/15)). Daily->hourly temperature for station
%   data uses the CLOCK base (clock day index + clock hour + clock sunrise/sunset)
%   so the daily Tmax/Tmin attach to the correct local-standard day.
%
%   Stages (in build order):
%     F1 raw      named columns + raw units -> legacy measurement layer
%     --          buildTimeAxis: per-step clock/UTC axis
%     F4 solar    UTC + lon/lat -> tts/azi/hour_angle/solarHour + sunrise/sunset
%                 (both clock and solar time) + day_length
%     F2 hourly   init.buildHourlyForcing: daily->hourly expansion (legacy_sine or
%                 sunrise_anchored) + resolved atmospheric drivers Pg/Ca/Oa
%     F3 derived  init.buildDerivedForcing: es/ea/delta/qs/qa/RHOV_A/Theta
%     gas         Pg/Ca/Oa/ea -> dry-air pressure and CO2/O2 partial pressures
%
%   init.buildForcingDrivers builds F2/F4 on the full source timeline, selects
%   the simulation-start suffix for all three driver containers, then builds
%   F3/gas. The groundwater adapter uses that same entry after changing raw rain.
%
%   Inputs:
%     rawForcing  table   raw.forcing (named columns, RAW Excel units)
%     built       struct  needs control.time, domain.site (lat/lon/elev),
%                         forcing.config (wind/radiation/CO2/O2 scalars).
%                         Optional: if omitted, only F1 (raw) is built.
%
%   Returns forcing with:
%     raw        struct  raw measurement layer, original cadence after unit conversion
%     timeAxis   struct  forcing-clock and UTC timestamps for each model hour
%     solar      struct  solar geometry from UTC + site latitude/longitude
%     hourly     struct  hourly meteorological drivers + resolved Pg/Ca/Oa
%     derived    struct  thermodynamic derivatives (es/ea/qs/qa/RHOV_A/Theta)
%     gas        struct  atmospheric gas partial pressures and O2 molar density
%
%   raw and to_hourly non-temperature fields still match legacy MeteoData;
%   solar/temperature use the corrected UTC time base (verified physically).
%
%   Editor: Jizixing Zhu; Date: 2026-06-26.
%   Change: 2026-10-28 Wenjie Liu -- align hourly, solar and time-axis rows
%     together through buildForcingDrivers; retain the complete raw source.
%   Change: 2026-07-01 Jizixing Zhu — split F2 buildHourly -> init.buildHourlyForcing
%     and F3 buildDerived -> init.buildDerivedForcing (public, adapter-rebuildable);
%     getForcingConfig moved to +init/private. buildForcing keeps F1 raw + F4 solar.

    forcing = struct();
    forcing.raw = buildRawLayer(rawForcing);              % F1

    if nargin >= 2 && ~isempty(built)
        forcing.config = getForcingConfig(built);
        forcing = init.buildForcingDrivers(forcing, built);
    end
end

% ===========================================================================
function raw = buildRawLayer(t)
%BUILDRAWLAYER  Convert named forcing columns to the raw measurement layer.
%
%   Inputs:
%     t  table  raw.forcing from io.loadSite, original file cadence and units
%
%   Returns raw with daily/hourly measurement fields. Conversions match
%   legacy MeteoData.load: RH percent -> fraction, Pg hPa -> Pa. Absent nitrogen
%   columns default to zero (legacy behavior), while absent meteorological
%   columns default to NaN so later checks can detect missing drivers.
%
%   Editor: Jizixing Zhu; Date: 2026-06-26.
    nrow = height(t);
    raw = struct();
    % temperature [degC]
    raw.Tmax = getcol(t, 'tmax_c',  NaN, nrow);
    raw.Tmin = getcol(t, 'tmin_c',  NaN, nrow);
    raw.Ta   = getcol(t, 'tmean_c', NaN, nrow);   % mean air temperature
    raw.Tdew = getcol(t, 'tdew_c',  NaN, nrow);
    % relative humidity: percent -> fraction
    raw.RH    = getcol(t, 'rh_mean_pct', NaN, nrow) * 0.01;
    raw.RHmax = getcol(t, 'rh_max_pct',  NaN, nrow) * 0.01;
    raw.RHmin = getcol(t, 'rh_min_pct',  NaN, nrow) * 0.01;
    % wind [m/s]
    raw.Uu = getcol(t, 'wind_u_ms', NaN, nrow);
    raw.Uv = getcol(t, 'wind_v_ms', NaN, nrow);
    raw.U  = getcol(t, 'wind_ms',   NaN, nrow);
    % radiation / sunshine
    raw.n   = getcol(t, 'sunshine_h', NaN, nrow);  % actual sunshine hours [h]
    raw.Rn  = getcol(t, 'rn_mjm2',  NaN, nrow);    % net radiation [MJ/m2]
    raw.Rin = getcol(t, 'rin_wm2',  NaN, nrow);    % incoming shortwave [W/m2]
    raw.Rli = getcol(t, 'rli_wm2',  NaN, nrow);    % incoming longwave [W/m2]
    % precipitation [mm]
    raw.Precip = getcol(t, 'precip_mm', NaN, nrow);
    % soil temperature observations [degC]
    raw.Ts = getcol(t, 'ts_c', NaN, nrow);
    raw.Tb = getcol(t, 'tb_c', NaN, nrow);
    % surface matric head observation [cm]
    raw.h_SUR = getcol(t, 'h_sur_cm', NaN, nrow);
    % air pressure: hPa -> Pa
    raw.Pg = getcol(t, 'pg_hpa', NaN, nrow) * 100;
    % atmospheric composition
    raw.Ca = getcol(t, 'ca_ppm', NaN, nrow);       % CO2 [ppm]
    % nitrogen deposition + bottom BC: absent columns default to 0 (legacy)
    raw.NH4_dep = getcol(t, 'nh4_dep_mgm2d', 0, nrow);
    raw.NO3_dep = getcol(t, 'no3_dep_mgm2d', 0, nrow);
    raw.NH4_bot = getcol(t, 'nh4_bot_mgl',   0, nrow);
    raw.NO3_bot = getcol(t, 'no3_bot_mgl',   0, nrow);

    raw.n_raw = nrow;
end

% ===========================================================================
function v = getcol(t, name, default, nrow)
%GETCOL  Read a table column by name with a scalar default for missing columns.
%
%   Inputs:
%     t        table  source table
%     name     char   column name
%     default  scalar value used when the column is absent
%     nrow     [-]    expected output row count
%
%   Returns:
%     v        [nrow x 1] double column vector
%
%   Editor: Jizixing Zhu; Date: 2026-06-26.
    if ismember(name, t.Properties.VariableNames)
        v = double(t.(name));
        v = v(:);
    else
        v = default * ones(nrow, 1);
    end
end
