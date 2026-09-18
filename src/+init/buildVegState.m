function vegState = buildVegState(built, veg, raw)
%BUILDVEGSTATE  Prognostic vegetation state initial values (built.vegetation_state).
%
%   vegState = init.buildVegState(built, veg) builds the per-day vegetation state
%   that the main loop advances: the phenology driver (V3) and the LAI/hc/SAI
%   series + per-day LAI state init (V4).
%
%   The driver holds the forcing-only quantities that can be precomputed once
%   (the GDD heat-unit curve from daily temperature). The LAI itself is NOT
%   precomputed for model mode -- it is an ALAI_MIN placeholder here and the main
%   loop overwrites it each day via the hydraulic-LAI update (not done here).
%
%   Inputs:
%     built  struct  needs forcing.raw.Tmax/Tmin or hourly forcing.Ta + control.time
%     veg    struct  built.parameters.veg (uses canopy.* + phenology.*)
%     raw    struct  optional raw.vegetation table with DOY, LAI, hc, SAI data
%
%   Returns vegState:
%     vegState.phen_driver  struct  phen_CHU/HU/HUI/dHUF + phen_Tavg/Tmin/Tmax_daily
%                                   (empty when LAI_mode ~= 'model')      [V3]
%     vegState.lai_state    struct  LAI/hc/SAI step series + per-day LAI state
%                                   (lai_day_prev/dormant/LAI_D/HUI_D/yr,
%                                   lai_model_current, wt_day_min)         [V4]
%
%   Reference: legacy VegetationData.init_phenology / to_hourly / init_day_state.
% Editor: Jizixing Zhu; Date: 2026-06-27.
% Change: 2026-07-16 Jizixing Zhu -- name the measured canopy input explicitly.
% Change: 2026-08-24 Jizixing Zhu -- retain the temperature phenology driver
%   for PlantCN even when canopy LAI is observation-driven.

    if nargin < 3
        raw = struct();
    end
    vegState = struct();
    vegState.phen_driver = buildPhenDriver(built, veg);        % V3
    vegState.lai_state   = buildLaiState(built, veg, raw);     % V4
end

% ===========================================================================
function ls = buildLaiState(built, veg, raw)
% V4: LAI/hc/SAI step series + per-day LAI state init. Reproduces
% VegetationData.to_hourly (series, length Nmsrmn) and init_day_state (state),
% per the agreed boundary: model LAI is only an ALAI_MIN placeholder (the main
% loop updates it daily); data mode interpolates measured DOY data from
% raw.vegetation (canopy_measurements.csv: DOY, LAI, hc, Vcmax, Cab, SAI).
% Ref: legacy/src/+core/VegetationData.m:372-482, :730-745.
% Editor: Jizixing Zhu; Date: 2026-06-27.
    Nh = built.control.time.Nmsrmn;
    c  = veg.canopy;
    ls = struct();
    obs = parseVegetationObservations(raw);
    if needsData(c) && isempty(obs.doy)
        error('buildVegState:missingVegetationData', ...
            ['LAI/hc/SAI data mode needs a vegetation observation table. ', ...
             'Set vegetation.vegetation_data_file or provide canopy_measurements.csv in the site folder.']);
    end
    doyClamped = [];
    if needsData(c)
        doyClamped = simulationDoyClamped(built.control.time, obs.doy);
    end

    % --- LAI series (constant value / model placeholder / data interp) ---
    switch lower(c.LAI_mode)
        case {'constant', 'function'}
            ls.LAI = c.LAI_constant * ones(Nh, 1);
        case 'model'
            ls.LAI = veg.phenology.ALAI_MIN * ones(Nh, 1);   % placeholder; main loop overwrites
        case 'data'
            ls.LAI = max(0, interpolateObservedSeries(obs, 'LAI', doyClamped));
        otherwise
            error('buildVegState:badLAImode', 'unknown LAI_mode "%s"', c.LAI_mode);
    end

    % --- hc series ---
    switch lower(c.hc_mode)
        case {'constant', 'function', 'model'}
            ls.hc = c.canopy_height * ones(Nh, 1);
        case 'data'
            ls.hc = max(0.01, interpolateObservedSeries(obs, 'hc', doyClamped));
        otherwise
            error('buildVegState:badHCmode', 'unknown hc_mode "%s"', c.hc_mode);
    end

    % --- SAI series ---
    switch lower(c.SAI_mode)
        case {'constant', 'function', 'model'}
            ls.SAI = c.SAI_constant * ones(Nh, 1);
        case 'data'
            ls.SAI = max(0, interpolateObservedSeries(obs, 'SAI', doyClamped));
        otherwise
            error('buildVegState:badSAImode', 'unknown SAI_mode "%s"', c.SAI_mode);
    end

    if ~isempty(obs.doy)
        ls.observed = obs;
    end

    % --- per-day LAI state init (== init_day_state) ---
    ls.lai_day_dormant = false;
    ls.lai_day_LAI_D   = 0;
    ls.lai_day_HUI_D   = veg.phenology.DLAI;
    ls.lai_day_yr      = 0;
    if ~isempty(ls.LAI)
        ls.lai_day_prev = ls.LAI(1);
    else
        ls.lai_day_prev = veg.phenology.ALAI_MIN;
    end
    ls.lai_model_current = ls.lai_day_prev;
    ls.wt_day_min = 1.0;   % [-] daily min hydraulic stress, reset after each daily LAI update
end

% ===========================================================================
function tf = needsData(canopy)
% True when any vegetation-structure series is requested from observations.
    tf = strcmpi(canopy.LAI_mode, 'data') || strcmpi(canopy.hc_mode, 'data') || ...
        strcmpi(canopy.SAI_mode, 'data');
end

% ===========================================================================
function obs = parseVegetationObservations(raw)
% Parse optional DOY-indexed LAI/hc/SAI observations from raw.vegetation.
    obs = struct('doy', [], 'LAI', [], 'hc', [], 'SAI', []);
    if ~isstruct(raw) || ~isfield(raw, 'vegetation') || ~istable(raw.vegetation) || isempty(raw.vegetation)
        return;
    end

    tbl = raw.vegetation;
    obs.doy = readNumericColumn(tbl, {'doy','DOY','day_of_year','DayOfYear'}, 1, true);
    obs.LAI = readNumericColumn(tbl, {'LAI','lai'}, 2, false);
    obs.hc  = readNumericColumn(tbl, {'hc','HC','canopy_height','canopyHeight'}, 3, false);
    obs.SAI = readNumericColumn(tbl, {'SAI','sai'}, 6, false);

    ok = isfinite(obs.doy);
    obs.doy = obs.doy(ok);
    for name = {'LAI','hc','SAI'}
        v = obs.(name{1});
        if ~isempty(v)
            obs.(name{1}) = v(ok);
        end
    end
    [obs.doy, order] = sort(obs.doy);
    for name = {'LAI','hc','SAI'}
        v = obs.(name{1});
        if ~isempty(v)
            obs.(name{1}) = v(order);
        end
    end
    [obs.doy, ia] = unique(obs.doy, 'stable');
    for name = {'LAI','hc','SAI'}
        v = obs.(name{1});
        if ~isempty(v)
            obs.(name{1}) = v(ia);
        end
    end
end

% ===========================================================================
function values = readNumericColumn(tbl, names, fallbackIndex, required)
% Read a numeric table column by name, with a legacy positional fallback.
    values = [];
    vars = tbl.Properties.VariableNames;
    idx = [];
    for i = 1:numel(names)
        hit = find(strcmpi(vars, names{i}), 1);
        if ~isempty(hit)
            idx = hit;
            break;
        end
    end
    if isempty(idx) && fallbackIndex <= width(tbl)
        idx = fallbackIndex;
    end
    if isempty(idx)
        if required
            error('buildVegState:missingVegetationColumn', ...
                'Vegetation observation table is missing required DOY column.');
        end
        return;
    end
    values = tbl.(vars{idx});
    if iscell(values) || isstring(values)
        values = str2double(string(values));
    end
    values = double(values(:));
end

% ===========================================================================
function doyClamped = simulationDoyClamped(time, obsDoy)
% Build the simulation DOY vector and clamp to the observed seasonal support.
    Nh = time.Nmsrmn;
    timestamps = time.start_time + seconds((0:Nh - 1)' * time.Delt_t);
    doySim = day(timestamps, 'dayofyear') + ...
        (hour(timestamps) + minute(timestamps) / 60 + second(timestamps) / 3600) / 24;
    doyClamped = max(min(obsDoy), min(max(obsDoy), doySim));
end

% ===========================================================================
function series = interpolateObservedSeries(obs, name, doyClamped)
% Interpolate one observed vegetation-structure series to model time steps.
    values = obs.(name);
    if isempty(values) || numel(values) ~= numel(obs.doy)
        error('buildVegState:missingVegetationSeries', ...
            'Vegetation observation table needs at least two finite %s points.', name);
    end
    ok = isfinite(obs.doy) & isfinite(values);
    if nnz(ok) < 2
        error('buildVegState:missingVegetationSeries', ...
            'Vegetation observation table needs at least two finite %s points.', name);
    end
    series = interp1(obs.doy(ok), values(ok), doyClamped, 'linear');
end

% ===========================================================================
function pd = buildPhenDriver(built, veg)
% V3: daily temperature and GDD driver for canopy and PlantCN phenology. Reuses the same static
% method as VegetationData.init_phenology so the result matches element-for-
% element when daily Tmax/Tmin are supplied. Observation-driven LAI still needs
% this driver when PlantCN is active: organ onset, allocation, turnover, and
% litter are prognostic even though the canopy area is prescribed. Daily
% Tmax/Tmin are read from the raw forcing layer when available; otherwise they
% are derived from hourly air temperature by clock day.
% Editor: Jizixing Zhu; Date: 2026-06-27.
% Change: 2026-08-24 Jizixing Zhu -- decouple PlantCN phenology forcing from LAI mode.
    pd = struct();
    plantCnOn = isfield(built, 'control') && isfield(built.control, 'switches') ...
        && isfield(built.control.switches, 'lPlantCN') ...
        && built.control.switches.lPlantCN == 1;
    if ~strcmpi(veg.canopy.LAI_mode, 'model') && ~plantCnOn
        return;
    end
    [Tmax, Tmin] = phenologyDailyTemperatureRange(built);
    ph   = veg.phenology;

    [pd.phen_CHU, pd.phen_HU, pd.phen_HUI, ~, ~, ~, pd.phen_dHUF] = ...
        vegetation.phenology.CalculateHeatUnit(Tmax, Tmin, ph.sowing_day, ph.T_BASE, ph.PHU, ...
            ph.FRGW1, ph.LAIMX1, ph.FRGW2, ph.LAIMX2);

    pd.phen_Tavg_daily = (Tmax + Tmin) / 2;     % daily mean T (phenology dispatcher)
    pd.phen_Tmin_daily = Tmin;                  % daily Tmin (Cesaraccio chilling)
    pd.phen_Tmax_daily = Tmax;                  % daily Tmax
end

% ===========================================================================
function [Tmax, Tmin] = phenologyDailyTemperatureRange(built)
% Return daily Tmax/Tmin used by the heat-unit LAI driver.
%
%   Inputs:
%     built [in] struct  forcing container with raw daily Tmax/Tmin or hourly Ta
%
%   Outputs:
%     Tmax [degC, nday x 1] daily maximum air temperature
%     Tmin [degC, nday x 1] daily minimum air temperature
%
%   Method:
%     Prefer measured daily Tmax/Tmin when they are present. For land-model
%     forcing at hourly/sub-hourly cadence, derive Tmax/Tmin from forcing.hourly.Ta
%     grouped by forcing.timeAxis.clockDayIndex. The derived daily range is a
%     forcing auxiliary cache, not a prognostic state.
%
%   Editor: Jizixing Zhu; Date: 2026-07-08.
    rawTmax = [];
    rawTmin = [];
    if isfield(built, 'forcing') && isfield(built.forcing, 'raw') ...
            && isfield(built.forcing.raw, 'Tmax') && isfield(built.forcing.raw, 'Tmin')
        rawTmax = built.forcing.raw.Tmax(:);
        rawTmin = built.forcing.raw.Tmin(:);
    end

    [hourlyTmax, hourlyTmin] = deriveDailyRangeFromHourlyTa(built);
    [Tmax, Tmin] = mergeDailyTemperatureRange(rawTmax, rawTmin, hourlyTmax, hourlyTmin);
    ok = isfinite(Tmax) & isfinite(Tmin);
    if ~any(ok)
        error('buildVegState:missingPhenologyTemperature', ...
            ['LAI_mode="model" needs daily Tmax/Tmin or hourly air temperature Ta ', ...
             'to build the heat-unit LAI driver.']);
    end

    Tmax = fillMissingDailyTemperature(Tmax);
    Tmin = fillMissingDailyTemperature(Tmin);
end

% ===========================================================================
function [Tmax, Tmin] = deriveDailyRangeFromHourlyTa(built)
% Aggregate native-cadence air temperature to clock-day Tmax/Tmin.
    Tmax = [];
    Tmin = [];
    if ~isfield(built, 'forcing') || ~isfield(built.forcing, 'hourly') ...
            || ~isfield(built.forcing.hourly, 'Ta') || ~isfield(built.forcing, 'timeAxis') ...
            || ~isfield(built.forcing.timeAxis, 'clockDayIndex')
        return
    end

    Ta = built.forcing.hourly.Ta(:);
    dayIndex = built.forcing.timeAxis.clockDayIndex(:);
    n = min(numel(Ta), numel(dayIndex));
    Ta = Ta(1:n);
    dayIndex = dayIndex(1:n);
    ok = isfinite(Ta) & isfinite(dayIndex) & dayIndex >= 1;
    if ~any(ok)
        return
    end

    dayIndex = round(dayIndex(ok));
    Ta = Ta(ok);
    nDay = max(dayIndex);
    Tmax = accumarray(dayIndex, Ta, [nDay, 1], @max, NaN);
    Tmin = accumarray(dayIndex, Ta, [nDay, 1], @min, NaN);
end

% ===========================================================================
function [Tmax, Tmin] = mergeDailyTemperatureRange(rawTmax, rawTmin, hourlyTmax, hourlyTmin)
% Prefer daily measured extrema, filling gaps from hourly-derived extrema.
    rawIdx = [];
    nRaw = min(numel(rawTmax), numel(rawTmin));
    if nRaw > 0
        rawOk = isfinite(rawTmax(1:nRaw)) & isfinite(rawTmin(1:nRaw));
        rawIdx = find(rawOk);
    end

    nRawDay = 0;
    if ~isempty(rawIdx)
        nRawDay = max(rawIdx);
    end
    nDay = max([nRawDay, numel(hourlyTmax), numel(hourlyTmin), 0]);
    Tmax = NaN(nDay, 1);
    Tmin = NaN(nDay, 1);

    if ~isempty(hourlyTmax)
        Tmax(1:numel(hourlyTmax)) = hourlyTmax(:);
    end
    if ~isempty(hourlyTmin)
        Tmin(1:numel(hourlyTmin)) = hourlyTmin(:);
    end

    if ~isempty(rawIdx)
        Tmax(rawIdx) = rawTmax(rawIdx);
        Tmin(rawIdx) = rawTmin(rawIdx);
    end
end

% ===========================================================================
function x = fillMissingDailyTemperature(x)
% Fill isolated missing daily temperatures without inventing a seasonal curve.
    good = isfinite(x);
    if all(good)
        return
    end
    if ~any(good)
        return
    end

    idx = (1:numel(x))';
    x(~good) = interp1(idx(good), x(good), idx(~good), 'linear', 'extrap');
end
