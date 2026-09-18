function phen = advanceDay(phen, phenDriver, pftType, ts, startTime, latDeg, sitePhenology)
%ADVANCEDAY  Advance the CLM5/LPJ-GUESS phenology one day (leaf-on fraction phen.now).
%
%   phen = vegetation.integrated.phenology.advanceDay(phen, phenDriver, pftType, ts, ...
%       startTime, latDeg) advances the persistent phenology state by ONE day: it
%   assembles the daily forcing (mean/min/max temperature, 30-day running mean,
%   daylength) and PFT phenology parameters, then calls phenologyDispatch to update
%   the CLM5 leaf-on/off accumulators and phen.now in [0,1]. It is the pure per-day
%   core migrated from COMPASS.m:1298-1367 (the caller gates it to the first outer
%   step of each day and commits the result only on an accepted step).
%
%   Inputs:
%     phen        [in]  phenology state struct (Y.vegetation.phenology): .state (PhenologyState), .now [-],
%                       .year [-], .year_mtemp_max [degC] (running, reset per year)
%     phenDriver  [in]  built.vegetation_state.phen_driver: .phen_Tavg_daily /
%                       .phen_Tmin_daily / .phen_Tmax_daily [degC, per model day]
%     pftType     [in]  PFT name (built.parameters.veg.PFT_type, e.g. 'Apple')
%     ts          [in]  datetime  end-of-step clock time (day-of-year / year source)
%     startTime   [in]  datetime  simulation start (multi-year day-index anchor)
%     latDeg      [in]  [degN]    site latitude (daylength via solar declination)
%
%   Returns:
%     phen        updated struct (.state/.now/.year/.year_mtemp_max)
%
%   Note (no stage-4 dependency): water_stress is passed as 0, matching the legacy
%   call (COMPASS.m:1365); the stress_decid branch's soil-water coupling is not wired.
%
%   Reference: COMPASS.m:1298-1367; phenologyDispatch (CLM5 CNPhenologyMod.F90:241-342).
%   Editor: Jizixing Zhu; Date: 2026-07-02.
%   Change: 2026-07-05 Jizixing Zhu — doc: phen arg is Y.vegetation.phenology (was runtime.phenology); ownership migration increment 5. Body unchanged.
%   Change: 2026-08-24 Jizixing Zhu -- index daily forcing by elapsed model
%           day so simulations starting after January 1 remain clock-aligned.
    if nargin < 7 || isempty(sitePhenology), sitePhenology = struct(); end
    doy    = day(ts, 'dayofyear');
    yr     = year(ts);
    dayIdx = floor(days(ts - startTime)) + 1;       % one-based elapsed forcing day

    % Year boundary: reset the running max monthly-mean temperature (season_decid).
    if yr ~= phen.year
        phen.year = yr;
        phen.year_mtemp_max = -999;
    end

    [Tavg, Tmin, Tmax] = dailyTemps(phenDriver, dayIdx);
    mtemp = runningMonthlyMean(phenDriver, dayIdx, Tavg);   % 30-day running mean
    phen.year_mtemp_max = max(phen.year_mtemp_max, mtemp);

    forcing = struct('T_avg', Tavg, 'Tmin', Tmin, 'Tmax', Tmax, ...
                     'dayl_h', dayLength(latDeg, doy), 'mtemp', mtemp, ...
                     'latitude_deg', latDeg, 'water_stress', 0);

    [phen.state, ~] = vegetation.integrated.phenology.phenologyDispatch( ...
        phen.state, phenPftParams(pftType, sitePhenology, phenDriver, startTime, latDeg), ...
        forcing, doy, phen.year_mtemp_max);
    phen.now = phen.state.phen;
end

% ---------------------------------------------------------------------------
function [Tavg, Tmin, Tmax] = dailyTemps(pd, dayIdx)
% Daily mean/min/max temperature for absolute model day dayIdx (legacy 1309-1321).
% Min/Max default to the mean when only a mean series is present.
% Editor: Jizixing Zhu; Date: 2026-07-02.
    Tavg = 15;                                   % [degC] fallback (legacy 1313)
    if hasDay(pd, 'phen_Tavg_daily', dayIdx)
        Tavg = pd.phen_Tavg_daily(dayIdx);
    end
    Tmin = Tavg;  Tmax = Tavg;
    if hasDay(pd, 'phen_Tmin_daily', dayIdx) && hasDay(pd, 'phen_Tmax_daily', dayIdx)
        Tmin = pd.phen_Tmin_daily(dayIdx);
        Tmax = pd.phen_Tmax_daily(dayIdx);
    end
end

% ---------------------------------------------------------------------------
function mtemp = runningMonthlyMean(pd, dayIdx, Tavg)
% 30-day trailing mean of daily mean temperature (legacy 1322-1327); Tavg fallback.
% Editor: Jizixing Zhu; Date: 2026-07-02.
    mtemp = Tavg;
    if isfield(pd, 'phen_Tavg_daily') && ~isempty(pd.phen_Tavg_daily)
        lastDay = min(dayIdx, numel(pd.phen_Tavg_daily));
        windowStart = max(1, lastDay - 29);
        mtemp = mean(pd.phen_Tavg_daily(windowStart:lastDay));
    end
end

% ---------------------------------------------------------------------------
function h = dayLength(latDeg, doy)
% Daylength [h] from solar declination (legacy 1331-1335).
% Editor: Jizixing Zhu; Date: 2026-07-02.
    latRad = deg2rad(latDeg);
    decl   = deg2rad(23.45) * sin(2*pi/365 * (doy - 81));
    cosOmega = min(1, max(-1, -tan(latRad) * tan(decl)));
    h = 24/pi * acos(cosOmega);
end

% ---------------------------------------------------------------------------
function p = phenPftParams(pftType, sitePhenology, phenDriver, startTime, latDeg)
% PFT phenology parameters from the PFT table (legacy 1340-1357).
% Editor: Jizixing Zhu; Date: 2026-07-02.
    g = @(name, varargin) vegetation.getPFTValue(pftType, name, varargin{:});
    p = struct( ...
        'evergreen',                g('evergreen'), ...
        'season_decid',             g('season_decid'), ...
        'stress_decid',             g('stress_decid'), ...
        'crop_phenology',           g('crop_phenology'), ...
        'k_chilla',                 g('k_chilla'), ...
        'k_chillb',                 g('k_chillb'), ...
        'k_chillk',                 g('k_chillk'), ...
        'phengdd5ramp',             g('phengdd5ramp'), ...
        'crit_dayl_h',              g('crit_dayl_h'), ...
        'ndays_on',                 g('ndays_on'), ...
        'ndays_off',                g('ndays_off'), ...
        'GDD_bloom_after_budburst', g('GDD_bloom_after_budburst'), ...
        'TSum_fruit',               g('TSum_fruit'), ...
        'harvest_DS',               g('harvest_DS'), ...
        'use_cesaraccio_chill',     g('use_cesaraccio_chill', 0), ...
        'crequ',                    g('crequ', 0), ...
        'baset',                    g('baset', 4), ...
        'winter_wheat_cn',          g('winter_wheat_cn', 0), ...
        'arooti',                   g('arooti', 0), ...
        'arootf',                   g('arootf', 0), ...
        'astemf',                   g('astemf', 0), ...
        'allconss',                 g('allconss', 1), ...
        'bfact',                    g('bfact', 0), ...
        'fleafi',                   g('fleafi', 0), ...
        'lfemerg',                  g('lfemerg', 0), ...
        'grnfill',                  g('grnfill', 1), ...
        'mxmat',                    g('mxmat', 365), ...
        'sowing_day',               readField(sitePhenology, 'sowing_day', 1), ...
        'gdd1020',                  cropGdd1020(phenDriver, startTime, latDeg));
end

function value = readField(s, name, fallback)
value = fallback;
if isstruct(s) && isfield(s, name) && ~isempty(s.(name)), value = s.(name); end
end

function gdd1020 = cropGdd1020(pd, startTime, latDeg)
% CoLM climate accumulator: mean April-September (NH) or Oct-March (SH)
% sum of max(Tref-10 degC,0), expressed in degC d.
gdd1020 = 0;
if ~isfield(pd, 'phen_Tavg_daily') || isempty(pd.phen_Tavg_daily), return; end
t = pd.phen_Tavg_daily(:);
dates = startTime + days((0:numel(t)-1)');
if latDeg >= 0
    active = month(dates) >= 4 & month(dates) <= 9;
else
    active = month(dates) > 9 | month(dates) < 4;
end
years = unique(year(dates(active)));
annual = zeros(0,1);
for k = 1:numel(years)
    take = active & year(dates) == years(k);
    if nnz(take) >= 180
        annual(end+1,1) = sum(max(0, t(take) - 10)); %#ok<AGROW>
    end
end
if isempty(annual)
    gdd1020 = sum(max(0, t(active) - 10));
else
    gdd1020 = mean(annual);
end
end

% ---------------------------------------------------------------------------
function tf = hasDay(pd, name, dayIdx)
% True when field `name` exists, is non-empty, and covers absolute day dayIdx.
% Editor: Jizixing Zhu; Date: 2026-07-02.
    tf = isfield(pd, name) && ~isempty(pd.(name)) && dayIdx >= 1 && dayIdx <= numel(pd.(name));
end
