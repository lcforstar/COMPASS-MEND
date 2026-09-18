classdef Climate10dBuffer < handle
%CLIMATE10DBUFFER  10-day rolling-mean climate buffer (auxiliary, shared).
%   Editor: Jizixing Zhu; Date: 2026-06-27.
%   Change: 2026-06-27 Jizixing Zhu -- moved from +core to
%     +vegetation/+integrated/+LUNA (LUNA-owned; decouples the new init chain from
%     +core). Construct as vegetation.integrated.LUNA.Climate10dBuffer(...).
%   Change: 2026-07-01 Jizixing Zhu -- add toStruct/fromStruct for checkpoint/resume
%     (handle class serialized explicitly instead of relying on MATLAB save(handle)).
%
%   Variable class: AUXILIARY (per compass-coding-standard §1).
%     - 不是 prognostic state (可以从 KT-1 forcing 重新累加重建)
%     - 不是 forcing (内部从 KT forcing 计算的派生量)
%     - 不是 diagnostic (会被 LUNA 反向读取作输入, 非只写)
%
%   Encapsulates rolling buffers for Ta, leaf temperature, PAR, RH, CO2, Pg,
%   Oa, boundary-layer resistance and GPP
%   for any downstream module needing 10-day mean climate driver:
%     - LUNA (M2-B.3): Vcmax25/Jmax25 优化输入 (Ali 2016 Eq.6 climate driver)
%     - M2-B.4 turnover + N_storage: 物候触发条件 (future)
%     - SOM (future): 凋落物分解温度调节 (Lloyd & Taylor 1994)
%
%   Replaces ad-hoc Ta_history array in COMPASS.m (line 144). Old:
%     Ta_history(KT) = meteo.Ta;
%     canopy.Tgrowth = mean(Ta_history(max(1, KT - 240 + 1):KT), 'omitnan');
%   New:
%     clim.push(KT, struct('Ta', meteo.Ta, ...));
%     canopy.Tgrowth = clim.getMeans().Ta;
%
%   Usage (in COMPASS main loop):
%       % --- init (before main loop) ---
%       clim = vegetation.integrated.LUNA.Climate10dBuffer(steps_per_day, 10);
%
%       % --- each KT (after forcing read, before ebal) ---
%       clim.push(KT, struct( ...
%           'Ta',      meteo.Ta,                   ...   % [°C]
%           'Tleaf_d', T_leaf_day,                 ...   % [°C] only daytime KT
%           'Tleaf_n', T_leaf_night,               ...   % [°C] only nighttime KT
%           'PAR',     rad.PAR,                    ...   % [mol photons/m2/s]
%           'RH',      meteo.RH,                   ...   % [-] 0..1
%           'CO2',     meteo.CO2,                  ...   % [ppm]
%           'rb',      thermal_rb_prev,            ...   % [s/m]
%           'GPP',     fluxes.GPP,                 ...   % [kg C/m²/s]
%           'is_day',  zenith_sin > 0));                 % [bool]
%
%       % --- when triggering LUNA (every 24 KT) ---
%       if clim.isWarm() && mod(KT, 24) == 0
%           lunaIn  = clim.assembleLunaInput(N_leaf_prev, canopy, dayl_sec, 1);
%           lunaOut = vegetation.integrated.LUNA.computeLuna(lunaIn, lunaState, p);
%       end
%
%   Reference:
%     Ali AA et al. (2016) GMD 9:587-606 — LUNA v1.0, §2.2 "10-day mean climate"
%     CLM5 LunaMod.F90 line 280-310 — par240d_z, t_veg10_day_p, etc. accumulators
%     CLM5 Tech Note Lawrence et al. (2019) Ch.9 §9.3
%
%   Author: Zhu Jizixing 2026-05-20  (M2-B.3 S5.1)
%
%       % --- each KT (after forcing read) ---
%       clim.push(KT, struct( ...
%           'Ta',      meteo.Ta, ...
%           'Tleaf_d', T_leaf_day, ...    % only daytime KT
%           'Tleaf_n', T_leaf_night, ...  % only nighttime KT
%           'PAR',     rad.PAR, ...
%           'RH',      rh_now, ...
%           'CO2',     co2_now, ...
%           'rb',      thermal_rb_prev, ...
%           'GPP',     fluxes.GPP, ...
%           'is_day',  zenith_sin > 0));
%
%       % --- when triggering LUNA (every 24 KT) ---
%       if clim.isWarm() && mod(KT, 24) == 0
%           lunaIn = clim.assembleLunaInput(N_leaf_prev, canopy, dayl_KT);
%           lunaOut = vegetation.integrated.LUNA.computeLuna(lunaIn, lunaState, p);
%       end
%
%   Author: Zhu Jizixing 2026-05-20  (M2-B.3 S5.1)
%
%   Reference:
%     Ali et al. 2016 GMD 9:587 — LUNA 10-day mean climate driver
%     CLM5 LunaMod.F90 line 280-310 — par240d_z, t_veg10_day, etc. accum

    properties
        % --- ring buffers, length = window (= steps_per_day × window_days) ---
        Ta              % [°C, window×1]            air temperature
        Tleaf_d         % [°C, window×1]            daytime leaf temperature (NaN on night KT)
        Tleaf_n         % [°C, window×1]            nighttime leaf temperature (NaN on day KT)
        PAR             % [mol photons/m2/s, window×1]  top-canopy PAR (visible band)
        PAR_max_today   % [mol photons/m2/s, window_days×1]  daily-max PAR (10-day ring)
        RH              % [-, window×1]             relative humidity 0..1
        ea              % [hPa, window×1]           atmospheric vapour pressure
        CO2             % [ppm, window×1]           atmospheric CO2
        Pg              % [Pa, window×1]            atmospheric pressure
        Oa              % [mmol/mol, window×1]      atmospheric O2 mole fraction
        rb              % [s/m, window×1]           boundary-layer resistance
        GPP             % [kg C/(m²·s), window×1]   instantaneous GPP per KT

        % --- ring pointers ---
        head            % [-]   next write index 1..window
        filled          % [-]   # entries written so far (≤ window)

        % --- daily accumulation ---
        steps_per_day   % [-]                       e.g. 24 for hourly forcing
        gpp_today       % [kg C/m²]                 accumulating since 00:00 today
        gpp_yesterday   % [kg C/m²/day]             total from last full day (LUNA decay test)
        par_max_today_buf  % [mol photons/m2/s]      running max within current day
        day_count       % [-]                       days elapsed since buffer init
        kt_within_day   % [-]                       KT index within current day 0..(steps_per_day-1)

        % --- config (immutable after construction) ---
        window          % [-]                       = window_days × steps_per_day, typically 240
        window_days     % [-]                       default 10 days (CLM5 LUNA convention)
    end

    methods
        function obj = Climate10dBuffer(steps_per_day, window_days)
            %CLIMATE10DBUFFER  Construct empty buffer.
            if nargin < 1 || isempty(steps_per_day), steps_per_day = 24; end
            if nargin < 2 || isempty(window_days),   window_days   = 10; end
            obj.steps_per_day = steps_per_day;
            obj.window_days   = window_days;
            obj.window        = steps_per_day * window_days;

            obj.Ta            = nan(obj.window, 1);
            obj.Tleaf_d       = nan(obj.window, 1);
            obj.Tleaf_n       = nan(obj.window, 1);
            obj.PAR           = nan(obj.window, 1);
            obj.RH            = nan(obj.window, 1);
            obj.ea            = nan(obj.window, 1);
            obj.CO2           = nan(obj.window, 1);
            obj.Pg            = nan(obj.window, 1);
            obj.Oa            = nan(obj.window, 1);
            obj.rb            = nan(obj.window, 1);
            obj.GPP           = nan(obj.window, 1);
            obj.PAR_max_today = nan(window_days, 1);

            obj.head = 1;
            obj.filled = 0;
            obj.gpp_today = 0;
            obj.gpp_yesterday = 0;
            obj.par_max_today_buf = 0;
            obj.day_count = 0;
            obj.kt_within_day = 0;
        end

        function push(obj, ~, f)
            %PUSH  Append one KT of climate.
            %   f struct, fields (use NaN if not applicable):
            %     .Ta .Tleaf_d .Tleaf_n .PAR .RH .ea .CO2 .rb .GPP .is_day
            obj.Ta(obj.head)      = getOr(f, 'Ta', NaN);
            obj.PAR(obj.head)     = getOr(f, 'PAR', NaN);
            obj.RH(obj.head)      = getOr(f, 'RH', NaN);
            obj.ea(obj.head)      = getOr(f, 'ea', NaN);
            obj.CO2(obj.head)     = getOr(f, 'CO2', NaN);
            obj.Pg(obj.head)      = getOr(f, 'Pg', NaN);
            obj.Oa(obj.head)      = getOr(f, 'Oa', NaN);
            obj.rb(obj.head)      = getOr(f, 'rb', NaN);
            obj.GPP(obj.head)     = getOr(f, 'GPP', 0);

            % Day/night split for leaf T
            is_day = getOr(f, 'is_day', true);
            if is_day
                obj.Tleaf_d(obj.head) = getOr(f, 'Tleaf_d', NaN);
                obj.Tleaf_n(obj.head) = NaN;
            else
                obj.Tleaf_d(obj.head) = NaN;
                obj.Tleaf_n(obj.head) = getOr(f, 'Tleaf_n', NaN);
            end

            % Daily accumulation
            gpp_kt = getOr(f, 'GPP', 0);
            if isfinite(gpp_kt) && gpp_kt > 0
                obj.gpp_today = obj.gpp_today + gpp_kt * 3600;   % kg C/m²/s × s = kg C/m²
            end
            par_kt = getOr(f, 'PAR', 0);
            if isfinite(par_kt) && par_kt > obj.par_max_today_buf
                obj.par_max_today_buf = par_kt;
            end

            obj.kt_within_day = obj.kt_within_day + 1;
            if obj.kt_within_day >= obj.steps_per_day
                obj.day_count = obj.day_count + 1;
                obj.kt_within_day = 0;

                % Shift PAR_max_today ring (window_days entries)
                idx = mod(obj.day_count - 1, obj.window_days) + 1;
                obj.PAR_max_today(idx) = obj.par_max_today_buf;
                obj.par_max_today_buf = 0;

                obj.gpp_yesterday = obj.gpp_today;
                obj.gpp_today = 0;
            end

            % Advance ring
            obj.head = mod(obj.head, obj.window) + 1;
            obj.filled = min(obj.filled + 1, obj.window);
        end

        function m = getMeans(obj)
            %GETMEANS  Return struct of 10-day means (uses available data, NaN if 0 filled).
            valid = obj.values_so_far();   % indices to use
            m.Ta      = nanmeanSafe(obj.Ta(valid));
            m.Tleaf_d = nanmeanSafe(obj.Tleaf_d(valid));
            m.Tleaf_n = nanmeanSafe(obj.Tleaf_n(valid));
            m.PAR     = nanmeanSafe(obj.PAR(valid));
            m.RH      = nanmeanSafe(obj.RH(valid));
            m.ea      = nanmeanSafe(obj.ea(valid));
            m.CO2     = nanmeanSafe(obj.CO2(valid));
            m.Pg      = nanmeanSafe(obj.Pg(valid));
            m.Oa      = nanmeanSafe(obj.Oa(valid));
            m.rb      = nanmeanSafe(obj.rb(valid));
            m.PAR_max = nanmeanSafe(obj.PAR_max_today);   % 10-day mean of daily maxes
            m.n_filled = obj.filled;
        end

        function gpp = getGPPYesterday(obj)
            %GETGPPYESTERDAY  Total accumulated GPP from previous full day [kg C/m²/day].
            gpp = obj.gpp_yesterday;
        end

        function tf = isWarm(obj)
            %ISWARM  Returns true once buffer has at least 1 full day of data.
            %   LUNA falls back to PFT-table defaults if isWarm() == false.
            tf = obj.filled >= obj.steps_per_day;
        end

        function tf = isFull(obj)
            %ISFULL  Returns true once buffer has 10 full days (window saturated).
            tf = obj.filled >= obj.window;
        end

        function lunaIn = assembleLunaInput(obj, ~, canopy, dayl_sec, c3psn)
            %ASSEMBLELUNAINPUT  Build LUNA input struct from 10-day means.
            %   canopy struct: .LAI .slatop .nrad .tlai_z .rabsorb .rhol .taul
            %   dayl_sec   [s]   current daylight length
            %   c3psn      [-]   1 for C3, 0 for C4
            m = obj.getMeans();

            % [2026-05-20 lnc_opt=.false. fix]
            % CLM5 default (PhotosynthesisMod.F90:1321):
            %   lnc(p) = 1 / (slatop × leafcn)        ← PFT constant
            % This decouples LUNA from prognostic N_leaf pool, avoiding the
            % mid-year start problem (CLM5 cold-start leafn=0, takes a full
            % year of phenology to ramp up to summer values).
            % The prognostic N_leaf pool is still used for root uptake demand
            % and NEE, just not as LUNA's lnc input.
            %
            % Original (lnc_opt=.true. path, kept for reference):
            %   LAI_safe = max(canopy.LAI, 0.1);
            %   lunaIn.lnc = N_leaf_per_m2_ground / LAI_safe;
            slatop = getOr2(canopy, 'slatop', 0.0125);            % [m²/g C]
            leafcn = getOr2(canopy, 'leafcn', 25.0);              % [g C/g N]
            lunaIn.lnc = 1.0 / (slatop * leafcn);                 % [g N/m² leaf]

            lunaIn.tlai = canopy.LAI;
            lunaIn.nrad = getOr2(canopy, 'nrad', 1);
            % [2026-05-20 PAR unit fix] rad.PAR in COMPASS is mol/m²/s
            % (RTMo.m line 510); LUNA expects W/m² visible band (LunaMod.F90:312
            % `PARTop = RadTop * 4.6` converts W/m² → μmol photons/m²/s).
            % So: rad.PAR [mol/m²/s] × 1e6 → μmol/m²/s, then /4.6 → W/m²
            par_W_m2 = m.PAR     * 1.0e6 / 4.6;   % [W/m² visible]
            parmx_W_m2 = m.PAR_max * 1.0e6 / 4.6;
            if lunaIn.nrad == 1
                lunaIn.tlai_z   = canopy.LAI;
                lunaIn.par240d_z = par_W_m2;
                lunaIn.par240x_z = parmx_W_m2;
            else
                lunaIn.tlai_z   = getOr2(canopy, 'tlai_z', ones(lunaIn.nrad,1) * canopy.LAI/lunaIn.nrad);
                lunaIn.par240d_z = ones(lunaIn.nrad,1) * par_W_m2;
                lunaIn.par240x_z = ones(lunaIn.nrad,1) * parmx_W_m2;
            end
            lunaIn.slatop  = getOr2(canopy, 'slatop', 0.0125);   % BDT temperate default
            lunaIn.rabsorb = getOr2(canopy, 'rabsorb', 0.85);

            % climate 10-day means (already filled by getMeans)
            phys = equations.physicalConstants();
            pbot10 = m.Pg;
            if ~isfinite(pbot10) || pbot10 <= 0
                pbot10 = phys.Pstd;
            end
            Oa10 = m.Oa;
            if ~isfinite(Oa10) || Oa10 <= 0
                Oa10 = phys.Oa_mmolmol_default;
            end

            lunaIn.forc_pbot10 = pbot10;        % [Pa] atmospheric pressure
            lunaIn.relh10    = max(0.25, min(1.0, m.RH));
            ea10_Pa = m.ea * 100.0;                    % [Pa] water vapour pressure
            if ~isfinite(ea10_Pa) || ea10_Pa < 0
                ea10_Pa = 0.0;
            end
            pdry10 = max(pbot10 - ea10_Pa, 1.0);       % [Pa] dry-air pressure
            lunaIn.CO2a10    = m.CO2 * 1e-6 * pdry10;  % [Pa] CO2 partial pressure
            lunaIn.O2a10     = Oa10 * 1e-3 * pdry10;   % [Pa] O2 partial pressure
            lunaIn.rb10      = max(m.rb, 0.1);
            lunaIn.hourpd    = dayl_sec / 3600.0;        % [hr]
            lunaIn.tair10    = m.Ta;
            lunaIn.tleafd10  = nanmeanSafe([m.Tleaf_d, m.Ta]);  % fallback to Ta if no daytime data
            lunaIn.tleafn10  = nanmeanSafe([m.Tleaf_n, m.Ta]);
            % NOTE: CLM5 LunaMod.F90 line 300 has tleaf10 = tleafd10 (apparent bug
            % in CLM5 but we follow 1:1 for max_daily_pchg = EnzTurnoverTFactor consistency).
            % Audit 2026-05-20 Bug #1.
            lunaIn.tleaf10   = lunaIn.tleafd10;
            lunaIn.gpp_day   = obj.getGPPYesterday();
            lunaIn.c3psn     = c3psn;
        end

        function s = toStruct(obj)
            %TOSTRUCT  Serialize all buffer state to a plain struct (checkpoint/resume).
            %   Copies every public property so a saved struct round-trips through
            %   fromStruct without depending on MATLAB's save(handle) behavior.
            % Editor: Jizixing Zhu; Date: 2026-07-01.
            s = struct();
            p = properties(obj);
            for i = 1:numel(p)
                s.(p{i}) = obj.(p{i});
            end
        end
    end

    methods (Static)
        function obj = fromStruct(s)
            %FROMSTRUCT  Rebuild a Climate10dBuffer from a toStruct snapshot.
            %   Reconstructs the object with its saved window sizing, then restores
            %   every saved property (ring buffers, pointers, daily accumulators).
            % Editor: Jizixing Zhu; Date: 2026-07-01.
            obj = vegetation.integrated.LUNA.Climate10dBuffer(s.steps_per_day, s.window_days);
            p = properties(obj);
            for i = 1:numel(p)
                if isfield(s, p{i})
                    obj.(p{i}) = s.(p{i});
                end
            end
        end
    end

    methods (Access = private)
        function idx = values_so_far(obj)
            if obj.filled >= obj.window
                idx = 1:obj.window;
            else
                idx = 1:obj.filled;
            end
        end
    end
end

% --- local helpers ---
function v = getOr(s, fname, default)
    if isstruct(s) && isfield(s, fname) && ~isempty(s.(fname))
        v = s.(fname);
    else
        v = default;
    end
end

function v = getOr2(s, fname, default)
    if isstruct(s) && isfield(s, fname) && ~isempty(s.(fname))
        v = s.(fname);
    elseif isobject(s) && isprop(s, fname) && ~isempty(s.(fname))
        v = s.(fname);
    else
        v = default;
    end
end

function m = nanmeanSafe(x)
    if isempty(x), m = NaN; return; end
    x = x(:);
    valid = isfinite(x);
    if ~any(valid)
        m = NaN;
    else
        m = mean(x(valid));
    end
end
