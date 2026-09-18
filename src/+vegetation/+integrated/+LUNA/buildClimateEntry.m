function climate = buildClimateEntry(model, kt)
%BUILDCLIMATEENTRY  Build one accept-only LUNA climate-buffer entry.
%
%   Inputs:
%     model [in] Model shell after the canopy energy balance (scope_step has
%                meteo, thermal, rad, fluxes, angles).
%     kt    [-]  One-based outer-step index into forcing.hourly.
%
%   Return:
%     climate [struct] One hourly sample for Climate10dBuffer:
%       .Ta       [degC]  Air temperature
%       .Tleaf_d  [degC]  Daytime leaf temperature (NaN at night)
%       .Tleaf_n  [degC]  Nighttime leaf temperature (NaN by day)
%       .PAR      [umol m-2 s-1 or model rad unit] Instantaneous PAR
%       .RH       [-]     Relative humidity in [0,1]
%       .ea       [hPa]   Air vapor pressure (if available)
%       .CO2      [ppm]   Atmospheric CO2 concentration
%       .rb       [s m-1] Leaf boundary-layer resistance (thermal.rbc)
%       .GPP      [umol m-2 s-1 or model flux unit] Canopy GPP this step
%       .is_day   [0/1]   Daytime flag
%
%   Method:
%     LUNA (CLM5-style) acclimates leaf N / Vcmax using ~10-day means of
%     climate and productivity. This function does NOT run LUNA. It only
%     packages the current hour so commitPending can push it into
%     Climate10dBuffer after outer-step acceptance (reject => drop).
%
%     Day/night leaf-temperature slots are mutually exclusive so the buffer
%     can form separate daytime and nighttime leaf-T averages used by LUNA
%     temperature-response functions (e.g. Bernacchi / Kattge forms).
%
%   Reference:
%     vegetation.integrated.LUNA.Climate10dBuffer / runDailyLuna / computeLuna
%     CLM5 LUNA (Ali et al. style N optimization; CLM5 photosynthesis docs)
%
%   Editor: Jizixing Zhu; Date: 2026-07-16.
%   Change: 2026-07-17 Jizixing Zhu — expand file and inline comments for
%     formula meaning and day/night leaf-T slots.
    scope_step = model.built.workspace.scope_step;
    meteo = scope_step.meteo;

    % --- Daytime flag -------------------------------------------------------
    % Solar zenith angle tts [deg]: cos(tts) > 0  <=>  sun above horizon
    %   (tts < 90 deg). Used by SCOPE geometry when available.
    % Fallback: Rin > 0  (downwelling shortwave present).
    is_day = false;
    if isfield(scope_step, 'angles') ...
            && isfield(scope_step.angles, 'tts') ...
            && ~isempty(scope_step.angles.tts)
        is_day = cosd(scope_step.angles.tts) > 0;   % cosd: cos of degrees
    elseif isfield(meteo, 'Rin')
        is_day = meteo.Rin > 0;                     % [W m-2] incident SW
    end

    % Instantaneous PAR for LUNA light-environment statistics.
    par = readField(scope_step.rad, 'PAR', 0);
    if ~isscalar(par)
        par = 0;   % non-scalar / spectral vector not used here
    end

    fluxes = readField(scope_step, 'fluxes', struct());
    thermal = readField(scope_step, 'thermal', struct());

    % GPP this outer step: drives productivity memory in LUNA N investment.
    gpp = readField(fluxes, 'GPP', 0);

    % Leaf boundary resistance rb [s m-1]: aerodynamic limitation on CO2/H2O
    % exchange at the leaf surface (LUNA / photosynthesis environment).
    % Default 60 s/m if ebal did not write thermal.rbc.
    boundary_resistance = readField(thermal, 'rbc', 60);

    % Canopy-mean leaf temperature [degC] for acclimation.
    % Prefer ebal Tcave (big-leaf / two-leaf average); else air temperature.
    leaf_temperature = meteo.Ta;   % [degC] fallback
    if isfield(thermal, 'Tcave') && isfinite(thermal.Tcave)
        leaf_temperature = thermal.Tcave;   % [degC] canopy-averaged leaf T
    end

    % Mutual exclusion: one hour contributes to either day or night leaf-T mean.
    %   Tleaf_d  used in daytime temperature responses / day hours
    %   Tleaf_n  used in nighttime respiration temperature responses
    leaf_temperature_day = NaN;
    leaf_temperature_night = NaN;
    if is_day
        leaf_temperature_day = leaf_temperature;
    else
        leaf_temperature_night = leaf_temperature;
    end

    climate = struct( ...
        'Ta', meteo.Ta, ...                         % [degC] air T
        'Tleaf_d', leaf_temperature_day, ...        % [degC] day leaf T or NaN
        'Tleaf_n', leaf_temperature_night, ...      % [degC] night leaf T or NaN
        'PAR', par, ...                             % light for 10-d mean
        'RH', forcingRelativeHumidity(model, kt), ... % [-] vapor demand
        'ea', readField(meteo, 'ea', NaN), ...      % [hPa] vapor pressure
        'CO2', readField(meteo, 'Ca', 400), ...     % [ppm] Ca (default 400)
        'rb', boundary_resistance, ...              % [s m-1] leaf boundary layer
        'GPP', gpp, ...                             % productivity this hour
        'is_day', is_day);                          % [-] day/night mask
end

function relative_humidity = forcingRelativeHumidity(model, kt)
% Read forcing RH and clamp to [0,1] (accepts 0-1 or percent-like values only if already fraction).
% Editor: Jizixing Zhu; Date: 2026-07-16.
% Change: 2026-07-17 Jizixing Zhu — clarify clamp as physical RH bounds.
    relative_humidity = NaN;
    hourly = model.built.forcing.hourly;
    if isfield(hourly, 'RH') && numel(hourly.RH) >= kt
        % RH in site forcing is stored as a fraction in [0,1].
        relative_humidity = max(0, min(1, hourly.RH(kt)));
    end
end

function value = readField(container, name, fallback)
% Read a nonempty struct field or return a local fallback.
% Editor: Jizixing Zhu; Date: 2026-07-16.
    value = fallback;
    if isfield(container, name) && ~isempty(container.(name))
        value = container.(name);
    end
end
