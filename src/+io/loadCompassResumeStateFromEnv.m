function loaded = loadCompassResumeStateFromEnv()
% loadCompassResumeStateFromEnv Restore a COMPASS resume checkpoint.
%
% COMPASS.m is a script with many live state variables. Keeping the restore
% field list here avoids spreading checkpoint bookkeeping through the main
% physical loop.

loaded = false;
resumeFile = strtrim(getenv('STEMMUS_RESUME_FILE'));
if isempty(resumeFile)
    return
end

if ~exist(resumeFile, 'file')
    warning('COMPASS:ResumeFileNotFound', ...
        'STEMMUS_RESUME_FILE does not exist: %s', resumeFile);
    return
end

resumeLoaded = load(resumeFile, 'ResumeState');
if ~isfield(resumeLoaded, 'ResumeState')
    warning('COMPASS:InvalidResumeFile', ...
        'Resume file does not contain ResumeState: %s', resumeFile);
    return
end

R = resumeLoaded.ResumeState;
% Sim_output / flux_vars / Scope_output dropped (2026-05-15): the new
% diagnostics.Registry handles output buffers via NetCDF flush, not via
% in-memory structs restored on resume. The first KTs after resume re-fill
% the registry from scratch starting at R.KT.
restoreFields = { ...
    'SoilVariables', 'TransportCoefficient', 'BoundaryCondition', ...
    'SAVE', 'DSTOR', 'DSTOR0', 'NBChh', ...
    'AirVars', 'EnergyVars', 'Sa', ...
    'TimeStep', 'SUMTIME', 'NoTime', 'DayHour', 'StepFailed', ...
    'current_day', 'current_hour', 'current_date', 'current_year', ...
    'EVAP', 'CanopyInterception', 'Throughfall', 'Infiltration', ...
    'Infiltration_act', 'TopFlux_potential', 'SWE', ...
    'QET', 'QEB', 'QMT', 'QMB', 'RS', ...
    'PicardSummary', 'PicardHistory', ...
    'KT', 'TIME', 'Delt_old', 'Delt_t_outer', ...
    'wt_day_min', 'lai_day_prev_idx', ...
    'canopy', 'Soil', 'fluxes', 'thermal', 'rad', 'profiles', 'hydraulic', ...
    'Ta_history', ...
    'NitrogenVariables', 'NitrogenFluxes', 'NitrogenDiagnostics', 'NITTT'};

for i = 1:numel(restoreFields)
    name = restoreFields{i};
    if isfield(R, name)
        assignin('caller', name, R.(name));
    end
end

if isfield(R, 'Delt_t_outer')
    assignin('caller', 'Delt_t', R.Delt_t_outer);
end

% Note: KT and TIME are already restored as part of restoreFields above.
% If the caller's cfg.Time.KT / cfg.Time.TIME copies need refreshing, the
% main loop (COMPASS.m) does that explicitly — we no longer evalin('caller','cfg').

printResumeMessage(resumeFile, R);
loaded = true;
end

function printResumeMessage(resumeFile, R)
kt = NaN;
timeValue = NaN;
if isfield(R, 'KT'), kt = R.KT; end
if isfield(R, 'TIME'), timeValue = R.TIME; end

if isfield(R, 'failed_KT')
    fprintf('[RESUME] Loaded %s at KT=%d TIME=%.1f s; next failed KT was %d.\n', ...
        resumeFile, kt, timeValue, R.failed_KT);
else
    fprintf('[RESUME] Loaded %s at KT=%d TIME=%.1f s.\n', ...
        resumeFile, kt, timeValue);
end
end
