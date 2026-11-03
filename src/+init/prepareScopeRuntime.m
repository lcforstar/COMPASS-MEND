function scopeRun = prepareScopeRuntime(built, scopeStatic)
%PREPARESCOPERUNTIME  Build pre-loop SCOPE runtime scaffolding.
%
%   scopeRun = init.prepareScopeRuntime(built, scopeStatic) separates the
%   once-per-run SCOPE runtime preparation from the per-time-step SCOPE input
%   slicing. It replaces the non-physical parts of the old
%   run_stemmus_initialize block [16]: xyt/telmax bookkeeping, output-structure
%   templates, and the pointer to the immutable optical library already loaded
%   by init.prepareRuntime.
%
%   Inputs:
%     built        [in] struct  initialized build container with forcing,
%                              control.time and parameters.scope.
%     scopeStatic  [in] struct  built.p.scope_static (the init static-optics cache), including
%                              spectral grids and immutable optical-library data.
%
%   Boundary: runtime.scope (this output) holds the runtime SCOPE workspace -- output templates,
%   time window, and a path reference to the static cache. It does NOT own the static optics cache;
%   that is built.p.scope_static (increment 4, 2026-07-05).
%
%   Returns:
%     scopeRun [struct] SCOPE runtime scaffold:
%       .xyt                 time/year vectors in SCOPE DOY convention
%       .time_window         start/end DOY and nearest indices
%       .telmax              number of SCOPE time elements for simulation=1
%       .vi                  legacy index vector placeholder (kept all ones)
%       .output_template     rad/thermal/fluxes/profiles templates
%       .step_entry          "init.prepareScopeStep(model, kt)"
%
%   Closure:
%     No state is advanced. Dynamic per-step objects (meteo, angles, leafbio,
%     canopy and Soil) are deliberately built by init.prepareScopeStep, because
%     they depend on kt, LAI, soil moisture, LUNA and solar geometry.
%
%   Editor: Jizixing Zhu; Date: 2026-06-28.
%   Change: 2026-07-05 Jizixing Zhu — scopeStatic now sourced from built.p.scope_static (was runtime.scope_static); doc the runtime.scope/built.p ownership boundary; update scope_static_path; increment 4.

    if nargin < 2 || isempty(scopeStatic)
        error('prepareScopeRuntime:missingScopeStatic', ...
            'scopeStatic (built.p.scope_static) is required.');
    end

    scope = built.parameters.scope;
    timeAxis = built.forcing.timeAxis;
    nTime = min(numel(timeAxis.clockTime), built.control.time.tS);
    clockTime = timeAxis.clockTime(1:nTime);

    scopeRun = struct();
    scopeRun.status = "prepared";
    scopeRun.options = scope.options;
    scopeRun.xyt = buildScopeXyt(clockTime);
    scopeRun.time_window = buildTimeWindow(scopeRun.xyt, built.control.time);
    scopeRun.telmax = determineTelmax(scope.options, scopeRun.xyt);
    scopeRun.vi = ones(1, 1);
    scopeRun.vmax = ones(1, 1);
    scopeRun.scope_static_path = "model.built.p.scope_static";
    scopeRun.step_entry = "init.prepareScopeStep(model, kt)";
    scopeRun.output_template = buildOutputTemplates(scopeStatic.spectral);
end

% -------------------------------------------------------------------------
function xyt = buildScopeXyt(clockTime)
% Build SCOPE time vectors: t = DOY-1 + clock-hour/24 and calendar year.
% Editor: Jizixing Zhu; Date: 2026-06-28.
    hr = hour(clockTime) + minute(clockTime) / 60 + second(clockTime) / 3600;
    xyt = struct();
    xyt.t = (day(clockTime, 'dayofyear') - 1) + hr / 24;
    xyt.year = year(clockTime);
end

% -------------------------------------------------------------------------
function tw = buildTimeWindow(xyt, timeCfg)
% Compute SCOPE start/end DOY and nearest indices from buildTime's clock time.
% Editor: Jizixing Zhu; Date: 2026-06-28.
    tw = struct();
    tw.startDOY = doyFraction(timeCfg.start_clock);
    tw.endDOY   = doyFraction(timeCfg.end_clock);
    [~, tw.I_tmin] = min(abs(xyt.t - tw.startDOY));
    [~, tw.I_tmax] = min(abs(xyt.t - tw.endDOY));
end

% -------------------------------------------------------------------------
function d = doyFraction(t)
% Convert a datetime to the SCOPE DOY convention: Jan 1 00:00 = 0.
% Editor: Jizixing Zhu; Date: 2026-06-28.
    d = (day(t, 'dayofyear') - 1) + ...
        hour(t) / 24 + minute(t) / 1440 + second(t) / 86400;
end

% -------------------------------------------------------------------------
function telmax = determineTelmax(options, xyt)
% Determine SCOPE loop length for the time-series mode used by COMPASS.
% Editor: Jizixing Zhu; Date: 2026-06-28.
    switch options.simulation
        case 0
            telmax = 1;
        case 1
            telmax = numel(xyt.t);
        case 2
            telmax = 1;
        otherwise
            error('prepareScopeRuntime:badSimulationMode', ...
                'scope.options.simulation must be 0, 1 or 2; got %d.', options.simulation);
    end
end

% -------------------------------------------------------------------------
function out = buildOutputTemplates(spectral)
% Allocate SCOPE output templates without running a radiative-transfer step.
% Editor: Jizixing Zhu; Date: 2026-06-28.
    [rad, thermal, fluxes] = io.initialize_output_structures(spectral);
    out = struct();
    out.rad = rad;
    out.thermal = thermal;
    out.fluxes = fluxes;
    out.profiles = struct();
end
