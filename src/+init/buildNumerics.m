function num = buildNumerics(rawConfig)
%BUILDNUMERICS  Normalize [numerics] into legacy solver settings.
%
%   num = init.buildNumerics(rawConfig) maps the friendly site.toml numerics
%   names onto the legacy solver fields (NIT, lUpW, kernels, ...). Pure config
%   translation -- no physics.
%
%   Inputs:
%     rawConfig  struct  raw.config (uses .numerics)
%
%   Returns num  struct  legacy numeric solver settings:
%     NIT/Max_Iter_EB/Tol_EB/tol_theta/tol_head_saturated_cm/
%     tol_water_mass_mm/tol_surface_water_cm/tol_energy/
%     tol_energy_balance_J_cm2/tol_air_relative/tol_air_mass_g_cm2/
%     tol_air_mass_relative,
%     HydraulicKernel/MatricKernel,
%     lExtrap/lUpW/lSat/lArtD, checkpointInterval/diagnosticInterval,
%     substep.{enabled,grow_factor,shrink_factor,successes_before_grow,hold_after_reject,grow_when_iters_below}.
%
%   Reference: legacy ModelSettings (field names).
% Editor: Jizixing Zhu; Date: 2026-06-27.
% Change: 2026-07-11 Jizixing Zhu -- allow water-only sites to omit inactive energy and transport numerics.
% Change: 2026-07-12 Jizixing Zhu -- make HYDRUS increment and actual-mass
%           acceptance controls fixed solver behavior.
% Change: 2026-07-13 Jizixing Zhu -- add the relative dry-air pressure
%           increment tolerance used by the unified soil Picard solve.
% Change: 2026-07-15 Jizixing Zhu -- remove the completed validation switch;
%           h/U/P_g is now the sole production soil coordinate path.
% Change: 2026-07-17 Jizixing Zhu -- normalize explicit conservation gates
%           for surface water, internal energy, and dry-air mass.
% Change: 2026-07-18 Jizixing Zhu -- remove the discarded A1 surface/soil reaction
%           correction controls after the bounded experiment showed no gain.
% Change: 2026-06-28 Jizixing Zhu — add control.numerics.substep adaptive-step strategy (from COMPASS.m [6] debug_substep_*).
% Change: 2026-06-28 Jizixing Zhu — rename substep fields to self-documenting names (grow_factor/shrink_factor/successes_before_grow/grow_when_iters_below).
% Change: 2026-07-19 Jizixing Zhu -- normalize nested Dirichlet soil-interface
%           temperature closure controls.
% Change: 2026-07-20 Jizixing Zhu -- normalize atmospheric heat-closure temporal
%           retry and bounded-residual fallback controls.
% Change: 2026-08-16 Jizixing Zhu -- expose the dry-air column inventory-relative
%           conservation tolerance as an explicit SITE numerical control.

    nm  = rawConfig.numerics;
    num = struct();
    num.NIT             = nm.NIT;
    num.Max_Iter_EB     = readOptionalNum(nm, 'max_iter_eb', 1);
    num.Tol_EB          = readOptionalNum(nm, 'tol_eb', 1);
    num.tol_theta = readOptionalNum(nm, 'tol_theta', 1e-4);                 % [-]
    num.tol_head_saturated_cm = readOptionalNum( ...
        nm, 'tol_head_saturated_cm', 0.1);                                  % [cm]
    num.tol_water_mass_mm = readOptionalNum( ...
        nm, 'tol_water_mass_mm', 1e-4);                                     % [mm]
    num.tol_surface_water_cm = readOptionalNum( ...
        nm, 'tol_surface_water_cm', max(1e-10, num.tol_water_mass_mm / 10)); % [cm]
    num.tol_energy      = readOptionalNum(nm, 'tol_energy', 1);   % [degC-equivalent] inactive in water-only sites
    num.tol_energy_balance_J_cm2 = readOptionalNum( ...
        nm, 'tol_energy_balance_J_cm2', 1e-2);                    % [J cm-2]
    num.tol_air_relative = readOptionalNum(nm, 'tol_air_relative', 1e-3); % [-]
    num.tol_air_mass_g_cm2 = readOptionalNum( ...
        nm, 'tol_air_mass_g_cm2', 1e-8);                          % [g cm-2]
    num.tol_air_mass_relative = readOptionalNum( ...
        nm, 'tol_air_mass_relative', 1e-4);                        % [-]
    num.HydraulicKernel = char(nm.hydraulic_kernel);
    num.MatricKernel    = char(nm.matric_kernel);
    num.lExtrap         = double(readOptionalBool(nm, 'extrapolate', true));
    num.lUpW            = double(readOptionalBool(nm, 'upwind', false));
    num.lSat            = double(readOptionalBool(nm, 'saturation_enforce', true));
    num.lArtD           = double(readOptionalBool(nm, 'artificial_disp', false));
    num.checkpointInterval = readOptionalNum(nm, 'checkpoint_interval', 0);
    num.diagnosticInterval = readOptionalNum(nm, 'diagnostic_interval', 0);
    num.soil_interface_closure = struct( ...
        'enabled', readOptionalBool( ...
            nm, 'soil_interface_closure_enable', false), ...
        'tolerance_W_m2', readOptionalNum( ...
            nm, 'soil_interface_closure_tol_W_m2', 5), ...
        'max_iterations', max(1, round(readOptionalNum( ...
            nm, 'soil_interface_closure_max_iters', 6))), ...
        'substep_retry_enabled', readOptionalBool( ...
            nm, 'soil_interface_closure_substep_retry_enable', true), ...
        'substep_retry_factor', readOptionalNum( ...
            nm, 'soil_interface_closure_substep_retry_factor', 0.5), ...
        'substep_retry_max_levels', max(0, round(readOptionalNum( ...
            nm, 'soil_interface_closure_substep_retry_max_levels', 3))), ...
        'bounded_acceptance_enabled', readOptionalBool( ...
            nm, 'soil_interface_closure_bounded_acceptance_enable', true), ...
        'bounded_tolerance_W_m2', readOptionalNum( ...
            nm, 'soil_interface_closure_bounded_tol_W_m2', 15), ...
        'bounded_relative_tolerance', readOptionalNum( ...
            nm, 'soil_interface_closure_bounded_relative_tol', 0.20));
    % adaptive substep controller strategy (legacy COMPASS.m [6] debug_substep_*).
    % Static rules for shrinking/growing the outer step under non-convergence; the
    % runtime state (current dt / previous dt / retry count) lives in runOneStep,
    % and the floor min-dt is control.time.Delt_t_min (site.toml [time] dt_min).
    sub = struct();
    sub.enabled               = readOptionalBool(nm, 'substep_enabled', true);
    sub.grow_factor           = readOptionalNum(nm, 'substep_grow_factor', 2.0);            % [-] multiply dt on sustained success
    sub.shrink_factor         = readOptionalNum(nm, 'substep_shrink_factor', 0.5);          % [-] multiply dt on rejection
    sub.successes_before_grow = readOptionalNum(nm, 'substep_successes_before_grow', 4);    % [-] consecutive successes before growing dt
    sub.hold_after_reject     = readOptionalNum(nm, 'substep_hold_after_reject', 6);        % [-] steps to hold small dt after a reject
    sub.grow_when_iters_below = max(8, ceil(0.2 * num.NIT));                                % [-] grow dt only when Picard iters below this
    num.substep = sub;
end

function v = readOptionalBool(s, name, defaultValue)
% Read an optional boolean numerics key with a default.
% Editor: Jizixing Zhu; Date: 2026-06-28.
    if isfield(s, name)
        v = logical(s.(name));
    else
        v = logical(defaultValue);
    end
end

% ---------------------------------------------------------------------------
function v = readOptionalNum(s, name, defaultValue)
% Read an optional numeric numerics key with a finite-value default.
% Editor: Jizixing Zhu; Date: 2026-06-28.
    if isfield(s, name) && ~isempty(s.(name)) && isfinite(s.(name))
        v = double(s.(name));
    else
        v = double(defaultValue);
    end
end
