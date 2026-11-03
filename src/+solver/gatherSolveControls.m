function solve_control = gatherSolveControls(built)
%GATHERSOLVECONTROLS  Build current canopy and soil iteration controls.
%
%   Input:
%     built [in] Built model containing immutable numerical configuration.
%
%   Return:
%     solve_control [struct] Energy, water, heat, and dry-air convergence
%                   controls for the current outer-step scheduler request.
%
%   Editor: Jizixing Zhu; Date: 2026-07-16.
%   Change: 2026-07-17 Jizixing Zhu -- pass every coupled conservation gate
%           explicitly into the trial solver.
%   Change: 2026-08-16 Jizixing Zhu -- forward the SITE-owned dry-air
%           inventory-relative conservation tolerance to every trial.
    numerics = built.control.numerics;
    solve_control = struct();
    solve_control.maxit = numerics.Max_Iter_EB;
    solve_control.maxEBer = numerics.Tol_EB;
    solve_control.Wc = 1;
    solve_control.NIT = numerics.NIT;
    solve_control.tol_theta = readNumeric(numerics, 'tol_theta', 1e-4);
    solve_control.tol_head_saturated_cm = readNumeric( ...
        numerics, 'tol_head_saturated_cm', 0.1);
    solve_control.tol_water_mass_mm = readNumeric( ...
        numerics, 'tol_water_mass_mm', 1e-4);
    solve_control.tol_surface_water_cm = readNumeric( ...
        numerics, 'tol_surface_water_cm', ...
        max(1e-10, solve_control.tol_water_mass_mm / 10));
    solve_control.tol_energy = numerics.tol_energy;
    solve_control.tol_energy_balance_J_cm2 = readNumeric( ...
        numerics, 'tol_energy_balance_J_cm2', 1e-2);
    solve_control.tol_air_relative = readNumeric( ...
        numerics, 'tol_air_relative', 1e-3);
    solve_control.tol_air_mass_g_cm2 = readNumeric( ...
        numerics, 'tol_air_mass_g_cm2', 1e-8);
    solve_control.tol_air_mass_relative = readNumeric( ...
        numerics, 'tol_air_mass_relative', 1e-4);
end

function value = readNumeric(container, name, fallback)
% Read a nonempty numerical control or return its documented fallback.
    value = fallback;
    if isfield(container, name) && ~isempty(container.(name))
        value = container.(name);
    end
end
