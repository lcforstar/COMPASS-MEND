function [model, advance, storage_before_mm] = advance( ...
        model, kt, dt_outer_s, dt_min_s, strategy, solve_control, ...
        surface_water_request, t_start_s, solve_function)
%ADVANCE Advance the coupled unsaturated-soil candidate for one outer interval.
%
%   Inputs:
%     model                 [inout] Accepted shell plus disposable trial state.
%     kt                    [-] Atmospheric forcing index.
%     dt_outer_s, dt_min_s  [s] Outer interval and adaptive floor.
%     strategy              [struct] Adaptive substep controller settings.
%     solve_control         [struct] Coupled h/U/P_g controls and root request.
%     surface_water_request [struct] Typed surface-to-soil water request.
%     t_start_s             [s] Accepted physical start time.
%     solve_function        [function_handle] One coupled soil attempt.
%
%   Returns:
%     model             [inout] Best soil trial only; not yet accepted state.
%     advance           [struct] Adaptive solve and interface-closure result.
%     storage_before_mm [mm] Pre-solve soil water storage for the ledger.
%
%   Physical contract:
%     Soil owns the h/U/P_g solve, internal phase partition, adaptive retries,
%     and the actual top/bottom flux result.  The surface request is immutable
%     during each trial.  The facade preserves the existing dynamic heat
%     closure and never reruns ebal inside the soil Picard iteration.
%
%   Editor: Jizixing Zhu; Date: 2026-07-20.
    [model, advance, storage_before_mm] = ...
        solver.advanceSoilWithInterfaceTemperatureClosure( ...
        model, kt, dt_outer_s, dt_min_s, strategy, solve_control, ...
        surface_water_request, t_start_s, solve_function);
end
