function solve_function = resolveAttemptFunction(model)
%RESOLVEATTEMPTFUNCTION  Select the production or injected soil trial solve.
%
%   Input:
%     model [in] Model runtime containing an optional debug attempt_solve.
%
%   Return:
%     solve_function [function_handle] Function with signature
%                    [model,outcome]=f(model,kt,dt_s,solve_control).
%
%   Dependency injection supports transaction tests without adding test
%   branches to the production scheduler and maps to future Julia dispatch.
%   Editor: Jizixing Zhu; Date: 2026-07-16.
    solve_function = @solver.attemptSolveSoilEquations;
    if isfield(model.runtime, 'debug') ...
            && isfield(model.runtime.debug, 'attempt_solve') ...
            && isa(model.runtime.debug.attempt_solve, 'function_handle')
        solve_function = model.runtime.debug.attempt_solve;
    end
end
