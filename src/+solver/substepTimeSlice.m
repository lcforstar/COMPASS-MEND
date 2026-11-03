function dt = substepTimeSlice(sub_dt_next, sub_elapsed, Delt_t_outer)
%SUBSTEPTIMESLICE  Clip the next trial substep length to the remaining outer step (Layer A / A1).
%
%   dt = solver.substepTimeSlice(sub_dt_next, sub_elapsed, Delt_t_outer) returns the length of the
%   current trial substep as min(sub_dt_next, Delt_t_outer - sub_elapsed). Ports COMPASS.m:1566-1569
%   (the body's slice line). It is a PURE function: it computes only the current trial dt. It does
%   NOT advance sub_elapsed, does NOT touch sub_dt_next, does NOT call solver.decideSubstep, does not
%   mutate runtime.substep, opens no loop, and wires nothing into runOneStep/attemptSolve. The
%   multi-substep accumulator loop is a later loop-wiring increment.
%
%   Inputs:
%     sub_dt_next   [in]  [s]  next trial length proposed by solver.decideSubstep as its out.dt field,
%                              capped at the full outer step but NOT at the remaining hour. Must be > 0.
%     sub_elapsed   [in]  [s]  outer-step time already consumed by accepted substeps. Must be >= 0.
%     Delt_t_outer  [in]  [s]  configured outer step length. Must be > 0.
%
%   Returns:
%     dt  [s]  min(sub_dt_next, Delt_t_outer - sub_elapsed); 0 when the outer step is complete
%              (remaining within +/- tol). A zero slice is the future loop's stop signal.
%
%   Guards (fail-fast): sub_dt_next <= 0 -> badDt; Delt_t_outer <= 0 -> badOuterDt;
%     sub_elapsed < 0 -> badElapsed; sub_elapsed beyond Delt_t_outer by more than tol -> overElapsed.
%     A remaining within [-tol, +tol] collapses to dt = 0 (tol = 1e-9 s, matching the legacy loop
%     guard COMPASS.m:1566, so a tiny-positive remaining is not emitted as a micro-substep and a
%     float-noise overshoot is not misread as an over-run).
%
%   Reference: COMPASS.m:1566-1569 (sub_remaining / Delt_t = min(sub_dt_next, sub_remaining)).
%   Editor: Jizixing Zhu; Date: 2026-07-07.
    if ~(isscalar(sub_dt_next) && isfinite(sub_dt_next) && sub_dt_next > 0)
        error('solver:substepTimeSlice:badDt', ...
            'sub_dt_next must be a finite positive scalar; got %g.', sub_dt_next);
    end
    if ~(isscalar(Delt_t_outer) && isfinite(Delt_t_outer) && Delt_t_outer > 0)
        error('solver:substepTimeSlice:badOuterDt', ...
            'Delt_t_outer must be a finite positive scalar; got %g.', Delt_t_outer);
    end
    if ~(isscalar(sub_elapsed) && isfinite(sub_elapsed) && sub_elapsed >= 0)
        error('solver:substepTimeSlice:badElapsed', ...
            'sub_elapsed must be a finite non-negative scalar; got %g.', sub_elapsed);
    end

    tol = 1e-9;                                    % [s] matches the legacy loop guard COMPASS.m:1566
    remaining = Delt_t_outer - sub_elapsed;        % [s] outer step not yet consumed
    if remaining < -tol
        error('solver:substepTimeSlice:overElapsed', ...
            'sub_elapsed (%g) exceeds Delt_t_outer (%g) beyond tolerance.', sub_elapsed, Delt_t_outer);
    end
    if remaining <= tol
        dt = 0;                                    % outer step complete: the future loop must not solve
        return;
    end
    dt = min(sub_dt_next, remaining);              % [s] clip the trial length to the remaining hour
end
