function [sub, out] = decideSubstep(sub, strat, dtOuter, dtMin, phase, info)
%DECIDESUBSTEP  Adaptive-substep controller: a pure start/accept/reject state machine.
%
%   [sub, out] = solver.decideSubstep(sub, strat, dtOuter, dtMin, phase, info)
%   advances the adaptive-substep runtime counters and returns the next step
%   decision. It is PURE: no I/O, no global state. The fixed strategy lives in
%   control.numerics.substep; the mutable counters live in runtime.substep. This
%   is the variable-step "shrink on non-convergence / grow on sustained success"
%   logic, factored out of runOneStep so it can be unit-tested without a solver.
%
%   Inputs:
%     sub     [inout] struct  runtime.substep counters:
%                       Delt_t [s] next suggested trial length,
%                       Delt_old [s] last accepted actual trial length,
%                       retry_count/success_count/hold_count [-], last_iters [-]
%     strat   [in]    struct  control.numerics.substep strategy:
%                       grow_factor [-], shrink_factor [-],
%                       successes_before_grow [-], hold_after_reject [-],
%                       grow_when_iters_below [-]
%     dtOuter [in]    [s]     configured outer step length (grow ceiling)
%     dtMin   [in]    [s]     minimum substep length (shrink floor)
%     phase   [in]    char    'start' | 'accept' | 'reject'
%     info    [in]    struct  phase payload: .dt [s] actual trial length
%                             (accept/reject), .iters [-] iteration count
%                             (accept)
%
%   Returns:
%     sub  [out] struct  updated counters
%     out  [out] struct  decision: .dt [s] length to attempt next (all phases),
%                        .give_up [logical] true when a reject reaches dtMin
%
%   State machine:
%     disabled: attempt dtOuter once; a rejection gives up without shrinking.
%     start : grow Delt_t by grow_factor (capped at dtOuter) only when a speed-up
%             is earned -- success_count >= successes_before_grow AND hold_count==0
%             AND last_iters < grow_when_iters_below -- then reset success_count;
%             otherwise keep Delt_t.
%     accept: store actual info.dt in Delt_old, then success_count++,
%             retry_count=0, decrement hold_count, and store last_iters;
%             then check earned-growth: if success_count >= successes_before_grow AND
%             hold_count==0 AND last_iters < grow_when_iters_below, grow Delt_t and
%             reset success_count.  Always returns the next suggested dt via out.dt.
%     reject: retry_count++, success_count=0, Delt_t = max(shrink_factor*dt,dtMin),
%             hold_count = hold_after_reject. The controller gives up only after
%             an attempted step at dtMin has itself failed.
%
%   Editor: Jizixing Zhu; Date: 2026-06-29.
%   Change: 2026-07-12 Jizixing Zhu -- honor strategy.enabled and require a real
%           solve attempt at dtMin before reporting give-up.
%   Change: 2026-07-14 Jizixing Zhu -- allow earned timestep growth between
%           accepted tiles while preserving accepted-dt and rollback semantics.
%   Change: 2026-07-14 Jizixing Zhu -- keep a clipped final tile in Delt_old
%           without overwriting the already earned next-step suggestion.

    out = struct('dt', sub.Delt_t, 'give_up', false);
    enabled = getFieldDef(strat, 'enabled', true);
    assert((islogical(enabled) || isnumeric(enabled)) && isscalar(enabled), ...
        'decideSubstep:badEnabled', 'strategy.enabled must be a logical scalar.');

    if ~logical(enabled)
        switch lower(char(phase))
            case 'start'
                sub.Delt_t = dtOuter;
                out.dt = dtOuter;
            case 'accept'
                sub.Delt_old = getFieldDef(info, 'dt', sub.Delt_t);
                sub.retry_count = 0;
                sub.success_count = 0;
                sub.hold_count = 0;
                sub.last_iters = getFieldDef(info, 'iters', 0);
            case 'reject'
                failedDt = getFieldDef(info, 'dt', sub.Delt_t);
                sub.Delt_t = failedDt;
                sub.retry_count = sub.retry_count + 1;
                sub.success_count = 0;
                sub.hold_count = 0;
                out.dt = failedDt;
                out.give_up = true;
            otherwise
                error('decideSubstep:badPhase', ...
                    'phase must be start|accept|reject, got "%s".', char(phase));
        end
        return;
    end

    switch lower(char(phase))
        case 'start'
            [sub, dt] = applyEarnedGrowth(sub, strat, dtOuter);
            out.dt = dt;

        case 'accept'
            sub.Delt_old = getFieldDef(info, 'dt', sub.Delt_t);
            sub.retry_count = 0;
            sub.success_count = sub.success_count + 1;
            if sub.hold_count > 0
                sub.hold_count = sub.hold_count - 1;
            end
            sub.last_iters = getFieldDef(info, 'iters', 0);
            [sub, dt] = applyEarnedGrowth(sub, strat, dtOuter);
            out.dt = dt;

        case 'reject'
            failedDt = getFieldDef(info, 'dt', sub.Delt_t);
            sub.retry_count = sub.retry_count + 1;
            sub.success_count = 0;
            sub.hold_count = strat.hold_after_reject;
            if failedDt <= dtMin
                out.give_up = true;
                dt = dtMin;
            else
                dt = max(failedDt * strat.shrink_factor, dtMin);
            end
            sub.Delt_t = dt;
            out.dt = dt;

        otherwise
            error('decideSubstep:badPhase', ...
                'phase must be start|accept|reject, got "%s".', char(phase));
    end
end

% ---------------------------------------------------------------------------
function [sub, dt] = applyEarnedGrowth(sub, strat, dtOuter)
% Shared growth gate used by both 'start' and 'accept'.  Grows Delt_t by
% grow_factor (capped at dtOuter) only when the earned-growth conditions are
% met; otherwise leaves Delt_t unchanged.  On growth, success_count is reset
% to zero so the next call cannot double-grow.
% Editor: Jizixing Zhu; Date: 2026-07-14.
    % The suggestion is always bounded by the active outer-step ceiling,
    % even when no new growth is earned (for example after a shorter final
    % calendar interval or a restored legacy checkpoint).
    dt = min(sub.Delt_t, dtOuter);
    sub.Delt_t = dt;
    lastIters = getFieldDef(sub, 'last_iters', 0);
    earnedGrow = sub.success_count >= strat.successes_before_grow ...
        && sub.hold_count == 0 ...
        && lastIters < strat.grow_when_iters_below;
    if earnedGrow
        dt = min(dt * strat.grow_factor, dtOuter);
        sub.success_count = 0;
    end
    sub.Delt_t = dt;
end

% ---------------------------------------------------------------------------
function v = getFieldDef(s, name, defaultValue)
% Read a struct field with a default (keeps the controller self-contained).
% Editor: Jizixing Zhu; Date: 2026-07-14.
    if isfield(s, name) && ~isempty(s.(name))
        v = s.(name);
    else
        v = defaultValue;
    end
end
