function built = extrapolateSubstepGuess(built, kt, dt_s, prev_dt_s)
%EXTRAPOLATESUBSTEPGUESS  Variable-dt Picard initial guess for a substep (Layer A / A3).
%
%   built = soil.water.extrapolateSubstepGuess(built, kt, dt_s, prev_dt_s) writes THIS substep's
%   Picard initial guess into built.workspace.step.trial.{h,T,p_gas} by variable-dt extrapolation
%   from the accepted state Y and the previous-substep-start anchor. Ports COMPASS.m:1586-1608.
%   It is the Picard initial GUESS only -- NOT accept, NOT carry, NOT a previous-update. It writes
%   trial ONLY (never Y) and calls no refresh/solve.
%
%   Inputs:
%     built     [inout] container; reads Y.{h,T,p_gas}, workspace.step.previous.{h,T,p_gas},
%                       workspace.step.trial.{h,T,p_gas} (the warm-start), control.numerics.lExtrap,
%                       control.switches.{lThmrlefc,lSoilairefc}
%     kt        [in]    [-] outer-step index (gate: extrapolate only when kt>1)
%     dt_s      [in]    [s] this substep's dt
%     prev_dt_s [in]    [s] the previous substep dt (alpha_sub = dt_s/prev_dt_s). EXPLICIT INPUT --
%                       the future loop supplies Delt_old (first substep) or sub_last_accepted_dt;
%                       A3 does NOT read any loop state. Must be > 0 when the gate is on (else error).
%
%   Returns:
%     built  with trial.{h,T,p_gas} overwritten by the extrapolated guess (when kt>1 && lExtrap==1),
%            else unchanged (the warm-start trial = Y stays).
%
%   Method (COMPASS.m:1596-1606, elementwise / regional-safe):
%     trial.h = Y.h + alpha*(Y.h - previous.h) where trial.h<-0.01 & Y.h<-0.01, else trial.h = Y.h;
%     trial.T = Y.T + alpha*(Y.T - previous.T) (if lThmrlefc==1);
%     trial.p_gas = Y.p_gas + alpha*(Y.p_gas - previous.p_gas) (if lSoilairefc==1).
%
%   Reference: docs/superpowers/specs/2026-07-07-stage4-layerA-substep-prepare-design.md sec 4 (A3).
%   Oracle: tests/test_extrapolateSubstepGuess.m (golden reference vs COMPASS.m:1596-1606).
%   Editor: Jizixing Zhu; Date: 2026-07-07.
%   Change: 2026-07-14 Jizixing Zhu -- extrapolate U when heat equation is enabled.
    if ~(kt > 1 && built.control.numerics.lExtrap == 1)
        return;   % gate off: keep the warm-start trial = Y (prev_dt_s unused, not validated)
    end
    if prev_dt_s <= 0
        error('soil:extrapolateSubstepGuess:badPreviousDt', ...
            'prev_dt_s must be > 0 when extrapolating (alpha = dt_s/prev_dt_s); got %g.', prev_dt_s);
    end
    alpha = dt_s / prev_dt_s;
    sw = built.control.switches;
    Yv = built.Y;
    pv = built.workspace.step.previous;
    tr = built.workspace.step.trial;

    % h: extrapolate where the incoming trial.h (=hh) and Y.h (=h) are both < -0.01; else fall back.
    coldMask = (tr.h < -0.01) & (Yv.h < -0.01);
    newh = Yv.h;
    newh(coldMask) = Yv.h(coldMask) + alpha .* (Yv.h(coldMask) - pv.h(coldMask));
    tr.h = newh;

    if isfield(sw, 'lThmrlefc') && sw.lThmrlefc == 1
        tr.T = Yv.T + alpha .* (Yv.T - pv.T);
    end
    if isfield(sw, 'lSoilairefc') && sw.lSoilairefc == 1
        tr.p_gas = Yv.p_gas + alpha .* (Yv.p_gas - pv.p_gas);
    end
    if isfield(sw, 'lThmrlefc') && sw.lThmrlefc == 1 ...
            && isfield(Yv, 'U') && isfield(pv, 'U')
        tr.U = Yv.U + alpha .* (Yv.U - pv.U);
    end

    built.workspace.step.trial = tr;
end
