function built = refreshWaterAux(built, dt_s, L_f, Tmin)
%REFRESHWATERAUX  B1-complete water-aux refresh: freeze-point head + water-content aux.
%
%   built = soil.refreshWaterAux(built, dt_s, L_f, Tmin) is the B1 refresh entry of the
%   Stage-4 pre-solve refresh layer. Internal order is PINNED to legacy COMPASS.m:1677-1684:
%     1. (if lFreezing==1) recompute the freeze-point TT_CRIT / trial.h_frez via
%        soil.water.computeFreezingPoint (the HT_frez head) and write them into the container;
%     2. refresh the water-content auxiliary cache via soil.water.updateWaterContentAuxiliary.
%   The freeze-point step is the FIRST step, not an independent bypass. This is a cross-domain
%   soil pre-solve orchestration adapter, so it lives at top-level +soil (not +water).
%
%   Inputs:
%     built [inout] canonical model container (built.Y, built.workspace, built.parameters,
%                   built.control)
%     dt_s  [in]    [s] outer update interval (0 during init/refresh-only use)
%     L_f   [in]    [J/kg] fusion latent heat; empty -> resolved as:
%                   built.workspace.step.L_f (set by runOneStep:298 = computeFreezeLatent)
%                   if present, else lambda_fus when lFreezing==1, else 0. lFreezing==1
%                   REQUIRES L_f~=0 (the freeze-point formula divides by L_f) -> errors otherwise.
%     Tmin  [in]    [degC] air temperature Ta at the step (freeze-point clamp floor); default 0
%
%   Returns:
%     built with derived.TT_CRIT, derived.DhFrezDh/DhFrezDT, and trial.h_frez
%     refreshed (when lFreezing==1) and the
%     water-aux set (Theta_LL/UU/II, Se, KL_h/KfL_h, DTheta_LLh/UUh, Ratio_ice,
%     Gamma_hh, COR/CORh) refreshed.
%
%   RULE: new stage-4 code MUST NOT call legacy HT_frez / UpdateSoilWaterContent (read-only
%   oracles). The new refresh path is computeFreezingPoint + updateWaterContentAuxiliary.
%
%   Reference: docs/superpowers/specs/2026-07-07-stage4-presolve-refresh-roadmap-design.md sec 5
%   (B1-tail). Oracle: tests/test_refreshWaterAux.m, tests/test_computeFreezingPoint_oracle.m.
%   Editor: Jizixing Zhu; Date: 2026-07-07.
%   Change: 2026-07-11 Jizixing Zhu -- retain the active-branch freezing-head
%           derivatives needed by the conservative phase-storage Jacobian.
%   Change: 2026-07-14 Jizixing Zhu -- remove the retired KfL_T field from the
%           documented strict-SFCC refresh contract.
    if nargin < 2 || isempty(dt_s), dt_s = 0; end
    if nargin < 3, L_f = []; end
    if nargin < 4 || isempty(Tmin), Tmin = 0; end

    sw = built.control.switches;
    lFreezing = isfield(sw, 'lFreezing') && sw.lFreezing == 1;

    % L_f resolution -- match the runOneStep step-state source (runOneStep.m:298 sets
    % built.workspace.step.L_f = computeFreezeLatent(switches)). Order: explicit arg >
    % step-state L_f > freeze-switch default. Do NOT bypass the step-state override.
    if isempty(L_f)
        if isfield(built.workspace, 'step') && isfield(built.workspace.step, 'L_f') ...
                && ~isempty(built.workspace.step.L_f)
            L_f = built.workspace.step.L_f;
        elseif lFreezing
            L_f = equations.getPhysicalConstant('lambda_fus');
        else
            L_f = 0;
        end
    end
    if lFreezing && L_f == 0
        error('soil:refreshWaterAux:zeroLatentHeatWithFreezing', ...
            'lFreezing=1 requires L_f~=0 (the freeze-point formula divides by L_f).');
    end

    if lFreezing
        T0 = equations.getPhysicalConstant('T_freeze');
        g  = equations.getPhysicalConstant('g');
        hd = built.parameters.soil.constants.hd;
        nCell = 1;
        if isfield(built.domain, 'nCell') && ~isempty(built.domain.nCell)
            nCell = built.domain.nCell;
        end
        % Regional-safe: loop soil columns; trial/derived may be [NN x nCell].
        for c = 1:nCell
            hh = soil.water.getNodeColumn(built.workspace.step.trial.h, c);
            TT = soil.water.getNodeColumn(built.workspace.step.trial.T, c);
            [TT_CRIT, h_frez, freezePartials] = ...
                soil.water.computeFreezingPoint(hh, TT, T0, g, L_f, hd, Tmin);
            built.workspace.derived.TT_CRIT = soil.water.setNodeColumn( ...
                built.workspace.derived.TT_CRIT, c, TT_CRIT);
            built.workspace.derived.DhFrezDh = soil.water.setNodeColumn( ...
                built.workspace.derived.DhFrezDh, c, freezePartials.dh_frez_dh);
            built.workspace.derived.DhFrezDT = soil.water.setNodeColumn( ...
                built.workspace.derived.DhFrezDT, c, freezePartials.dh_frez_dT);
            built.workspace.step.trial.h_frez = soil.water.setNodeColumn( ...
                built.workspace.step.trial.h_frez, c, h_frez);
        end
    end

    built = soil.water.updateWaterContentAuxiliary(built, dt_s, L_f);
end
