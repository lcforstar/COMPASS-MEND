function [model, pending, context] = stageSurfaceVegetation( ...
        model, kt, dt_s, fertilizer_applied)
%STAGESURFACEVEGETATION Stage the coupled canopy--surface candidate graph.
%
%   Inputs:
%     model              [inout] Accepted model shell with trial auxiliaries.
%     kt                 [-] Forcing index.
%     dt_s               [s] Actual outer-step duration.
%     fertilizer_applied [struct] Accept-only management event state.
%
%   Returns:
%     model   [inout] Shell with disposable SCOPE/snow/canopy auxiliaries.
%     pending [struct] Uncommitted domain candidates and interface requests.
%     context [struct] Accepted-snow/canopy input needed only for one bounded
%             snow-event restage.
%
%   Physical contract:
%     This is the ordered coupling of atmospheric forcing, canopy interception,
%     snow surface properties, radiation, ebal, canopy water, root uptake, and
%     daily vegetation.  It does not commit a physical state.  The order is
%     load-bearing because the same current-step snow candidate must be seen by
%     radiation, energy, and the soil boundary request.
%
%   Editor: Jizixing Zhu; Date: 2026-07-20.
%   Change: 2026-07-20 Jizixing Zhu -- seed the detached surface-water candidate
%     from accepted surface state before the soil active-set solve.
    [model, canopyInput] = vegetation.integrated.stageCanopyInputs( ...
        model, kt, dt_s);
    model = vegetation.integrated.stageCanopyEnergy(model, kt, dt_s);
    [model, postEbal] = vegetation.integrated.stagePostEbal( ...
        model, kt, dt_s, canopyInput.lai_state);
    postEbal.surface_water = model.built.Y.surface_water;
    dailyVegetation = vegetation.integrated.stageDailyVegetation(model, kt, dt_s);
    pending = step.buildPending( ...
        canopyInput, postEbal, dailyVegetation, fertilizer_applied);
    context = struct('canopy_input', canopyInput);
end
