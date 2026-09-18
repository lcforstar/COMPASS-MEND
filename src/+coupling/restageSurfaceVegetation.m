function [model, pending] = restageSurfaceVegetation( ...
        model, kt, dt_s, context, fertilizer_applied, snow_surface, snow_override)
%RESTAGESURFACEVEGETATION Rebuild candidates once after a snow regime event.
%
%   Inputs:
%     model              [inout] Restored accepted-state shell.
%     kt                 [-] Forcing index.
%     dt_s               [s] Outer-step duration.
%     context            [struct] Context returned by stageSurfaceVegetation.
%     fertilizer_applied [struct] Accept-only management event state.
%     snow_surface       [struct] Completed snow candidate for the corrected
%                        radiation/energy pass.
%     snow_override      [struct] Optional conservative terminal-meltout reuse.
%
%   Returns:
%     model   [inout] Restaged disposable coupling auxiliaries.
%     pending [struct] Rebuilt uncommitted domain candidates/requests.
%
%   Physical contract:
%     A snow regime event changes optical, turbulent, and soil-boundary terms
%     together.  This bounded correction therefore rebuilds the complete
%     surface--vegetation graph from the same accepted state; it never swaps
%     only snow storage onto stale ebal fluxes and never advances snow twice.
%
%   Editor: Jizixing Zhu; Date: 2026-07-20.
%   Change: 2026-07-20 Jizixing Zhu -- preserve the detached surface-water candidate
%     during snow-event restaging so the transaction remains complete.
    assert(isstruct(context) && isfield(context, 'canopy_input') ...
            && isfield(context.canopy_input, 'lai_state'), ...
        'coupling:restageSurfaceVegetation:badContext', ...
        'The original canopy input and LAI candidate are required for restage.');
    model = vegetation.integrated.stageCanopyEnergy( ...
        model, kt, dt_s, snow_surface, false);
    [model, postEbal] = vegetation.integrated.stagePostEbal( ...
        model, kt, dt_s, context.canopy_input.lai_state, snow_override);
    % Keep the detached surface-water candidate in the same ownership path as
    % the ordinary staging pass; the snow restage must be a complete candidate
    % graph, not a partial replacement of snow/ebal fields.
    postEbal.surface_water = model.built.Y.surface_water;
    dailyVegetation = vegetation.integrated.stageDailyVegetation(model, kt, dt_s);
    pending = step.buildPending( ...
        context.canopy_input, postEbal, dailyVegetation, fertilizer_applied);
end
