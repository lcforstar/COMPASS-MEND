function surface = stageSurfaceWaterBoundary(built, kt) %#ok<INUSD>
%STAGESURFACEWATERBOUNDARY  Reset top-water boundary compatibility scratch.
%
%   Inputs:
%     built [in] Built model after the soil trial is warm-started.
%     kt    [-] One-based outer-step index.
%
%   Return:
%     surface [struct] Current surface boundary auxiliary with NBChh,
%             h_target, and hN staged for the water solve.
%
%   Atmospheric wet, ponded, runoff, and residual-dry regimes are resolved
%   transactionally from the typed request by attachSurfaceState and
%   evaluateSurfaceTrial. This compatibility stage must not preselect any of
%   those regimes from precipitation or pressure head. It only resets legacy
%   kernel scratch and retains an explicitly configured Dirichlet top code.
%   Editor: Jizixing Zhu; Date: 2026-07-16.
%   Change: 2026-07-17 Jizixing Zhu -- stop pinning a stored pond to its old depth;
%           the implicit surface-reservoir active set now owns ponding.
%   Change: 2026-07-17 Jizixing Zhu -- reset the staged liquid-to-snow
%           transfer at every new outer-step boundary stage.
%   Change: 2026-07-17 Jizixing Zhu -- remove the competing no-rain/hcrit switch;
%           the converged typed active set is the sole atmospheric policy.
    surface = built.workspace.surface;
    surface.frozen_transfer_cm = 0;
    surface.NBChh = defaultTopBoundaryCode(built);
    surface.h_target = 0;
    surface.hN = 0;
end

function code = defaultTopBoundaryCode(built)
% Return the configured top water-BC code, defaulting to flux mode.
    code = 2;
    if isfield(built, 'boundary') && isfield(built.boundary, 'NBCh_top') ...
            && built.boundary.NBCh_top == 1
        code = 1;
    end
end
