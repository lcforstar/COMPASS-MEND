function request = attachSurfaceState(request, built, regime, predictWetStart)
%ATTACHSURFACESTATE  Add accepted liquid-surface state to an atmospheric request.
%
%   The materialized atmospheric forcing remains immutable. This adapter adds
%   the accepted substep-start depression storage and one trial-local active
%   regime. The resulting request is still plain data and may be discarded or
%   rebuilt without changing model state.
%
%   Regimes:
%     supply_limited   all available liquid is prescribed; endpoint h <= 0
%     storage          endpoint positive head is implicit liquid depth
%     capacity_limited endpoint head equals capacity; excess is runoff
%     dry_evaporation_limited endpoint h is held at the accepted residual-
%                        water state; the solved total-water reaction is
%                        realized evaporation below demand
%
%   Frozen retained ponding is transferred transactionally into the bulk-snow
%   S/S_l/U candidate before this request is materialized. Only the remaining
%   liquid depth supplies pressure and infiltration here. A failed outer step
%   commits neither side of that internal transfer.
%
%   predictWetStart is true only for the first trial at an outer-step start.
%   On an abrupt positive liquid-supply jump, the first trial then uses the
%   wet feasible set (h=0 for zero capacity, implicit storage otherwise).
%   This is a predictor only: the solved reaction decides whether the same-dt
%   active set remains wet or returns to supply limitation.
%
%   Editor: Jizixing Zhu; Date: 2026-07-17.
%   Change: 2026-07-17 Jizixing Zhu -- start abrupt wetting events from the
%           head/storage feasible set without using a P>Ksat heuristic.
%   Change: 2026-07-17 Jizixing Zhu -- consume the staged liquid-to-snow
%           transfer instead of aborting long cold-season simulations.
%   Change: 2026-07-17 Jizixing Zhu -- materialize the pre-existing hcrit
%           dry-atmosphere constraint as an explicit active set instead of
%           overriding it with the typed atmospheric flux request.
%   Change: 2026-07-17 Jizixing Zhu -- enter the dry set from accepted residual-water
%           state rather than a standalone hcrit crossing, and preserve the
%           accepted head so activation cannot inject liquid water.
%   Change: 2026-07-20 Jizixing Zhu -- read the surface-owned substep candidate
%     instead of requiring an intermediate write to accepted Y.
    if nargin < 3, regime = ''; end
    if nargin < 4, predictWetStart = false; end
    if ~isstruct(request) || ~isfield(request, 'top') ...
            || ~strcmp(request.top.kind, 'atmospheric')
        return;
    end

    surfaceWater = surface.currentWaterState(built);
    depth0 = double(surfaceWater.DSTOR0);
    frozenTransfer = stagedFrozenTransfer(built);
    assert(frozenTransfer <= depth0 + 1e-10, ...
        'attachSurfaceState:frozenTransferAboveStorage', ...
        'Staged frozen transfer %.6g cm exceeds liquid storage %.6g cm.', ...
        frozenTransfer, depth0);
    depth0 = max(depth0 - frozenTransfer, 0);
    capacity = double(built.workspace.surface.DSTMAX);
    assert(isfinite(depth0) && depth0 >= 0 && isfinite(capacity) ...
            && capacity >= 0, ...
        'attachSurfaceState:badState', ...
        'Surface depth and capacity must be finite and nonnegative [cm].');
    tol = 1e-10;
    assert(depth0 <= capacity + tol, ...
        'attachSurfaceState:depthAboveCapacity', ...
        'Accepted surface depth %.6g cm exceeds capacity %.6g cm.', ...
        depth0, capacity);
    depth0 = min(depth0, capacity);
    dryHeadReference = double(built.workspace.surface.hcrit);
    assert(isscalar(dryHeadReference) && isfinite(dryHeadReference) ...
            && dryHeadReference < 0, ...
        'attachSurfaceState:badDryHead', ...
        'The atmospheric dry-head limit must be one finite negative scalar [cm].');
    acceptedDry = acceptedSurfaceIsResidual(built);
    dryHead = dryHeadReference;
    if acceptedDry
        dryHead = min(dryHeadReference, double(built.Y.h(end)));
    end

    initialization = 'explicit_active_set';
    if isempty(regime)
        if depth0 > tol && capacity > tol
            regime = 'storage';
            initialization = 'retained_storage';
        elseif logical(predictWetStart) && hasAbruptWetting(request, built)
            regime = wetRegime(capacity, tol);
            initialization = 'wet_feasibility_predictor';
        elseif shouldStartDryLimited(request, acceptedDry, tol)
            regime = 'dry_evaporation_limited';
            initialization = 'accepted_residual_state';
        else
            regime = 'supply_limited';
            initialization = 'accepted_unponded_state';
        end
    end
    allowed = {'supply_limited','storage','capacity_limited', ...
        'dry_evaporation_limited'};
    assert(any(strcmp(regime, allowed)), ...
        'attachSurfaceState:badRegime', ...
        'Unknown atmospheric surface regime "%s".', regime);
    if strcmp(regime, 'storage')
        assert(capacity > tol, 'attachSurfaceState:storageWithoutCapacity', ...
            'The storage regime requires positive surface capacity.');
    end

    request.top.surface = struct( ...
        'regime', regime, ...
        'storage_start_cm', depth0, ...
        'capacity_cm', capacity, ...
        'dry_head_limit_cm', dryHead, ...
        'initialization', initialization);
end

function tf = shouldStartDryLimited(request, acceptedDry, tol)
% Start from a previously accepted residual-water state with no liquid supply.
% Liquid supply always releases this predictor so rainfall can infiltrate.
% Editor: Jizixing Zhu; Date: 2026-07-17.
    supplyTol = max(1e-14, tol / max(request.dt_s, 1));
    tf = request.top.liquid_supply_cm_s <= supplyTol ...
        && request.top.evap_demand_cm_s > supplyTol ...
        && logical(acceptedDry);
end

function tf = acceptedSurfaceIsResidual(built)
% Test the top endpoint against the retention model's accepted residual state.
% Requiring both total and liquid water prevents freezing alone from being
% mistaken for atmospheric desiccation. The scale-aware tolerance only
% recognizes the kernel's residual clamp; it is not a tunable moisture band.
% Editor: Jizixing Zhu; Date: 2026-07-17.
    thetaR = double(built.parameters.soil.theta_r(end));
    thetaS = double(built.parameters.soil.theta_s(end));
    thetaL = double(built.Y.Theta_L(end, end));
    thetaU = thetaL;
    if isfield(built.Y, 'Theta_U') && ~isempty(built.Y.Theta_U)
        thetaU = double(built.Y.Theta_U(end, end));
    end
    thetaTol = max(1e-10, 1e-8 * max(thetaS - thetaR, 1e-3));
    tf = isfinite(thetaL) && isfinite(thetaU) ...
        && thetaL <= thetaR + thetaTol ...
        && thetaU <= thetaR + thetaTol;
end

function amount = stagedFrozenTransfer(built)
% Read the accept-only internal liquid-surface to snow transfer.
    amount = 0;
    if isfield(built.workspace, 'surface') ...
            && isfield(built.workspace.surface, 'frozen_transfer_cm') ...
            && ~isempty(built.workspace.surface.frozen_transfer_cm)
        amount = double(built.workspace.surface.frozen_transfer_cm);
    end
    assert(isscalar(amount) && isfinite(amount) && amount >= 0, ...
        'attachSurfaceState:badFrozenTransfer', ...
        'Frozen surface-water transfer must be finite and nonnegative [cm].');
end

function tf = hasAbruptWetting(request, built)
% Compare the current post-canopy/snow liquid supply with the last accepted one.
% Editor: Jizixing Zhu; Date: 2026-07-17.
    currentNet = max(request.top.liquid_supply_cm_s ...
        - request.top.evap_demand_cm_s, 0);
    previousSupply = 0;
    previousEvap = 0;
    if isfield(built.workspace, 'flux')
        if isfield(built.workspace.flux, 'Throughfall') ...
                && isscalar(built.workspace.flux.Throughfall) ...
                && isfinite(built.workspace.flux.Throughfall)
            previousSupply = built.workspace.flux.Throughfall;
        end
        if isfield(built.workspace.flux, 'EVAP') ...
                && isscalar(built.workspace.flux.EVAP) ...
                && isfinite(built.workspace.flux.EVAP)
            previousEvap = built.workspace.flux.EVAP;
        end
    end
    previousNet = max(previousSupply - previousEvap, 0);
    jumpTol = max(1e-10, 1e-6 * max(currentNet, previousNet)); % [cm s-1]
    tf = currentNet > previousNet + jumpTol;
end

function regime = wetRegime(capacity, tol)
% Select the wet feasible set without predicting the eventual runoff amount.
% Editor: Jizixing Zhu; Date: 2026-07-17.
    if capacity > tol
        regime = 'storage';
    else
        regime = 'capacity_limited';
    end
end

function temperature = surfaceTemperatureC(built)
% Read the staged physical surface temperature, falling back to topsoil T.
% Editor: Jizixing Zhu; Date: 2026-07-17.
    temperature = built.Y.T(end) - 273.15;
    if isfield(built.workspace, 'scope_step') ...
            && isfield(built.workspace.scope_step, 'thermal') ...
            && isfield(built.workspace.scope_step.thermal, 'Tss') ...
            && isscalar(built.workspace.scope_step.thermal.Tss) ...
            && isfinite(built.workspace.scope_step.thermal.Tss)
        temperature = built.workspace.scope_step.thermal.Tss;
    end
end
