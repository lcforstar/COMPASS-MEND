function [built, CHKWATER, waterFlux, boundaryApplied, linearSystem] = solveWaterTrial( ...
        built, dt, kt, waterBoundaryFlux, thetaEqN, waterRootSink, ...
        boundaryRequest, waterCoefficients, waterMatrix, assembleOnly)  %#ok<INUSD>
%SOLVEWATERTRIAL  One trial-substep Richards (water) solve on the built's trial state (Layer 4 / I3).
%
%   [built, CHKWATER] = soil.water.solveWaterTrial(
%   built, dt, kt, waterBoundaryFlux, thetaEqN) assembles the legacy
%   SoilVariables/TransportCoefficient/cfg/BoundaryCondition frame from built (trial scratch + the
%   refreshed aux in workspace.derived), calls the EXISTING soil.water.solveSoilMoistureBalance
%   kernel (READ-ONLY reuse), and scatters the updated matric head and the boundary-flux snapshot
%   back to the trial scratch. It is a marshaling adapter: it rewrites no physics and enters no other
%   solve. The kernel updates ONLY hh (Theta_LL is recomputed by the next Picard refresh, not here).
%
%   Inputs:
%     built  [inout] container. Reads workspace.step.trial.{h [cm],T [K],h_frez,p_gas} (the Picard
%                    iterate), Y.{h,T} (the current accepted substep-start time layer), Y.Theta_L
%                    (substep-start liquid theta), workspace.derived.* (refreshed aux: KL_h/KfL_h/
%                    DTheta_LLh/RHOV/COR/CORh/Theta_LL/Theta_UU/Theta_V/D_V/D_Vg/...), workspace.flux.
%                    {Srt,V_A,Throughfall,EVAP}, workspace.solver.SAVE, workspace.surface.{DSTMAX,
%                    NBChh,hcrit,h_target,hN}, boundary.{NBCh_top,NBCh_bottom,h_top,h_bottom},
%                    Y.surface_water.DSTOR0, parameters.soil.ks, control.switches/numerics/domain.
%     dt     [in]    [s] trial substep length (Delt_t).
%     kt     [in]    [-] outer-step index. NOT forwarded into the current-step legacy frame (which is
%                    driven at localKt=1, see the body); retained for the caller contract/diagnostics.
%     thetaEqN [in]   [-, NL x 2] accepted-substep total water-equivalent
%                    endpoint content; fixed through every Picard iteration.
%
%   Returns:
%     built     with workspace.step.trial.h <- SoilVariables.hh and workspace.solver.SAVE <- SAVE.
%               Writes ONLY trial scratch; NOT Y, NOT pending, NOT the runtime clock.
%     CHKWATER  [cm, per node] the Picard water residual; the caller converges on max(CHKWATER)
%               (legacy Error_Water = max(CHKWATER), COMPASS.m:1836).
%     waterFlux struct of the kernel's water-boundary DIRECT outputs for this solve (caller-facing
%               scalars): QMT/QMB [cm/s] top/bottom matric flux, AVAIL0 [cm/s] available surface
%               supply, Throughfall [cm/s], Infiltration/Infiltration_act [cm/s] potential/actual
%               infiltration, TopFlux_potential [cm/s] top BC potential, runoffHort [cm/s]. For a
%               transactional atmospheric request runoffHort is the realized capacity overflow;
%               direct legacy calls retain the historical Ks pre-clip. RETURN ONLY -- not scattered
%               to built (the WBC-2 accumulator consumes it).
%     boundaryApplied struct describing the active mathematical top/bottom BC
%               kind used by the kernel for this trial; empty on the legacy
%               no-request path. For atmospheric requests, the top kind records
%               the current legacy switch as fixed_head or fixed_flux.
%
%   Scope: the water-STATE outputs hh + SAVE are SCATTERED to built (trial.h, solver.SAVE); the
%   water-BOUNDARY fluxes the kernel returns are packed into waterFlux and RETURNED (WBC-1) for the
%   caller-level water balance -- they are NOT scattered to built. This adapter does NOT update
%   Theta_LL (the kernel writes only hh; Theta_LL is recomputed by the next Picard refreshWaterAux).
%
%   Reference: legacy Picard water solve COMPASS.m:1731-1735; kernel soil.water.solveSoilMoistureBalance.
%   Editor: Jizixing Zhu; Date: 2026-07-07.
%   Change: 2026-07-08 Jizixing Zhu -- WBC-1: return the water-boundary direct outputs as waterFlux (return-only; scatter unchanged).
%   Change: 2026-07-09 Jizixing Zhu -- fix kt>=2 crash (soak-discovered): drive the legacy kernel at localKt=1 (the current-step frame is length-1), not the global kt; read the KT-indexed waterFlux outputs back at localKt; add a single-column guard.
%   Change: 2026-07-09 Jizixing Zhu -- take the this-step water-boundary forcing (waterBoundaryFlux.Throughfall/EVAP) as an optional 4th arg (default built.workspace.flux); the solve top BC now reads the current-step pending forcing, not the last-committed built.workspace.flux.
%   Change: 2026-07-10 Jizixing Zhu -- forward accepted storage to the water time term.
%   Change: 2026-07-10 Jizixing Zhu -- map legacy h/T to the current accepted substep-start Y.h/Y.T; workspace.step.previous is hOLD/TOLD for A3 extrapolation, not the Richards time-old layer.
%   Change: 2026-07-10 Jizixing Zhu -- accept optional 6th arg waterRootSink (this-step root sink) -> SV.Srt; default built.workspace.flux.Srt (backward-compat).
%   Change: 2026-07-11 Jizixing Zhu -- fail fast when a typed series BC reaches the legacy adapter without a per-trial materialized request.
%   Change: 2026-07-11 Jizixing Zhu -- W1b: accept an optional immutable BoundaryRequest and translate canonical top/bottom policies into the existing Richards-kernel boundary frame; the no-request path remains backward-compatible.
%   Change: 2026-07-11 Jizixing Zhu -- W1b: return the applied mathematical BC kind separately from the request so atmospheric fixed-head/fixed-flux switching is not mislabeled in audit results.
%   Change: 2026-07-11 Jizixing Zhu -- route the new transactional solve through the
%           conservative phase-storage coefficient assembler; legacy entry
%           points retain calculateMatricCoefficients unchanged.
%   Change: 2026-07-12 Jizixing Zhu -- permit positive pressure heads on the new
%           variably saturated path; legacy direct callers retain their clamp.
%   Change: 2026-07-12 Jizixing Zhu -- read the canonical liquid-temperature
%           partial after removal of the dynamic legacy DTheta_LLT cache.
%   Change: 2026-07-12 Jizixing Zhu -- replace the head/mixed control struct
%           with required thetaEqN; mixed form is now the sole formulation.
%   Change: 2026-07-13 Jizixing Zhu -- select the shared scalar Thomas core
%           for the new transactional path; legacy direct callers remain on
%           their unchanged compatibility branch.
%   Change: 2026-07-13 Jizixing Zhu -- consume the coefficient and common-FEM
%           matrix banks already built at Picard refresh; omitted optional
%           inputs retain the direct-kernel compatibility behavior.
%   Change: 2026-07-13 Jizixing Zhu -- marshal gas pressure mobility through
%           the Kg field without changing the water solve.
%   Change: 2026-07-14 Jizixing Zhu -- expose an assembly-only mode and the
%           post-boundary water linear system for correction-coordinate solves.
%   Change: 2026-07-17 Jizixing Zhu -- translate the atmospheric surface-water
%           active set without collapsing its supply and evaporation drivers.
%   Change: 2026-07-17 Jizixing Zhu -- retain legacy infiltration candidate fields
%           only inside this kernel adapter; the coupled solver accepts and
%           publishes the realized QMT/QMB and surface-reservoir ledger.
%   Change: 2026-07-20 Jizixing Zhu -- marshal the detached surface-water candidate
%     rather than an intermediate accepted-state write.
    if nargin < 4 || isempty(waterBoundaryFlux); waterBoundaryFlux = built.workspace.flux; end
    assert(nargin >= 5 && ~isempty(thetaEqN), ...
        'solveWaterTrial:missingThetaEqN', ...
        'Accepted-state thetaEqN is required by the mixed water equation.');
    if nargin < 6 || isempty(waterRootSink); waterRootSink = built.workspace.flux.Srt; end        % this-step root sink; default = committed Srt (backward-compat)
    if nargin < 7; boundaryRequest = []; end
    if nargin < 8; waterCoefficients = []; end
    if nargin < 9; waterMatrix = []; end
    if nargin < 10 || isempty(assembleOnly); assembleOnly = false; end
    [SV, TC, cfg, BC, ex, boundaryApplied] = assembleWaterSolveFrame( ...
        built, waterBoundaryFlux, waterRootSink, boundaryRequest);
    % The legacy kernel indexes its per-step surface forcing (Throughfall/EVAP) AND its per-step outputs
    % (QMT/QMB/Infiltration/Infiltration_act/TopFlux_potential) by the positional index KT. The new chain
    % marshals a CURRENT-STEP frame: ex.Throughfall/ex.EVAP are length-1 [1x1] scalars (allocateWorkspace
    % fx.Throughfall/EVAP = zeros(1,nCell)), NOT the full tS time series. So the legacy KT here is a LOCAL
    % index of 1, NOT the global outer-step kt (which stays with runOneStep/Registry/diagnostics and is
    % never forwarded into this frame). Passing the global kt read Throughfall(kt) out of bounds at
    % kt>=2. KT is purely positional in the kernel (no physical/first-step meaning; the only KT==1
    % branches are lDEBUG_NAN_TRACE-gated traces, off here), so mapping it to 1 is exact.
    localKt = 1;                                 % [-] positional index into the current-step legacy frame
    cfg.AssembleOnly = logical(assembleOnly);
    assert(isscalar(ex.Throughfall) && isscalar(ex.EVAP), ...
        ['solveWaterTrial: single-column only -- ex.Throughfall/EVAP must be [1x1]. A regional ' ...
         '[1 x nCell>1] frame would silently read cell-1 as the current step; per-cell looping is not wired.']);
    [SV, ~, CHKWATER, QMT, QMB, AVAIL0, Throughfall, Infiltration, Infiltration_act, ...
        TopFlux_potential, SAVE, runoffHort, linearSystem] = soil.water.solveSoilMoistureBalance( ...
        SV, TC, cfg, BC, dt, localKt, ex.SAVE, ex.Throughfall, ex.EVAP, ex.hN, ex.DSTMAX, ...
        ex.h_target, ex.NBChh, ex.hcrit, thetaEqN, waterCoefficients, waterMatrix);
    if ~assembleOnly
        built.workspace.step.trial.h  = SV.hh;  % [cm] the solved matric-head iterate
    end
    built.workspace.solver.SAVE   = SAVE;       % [-]  shared water/thermal/gas boundary-flux snapshot
    % WBC-1: pack the caller-facing scalar of each water-boundary direct output. The KT-indexed outputs
    % (QMT/QMB/Throughfall/Infiltration/Infiltration_act/TopFlux_potential) are read back at localKt (=1),
    % the SAME index the kernel wrote; AVAIL0 and runoffHort are already scalars. RETURN ONLY -- nothing
    % here is scattered to built.
    waterFlux = struct( ...
        'QMT',               QMT(localKt), ...            % [cm/s] top matric flux
        'QMB',               QMB(localKt), ...            % [cm/s] bottom matric flux
        'AVAIL0',            AVAIL0, ...                  % [cm/s] available surface supply (scalar)
        'Throughfall',       Throughfall(localKt), ...    % [cm/s] canopy throughfall to soil top
        'Infiltration',      Infiltration(localKt), ...   % [cm/s] potential infiltration
        'Infiltration_act',  Infiltration_act(localKt), ...% [cm/s] actual infiltration (net of runoff)
        'TopFlux_potential', TopFlux_potential(localKt), ...% [cm/s] top BC potential flux
        'runoffHort',        runoffHort);                 % [cm/s] Hortonian surface runoff (scalar)
end

% -----------------------------------------------------------------------------
function [SV, TC, cfg, BC, ex, boundaryApplied] = assembleWaterSolveFrame(built, wbFlux, waterRootSink, boundaryRequest)
% Assemble the legacy water-solve frame FROM built. wbFlux supplies the this-step water-boundary forcing
% (Throughfall/EVAP), decoupled from the last-committed built.workspace.flux. Private to the adapter; the oracle uses a SEPARATE
% assembler (wofixture.waterSolveLegacyFrame) so the bit-exact oracle catches a marshaling bug. Field
% list = the transitive reads of solveSoilMoistureBalance (calculateMatricCoefficients scalar +
% assembleCoefficientMatrices + calculateTimeDerivatives + calculateBoundaryConditions +
% solveTridiagonalMatrixEquations). Node-vector fields use (:); layer fields stay [NL x nD]. dt/kt are
% NOT needed here (solveWaterTrial passes them to the kernel directly).
    tr = built.workspace.step.trial;
    d  = built.workspace.derived;      fl = built.workspace.flux;
    Y  = built.Y;                      S  = built.parameters.soil;
    Y.surface_water = surface.currentWaterState(built);
    sw = built.control.switches;

    SV = struct();
    % state / iterate (node vectors)
    % Legacy time levels: hOLD/TOLD <- workspace.step.previous (A3 only), h/T <- accepted Y,
    % hh/TT <- Picard trial. The Richards time term must use Y as the current substep start.
    SV.hh = tr.h(:);  SV.h = Y.h(:);  SV.TT = tr.T(:);  SV.T = Y.T(:);
    SV.h_frez = tr.h_frez(:);  SV.hh_frez = tr.h_frez(:);  SV.P_gg = tr.p_gas(:);
    % vapor / correction aux (node vectors from the refresh)
    SV.RHOV = d.RHOV(:);  SV.DRHOVh = d.DRHOVh(:);  SV.DRHOVT = d.DRHOVT(:);
    SV.COR = d.COR(:);    SV.CORh = d.CORh(:);
    % water-content layers ([NL x nD]): *_L/_U = step-start (carried), *_LL/_UU = current (refreshed)
    SV.Theta_LL = d.Theta_LL;  SV.Theta_L = Y.Theta_L;
    SV.Theta_UU = d.Theta_UU;  SV.Theta_U = d.Theta_UU;   % Y carries no Theta_U -> use current (upper secant degenerate under no-freeze)
    SV.DTheta_LLh = d.DTheta_LLh;  SV.DTheta_UUh = d.DTheta_UUh;
    SV.Theta_V = d.Theta_V;  SV.Theta_II = d.Theta_II;
    SV.KfL_h = d.KfL_h;  SV.KL_h = d.KL_h;
    SV.DThetaLiquidDh = d.DThetaLiquidDh;
    SV.DThetaLiquidDT = d.DThetaLiquidDT;
    SV.DThetaWaterEquivDh = d.DThetaWaterEquivDh;
    SV.DThetaWaterEquivDT = d.DThetaWaterEquivDT;
    SV.DThetaIceDh = d.DThetaIceDh;  SV.DThetaIceDT = d.DThetaIceDT;
    SV.DThetaGasDh = d.DThetaGasDh;  SV.DThetaGasDT = d.DThetaGasDT;
    SV.Srt = waterRootSink;   % this-step root sink (caller defaults it to fl.Srt for backward-compat)

    TC = struct('D_V', d.D_V, 'D_Vg', d.D_Vg, 'DTheta_LLT', d.DThetaLiquidDT, ...
        'Eta', d.Eta, 'KL_T', d.KL_T, 'D_Ta', d.D_Ta, 'Kg', d.Kg, 'V_A', fl.V_A);

    cfg = struct();
    activeBottomNode = 1;
    externalModflow = struct('GroundwaterCoupling', 0);
    if ~isempty(boundaryRequest) && isfield(boundaryRequest, 'bottom') ...
            && isfield(boundaryRequest.bottom, 'kind') ...
            && strcmp(boundaryRequest.bottom.kind, 'modflow_fixed_head')
        modflow = boundaryRequest.bottom.modflow;
        activeBottomNode = modflow.active_node;
        externalModflow = struct( ...
            'GroundwaterCoupling', 1, ...
            'indxBotmLayer', modflow.active_node, ...
            'indxBotmLayer_R', modflow.active_layer_top_down, ...
            'headBotmLayer', modflow.head_bottom_layer_cm, ...
            'topLevel', modflow.aquifer_top_elevation_cm, ...
            'soilLayerThickness', modflow.soil_layer_thickness_cm);
    end
    cfg.Grid = struct('NN', built.domain.NN, 'NL', built.domain.NL, ...
        'nD', built.domain.grid.nD, 'DeltZ', built.domain.dz, ...
        'indxBotmLayer', activeBottomNode);
    cfg.MatricKernel = 'conservative';
    cfg.AllowPositiveHead = true;
    cfg.UseCommonThomas = true;
    cfg.lDualPor   = getFlag(sw, 'lDualPor', 0);
    cfg.lhThmrl    = getFlag(sw, 'lhThmrl', 1);
    cfg.lThmrlefc  = getFlag(sw, 'lThmrlefc', 0);
    cfg.lSoilairefc = getFlag(sw, 'lSoilairefc', 0);
    cfg.lKLa_Switch = getFlag(sw, 'lKLa_Switch', 0);
    cfg.lDVa_Switch = getFlag(sw, 'lDVa_Switch', 0);
    cfg.lDEBUG_NAN_TRACE = false;  cfg.lDEBUG_RAIN_BC = false;

    b = built.boundary;
    BC = struct('NBCh', b.NBCh_top, 'NBChB', b.NBCh_bottom, 'BCh', b.h_top, 'BChB', b.h_bottom, ...
        'Ks0', S.ks(end), 'DSTOR0', Y.surface_water.DSTOR0, ...
        'SurfaceRegime', 'legacy', ...
        'GroundwaterSettings', externalModflow);

    su = built.workspace.surface;
    ex = struct('SAVE', built.workspace.solver.SAVE, 'Throughfall', wbFlux.Throughfall, 'EVAP', wbFlux.EVAP, ...
        'hN', su.hN, 'DSTMAX', su.DSTMAX, 'h_target', su.h_target, 'NBChh', su.NBChh, 'hcrit', su.hcrit);
    boundaryApplied = struct();
    if isempty(boundaryRequest)
        assertNoUnmaterializedSeries(b);
    else
        [BC, ex, boundaryApplied] = applyBoundaryRequest(BC, ex, boundaryRequest);
    end
end

% -----------------------------------------------------------------------------
function [BC, ex, applied] = applyBoundaryRequest(BC, ex, request)
% Translate the canonical public request into the legacy kernel's sign/codes.
    validateBoundaryRequest(request);
    topKind = request.top.kind;

    switch request.top.kind
        case 'fixed_head'
            BC.NBCh = 1;
            BC.BCh = request.top.head_cm;
            ex.hN = request.top.head_cm;   % calculateBoundaryConditions uses hN for top Dirichlet
            ex.Throughfall = 0;
            ex.EVAP = 0;
        case 'fixed_flux'
            BC.NBCh = 2;
            BC.BCh = -request.top.requested_flux_into_domain_cm_s;   % legacy top flux is upward-positive
            ex.NBChh = 2;   % prescribed flux must not be replaced by the natural dry/ponded switch
            ex.Throughfall = 0;
            ex.EVAP = 0;
        case 'zero_flux'
            BC.NBCh = 4;
            BC.BCh = 0;
            ex.Throughfall = 0;
            ex.EVAP = 0;
        case 'atmospheric'
            BC.NBCh = 3;
            BC.BCh = 0;
            ex.Throughfall = request.top.liquid_supply_cm_s;
            ex.EVAP = request.top.evap_demand_cm_s;
            if isfield(request.top, 'surface')
                active = request.top.surface;
                BC.SurfaceRegime = active.regime;
                BC.DSTOR0 = active.storage_start_cm;
                ex.DSTMAX = active.capacity_cm;
                switch active.regime
                    case 'supply_limited'
                        ex.NBChh = 2;
                        topKind = 'fixed_flux';
                    case 'storage'
                        ex.NBChh = 2;
                        topKind = 'surface_storage';
                    case 'capacity_limited'
                        ex.NBChh = 1;
                        ex.h_target = active.capacity_cm;
                        ex.hN = active.capacity_cm;
                        topKind = 'fixed_head';
                    case 'dry_evaporation_limited'
                        ex.NBChh = 1;
                        ex.h_target = active.dry_head_limit_cm;
                        ex.hN = active.dry_head_limit_cm;
                        topKind = 'fixed_head';
                    otherwise
                        error('solveWaterTrial:badSurfaceRegime', ...
                            'Unknown atmospheric surface regime "%s".', ...
                            active.regime);
                end
            elseif ex.NBChh == 1
                topKind = 'fixed_head';
            else
                topKind = 'fixed_flux';
            end
        otherwise
            unsupportedRequest('top', request.top.kind);
    end

    switch request.bottom.kind
        case 'modflow_fixed_head'
            BC.NBChB = 1;
        case 'fixed_head'
            BC.NBChB = 1;
            BC.BChB = request.bottom.head_cm;
        case 'fixed_flux'
            BC.NBChB = 2;
            BC.BChB = request.bottom.requested_flux_into_domain_cm_s;   % legacy bottom sign already points into domain
        case 'zero_flux'
            BC.NBChB = 4;
            BC.BChB = 0;
        case 'free_drainage'
            BC.NBChB = 3;
            BC.BChB = 0;
        otherwise
            unsupportedRequest('bottom', request.bottom.kind);
    end
    applied = struct('top_kind', topKind, 'bottom_kind', request.bottom.kind);
end

% -----------------------------------------------------------------------------
function validateBoundaryRequest(request)
% Fail fast at the solver boundary; materialization bugs must not fall back.
    ok = isstruct(request) && isscalar(request) ...
        && isfield(request, 'schema_version') && isequal(request.schema_version, 1) ...
        && isfield(request, 'sign_convention') ...
        && strcmp(request.sign_convention, 'flux_into_domain_positive') ...
        && isfield(request, 'top') && isfield(request, 'bottom');
    assert(ok, 'solveWaterTrial:badWaterBoundaryRequest', ...
        'boundaryRequest must use canonical water-boundary request schema version 1.');
end

% -----------------------------------------------------------------------------
function unsupportedRequest(position, kind)
% The solver adapter supports only endpoint kinds represented by the kernel.
    error('solveWaterTrial:unsupportedWaterBoundaryRequest', ...
        'Unsupported %s water boundary request kind "%s".', position, kind);
end

% -----------------------------------------------------------------------------
function assertNoUnmaterializedSeries(boundary)
% Series policies require the W1b request-to-solver adapter; never use NaN static aliases.
    if ~isfield(boundary, 'water') || ~isstruct(boundary.water)
        return;
    end
    positions = {'top','bottom'};
    for i = 1:numel(positions)
        position = positions{i};
        if isfield(boundary.water, position) ...
                && isfield(boundary.water.(position), 'source') ...
                && strcmp(boundary.water.(position).source, 'series')
            error('solveWaterTrial:unmaterializedWaterBoundary', ...
                ['The %s water boundary uses a series source, but no per-trial ' ...
                 'BoundaryRequest is wired into solveWaterTrial yet. Materialize ' ...
                 'the request in the W1b scheduler before calling the solver.'], position);
        end
    end
end

% -----------------------------------------------------------------------------
function v = getFlag(sw, name, dflt)
    v = dflt;
    if isfield(sw, name) && ~isempty(sw.(name)), v = sw.(name); end
end
