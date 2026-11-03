function decision = evaluateSurfaceTrial( ...
        request, surfaceState, massToleranceCm, trialFlags)
%EVALUATESURFACETRIAL  Test atmospheric complementarity and close surface water.
%
%   QMT is positive outward from the soil. Let d be retained liquid depth,
%   h the top-node pressure head, r the overflow reaction, and qin=-QMT the
%   signed inward soil flux. A valid accepted tile satisfies
%
%     d1 - d0 = dt * (P_liq - E - qin - r),
%     0 <= d1  _|_  d1-h >= 0,
%     0 <= r   _|_  dmax-d1 >= 0.
%
%   The two complementarity pairs distinguish dry supply limitation, retained
%   ponding below capacity, and capacity-limited overflow. An invalid active
%   set returns the next regime and no state is committed. A valid result
%   returns endpoint liquid depth, runoff, mass residual, and KKT diagnostics.
%
%   Inputs:
%     request          immutable atmospheric boundary request for one trial
%     surfaceState     exact candidate record with head_cm [cm],
%                      dry_state_depleted [-], and water_flux.QMT [cm s-1,
%                      positive outward from soil]
%     massToleranceCm  accepted surface ledger tolerance [cm per trial]
%     trialFlags       projection and dry-head target audit flags
%
%   The function returns a decision only. It never mutates trial or accepted
%   state; the outer scheduler owns re-solve, commit, and rollback.
%
%   For capacity_cm=0 this reduces to the no-ponding pressure/flux
%   Signorini condition in Gatti et al. (2024), arXiv:2407.07865,
%   Sect. 2.2 Eq. (2). Positive capacity generalizes that condition with the
%   prognostic COMPASS DSTOR reservoir; no Nitsche or mixed-FEM terms enter.
%
%   Editor: Jizixing Zhu; Date: 2026-07-17.
%   Change: 2026-07-17 Jizixing Zhu -- document the finite-storage
%           generalization of the no-ponding Signorini boundary.
%   Change: 2026-07-17 Jizixing Zhu -- express finite ponding as two
%           complementarity pairs and expose their accepted-trial residuals.
%   Change: 2026-07-17 Jizixing Zhu -- expose the external liquid-supply and
%           evaporation terms needed by the combined surface-soil ledger.
%   Change: 2026-07-17 Jizixing Zhu -- scale reaction feasibility with the
%           declared surface-water amount tolerance instead of an unrelated
%           machine-scale head tolerance.
%   Change: 2026-07-17 Jizixing Zhu -- make the finite hcrit atmospheric
%           drying constraint an explicit active set and distinguish realized
%           evaporation from the upstream SCOPE evaporation demand.
%   Change: 2026-07-17 Jizixing Zhu -- retain the flux set at the dry-head corner
%           when its supply deficit is already within the surface mass
%           tolerance, preventing tolerance-inconsistent active-set chatter.
%   Change: 2026-07-17 Jizixing Zhu -- activate the dry set from solved residual-
%           water state, not from pressure head alone, and return the
%           state-consistent head that the dry branch must retain.
%   Change: 2026-08-14 Jizixing Zhu -- require the
%           prescribed-flux solution to cross below the dry-head lower bound
%           before a residual-water endpoint may activate the dry set.
%   Change: 2026-08-14 Jizixing Zhu -- separate
%           the dry-head tolerance from the ponding-head tolerance, require
%           the residual-state condition on both entry and retention of the
%           dry set, and expose the signed dry-head gap for audit.
%   Change: 2026-08-15 Jizixing Zhu -- consume the exact coupled surface state and
%           water flux as one inseparable candidate record.
    if nargin < 3 || isempty(massToleranceCm)
        massToleranceCm = 1e-5;
    end
    if nargin < 4 || isempty(trialFlags)
        trialFlags = struct();
    end
    assert(isscalar(massToleranceCm) && isfinite(massToleranceCm) ...
            && massToleranceCm > 0, ...
        'evaluateSurfaceTrial:badTolerance', ...
        'massToleranceCm must be a positive finite scalar [cm].');
    decision = emptyDecision();
    if ~isstruct(request) || ~isfield(request, 'top') ...
            || ~strcmp(request.top.kind, 'atmospheric') ...
            || ~isfield(request.top, 'surface')
        decision.valid = true;
        return;
    end
    assert(isstruct(surfaceState) ...
            && isfield(surfaceState, 'head_cm') ...
            && isfield(surfaceState, 'dry_state_depleted') ...
            && isfield(surfaceState, 'water_flux') ...
            && isstruct(surfaceState.water_flux) ...
            && isfield(surfaceState.water_flux, 'QMT') ...
            && isfinite(surfaceState.head_cm) ...
            && isfinite(surfaceState.water_flux.QMT), ...
        'evaluateSurfaceTrial:badTrial', ...
        'A complete finite coupled surface state and actual QMT are required.');
    surfaceHeadCm = double(surfaceState.head_cm);
    waterFlux = surfaceState.water_flux;

    surface = request.top.surface;
    d0 = surface.storage_start_cm;
    dmax = surface.capacity_cm;
    dt = request.dt_s;
    supply = request.top.liquid_supply_cm_s;
    evap = request.top.evap_demand_cm_s;
    qSurface = supply - evap + waterFlux.QMT;
    headTol = max(1e-10, 1e-8 * max(1, dmax));
    fluxTol = max(1e-14, massToleranceCm / dt);
    activeDryHead = -Inf;
    if isfield(surface, 'dry_head_limit_cm') ...
            && ~isempty(surface.dry_head_limit_cm)
        activeDryHead = double(surface.dry_head_limit_cm);
    end
    dryHeadTol = dryHeadToleranceCm(surface, activeDryHead);
    dryHeadProjected = isfield(trialFlags, 'dry_head_projected') ...
        && logical(trialFlags.dry_head_projected);
    dryStateDepleted = logical(surfaceState.dry_state_depleted);
    dryHeadTarget = surfaceHeadCm;
    if isfield(trialFlags, 'dry_head_target_cm') ...
            && isfinite(trialFlags.dry_head_target_cm)
        dryHeadTarget = double(trialFlags.dry_head_target_cm);
    end

    decision.regime = surface.regime;
    decision.next_regime = surface.regime;
    decision.psi_top_cm = surfaceHeadCm;
    decision.liquid_supply_cm_s = supply;
    decision.evap_demand_cm_s = evap;
    decision.q_in_actual_cm_s = -waterFlux.QMT;
    decision.trial_balance_rate_cm_s = qSurface;
    decision.trial_dry_head_projected = dryHeadProjected;
    decision.trial_dry_state_depleted = dryStateDepleted;
    decision.dry_head_limit_cm = activeDryHead;
    decision.dry_head_tolerance_cm = dryHeadTol;
    if isfinite(activeDryHead)
        decision.dry_head_gap_cm = surfaceHeadCm - activeDryHead;
    end
    evapActual = evap;
    switch surface.regime
        case 'supply_limited'
            if surfaceHeadCm > headTol
                decision.next_regime = nextWetRegime(dmax, headTol);
                return;
            elseif isfinite(activeDryHead) && supply <= fluxTol ...
                    && evap > fluxTol && dryStateDepleted ...
                    && surfaceHeadCm < activeDryHead - dryHeadTol
                decision.next_regime = 'dry_evaporation_limited';
                decision.next_dry_head_limit_cm = ...
                    min(dryHeadTarget, activeDryHead);
                return;
            end
            d1 = 0;
            runoff = 0;
        case 'storage'
            if surfaceHeadCm < -headTol
                decision.next_regime = 'supply_limited';
                return;
            elseif surfaceHeadCm > dmax + headTol
                decision.next_regime = 'capacity_limited';
                return;
            end
            d1 = min(max(surfaceHeadCm, 0), dmax);
            runoff = 0;
        case 'capacity_limited'
            d1 = dmax;
            runoff = qSurface - (d1 - d0) / dt;
            if runoff < -fluxTol
                if dmax > headTol
                    decision.next_regime = 'storage';
                else
                    decision.next_regime = 'supply_limited';
                end
                return;
            end
            runoff = max(runoff, 0);
        case 'dry_evaporation_limited'
            if ~isfinite(activeDryHead)
                error('evaluateSurfaceTrial:missingDryHead', ...
                    ['The dry_evaporation_limited regime requires ', ...
                     'surface.dry_head_limit_cm.']);
            end
            if supply > fluxTol || ~dryStateDepleted ...
                    || abs(surfaceHeadCm - activeDryHead) > dryHeadTol
                decision.next_regime = 'supply_limited';
                return;
            end
            % QMT is the solved total-water reaction, positive outward.
            % With no retained liquid, surface conservation gives the
            % realized atmospheric evaporation E_actual=P_liq+QMT.  The
            % internal water operator still contains liquid capillary and
            % vapor transport; this active set does not make the top element
            % impermeable.  Rainfall releases the fixed dry head above.
            evapActual = supply + waterFlux.QMT;
            if evapActual < -fluxTol || evapActual > evap + fluxTol
                decision.next_regime = 'supply_limited';
                return;
            end
            evapActual = min(max(evapActual, 0), evap);
            d1 = 0;
            runoff = 0;
        otherwise
            error('evaluateSurfaceTrial:badRegime', ...
                'Unknown atmospheric surface regime "%s".', surface.regime);
    end

    decision.valid = true;
    decision.storage_end_cm = d1;
    decision.runoff_cm_s = runoff;
    decision.residual_cm = d1 - d0 ...
        - dt * (supply - evapActual + waterFlux.QMT - runoff);
    decision.q_available_cm_s = supply + d0 / dt;
    decision.evap_actual_cm_s = evapActual;
    decision.evap_unmet_cm_s = max(evap - evapActual, 0);
    decision.storage_tendency_cm_s = (d1 - d0) / dt;
    decision.excess_reaction_cm_s = supply - evapActual ...
        - decision.q_in_actual_cm_s - decision.storage_tendency_cm_s;
    decision.head_gap_cm = d1 - surfaceHeadCm;
    decision.capacity_gap_cm = dmax - d1;
    decision.lower_complementarity_cm2 = ...
        abs(d1 * decision.head_gap_cm);
    decision.upper_complementarity_cm2_s = ...
        abs(decision.excess_reaction_cm_s * decision.capacity_gap_cm);
    decision.mode_code = regimeCode(surface.regime);

    headProductTol = headTol * max(1, dmax);
    fluxProductTol = fluxTol * max(1, dmax);
    decision.complementarity_valid = d1 >= -headTol ...
        && decision.head_gap_cm >= -headTol ...
        && decision.capacity_gap_cm >= -headTol ...
        && decision.excess_reaction_cm_s >= -fluxTol ...
        && decision.lower_complementarity_cm2 <= headProductTol ...
        && decision.upper_complementarity_cm2_s <= fluxProductTol;
end

function regime = nextWetRegime(capacity, tol)
% A zero-capacity surface goes directly to the runoff constraint.
% Editor: Jizixing Zhu; Date: 2026-07-17.
    if capacity > tol
        regime = 'storage';
    else
        regime = 'capacity_limited';
    end
end

function tol = dryHeadToleranceCm(surface, activeDryHead)
% Scale the dry lower-bound test on its own pressure-head magnitude.
%
% The ponding tolerance is O(surface storage), whereas an accepted dry head
% may be O(1e6 cm). Reusing the ponding ruler makes the dry inequality depend
% on floating-point noise. The optional relative field exists only for this
% dry-head sensitivity gate; requests may override it only for controlled tests.
% Editor: Jizixing Zhu; Date: 2026-08-14.
    relativeTol = 1e-8;
    if isfield(surface, 'dry_head_relative_tolerance') ...
            && ~isempty(surface.dry_head_relative_tolerance)
        relativeTol = double(surface.dry_head_relative_tolerance);
    end
    assert(isscalar(relativeTol) && isfinite(relativeTol) ...
            && relativeTol > 0, ...
        'evaluateSurfaceTrial:badDryHeadTolerance', ...
        'dry_head_relative_tolerance must be positive and finite.');
    if isfinite(activeDryHead)
        tol = max(1e-10, relativeTol * max(1, abs(activeDryHead)));
    else
        tol = NaN;
    end
end

function decision = emptyDecision()
% Fixed audit schema for atmospheric and non-atmospheric requests.
% Editor: Jizixing Zhu; Date: 2026-07-17.
    decision = struct( ...
        'valid', false, ...
        'regime', '', ...
        'next_regime', '', ...
        'storage_end_cm', NaN, ...
        'runoff_cm_s', NaN, ...
        'residual_cm', NaN, ...
        'psi_top_cm', NaN, ...
        'q_available_cm_s', NaN, ...
        'liquid_supply_cm_s', NaN, ...
        'evap_demand_cm_s', NaN, ...
        'evap_actual_cm_s', NaN, ...
        'evap_unmet_cm_s', NaN, ...
        'q_in_actual_cm_s', NaN, ...
        'trial_balance_rate_cm_s', NaN, ...
        'trial_dry_head_projected', false, ...
        'trial_dry_state_depleted', false, ...
        'dry_head_limit_cm', NaN, ...
        'dry_head_tolerance_cm', NaN, ...
        'dry_head_gap_cm', NaN, ...
        'next_dry_head_limit_cm', NaN, ...
        'storage_tendency_cm_s', NaN, ...
        'excess_reaction_cm_s', NaN, ...
        'head_gap_cm', NaN, ...
        'capacity_gap_cm', NaN, ...
        'lower_complementarity_cm2', NaN, ...
        'upper_complementarity_cm2_s', NaN, ...
        'mode_code', NaN, ...
        'complementarity_valid', false);
end

function code = regimeCode(regime)
% Stable public code: 0 supply, 1 storage, 2 overflow, 3 residual dry state.
% Editor: Jizixing Zhu; Date: 2026-07-17.
    codes = struct('supply_limited', 0, 'storage', 1, ...
        'capacity_limited', 2, 'dry_evaporation_limited', 3);
    code = codes.(regime);
end
