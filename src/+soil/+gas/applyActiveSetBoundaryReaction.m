function [QAT, QAB, audit] = applyActiveSetBoundaryReaction( ...
        QAT, QAB, reaction, active, boundaryApplied)
%APPLYACTIVESETBOUNDARYREACTION Close a top-connected disappearance front.
%
%   A dry-air lower-bound row replaces the unconstrained conservation row
%   when its nodal mass would become negative. The positive KKT reaction is
%   the outward dry-air flux that the zero-mass state can no longer supply.
%   For an active front in the upper half of a column connected to a fixed-
%   pressure atmosphere, subtracting that rate from legacy outward-positive
%   QAT makes the applied boundary flux consistent with the constrained
%   state. QAB is unchanged. The correction is derived from the frozen
%   conservative row; it is never inferred from the later ledger residual.
%
%   This first production contract intentionally supports only the observed
%   upper, atmosphere-connected disappearance branch. A lower or disconnected
%   active set fails explicitly instead of silently choosing a flux route.
%
%   Inputs/outputs use [g cm^-2 s^-1].
%   Editor: Jizixing Zhu; Date: 2026-08-12.

    reaction = reaction(:);
    active = logical(active(:));
    n = numel(active);
    assert(numel(reaction) == n && all(isfinite(reaction)) ...
            && isscalar(QAT) && isfinite(QAT) ...
            && isscalar(QAB) && isfinite(QAB), ...
        'applyActiveSetBoundaryReaction:badInput', ...
        'Fluxes/reaction must be finite and active shape-compatible.');
    audit = struct('applied', false, 'active_node_count', nnz(active), ...
        'suppressed_top_outflow_g_cm2_s', 0);
    if ~any(active)
        return;
    end
    assert(isstruct(boundaryApplied) ...
            && isfield(boundaryApplied, 'top_kind') ...
            && strcmp(boundaryApplied.top_kind, 'fixed_pressure') ...
            && all(find(active) > floor(n/2)), ...
        'applyActiveSetBoundaryReaction:unsupportedConnectivity', ...
        ['Dry-air disappearance flux routing currently requires every ', ...
         'active node in the upper half and a fixed-pressure top boundary.']);
    activeReaction = reaction(active);
    reactionScale = max(max(abs(reaction)), 1);
    assert(all(activeReaction >= -128 * eps(reactionScale)), ...
        'applyActiveSetBoundaryReaction:wrongDualSign', ...
        ['A bound-active dry-air row must have a nonnegative outward-flow ', ...
         'suppression reaction.']);
    suppressedOutflow = sum(activeReaction);
    QAT = QAT - suppressedOutflow;
    audit.applied = true;
    audit.suppressed_top_outflow_g_cm2_s = suppressedOutflow;
end
