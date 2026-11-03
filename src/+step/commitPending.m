function model = commitPending(model, pending, kt, soil_candidate)
%COMMITPENDING  Commit all accepted outer-step candidates in fixed order.
%
%   Inputs:
%     model   [inout] Model transaction shell before the accepted commit.
%     pending [in] Destination-keyed candidate state from step.emptyPending.
%     kt      [-] Accepted one-based outer-step index.
%     soil_candidate [struct, optional] Explicit coupled-soil endpoint.  When
%             omitted, it is captured from the legacy solver workspace.
%
%   Return:
%     model   [out] Model with prognostic, flux, diagnostic, and management
%             candidates committed.
%
%   Commit order is load-bearing. The function is called only after the soil
%   scheduler completes the full outer step. A rejected step discards pending
%   and restores its snapshot instead of calling this function.
%   Editor: Jizixing Zhu; Date: 2026-07-16.
%   Change: 2026-07-17 Jizixing Zhu -- commit snow prognostics and diagnostics
%     only after the coupled outer step is accepted.
%   Change: 2026-07-18 Jizixing Zhu -- commit SCOPE force-restore history only after
%     the complete outer step is accepted.
%   Change: 2026-07-20 Jizixing Zhu -- accept the soil-owned candidate through
%     the transaction schema instead of implicitly reading mutable trial data.
%   Change: 2026-07-20 Jizixing Zhu -- retain this legacy helper as a typed-ledger
%     adapter; production runOneStep uses step.commitTransaction directly.
    step.validatePending(pending);

    if nargin < 4 || isempty(soil_candidate)
        soil_candidate = interfaces.captureSoilCandidate(model.built);
    end
    model.built = soil.commitCandidate(model.built, soil_candidate);
    model.built = surface.commitCandidate(model.built, pending.candidate.surface);
    model.built = vegetation.commitCandidate( ...
        model.built, pending.candidate.vegetation, kt);
    legacyLedger = struct( ...
        'surface_to_soil', struct( ...
            'water_request', pending.request.surface_to_soil.water), ...
        'vegetation_to_soil', struct( ...
            'water_request', pending.request.vegetation_to_soil.water, ...
            'nitrogen_request', pending.request.vegetation_to_soil.nitrogen), ...
        'vegetation_to_atmosphere', ...
            pending.ledger.vegetation_to_atmosphere);
    model = step.commitAncillary(model, legacyLedger, struct( ...
        'diagnostics', pending.diagnostics, ...
        'management', pending.management), kt);
end
