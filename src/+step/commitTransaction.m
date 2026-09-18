function model = commitTransaction(model, transaction, kt)
%COMMITTRANSACTION Atomically install validated multi-domain candidates.
%
%   Inputs:
%     model       [inout] Model shell after a completed outer solve.
%     transaction [struct] Destination-owned candidate and ledger payload.
%     kt          [-] Accepted one-based outer-step index.
%
%   Returns:
%     model       [out] Model with soil, vegetation, and surface candidates
%                 installed together.
%
%   Physical contract:
%     The routine changes no equation and recomputes no flux.  It performs the
%     discrete acceptance operation only after every candidate is present:
%     soil h/U/P_g endpoint first, then surface and vegetation candidates.
%     Diagnostics and management are write-only ancillary records.  A rejected
%     outer step never calls here.
%
%   Editor: Jizixing Zhu; Date: 2026-07-20.
%   Change: 2026-07-20 Jizixing Zhu -- project compatibility outputs from the
%     validated typed ledger after all domain candidates commit.
    interfaces.validateTransaction(transaction);
    model.built = soil.commitCandidate(model.built, transaction.candidate.soil);
    model.built = surface.commitCandidate( ...
        model.built, transaction.candidate.surface);
    model.built = vegetation.commitCandidate( ...
        model.built, transaction.candidate.vegetation, kt);
    model = step.commitAncillary( ...
        model, transaction.ledger, transaction.ancillary, kt);
end
