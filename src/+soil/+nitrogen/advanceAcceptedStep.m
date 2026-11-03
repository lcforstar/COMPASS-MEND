function model = advanceAcceptedStep(model, kt, dt_outer_s) %#ok<INUSD>
%ADVANCEACCEPTEDSTEP  Reject the retired direct-mutating soil-N entrypoint.
%
%   Inputs:
%     model      [inout] Model with accepted water, heat, and dry-air state.
%     kt         [-] One-based outer-step index.
%     dt_outer_s [s] Accepted outer-step duration.
%
%   The production path is advanceAcceptedCandidate -> soil.commitCandidate,
%   so N pools, O2, history, and the ledger commit atomically with h/U/P_g.
%   This legacy mutating signature has no transaction argument and is kept
%   only to prevent a caller from bypassing that contract.
%   Editor: Jizixing Zhu; Date: 2026-07-16.
% Change: 2026-07-20 Jizixing Zhu -- supersede the temporary migration guard with
%   the candidate/transaction N bridge.
    error('soil:nitrogen:legacyDirectAdvanceForbidden', ...
        ['advanceAcceptedStep is retired. Use advanceAcceptedCandidate and ' ...
         'commit it through step.commitTransaction.']);
end
