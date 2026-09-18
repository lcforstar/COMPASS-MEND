function hash = hashTrialState(built)
%HASHTRIALSTATE Return a stable identity for one h/U/T/P_g trial state.
%
%   The hash is used only to prove that a residual/ledger record belongs to
%   the candidate currently under review. It is not a physical diagnostic
%   and never changes solver state.
%
%   Editor: Jizixing Zhu; Date: 2026-08-15.
%   Change: 2026-08-20 Jizixing Zhu -- hash the exact IEEE-754 state bytes
%     instead of formatting every value as decimal text. The identity remains
%     bit-exact within a run while avoiding allocation-heavy sprintf work.
    trial = built.workspace.step.trial;
    state = [double(trial.h(:)); double(trial.U(:)); ...
        double(trial.T(:)); double(trial.p_gas(:))];
    payload = typecast(state, 'uint8');
    digest = java.security.MessageDigest.getInstance('SHA-256');
    digest.update(payload);
    hash = lower(sprintf('%02x', typecast(digest.digest, 'uint8')));
end
