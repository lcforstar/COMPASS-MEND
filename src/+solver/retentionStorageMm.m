function storage_mm = retentionStorageMm(built, dt_s, kt, firstNode)
%RETENTIONSTORAGEMM  Evaluate accepted soil water with solver constitutive law.
%
%   Inputs:
%     built [in] Built model at the accepted outer-step start.
%     dt_s  [s] Actual outer-step duration used for coefficient refresh.
%     kt    [-] One-based outer-step index.
%
%   Return:
%     firstNode [-, optional] First active Richards element; default 1.
%     storage_mm [mm] Liquid plus vapor-equivalent Richards-domain storage.
%
%   A value-semantics clone is refreshed at accepted Y so the diagnostic
%   ruler matches the mixed Richards endpoint relation without mutating model.
%   Editor: Jizixing Zhu; Date: 2026-07-16.
%   Change: 2026-07-20 Jizixing Zhu — restrict MODFLOW storage to the
%           externally bounded Richards control volume.
    if nargin < 4 || isempty(firstNode), firstNode = 1; end
    probe = built;
    probe.workspace.step.trial.h = built.Y.h;
    probe.workspace.step.trial.T = built.Y.T;
    if isfield(built.Y, 'p_gas') ...
            && isfield(probe.workspace.step.trial, 'p_gas')
        probe.workspace.step.trial.p_gas = built.Y.p_gas;
    end
    probe = soil.refreshSoilCoefficients(probe, struct( ...
        'L_f', built.workspace.step.L_f, 'dt_s', dt_s, 'KT', kt));
    derived = probe.workspace.derived;
    storage_mm = solver.richardsWaterStorageMm( ...
        derived.Theta_UU, derived.Theta_V, derived.RHOV, ...
        built.domain.dz, firstNode);
end
