function candidate = reconcileRootNitrogenCandidate(candidate, requested_g_m2_s, actual_g_m2_s, dt_s)
%RECONCILEROOTNITROGENCANDIDATE Reconcile plant N with accepted soil uptake.
%   The plant candidate is generated before the soil N solve.  This adapter
%   applies the soil-side actual-minus-requested correction before the common
%   transaction is created, so plant and soil commit one N amount.
%   Editor: Jizixing Zhu; Date: 2026-09-12.
    assert(isstruct(candidate) && isfield(candidate,'plant_cn') ...
        && ~isempty(candidate.plant_cn), ...
        'vegetation:rootNReconcileCandidate','Plant-CN candidate is required.');
    assert(isfinite(requested_g_m2_s)&&isfinite(actual_g_m2_s) ...
        &&requested_g_m2_s>=0&&actual_g_m2_s>=0&&isfinite(dt_s)&&dt_s>0, ...
        'vegetation:rootNReconcileValues','Root-N reconciliation values are invalid.');
    delta=(actual_g_m2_s-requested_g_m2_s)*dt_s;
    nv=candidate.plant_cn;
    assert(isfield(nv,'N_labile')&&numel(nv.N_labile)>=2 ...
        &&isfield(nv,'cum_U_N'), ...
        'vegetation:rootNReconcileSchema','Plant candidate lacks N reconciliation fields.');
    nv.N_labile(2)=nv.N_labile(2)+delta;
    nv.cum_U_N=nv.cum_U_N+delta;
    assert(nv.N_labile(2)>=-1e-12&&nv.cum_U_N>=-1e-12, ...
        'vegetation:rootNReconcileNegative','Root-N reconciliation made plant N negative.');
    candidate.plant_cn=nv;
    candidate.root_nitrogen.actual_U_N_total=actual_g_m2_s;
    candidate.root_nitrogen.requested_U_N_total=requested_g_m2_s;
    candidate.root_nitrogen.unmet_U_N_total=max(requested_g_m2_s-actual_g_m2_s,0);
end
