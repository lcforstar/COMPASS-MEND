function candidate = reconcileRootNitrogenCandidate(candidate, requested_g_m2_s, actual_g_m2_s, dt_s)
%RECONCILEROOTNITROGENCANDIDATE Reconcile plant N with accepted soil uptake.
%   The plant candidate is generated before the soil N solve.  This adapter
%   applies the soil-side actual-minus-requested correction before the common
%   transaction is created, so plant and soil commit one N amount.
%   Winter wheat is replayed from accepted C/N input with the same frozen
%   demand and events and actual U_N [g N m^-2 s^-1]. This recomputes N
%   allocation, C:N limits and annual WSC supply together. A failed retry
%   repeats that accepted input; it never advances the requested candidate.
%   External litter/harvest amounts already consumed by soil must be identical.
%   Editor: Jizixing Zhu; Date: 2026-09-12.
% Change: 2026-09-30 Wenjie Liu -- actual crop uptake rematerializes the
% whole PlantCN candidate instead of withdrawing unavailable mobile N.
    assert(isstruct(candidate) && isfield(candidate,'plant_cn') ...
        && ~isempty(candidate.plant_cn), ...
        'vegetation:rootNReconcileCandidate','Plant-CN candidate is required.');
    assert(isscalar(requested_g_m2_s)&&isscalar(actual_g_m2_s)&&isscalar(dt_s) ...
        &&isfinite(requested_g_m2_s)&&isfinite(actual_g_m2_s) ...
        &&requested_g_m2_s>=0&&actual_g_m2_s>=0&&isfinite(dt_s)&&dt_s>0, ...
        'vegetation:rootNReconcileValues','Root-N reconciliation values are invalid.');
    delta=(actual_g_m2_s-requested_g_m2_s)*dt_s;
    nv=candidate.plant_cn;
    assert(isfield(nv,'N_labile')&&numel(nv.N_labile)>=2 ...
        &&isfield(nv,'cum_U_N'), ...
        'vegetation:rootNReconcileSchema','Plant candidate lacks N reconciliation fields.');
    replay = [];
    if isfield(candidate,'plant_cn_replay'),replay=candidate.plant_cn_replay;end
    isCropReplay = isstruct(replay)&&isfield(replay,'demand_ctx') ...
        &&replay.demand_ctx.plant_cn_params.resolved.winter_wheat_cn==1;
    if delta~=0 && isCropReplay
        required={'accepted_nv','env_turn','dt_s'};
        assert(all(isfield(replay,required))&&replay.dt_s==dt_s, ...
            'vegetation:rootNReplayContract', ...
            'Winter-wheat replay requires the same accepted input and time interval.');
        nvActual = replay.accepted_nv;
        nvActual.cum_U_N = nvActual.cum_U_N+actual_g_m2_s*dt_s;
        [nvActual,~] = vegetation.integrated.runPlantCNStep( ...
            nvActual,replay.demand_ctx,actual_g_m2_s,replay.env_turn,dt_s,2);
        external = {'cum_C_harvest','cum_N_harvest', ...
            'cum_C_leaf_litter','cum_N_leaf_litter','cum_C_root_litter','cum_N_root_litter', ...
            'cum_C_prune_litter','cum_N_prune_litter','cum_C_mortality','cum_N_mortality'};
        for j=1:numel(external)
            name=external{j};
            assert(isfield(nv,name)&&isfield(nvActual,name) ...
                &&isequal(nv.(name),nvActual.(name)), ...
                'vegetation:rootNReplayExternalTransfer', ...
                'Actual root-N replay changed already-consumed external plant ledger %s.',name);
        end
        nv=nvActual;
    else
        % Preserve the existing non-crop correction. For matching uptake,
        % leave the staged crop computation exactly intact.
        nv.N_labile(2)=nv.N_labile(2)+delta;
        nv.cum_U_N=nv.cum_U_N+delta;
    end
    assert(nv.N_labile(2)>=-1e-12&&nv.cum_U_N>=-1e-12, ...
        'vegetation:rootNReconcileNegative','Root-N reconciliation made plant N negative.');
    candidate.plant_cn=nv;
    candidate.root_nitrogen.actual_U_N_total=actual_g_m2_s;
    candidate.root_nitrogen.requested_U_N_total=requested_g_m2_s;
    candidate.root_nitrogen.unmet_U_N_total=max(requested_g_m2_s-actual_g_m2_s,0);
end
