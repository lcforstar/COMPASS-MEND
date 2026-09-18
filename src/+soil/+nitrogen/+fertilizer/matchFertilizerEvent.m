function [should_apply, F_conc, fR, applied, pools] = matchFertilizerEvent(schedule, current_date, applied, dt_seconds)
%MATCHFERTILIZEREVENT  Pure fertilizer-event matcher with explicit once-only state.
%
%   [should_apply, F_conc, fR, applied, pools] = ...
%       soil.nitrogen.fertilizer.matchFertilizerEvent(schedule, current_date, applied)
%   checks whether any scheduled fertilizer event matches current_date (within a
%   1-hour window and not yet applied) and returns the UPDATED once-only flags plus
%   TWO views of the matched N:
%     - legacy aggregate: F_conc (summed concentration) + fR (summed form ratios),
%       consumed by the legacy soil.nitrogen.applyFertilizer / initFertilizerController;
%     - per-pool `pools`: strict concentration accumulated into each model state pool
%       (nh4_granule / nh4_liquid / no3_liquid / urea), used by the new-chain injector
%       soil.nitrogen.fertilizer.injectFertilizerPools. Per-pool accumulation stays
%       exact even when several events coincide with different pools/ratios (the
%       aggregate F_conc+fR cannot express that split -- hence the per-pool view).
%
%   Unlike the legacy closure (initFertilizerController), the once-only state
%   (`applied`) is an explicit in/out argument -- no hidden mutable closure. The
%   caller stores it in serializable runtime state (runtime.management.fertilizer),
%   so a checkpoint/resume does not re-fire past events. Maps directly to a Julia
%   pure function on the management state.
%
%   Inputs:
%     schedule     struct    .t_fert [datetime 1xn], .fert_kg_ha [1xn], .NH4_ratio/
%                            .NO3_ratio/.urea_ratio [1xn], .depth_cm [scalar cm],
%                            .pool_entry [1xn cellstr] target pool per event
%     current_date datetime  current outer-step timestamp
%     applied      [1xn logical]  which events have already fired
%     dt_seconds   optional positive step duration [s]. When supplied,
%                  current_date is the step start and events match [start,end).
%                  Omitted retains the legacy one-hour timestamp matcher.
%
%   Returns:
%     should_apply logical
%     F_conc       [scalar]  combined event concentration (fertilizerFlux sum) -- legacy
%     fR           struct(.NH4_ratio,.NO3_ratio,.urea_ratio) summed form ratios -- legacy
%     applied      [1xn logical]  updated once-only flags
%     pools        struct(.nh4_granule,.nh4_liquid,.no3_liquid,.urea) per-pool
%                  concentration [ug/cm3] to inject; all zero when no pool_entry field
%                  (legacy schedule) or no match.
%
%   Reference: extracted (pure) from soil.nitrogen.fertilizer.initFertilizerController.
%
%   Editor: Jizixing Zhu; Date: 2026-07-02.
%   Change: 2026-07-02 Jizixing Zhu — return pool_entry (nh4_granule vs nh4_liquid routing).
%   Change: 2026-07-02 Jizixing Zhu — replace pool_entry 5th output with per-pool `pools` (strict multi-event allocation; aggregate F_conc/fR kept for legacy).
% Change: 2026-09-16 Jizixing Zhu -- explicit half-open production event interval.
    use_interval=nargin>=4;
    if use_interval
        validateattributes(dt_seconds,{'numeric'},{'scalar','real','finite','positive'});
    end
    F_conc = 0;
    should_apply = false;
    fR = struct('NH4_ratio', 0, 'NO3_ratio', 0, 'urea_ratio', 0);
    pools = struct('nh4_granule', 0, 'nh4_liquid', 0, 'no3_liquid', 0, 'urea', 0);
    has_pool = isfield(schedule, 'pool_entry');
    eps_tol = hours(1);
    for i = 1:numel(schedule.t_fert)
        if use_interval
            matches=schedule.t_fert(i)>=current_date && ...
                schedule.t_fert(i)<current_date+seconds(dt_seconds);
        else
            matches=current_date>=schedule.t_fert(i) && ...
                current_date<schedule.t_fert(i)+eps_tol;
        end
        if matches && ~applied(i)
            Fi = soil.nitrogen.fertilizer.fertilizerFlux(schedule.fert_kg_ha(i), schedule.depth_cm);
            % legacy aggregate view: F_conc = sum of fluxes; fR accumulates the
            % FLUX-WEIGHTED ratio numerator (Fi*ratio) so applyFertilizer's F_conc*fR
            % reproduces sum_i(Fi*ratio_i) even for coincident events. Normalised below.
            F_conc = F_conc + Fi;
            fR.NH4_ratio  = fR.NH4_ratio  + Fi * schedule.NH4_ratio(i);
            fR.NO3_ratio  = fR.NO3_ratio  + Fi * schedule.NO3_ratio(i);
            fR.urea_ratio = fR.urea_ratio + Fi * schedule.urea_ratio(i);
            % per-pool view (strict): the NH4 fraction routes by this event's pool_entry
            % (nh4_granule vs nh4_liquid); NO3/urea fractions go to their fixed pools.
            if has_pool && numel(schedule.pool_entry) >= i
                pools = accumulatePools(pools, Fi, schedule, i);
            end
            should_apply  = true;
            applied(i)    = true;
        end
    end
    % Convert the weighted numerator to a flux-weighted mean ratio. Single event ->
    % fR == that event's ratios (F_conc*fR = Fi*ratio); coincident events no longer
    % double-count. F_conc = 0 leaves fR at zero (no match).
    if F_conc > 0
        fR.NH4_ratio  = fR.NH4_ratio  / F_conc;
        fR.NO3_ratio  = fR.NO3_ratio  / F_conc;
        fR.urea_ratio = fR.urea_ratio / F_conc;
    end
end

% ---------------------------------------------------------------------------
function pools = accumulatePools(pools, Fi, schedule, i)
% Route event i's N into per-pool concentrations. buildSources guarantees an
% NH4-bearing event (NH4_ratio>0) has an NH4 pool_entry; the otherwise-branch is a
% defensive fail-fast for an unvalidated schedule.
% Editor: Jizixing Zhu; Date: 2026-07-02.
    nh4_amt = Fi * schedule.NH4_ratio(i);
    if nh4_amt > 0
        switch lower(strtrim(schedule.pool_entry{i}))
            case 'nh4_granule', pools.nh4_granule = pools.nh4_granule + nh4_amt;
            case 'nh4_liquid',  pools.nh4_liquid  = pools.nh4_liquid  + nh4_amt;
            otherwise
                error('matchFertilizerEvent:nh4PoolMismatch', ...
                    ['event %d carries NH4 (NH4_ratio=%.3g) but pool_entry="%s" is not an ', ...
                     'NH4 pool (nh4_granule/nh4_liquid).'], i, schedule.NH4_ratio(i), schedule.pool_entry{i});
        end
    end
    pools.no3_liquid = pools.no3_liquid + Fi * schedule.NO3_ratio(i);
    pools.urea       = pools.urea       + Fi * schedule.urea_ratio(i);
end
