function [pools_out, diag, state_out] = turnoverAndHarvest( ...
    pools_in, state_in, pftParams, env, dt)
%TURNOVERANDHARVEST  Plant N pool turnover + Apple harvest + retransl + spring release (M2-B.4 主调度).
%
%   每 KT 调一次. 按 frPHU 触发 4 个事件:
%     (a) Apple fruit harvest (一次性, frPHU > harvest_threshold)
%     (b) Leaf senescence + retransl (每日, frPHU > DLAI_threshold)
%     (c) Fine-root turnover (每日, 全年)
%     (d) Spring N_storage release (每日, 0 < frPHU < spring_release_endphu)
%
%   Returns:
%     pools_out  struct  updated 6-pool plant N: N_labile, N_leaf, N_root_fine, N_livestem, N_deadstem, N_fruit, N_storage
%     diag       struct  this-KT diagnostic: harvest_today, leaf_litter_today, root_litter_today,
%                        storage_release_today
%     state_out  struct  updated state: harvested_this_year (1=already done), year_now
%
%   Inputs:
%     pools_in   struct  6 N pools (scalar per pool):
%                  .N_labile, .N_leaf, .N_root_fine, .N_stem_wood, .N_fruit, .N_storage  [g N/m² ground]
%     state_in   struct  state carried across KT:
%                  .harvested_this_year [-,{0,1}]  fruit-harvest event flag for this growing season
%                  .year_now             [int]      to detect year rollover (reset harvested flag)
%     pftParams  struct  from PFT table (M2-B.4 [M] group + [J] + [L]):
%                  .harvest_threshold, .harvest_fraction
%                  .root_turnover_rate
%                  .spring_release_endphu, .spring_release_days
%                  .leaf_long_days, .retransl_fraction  ([J] phenology)
%                  .is_perennial, .has_fruit, .has_storage_pool  ([B] attributes)
%     env        struct  current state:
%                  .frPHU  [-]   accumulated heat unit fraction
%                  .LAI    [-]   current LAI
%                  .DLAI   [-]   EPIC senescence threshold (= cfg.Veg.DLAI = 0.45)
%                  .year   [int] current calendar year (for harvested_this_year reset)
%     dt         [s]            outer time step (3600)
%
%   Variable class (per compass-coding-standard §1):
%     pools_in/pools_out — prognostic state (caller owns)
%     state_in/state_out — auxiliary (can be reconstructed from frPHU history)
%     pftParams          — parameters (immutable)
%     env                — forcing slice
%     diag               — diagnostics (only output)
%
%   Reference:
%     EPIC harvest index (Williams 1989 Trans.ASAE 32:497)
%     Aerts R (1996) Oikos 67:469 — leaf N resorption
%     Wells & Eissenstat (2002) JPGR 21:324 — Apple fine root turnover
%     cn_coupling_path_and_literature_audit.md §7 (M2-B.4 plan)
%
%   Author: Zhu Jizixing 2026-05-21 (M2-B.4 S2.5)
% Editor: Wenjie Liu; Date: 2026-10-22.
% Change: 2026-10-22 Wenjie Liu -- route 15% of actual winter-wheat straw N
% to the paired stubble AOM return; preserve seed and one-shot harvest latches.
% Change: 2026-09-30 Wenjie Liu -- annual crop harvest is one-shot, pairs
%   retained seed C/N, and accounts for remaining mobile/storage N in the
%   existing Changwu shoot-removal/root-litter management ledger.

    % --- Init outputs ---
    pools_out = pools_in;
    diag = struct('harvest_today', 0, 'leaf_litter_today', 0, ...
                  'root_litter_today', 0, ...
                  'storage_release_today', 0, ...        % total reserve->labile (0 since springStorageRelease removed)
                  'storage_release_root_today', 0, ...   % [dual-storage] root reserve -> labile (C-mirror)
                  'storage_release_stem_today', 0, ...   % [dual-storage] stem reserve -> labile (C-mirror)
                  'stem_turnover_today', 0);   % [13c] live->dead stem turnover
    state_out = state_in;

    if isfield(pftParams, 'winter_wheat_cn') && pftParams.winter_wheat_cn == 1
        [pools_out, diag, state_out] = winterWheatTurnoverN( ...
            pools_out, diag, state_out, pftParams, env);
        return;
    end

    % --- Re-arm the one-shot harvest each fruit season ---
    % [13b fix 2026-06-04] Re-arm via the phenology's own harvested_flag (reset to
    % false at COLDEST_DAY in seasonDecidPhenology), NOT the legacy EPIC frPHU.
    % frPHU is a MULTI-YEAR MONOTONIC accumulator (Changwu y1 end 0.648 -> y2 stays
    % 0.648 -> 1.0, never < 0.01), so the old frPHU<0.01 reset never fired in year-2+
    % and blocked every harvest after year 1. harvested_flag is true only from the
    % harvest day until the COLDEST_DAY reset, so resetting on ~harvested_flag re-arms
    % each new season AND keeps the flag set across the calendar year boundary (no
    % double-harvest).
    if isfield(env, 'phen_state') && isfield(env.phen_state, 'harvested_flag') ...
            && ~env.phen_state.harvested_flag
        state_out.harvested_this_year = 0;
    end

    % --- Skip if not perennial (crop C/N requires its own CoLM crop state machine) ---
    if ~isfield(pftParams, 'is_perennial') || pftParams.is_perennial ~= 1
        return;
    end

    % --- [Phase D Phenology 2026-05-26] Build trigger signals from PhenologyState ---
    % env.phen_state (PhenologyState struct) fields (CLM5+LPJ-GUESS schema):
    %   .offset_flag      true during CLM5 ndays_off=15 day senescence ramp
    %   .budburst_doy     set on Sykes 1996 bud break event
    %   .DS_fruit         [0,1] fruit development stage (apple-class only)
    %   .harvested_flag   true after DS_fruit > harvest_DS (one-shot per year)
    %   .harvest_doy      set on harvest event day
    %   .aphen            cumulative leaf-on days (LPJ-GUESS APHEN_MAX=210)
    % Source: src/+vegetation/+integrated/+phenology/seasonDecidPhenology.m
    %
    % [P5 fix] harvest uses DS_fruit (independent fruit clock), NOT aphen proxy.
    % [P6 fix] spring uses doy − budburst_doy window (day-based), NOT aphen window.
    if isfield(env, 'phen_state')
        ps = env.phen_state;
        % Harvest signal: triggered ONCE per year by independent DS_fruit clock
        % (P5 fix: separate fruit phenology from leaf phenology)
        harvest_signal = ps.harvested_flag * 1.0;
        % Senescence signal: CLM5 daylength offset flag (unchanged)
        sene_signal    = ps.offset_flag * 1.0;
        % [unify-spring-release 2026-06-06] spring_signal removed: spring remobilisation
        %   is now springOnsetXfer (budburst-window-driven in runPlantCNStep), not the
        %   deleted springStorageRelease that consumed this signal.
    else
        % Fallback: legacy frPHU path (preserves backward compat if dispatcher absent)
        harvest_signal = env.frPHU;
        sene_signal    = env.frPHU;
    end

    % --- (a) Apple fruit harvest (Phase D: aphen-driven) ---
    if isfield(pftParams, 'has_fruit') && pftParams.has_fruit == 1
        % Threshold 0.5: harvest_signal is 0 or 1, so threshold=0.5 triggers exactly once
        [pools_out.N_fruit, diag.harvest_today, harvested_now] = ...
            vegetation.integrated.turnover.appleFruitHarvest( ...
                pools_out.N_fruit, harvest_signal, ...
                0.5, pftParams.harvest_fraction, ...
                state_out.harvested_this_year);
        if harvested_now
            state_out.harvested_this_year = 1;
        end
    end

    % --- (b) Leaf senescence + retransl (CLM5 strict: ramp + last-day flush) ---
    % Source: CLM5 CNPhenologyMod.F90:808-820 (ramp) + 2448-2453 (flush residual)
    if pftParams.has_storage_pool == 1
        if isfield(pftParams, 'ndays_off') && ~isnan(pftParams.ndays_off)
            ndays_off_use = pftParams.ndays_off;
        else
            ndays_off_use = 15;
        end
        % Pass phen_state.offset_counter so last-day-flush trigger works
        offset_counter_use = [];
        if isfield(env, 'phen_state') && isfield(env.phen_state, 'offset_counter')
            offset_counter_use = env.phen_state.offset_counter;
        end
        % [leaf-flush 2026-06-08] If the cold-month reset interrupted the leaf-fall ramp
        %   (seasonDecidPhenology sets leaf_flush_pending and clears offset_flag in the same
        %   step), force the last-day flush: override sene_signal to 1 and offset_counter to 1
        %   so leafSenescenceRetransl drains ALL residual leaf C/N to litter this step instead
        %   of returning early on sene_signal=offset_flag=0 (which stranded ~24 gC over winter
        %   in v62). Physically = a hard freeze defoliating the tree in one step.
        if isfield(env, 'phen_state') && isfield(env.phen_state, 'leaf_flush_pending') ...
                && env.phen_state.leaf_flush_pending
            sene_signal        = 1.0;
            offset_counter_use = 1;     % triggers leafSenescenceRetransl drain-all path
        end
        % [resorbed-N split 2026-06-11] Leaf-resorbed N is stored across woody sites, not only stem.
        %   For apple, roots and bark are the MAIN sinks of foliar/autumn-resorbed N (Dong, Cheng,
        %   Scagel & Fuchigami 2002 Tree Physiol 22:1305-1310, p.1308); over-wintering N is stored in
        %   roots, branches and trunk (Millard & Grelet 2010 Tree Physiol 30:1083-1095, p.1085).
        %   leafSenescenceRetransl adds the resorbed N to the stem pool; we then move a root share
        %   (f_storage_root, the reserve root:stem partition) into the coarse-root reserve. This gives
        %   the root N reserve a physical inflow -- without it, since the old C-derived refill was
        %   removed (see runPlantCNStep Step 6), the root reserve would starve. The root share is a
        %   calibration parameter: Dong 2002 implies roots may warrant more than the wood-mass-based
        %   f_storage_root; revisit with a quantitative apple N-partition measurement.
        N_stem_storage_pre = pools_out.N_stem_storage;
        [pools_out.N_leaf, pools_out.N_stem_storage, diag.leaf_litter_today] = ...
            vegetation.integrated.turnover.leafSenescenceRetransl( ...
                pools_out.N_leaf, pools_out.N_stem_storage, sene_signal, env.LAI, ...
                ndays_off_use, pftParams.retransl_fraction, ...
                0.5, dt, offset_counter_use);
        N_retransl_to_store = pools_out.N_stem_storage - N_stem_storage_pre;   % [g N/m2] this-step leaf-resorbed N
        if N_retransl_to_store > 0 && isfield(pftParams, 'f_storage_root')
            N_retransl_to_root = pftParams.f_storage_root * N_retransl_to_store;   % [g N/m2] root share
            pools_out.N_stem_storage = pools_out.N_stem_storage - N_retransl_to_root;
            pools_out.N_root_storage = pools_out.N_root_storage + N_retransl_to_root;
        end
    end

    % --- (c) Fine-root turnover [F5 2026-05-31]: phenology-gated (CLM5 bglfr) ---
    % CLM5 CNPhenologyMod.F90:2623 ties frootc_to_litter to leaf bglfr (=0 in dormancy).
    % Proxy: phen_gate = ps.phen (leaf activity 0..1); dormant phen~0 -> turnover slows,
    %   growing season phen~1 -> full rate. Avoids year-round root shedding in winter.
    % Fallback (no phen_state): phen_gate=1 (year-round, backward compatible).
    if isfield(env, 'phen_state') && isfield(env.phen_state, 'phen')
        froot_phen_gate = max(0, min(1, env.phen_state.phen));
    else
        froot_phen_gate = 1.0;
    end
    [pools_out.N_root_fine, diag.root_litter_today] = ...
        vegetation.integrated.turnover.fineRootTurnover( ...
            pools_out.N_root_fine, pftParams.root_turnover_rate, dt, froot_phen_gate);

    % --- (c2) Fine-root OFFSET senescence [CLM5-exact 2026-06-13] -- die back with the leaf ---
    % CLM5 CNPhenologyMod.F90:3334 sheds fine root in LOCK-STEP with the leaf during the autumn offset
    %   (frootc_to_litter uses the SAME quadratic ramp as leafc_to_litter). For a season-deciduous tree
    %   the fine root dies back with the leaf -> ~0 in winter (Dombrowski 2022 Fig 4a; FruitTree fork
    %   season_decid path). COMPASS ties it to the leaf shed FRACTION this step: any leaf shed (b) carries
    %   the same fraction of fine root to litter. This makes the fine root SEASONAL (grow with leaf via
    %   afroot allocation, shed with leaf at senescence). Gated by has_storage_pool (deciduous tree with
    %   the leaf-offset machinery); N to litter (no retransl, conservative). C mirrors via root_litter_today.
    if pftParams.has_storage_pool == 1 && diag.leaf_litter_today > 0 && pools_out.N_root_fine > 0
        N_leaf_pre_sen = pools_out.N_leaf + diag.leaf_litter_today;     % leaf N before this-step senescence
        if N_leaf_pre_sen > 1e-12
            leaf_shed_frac          = min(1, diag.leaf_litter_today / N_leaf_pre_sen);
            dN_froot_off            = leaf_shed_frac * pools_out.N_root_fine;
            pools_out.N_root_fine   = pools_out.N_root_fine - dN_froot_off;
            diag.root_litter_today  = diag.root_litter_today + dN_froot_off;   % C follows via turnoverCarbonMirror
        end
    end

    % --- (c2) Live->dead stem turnover [13c, fixes N3]: livestem ages into deadstem ---
    % CLM5 CNPhenologyMod.F90:2754 ctovr = livestemc * lwtop (annual fraction per second).
    % N: dead wood only keeps CN_lw/CN_dw share (CN_dw>>CN_lw); rest retransl -> N_storage
    %   (CLM5 Veg_Phenology Tech Note:327 N retransl out of live wood during turnover).
    if isfield(pftParams,'lwtop') && pftParams.lwtop > 0 && pools_out.N_livestem > 0
        lwtop_frac = pftParams.lwtop * dt / 31536000;          % [1/yr] annual frac -> per dt
        dN_lw      = pools_out.N_livestem * lwtop_frac;        % [g N/m2] N leaving live stem
        dN_to_dead = dN_lw * (pftParams.CN_lw / pftParams.CN_dw);   % dead wood keeps fewer N (high C:N)
        dN_retrans = dN_lw - dN_to_dead;                       % rest recycled to storage
        pools_out.N_livestem = pools_out.N_livestem - dN_lw;
        pools_out.N_deadstem = pools_out.N_deadstem + dN_to_dead;
        % [dual-storage 2026-06-06] live-wood retransl N -> STEM reserve (stem source)
        pools_out.N_stem_storage = pools_out.N_stem_storage + dN_retrans;
        diag.stem_turnover_today = dN_lw;                      % [g N/m2] for C mirror + diag
    else
        diag.stem_turnover_today = 0;
    end

    % --- (c3) Live->dead coarse-root turnover [coarse-root 2026-06-12]: livecroot ages into deadcroot ---
    % Belowground mirror of (c2) stem, SAME lwtop. The C side (turnoverCarbonMirror c3) moves ALL
    % livecroot C -> deadcroot; here the N follows the SAME physics as stem: dead coarse root keeps
    % only CN_lw/CN_dw of the N (heartwood N-poor), the rest retransl -> ROOT reserve (belowground
    % source). Without this, livecroot C:N would drift down -> N_livecroot (hence R_m_livecroot)
    % inflates. CLM5 CNPhenologyMod.F90:2773 + N retransl (Veg_Phenology TechNote:327).
    if isfield(pftParams,'lwtop') && pftParams.lwtop > 0 && isfield(pools_out,'N_livecroot') && pools_out.N_livecroot > 0
        lwtop_frac_cr = pftParams.lwtop * dt / 31536000;          % [1/yr] annual frac -> per dt
        dN_lw_cr      = pools_out.N_livecroot * lwtop_frac_cr;    % [g N/m2] N leaving live coarse root
        dN_to_dead_cr = dN_lw_cr * (pftParams.CN_lw / pftParams.CN_dw);  % dead coarse root keeps fewer N
        dN_retrans_cr = dN_lw_cr - dN_to_dead_cr;                % rest recycled to ROOT reserve
        pools_out.N_livecroot    = pools_out.N_livecroot    - dN_lw_cr;
        pools_out.N_deadcroot    = pools_out.N_deadcroot    + dN_to_dead_cr;
        pools_out.N_root_storage = pools_out.N_root_storage + dN_retrans_cr;
    end

    % --- (d) Spring reserve release: REMOVED 2026-06-06 (unify-spring-release) ---
    % The old springStorageRelease (reserve -> N_labile during the budburst window) is
    % deleted: spring remobilisation is now done solely by springOnsetXfer (Step 8 in
    % runPlantCNStep), the CLM5-faithful storage->xfer->display path, for ALL storage PFTs
    % (both LAI engines). This removes the E2 double-release (two mechanisms draining the
    % same reserve in overlapping spring windows). diag fields stay 0 (no C-mirror release).
end

function [p, d, state] = winterWheatTurnoverN(p, d, state, prm, env)
% CoLM crop onset/offset adapter. Grain retained for the next 3 gC m-2 seed;
% Change: 2026-09-09 Jizixing Zhu -- Changwu NM/NPM management (Su Fuyuan
% thesis pp.15-16): remaining grain and all standing shoots leave the field;
% roots enter belowground litter. Harvest ledger includes straw, not grain only.
ps = env.phen_state;
% Daily phenology flags can persist over several hourly steps. Keep the
% applied harvest latch in the existing accepted turnover state, and rearm
% only in the next live pre-grain crop (not at the calendar-year boundary).
if (isfield(ps,'crop_planting_event') && ps.crop_planting_event) ...
        || (isfield(ps,'crop_live') && ps.crop_live ...
            && isfield(ps,'crop_phase') && ps.crop_phase < 3)
    state.harvested_this_year = 0;
end
already_harvested = isfield(state,'harvested_this_year') && state.harvested_this_year;
d.crop_emergence_event = isfield(ps, 'crop_emergence_event') && ps.crop_emergence_event;
d.crop_harvest_event = isfield(ps, 'crop_harvest_event') && ps.crop_harvest_event ...
    && ~already_harvested;
if isfield(ps, 'crop_emergence_event') && ps.crop_emergence_event
    p.N_leaf = p.N_leaf + p.N_fruit;
    p.N_fruit = 0;
end
if ~d.crop_harvest_event
    return;
end
assert(isfield(env,'C_fruit_after_allocation') ...
    && isscalar(env.C_fruit_after_allocation) ...
    && isfinite(env.C_fruit_after_allocation) && env.C_fruit_after_allocation >= 0, ...
    'turnoverAndHarvest:missingCropGrainCarbon', ...
    'Crop seed retention requires the same-step physical grain carbon.');
assert(isscalar(prm.CN_L) && isfinite(prm.CN_L) && prm.CN_L > 0 ...
    && isscalar(prm.initial_seed_c) && isfinite(prm.initial_seed_c) ...
    && prm.initial_seed_c >= 0, 'turnoverAndHarvest:invalidCropSeed', ...
    'Crop seed C amount and C:N must be finite and physically valid.');
d.crop_seed_carbon = min([prm.initial_seed_c, ...
    env.C_fruit_after_allocation, p.N_fruit * prm.CN_L]);
seedN = d.crop_seed_carbon / prm.CN_L;
% COMPASS management adaptation, NOT a claimed CoLM retransn reset:
% the crop's unallocated mobile N and aboveground reserve leave with the
% removed shoots; any belowground reserve joins root residue. Every amount
% is debited and entered in an existing physical export/litter ledger.
% Change: 2026-10-22 Wenjie Liu -- prescribed stubble management supersedes
% complete standing-shoot export; use the same fraction as the carbon mirror.
f_stubble = 0.15;
straw_N = p.N_leaf + p.N_livestem + p.N_deadstem + p.N_stem_storage;
stubble_N = f_stubble * straw_N;
d.harvest_today = p.N_fruit - seedN + (1 - f_stubble) * straw_N + p.N_labile;
d.leaf_litter_today = stubble_N;
d.root_litter_today = p.N_root_fine + p.N_livecroot + p.N_deadcroot + p.N_root_storage;
p.N_labile = 0;
p.N_root_storage = 0;
p.N_stem_storage = 0;
p.N_leaf = 0;
p.N_root_fine = 0;
p.N_livestem = 0;
p.N_deadstem = 0;
p.N_livecroot = 0;
p.N_deadcroot = 0;
p.N_fruit = seedN;
state.harvested_this_year = 1;
end

% Change: 2026-09-30 Wenjie Liu -- harvest explicit standing senescent
% leaf N through the existing standing-shoot export ledger, once per crop.
% CoLM crop onset/offset adapter. Grain retained for the next 3 gC m-2 seed;
% Change: 2026-09-30 Wenjie Liu -- withdraw this round's APSIM stem WSC
%   and standing-senescent-leaf extension; retain the original COMPASS litter
%   and CoLM crop C/N allocation/harvest contracts. Earlier Change records are historical.
