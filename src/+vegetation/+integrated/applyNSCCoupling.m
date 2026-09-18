function [C_pools_out, nsc_diag] = applyNSCCoupling( ...
    C_pools_prev, GPP_kg, R_m_total_kg, R_g_kg, NPP_kg, current_DOY, pftParams, dt, A_N_g, phen_alloc, C_xsmr_prev_g)
%APPLYNSCCOUPLING  NSC-limited respiration + C_storage overdraft (M2-B.8b 全年路径).
%
%   [v12 2026-05-29] 可选第 9 参 A_N_g (struct, 实际 N 分配 [g N/m²/dt], 字段
%     .leaf/.root/.stem/.fruit)。透传给 computeCarbonAllocation 触发 C-follows-N
%     (A_C=A_N×CN_target)。不传则 C 分配走旧式 NPP×f_alloc。
%
%   Wraps computeCarbonAllocation with two physically-required couplings that the
%   raw alloc cannot enforce on its own:
%
%     (a) Respiration cap: R_m_actual = min(R_a_potential, GPP×dt + NSC_avail).
%         当 NSC (= C_labile + C_storage) 不够维持理论呼吸时, 下调 R_m_actual.
%         物理对标 CLM5 curmr_ratio (NutrientCompetitionCLM45defaultMod.F90:709-720):
%           availc = gpp - mr; if mr>0 && availc<0: curmr_ratio = gpp/mr (GPP 按比例 cover 呼吸).
%         触发即记一次 storage_drawdown_event (= 本步动用了储备, 非致死).
%         [命名澄清 2026-05-29] 旧名 starvation_event 易误解为"碳饥饿致死"。本 flag 仅标记
%         "本步 R_m 需求超 GPP, 动用 NSC 储备"——是苹果树正常 NSC 周转 (白天存糖/夜间果实膨大
%         期取储备), 不触发任何死亡。
%         对比真碳饥饿致死 (FATES EDMortalityFunctionsMod.F90:201-228, D:/OneDrive/fates/):
%           cmort 基于 frac = store_c / target_leaf_c (储备碳/满冠目标叶碳), frac 持续低于
%           mort_upthresh_cstarvation 阈值才按线性/指数致死 — 看的是储备持续耗尽, 不是单步呼吸>GPP。
%         CLM5 大叶模型 (CNGapMortalityMod.F90) 本身不做碳饥饿致死 (用固定 gap mortality 速率)。
%         COMPASS 当前不接致死: storage_drawdown_event 仅诊断计数。
%
%     (b) C_storage overdraft: alloc 后若 C_labile < 0, 从 C_storage 借满补 0.
%         物理: NSC 池统一调用, 不分 "现金/存款" 等级. C_storage 全年波动 (秋充冬消春释),
%         不再只在春萌 (springStorageRelease.m) 时单向流出.
%
%   Returns:
%     C_pools_out  struct  6 updated C pools (C_labile ≥ 0 严格保证)
%     nsc_diag     struct  per-KT diagnostics:
%       .R_a_potential_g     [g C/m²/dt]  理论呼吸 (温度+N驱动, 上调前)
%       .R_a_actual_g        [g C/m²/dt]  实际呼吸 (NSC cap 后)
%       .R_a_actual_kg       [kg C/(m²·s)] 同上, kg/s 单位 (供 cum_R_a 累加)
%       .storage_drawdown_event [-]       1 = 本步 R_m>GPP 动用 NSC 储备, 0 = 未动用 (非致死)
%       .nsc_total_g         [g C/m²]     physical NSC = C_labile + C_root_storage + C_stem_storage (sugar+starch, entry value)
%       .nsc_net_g           [g C/m²]     net NSC = nsc_total_g + C_xsmr (debt-inclusive, <= nsc_total_g)
%       .c_storage_drawdown_g [g C/m²/dt] C_storage 本步净变化 (+ = 充值, − = 取款)
%
%   Inputs:
%     C_pools_prev   struct  6 C pools 上 KT 末值 [g C/m²]
%     GPP_kg         [kg C/(m²·s)]  fluxes.GPP
%     R_m_total_kg   [kg C/(m²·s)]  R_m_leaf+R_m_root+R_m_stem (computeNPP 输出, "potential")
%     R_g_kg         [kg C/(m²·s)]  R_g potential (computeNPP)
%     NPP_kg         [kg C/(m²·s)]  NPP potential (= max(GPP - R_m_pot - R_g_pot, 0))
%     current_DOY    [-]            fruit gate
%     pftParams      struct         f_alloc_*, has_fruit, CN_*, retransl_fraction,
%                                    storage_cushion [-] (C_storage 目标=cushion×C_leaf, FATES)
%     dt             [s]            outer time step (3600)
%
%   Method:
%     1. R_a_potential_kg  = R_m_total_kg + R_g_kg
%     2. NSC_avail_g       = C_labile_prev + C_storage_prev
%     3. cap_g             = GPP_kg × 1000 × dt + NSC_avail_g
%     4. R_a_actual_g      = min(R_a_potential_kg × 1000 × dt, cap_g)
%     5. R_a_actual_kg     = R_a_actual_g / (1000 × dt)   [kg/s, equivalent rate]
%     6. 调 computeCarbonAllocation 传 (GPP_kg, R_a_actual_kg, 0, NPP_kg_capped, ...)
%        (R_g 已并入 R_a_actual, 第二参数置零避免重复扣)
%     7. Overdraft: if C_pools.C_labile < 0 → 从 C_storage 借 min(|deficit|, C_storage)
%        剩余 deficit (storage 也耗尽) → C_labile clamp 0 + 守恒由 R_a cap 保证 (此时才接近真饥饿, 但不致死)
%
%   Conservation:
%     ΔC_labile + ΔC_storage + Σ ΔC_organ = GPP×dt − R_a_actual×dt
%     即 cap 后 GPP - R_a 严格等于 6 池 net change. 残差 → 机器精度.
%
%   Reference:
%     Sala A, Woodruff DR, Meinzer FC (2012) Carbon dynamics in trees: feast or famine?
%       Tree Physiol 32:764-775
%     Hartmann H, Trumbore S (2016) Understanding the roles of NSC in tree physiology.
%       New Phytol 211:386-403
%     DAISY (Hansen 1990) production.C — single NSC pool drives R_m
%     FATES EDMortalityFunctionsMod.F90:201-228 (D:/OneDrive/fates/) — 碳饥饿致死 cmort
%       基于 frac=store_c/target_leaf_c (储备/目标叶碳), 与本函数的 R_a cap 不同 (见上 (a) 命名澄清)
%     cn_coupling_path_and_literature_audit.md §12.26 (M2-B.8b NSC 联动设计)
%
%   Variable class (per compass-coding-standard §1):
%     C_pools_prev/out  — prognostic state (caller owns 6 池)
%     pftParams         — parameters (immutable)
%     GPP/R_m/R_g/NPP   — forcing slice this dt
%     nsc_diag          — diagnostics only
%
%   Author: Zhu Jizixing 2026-05-22 (M2-B.8b S1)
%   Editor: Jizixing Zhu; Date: 2026-08-24.
%   Change: 2026-08-24 Jizixing Zhu -- enforce one no-debt affordability
%           budget across maintenance and growth respiration.

    % --- Step 1: respiration fluxes this dt ---
    R_m_g       = max(R_m_total_kg, 0) * 1000 * dt;     % [g C/m2/dt] maintenance respiration
    R_g_g       = max(R_g_kg, 0)       * 1000 * dt;     % [g C/m2/dt] growth respiration
    R_g_g_pot   = R_g_g;                                  % preserve pre-affordability potential
    GPP_in_KT_g = max(GPP_kg, 0)       * 1000 * dt;     % [g C/m2/dt]
    if nargin < 11 || isempty(C_xsmr_prev_g), C_xsmr_prev_g = 0; end   % [g C/m2] xsmrpool prev (<=0 = debt)

    % --- Step 1b: maintenance-respiration substrate scheme (EXPERIMENT switch, 2026-06-12) ---
    % [resp-substrate EXP 2026-06-12] Maintenance respiration is C-substrate limited, NOT a fixed
    %   coefficient. Four candidate schemes are switchable via pftParams.resp_scheme so they can be
    %   compared in controlled runs (only the scheme differs). Once the winner is chosen, prune the
    %   losing branches and promote the choice to config (remove this EXPERIMENT block / numeric switch).
    %     0 = none  : beta_sub = 1 (baseline; fixed-coefficient R_m, unchanged).
    %     1 = bruhn : nocturnal decline. Bruhn et al. 2022 Nat Commun 13:5650 Eq 3-4
    %                 R/R0 = 1 - 0.08*h^0.54, h = hours since sunset (night only; day beta=1).
    %     2 = tc_mm : Michaelis-Menten substrate turn-down. Thornley & Cannell 2000 Ann Bot 85:55
    %                 Eq 17 beta = C_labile/(K_m + C_labile); -> 1 when ample (GMR limit), -> 0 as
    %                 C_labile depletes; prevents carbohydrate/xsmr from going chronically negative.
    %     3 = jones : fast-pool respiration. Jones et al. 2024 Commun Earth Environ 5:148 Eq 8
    %                 R = (1-Yg)*kf*S_f, with S_f = C_labile (slow pool = COMPASS storage remob).
    %                 REPLACES the temperature-based R_m with a substrate-proportional rate.
    %   Principle: Cannell & Thornley 2000 Ann Bot 85:45. Empirical: O'Leary et al. 2017 Plant Physiol
    %   174:2261 (nighttime leaf R ~ carbohydrate store, r2=0.73).
    %   K_m / kf / Yg: form from cited papers; absolute values site-tuned for Changwu apple.
    %   [SOURCE UNVERIFIED 2026-06-12] no published apple K_m located; tuned to CUE 0.4-0.6 + debt bound.
    resp_scheme = 0; resp_Km = 2; resp_kf = 0.227; resp_Yg = 0.75; night_h = 0;
    if isfield(pftParams,'resp_scheme')      && ~isempty(pftParams.resp_scheme),      resp_scheme = pftParams.resp_scheme; end
    if isfield(pftParams,'resp_Km')          && ~isempty(pftParams.resp_Km),          resp_Km     = pftParams.resp_Km; end
    if isfield(pftParams,'resp_jones_kf')    && ~isempty(pftParams.resp_jones_kf),    resp_kf     = pftParams.resp_jones_kf; end
    if isfield(pftParams,'resp_jones_Yg')    && ~isempty(pftParams.resp_jones_Yg),    resp_Yg     = pftParams.resp_jones_Yg; end
    if isfield(pftParams,'resp_night_hours') && ~isempty(pftParams.resp_night_hours), night_h     = pftParams.resp_night_hours; end
    Cl_prev   = max(C_pools_prev.C_labile, 0);                      % [g C/m2] free soluble-sugar pool (jones fast pool)
    % [resp substrate fix 2026-06-14] The Michaelis-Menten substrate turn-down (schemes 2/4/5) reads the
    %   TOTAL remobilisable NSC (soluble sugar C_labile + remobilisable starch C_root_storage +
    %   C_stem_storage), NOT the free-sugar pool alone. C_labile is swept toward storage every step
    %   (computeCarbonAllocation Step 3c), so keying the turn-down on C_labile fired the down-regulation
    %   even when starch reserves were abundant (spurious starvation signal). Starch remobilises to sugar
    %   on demand to fund maintenance respiration, so the respiratory substrate is the whole mobile pool.
    %   Martinez-Vilalta et al. 2016 Ecol Monogr 86:495 Fig.8 + Discussion (soluble-sugar depletion rare,
    %     starch depletion common; sugars kept above a critical threshold via starch<->sugar interconversion)
    %   O'Leary et al. 2017 Plant Physiol 174:2261 (night leaf R scales with total carbohydrate store)
    %   The jones scheme (3) is explicitly a sugar-only fast pool, so it keeps Cl_prev.
    C_resp_substrate = Cl_prev + max(C_pools_prev.C_root_storage, 0) + ...
                       max(C_pools_prev.C_stem_storage, 0);         % [g C/m2] total remobilisable NSC
    R_m_g_pot = R_m_g;                                              % [g C/m2/dt] potential (temp+N) R_m before scheme
    switch resp_scheme
        case 1   % bruhn nocturnal decline (night only)
            if night_h > 0
                R_m_g = R_m_g * max(1 - 0.08 * night_h^0.54, 0.30);
            end
        case {2, 4, 5}   % tc_mm (2), tc_mm_nodebt (4), tc_mm_fastrec (5): M-M turn-down on total NSC (Step 2 differs)
            if resp_Km > 0
                R_m_g = R_m_g * (C_resp_substrate / (resp_Km + C_resp_substrate));
            end
        case 3   % jones fast-pool: R = (1-Yg)*kf*S_f  (S_f = C_labile, sugar-only fast pool)
            R_m_g = (1 - resp_Yg) * resp_kf * Cl_prev * (dt / 3600);
        otherwise   % 0 = none (baseline)
            % R_m_g unchanged
    end
    beta_sub = R_m_g / max(R_m_g_pot, eps);                         % [-] effective turn-down (diagnostic)

    % --- Step 2: maintenance-respiration funding (xsmr debt pool vs nodebt) ---
    if resp_scheme == 4
        % [resp tc_mm_nodebt EXP 2026-06-12] Fund the (beta-throttled) R_m DIRECTLY from C_labile,
        %   with NO xsmrpool debt. Cap R_m at what GPP + C_labile can fund this step, so the plant
        %   respires only what it can afford (graded down-regulation), never carrying a deficit. The
        %   Step 1b beta turn-down already drives R_m -> 0 as TOTAL NSC (sugar+starch) depletes; this branch lets the
        %   labile sugar (which the GPP-only xsmr funding leaves unused) pay night respiration,
        %   eliminating the summer debt blow-up seen with the xsmrpool. Thornley & Cannell 2000 Ann
        %   Bot 85:55 Eq 17 (substrate turn-down prevents carbohydrate < 0 -> no debt pool needed).
        % [nodebt+storage] Cap total R_m+R_g at GPP + C_labile + C_storage (sugar + remobilisable
        %   starch). Restores the pre-xsmr resp-fix-A NSC_avail = C_labile + C_storage cap (2026-06-07):
        %   when sugar is short the plant remobilises its starch reserve to fund maintenance, rather
        %   than carrying debt. The Step 1b beta turn-down (on total NSC = C_labile+C_storage) still grades R_m down.
        nsc_avail_g = C_pools_prev.C_labile + C_pools_prev.C_root_storage + C_pools_prev.C_stem_storage;
        total_affordable_g = max(GPP_in_KT_g + nsc_avail_g, 0);
        R_m_g       = min(R_m_g, total_affordable_g);                  % maintenance has first claim
        R_g_g       = min(R_g_g, max(total_affordable_g-R_m_g,0));    % growth respiration shares the same C budget
        R_m_cur_g   = R_m_g;                                         % full (throttled+capped) R_m paid from cpool
        R_m_xs_g    = 0;                                             % no xsmr debt
        availc_g    = max(GPP_in_KT_g - R_m_g, 0);
        recover_g   = 0;
    else
        % --- CLM5 cpool + xsmrpool (NutrientCompetitionCLM45defaultMod.F90:704-744) ---
        % [xsmrpool 2026-06-11] Maintenance respiration is paid in FULL. GPP funds fraction curmr_ratio
        %   (drawn from cpool = C_labile); the shortfall R_m_xs is carried as xsmrpool (C_xsmr) DEBT
        %   (goes negative) and repaid from future GPP surplus, with priority BEFORE growth, over
        %   dayscrecover days. REPLACES the old "cap R_a at GPP+C_labile" (resp-fix-A 2026-06-07), which
        %   under-paid respiration under carbon stress. cpool (C_labile) and xsmrpool (C_xsmr) are SEPARATE
        %   pools, ported one-to-one from CLM5; CO2 to atmosphere = full R_a regardless of funding pool.
        %   State updates mirror CNCStateUpdate1Mod.F90:302-327.
        dayscrecover = 30;                                              % [days] CLM5 default
        if isfield(pftParams,'xsmr_recovery_days') && ~isempty(pftParams.xsmr_recovery_days)
            dayscrecover = pftParams.xsmr_recovery_days;
        end
        availc_g = GPP_in_KT_g - R_m_g;                                 % CLM5:704
        if R_m_g > 0 && availc_g < 0
            curmr_ratio = GPP_in_KT_g / R_m_g;                          % CLM5:711  GPP-funded fraction of MR
        else
            curmr_ratio = 1;                                            % CLM5:713
        end
        R_m_cur_g = R_m_g * curmr_ratio;                                % CLM5:715  MR paid now (from cpool)
        R_m_xs_g  = R_m_g - R_m_cur_g;                                  % CLM5:716  MR shortfall -> xsmr debt
        availc_g  = max(availc_g, 0);                                   % CLM5:727
        recover_g = 0;                                                  % [g C/m2/dt] cpool -> xsmr deficit recovery
        if C_xsmr_prev_g < 0                                            % CLM5:730  running an xsmr deficit
            if resp_scheme == 5
                % [tc_mm_fastrec EXP 2026-06-12] Repay the debt from ALL leftover GPP this step (capped
                %   at availc below), instead of CLM5's 30-day rate-limit. CLM5 dayscrecover=30 is sized
                %   for SEASONAL recovery, but the deficit here is driven by the DAILY night/day cycle;
                %   throttling daily repayment to 1/30-of-debt lets the nightly deficit outpace it ->
                %   debt creeps to -90. Daytime leftover GPP (verified ~1.9 vs 1.86/night) is enough to
                %   clear it if not throttled. Keeps the CLM5 xsmrpool structure; only the recovery rate.
                recover_g = -C_xsmr_prev_g;                            % full repayment (capped at availc)
            else
                recover_g = -C_xsmr_prev_g / (dayscrecover * 86400) * dt;   % CLM5:735  repay debt over tau (30 d)
            end
            if recover_g < availc_g                                     % CLM5:736
                availc_g = availc_g - recover_g;                        % CLM5:738  recovery taken before growth
            else
                recover_g = availc_g;                                   % CLM5:741
                availc_g  = 0;                                          % CLM5:742
            end
        end
    end

    % --- Step 3: growth budget (growth respiration stays with growth) ---
    NPP_actual_g  = max(availc_g - R_g_g, 0);                       % [g C/m2/dt] -> display growth
    NPP_actual_kg = NPP_actual_g / (1000 * dt);
    % cpool funds R_m_cur + R_g; xsmr funds R_m_xs. Pass cpool's respiration share to allocation so
    %   computeCarbonAllocation debits C_labile by R_m_cur+R_g (NOT the full R_m).
    R_a_cpool_kg    = (R_m_cur_g + R_g_g) / (1000 * dt);            % [kg C/(m2 s)] cpool-funded respiration
    R_a_actual_g    = R_m_g + R_g_g;                               % [g C/m2/dt] ACTUAL respiration after substrate turn-down (real CO2 flux; funds pools)
    R_a_actual_kg   = R_a_actual_g / (1000 * dt);
    % [Ra potential/actual split 2026-06-15] r_a_potential must hold the TRUE pre-turn-down value:
    %   R_m_g_pot (= temp+N maintenance R_m BEFORE the resp_scheme M-M turn-down, captured line 128) + R_g_g.
    %   Was 'R_a_potential_g = R_a_actual_g', which silently equated the two so r_a_potential /
    %   cum_r_a_potential lost the un-throttled value. R_g_g is not throttled, so it is common to both.
    R_a_potential_g = R_m_g_pot + R_g_g_pot;                      % [g C/m2/dt] potential (before substrate/affordability limits)
    storage_drawdown_event = double(R_m_xs_g > 1e-9);             % [-] 1 = ran on xsmr debt this step (diag)

    resp_only = nargin >= 9 && isstruct(A_N_g) && isfield(A_N_g, 'mode') && strcmp(A_N_g.mode, 'resp_only');
    if resp_only
        C_pools_out = C_pools_prev;
        C_pools_out.C_xsmr = C_xsmr_prev_g;
        nsc_diag.R_a_potential_g        = R_a_potential_g;
        nsc_diag.resp_theta_actual      = max(0, R_m_g / max(R_m_g_pot, eps));
        nsc_diag.R_a_actual_g           = R_a_actual_g;
        nsc_diag.R_a_actual_kg          = R_a_actual_kg;
        nsc_diag.NPP_actual_kg          = NPP_actual_kg;
        nsc_diag.storage_drawdown_event = storage_drawdown_event;
        nsc_diag.nsc_total_g            = C_pools_prev.C_labile + C_pools_prev.C_root_storage + ...
                                          C_pools_prev.C_stem_storage;
        nsc_diag.nsc_net_g              = nsc_diag.nsc_total_g + C_xsmr_prev_g;
        nsc_diag.c_storage_drawdown_g   = 0;
        nsc_diag.C_xsmr                 = C_xsmr_prev_g;
        nsc_diag.R_m_xs_g               = R_m_xs_g;
        nsc_diag.xsmr_recover_g         = recover_g;
        nsc_diag.resp_beta_sub          = beta_sub;
        nsc_diag.resp_substrate_g       = C_resp_substrate;
        return
    end

    % --- Step 4: allocation (cpool side); reserve the xsmr recovery from cpool BEFORE alloc+storage ---
    % [xsmrpool] recovery (cpool->xsmr) has priority over growth AND over the storage refill. Removing
    %   it from C_labile before computeCarbonAllocation stops the storage surplus from double-claiming
    %   it (which would push C_labile negative). Conservation verified:
    %   ΔC_labile + ΔC_xsmr + Σorgan + Δstorage = GPP - R_m - R_g.
    if nargin < 10, phen_alloc = []; end
    c_storage_in = C_pools_prev.C_root_storage + C_pools_prev.C_stem_storage;   % [g C/m2] for c_storage_drawdown diag
    C_pools_for_alloc = C_pools_prev;
    C_pools_for_alloc.C_labile = C_pools_prev.C_labile - recover_g;              % cpool hands recovery to xsmrpool
    % [labile sugar floor 2026-06-14] Soluble-sugar working-pool floor protected during allocation
    %   (computeCarbonAllocation Step 3c). Sized from POTENTIAL maintenance respiration R_m_total_kg
    %   (temperature/N-driven, BEFORE the substrate turn-down above) so a low substrate signal cannot
    %   shrink its own target. See docs/nsc_labile_sugar_floor_design.md.
    C_labile_target_g = vegetation.carbon.computeLabileSugarTarget(R_m_total_kg, pftParams, dt);   % [g C/m2]
    if nargin >= 9 && ~isempty(A_N_g)
        [C_pools_out, c_alloc_diag] = vegetation.carbon.computeCarbonAllocation( ...
            C_pools_for_alloc, GPP_kg, R_a_cpool_kg, 0, NPP_actual_kg, ...
            current_DOY, pftParams, dt, A_N_g, phen_alloc, C_labile_target_g);
    else
        [C_pools_out, c_alloc_diag] = vegetation.carbon.computeCarbonAllocation( ...
            C_pools_for_alloc, GPP_kg, R_a_cpool_kg, 0, NPP_actual_kg, ...
            current_DOY, pftParams, dt, [], phen_alloc, C_labile_target_g);
    end

    % --- Step 5: C_labile non-negativity + xsmrpool update (CNCStateUpdate1Mod.F90:302-327) ---
    % [conservation fix 2026-06-12] If cpool dips below 0 (the recovery debit + storage refill drew
    %   more than was available this step), do NOT clamp to 0 — that would CREATE carbon (the v1 run
    %   leaked 0.38 gC/yr this way). Instead carry the shortfall as additional xsmr debt: it is
    %   deferred carbon the cpool could not fund, exactly what xsmrpool represents. Conserves
    %   (C_labile += deficit, C_xsmr -= deficit; net 0).
    cpool_deficit_g    = max(-C_pools_out.C_labile, 0);            % [g C/m2] cpool shortfall this step
    C_pools_out.C_labile = C_pools_out.C_labile + cpool_deficit_g; % -> 0 (no carbon created)
    if resp_scheme == 4
        % [nodebt 2026-06-12] Cover the C_labile deficit from C_storage (starch remobilisation), root
        %   reserve first then stem. The shared total-respiration cap above guarantees no residual
        %   demand remains after both reserves are exhausted; any residual is an internal error.
        draw_r = min(cpool_deficit_g, max(C_pools_out.C_root_storage, 0));
        C_pools_out.C_root_storage = C_pools_out.C_root_storage - draw_r;
        draw_s = min(cpool_deficit_g - draw_r, max(C_pools_out.C_stem_storage, 0));
        C_pools_out.C_stem_storage = C_pools_out.C_stem_storage - draw_s;
        residual_g   = cpool_deficit_g - draw_r - draw_s;
        assert(residual_g<=1e-9, ...
            'applyNSCCoupling:nodebtAffordabilityFailure', ...
            ['tc_mm_nodebt produced an unfunded carbon demand of %g g C m^-2 ', ...
             'after the total-respiration affordability cap.'],residual_g);
        C_xsmr_new_g = C_xsmr_prev_g;                              % nodebt never creates new debt
    else
        C_xsmr_new_g = C_xsmr_prev_g - R_m_xs_g + recover_g - cpool_deficit_g;   % [xsmrpool] MR shortfall + cpool deficit as debt; recovery repays
    end
    C_pools_out.C_xsmr = C_xsmr_new_g;                             % [g C/m2] xsmrpool (<=0 = debt)

    % --- Step 6: Diagnostics ---
    nsc_diag.R_a_potential_g     = R_a_potential_g;
    % [actual organ resp 2026-06-15] Effective maintenance throttle INCLUDING the nodebt afford-cap:
    %   R_m_g is the FINAL maintenance respiration (after substrate turn-down + cap); R_m_g_pot is the
    %   pre-turn-down temp+N R_m (line 128). theta = R_m_g/R_m_g_pot is the true actual/potential ratio,
    %   stricter than resp_beta_sub (which is computed BEFORE the cap). Used by runPlantCNStep to split
    %   actual Ra into organ components (actual organ Rm = theta * potential organ Rm).
    % Unclamped above 1 ON PURPOSE: the organ split (runPlantCNStep) needs sum(theta*organ_pot) = R_m_g to
    %   close. For turn-down schemes theta in [0,1]; a substitution scheme (jones) can raise actual R_m above
    %   potential -> theta>1, and clamping to 1 would make the organ-actual sum < R_m_g. Only the lower bound
    %   (theta>=0) and the /0 guard are kept.
    nsc_diag.resp_theta_actual   = max(0, R_m_g / max(R_m_g_pot, eps));   % [-] actual R_m / potential R_m (1=none, <1 turn-down, >1 under substitution e.g. jones)
    nsc_diag.R_a_actual_g        = R_a_actual_g;
    nsc_diag.R_a_actual_kg       = R_a_actual_kg;
    nsc_diag.NPP_actual_kg       = NPP_actual_kg;
    nsc_diag.storage_drawdown_event = storage_drawdown_event;
    % [NSC naming split 2026-06-14] Separate the PHYSICAL carbohydrate pool from the debt-inclusive
    %   net. nsc_total_g = real remobilisable carbohydrate (soluble sugar C_labile + starch
    %   C_root_storage + C_stem_storage) — the SAME pool the M-M maintenance-resp turn-down reads as
    %   substrate. nsc_net_g additionally folds in the xsmrpool maintenance-respiration debt
    %   (C_xsmr <= 0), so nsc_net_g <= nsc_total_g. Keeping them apart stops the debt from
    %   masquerading as available carbohydrate when reading "NSC total".
    %   (xsmrpool: CLM5 CNVegCarbonStateType.F90 xsmrpool; starch as buffered mobilisable reserve:
    %    Martinez-Vilalta et al. 2016 Ecol Monogr 86:495.)
    nsc_diag.nsc_total_g         = C_pools_prev.C_labile + C_pools_prev.C_root_storage + ...
                                   C_pools_prev.C_stem_storage;       % [g C/m2] physical NSC (sugar + starch)
    nsc_diag.nsc_net_g           = nsc_diag.nsc_total_g + C_xsmr_prev_g;  % [g C/m2] net NSC incl. xsmr debt (<= nsc_total_g)
    nsc_diag.c_storage_drawdown_g = (C_pools_out.C_root_storage + C_pools_out.C_stem_storage) - c_storage_in;
    % [xsmrpool 2026-06-11] xsmr pool + fluxes for write-back (runPlantCNStep) and diagnostics
    nsc_diag.C_xsmr        = C_xsmr_new_g;          % [g C/m2] xsmrpool after this step (<=0 = debt)
    nsc_diag.R_m_xs_g      = R_m_xs_g;              % [g C/m2/dt] maintenance respiration funded by xsmr debt
    nsc_diag.xsmr_recover_g = recover_g;            % [g C/m2/dt] debt repaid from GPP this step
    nsc_diag.resp_beta_sub  = beta_sub;             % [-] [resp-substrate EXP] effective R_m turn-down (1=none)
    nsc_diag.resp_substrate_g  = C_resp_substrate;  % [g C/m2] [resp substrate fix] total remobilisable NSC seen by M-M turn-down
    nsc_diag.C_labile_target_g = C_labile_target_g; % [g C/m2] [labile floor] soluble-sugar floor protected this step
    nsc_diag.C_storage_to_labile = c_alloc_diag.C_storage_to_labile;   % [g C/m2/dt] [labile floor] starch->sugar remobilisation this step
    % [v8] C 各器官分配诊断 (computeCarbonAllocation diag, g C/m²/dt; NPP×f_alloc, fruit gate)
    nsc_diag.A_C_leaf  = c_alloc_diag.A_C_leaf;
    nsc_diag.A_C_root  = c_alloc_diag.A_C_root;
    nsc_diag.A_C_stem  = c_alloc_diag.A_C_stem;
    nsc_diag.A_C_fruit = c_alloc_diag.A_C_fruit;
    % [mass-balance N reserve 2026-06-11] C-side reserve-refill diagnostics only. The N reserve is no
    %   longer derived from these (former N = C_to_storage / CN_storage_target removed; the N reserve
    %   is mass-balanced from resorption in runPlantCNStep Step 6). Kept as C-flux diagnostics.
    nsc_diag.C_to_root_storage = c_alloc_diag.C_to_root_storage;   % [g C/m2/dt]
    nsc_diag.C_to_stem_storage = c_alloc_diag.C_to_stem_storage;   % [g C/m2/dt]
end
