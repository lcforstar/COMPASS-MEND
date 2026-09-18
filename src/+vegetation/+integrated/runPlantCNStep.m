function [NitrogenVariables, R_a_plant_kg] = runPlantCNStep( ...
    NitrogenVariables, demand_ctx, U_N_total_g_m2_s, env_turn, Delt_t, KT)
%RUNPLANTCNSTEP  Per-KT plant C/N dispatcher (M2-B.1b + M2-B.4 + M2-B.8 + M2-B.8b).
%
%   主调度器: 把"上 KT pools → NPP/R_a → N alloc → C alloc (NSC cap) → N turnover →
%   C turnover mirror → 写回"整套链路集中到一个函数, 让 COMPASS.m 主循环只剩 1 行.
%
%   Returns:
%     NitrogenVariables  struct  完整更新后的 plant C/N state (6 N 池 + 6 C 池 + cum_* + diag)
%     R_a_plant_kg       [kg C/(m²·s)]  当 KT 整株自养呼吸 potential (cap 前; 仅作 r_a_potential 诊断. NEE/Reco 改用 nsc_diag.R_a_actual_kg)
%
%   Inputs:
%     NitrogenVariables  struct  在调用前的 plant C/N state (KT-1 末值 + 累积器)
%     demand_ctx         struct  来自 runPlantNDemand 输出, 字段:
%                          .pftType            [-]            PFT 标签
%                          .pftP_demand        struct         demand-side PFT params (r_m10, R_g_fraction, f_alloc_*, CN_*)
%                          .GPP_for_demand     [kg C/(m²·s)]  fluxes.GPP (KT)
%                          .Rd_leaf_for_demand [kg C/(m²·s)]
%                          .T_air_demand       [°C]
%                          .T_soil_demand      [°C]
%                          .D_N_for_demand     struct         organ demand (复用)
%                          .DOY_for_demand     [-]            day-of-year
%     U_N_total_g_m2_s   [g N/(m²·s)]  来自 applyRootUptake 的根吸总量
%     env_turn           struct  turnover env: .frPHU .LAI .DLAI .year
%     Delt_t             [s]
%     KT                 [-]
%
%   Method (内部串调 6 步, 详见 cn_coupling §12.26):
%     1. 装填 N pools_prev (KT-1 末值)
%     2. computeNPP → NPP_kg, R_a_plant_kg, npp_diag (R_m_leaf/root/stem + R_g)
%     3. allocatePlantN → pools_new (5 N pools updated)
%     4. 装填 C pools_prev + pftP_carbon (mirror demand 的 f_alloc/CN_target)
%     5. applyNSCCoupling → C_pools_after_alloc, nsc_diag (NSC cap + overdraft)
%     6. turnoverAndHarvest → pools_turn_out (6 N pools after turnover), turn_diag
%     7. turnoverCarbonMirror → C_pools_final (6 C pools after mirror), C_turn_diag
%     8. 写回 NitrogenVariables + 累积 cum_* + plantN_diag/nsc_diag/turnover_diag
%
%   Conservation:
%     Plant N: ΔΣpoolN + cum_lost = ΣU_N×dt + init (机器精度)
%     Plant C: ΔΣpoolC + cum_lost = ΣGPP×dt - ΣR_a_actual×dt (机器精度, NSC cap 后)
%     全年 smoke 验证残差 < 1e-10 g/m².
%
%   Reference:
%     cn_coupling §12.26 (M2-B.8b NSC 联动 + 二期重构)
%     +vegetation/+nitrogen/computeNPP.m, allocatePlantN.m
%     +vegetation/+integrated/applyNSCCoupling.m, turnoverAndHarvest.m, turnoverCarbonMirror.m
%
%   Variable class (per compass-coding-standard §1):
%     NitrogenVariables  — prognostic state (inout, caller owns)
%     demand_ctx         — forcing slice (in) 上游 runPlantNDemand 算出
%     U_N_total_g_m2_s   — forcing (in) 上游 applyRootUptake 算出
%     env_turn           — forcing (in)
%     R_a_plant_kg       — diagnostic (out) 供 NEE 计算
%
%   Author: Zhu Jizixing 2026-05-22 (M2-B.8c 二期重构 S2)
%   Editor: Jizixing Zhu; Date: 2026-08-24.
% Change: 2026-06-28 Jizixing Zhu — spring onset-xfer reserve split routed by PFT class
%   (fruit tree -> FruitTree GDD phase; natural woody -> c_allometry; crop/legacy -> f_alloc).
%   Previously always f_alloc_*, which mis-split Apple's spring remobilisation.
% Change: 2026-08-24 Jizixing Zhu -- diagnose realized growth respiration
%   from the shared no-debt respiration budget.

    % ------------------------------------------------------------
    % Step 1: N pools_prev (KT-1 末值, KT=1 时用 init)
    % ------------------------------------------------------------
    idx_prev = max(1, KT - 1);
    pools_prev.N_labile     = NitrogenVariables.N_labile(idx_prev);
    pools_prev.N_leaf       = NitrogenVariables.N_leaf(idx_prev);
    pools_prev.N_root_fine  = NitrogenVariables.N_root_fine(idx_prev);
    pools_prev.N_livestem   = NitrogenVariables.N_livestem(idx_prev);   % [13c]
    pools_prev.N_deadstem   = NitrogenVariables.N_deadstem(idx_prev);   % [13c]
    pools_prev.N_livecroot  = NitrogenVariables.N_livecroot(idx_prev);  % [coarse-root 2026-06-12]
    pools_prev.N_deadcroot  = NitrogenVariables.N_deadcroot(idx_prev);  % [coarse-root 2026-06-12]
    pools_prev.N_fruit      = NitrogenVariables.N_fruit(idx_prev);

    % ------------------------------------------------------------
    % Step 2: NPP / R_a 完整分解 (复用 demand 块的 GPP/Rd/T, 用 prev pools)
    % ------------------------------------------------------------
    N_fruit_for_resp = pools_prev.N_fruit;
    if isfield(demand_ctx, 'phen_alloc') && isfield(demand_ctx.phen_alloc, 'is_crop') ...
            && demand_ctx.phen_alloc.is_crop == 1 && ~demand_ctx.phen_alloc.crop_live
        N_fruit_for_resp = 0;
    end
    [NPP_kg, R_a_plant_kg, npp_diag] = vegetation.nitrogen.computeNPP( ...
        demand_ctx.GPP_for_demand, demand_ctx.Rd_leaf_for_demand, ...
        pools_prev.N_root_fine, pools_prev.N_livestem, pools_prev.N_livecroot, N_fruit_for_resp, ...
        demand_ctx.T_air_demand, demand_ctx.T_soil_demand, demand_ctx.pftP_demand);

    % ------------------------------------------------------------
    % Step 3: N alloc (5 organ pool)
    % ------------------------------------------------------------
    [pools_new, alloc_diag] = vegetation.nitrogen.allocatePlantN( ...
        pools_prev, U_N_total_g_m2_s, demand_ctx.D_N_for_demand, Delt_t, demand_ctx.pftP_demand);

    NitrogenVariables.N_labile(KT)     = pools_new.N_labile;
    NitrogenVariables.N_leaf(KT)       = pools_new.N_leaf;
    NitrogenVariables.N_root_fine(KT)  = pools_new.N_root_fine;
    NitrogenVariables.N_livestem(KT)   = pools_new.N_livestem;   % [13c]
    NitrogenVariables.N_deadstem(KT)   = pools_new.N_deadstem;   % [13c]
    NitrogenVariables.N_livecroot(KT)  = pools_new.N_livecroot;  % [coarse-root 2026-06-12]
    NitrogenVariables.N_deadcroot(KT)  = pools_new.N_deadcroot;  % [coarse-root 2026-06-12]
    NitrogenVariables.N_fruit(KT)      = pools_new.N_fruit;

    NitrogenVariables.plantN_diag(KT).NPP         = NPP_kg;
    NitrogenVariables.plantN_diag(KT).R_a_plant   = R_a_plant_kg;
    NitrogenVariables.plantN_diag(KT).R_m_leaf    = npp_diag.R_m_leaf;
    NitrogenVariables.plantN_diag(KT).R_m_root    = npp_diag.R_m_root;
    NitrogenVariables.plantN_diag(KT).R_m_stem    = npp_diag.R_m_stem;
    NitrogenVariables.plantN_diag(KT).R_m_livecroot = npp_diag.R_m_livecroot;   % [coarse-root 2026-06-12]
    NitrogenVariables.plantN_diag(KT).R_m_fruit   = npp_diag.R_m_fruit;   % [FruitTree] perennial fruit maintenance resp
    NitrogenVariables.plantN_diag(KT).R_g         = npp_diag.R_g;
    NitrogenVariables.plantN_diag(KT).NPP_net     = npp_diag.NPP_net;   % unclamped, for true CUE
    NitrogenVariables.plantN_diag(KT).f_N_growth  = alloc_diag.f_N_growth;
    NitrogenVariables.plantN_diag(KT).A_leaf      = alloc_diag.A_leaf;
    NitrogenVariables.plantN_diag(KT).A_root      = alloc_diag.A_root;
    NitrogenVariables.plantN_diag(KT).A_stem      = alloc_diag.A_stem;
    NitrogenVariables.plantN_diag(KT).A_croot     = alloc_diag.A_croot;   % [coarse-root 2026-06-12]
    NitrogenVariables.plantN_diag(KT).A_fruit     = alloc_diag.A_fruit;

    % ------------------------------------------------------------
    % Step 4: C pools_prev + pftP_carbon (mirror demand 的 alloc/CN)
    % ------------------------------------------------------------
    % [resp-substrate EXP 2026-06-12] track consecutive dark hours since sunset (bruhn scheme h).
    %   night = no gross photosynthesis this step (GPP ~ 0); resets to 0 each daylight step.
    if demand_ctx.GPP_for_demand <= 1e-9
        NitrogenVariables.resp_night_hours = NitrogenVariables.resp_night_hours + 1;
    else
        NitrogenVariables.resp_night_hours = 0;
    end
    C_pools_prev = struct( ...
        'C_labile',     NitrogenVariables.C_labile(idx_prev), ...
        'C_leaf',       NitrogenVariables.C_leaf(idx_prev), ...
        'C_root_fine',  NitrogenVariables.C_root_fine(idx_prev), ...
        'C_livestem',   NitrogenVariables.C_livestem(idx_prev), ...   % [13c]
        'C_deadstem',   NitrogenVariables.C_deadstem(idx_prev), ...   % [13c]
        'C_livecroot',  NitrogenVariables.C_livecroot(idx_prev), ...  % [coarse-root 2026-06-12]
        'C_deadcroot',  NitrogenVariables.C_deadcroot(idx_prev), ...  % [coarse-root 2026-06-12]
        'C_fruit',      NitrogenVariables.C_fruit(idx_prev), ...
        'C_root_storage', NitrogenVariables.C_root_storage(idx_prev), ...   % [dual-storage]
        'C_stem_storage', NitrogenVariables.C_stem_storage(idx_prev));      % [dual-storage]
    pftP_carbon = struct( ...
        'f_alloc_leaf',     demand_ctx.pftP_demand.f_alloc_leaf, ...
        'f_alloc_root',     demand_ctx.pftP_demand.f_alloc_root, ...
        'f_alloc_stem',     demand_ctx.pftP_demand.f_alloc_stem, ...
        'f_alloc_fruit',    demand_ctx.pftP_demand.f_alloc_fruit, ...
        'use_fruittree_alloc', demand_ctx.pftP_demand.use_fruittree_alloc, ... % [F2/13a] forward apple GDD-alloc gate to C side
        'ft_GDDlfmat',      demand_ctx.pftP_demand.ft_GDDlfmat, ...  % [F2/13a] leaf-taper GDD endpoint
        'ft_bfact',         demand_ctx.pftP_demand.ft_bfact, ...     % [F2/13a] leaf-taper shape factor
        'harvest_DS',       demand_ctx.pftP_demand.harvest_DS, ...   % [#1] fruit-gate upper bound (DS_fruit)
        'fleafi',           demand_ctx.pftP_demand.fleafi, ...       % [13a] 4-phase allocation params (forward to C side)
        'arooti',           demand_ctx.pftP_demand.arooti, ...
        'arootf',           demand_ctx.pftP_demand.arootf, ...
        'arootf2',          demand_ctx.pftP_demand.arootf2, ...
        'astemf',           demand_ctx.pftP_demand.astemf, ...
        'aleafstor',        demand_ctx.pftP_demand.aleafstor, ...
        'grnfill',          demand_ctx.pftP_demand.grnfill, ...
        'grnrp',            demand_ctx.pftP_demand.grnrp, ...   % [5-phase 2026-06-13] huiripe offset -> fruitTreeAllocation P2/P3 boundary
        'hybgdd',           demand_ctx.pftP_demand.hybgdd, ...
        'has_fruit',        demand_ctx.pftP_demand.has_fruit, ...
        'CN_L',             demand_ctx.pftP_demand.CN_L, ...
        'CN_target_root',   demand_ctx.pftP_demand.CN_target_root, ...
        'CN_target_stem',   demand_ctx.pftP_demand.CN_target_stem, ...
        'CN_target_fruit',  demand_ctx.pftP_demand.CN_target_fruit, ...
        'CN_lw',            vegetation.getPFTValue(demand_ctx.pftType, 'CN_lw'), ...   % [13c per-pool clamp 2026-06-12] live wood C:N
        'CN_dw',            vegetation.getPFTValue(demand_ctx.pftType, 'CN_dw'), ...   % [13c per-pool clamp] dead wood C:N
        'croot_stem_alloc', vegetation.getPFTValue(demand_ctx.pftType, 'croot_stem_alloc', 0), ... % [coarse-root] NPP coarse:stem ratio
        'froot_leaf',       vegetation.getPFTValue(demand_ctx.pftType, 'froot_leaf', 0), ...       % [fine-root] fine-root:leaf allometry
        'stem_leaf',        vegetation.getPFTValue(demand_ctx.pftType, 'stem_leaf', 0), ...        % [CLM5 all-PFT 2026-06-13] f3 stem:leaf (natural woody; apple uses GDD-phase astem instead)
        'retransl_fraction', vegetation.getPFTValue(demand_ctx.pftType, 'retransl_fraction'), ...
        'storage_cushion',   demand_ctx.pftP_demand.storage_cushion, ...
        'cn_flex_offset',    demand_ctx.pftP_demand.cn_flex_offset, ...   % [gC/gN] [FlexibleCN] 池 C:N 夹紧半幅
        'fcur',              demand_ctx.pftP_demand.fcur, ...             % [-] [A] 当年显示比例, (1−fcur)→storage
        'ndays_on',          vegetation.getPFTValue(demand_ctx.pftType, 'ndays_on', 30), ...   % [E] onset-xfer duration (CLM5 30 d, general)
        'fstor2tran',        vegetation.getPFTValue(demand_ctx.pftType, 'fstor2tran', 0), ...   % [E] storage->leaf onset frac (0=off unless set in PFT table; deciduous storage PFTs=0.5)
        'has_storage_pool',  vegetation.getPFTValue(demand_ctx.pftType, 'has_storage_pool', 0), ...   % [E] onset-xfer needs a reserve pool
        'flivewd',           demand_ctx.pftP_demand.flivewd, ...          % [13c] new stem C fraction to live wood (computeCarbonAllocation)
        'croot_stem',        vegetation.getPFTValue(demand_ctx.pftType, 'croot_stem', 0), ...   % [coarse-root-refill] coarse-root:stem ratio -> belowground reserve refill in computeCarbonAllocation (apple 0.225 Zanotelli 2013; 0=non-woody). MISSING before 2026-06-10 -> croot_stem silently 0 in v63/v64/v65.
        'nsc_frac',          vegetation.getPFTValue(demand_ctx.pftType, 'nsc_frac', 0.07), ...   % [target-reserve 2026-06-10] woody NSC C-fraction; reserve target = nsc_frac x own wood (computeReserveTarget). Apple 0.131; forwarded here to avoid the croot_stem-trap (3-segment chain).
        'lwtop',             vegetation.getPFTValue(demand_ctx.pftType, 'lwtop'), ...   % [13c] live->dead stem C turnover (mirror)
        'SLA0',              demand_ctx.pftP_demand.SLA0, ...              % [peaklai] LAI=f(C_leaf) cap check (computeCarbonAllocation)
        'SLAm',              demand_ctx.pftP_demand.SLAm, ...              % [peaklai] SLA canopy gradient
        'laimx',             demand_ctx.pftP_demand.laimx, ...            % [peaklai] physiological peak LAI -> leaf alloc cap
        'harvest_fraction',  vegetation.getPFTValue(demand_ctx.pftType, 'harvest_fraction'), ...   % [-] [v11] 采收 C 直接按比例扣 C_fruit (不用 CN 反推)
        'winter_wheat_cn',   vegetation.getPFTValue(demand_ctx.pftType, 'winter_wheat_cn', 0), ...
        'initial_seed_c',    vegetation.getPFTValue(demand_ctx.pftType, 'initial_seed_c', 3), ...
        'resp_scheme',       vegetation.getPFTValue(demand_ctx.pftType, 'resp_scheme', 0), ...      % [resp-substrate EXP 2026-06-12] 0=none 1=bruhn 2=tc_mm 3=jones
        'resp_Km',           vegetation.getPFTValue(demand_ctx.pftType, 'resp_Km', 2), ...          % [tc_mm] M-M half-saturation [g C/m2]
        'resp_jones_kf',     vegetation.getPFTValue(demand_ctx.pftType, 'resp_jones_kf', 0.227), ...% [jones] fast-pool turnover [1/hr]
        'resp_jones_Yg',     vegetation.getPFTValue(demand_ctx.pftType, 'resp_jones_Yg', 0.75), ... % [jones] growth yield
        'xsmr_recovery_days',vegetation.getPFTValue(demand_ctx.pftType, 'xsmr_recovery_days', 30));
    pftP_carbon.resp_night_hours = NitrogenVariables.resp_night_hours;   % [resp-substrate EXP] bruhn h (this step)
    pftP_carbon.labile_resp_days = vegetation.getPFTValue(demand_ctx.pftType, 'labile_resp_days', 1.0);   % [labile floor 2026-06-14] [d] soluble-sugar working-pool residence (computeLabileSugarTarget)
    pftP_carbon.C_labile_abs_min = vegetation.getPFTValue(demand_ctx.pftType, 'C_labile_abs_min', 0.01);  % [labile floor 2026-06-14] [g C/m2] small nonzero floor

    % ------------------------------------------------------------
    % Step 5: NSC 联动 C alloc (cap + overdraft + storage_drawdown 诊断)
    % ------------------------------------------------------------
    R_m_total_kg = npp_diag.R_m_leaf + npp_diag.R_m_root + npp_diag.R_m_stem + npp_diag.R_m_livecroot + npp_diag.R_m_fruit;
    % C/N 耦合模式 (cfg.Veg.cn_coupling_mode), 决定怎么把 N 分配信息传给 C alloc:
    %   flexible_cn (默认, CLM5 FlexibleCN): 传本步分配后器官 N 池 (pools_new.N_*), C 侧 NPP×f_alloc
    %     后按池 C:N 夹紧到 [target, target+offset] → C:N 带内浮动。
    %   c_follows_n (刚性): 传 A_N + actual_CN, A_C=A_N×actual_CN (弹性 inert → 池 C:N≡target)。
    %   parallel (旧): 不传 → A_C=NPP×f_alloc (C:N 无界漂移)。
    if strcmp(demand_ctx.cn_coupling_mode, 'flexible_cn')
        A_N_g = struct('mode', 'flexible_cn', ...
                       'N_leaf', pools_new.N_leaf,      'N_root',  pools_new.N_root_fine, ...
                       'N_stem', pools_new.N_livestem + pools_new.N_deadstem, ...   % [13c] total stem N for C-side clamp
                       'N_livestem', pools_new.N_livestem, 'N_deadstem', pools_new.N_deadstem, ...   % [13c per-pool clamp 2026-06-12]
                       'N_livecroot', pools_new.N_livecroot, 'N_deadcroot', pools_new.N_deadcroot, ...   % [coarse-root 2026-06-12]
                       'N_fruit', pools_new.N_fruit);
        [C_pools_after_alloc, nsc_diag] = vegetation.integrated.applyNSCCoupling( ...
            C_pools_prev, demand_ctx.GPP_for_demand, R_m_total_kg, npp_diag.R_g, ...
            NPP_kg, demand_ctx.DOY_for_demand, pftP_carbon, Delt_t, A_N_g, demand_ctx.phen_alloc, ...
            NitrogenVariables.C_xsmr(idx_prev));   % [xsmrpool] prev maintenance-resp debt pool
    else   % [cn-mode 2026-06-11] c_follows_n branch REMOVED (C-dominant only); this else = parallel (legacy, no clamp, A_N_g=[])
        [C_pools_after_alloc, nsc_diag] = vegetation.integrated.applyNSCCoupling( ...
            C_pools_prev, demand_ctx.GPP_for_demand, R_m_total_kg, npp_diag.R_g, ...
            NPP_kg, demand_ctx.DOY_for_demand, pftP_carbon, Delt_t, [], demand_ctx.phen_alloc, ...
            NitrogenVariables.C_xsmr(idx_prev));   % [xsmrpool] prev maintenance-resp debt pool
    end
    NitrogenVariables.C_xsmr(KT) = nsc_diag.C_xsmr;   % [xsmrpool] write-back this-step pool (<=0 = debt)

    NitrogenVariables.cum_GPP              = NitrogenVariables.cum_GPP + ...
                                             max(demand_ctx.GPP_for_demand, 0) * 1000 * Delt_t;
    NitrogenVariables.cum_R_a              = NitrogenVariables.cum_R_a + nsc_diag.R_a_actual_g;
    NitrogenVariables.cum_R_a_potential    = NitrogenVariables.cum_R_a_potential + nsc_diag.R_a_potential_g;
    NitrogenVariables.cum_storage_drawdown_count = NitrogenVariables.cum_storage_drawdown_count + nsc_diag.storage_drawdown_event;
    NitrogenVariables.nsc_diag(KT) = nsc_diag;

    % [actual organ resp 2026-06-15] Split the substrate-throttled (actual) Ra into organ components so
    %   downstream analysis reads actual CO2 directly instead of reconstructing potential*beta. theta
    %   (= actual R_m / potential R_m, INCLUDING the nodebt afford-cap, from applyNSCCoupling) applies
    %   uniformly to every organ's maintenance R_m; growth respiration R_g is NOT throttled. All in
    %   [kg C/(m2 s)], same as the potential organ R_m diagnostics. ra_above = leaf+stem+fruit Rm + R_g.
    theta_actual = nsc_diag.resp_theta_actual;
    NitrogenVariables.plantN_diag(KT).r_m_leaf_actual      = theta_actual * npp_diag.R_m_leaf;
    NitrogenVariables.plantN_diag(KT).r_m_root_actual      = theta_actual * npp_diag.R_m_root;
    NitrogenVariables.plantN_diag(KT).r_m_stem_actual      = theta_actual * npp_diag.R_m_stem;
    NitrogenVariables.plantN_diag(KT).r_m_livecroot_actual = theta_actual * npp_diag.R_m_livecroot;
    NitrogenVariables.plantN_diag(KT).r_m_fruit_actual     = theta_actual * npp_diag.R_m_fruit;
    r_g_actual = max(nsc_diag.R_a_actual_kg-theta_actual*R_m_total_kg,0);
    NitrogenVariables.plantN_diag(KT).r_g_actual           = r_g_actual;
    NitrogenVariables.plantN_diag(KT).ra_above_actual      = theta_actual * (npp_diag.R_m_leaf + npp_diag.R_m_stem + npp_diag.R_m_fruit) + r_g_actual;
    NitrogenVariables.plantN_diag(KT).ra_below_actual      = theta_actual * (npp_diag.R_m_root + npp_diag.R_m_livecroot);

    % [actual NPP basis 2026-06-15] NPP/NPP_net outputs must match the actual-Ra basis (autotrophic_respiration,
    %   NEE/Reco, plant C conservation all use actual Ra). computeNPP's NPP/NPP_net use the PRE-turn-down
    %   R_a_plant (potential), which would make npp_net + autotrophic_respiration(actual) != GPP and pull
    %   CUE back to the old potential basis (CUE 0.37 vs the actual 0.58). Overwrite with GPP - R_a_actual.
    NitrogenVariables.plantN_diag(KT).NPP_net = demand_ctx.GPP_for_demand - nsc_diag.R_a_actual_kg;          % actual, unclamped (true CUE)
    NitrogenVariables.plantN_diag(KT).NPP     = max(demand_ctx.GPP_for_demand - nsc_diag.R_a_actual_kg, 0);  % actual, clamped >=0

    % ------------------------------------------------------------
    % Step 6: N turnover (4 events: harvest/leaf_senesce/root_turnover/spring_release)
    % ------------------------------------------------------------
    % [mass-balance N reserve 2026-06-11] The N reserve evolves by mass balance from real N fluxes
    %   only -- it is NOT derived from the C reserve via a fixed storage C:N. The former
    %   CN_storage_target = 20 (N_to_storage = C_to_storage / 20) had no physical/literature source:
    %   "20" = 135/6.70, where N_storage_init = 6.70 had itself been set to make that ratio 20
    %   (circular), and the only CLM5 line cited for it (NutrientCompetitionFlexibleCNMod.F90:1660) is
    %   the LEAF C:N +/-10 band, not a reserve C:N. CLM5 holds retranslocated N in a SEPARATE
    %   mass-balanced pool (retransn_patch, CNVegNitrogenStateType.F90:60) filled by N resorption
    %   (leafn_to_retransn, CNPhenologyMod.F90:2433), never by dividing C storage by a C:N.
    %   Reserve N here = prev reserve + leaf/live-wood N resorption (added in turnoverAndHarvest)
    %   - spring remobilisation (springOnsetXfer). Storage C:N is therefore an EMERGENT diagnostic
    %   (higher after summer starch-only C refill, lower after autumn N resorption), not a target.
    %   The growing-season C reserve refill carries no matching N inflow on purpose: it is
    %   carbohydrate (starch), which is N-free.
    pools_turn_in = struct( ...
        'N_labile',     pools_new.N_labile, ...   % no fabricated C:N skim from labile into storage
        'N_leaf',       pools_new.N_leaf, ...
        'N_root_fine',  pools_new.N_root_fine, ...
        'N_livestem',   pools_new.N_livestem, ...   % [13c]
        'N_deadstem',   pools_new.N_deadstem, ...   % [13c]
        'N_livecroot',  pools_new.N_livecroot, ...  % [coarse-root 2026-06-12] for N live->dead mirror
        'N_deadcroot',  pools_new.N_deadcroot, ...  % [coarse-root 2026-06-12]
        'N_fruit',      pools_new.N_fruit, ...
        ... % [mass-balance N reserve 2026-06-11] prev reserve only; resorption added in turnoverAndHarvest
        'N_root_storage', NitrogenVariables.N_root_storage(idx_prev), ...
        'N_stem_storage', NitrogenVariables.N_stem_storage(idx_prev));
    pftP_turn = struct( ...
        'harvest_fraction',      vegetation.getPFTValue(demand_ctx.pftType, 'harvest_fraction'), ...
        'root_turnover_rate',    vegetation.getPFTValue(demand_ctx.pftType, 'root_turnover_rate'), ...
        'spring_release_days',   vegetation.getPFTValue(demand_ctx.pftType, 'spring_release_days'), ...
        'leaf_long_days',        vegetation.getPFTValue(demand_ctx.pftType, 'leaf_long_days'), ...
        'ndays_off',             vegetation.getPFTValue(demand_ctx.pftType, 'ndays_off'), ...   % Phase D leaf senescence ramp
        'retransl_fraction',     pftP_carbon.retransl_fraction, ...
        'f_storage_root',        vegetation.getPFTValue(demand_ctx.pftType, 'f_storage_root', 0.18), ...   % [resorbed-N split] root share of resorbed leaf N (Dong 2002; reserve root:stem partition)
        'is_perennial',          vegetation.getPFTValue(demand_ctx.pftType, 'is_perennial'), ...
        'has_fruit',             vegetation.getPFTValue(demand_ctx.pftType, 'has_fruit'), ...
        'has_storage_pool',      vegetation.getPFTValue(demand_ctx.pftType, 'has_storage_pool'), ...
        'lwtop',                 vegetation.getPFTValue(demand_ctx.pftType, 'lwtop'), ...   % [13c] live->dead stem turnover
        'CN_lw',                 vegetation.getPFTValue(demand_ctx.pftType, 'CN_lw'), ...   % [13c]
        'CN_dw',                 vegetation.getPFTValue(demand_ctx.pftType, 'CN_dw'));      % [13c]
    pftP_turn.winter_wheat_cn = vegetation.getPFTValue(demand_ctx.pftType, 'winter_wheat_cn', 0);
    pftP_turn.initial_seed_c = vegetation.getPFTValue(demand_ctx.pftType, 'initial_seed_c', 3);
    pftP_turn.CN_L = vegetation.getPFTValue(demand_ctx.pftType, 'CN_L');

    [pools_turn_out, turn_diag, turn_state_out] = vegetation.integrated.turnoverAndHarvest( ...
        pools_turn_in, NitrogenVariables.turnover_state, pftP_turn, env_turn, Delt_t);

    NitrogenVariables.N_labile(KT)    = pools_turn_out.N_labile;
    NitrogenVariables.N_leaf(KT)      = pools_turn_out.N_leaf;
    NitrogenVariables.N_root_fine(KT) = pools_turn_out.N_root_fine;
    NitrogenVariables.N_livestem(KT)  = pools_turn_out.N_livestem;   % [13c]
    NitrogenVariables.N_deadstem(KT)  = pools_turn_out.N_deadstem;   % [13c]
    NitrogenVariables.N_livecroot(KT) = pools_turn_out.N_livecroot;  % [coarse-root 2026-06-12]
    NitrogenVariables.N_deadcroot(KT) = pools_turn_out.N_deadcroot;  % [coarse-root 2026-06-12]
    NitrogenVariables.N_fruit(KT)     = pools_turn_out.N_fruit;
    NitrogenVariables.N_root_storage(KT) = pools_turn_out.N_root_storage;   % [dual-storage]
    NitrogenVariables.N_stem_storage(KT) = pools_turn_out.N_stem_storage;   % [dual-storage]
    NitrogenVariables.turnover_state  = turn_state_out;

    NitrogenVariables.cum_N_harvest     = NitrogenVariables.cum_N_harvest     + turn_diag.harvest_today;
    NitrogenVariables.cum_N_leaf_litter = NitrogenVariables.cum_N_leaf_litter + turn_diag.leaf_litter_today;
    NitrogenVariables.cum_N_root_litter = NitrogenVariables.cum_N_root_litter + turn_diag.root_litter_today;

    NitrogenVariables.turnover_diag(KT).harvest_today         = turn_diag.harvest_today;
    NitrogenVariables.turnover_diag(KT).leaf_litter_today     = turn_diag.leaf_litter_today;
    NitrogenVariables.turnover_diag(KT).root_litter_today     = turn_diag.root_litter_today;
    NitrogenVariables.turnover_diag(KT).storage_release_today = turn_diag.storage_release_today;
    NitrogenVariables.frPHU_diag(KT) = env_turn.frPHU;

    % ------------------------------------------------------------
    % Step 7: C turnover mirror (N→C via 实际池 C:N, FlexibleCN 下保浮动不回弹)
    % ------------------------------------------------------------
    % 传 pools_turn_in (周转前器官 N 池) → 凋落/根周转按实际池 C:N=C_pool/N_pool 扣碳,
    %   而非 CN_target, 否则周转把池 C:N 拉回 target (解第 3 把锁)。池在 target 时与旧式等价。
    % [E] onset-xfer gate context (consumed by Step 8 below)
    if isfield(env_turn, 'phen_state') && isfield(env_turn.phen_state, 'budburst_doy')
        bb_doy = env_turn.phen_state.budburst_doy;
    else
        bb_doy = NaN;
    end
    % [unify-spring-release 2026-06-06] use_leafc_sla no longer gates onset-xfer (removed);
    %   onset-xfer now runs for all storage PFTs regardless of LAI engine (see Step 8).
    C_pools_after_alloc.C_xsmr = NitrogenVariables.C_xsmr(KT);
    [C_pools_final, C_turn_diag] = vegetation.integrated.turnoverCarbonMirror( ...
        C_pools_after_alloc, turn_diag, pftP_carbon, Delt_t, pools_turn_in);

    % Crop cpool cleanup is a nonnegative atmospheric C flux. Negative xsmr
    % is a maintenance-debt state, so clearing it is recorded separately as
    % a bookkeeping source and never allowed to make respiration negative.
    crop_cleanup_g = C_turn_diag.C_crop_cleanup_to_atm_today;
    xsmr_writeoff_g = C_turn_diag.C_xsmr_debt_writeoff_today;
    if isfield(C_pools_final, 'C_xsmr')
        NitrogenVariables.C_xsmr(KT) = C_pools_final.C_xsmr;
    end
    if crop_cleanup_g ~= 0
        crop_cleanup_kg_s = crop_cleanup_g / (1000 * Delt_t);
        NitrogenVariables.cum_R_a = NitrogenVariables.cum_R_a + crop_cleanup_g;
        nsc_diag.C_xsmr = NitrogenVariables.C_xsmr(KT);
        nsc_diag.R_a_actual_g = nsc_diag.R_a_actual_g + crop_cleanup_g;
        nsc_diag.R_a_actual_kg = nsc_diag.R_a_actual_kg + crop_cleanup_kg_s;
        NitrogenVariables.nsc_diag(KT) = nsc_diag;
        NitrogenVariables.plantN_diag(KT).ra_above_actual = ...
            NitrogenVariables.plantN_diag(KT).ra_above_actual + crop_cleanup_kg_s;
        NitrogenVariables.plantN_diag(KT).NPP_net = ...
            demand_ctx.GPP_for_demand - nsc_diag.R_a_actual_kg;
        NitrogenVariables.plantN_diag(KT).NPP = max( ...
            NitrogenVariables.plantN_diag(KT).NPP_net, 0);
    end
    if xsmr_writeoff_g ~= 0
        NitrogenVariables.cum_C_xsmr_writeoff = ...
            NitrogenVariables.cum_C_xsmr_writeoff + xsmr_writeoff_g;
    end

    NitrogenVariables.C_labile(KT)     = C_pools_final.C_labile;
    NitrogenVariables.C_leaf(KT)       = C_pools_final.C_leaf;
    NitrogenVariables.C_root_fine(KT)  = C_pools_final.C_root_fine;
    NitrogenVariables.C_livestem(KT)   = C_pools_final.C_livestem;   % [13c]
    NitrogenVariables.C_deadstem(KT)   = C_pools_final.C_deadstem;   % [13c]
    NitrogenVariables.C_livecroot(KT)  = C_pools_final.C_livecroot;  % [coarse-root 2026-06-12]
    NitrogenVariables.C_deadcroot(KT)  = C_pools_final.C_deadcroot;  % [coarse-root 2026-06-12]
    NitrogenVariables.C_fruit(KT)      = C_pools_final.C_fruit;
    NitrogenVariables.C_root_storage(KT) = C_pools_final.C_root_storage;   % [dual-storage]
    NitrogenVariables.C_stem_storage(KT) = C_pools_final.C_stem_storage;   % [dual-storage]

    NitrogenVariables.cum_C_harvest     = NitrogenVariables.cum_C_harvest     + C_turn_diag.C_harvest_today;
    NitrogenVariables.cum_C_leaf_litter = NitrogenVariables.cum_C_leaf_litter + C_turn_diag.C_leaf_litter_today;
    NitrogenVariables.cum_C_root_litter = NitrogenVariables.cum_C_root_litter + C_turn_diag.C_root_litter_today;

    % ------------------------------------------------------------
    % Step 8: spring onset-xfer (storage->xfer->display, faithful CLM5 FruitTree perennial)
    % ------------------------------------------------------------
    % [unify-spring-release 2026-06-06] onset-xfer is now the SINGLE spring storage->growth
    % mechanism for ALL storage PFTs (gate decoupled from the LAI engine). At budburst it
    % mobilises fstor2tran of the NSC reserve into per-organ transfer pools (split by the
    % budburst allocation), then drains to display over ndays_stor. C and N moved at the
    % same fraction -> conservation closes + C:N preserved. CLM5 CNPhenologyMod.F90:1946-1996.
    % Previously gated on use_leafc_sla because onset-xfer was added to break the leafc_sla
    % death-spiral; but EPIC-engine storage PFTs need the same spring remobilisation, and
    % under EPIC the extra C/N into C_leaf is harmless (EPIC LAI = S-curve, not f(C_leaf)).
    % This replaces the old springStorageRelease (storage->labile), which is now removed.
    % See springOnsetXfer.m.
    use_onset = isfield(pftP_carbon,'has_storage_pool') ...
        && pftP_carbon.has_storage_pool == 1 && pftP_carbon.fstor2tran > 0;
    if use_onset
        DOY_o = demand_ctx.DOY_for_demand;
        ndays_stor = vegetation.getPFTValue(demand_ctx.pftType, 'ndays_stor', 30);
        if isnan(bb_doy) || DOY_o < bb_doy
            NitrogenVariables.onset_mobilized_flag = false;   % re-arm one-shot before budburst
        end
        if ~isnan(bb_doy) && DOY_o >= bb_doy && DOY_o <= bb_doy + ndays_stor
            onset = struct('mobilize', ~NitrogenVariables.onset_mobilized_flag, ...
                'counter', max(Delt_t, (bb_doy + ndays_stor - DOY_o) * 86400));
            % [onset-alloc routing 2026-06-28 Jizixing Zhu] DESIGN: this is the TWO-POOL
            % equivalent of CLM5's multi-pool spring onset. CLM5 keeps a separate storage pool
            % per organ (leafc_storage/livestemc_storage/...) and at budburst drains each into
            % its own organ (CNPhenologyMod.F90:626-643, strictly diagonal). COMPASS collapses
            % the aboveground storages into ONE wood-based reserve (C_stem_storage = nsc_frac x
            % stem wood; Step 3c), which is physiologically apt for deciduous trees (no leaves in
            % winter -> winter NSC lives in wood+root, and spring leaf flush is remobilised from
            % proximal woody NSC: Loescher 1990 HortScience 25:274; Chapin 1990 ARES 21:423).
            % So the single aboveground reserve must be SPLIT between leaf flush and new stem.
            % We split it by THIS PFT's spring growth allocation (same rule as normal growth) so
            % the reserve funds spring tissue in the proportions that tissue is being grown:
            %   * fruit tree (use_fruittree_alloc==1): FruitTree GDD-phase split
            %     (fruitTreeAllocation, fed the SAME demand_ctx.phen_alloc as Step 5);
            %   * natural woody (allometry, froot_leaf/stem_leaf>0): c_allometry ratios
            %     [1, froot_leaf, stem_leaf] -> normalised below to a_leaf:a_root:a_stem;
            %   * crops / legacy: fixed f_alloc_*.
            % DEVIATION from CLM5 (accepted 2026-06-28, kept by design): CLM5 sizes the leaf
            % flush from LAST year's accumulated leafc_storage; COMPASS uses THIS spring's
            % allocation ratio x the current wood reserve. Equivalent when allocation is stable
            % year-to-year; not a faithful per-organ port. Do NOT describe as "exactly CLM5".
            % Before this fix the split was ALWAYS f_alloc_*, so Apple's spring remobilisation
            % was pulled back to the Zanotelli fallback (0.12/0.08) instead of its FruitTree
            % allocation, even though normal NPP allocation already routes through FruitTree.
            if isfield(pftP_carbon,'use_fruittree_alloc') && pftP_carbon.use_fruittree_alloc == 1
                pa = demand_ctx.phen_alloc;
                gdd_on  = NaN; harv_on = 0;                        % defensive: match computeCarbonAllocation
                if isfield(pa,'gdd_budburst'), gdd_on  = pa.gdd_budburst; end
                if isfield(pa,'harvested'),    harv_on = pa.harvested;    end
                a_on = vegetation.carbon.fruitTreeAllocation(gdd_on, harv_on, pftP_carbon);
                fa = [a_on.leaf, a_on.root, a_on.stem];
            elseif pftP_carbon.froot_leaf > 0 || pftP_carbon.stem_leaf > 0
                fa = [1, pftP_carbon.froot_leaf, pftP_carbon.stem_leaf];   % c_allometry a_leaf:a_root:a_stem
            else
                fa = [pftP_carbon.f_alloc_leaf, pftP_carbon.f_alloc_root, pftP_carbon.f_alloc_stem];
            end
            fa = max(fa, 0);
            if sum(fa) <= 0, fa = [1 0 0]; end
            fa = fa / sum(fa);
            alloc = struct('leaf', fa(1), 'root', fa(2), 'stem', fa(3));
            % [diagonal dual-storage 2026-06-07] CLM5-faithful per-organ onset: each reserve
            %   is mobilised into its OWN organ's xfer (root reserve -> root; stem reserve ->
            %   leaf+livestem), with NO combined-total + fr_C re-split. The v59 fr_C re-split
            %   was a one-way zero-ratchet (once C_root_storage hit 0, fr_C=0 pinned root at 0
            %   forever). springOnsetXfer now reads/writes the dual reserve fields directly and
            %   debits each pool itself, so each reserve refills from its own allocation share
            %   and is never re-split. See springOnsetXfer.m + CLM5 CNCStateUpdate1Mod.F90:440-461.
            op = struct( ...
                'C_root_storage', NitrogenVariables.C_root_storage(KT), 'N_root_storage', NitrogenVariables.N_root_storage(KT), ...
                'C_stem_storage', NitrogenVariables.C_stem_storage(KT), 'N_stem_storage', NitrogenVariables.N_stem_storage(KT), ...
                'C_leaf',    NitrogenVariables.C_leaf(KT),          'N_leaf',    NitrogenVariables.N_leaf(KT), ...
                'C_root',    NitrogenVariables.C_root_fine(KT),     'N_root',    NitrogenVariables.N_root_fine(KT), ...
                'C_livestem',NitrogenVariables.C_livestem(KT),      'N_livestem',NitrogenVariables.N_livestem(KT), ...
                'C_leaf_xfer',    NitrogenVariables.C_leaf_xfer(idx_prev),     'N_leaf_xfer',    NitrogenVariables.N_leaf_xfer(idx_prev), ...
                'C_root_xfer',    NitrogenVariables.C_root_xfer(idx_prev),     'N_root_xfer',    NitrogenVariables.N_root_xfer(idx_prev), ...
                'C_livestem_xfer',NitrogenVariables.C_livestem_xfer(idx_prev), 'N_livestem_xfer',NitrogenVariables.N_livestem_xfer(idx_prev));
            op = vegetation.integrated.turnover.springOnsetXfer(op, onset, pftP_carbon.fstor2tran, alloc, Delt_t);
            % Diagonal write-back: springOnsetXfer already debited each reserve pool in place.
            NitrogenVariables.C_root_storage(KT) = op.C_root_storage;  NitrogenVariables.N_root_storage(KT) = op.N_root_storage;
            NitrogenVariables.C_stem_storage(KT) = op.C_stem_storage;  NitrogenVariables.N_stem_storage(KT) = op.N_stem_storage;
            NitrogenVariables.C_leaf(KT)=op.C_leaf;                NitrogenVariables.N_leaf(KT)=op.N_leaf;
            NitrogenVariables.C_root_fine(KT)=op.C_root;           NitrogenVariables.N_root_fine(KT)=op.N_root;
            NitrogenVariables.C_livestem(KT)=op.C_livestem;        NitrogenVariables.N_livestem(KT)=op.N_livestem;
            NitrogenVariables.C_leaf_xfer(KT)=op.C_leaf_xfer;          NitrogenVariables.N_leaf_xfer(KT)=op.N_leaf_xfer;
            NitrogenVariables.C_root_xfer(KT)=op.C_root_xfer;          NitrogenVariables.N_root_xfer(KT)=op.N_root_xfer;
            NitrogenVariables.C_livestem_xfer(KT)=op.C_livestem_xfer;  NitrogenVariables.N_livestem_xfer(KT)=op.N_livestem_xfer;
            if onset.mobilize, NitrogenVariables.onset_mobilized_flag = true; end
        end
    end

    % ------------------------------------------------------------
    % Step 9: autumn pruning (Group G) — prune annual wood increment to litter
    % ------------------------------------------------------------
    % CLM5 FruitTree: deadstemc_soy captured Jan 1; at autumn offset prune prune_fr of the
    %   annual deadstem increment (deadstemc - deadstemc_soy) to litter (mulch_pruning=1).
    %   CNPhenologyMod.F90:1675 (soy capture) + CNCStateUpdate1Mod.F90:304-307 (deadstem->litter).
    prune_fr = vegetation.getPFTValue(demand_ctx.pftType, 'prune_fr', 0);
    if prune_fr > 0
        DOY_p = demand_ctx.DOY_for_demand;
        if DOY_p <= 1   % year start: capture deadstem at start of year, re-arm one-shot prune
            NitrogenVariables.C_deadstem_soy = NitrogenVariables.C_deadstem(KT);
            NitrogenVariables.N_deadstem_soy = NitrogenVariables.N_deadstem(KT);
            NitrogenVariables.pruned_this_year = false;
        end
        % autumn senescence (leaf shedding, 2nd half of year): one-shot prune the wood increment
        if ~NitrogenVariables.pruned_this_year && DOY_p > 180 ...
                && isfield(turn_diag,'leaf_litter_today') && turn_diag.leaf_litter_today > 0
            dC_prune = prune_fr * max(0, NitrogenVariables.C_deadstem(KT) - NitrogenVariables.C_deadstem_soy);
            dN_prune = prune_fr * max(0, NitrogenVariables.N_deadstem(KT) - NitrogenVariables.N_deadstem_soy);
            NitrogenVariables.C_deadstem(KT)     = NitrogenVariables.C_deadstem(KT) - dC_prune;
            NitrogenVariables.N_deadstem(KT)     = NitrogenVariables.N_deadstem(KT) - dN_prune;
            NitrogenVariables.cum_C_prune_litter = NitrogenVariables.cum_C_prune_litter + dC_prune;
            NitrogenVariables.cum_N_prune_litter = NitrogenVariables.cum_N_prune_litter + dN_prune;
            NitrogenVariables.pruned_this_year   = true;
        end
    end

    % ------------------------------------------------------------
    % Step 10: background gap mortality [source-verified 2026-06-13] -- the MISSING turnover
    % ------------------------------------------------------------
    % CLM5 CNGapMortalityMod.F90:196-200: every woody pool loses r_mort fraction per year to
    %   litter/CWD: m_<pool>c_to_litter = <pool>c * m, m = r_mort/(days_per_yr*secsperday)
    %   (CNGapMortalityMod.F90:146). FruitTree keeps this for the apple PFT (r_mort=0.02,
    %   clm5_params.c171117__FruitTree.nc citrus idx 35). COMPASS was MISSING this entirely, so
    %   deadstem/deadcroot accumulated indefinitely (+49%/yr) -- only the deadstem-increment prune
    %   bounded stem, and coarse root had NO sink at all. This pool-proportional sink bounds the
    %   perennial woody pools at steady state. Applied to live+dead stem and live+dead coarse root;
    %   C and N at the SAME fraction (C:N preserved). Removed C/N -> cum_C_mortality/cum_N_mortality,
    %   summed into plant_c/n_residual in COMPASS.m like the other litter outflows (conservation).
    r_mort_p = vegetation.getPFTValue(demand_ctx.pftType, 'r_mort', 0);
    if r_mort_p > 0
        m_frac = r_mort_p * Delt_t / 31536000;   % [-] annual fraction -> per dt
        woodC = {'C_livestem','C_deadstem','C_livecroot','C_deadcroot'};
        woodN = {'N_livestem','N_deadstem','N_livecroot','N_deadcroot'};
        dC_mort = 0;  dN_mort = 0;
        for iw = 1:numel(woodC)
            dC = NitrogenVariables.(woodC{iw})(KT) * m_frac;
            dN = NitrogenVariables.(woodN{iw})(KT) * m_frac;
            NitrogenVariables.(woodC{iw})(KT) = NitrogenVariables.(woodC{iw})(KT) - dC;
            NitrogenVariables.(woodN{iw})(KT) = NitrogenVariables.(woodN{iw})(KT) - dN;
            dC_mort = dC_mort + dC;  dN_mort = dN_mort + dN;
        end
        NitrogenVariables.cum_C_mortality = NitrogenVariables.cum_C_mortality + dC_mort;
        NitrogenVariables.cum_N_mortality = NitrogenVariables.cum_N_mortality + dN_mort;
    end
end

% Change: 2026-09-02 Jizixing Zhu -- keep crop cleanup respiration
% nonnegative and ledger xsmr debt writeoff separately for exact closure.
