function [demand_ctx, root_uptake_ctx] = runPlantNDemand( ...
    KT, cfg, MeteoData, SoilVariables, NitrogenVariables, plantHydro, fluxes, canopy, current_date_KT, Delt_t, phen_alloc)
%RUNPLANTNDEMAND  Per-KT plant N demand computation (M2-B.2 DAISY = demand).
%
%   把 COMPASS.m 主循环 ~1085-1217 行的 demand 计算块抽出来. 输出两个 struct:
%     - demand_ctx     供 runPlantCNStep 复用 (NPP/D_N/T/GPP/PFT params)
%     - root_uptake_ctx 供 applyRootUptake 调用 (daisyParams + nacryParams)
%
%   Returns:
%     demand_ctx   struct, 字段:
%       .pftType            [string]       PFT 标签 (e.g. 'Apple')
%       .pftP_demand        struct         demand-side PFT params (r_m10, R_g_fraction, f_alloc_*, CN_*, has_fruit)
%       .GPP_for_demand     [kg C/(m²·s)]  fluxes.GPP (KT 值)
%       .Rd_leaf_for_demand [kg C/(m²·s)]  fluxes.Rd_leaf
%       .T_air_demand       [°C]           MeteoData.Ta(KT)
%       .T_soil_demand      [°C]           layer-mean SoilVariables.TT
%       .D_N_for_demand     struct         organ demand (M2-B.1b computeOrganDemand 输出)
%       .DOY_for_demand     [-]            day-of-year (1-365/366)
%     root_uptake_ctx  struct, 字段:
%       .daisyParams        struct         (r_root_m, SRL)
%       .nacryParams        struct         (Vmax_*, Km_*, Cmin_*, solver_*, demand_NH4/NO3)
%       .demand_NH4_g_m2_s  [g N/(m²·s)]   缓存供 reg.write
%       .demand_NO3_g_m2_s  [g N/(m²·s)]
%
%   Inputs:
%     KT                 [-]
%     cfg                ModelSettings (含 cfg.Veg, cfg.Grid)
%     MeteoData          forcing handle
%     SoilVariables      struct (含 TT)
%     NitrogenVariables  struct (含 N pool 时间序列, 上 KT 末值要的)
%     plantHydro         struct (含 r_root_m)
%     fluxes             struct (含 GPP, Rd_leaf, 可能为空/无字段 → 默认 0)
%     canopy             struct (含 LAI) — [F3后保留供接口兼容; demand 改 D_N×kstore 后不再用 LAI]
%     current_date_KT    datetime  当前 KT 对应日期 (供 DOY 计算)
%     Delt_t             [s]
%
%   Method:
%     1. 从 PFT 表加载所有 demand-side 参数 (避免主循环 4 处 vegetation.getPFTValue 散调)
%     2. 装填 LAI/N_target 推 N_labile_target = N_target_per_LAI × LAI
%     3. computeNPP (用 KT-1 N pool + 当前 GPP/Rd/T) → NPP_for_demand
%     4. computeOrganDemand (NPP × f_alloc / CN) → D_N
%     5. demand 公式: labile_deficit = max(0, N_target - N_labile_prev) → /dt
%        (M2-B.1b 修正: 不加 D_N_growth, 让根只补 labile 缺口, 见 cn_coupling §12.5)
%     6. NH4/NO3 拆分 (Vmax_NH4/(Vmax_NH4+Vmax_NO3) 比例)
%
%   Reference:
%     cn_coupling §12.5 (M2-B.2 demand 公式修正)
%     cn_coupling §7.2 (M2-B.4 N pool 演化公式)
%     COMPASS.m 主循环 1085-1217 (本函数抽取来源)
%
%   Variable class (per compass-coding-standard §1):
%     KT, cfg, MeteoData, SoilVariables, NitrogenVariables, plantHydro, fluxes — forcing/state (in)
%     demand_ctx, root_uptake_ctx — computed context (out, 给下游函数复用)
%
%   Author: Zhu Jizixing 2026-05-22 (M2-B.8c 二期重构 S1)
% Change: 2026-09-30 Wenjie Liu -- consume one compiled parameter view so
% site C:N overrides govern demand and downstream allocation consistently.
%   Change: 2026-09-30 Wenjie Liu -- forward the crop physical-NSC funding
%   gate to resp_only so N demand and applied C allocation share one budget.
%   Change: 2026-09-30 Wenjie Liu -- annual winter wheat takes only the
%   current growth N deficit after existing mobile N; no woody/next-year
%   reserve demand. Other PFT demand equations are unchanged.

    % ---------- PFT params (一次性加载, 后续函数复用) ----------
    pftType = cfg.Veg.PFT_type;

    assert(isfield(cfg.Veg, 'plant_cn') && ...
        isfield(cfg.Veg.plant_cn, 'resolved'), ...
        'runPlantNDemand:missingCompiledParameters', ...
        'Demand requires the immutable compiled PlantCN parameter view.');
    pftP_demand = cfg.Veg.plant_cn.resolved;
    pftP_demand.cn_flex_offset = cfg.Veg.cn_flex_offset;
    pftP_demand.fcur = cfg.Veg.fcur;
    pftP_demand.storage_cushion = cfg.Veg.storage_cushion;
    pftP_demand.laimx = cfg.Veg.laimx;

    % Nacry root-membrane PFT params (M1)
    Vmax_NH4_pft = pftP_demand.Vmax_NH4_umol_g_h;
    Vmax_NO3_pft = pftP_demand.Vmax_NO3_umol_g_h;
    Km_NH4_pft   = pftP_demand.Km_NH4;
    Km_NO3_pft   = pftP_demand.Km_NO3;
    Cmin_NH4_pft = pftP_demand.Cmin_NH4;
    Cmin_NO3_pft = pftP_demand.Cmin_NO3;
    SRL_pft      = pftP_demand.SRL;

    daisyParams = struct( ...
        'r_root_m', plantHydro.r_root_m, ...
        'SRL',      SRL_pft);

    % ---------- KT-1 末值 (KT=1 时用 init; 供 computeNPP 的滞后 organ N + F3 储备目标) ----------
    idx_prev = max(1, KT - 1);
    N_labile_for_demand   = NitrogenVariables.N_labile(idx_prev);     % [g N/m²] [F3] 储备缺口用
    N_root_prev_demand    = NitrogenVariables.N_root_fine(idx_prev);
    N_livestem_prev_demand = NitrogenVariables.N_livestem(idx_prev);  % [g N/m²] [13c] live stem N (respires)
    N_livecroot_prev_demand = NitrogenVariables.N_livecroot(idx_prev); % [g N/m²] [coarse-root] live coarse-root N (respires)
    N_fruit_prev_demand   = NitrogenVariables.N_fruit(idx_prev);      % [g N/m²] [FruitTree] fruit N (respires; non-fruit = 0)
    C_leaf_prev_demand    = NitrogenVariables.C_leaf(idx_prev);       % [g C/m²] [F3] 储备目标(萌发项)
    C_root_prev_demand    = NitrogenVariables.C_root_fine(idx_prev);  % [g C/m²] [F3]
    C_livestem_prev_demand = NitrogenVariables.C_livestem(idx_prev);  % [g C/m²] [F3/13c] sapwood-capacity term

    % ---------- 当前 KT 的 GPP/Rd/T ----------
    GPP_for_demand     = 0;
    Rd_leaf_for_demand = 0;
    if ~isempty(fluxes) && isstruct(fluxes) && isfield(fluxes, 'GPP')
        GPP_for_demand = fluxes.GPP;
        if isfield(fluxes, 'Rd_leaf'), Rd_leaf_for_demand = fluxes.Rd_leaf; end
    end
    if nargin < 11 || isempty(phen_alloc), phen_alloc = struct(); end
    crop_is_inactive = isfield(phen_alloc, 'is_crop') && phen_alloc.is_crop == 1 ...
        && (~isfield(phen_alloc, 'crop_live') || ~phen_alloc.crop_live);
    if crop_is_inactive
        % CLM5/CoLM clears crop C state at harvest. A lagged canopy LAI may
        % still produce fluxes in this adapter, but those fluxes must not
        % refill cpool or create new xsmr debt after crop_live becomes false.
        GPP_for_demand = 0;
        Rd_leaf_for_demand = 0;
    end
    T_air_demand   = MeteoData.Ta(KT);
    T_soil_demand  = mean(SoilVariables.TT(:));
    DOY_for_demand = day(current_date_KT, 'dayofyear');
    % ---------- Actual-NPP basis for N demand ----------
    % Compute respiration and substrate limitation first, without advancing C pools,
    % so root N demand and plant N allocation use the same NPP basis as C allocation.
    N_fruit_for_resp = N_fruit_prev_demand;
    if isfield(phen_alloc, 'is_crop') && phen_alloc.is_crop == 1 ...
            && (~isfield(phen_alloc, 'crop_live') || ~phen_alloc.crop_live)
        N_fruit_for_resp = 0; % retained grain seed is not a respiring fruit organ
    end
    [NPP_potential_for_demand, ~, npp_potential_diag] = vegetation.nitrogen.computeNPP( ...
        GPP_for_demand, Rd_leaf_for_demand, ...
        N_root_prev_demand, N_livestem_prev_demand, N_livecroot_prev_demand, N_fruit_for_resp, ...
        T_air_demand, T_soil_demand, pftP_demand);
    C_pools_prev_for_demand = struct( ...
        'C_labile',       NitrogenVariables.C_labile(idx_prev), ...
        'C_leaf',         NitrogenVariables.C_leaf(idx_prev), ...
        'C_root_fine',    NitrogenVariables.C_root_fine(idx_prev), ...
        'C_livestem',     NitrogenVariables.C_livestem(idx_prev), ...
        'C_deadstem',     NitrogenVariables.C_deadstem(idx_prev), ...
        'C_livecroot',    NitrogenVariables.C_livecroot(idx_prev), ...
        'C_deadcroot',    NitrogenVariables.C_deadcroot(idx_prev), ...
        'C_fruit',        NitrogenVariables.C_fruit(idx_prev), ...
        'C_root_storage', NitrogenVariables.C_root_storage(idx_prev), ...
        'C_stem_storage', NitrogenVariables.C_stem_storage(idx_prev));
    resp_night_hours_current = 0;
    if isfield(NitrogenVariables, 'resp_night_hours') && ~isempty(NitrogenVariables.resp_night_hours)
        resp_night_hours_current = NitrogenVariables.resp_night_hours;
    end
    if GPP_for_demand <= 1e-9
        resp_night_hours_current = resp_night_hours_current + 1;
    else
        resp_night_hours_current = 0;
    end
    pftP_demand.resp_night_hours = resp_night_hours_current;
    C_xsmr_prev_for_demand = 0;
    if isfield(NitrogenVariables, 'C_xsmr')
        C_xsmr_prev_for_demand = NitrogenVariables.C_xsmr(idx_prev);
    end
    R_m_total_for_demand = npp_potential_diag.R_m_leaf + npp_potential_diag.R_m_root + ...
                           npp_potential_diag.R_m_stem + npp_potential_diag.R_m_livecroot + ...
                           npp_potential_diag.R_m_fruit;
    [~, nsc_basis] = vegetation.integrated.applyNSCCoupling( ...
        C_pools_prev_for_demand, GPP_for_demand, R_m_total_for_demand, npp_potential_diag.R_g, ...
        NPP_potential_for_demand, DOY_for_demand, pftP_demand, Delt_t, ...
        struct('mode', 'resp_only'), phen_alloc, C_xsmr_prev_for_demand);
    NPP_for_demand = nsc_basis.NPP_actual_kg;
    % [Group D 2026-06-05] CLM5 FruitTree storage_flag: after harvest fcur=0 so ALL growth
    %   (both C and N, same fcur) routes to storage to refill reserves. Set here on
    %   pftP_demand so computeOrganDemand (N demand) + allocatePlantN (N alloc) + the
    %   forwarded pftP_carbon (C alloc) all see fcur=0. NutrientCompetitionFlexibleCNMod.F90:548-554.
    if isfield(phen_alloc, 'harvested') && phen_alloc.harvested == 1
        pftP_demand.fcur = 0;
    end
    % [peaklai 2026-06-05] Pass previous-step leaf C so computeOrganDemand can apply
    % the same LAI>=laimx leaf-cap as the C side (keeps C/N demand consistent).
    phen_alloc.C_leaf_prev = C_leaf_prev_demand;
    [D_N_for_demand, ~] = vegetation.nitrogen.computeOrganDemand( ...
        NPP_for_demand, DOY_for_demand, pftP_demand, phen_alloc);

    % --- N demand [F3 2026-05-31, 13c sapwood term 2026-06-01]: growth + labile refill ---
    % LPJ-GUESS need-based uptake (canexch.cpp:1163 ndemand = growth + storendemand;
    %   Smith et al. 2014 GMD 7:2261). Demand = current-step organ growth (D_N) plus
    %   topping the labile N pool up to a storage target; refill stops when labile is
    %   full (self-limiting), so N-rich sites do not accumulate luxury N unbounded.
    %   (Replaces D_N x kstore which kept over-absorbing 15% and let N_labile reach ~22.)
    % Storage target = max(sapwood-capacity term, next-year regrow term)
    %   LPJ-GUESS growth.cpp:1491 (tree branch):
    %     max_n_storage = max(cmass_sap * fnstorage / cton_leaf , retransn_nextyear)
    %   Both terms expressed in leaf C:N (cton_leaf_bg), not stem/root C:N.
    %   13c: sapwood-capacity term now uses real C_livestem (live wood = phloem+ray
    %   parenchyma, CLM5 CNGRespMod.F90:64), stable across season (was omitted pre-13c,
    %   leaving only the regrow term that tracked seasonal-inflated C_leaf -> N_labile ~22).
    %   fnstorage = 0.3 (LPJ-GUESS europe.ins:288, boreal/temperate tree default;
    %   0.05-0.15 herb/grass — apple is a tree, 0.3 applies).
    if pftP_demand.winter_wheat_cn == 1
        % CoLM NutrientCompetition (commit 626347a, lines 506-526) deploys
        % retranslocated N before requesting soil N. COMPASS uses N_labile
        % as the within-season mobile pool (including grain-fill retrans N).
        % This is not a perennial reserve or an additional uptake source.
        demand_total_g_m2_s = max(0, D_N_for_demand.total - ...
            N_labile_for_demand / Delt_t);
        if crop_is_inactive, demand_total_g_m2_s = 0; end
    else
    fnstorage = 0.3;                                                    % [-] LPJ-GUESS tree default
    sap_term = C_livestem_prev_demand / max(pftP_demand.CN_L, 1) * fnstorage;   % [g N/m²] sapwood capacity
    regrow_term = (C_leaf_prev_demand / max(pftP_demand.CN_L, 1) ...
                 + C_root_prev_demand / max(pftP_demand.CN_target_root, 1)) ...
                 * pftP_demand.nrelocfrac;                             % [g N/m²] next-year leaf+root regrow
    N_labile_target = max(sap_term, regrow_term);                      % [g N/m²] LPJ max(...)
    % [F3-combined 2026-06-02] LPJ-GUESS gates uptake on the COMBINED store
    % nstore() = nstore_labile + nstore_longterm, not labile alone
    % (LPJ-GUESS v4.1.1 guess.cpp:1679 storendemand = max(0, target - nstore())).
    % COMPASS holds the reserve in two pools (N_labile + N_storage). Comparing the
    % target against N_labile alone let root uptake keep pulling N while N_storage
    % was already full from autumn retranslocation, so the combined reserve grew to
    % ~2x target (N_labile reached 24.9 vs target 11.6). Internal transfers
    % (retransl, fcur->storage, spring release) only shuffle N between the two
    % pools; root uptake is the sole EXTERNAL inflow, so gating it on the combined
    % reserve bounds the total.
    % [dual-storage 2026-06-06] combined reserve = root + stem (single N_storage removed)
    N_storage_for_demand = NitrogenVariables.N_root_storage(idx_prev) + ...
                           NitrogenVariables.N_stem_storage(idx_prev);  % [g N/m2] winter reserve (root+stem)
    N_reserve_for_demand = N_labile_for_demand + N_storage_for_demand;  % [g N/m2] LPJ nstore()
    storage_deficit = max(0, N_labile_target - N_reserve_for_demand);   % [g N/m2], 0 when reserve full
    demand_total_g_m2_s = D_N_for_demand.total + storage_deficit / Delt_t;
    end

    % ---------- NH4/NO3 拆分 (生理偏好) ----------
    vmax_total = Vmax_NH4_pft + Vmax_NO3_pft;
    if vmax_total > 1e-6
        frac_NH4 = Vmax_NH4_pft / vmax_total;
        frac_NO3 = Vmax_NO3_pft / vmax_total;
    else
        frac_NH4 = 0.5; frac_NO3 = 0.5;
    end
    demand_NH4_g_m2_s = demand_total_g_m2_s * frac_NH4;
    demand_NO3_g_m2_s = demand_total_g_m2_s * frac_NO3;

    nacryParams = struct( ...
        'Vmax_NH4_umol_g_h', Vmax_NH4_pft, ...
        'Vmax_NO3_umol_g_h', Vmax_NO3_pft, ...
        'Km_NH4',            Km_NH4_pft, ...
        'Km_NO3',            Km_NO3_pft, ...
        'Cmin_NH4',          Cmin_NH4_pft, ...
        'Cmin_NO3',          Cmin_NO3_pft, ...
        'solver_tol',        cfg.Veg.root_solver_tol, ...
        'solver_max_iter',   cfg.Veg.root_solver_max_iter, ...
        'demand_NH4_g_m2_s', demand_NH4_g_m2_s, ...
        'demand_NO3_g_m2_s', demand_NO3_g_m2_s);

    % ---------- 组装 output structs ----------
    % [cn-mode 2026-06-11] C/N coupling mode. Model is C-DOMINANT only: default = flexible_cn
    %   (A_C=NPP*f_alloc + FlexibleCN clamp). The N-dominant 'c_follows_n' was REMOVED. NOTE: the
    %   old fallback here was 'c_follows_n', so whenever cfg.Veg.cn_coupling_mode was unset the model
    %   silently ran N-dominant — that latent bug is fixed by defaulting to flexible_cn.
    if isfield(cfg.Veg, 'cn_coupling_mode') && ~isempty(cfg.Veg.cn_coupling_mode)
        cn_coupling_mode = cfg.Veg.cn_coupling_mode;
    else
        cn_coupling_mode = 'flexible_cn';
    end
    if strcmp(cn_coupling_mode, 'c_follows_n')
        cn_coupling_mode = 'flexible_cn';   % c_follows_n removed -> force C-dominant
    end

    demand_ctx = struct( ...
        'pftType',            pftType, ...
        'pftP_demand',        pftP_demand, ...
        'plant_cn_params',    cfg.Veg.plant_cn, ...
        'GPP_for_demand',     GPP_for_demand, ...
        'Rd_leaf_for_demand', Rd_leaf_for_demand, ...
        'T_air_demand',       T_air_demand, ...
        'T_soil_demand',      T_soil_demand, ...
        'D_N_for_demand',     D_N_for_demand, ...
        'DOY_for_demand',     DOY_for_demand, ...
        'phen_alloc',         phen_alloc, ...
        'cn_coupling_mode',   cn_coupling_mode);

    root_uptake_ctx = struct( ...
        'daisyParams',        daisyParams, ...
        'nacryParams',        nacryParams, ...
        'demand_NH4_g_m2_s',  demand_NH4_g_m2_s, ...
        'demand_NO3_g_m2_s',  demand_NO3_g_m2_s);
end

% Change: 2026-08-24 Jizixing Zhu -- forward fixed CLM5/CoLM crop
% respiration and xsmr recovery parameters, and gate post-harvest crop
% assimilation out of the PlantCN ledger.

% Change: 2026-09-30 Wenjie Liu -- withdraw this round's APSIM stem WSC
%   and standing-senescent-leaf extension; retain the original COMPASS litter
%   and CoLM crop C/N allocation/harvest contracts. Earlier Change records are historical.
