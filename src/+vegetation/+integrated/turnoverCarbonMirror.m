function [C_pools_out, C_diag] = turnoverCarbonMirror( ...
    C_pools_in, N_diag, pftParams, dt, N_pools)
%TURNOVERCARBONMIRROR  Mirror N turnover to C pools using actual pool C:N (M2-B.8 S3 + FlexibleCN).
%
%   [FlexibleCN 2026-05-30] 第 5 参 N_pools (周转前器官 N 池): 叶凋落/根周转/储备释放按
%     "实际池 C:N = C_pool/N_pool" 扣碳, 而非 CN_target — 否则周转把 FlexibleCN 浮动的池
%     C:N 拉回 target。池在 target 时与旧式 (×CN_target) 等价, 故向后兼容。不传 N_pools → 回退 CN_target。
%
%   M2-B.4 turnoverAndHarvest 在 N path 算 4 个 event 的 N 转移 (harvest_today,
%   leaf_litter_today, root_litter_today, storage_release_today). 这里用 CN_target
%   把 N 转换成 C, 同步更新 6 个 C 池.
%
%   Returns:
%     C_pools_out  struct  6 updated C pools
%     C_diag       struct  per-KT C diagnostic:
%       .C_harvest_today, .C_leaf_litter_today, .C_root_litter_today, .C_storage_release_today
%
%   Inputs:
%     C_pools_in   struct  6 C pools (上 KT 末值 already updated by computeCarbonAllocation)
%     N_diag       struct  M2-B.4 turnoverAndHarvest 输出的 4 字段 (g N/m²/KT):
%                          .harvest_today, .leaf_litter_today, .root_litter_today, .storage_release_today
%     pftParams    struct  CN_target_*, retransl_fraction
%     dt           [s]
%
%   Method (DAISY 简化路线, cn_coupling §7b.4):
%     - C harvest    = C_fruit × harvest_fraction   [v11: 直接按比例扣 C_fruit, 不用 CN 反推, 守恒]
%     - C leaf litter = N leaf litter × CN_L            (× (1-retransl_fraction) 已在 N path 算)
%     - C leaf retrans (隐) = leaf_senesce_total × CN_L × retransl_fraction → C_storage
%     - C root litter = N root litter × CN_target_root
%     - C spring release = N spring release × CN_L
%
%   重要: N path 的 leaf_litter_today 已经是 (1-retransl) × dN_leaf_senesce. 反推 dN_leaf_senesce
%   = N_diag.leaf_litter_today / (1-retransl). 总 leaf C 损失 = leaf_litter + retransl_to_storage.
%
%   守恒: ΔC_leaf = -leaf_litter - leaf_retransl_to_storage = -dN_leaf_senesce × CN_L (总损失)
%         ΔC_storage = +leaf_retransl - spring_release × CN_L
%
%   Reference:
%     cn_coupling §7b.4 (M2-B.8 C 池演化方程)
%     M2-B.4 turnover 子函数 4 个 (N path)
%
%   Variable class: C_pools — prognostic state, N_diag — input forcing, C_diag — diagnostic
%
%   Author: Zhu Jizixing 2026-05-21 (M2-B.8 S3)

    C_pools_out = C_pools_in;
    C_diag = struct('C_harvest_today',0, 'C_leaf_litter_today',0, ...
                    'C_root_litter_today',0, 'C_storage_release_today',0, ...
                    'C_cpool_to_atm_today',0, 'C_xsmr_to_atm_today',0, ...
                    'C_xsmr_debt_writeoff_today',0, ...
                    'C_crop_cleanup_to_atm_today',0);

    if isfield(pftParams, 'winter_wheat_cn') && pftParams.winter_wheat_cn == 1
        [C_pools_out, C_diag] = winterWheatTurnoverC( ...
            C_pools_out, C_diag, N_diag, pftParams);
        return;
    end

    CN_L     = pftParams.CN_L;
    CN_root  = pftParams.CN_target_root;
    retransl_frac = pftParams.retransl_fraction;
    % [v11] CN_fruit 不再用于采收 (改 C_fruit×harvest_fraction 直接扣, 见 (a))

    % [FlexibleCN 2026-05-30] 凋落/周转/释放按"实际池 C:N"= C_pool/N_pool 扣碳 (而非 CN_target),
    %   使 FlexibleCN 浮动的池不被周转拉回 target。无 N_pools 或 N≈0 → 回退 CN_target (向后兼容)。
    if nargin >= 5 && ~isempty(N_pools)
        cn_leaf         = poolCN(C_pools_in.C_leaf,         N_pools.N_leaf,         CN_L);
        cn_root         = poolCN(C_pools_in.C_root_fine,    N_pools.N_root_fine,    CN_root);
        % [dual-storage 2026-06-06] actual C:N of each reserve pool
        cn_root_storage = poolCN(C_pools_in.C_root_storage, N_pools.N_root_storage, CN_L);
        cn_stem_storage = poolCN(C_pools_in.C_stem_storage, N_pools.N_stem_storage, CN_L);
    else
        cn_leaf = CN_L;  cn_root = CN_root;
        cn_root_storage = CN_L;  cn_stem_storage = CN_L;
    end

    % --- (a) Apple harvest: C_fruit removed ---
    % [v11 2026-05-29] 物理: 采收摘走整个果实, C 和 N 各按同一采收比例 (harvest_fraction)
    %   从各自池扣 (N: N_fruit×hf in turnoverAndHarvest; C: C_fruit×hf 这里), 不互相 CN 换算。
    % 旧 bug (v10): C_harvest = N_harvest×CN_fruit 反推 — 当 grainCN=139 时反推 C 远超 C_fruit
    %   实际值 (C_fruit 按 f_alloc 碳流长, C:N≠139) → max(,0) clamp → 守恒账多记 ~29.6 gC。
    % 现: C 直接按比例扣 C_fruit, 自动守恒 (扣的就是池里有的)。
    if isfield(N_diag, 'harvest_today') && N_diag.harvest_today > 0 ...
            && isfield(pftParams, 'harvest_fraction') && pftParams.harvest_fraction > 0
        C_diag.C_harvest_today = C_pools_in.C_fruit * pftParams.harvest_fraction;   % [g C/m²]
        C_pools_out.C_fruit    = C_pools_in.C_fruit - C_diag.C_harvest_today;
        C_pools_out.C_fruit    = max(C_pools_out.C_fruit, 0);   % 数值安全 (理论不触发)
    end

    % --- (b) Leaf senescence + retransl: 总损失 = leaf_litter / (1-retransl) ---
    % [一致性临时件 2026-05-30] 这里叶 C 按 retransl_frac(0.7) 回收进 C_storage, 与 CLM5/FruitTree 不一致
    %   (CLM5+FruitTree: 叶 C → litter, 不回收; 仅 N 回收 leafn_to_retransn; storage 靠分配 fcur+aleafstor)。
    %   但它现役充当 FruitTree 采后 aleafstor=0.3 储备充值的功能 (秋季 storage 回补全靠它), 裸删会塌 storage。
    %   TODO Phase B: harvest/PHU 修好后, 改叶 C→litter + 接 FruitTree 采后 aleafstor=0.3 → C_storage 替代。
    %   (Dabbaghi 2022 GMD Table C1 aleafstor=0.3; CLM5 CNPhenologyMod.F90 leafc_to_litter)
    if isfield(N_diag, 'leaf_litter_today') && N_diag.leaf_litter_today > 0
        % Total N leaf senesce this KT (含 retrans 到 storage)
        if retransl_frac < 1
            dN_leaf_senesce_total = N_diag.leaf_litter_today / (1 - retransl_frac);
        else
            dN_leaf_senesce_total = 0;
        end
        dC_leaf_senesce_total = dN_leaf_senesce_total * cn_leaf;  % 总叶 C 损失 [g C/m²/KT] (实际池 C:N)

        % [Stage 2 followup 2026-06-10] Leaf C at senescence -> ALL to litter (CLM5/FruitTree-faithful:
        %   CNPhenologyMod.F90 leafc_to_litter; only N is retranslocated, C is shed). The previous
        %   leaf-C-retransl -> C_stem_storage was a documented TEMPORARY patch (comment above) to do the
        %   autumn storage refill; with the Stage 2 target-driven surplus-refill now doing that (capped
        %   at target), the uncapped leaf-C dump was the cause of stem_storage OVERSHOOTING its target
        %   (205 vs 149). Removing it fixes the overshoot AND restores CLM5 consistency. Leaf N retransl
        %   (N path) is unchanged — N is the scarce resource recovered.
        C_pools_out.C_leaf         = C_pools_in.C_leaf - dC_leaf_senesce_total;
        C_diag.C_leaf_litter_today = dC_leaf_senesce_total;   % all senesced leaf C -> litter (no storage retransl)
        C_pools_out.C_leaf         = max(C_pools_out.C_leaf, 0);
    end

    % --- (c) Root turnover: 100% to litter, no retrans ---
    if isfield(N_diag, 'root_litter_today') && N_diag.root_litter_today > 0
        C_diag.C_root_litter_today = N_diag.root_litter_today * cn_root;   % 实际池 C:N
        C_pools_out.C_root_fine    = C_pools_in.C_root_fine - C_diag.C_root_litter_today;
        C_pools_out.C_root_fine    = max(C_pools_out.C_root_fine, 0);
    end

    % --- (d) Spring storage release: root+stem reserve → C_labile ---
    % [dual-storage 2026-06-06] N side releases each reserve separately
    %   (storage_release_root_today / _stem_today); mirror C using each pool's actual C:N.
    if isfield(N_diag, 'storage_release_root_today') && N_diag.storage_release_root_today > 0
        dC_release_root = N_diag.storage_release_root_today * cn_root_storage;
        dC_release_root = min(dC_release_root, C_pools_out.C_root_storage);
        C_pools_out.C_root_storage = C_pools_out.C_root_storage - dC_release_root;
        C_pools_out.C_labile       = C_pools_out.C_labile       + dC_release_root;
        C_diag.C_storage_release_today = C_diag.C_storage_release_today + dC_release_root;
    end
    if isfield(N_diag, 'storage_release_stem_today') && N_diag.storage_release_stem_today > 0
        dC_release_stem = N_diag.storage_release_stem_today * cn_stem_storage;
        dC_release_stem = min(dC_release_stem, C_pools_out.C_stem_storage);
        C_pools_out.C_stem_storage = C_pools_out.C_stem_storage - dC_release_stem;
        C_pools_out.C_labile       = C_pools_out.C_labile       + dC_release_stem;
        C_diag.C_storage_release_today = C_diag.C_storage_release_today + dC_release_stem;
    end

    % --- (c2) Live->dead stem turnover C mirror [13c]: C_livestem -> C_deadstem (all C, no retrans) ---
    % CLM5 CNPhenologyMod.F90:2754: dC = livestemc * lwtop. C fully transfers (only N retransl, in N path).
    if isfield(pftParams,'lwtop') && pftParams.lwtop > 0 && C_pools_in.C_livestem > 0
        lwtop_frac = pftParams.lwtop * dt / 31536000;          % [1/yr] annual frac -> per dt
        dC_stem    = C_pools_in.C_livestem * lwtop_frac;       % [g C/m2]
        C_pools_out.C_livestem = C_pools_out.C_livestem - dC_stem;
        C_pools_out.C_deadstem = C_pools_out.C_deadstem + dC_stem;
        C_diag.C_stem_turnover_today = dC_stem;
    else
        C_diag.C_stem_turnover_today = 0;
    end

    % --- (c3) Live->dead coarse-root turnover C mirror [coarse-root 2026-06-12]: C_livecroot -> C_deadcroot ---
    % Belowground mirror of (c2) stem, SAME lwtop. CLM5 CNPhenologyMod.F90:2773 (livecrootc*lwtop ->
    % deadcrootc); state CNCStateUpdate1Mod.F90:285-286. C fully transfers (only N retransl, in N path).
    if isfield(pftParams,'lwtop') && pftParams.lwtop > 0 && isfield(C_pools_in,'C_livecroot') && C_pools_in.C_livecroot > 0
        lwtop_frac_cr = pftParams.lwtop * dt / 31536000;          % [1/yr] annual frac -> per dt
        dC_croot      = C_pools_in.C_livecroot * lwtop_frac_cr;   % [g C/m2]
        C_pools_out.C_livecroot = C_pools_out.C_livecroot - dC_croot;
        C_pools_out.C_deadcroot = C_pools_out.C_deadcroot + dC_croot;
        C_diag.C_croot_turnover_today = dC_croot;
    else
        C_diag.C_croot_turnover_today = 0;
    end

    % [E onset-xfer 2026-06-05] Spring storage->leaf bootstrap moved OUT of here into
    %   runPlantCNStep Step 8 (vegetation.integrated.turnover.springOnsetXfer): the faithful
    %   CLM5 FruitTree perennial onset moves ALL organs + N via per-organ transfer pools,
    %   which the C-only mirror here cannot express. See springOnsetXfer.m.
end

function [p, d] = winterWheatTurnoverC(p, d, nd, prm)
% Carbon mirror of the CoLM crop emergence/harvest event.
if isfield(nd, 'crop_emergence_event') && nd.crop_emergence_event
    p.C_leaf = p.C_leaf + p.C_fruit;
    p.C_fruit = 0;
end
if ~isfield(nd, 'crop_harvest_event') || ~nd.crop_harvest_event
    return;
end
seedC = min(p.C_fruit, prm.initial_seed_c);
% Change: 2026-09-09 Jizixing Zhu -- Changwu NM/NPM management: remove
% all standing shoot residue at harvest (Su Fuyuan thesis, pp.15-16).
% Harvest ledger is total off-field C, not grain-only yield. Seed reserve
% and root litter retain their existing ownership; apple branch unchanged.
d.C_harvest_today = p.C_fruit - seedC + p.C_leaf + p.C_livestem + p.C_deadstem;
d.C_leaf_litter_today = 0;
d.C_root_litter_today = p.C_root_fine + p.C_livecroot + p.C_deadcroot;
p.C_leaf = 0;
p.C_root_fine = 0;
p.C_livestem = 0;
p.C_deadstem = 0;
p.C_livecroot = 0;
p.C_deadcroot = 0;
p.C_fruit = seedC;

% CLM5 CNCStateUpdate1Mod.F90:435-461 and CoLM
% MOD_BGC_CNCStateUpdate1.F90:473-481 clear residual crop cpool and xsmrpool
% at harvest. C_xsmr<=0 is a signed maintenance-debt state, not physical
% carbon. Its reset is a bookkeeping source and must not be reported as
% negative atmospheric CO2. Positive xsmr, if supplied by a nonstandard
% caller, remains physical carbon and is released.
d.C_cpool_to_atm_today = p.C_labile;
if isfield(p, 'C_xsmr')
    d.C_xsmr_to_atm_today = max(p.C_xsmr, 0);
    d.C_xsmr_debt_writeoff_today = max(-p.C_xsmr, 0);
    p.C_xsmr = 0;
end
d.C_crop_cleanup_to_atm_today = d.C_cpool_to_atm_today + d.C_xsmr_to_atm_today;
p.C_labile = 0;
end

% Change: 2026-09-02 Jizixing Zhu -- separate nonnegative crop-pool CO2
% release from signed xsmr debt writeoff; preserve both in the C ledger.

function cn = poolCN(C_pool, N_pool, cn_fallback)
%POOLCN  实际池 C:N = C_pool/N_pool; N≈0 时回退 cn_fallback (=CN_target)。
    if N_pool > 1e-30
        cn = C_pool / N_pool;
    else
        cn = cn_fallback;
    end
end
