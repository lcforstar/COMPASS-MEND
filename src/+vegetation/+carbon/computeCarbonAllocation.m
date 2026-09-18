function [C_pools_new, diag] = computeCarbonAllocation( ...
    C_pools_prev, GPP_kg, R_m_total_kg, R_g_kg, NPP_kg, current_DOY, pftParams, dt, A_N_g, phen_alloc, C_labile_target_g)
%COMPUTECARBONALLOCATION  Plant C pool allocation (6 pools, M2-B.8 DAISY 简化路线).
%
%   [v12 2026-05-29 C 主导 + C:N 弹性] 可选第 9 参 A_N_g (struct):
%     .leaf/.root/.stem/.fruit  [g N/m²/dt]  实际 N 分配 (allocatePlantN, 已 ×f_N_growth)
%     .CN_leaf/.CN_root/.CN_stem/.CN_fruit [gC/gN]  弹性 actual_CN (N 不足时上浮至 CN_max)
%     传入时走 C 主导: A_C_organ = A_N_organ × actual_CN_organ。N 足→C:N=target→A_C=NPP·f_alloc;
%       N 不足→C:N 上浮 C 照长; 撞 CN_max→C 被砍余碳留 C_labile。不传走旧式 A_C=NPP×f_alloc (C/N 独立漂移)。
%     依据: CLM5 FlexibleCN nscalar (NutrientCompetitionFlexibleCNMod.F90:1660-1667);
%       FruitTree C 主导 (Dabbaghi 2022 GMD 15:5167 p.4); cn_coupling §F.9。
%
%   每 KT 更新 6 个 C 池 (跟 M2-B.1b allocatePlantN 的 N 路径 mirror).
%
%   Returns:
%     C_pools_new  struct  updated 6 C pools:
%       .C_labile, .C_leaf, .C_root_fine, .C_livestem, .C_deadstem, .C_fruit, .C_storage  [g C/m²]
%     diag         struct  per-KT allocation diagnostic (per organ):
%       .A_C_leaf, .A_C_root, .A_C_stem, .A_C_fruit  [g C/m²/KT, 整 dt 量]
%       .R_a_total_kg   [kg C/(m²·s)]   total autotrophic respiration (M2-B.1b output)
%       .GPP_in_KT      [g C/m²]         carbon coming in this KT
%
%   Inputs:
%     C_pools_prev   struct, fields: .C_labile, .C_leaf, .C_root_fine, .C_livestem, .C_deadstem, .C_fruit, .C_storage
%                                     all [g C/m²], 上 KT 末值
%     GPP_kg         [kg C/(m²·s)]   gross primary productivity (SCOPE biochem)
%     R_m_total_kg   [kg C/(m²·s)]   maintenance respiration sum (R_m_leaf+root+stem from computeNPP)
%     R_g_kg         [kg C/(m²·s)]   growth respiration (computeNPP)
%     NPP_kg         [kg C/(m²·s)]   NPP (= GPP − R_a, computeNPP output)
%     current_DOY    [-]              for fruit gate
%     pftParams      struct           f_alloc_leaf/root/stem/fruit, has_fruit, CN_*
%     dt             [s]              time step (3600)
%
%   Method (Sitch 2003 LPJ + DAISY simplification):
%     1. GPP × dt → C_labile (NSC pool 积累)
%     2. R_a_total × dt → 出 C_labile (维持+生长呼吸, M2-B.1b 已算)
%     3. NPP × dt × f_alloc_X → C_organ X (跟 N 同 f_alloc, 物候 fruit gate 同 N)
%     4. C_labile = C_labile_prev + GPP×dt − R_a×dt − Σ A_C_X
%
%   Conservation invariant:
%     dC_labile + Σ dC_organ = GPP×dt − R_a×dt
%   (turnover/harvest/retransl 在 turnoverAndHarvest 内处理, 不在本函数)
%
%   Reference:
%     Sitch S et al (2003) GCB 9:161 — LPJ alloc (Eq. 11-14, NPP×f_alloc)
%     DAISY production.C — AssG × E_X × f_X(DS) → W_X 积累 (单态 5 池)
%     CLM5 CNCStateUpdate1Mod.F90:133 — cpool += psn × dt (类比 C_labile)
%     cn_coupling_path_and_literature_audit.md §7b.4 (M2-B.8 公式)
%
%   Variable class (per compass-coding-standard §1):
%     C_pools_prev/new — prognostic state (caller owns)
%     pftParams        — parameters (immutable)
%     dt, GPP, R_m, R_g, NPP — forcing slice
%     diag             — diagnostics (only output)
%
%   Author: Zhu Jizixing 2026-05-21 (M2-B.8 S2)
% Change: 2026-06-28 Jizixing Zhu — allometry branch treats empty []  froot_leaf/stem_leaf
%   as 0 (was 1+f1+[] -> [] -> 1/[] dimension error when a PFT sets one ratio but not the other).

    % [F2 stage-13a / #1] phen_alloc optional (10th arg): struct of phenology drivers
    % for allocation. .gdd_budburst [°C·d] drives the leaf taper; .DS_fruit [-] gates
    % fruit allocation on the GDD bloom clock. Absent -> no taper + pre-bloom (no fruit).
    if nargin < 10 || isempty(phen_alloc), phen_alloc = struct(); end
    % [labile sugar floor 2026-06-14] Optional 11th arg: soluble-sugar floor [g C/m2] protected before
    %   storage refill (computeLabileSugarTarget). 0 = legacy behaviour (sweep all surplus to storage).
    if nargin < 11 || isempty(C_labile_target_g), C_labile_target_g = 0; end
    if ~isfield(phen_alloc, 'gdd_budburst'), phen_alloc.gdd_budburst = NaN; end
    if ~isfield(phen_alloc, 'phen'),         phen_alloc.phen         = 1;   end
    if ~isfield(phen_alloc, 'harvested'),    phen_alloc.harvested    = 0;   end
    gdd_budburst = phen_alloc.gdd_budburst;     % [°C·d] GDD since bud-break
    phen         = phen_alloc.phen;             % [-] phenology activity (0=dormant .. 1=active)
    harvested    = phen_alloc.harvested;        % [0/1] fruit harvested flag (-> reserve refill)

    % --- Constants ---
    kg_to_g = 1000;                          % [g/kg]
    cv_kgC_per_s = kg_to_g * dt;             % [kg C/(m²·s)] × s × 1000 = g C/m² in dt

    % --- Step 1: GPP 进 C_labile ---
    GPP_in_KT  = max(GPP_kg, 0) * cv_kgC_per_s;            % [g C/m²/dt]

    % --- Step 2: R_a 出 C_labile ---
    R_a_total_kg = R_m_total_kg + R_g_kg;                  % [kg C/(m²·s)]
    R_a_out_KT   = max(R_a_total_kg, 0) * cv_kgC_per_s;    % [g C/m²/dt]

    % --- Step 3: Allocation ---
    NPP_in_KT = max(NPP_kg, 0) * cv_kgC_per_s;             % [g C/m²/dt]

    % [cn-mode 2026-06-11] Model is C-DOMINANT ONLY. Mode (A_N_g.mode): 'flexible_cn' (default, CLM5
    %   FlexibleCN: A_C=NPP*f_alloc, then clamp organ C:N to [CN_target, CN_target+offset]) or legacy
    %   'parallel' (NPP*f_alloc, no clamp). The N-dominant 'c_follows_n' (A_C=A_N*CN) is REMOVED:
    %   carbon leads, nitrogen follows within the flexible C:N band (NutrientCompetitionFlexibleCNMod.F90:1660).
    cn_mode = 'parallel';
    if nargin >= 9 && ~isempty(A_N_g)
        if isfield(A_N_g, 'mode'), cn_mode = A_N_g.mode;
        else,                      cn_mode = 'flexible_cn'; end   % backward-compat default -> C-dominant
    end

    if strcmp(cn_mode, 'c_follows_n')
        % [REMOVED 2026-06-11] N-dominant mode deleted; model is C-dominant only.
        error('computeCarbonAllocation:cn_mode_removed', ...
              'c_follows_n (N-dominant) removed; cn_coupling_mode must be flexible_cn (C-dominant).');
    else
        % parallel + flexible_cn 共用 C 侧分配: A_C = NPP × f_alloc (C 主导)。
        use_crop_alloc = isfield(phen_alloc, 'is_crop') && phen_alloc.is_crop == 1;
        if use_crop_alloc
            % CoLM crop allocation coefficients are calculated by the dedicated
            % winter-wheat state machine and consumed unchanged on the C side.
            f_leaf  = phen_alloc.crop_a_leaf;
            f_root  = phen_alloc.crop_a_root;
            f_stem  = phen_alloc.crop_a_stem;
            f_fruit = phen_alloc.crop_a_repr;
            if ~phen_alloc.crop_live
                f_leaf = 0; f_root = 0; f_stem = 0; f_fruit = 0;
            end
        elseif isfield(pftParams, 'use_fruittree_alloc') && pftParams.use_fruittree_alloc == 1
            % [13a 2026-06-04] FruitTree 4-phase GDD allocation (apple): root decline +
            % leaf J-Norman taper + stem decay + fruit (= remainder) rise, keyed to GDD since
            % bud-break. Subsumes the constant f_alloc + DOY/GDD fruit gate + separate leaf
            % taper. CLM5-FruitTree NutrientCompetitionFlexibleCNMod.F90:1645-1754.
            a_ft    = vegetation.carbon.fruitTreeAllocation(gdd_budburst, harvested, pftParams);
            f_leaf  = a_ft.leaf;
            f_root  = a_ft.root;
            f_stem  = a_ft.stem;
            f_fruit = a_ft.repr;
        else
            % [CLM5 all-PFT allometry 2026-06-13] Natural woody PFTs (NET/BDT/shrub) use the SAME
            %   c_allometry as the fruit tree, but fed by the FIXED nc ratios (NOT the GDD phase logic):
            %   f1=froot_leaf, f3=stem_leaf (NutrientCompetitionFlexibleCNMod.F90:413/426). Convert the
            %   leaf-relative ratios to NPP shares -- with S=1+f1+f3: a_leaf=1/S, a_root=f1/S, a_stem=f3/S
            %   (no fruit). Fed to the SAME denom normalization below, these reproduce CLM5's nlc form
            %   exactly (A_leaf=nlc, A_root=nlc*f1, A_stem=nlc*f3, A_coarse=nlc*f2*f3; verified algebraically).
            %   PFTs WITHOUT the nc ratios (grass/crop: froot_leaf=stem_leaf=0) fall back to the legacy
            %   constant f_alloc split with the DOY 150-290 fruit window.
            froot_leaf_p = 0;  stem_leaf_p = 0;
            % Treat an unset (empty []) ratio as 0: a PFT may define froot_leaf but not
            % stem_leaf (e.g. apple). Without the ~isempty guard, 1+f1+[] -> [] -> 1/[]
            % "dimensions must agree". [robustness 2026-06-28 Jizixing Zhu]
            if isfield(pftParams, 'froot_leaf') && ~isempty(pftParams.froot_leaf), froot_leaf_p = pftParams.froot_leaf; end
            if isfield(pftParams, 'stem_leaf')  && ~isempty(pftParams.stem_leaf),  stem_leaf_p  = pftParams.stem_leaf;  end
            if froot_leaf_p > 0 || stem_leaf_p > 0
                S_allom = 1 + froot_leaf_p + stem_leaf_p;       % c_allometry numerator pre-coarse (shares sum to 1)
                f_leaf  = 1 / S_allom;                          % a_leaf
                f_root  = froot_leaf_p / S_allom;               % a_root = f1/S (afroot)
                f_stem  = stem_leaf_p  / S_allom;               % a_stem = f3/S
                f_fruit = 0;                                    % natural trees: no reproductive pool
            else
                % Legacy constant f_alloc (grass/crop without nc allometry), DOY 150-290 fruit window.
                fruit_active = (pftParams.has_fruit == 1) && ...
                               (current_DOY >= 150) && (current_DOY <= 290);
                if fruit_active
                    f_leaf  = pftParams.f_alloc_leaf;
                    f_root  = pftParams.f_alloc_root;
                    f_stem  = pftParams.f_alloc_stem;
                    f_fruit = pftParams.f_alloc_fruit;
                else
                    % Fruit window closed: fruit share -> stem (woody/storage), matching
                    % computeOrganDemand.m (N side); preserves C-dominant C/N consistency.
                    f_leaf  = pftParams.f_alloc_leaf;
                    f_root  = pftParams.f_alloc_root;
                    f_stem  = pftParams.f_alloc_stem + pftParams.f_alloc_fruit;
                    f_fruit = 0;
                end
            end
            if pftParams.CN_target_stem <= 0     % Annual PFT: 跳 stem alloc
                f_stem = 0;
            end
        end
        % [F2 dormancy gate 2026-06-04] Deciduous trees grow no leaves while dormant.
        % phen < 0.1 (winter / pre-budburst, canopy ~ ALAI_MIN floor) -> zero leaf
        % allocation; the unallocated carbon stays in the C_labile NSC reserve for next
        % spring's flush. Without this, f_leaf×NPP keeps feeding the leaf pool all winter
        % (ALAI_MIN floor -> winter GPP -> NPP) and C_leaf creeps up from 0. CLM5 sets
        % dormant aleaf~1e-5 with C routed to storage (CNPhenologyMod offset/dormancy).
        if ~use_crop_alloc && phen < 0.1
            f_leaf = 0;
        end

        % [peaklai 2026-06-05] LAI saturation cap: once canopy LAI reaches the
        % physiological peak (laimx), stop feeding the leaf pool and route that
        % carbon to the overwintering reserve (C_storage). Without this the leafc_sla
        % engine has no leaf-area negative feedback (constant f_alloc keeps inflating
        % C_leaf -> LAI 4.35 vs obs ~2.8).
        % WHY storage, not the current-year root pool: CLM5 redirects the capped leaf
        % share to root, but in CLM5 that root carbon is itself split fcur / (1-fcur)
        % into current root vs *frootc_storage* overwintering reserve
        % (NutrientCompetitionFlexibleCNMod.F90:522-523). COMPASS has no per-organ
        % storage sub-pools — its single C_storage IS the overwintering reserve. With
        % apple fcur=0.95 routing the capped C to the current-year root pool (a no-return
        % sink that still respires in winter) starved C_storage (yr1 end 41 vs 75 needed
        % for spring onset) and killed year-2 (LAI 0.06). Routing the capped leaf C to
        % C_storage is the architecture-equivalent of CLM5's (1-fcur)->frootc_storage
        % overwintering portion, restoring the multi-year cycle.
        % Source: CLM5-FruitTree CNVegStructUpdateMod.F90:173-181 (woody&perennial:
        %   dormant_flag==0 & tlai>=laimx -> peaklai=1; dormant -> peaklai=0)
        %         + NutrientCompetitionFlexibleCNMod.F90:1507-1511 (aleaf->~0) and
        %         522-523 (redirected C split into storage via 1-fcur).
        % Lagged by one step (uses previous C_leaf -> LAI), matching CLM5's
        % "enforce the cap in the following tstep" (FlexibleCN:1524-1527).
        if ~use_crop_alloc && phen >= 0.1 && isfield(pftParams, 'laimx') && pftParams.laimx > 0 ...
                && isfield(pftParams, 'SLA0') && pftParams.SLA0 > 0
            % LAI from previous-step leaf C (same form as computeLAIfromPhen leafc_sla)
            if isfield(pftParams, 'SLAm') && pftParams.SLAm > 1e-9
                LAI_prev = (pftParams.SLA0 / pftParams.SLAm) * ...
                           (exp(C_pools_prev.C_leaf * pftParams.SLAm) - 1);
            else
                LAI_prev = pftParams.SLA0 * C_pools_prev.C_leaf;
            end
            if LAI_prev >= pftParams.laimx
                % [③ peaklai no-redirect 2026-06-12] At peak LAI stop leaf growth (f_leaf=0).
                %   The OLD code redirected the freed share to ROOT (f_root += f_leaf), which
                %   inflated the fine root (root doubled 172->323 under nodebt; Zanotelli 2013
                %   apple fine-root NPP is ~1%). Fine root is now bound to leaf (A_C_root =
                %   froot_leaf*A_C_leaf below), so f_leaf=0 auto-drops fine-root growth at peak;
                %   the freed carbon flows to the FRUIT remainder (Step 3b ④). N side mirrors
                %   this (computeOrganDemand: f_leaf=0, no root redirect).
                f_leaf = 0;
            end
        end

        % [CLM5-exact allocation 2026-06-13] FruitTree fork NutrientCompetitionFlexibleCNMod.F90:
        %   cpool_to_<organ> = nlc*f_<organ>, nlc = plant_calloc/c_allometry, c_allometry(woody) =
        %   (1+g1)(1 + f1 + f3*(1+f2))  (:1897). For the crop/fruit tree f1=aroot/aleaf (:446),
        %   f3=astem/aleaf (:447), f2=croot_stem (:414 param, NO crop override). The 6 fluxes (:555-564):
        %   leaf=nlc, froot=nlc*f1, livestem=nlc*f3*f4, deadstem=nlc*f3*(1-f4),
        %   livecroot=nlc*f2*f3*f4, deadcroot=nlc*f2*f3*(1-f4). Since aleaf+aroot+astem+arepr=1, this
        %   reduces algebraically to dividing every organ by denom = 1 + a_stem*f2 (the coarse root is
        %   the extra f2*a_stem term in c_allometry). So A_leaf=a_leaf/denom, A_root=a_root/denom (afroot,
        %   Eq3), A_stem=a_stem/denom, A_coarse=f2*a_stem/denom (= f2*livestem, src :563), A_fruit=a_repr/
        %   denom (Eq7); total = (1+a_stem*f2)/denom * NPP = NPP (conserved). This REVERTS the 2026-06-12
        %   froot_leaf-bind + fruit-remainder, which deviated from the paper's Eq3/Eq7. [CLM5 all-PFT
        %   2026-06-13] now UNIVERSAL: natural woody PFTs feed a_leaf/a_root/a_stem from the nc f1/f3
        %   ratios (else branch above), grass/crop keep f_alloc; the SAME denom normalization serves all.
        use_ft = isfield(pftParams, 'use_fruittree_alloc') && pftParams.use_fruittree_alloc == 1;
        f2_croot = pftParams.croot_stem_alloc;             % [-] f2 coarse-root:stem (citrus 0.08; natural woody 0.30; grass/crop 0)
        denom = 1 + f_stem * f2_croot;                     % c_allometry normalization, UNIVERSAL (f2=0 -> denom=1, no coarse); apple denom value unchanged
        A_C_leaf  = f_leaf * NPP_in_KT / denom;            % [g C/m²/dt] leaf (peaklai/dormancy-capped)
        A_C_stem  = f_stem * NPP_in_KT / denom;            % aboveground wood
        A_C_livestem = A_C_stem * pftParams.flivewd;       % [13c] default live/dead split (clamp overrides)
        A_C_deadstem = A_C_stem * (1 - pftParams.flivewd);

        if strcmp(cn_mode, 'flexible_cn')
            % [FlexibleCN 2026-05-30] 池 C:N 带内浮动: 夹紧各器官 A_C 使新池 C:N ≤ CN_target+offset,
            %   超出 C 不分配 → 经 Step4 C_labile balance 留存。浮动来源 = N 充足期(满额=target) 与
            %   N 饿死期(C 照进/N≈0 → C:N 推高撞顶) 交替, 池 C:N 落 [target, target+offset]。
            %   新组织 C:N = A_C/A_N = CN_target/f_N ≥ target, 只顶上界(不破下界), 故只夹 CN_max。
            %   A_N_g.N_* = 本步分配后器官 N 池 (gN/m²)。
            %   对标 CLM5 NutrientCompetitionFlexibleCNMod.F90:1660-1664 (leafcn 在 [target±10] 浮动)。
            off = pftParams.cn_flex_offset;                 % [gC/gN] 带半幅
            A_C_leaf  = clampPoolCN(A_C_leaf,  C_pools_prev.C_leaf,      A_N_g.N_leaf,  pftParams.CN_L           + off);
            if pftParams.CN_target_stem > 0
                % [13c per-pool clamp 2026-06-12] Clamp livestem vs CN_lw and deadstem vs CN_dw
                %   SEPARATELY, then recombine. The old single clamp used TOTAL stem C:N vs
                %   CN_target_stem(=60); the dead-wood pool (C:N ~500) drove the total to ~198 >> 60,
                %   so ALL fresh stem C was clamped to 0 (A_C_stem=0 all year -> stem never grew ->
                %   coarse-root-via-stem = 0). Split A_C_stem by flivewd, clamp each wood type within
                %   its own C:N ceiling (live CN_lw, dead CN_dw). A_C_ls/A_C_ds reused at pool update.
                A_C_livestem = clampPoolCN(A_C_stem * pftParams.flivewd,       C_pools_prev.C_livestem, A_N_g.N_livestem, pftParams.CN_lw + off);
                A_C_deadstem = clampPoolCN(A_C_stem * (1 - pftParams.flivewd), C_pools_prev.C_deadstem, A_N_g.N_deadstem, pftParams.CN_dw + off);
                A_C_stem = A_C_livestem + A_C_deadstem;
            end
        end

        % [CLM5-exact 2026-06-13] fine root = a_root (afroot, NutrientCompetitionFlexibleCN:446 Eq3)
        %   * NPP/denom -- NOT froot_leaf*leaf (that was the 2026-06-12 deviation). Coarse root =
        %   f2*livestem (src :563). Derived AFTER leaf/stem clamps so bound organs scale with clamped parent.
        A_C_root = f_root * NPP_in_KT / denom;                     % afroot (or legacy f_alloc; denom=1 non-FT)
        A_C_livecroot = f2_croot * A_C_livestem;                   % = f2*livestem (src :563); 0 if f2=0 (non-woody)
        A_C_deadcroot = f2_croot * A_C_deadstem;
        if strcmp(cn_mode, 'flexible_cn')
            A_C_root = clampPoolCN(A_C_root, C_pools_prev.C_root_fine, A_N_g.N_root, pftParams.CN_target_root + off);
            if pftParams.croot_stem_alloc > 0
                A_C_livecroot = clampPoolCN(A_C_livecroot, C_pools_prev.C_livecroot, A_N_g.N_livecroot, pftParams.CN_lw + off);
                A_C_deadcroot = clampPoolCN(A_C_deadcroot, C_pools_prev.C_deadcroot, A_N_g.N_deadcroot, pftParams.CN_dw + off);
            end
        end

        % [CLM5-exact 2026-06-13] fruit = a_repr (arepr, NutrientCompetitionFlexibleCN Eq7) * NPP/denom
        %   -- NOT NPP-remainder (the 2026-06-12 deviation). fruitTreeAllocation already returns a_repr
        %   as the phase remainder (1 - aleaf - aroot - astem); here it is just scaled by 1/denom.
        A_C_fruit = f_fruit * NPP_in_KT / denom;
        if strcmp(cn_mode, 'flexible_cn') && pftParams.CN_target_fruit > 0
            A_C_fruit = clampPoolCN(A_C_fruit, C_pools_prev.C_fruit, A_N_g.N_fruit, pftParams.CN_target_fruit + off);
        end
    end

    % --- Step 3c: reserve loading (target-driven, Stage 2) ---
    % [STAGE 2 2026-06-10] The old (1-fcur) fixed-fraction loading + the post-harvest storage_flag
    %   (fcur=0) are REMOVED; the reserve now refills from C_labile surplus toward its wood-based
    %   target (below). fcur and A_C_organ_total were left dead by that change and are deleted here.
    % [dual-storage 2026-06-06] overwinter reserve split coarse-root(0.55)+stem(0.45).
    %   (1-fcur) organ growth + peaklai-capped leaf C both route to the two reserve pools
    %   by the same physiological fraction (Chapin 1990; Loescher 1990, see setInitialNitrogen).
    %   peaklai capped C: instead of orphaning to a pure-C single storage (old B5 C/N drift),
    %   it now enters the same root/stem reserves whose N is fed symmetrically on the N side.
    % [dual-storage 2026-06-06] Route (1-fcur) reserve BY SOURCE ORGAN (not a fixed split):
    %   root's (1-fcur) share -> root reserve; stem's -> stem reserve; leaf+fruit have no
    %   own reserve pool so their overwinter share + the peaklai-capped leaf C go to the
    %   stem (sapwood) reserve, the nearest woody store via phloem. The two reserves then
    %   evolve independently (filled by their own organ's source, drained by their own organ's
    %   spring growth); their ratio drifts naturally rather than being pinned to a constant.
    % [STAGE 2 2026-06-10] fc1=(1-fcur) fixed-fraction loading REMOVED; reserve is now target-driven
    %   (Step 3c below refills from surplus toward nsc_frac x wood). fcur/storage_flag superseded here.
    % [dual-storage] (1-fcur) reserve routed BY SOURCE ORGAN. peaklai already folded its
    %   freed leaf share into f_root above (CLM5-faithful), so the root reserve naturally
    %   receives that carbon via fc1*A_C_root — no separate peaklai-to-storage stream needed.
    % [coarse-root-refill 2026-06-07] CLM5 belowground reserve = frootc_storage(fine) +
    %   livecrootc_storage + deadcrootc_storage(coarse). The coarse-root storages are fed from
    %   the stem allocation stream scaled by f2=croot_stem (coarse-root:stem ratio):
    %   livecroot_storage = nlc*f2*f3*f4*(1-fcur) and deadcroot_storage = nlc*f2*f3*(1-f4)*(1-fcur)
    %   (CLM5 NutrientCompetitionFlexibleCNMod.F90:530,532), summing to nlc*f2*f3*(1-fcur).
    %   COMPASS has NO explicit coarse-root organ pool: A_C_root here is FINE root only (apple
    %   a_ft.root) and A_C_stem is ABOVEGROUND stem only (apple a_ft.stem). The coarse root is
    %   missing from the allocation split, so its reserve contribution is reconstructed by taking
    %   croot_stem*A_C_stem*(1-fcur) from the (aboveground) stem stream and routing it to the
    %   belowground reserve. This does NOT double-count: fine root (A_C_root) and reconstructed
    %   coarse root (croot_stem*A_C_stem) are distinct belowground parts. Without it, the root
    %   reserve only got the thin fine-root 5% flux and decayed year over year (112->61->38 over
    %   2yr) while stem reserve grew fat from autumn leaf retransl.
    %   croot_stem is PFT-table sourced: apple=0.225 (Zanotelli et al. 2013 Biogeosciences 10:3089
    %   woodBG/woodAG=16.3/72.3, coarse-root:aboveground-wood); natural tree/shrub=0.30 (CLM5).
    %   PFT-specific, so it MUST come from the PFT table, not be hardcoded. PFTs that do not set
    %   it default to 0 (no coarse-root split = legacy behaviour, safe for non-woody).
    if isfield(pftParams, 'croot_stem'), croot_stem = pftParams.croot_stem; else, croot_stem = 0; end  % [-] coarse-root:stem ratio
    % [target-reserve 2026-06-10] Wood-based reserve target for each pool (target = nsc_frac x own wood;
    %   coarse-root wood reconstructed via croot_stem). CONSUMED by the surplus-refill in Step 3c below.
    %   See docs/superpowers/plans/2026-06-10-target-driven-nsc-reserve.md
    if isfield(pftParams, 'nsc_frac'), nsc_frac = pftParams.nsc_frac; else, nsc_frac = 0.07; end  % [-] woody NSC C-fraction
    C_stem_wood_prev = C_pools_prev.C_livestem + C_pools_prev.C_deadstem;         % [g C/m2] perennial AG wood
    [target_stem_storage, target_root_storage] = ...
        vegetation.carbon.computeReserveTarget(C_stem_wood_prev, croot_stem, nsc_frac);
    % [target-driven 2026-06-10 STAGE 2] Reserve loading REPLACED. Instead of the fixed (1-fcur) skim of
    %   each organ's allocation, each woody reserve refills from the SURPLUS of C_labile (what is left
    %   after the current-year display growth is funded) TOWARD its wood-based target, with a FATES
    %   smooth taper (fills fast when depleted, tapers to 0 at the target -> neither starves nor bloats).
    %   Surplus-sourced -> never force-drains C_labile (the v10 cushion-failure mode). Display organs now
    %   receive the FULL allocation (no fcur skim). Autumn refill is automatic (post-harvest the display
    %   demand drops -> surplus rises -> storage tops up), so the old fcur=0 storage_flag is superseded.
    %   fates/parteh/PRTAllometricCarbonMod.F90:598-613 (smooth taper); Massoud et al. 2019 GMD 12:4133.
    display_total = A_C_leaf + A_C_root + A_C_stem + A_C_fruit + A_C_livecroot + A_C_deadcroot;  % [g C/m2] full current growth (incl coarse root)
    avail   = C_pools_prev.C_labile + GPP_in_KT - R_a_out_KT - display_total;     % [g C/m2] labile left after display
    % [labile sugar floor 2026-06-14] Only carbon ABOVE the soluble-sugar floor is surplus for starch
    %   storage; C_labile_target_g reserves a working sugar pool so C_labile is not swept to ~0 every step.
    %   See docs/nsc_labile_sugar_floor_design.md.
    surplus = max(avail - C_labile_target_g, 0);                                  % [g C/m2] excess above sugar floor refills storage
    refillTo = @(store, tgt, cap) min(max(tgt - store, 0), ...
                  cap * max(exp(-(store / max(tgt, 1e-9))^4) - exp(-1), 0));      % FATES smooth-taper toward target
    C_to_root_storage = refillTo(C_pools_prev.C_root_storage, target_root_storage, surplus);
    C_to_stem_storage = refillTo(C_pools_prev.C_stem_storage, target_stem_storage, surplus - C_to_root_storage);
    C_to_storage_tot  = C_to_root_storage + C_to_stem_storage;                    % total to reserve (labile balance)
    % [labile sugar floor 2026-06-14] Starch -> sugar remobilisation: when labile falls below its floor
    %   (deficit; mutually exclusive with surplus>0 since surplus>0 => avail>target => labile_pre>=target),
    %   refill it from the starch reserve toward the target, split by available reserve. This funds the
    %   shortfall BEFORE the hard zero/debt guard in applyNSCCoupling Step 5 (which then only handles any
    %   residual after BOTH reserves are exhausted). Martinez-Vilalta et al. 2016 Ecol Monogr 86:495
    %   (starch reservoir remobilises to keep soluble sugars above threshold).
    C_labile_pre   = avail - C_to_storage_tot;                                    % [g C/m2] labile before remob
    labile_deficit = max(C_labile_target_g - C_labile_pre, 0);                    % [g C/m2] shortfall below floor
    root_stor_av   = max(C_pools_prev.C_root_storage, 0);                         % [g C/m2] reserve available
    stem_stor_av   = max(C_pools_prev.C_stem_storage, 0);                         % [g C/m2] reserve available
    stor_av_tot    = root_stor_av + stem_stor_av;                                 % [g C/m2]
    C_storage_to_labile = min(labile_deficit, stor_av_tot);                       % [g C/m2] starch -> sugar this step
    draw_root      = C_storage_to_labile * (root_stor_av / max(stor_av_tot, 1e-12));   % split by available reserve
    draw_stem      = C_storage_to_labile - draw_root;
    % (fcur skim removed: A_C_* stay full; storage filled from surplus above, not an fcur fraction)

    % --- Step 4: C 池更新 (mass balance) ---
    % C_labile balance: dC_labile = GPP×dt − R_a×dt − (display_total + C_to_storage_tot)
    %   [STAGE 2] display_total = full organ allocation (no fcur skim); C_to_storage_tot = surplus-based
    %   refill toward target. Conservation: dC_labile + Σd(organ) + d(storage) = GPP − R_a (closes;
    %   storage just moves C between labile and reserve). C_labile stays ≥0 when avail≥0 (refill≤surplus);
    %   any deficit (avail<0) is covered from storage in applyNSCCoupling (Stage 3).
    % [M2-B.8b 2026-05-22] 移除 max(C_labile, 0) clamp (applyNSCCoupling 保 ≥0 不破守恒)。
    C_pools_new.C_labile     = C_pools_prev.C_labile + GPP_in_KT - R_a_out_KT - ...
                               (A_C_leaf + A_C_root + A_C_stem + A_C_livecroot + A_C_deadcroot + A_C_fruit + C_to_storage_tot) ...
                               + C_storage_to_labile;   % [labile floor 2026-06-14] + starch->sugar remob toward the sugar floor

    % Organ pools: 加 A_C_X 显示 (turnover/harvest 在外面处理)
    C_pools_new.C_leaf       = C_pools_prev.C_leaf      + A_C_leaf;
    C_pools_new.C_root_fine  = C_pools_prev.C_root_fine + A_C_root;
    % [13c per-pool clamp 2026-06-12] A_C_livestem/A_C_deadstem already set (default split + per-pool
    %   C:N clamp above); use them directly (re-splitting A_C_stem here would discard the per-pool clamp).
    C_pools_new.C_livestem   = C_pools_prev.C_livestem + A_C_livestem;
    C_pools_new.C_deadstem   = C_pools_prev.C_deadstem + A_C_deadstem;
    % [coarse-root 2026-06-12] explicit belowground-wood organ (mirror of livestem/deadstem)
    C_pools_new.C_livecroot  = C_pools_prev.C_livecroot + A_C_livecroot;
    C_pools_new.C_deadcroot  = C_pools_prev.C_deadcroot + A_C_deadcroot;
    C_pools_new.C_fruit      = C_pools_prev.C_fruit     + A_C_fruit;
    % [dual-storage 2026-06-06] root + stem reserve pools (+ 春萌/秋 retransl 在 turnover)
    C_pools_new.C_root_storage = C_pools_prev.C_root_storage + C_to_root_storage - draw_root;   % [labile floor] minus starch->sugar remob
    C_pools_new.C_stem_storage = C_pools_prev.C_stem_storage + C_to_stem_storage - draw_stem;   % [labile floor] minus starch->sugar remob

    % --- Diagnostic ---
    % [mass-balance N reserve 2026-06-11] C reserve refill exposed as a C-flux diagnostic only. The N
    %   reserve is NOT derived from it (former N = C_to_storage / CN_storage_target removed); the N
    %   reserve is mass-balanced from resorption in runPlantCNStep Step 6.
    diag.C_to_root_storage = C_to_root_storage;            % [g C/m2] this-step C refill -> root reserve
    diag.C_to_stem_storage = C_to_stem_storage;            % [g C/m2] this-step C refill -> stem reserve
    % [labile sugar floor 2026-06-14] sugar-floor + starch->sugar remobilisation diagnostics
    diag.C_storage_to_labile = C_storage_to_labile;        % [g C/m2] starch -> sugar remob this step
    diag.C_labile_target     = C_labile_target_g;          % [g C/m2] soluble-sugar floor this step
    diag.C_labile_pre_remob  = C_labile_pre;               % [g C/m2] labile before remob
    diag.A_C_leaf      = A_C_leaf;
    diag.A_C_root      = A_C_root;
    diag.A_C_stem      = A_C_stem;
    diag.A_C_fruit     = A_C_fruit;
    diag.GPP_in_KT     = GPP_in_KT;
    diag.R_a_out_KT    = R_a_out_KT;
    diag.NPP_in_KT     = NPP_in_KT;
    diag.cn_mode       = cnModeCode(cn_mode);              % 0=parallel, 1=c_follows_n, 2=flexible_cn
    % [target-reserve 2026-06-10 STAGE 1: DIAGNOSTIC] wood-based reserve targets (not yet consumed)
    diag.target_stem_storage = target_stem_storage;        % [g C/m2] nsc_frac x stem wood
    diag.target_root_storage = target_root_storage;        % [g C/m2] nsc_frac x coarse-root wood
end

function code = cnModeCode(m)
    switch m
        case 'flexible_cn', code = 2;
        case 'c_follows_n', code = 1;
        otherwise,          code = 0;   % parallel
    end
end

function A_C = clampPoolCN(A_C, C_prev, N_new, CN_max)
%CLAMPPOOLCN  夹紧器官 C 分配, 使分配后池 C:N ≤ CN_max。超出 C 不分配 (留 C_labile, 守恒由 balance 保证)。
%   A_C    [g C/m²]  本步拟分配器官 C
%   C_prev [g C/m²]  上 KT 末器官 C 池
%   N_new  [g N/m²]  本步分配后器官 N 池 (= N_prev + A_N)
%   CN_max [gC/gN]   池 C:N 上界 = CN_target + offset
    eps_n = 1e-30;
    if N_new > eps_n
        if (C_prev + A_C) / N_new > CN_max
            A_C = max(0, CN_max * N_new - C_prev);
        end
    else
        A_C = 0;   % 无 N → 该器官 C:N 无界, 不分配 C (全留 C_labile)
    end
end
