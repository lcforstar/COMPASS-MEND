function [pools_new, diag] = allocatePlantN(pools_prev, U_N_total_g_m2_s, D_N, dt, pftParams)
%ALLOCATEPLANTN  Update plant N pools (N_labile + 4 organ pools) per Δt, C-dominant + C:N elasticity.
%
%   [v12 2026-05-29 C 主导 + C:N 弹性] 第 5 参 pftParams (CN_target_* + cn_flex_offset)。
%     架构: C 主导 (computeOrganDemand 已 D_N = NPP·f_alloc / CN_target, 即 C 先定后反推 N 需求)。
%     N 不足时不立即砍 C, 而是让器官 C:N 朝 CN_max=CN_target+offset 上浮 (C 照长, 每单位 C 少配 N),
%     到 CN_max 仍不足才由 computeCarbonAllocation 砍 C。对标 CLM5 FlexibleCN nscalar
%     (NutrientCompetitionFlexibleCNMod.F90:1660-1667: leafcn 在 [target±10] 内浮动)。
%     diag 增加 .actual_CN_leaf/root/stem/fruit (供 C 模块 A_C=A_N×actual_CN 用)。
%
%   Returns:
%     pools_new  struct  Updated pools after dt [g N/m²]:
%                  .N_labile        single mobile pool
%                  .N_leaf          leaf N
%                  .N_root_fine     fine root N
%                  .N_livestem .N_deadstem  [13c] live/dead wood N (split by flivewd)
%                  .N_fruit         fruit N
%     diag       struct  Diagnostic [g N/(m²·s)]: .A_leaf, .A_root, .A_stem, .A_fruit
%                                                  .f_N_growth [-]
%                                                  .N_labile_supply [g N/m²] = N_labile_prev + U_N×dt
%
%   Inputs:
%     pools_prev  struct  Pools at start of Δt (5 prognostic fields, see Returns)
%     U_N_total_g_m2_s    [g N/(m²·s)]  Total N uptake from soil (NH4+NO3) via
%                                       SinkS depth-integrated, converted to g/m²/s
%     D_N         struct  From computeOrganDemand: .leaf, .root, .stem, .fruit, .total
%                          all [g N/(m²·s)]
%     dt          [s]     Outer time step (typically 3600 = 1 h)
%
%   Method (LPJ Sitch 2003 Eq. 41-43 + N balance):
%
%     1. N_labile supply:
%        N_labile_supply = N_labile_prev + U_N_total × dt  [g N/m²]
%
%     2. N limitation factor (你 md §8.7 Eq. 41):
%        f_N_growth = clamp(N_labile_supply / (D_N_total × dt + ε), 0, 1)
%
%        f_N_growth = 1: N_labile足够, 完整 demand 实现
%        f_N_growth < 1: N_labile不足, 各器官按比例缩减
%
%     3. Actual allocation per organ:
%        A_i = D_N,i × f_N_growth × dt        [g N/m²]
%        Note: A_i is integrated over dt (Δstate), not rate.
%
%     4. Pool updates (state = state + Δstate):
%        N_leaf_new  = N_leaf_prev  + A_leaf
%        N_root_new  = N_root_prev  + A_root
%        N_stem_new  = N_stem_prev  + A_stem
%        N_fruit_new = N_fruit_prev + A_fruit
%        N_labile_new = N_labile_supply - (A_leaf + A_root + A_stem + A_fruit)
%
%     5. Mass balance (debug check, not enforced):
%        ΔN_total = (N_organ_new - N_organ_prev) + (N_labile_new - N_labile_prev)
%                 = U_N_total × dt    (only input source in M2-B.1b)
%        M2-B.4 will add retransl and seasonal storage transfers.
%        M2-B.5 will add harvest output.
%
%   Notes:
%     - Pure function: pools_prev not mutated, returns pools_new.
%     - When D_N_total = 0 (winter, NPP=0): all A_i = 0, N_labile only accumulates U_N.
%     - When U_N = 0 AND N_labile_prev = 0: f_N_growth = 0, all A_i = 0 (starvation).
%
%   Reference:
%     Sitch S et al. (2003) GCB 9:161-185 Eq. 41-43 (N-limited growth scaling).
%     CLM5 Tech Note Ch.14 (allocation algorithm reference).
%     COMPASS Phase C M2-B.1b plan: cn_coupling_path_and_literature_audit.md §12.17.3.
%
%   Phase C M2-B.1b (2026-05-19): new module. Called by COMPASS.m AFTER root
%   uptake (SinkS → U_N_total) and computeOrganDemand. ClimaLand alignment.

% --- Constants ---
eps_safe = 1e-30;

% --- Step 1: N_labile supply (上 KT 末值 + 本 KT 根吸 × dt) ---
N_labile_supply = pools_prev.N_labile + U_N_total_g_m2_s * dt;

% --- Step 2: f_N_growth limitation factor (Eq. 41) ---
demand_integrated = D_N.total * dt;       % [g N/m²] (per dt)
if demand_integrated > eps_safe
    f_N_growth = max(0, min(1, N_labile_supply / (demand_integrated + eps_safe)));
else
    f_N_growth = 0;     % NPP=0 (winter) → no demand → no allocation
end

% --- Step 3: Actual organ N allocation [g N/m²] ---
% N 守恒不变: A_N = D_N × f_N_growth × dt (供给不足时等比例缩 N, Σ A_N ≤ N_labile_supply)。
A_leaf  = D_N.leaf  * f_N_growth * dt;
A_root  = D_N.root  * f_N_growth * dt;
A_stem  = D_N.stem  * f_N_growth * dt;
if isfield(D_N, 'croot'), A_croot = D_N.croot * f_N_growth * dt; else, A_croot = 0; end   % [coarse-root 2026-06-12]
if isfield(D_N, 'livestem')
    A_livestem  = D_N.livestem  * f_N_growth * dt;
    A_deadstem  = D_N.deadstem  * f_N_growth * dt;
    A_livecroot = D_N.livecroot * f_N_growth * dt;
    A_deadcroot = D_N.deadcroot * f_N_growth * dt;
else
    % Backward-compatible fallback for old D_N structs.
    A_livestem  = A_stem  * pftParams.flivewd;
    A_deadstem  = A_stem  * (1 - pftParams.flivewd);
    A_livecroot = A_croot * pftParams.flivewd;
    A_deadcroot = A_croot * (1 - pftParams.flivewd);
end
A_fruit = D_N.fruit * f_N_growth * dt;

% --- Step 3b: C:N 弹性 (C 主导, CLM5 nscalar) ---
% N 不足 (f_N_growth<1) 时, 器官 C:N 上浮 = CN_target/f_N_growth, clamp 到 [CN_target, CN_max]。
% C 模块用 actual_CN 算 A_C=A_N×actual_CN: f_N 小→actual_CN 大→C 照长 (少配 N/单位 C);
% 撞 CN_max 封顶后 C 不再随 f_N 增长, 余 C 由 computeCarbonAllocation 砍入 C_labile。
% 对标 CLM5 FlexibleCN:1660-1667 (actual_leafcn 在 [target±10] clamp + nscalar)。
offset = pftParams.cn_flex_offset;       % [gC/gN] 弹性半幅 (默认 10, config)
if f_N_growth > eps_safe
    cn_scale = 1 / f_N_growth;           % f_N<1 → >1, C:N 上浮倍率
else
    cn_scale = 1;                        % NPP=0/无供给, A_N=0, actual_CN 无意义, 置 target
end
actual_CN_leaf  = min(pftParams.CN_L              * cn_scale, pftParams.CN_L              + offset);
actual_CN_root  = min(pftParams.CN_target_root    * cn_scale, pftParams.CN_target_root    + offset);
actual_CN_stem  = min(pftParams.CN_target_stem    * cn_scale, pftParams.CN_target_stem    + offset);
actual_CN_fruit = min(pftParams.CN_target_fruit   * cn_scale, pftParams.CN_target_fruit   + offset);
% annual PFT: CN_target_stem=0 → actual_CN_stem=0 (C 模块跳 stem)
if pftParams.CN_target_stem <= 0, actual_CN_stem = 0; end

% --- Step 3c: REMOVED fcur split [N-side, 2026-06-11] ---
% The old (1-fcur) N storage split + the display fcur-skim are REMOVED: display organs receive the
%   FULL N allocation. The N RESERVE is updated in runPlantCNStep Step 6 by MASS BALANCE from real
%   N fluxes (leaf/live-wood resorption in, spring remobilisation out), NOT from a fixed storage C:N.
%   (The former C-derived N = C_to_storage / CN_storage_target = 20 had no source; see the note in
%   runPlantCNStep Step 6.)
N_to_root_storage = 0;   % [g N/m²] reserve N now updated by mass balance in runPlantCNStep, not here
N_to_stem_storage = 0;

% --- Step 4: Pool updates (state += Δstate) ---
pools_new.N_leaf      = pools_prev.N_leaf      + A_leaf;
pools_new.N_root_fine = pools_prev.N_root_fine + A_root;
% Live/dead woody N increments are computed with CN_lw/CN_dw demands upstream.
pools_new.N_livestem  = pools_prev.N_livestem  + A_livestem;
pools_new.N_deadstem  = pools_prev.N_deadstem  + A_deadstem;
pools_new.N_livecroot = pools_prev.N_livecroot + A_livecroot;
pools_new.N_deadcroot = pools_prev.N_deadcroot + A_deadcroot;
pools_new.N_fruit     = pools_prev.N_fruit     + A_fruit;
% N_labile = supply - display (full). Reserve N is subtracted later in runPlantCNStep (when the
%   C-side C_to_storage is known), preserving total-N conservation.
pools_new.N_labile    = N_labile_supply - (A_leaf + A_root + A_stem + A_croot + A_fruit);

% Clamp: numerical safety (shouldn't trigger if mass balance correct)
pools_new.N_labile = max(pools_new.N_labile, 0);

% --- Diagnostic ---
% A_i 保留 per-KT 增量 [g N/m²]，与 pool 单位一致，对账 Σ A_i ≡ Δ(pools)
diag.A_leaf          = A_leaf;       % [g N/m²] (per dt)
diag.A_root          = A_root;
diag.A_stem          = A_stem;
diag.A_croot         = A_croot;   % [coarse-root 2026-06-12] g N/m² (per dt)
diag.A_livestem      = A_livestem;
diag.A_deadstem      = A_deadstem;
diag.A_livecroot     = A_livecroot;
diag.A_deadcroot     = A_deadcroot;
diag.A_fruit         = A_fruit;
diag.f_N_growth      = f_N_growth;
diag.N_labile_supply = N_labile_supply;
diag.N_to_root_storage = N_to_root_storage;   % [g N/m²] [dual-storage] (1−fcur)→root reserve
diag.N_to_stem_storage = N_to_stem_storage;   % [g N/m²] [dual-storage] (1−fcur)→stem reserve, runPlantCNStep 累加
% actual_CN [gC/gN] 供 computeCarbonAllocation 算 A_C=A_N×actual_CN (C 主导 + 弹性)
diag.actual_CN_leaf  = actual_CN_leaf;
diag.actual_CN_root  = actual_CN_root;
diag.actual_CN_stem  = actual_CN_stem;
diag.actual_CN_fruit = actual_CN_fruit;

end
