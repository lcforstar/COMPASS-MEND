function Y = buildInitialPlantPools(raw, built)
%BUILDINITIALPLANTPOOLS  Plant C/N pool initial values (P2, carbon-primary).
%
%   Y = init.buildInitialPlantPools(raw, built) builds the plant C/N prognostic
%   pool initial values, faithfully reproducing the plant-pool block (block 4) of
%   the legacy init.setInitialNitrogen. Carbon is primary (from raw.plant /
%   plant_initial_pools.csv); nitrogen is derived as N = C / CN_target using the C:N
%   ratios in built.parameters.plant_cn. Caller stores it at built.Y.plant_cn.
%
%   Order: read C pool init (raw.plant) -> read plant_cn params (P3) -> build C
%   pools (stem live/dead by flivewd; coarse root = croot_stem x stem; storage
%   split by f_storage_root) -> derive N pools (live wood CN_lw, dead wood CN_dw)
%   -> frozen conservation baselines.
%
%   Inputs:
%     raw    struct  raw.plant (plant_initial_pools.csv: C_labile/C_leaf/C_root_fine/
%                    C_stem_wood/C_fruit/C_storage/N_storage); PFT fallback
%     built  struct  built.parameters.plant_cn (cn_ratio/allocation/storage) +
%                    built.parameters.veg.structure.PFT_type
%
%   Returns Y (grouped, all scalar init values [g C or N / m2]):
%     Y.C        labile/leaf/root_fine/livestem/deadstem/livecroot/deadcroot/
%                fruit/root_storage/stem_storage/xsmr  (xsmr = carbon debt, 0)
%     Y.N        same organs (no xsmr -- carbon-only debt)
%     Y.transfer C/N onset transfer pools (leaf/root/livestem), all 0
%     Y.baseline C_pool_init_const / N_pool_init_const (frozen, for conservation)
%
%   Reference: init.setInitialNitrogen block 4 (src/+init/setInitialNitrogen.m:110-291).
% Editor: Jizixing Zhu; Date: 2026-06-27.
% Change: 2026-07-04 Jizixing Zhu -- PR-D prereq: add Y.accum (15 cumulative flux accumulators incl. cum_U_N, init 0) + Y.state (turnover_state/onset_mobilized_flag/resp_night_hours/C_deadstem_soy/N_deadstem_soy/pruned_this_year), the prognostic plant C/N accumulators + one-shot flags runPlantCNStep advances; mirrors setInitialNitrogen.m:203-302. Carried via full-Y snapshot, committed accept-only by the coming scatterNitrogenVariables.
% Change: 2026-07-16 Jizixing Zhu -- document the descriptive plant initialization filename.

    pl  = raw.plant;
    pcn = built.parameters.plant_cn;
    PFT = built.parameters.veg.structure.PFT_type;

    % --- C:N targets + split params (from P3) ---
    CN_L     = pcn.cn_ratio.CN_L;
    CN_root  = pcn.cn_ratio.CN_target_root;
    CN_fruit = pcn.cn_ratio.CN_target_fruit;
    CN_lw    = pcn.cn_ratio.CN_lw;            % live wood (respires)
    CN_dw    = pcn.cn_ratio.CN_dw;            % dead wood (~N-free)
    flivewd    = pcn.allocation.flivewd;      % live-wood fraction (stem live/dead split)
    croot_stem = pcn.allocation.croot_stem;   % coarse-root:stem ratio
    f_sto_root = pcn.storage.f_storage_root;  % storage root/stem split fraction

    % --- C pools: carbon-primary from raw.plant (PFT fallback), >= 0 ---
    C_labile      = max(plantPool(pl, PFT, 'C_labile',    'C_labile_init',    0), 0);
    C_leaf        = max(plantPool(pl, PFT, 'C_leaf',      'C_leaf_init',      0), 0);
    C_root_fine   = max(plantPool(pl, PFT, 'C_root_fine', 'C_root_fine_init', 0), 0);
    C_stem_total  = max(plantPool(pl, PFT, 'C_stem_wood', 'C_stem_wood_init', 0), 0);
    C_fruit       = max(plantPool(pl, PFT, 'C_fruit',     'C_fruit_init',     0), 0);
    C_storage_tot = max(plantPool(pl, PFT, 'C_storage',   'C_storage_init',   0), 0);

    C_livestem  = C_stem_total * flivewd;            % live stem (respires)
    C_deadstem  = C_stem_total * (1 - flivewd);      % dead stem (heartwood)
    C_livecroot = croot_stem * C_livestem;           % coarse-root mirror of stem dual pool
    C_deadcroot = croot_stem * C_deadstem;
    C_root_storage = C_storage_tot * f_sto_root;     % coarse-root reserve
    C_stem_storage = C_storage_tot * (1 - f_sto_root); % stem sapwood reserve
    C_xsmr = 0;                                      % maintenance-respiration debt pool (CLM5 init 0)

    % --- N pools: derived N = C / CN_target; storage N-primary (own init) ---
    N_storage_tot = max(plantPool(pl, PFT, 'N_storage', 'N_storage_init', 0), 0);
    N_labile      = C_labile      / max(CN_L,     1);
    N_leaf        = C_leaf        / max(CN_L,     1);
    N_root_fine   = C_root_fine   / max(CN_root,  1);
    N_fruit       = max(plantPool(pl, PFT, 'N_fruit', '', ...
        C_fruit / max(CN_fruit, 1)), 0);
    N_livestem    = C_livestem    / max(CN_lw,    1);
    N_deadstem    = C_deadstem    / max(CN_dw,    1);
    N_livecroot   = C_livecroot   / max(CN_lw,    1);
    N_deadcroot   = C_deadcroot   / max(CN_dw,    1);
    N_root_storage = N_storage_tot * f_sto_root;
    N_stem_storage = N_storage_tot * (1 - f_sto_root);

    % --- grouped output ---
    Y.C = struct('labile', C_labile, 'leaf', C_leaf, 'root_fine', C_root_fine, ...
        'livestem', C_livestem, 'deadstem', C_deadstem, 'livecroot', C_livecroot, ...
        'deadcroot', C_deadcroot, 'fruit', C_fruit, 'root_storage', C_root_storage, ...
        'stem_storage', C_stem_storage, 'xsmr', C_xsmr);
    Y.N = struct('labile', N_labile, 'leaf', N_leaf, 'root_fine', N_root_fine, ...
        'livestem', N_livestem, 'deadstem', N_deadstem, 'livecroot', N_livecroot, ...
        'deadcroot', N_deadcroot, 'fruit', N_fruit, 'root_storage', N_root_storage, ...
        'stem_storage', N_stem_storage);

    % onset transfer pools (storage->display at budburst); init 0 (legacy: leaf/root/livestem)
    Y.transfer.C = struct('leaf', 0, 'root', 0, 'livestem', 0);
    Y.transfer.N = struct('leaf', 0, 'root', 0, 'livestem', 0);

    % frozen conservation baselines (xsmr=0 included to track the debt excursion)
    Y.baseline.C_pool_init_const = C_labile + C_leaf + C_root_fine + C_livestem + C_deadstem + ...
        C_fruit + C_livecroot + C_deadcroot + C_root_storage + C_stem_storage + C_xsmr;
    Y.baseline.N_pool_init_const = N_labile + N_leaf + N_root_fine + N_livestem + N_deadstem + ...
        N_fruit + N_livecroot + N_deadcroot + N_root_storage + N_stem_storage;

    % cumulative flux accumulators (prognostic; runPlantCNStep + root-uptake glue advance them,
    % committed accept-only; cum_U_N is the plant-N conservation inflow). All init 0. Mirrors
    % setInitialNitrogen.m:203-302. (PR-D: plant C/N + root N uptake migration.)
    Y.accum = struct( ...
        'cum_U_N', 0, 'cum_GPP', 0, 'cum_R_a', 0, 'cum_R_a_potential', 0, ...
        'cum_C_xsmr_writeoff', 0, ...
        'cum_storage_drawdown_count', 0, ...
        'cum_C_harvest', 0, 'cum_N_harvest', 0, ...
        'cum_C_leaf_litter', 0, 'cum_N_leaf_litter', 0, ...
        'cum_C_root_litter', 0, 'cum_N_root_litter', 0, ...
        'cum_C_prune_litter', 0, 'cum_N_prune_litter', 0, ...
        'cum_C_mortality', 0, 'cum_N_mortality', 0);

    % persistent one-shot flags + year-start snapshots (prognostic; committed accept-only). Mirrors
    % setInitialNitrogen.m:209/258-263. Field-by-field build (struct('year_now',[]) would make a
    % scalar struct with an empty field; kept explicit for clarity).
    Y.state.turnover_state       = struct('harvested_this_year', 0, 'year_now', [], 'frPHU_prev', []);
    Y.state.onset_mobilized_flag = false;   % [-] one-shot: storage mobilised this onset
    Y.state.resp_night_hours     = 0;       % [-] consecutive dark hours since sunset
    Y.state.C_deadstem_soy       = C_deadstem;   % [g C/m2] year-start deadstem C (Jan-1 capture)
    Y.state.N_deadstem_soy       = N_deadstem;   % [g N/m2] year-start deadstem N
    Y.state.pruned_this_year     = false;   % [-] one-shot: pruned this autumn
end

% ---------------------------------------------------------------------------
function v = plantPool(pl, PFT, poolName, pftName, default)
% Pool init value: raw.plant (plant_initial_pools.csv) by pool name > PFT table > default.
% Mirrors setInitialNitrogen.getInitVal (cfg.Veg override > PFT > default).
% Editor: Jizixing Zhu; Date: 2026-06-27.
    if istable(pl) && ismember('pool', pl.Properties.VariableNames) ...
            && ismember('value', pl.Properties.VariableNames)
        idx = strcmp(cellstr(string(pl.pool)), poolName);
        if any(idx)
            val = pl.value(find(idx, 1));
            if ~isempty(val) && isnumeric(val) && isfinite(val); v = val; return; end
        end
    end
    if isempty(pftName)
        v = default;
    else
        v = vegetation.getPFTValue(PFT, pftName, default);
    end
end

% Change: 2026-09-02 Jizixing Zhu -- initialize the explicit xsmr debt-
% writeoff ledger used by crop harvest without treating debt as atmospheric C.

% Change: 2026-09-30 Wenjie Liu -- initialize explicit annual standing-leaf
% C/N and stem WSC inventories and event memory only for enabled crop physiology.
% Change: 2026-09-30 Wenjie Liu -- withdraw this round's APSIM stem WSC
%   and standing-senescent-leaf extension; retain the original COMPASS litter
%   and CoLM crop C/N allocation/harvest contracts. Earlier Change records are historical.
