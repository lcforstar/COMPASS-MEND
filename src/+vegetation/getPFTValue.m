function value = getPFTValue(pft_name, fieldname, default_value)
% getPFTValue
% Retrieves a specific photosynthetic parameter from a given PFT.
%
% This utility accesses the internal PFT parameter table and extracts 
% a single field value for a specified PFT name.
%
% Inputs:
%   pft_name  - string, e.g. 'C3 grass', 'BET Tropical', etc.
%   fieldname - string, one of the field names listed below.
%
% Output:
%   value     - scalar value from the selected PFT and field.
%
% Example:
%   SLA = getPFTValue('C3 grass', 'SLA0');
%
% -------------------------------------------------------------------------
% Supported fields (all drawn from Tables 8.1 & 8.2):
%
%   'PFT'           - PFT short name (text)
%   'FullName'      - Full descriptive name (text)
%
%   'm'             - Ball–Berry stomatal slope [unitless]
%   'alpha'         - Quantum yield of electron transport [mol e⁻ mol⁻¹ photon]
%   'CN_L'          - Leaf carbon-to-nitrogen ratio [g C g⁻¹ N]
%
%   'SLA0'          - Specific leaf area at canopy top (x = 0) [m² leaf g⁻¹ C]
%   'SLAm'          - Linear slope of SLA with cumulative LAI [m² leaf g⁻¹ C]
%
%
%   'Vcmax25'       - Vcmax at 25 °C for top-leaf (no N-limit) [µmol CO₂ m⁻² s⁻¹]
%   'Ztop','Zbot'   - 冠层高度上/下界 [m]
%
% ------ Phase C M2-B.0 (2026-05-18) — 25 新增字段 (详 cn_coupling §12.13/12.14) ------
%
%   Group A — PFT 物种属性:
%     'is_perennial'        [0/1]   多年生开关 (annual=0)
%     'has_fruit'           [0/1]   有果实/种子收获池
%     'has_storage_pool'    [0/1]   做 N_storage 季节储备 (多年生落叶)
%     'photo_type'          'C3'/'C4'
%
%   Group B — 器官 C:N 目标 + C 分配比 (M2-B.1):
%     'CN_target_root'      [-]     细根 C:N (Yuan 2011, Aerts 1996)
%     'CN_target_stem'      [-]     木质部 C:N (Yang 2014; 草本=0)
%     'CN_target_fruit'     [-]     果实/种子 C:N (Mao 2018; 无收获=0)
%     'f_alloc_leaf'        [-]     C 分配到叶, 4 个 alloc 和 = 1 (Sitch 2003 LPJ)
%     'f_alloc_root'        [-]     C 分配到细根
%     'f_alloc_stem'        [-]     C 分配到木质部 (草本=0)
%     'f_alloc_fruit'       [-]     C 分配到果实 (生长期物候触发)
%       NOTE (2026-06-28): f_alloc_* is the LEGACY fixed-fraction allocation path.
%       It is ACTIVE only for PFTs without c_allometry ratios -- currently the crop
%       PFTs 15-22 (placeholder; pending a proper crop-phenology module). For PFTs
%       that set froot_leaf/stem_leaf (trees/shrubs 1-11,24 and grass 12-14) the
%       allocation goes through the CLM5 c_allometry branch instead, so f_alloc_*
%       there is a DEAD fallback (not read by computeCarbonAllocation/computeOrganDemand;
%       runPlantCNStep onset-mobilization still reads it -- see that file). Do NOT
%       read f_alloc_* as the allocation for an allometry PFT.
%
%   Group C — LUNA (Ali 2016, M2-B.3):
%     'Jmaxb1'              [-]    LUNA light response (CLM5 namelist default 0.093563; Apple [site] 0.1745)
%                                  Only LUNA field consumed by physics — see COMPASS.m:864.
%     % [REMOVED 2026-05-26] Rd_fraction / Jmaxb0 / N_struct_per_LMA / N_photo_ref /
%     % tc_j0 / H_LUNA / NUE_Vcmax25 / NUE_Jmax25 were dead fields — never read by
%     % physics. LUNA uses hardcoded CLM5 constants in
%     % src/+vegetation/+integrated/+LUNA/lunaParams.m (Jmaxb0=0.0311, SNC=0.004,
%     % Wc2Wjb0=0.8054, relhExp=6.0999). Rd_fraction had no consumer at all.
%
%   Group D — Nacry 2013 根膜参数 (从 cfg.Veg 迁移到 PFT 表):
%     'Vmax_NH4_umol_g_h'   [μmol N/(g-root·h)] AMT max (Nacry 2013)
%     'Vmax_NO3_umol_g_h'   [μmol N/(g-root·h)] NRT2 max
%     'Km_NH4'              [mg/L]              AMT 半饱和
%     'Km_NO3'              [mg/L]              NRT2 半饱和
%     'Cmin_NH4'            [mg/L]              根膜 NH4 floor
%     'Cmin_NO3'            [mg/L]              根膜 NO3 floor
%     'SRL'                 [m-root/g-root]     Specific root length (Wells 2002)
%     'Vmax_NO3_LATS_umol_g_h' [μmol N/g/h]     LATS Vmax (M3, 暂留 0)
%     'Km_LATS_NO3'         [mg/L]              LATS Km ≈ 70 (Tsay 2007)
%
%   Group E — 叶寿命 + retranslocation (M2-B.4):
%     'leaf_long_days'      [d]   叶寿命 (Niinemets 2010)
%     'retransl_fraction'   [-]   落叶 N 回收比例 (Aerts 1996: 阔叶 0.7, 针叶 0.5)
% -------------------------------------------------------------------------

     % Load parameter set
    PFTs = vegetation.getPFTParameterSet();

    % Find the index of the specified PFT
    idx = find(strcmpi({PFTs.PFT}, pft_name), 1);
    if isempty(idx)
        error('PFT "%s" not found.', pft_name);
    end

    % Extract the structure for the specified PFT
    pft = PFTs(idx);

    % Check whether the field exists AND is non-empty. A struct-array field set on only SOME PFTs
    % (e.g. croot_stem_alloc/stem_leaf on woody PFTs, use_fruittree_alloc on apple) EXISTS as [] on
    % every other PFT; treat that empty like a missing field so PFT-specific params fall back to the
    % caller's default (0) instead of returning [] (which then breaks scalar logic like denom = 1+a*[]).
    if ~isfield(pft, fieldname) || isempty(pft.(fieldname))
        % [F2/13a] Optional 3rd arg: return a default for PFT-specific fields that only some PFTs
        % define. Without it, the missing/empty-field error stands (strict default for required fields).
        if nargin >= 3
            value = default_value;
            return;
        end
        error('Field "%s" not found (or empty) in PFT structure.', fieldname);
    end

    % Return the specific value
    value = pft.(fieldname);
end
