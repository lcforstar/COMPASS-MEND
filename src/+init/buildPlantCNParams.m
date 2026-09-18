function pcn = buildPlantCNParams(raw, built)
%BUILDPLANTCNPARAMS  Plant C/N parameters from the PFT table + site overrides (P3).
%
%   pcn = init.buildPlantCNParams(raw, built) collects the plant carbon/nitrogen
%   parameters into a layered container. It ONLY organizes parameters (no pool
%   initialization -- that is init.buildInitialPlantPools, P2). Caller stores it
%   at built.parameters.plant_cn.
%
%   Priority for each field: site.toml [vegetation] explicit > PFT table.
%   Required plant-C/N fields are strict: a missing PFT-table value raises an
%   initialization error instead of falling through to a hidden code default.
%   Fields that are inactive for a PFT (for example, stem allometry in grass)
%   are set to zero explicitly. A few legacy compatibility constants are kept
%   isolated in compatPftValue so they do not masquerade as PFT-table values.
%
%   Inputs:
%     raw    struct  raw.config.vegetation (PFT_type + any explicit overrides)
%     built  struct  built.parameters.veg.structure.PFT_type
%
%   Returns pcn (layered):
%     .coupling     cn_coupling_mode, cn_flex_offset
%     .cn_ratio     CN_L, CN_target_root/stem/fruit, CN_lw, CN_dw
%     .allocation   f_alloc_leaf/root/stem/fruit, croot_stem, froot_leaf, stem_leaf, fcur
%     .respiration  r_m10, R_g_fraction
%     .turnover     leaf_long_days, ndays_on, ndays_off, lwtop, harvest_fraction, retransl_fraction
%     .storage      has_storage_pool, fstor2tran, f_storage_root, storage_cushion
%
%   Field names verified against src/+vegetation/getPFTParameterSet.m. Reference:
%   init.setInitialNitrogen (which reads the same PFT fields for pool init).
% Editor: Jizixing Zhu; Date: 2026-06-27.
% Change: 2026-06-28 Jizixing Zhu — f_alloc_* initialized by PFT class: allometry PFTs
%   without a storage pool (grass) carry the dead NaN sentinel; crops stay required.

    v   = raw.config.vegetation;
    PFT = built.parameters.veg.structure.PFT_type;

    % --- coupling (site.toml [switches]/[vegetation]) ---
    pcn.coupling.cn_coupling_mode = char(getf(v, 'cn_coupling_mode', 'flexible_cn'));
    pcn.coupling.cn_flex_offset   = getf(v, 'cn_flex_offset', 10.0);

    hasFruit       = logical(round(requiredPftValue(v, PFT, 'has_fruit')));
    hasStoragePool = logical(round(requiredPftValue(v, PFT, 'has_storage_pool')));

    % --- organ C:N targets (PFT) ---
    pcn.cn_ratio.CN_L            = requiredPftValue(v, PFT, 'CN_L');
    pcn.cn_ratio.CN_target_root  = requiredPftValue(v, PFT, 'CN_target_root');
    pcn.cn_ratio.CN_target_stem  = requiredPftValue(v, PFT, 'CN_target_stem');
    pcn.cn_ratio.CN_target_fruit = activeOrZero(v, PFT, 'CN_target_fruit', hasFruit);
    pcn.cn_ratio.CN_lw           = requiredPftValue(v, PFT, 'CN_lw');
    pcn.cn_ratio.CN_dw           = requiredPftValue(v, PFT, 'CN_dw');

    % --- carbon allocation (PFT) ---
    % The allocation engine differs by PFT class, so f_alloc_* is initialized per class:
    %   * c_allometry PFTs (froot_leaf/stem_leaf > 0): NPP is split via the allometry
    %     ratios in computeCarbonAllocation / computeOrganDemand. f_alloc_* is then read
    %     ONLY by the spring onset-xfer storage split (runPlantCNStep), which runs for
    %     storage PFTs. So an allometry PFT WITHOUT a storage pool (grass 12-14) never
    %     reads f_alloc_*: it is carried as the dead NaN sentinel set in
    %     getPFTParameterSet, not required-finite (which would reject NaN).
    %   * non-allometry PFTs (crops 15-22): f_alloc_* IS the active fixed-fraction
    %     split and stays strictly required.
    pcn.allocation.croot_stem    = optionalPftValue(v, PFT, 'croot_stem', 0);
    pcn.allocation.froot_leaf    = optionalPftValue(v, PFT, 'froot_leaf', 0);
    pcn.allocation.stem_leaf     = optionalPftValue(v, PFT, 'stem_leaf',  0);
    usesAllometry = (pcn.allocation.froot_leaf > 0) || (pcn.allocation.stem_leaf > 0);
    fAllocIsDead  = usesAllometry && ~hasStoragePool;   % grass: allometry + no onset-xfer
    if fAllocIsDead
        pcn.allocation.f_alloc_leaf = deadPftValue(v, PFT, 'f_alloc_leaf');
        pcn.allocation.f_alloc_root = deadPftValue(v, PFT, 'f_alloc_root');
        pcn.allocation.f_alloc_stem = deadPftValue(v, PFT, 'f_alloc_stem');
    else
        pcn.allocation.f_alloc_leaf = requiredPftValue(v, PFT, 'f_alloc_leaf');
        pcn.allocation.f_alloc_root = requiredPftValue(v, PFT, 'f_alloc_root');
        pcn.allocation.f_alloc_stem = requiredPftValue(v, PFT, 'f_alloc_stem');
    end
    pcn.allocation.f_alloc_fruit = activeOrZero(v, PFT, 'f_alloc_fruit', hasFruit);
    pcn.allocation.flivewd       = requiredPftValue(v, PFT, 'flivewd');   % live-wood fraction (live/dead stem split)
    pcn.allocation.fcur          = compatPftValue(v, PFT, 'fcur', 0.95);

    % --- respiration (PFT; Sitch 2003 + CLM5 grperc) ---
    pcn.respiration.r_m10        = requiredPftValue(v, PFT, 'r_m10');
    pcn.respiration.R_g_fraction = requiredPftValue(v, PFT, 'R_g_fraction');

    % --- turnover / phenology ramps (PFT) ---
    pcn.turnover.leaf_long_days    = requiredPftValue(v, PFT, 'leaf_long_days');
    pcn.turnover.ndays_on          = requiredPftValue(v, PFT, 'ndays_on');
    pcn.turnover.ndays_off         = requiredPftValue(v, PFT, 'ndays_off');
    pcn.turnover.lwtop             = requiredPftValue(v, PFT, 'lwtop');
    pcn.turnover.harvest_fraction  = activeOrZero(v, PFT, 'harvest_fraction', hasFruit);
    pcn.turnover.retransl_fraction = requiredPftValue(v, PFT, 'retransl_fraction');

    % --- storage / NSC reserve (PFT + site.toml) ---
    pcn.storage.has_storage_pool = double(hasStoragePool);
    pcn.storage.fstor2tran       = requiredPftValue(v, PFT, 'fstor2tran');
    pcn.storage.f_storage_root   = optionalPftValue(v, PFT, 'f_storage_root', 0.55);
    pcn.storage.storage_cushion  = compatPftValue(v, PFT, 'storage_cushion', 1.2);
end

% ---------------------------------------------------------------------------
function val = getf(s, name, default)
% site.toml [vegetation] field by name (unset/empty -> default).
% Editor: Jizixing Zhu; Date: 2026-06-27.
    val = default;
    if isfield(s, name) && ~isempty(s.(name)); val = s.(name); end
end

% ---------------------------------------------------------------------------
function val = requiredPftValue(rawVeg, PFT, name)
% Numeric required param: site.toml [vegetation] explicit > strict PFT table.
% Editor: Jizixing Zhu; Date: 2026-06-27.
    if hasNumericSiteOverride(rawVeg, name)
        val = rawVeg.(name);
        return;
    end

    try
        val = vegetation.getPFTValue(PFT, name);
    catch ME
        error('buildPlantCNParams:missingRequiredPFT', ...
            'PFT "%s" must define required plant C/N parameter "%s" or site.toml [vegetation] must override it. Original: %s', ...
            char(PFT), name, ME.message);
    end
    validateScalarFinite(val, PFT, name);
end

% ---------------------------------------------------------------------------
function val = optionalPftValue(rawVeg, PFT, name, inactiveDefault)
% Numeric optional param: site override > PFT table > explicit inactive value.
% Editor: Jizixing Zhu; Date: 2026-06-27.
    if hasNumericSiteOverride(rawVeg, name)
        val = rawVeg.(name);
    else
        val = vegetation.getPFTValue(PFT, name, inactiveDefault);
    end
    validateScalarFinite(val, PFT, name);
end

% ---------------------------------------------------------------------------
function val = compatPftValue(rawVeg, PFT, name, legacyDefault)
% Numeric compatibility param for legacy fields not yet complete in the PFT table.
% Editor: Jizixing Zhu; Date: 2026-06-27.
    if hasNumericSiteOverride(rawVeg, name)
        val = rawVeg.(name);
    else
        val = vegetation.getPFTValue(PFT, name, legacyDefault);
    end
    validateScalarFinite(val, PFT, name);
end

% ---------------------------------------------------------------------------
function val = deadPftValue(rawVeg, PFT, name)
% Allocation param that this PFT never reads (c_allometry PFT with no storage pool,
% e.g. grass): honour a site.toml [vegetation] override if present, otherwise pass the
% PFT-table value through WITHOUT the finite check -- the NaN dead-value sentinel set in
% getPFTParameterSet is expected here. Still requires a numeric scalar.
% Editor: Jizixing Zhu; Date: 2026-06-28.
    if hasNumericSiteOverride(rawVeg, name)
        val = rawVeg.(name);
    else
        val = vegetation.getPFTValue(PFT, name, NaN);
    end
    if ~(isnumeric(val) && isscalar(val))
        error('buildPlantCNParams:invalidPFTValue', ...
            'PFT "%s" parameter "%s" must be a numeric scalar.', char(PFT), name);
    end
end

% ---------------------------------------------------------------------------
function val = activeOrZero(rawVeg, PFT, name, isActive)
% Numeric param required only when the corresponding PFT process is active.
% Editor: Jizixing Zhu; Date: 2026-06-27.
    if isActive
        val = requiredPftValue(rawVeg, PFT, name);
    else
        val = optionalPftValue(rawVeg, PFT, name, 0);
    end
end

% ---------------------------------------------------------------------------
function tf = hasNumericSiteOverride(rawVeg, name)
% True when site.toml [vegetation] provides a usable scalar numeric override.
% Editor: Jizixing Zhu; Date: 2026-06-27.
    tf = isfield(rawVeg, name) && ~isempty(rawVeg.(name)) ...
        && isnumeric(rawVeg.(name)) && isscalar(rawVeg.(name)) && isfinite(rawVeg.(name));
end

% ---------------------------------------------------------------------------
function validateScalarFinite(val, PFT, name)
% Fail early when a PFT value is not a finite scalar numeric parameter.
% Editor: Jizixing Zhu; Date: 2026-06-27.
    if ~(isnumeric(val) && isscalar(val) && isfinite(val))
        error('buildPlantCNParams:invalidPFTValue', ...
            'PFT "%s" parameter "%s" must be a finite scalar numeric value.', char(PFT), name);
    end
end
