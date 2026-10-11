function veg = buildVeg(raw, built)
% Change: 2026-09-30 Wenjie Liu -- selectively restore COMPASS canopy
% behavior; retain winter-wheat phenology and CoLM C/N allocation.
%BUILDVEG  Build the static vegetation parameters from the [vegetation] config.
%
%   veg = init.buildVeg(raw, built) maps the friendly site.toml [vegetation]
%   names onto the legacy VegetationData parameter names the model consumes,
%   with enum/bool conversions, grouped into process sub-structs, and (when built
%   is supplied) the per-layer root distribution RAI. Caller stores it at
%   built.parameters.veg.
%
%   Stages: V1 static parameters (always); V2 root distribution RAI (needs
%   built.domain + raw.root). The phenology driver + LAI/hc/SAI state are later
%   stages (V3-V4) in init.buildVegState (-> built.vegetation_state).
%
%   Inputs:
%     raw    struct  io.loadSite output (uses .config.vegetation, .root)
%     built  struct  needs .domain (dz/NN/tot_depth) for V2 RAI. Optional: if
%                    omitted, only V1 (static params) is built.
%
%   Returns veg with sub-structs (legacy field names):
%     veg.structure      aerodynamic / Wallace-Verhoef canopy params
%     veg.canopy         LAI/hc/SAI modes + caps + vertical profile + interception
%     veg.phenology      GDD phenology params (model LAI)
%     veg.hydraulic      plant hydraulic scheme + dynamic-scheme + tree fallback
%     veg.root           root params (scalars) + RAI [NL x 1, bottom-up] (V2)
%     veg.photosynthesis leaf photosynthesis + stomatal (Medlyn) + LUNA
%
%   NOT here (per the buildVeg I/O contract): hydraulic_init_psi* and C-pool init
%   (-> buildInitialPlantPools); cn_coupling_mode/cn_flex_offset (-> plant_cn);
%   [scope]/[leaf_optical] (-> buildScope). Enum maps verified against
%   legacy VegetationData (Photo_Vcmax_flag: manual=2/pft_table=0; lLUNA 0/1;
%   Photo_LeafResp_Method: vcmax_proportional=0/atkin2015=2).
% Editor: Jizixing Zhu; Date: 2026-06-27.
% Change: 2026-07-01 Jizixing Zhu — restore Init_Plant validation gates (psil_min < 0, root_dist_mode in {0,1}, root-contact thresholds).
% Change: 2026-07-03 Jizixing Zhu — read veg.canopy.LIDFa/LIDFb from config (leaf inclination distribution for RTMo/RTMt), default SCOPE spherical -0.35/-0.15.
% Change: 2026-07-03 Jizixing Zhu — read veg.structure.CropCode from config (crop_code, default 1 = legacy VegetationData); feeds ebal resistances_clm45 leaf-width lookup (was hardwired 1 in runStepEbal).
% Change: 2026-07-05 Jizixing Zhu — add veg.canopy.Ps_max_diurnal (env PS_MAX via io.getEnvDouble) + Ps_min_diurnal (moved from runtime.canopy_water; state-ownership migration increment 3 stage 1).
% Change: 2026-07-16 Jizixing Zhu — resolve foliage clumping index from the selected PFT or an explicit site value and validate 0 < Omega <= 1.
% Change: 2026-07-20 Jizixing Zhu — add optional heat_gs_* under photosynthesis (default enable=false).
% Change: 2026-09-30 Wenjie Liu -- explicit annual-crop canopy/harvest policy;
% fail fast for unsupported prognostic height/SAI rather than constant proxies.

    v = raw.config.vegetation;
    veg = struct();

    % --- structure: aerodynamics / Wallace-Verhoef (friendly -> legacy) ---
    s = struct();
    s.PFT_type    = char(getf(v, 'PFT_type', 'Apple'));
    s.Method      = char(getf(v, 'method', 'verhoef'));
    s.SurfaceType = char(getf(v, 'surface_type', 'Grass'));
    s.CropCode    = getf(v, 'crop_code', 1);              % [-] leaf-width lookup for vegetation.energy.resistancesClm45(1=static leaf_width); legacy VegetationData default 1
    s.w_leaf_max  = getf(v, 'leaf_width', 0.05);          % friendly leaf_width -> w_leaf_max
    s.Cd     = getf(v, 'Cd', 0.2);
    s.rwc    = getf(v, 'rwc', 0.0);
    s.CR     = getf(v, 'CR', 0.35);
    s.CD1    = getf(v, 'CD1', 20.6);
    s.Psicor = getf(v, 'psicor', 0.2);                    % friendly psicor -> Psicor
    s.CSSOIL = getf(v, 'cssoil', 0.01);                   % friendly cssoil -> CSSOIL
    veg.structure = s;

    % --- canopy: LAI/hc/SAI modes + caps + vertical profile + interception ---
    c = struct();
    c.LAI_mode      = char(getf(v, 'LAI_mode', 'constant'));
    c.LAI_engine    = char(getf(v, 'LAI_engine', 'leafc_sla_phen'));
    c.hc_mode       = char(getf(v, 'hc_mode', 'constant'));
    c.SAI_mode      = char(getf(v, 'SAI_mode', 'constant'));
    c.LAI_max       = getf(v, 'LAI_max', 5.0);
    c.laimx         = getf(v, 'laimx', 0);
    c.LAI_constant  = getf(v, 'LAI_constant', 3.0);
    c.canopy_height = getf(v, 'canopy_height', 0.5);
    c.SAI_constant  = getf(v, 'SAI_constant', 0.5);
    c.LAI_profile_mode  = char(getf(v, 'LAI_profile_mode', 'uniform'));
    c.LAI_profile_peak  = getf(v, 'LAI_profile_peak', 0.50);
    c.LAI_profile_sigma = getf(v, 'LAI_profile_sigma', 0.24);
    c.LAI_profile_skew  = getf(v, 'LAI_profile_skew', 0.8);
    c.LAI_profile_k     = getf(v, 'LAI_profile_k', 1.5);
    c.vegetation_data_file = char(getf(v, 'vegetation_data_file', ''));
    c.canopy_storage_per_LAI = getf(v, 'canopy_storage_per_LAI', 0.1);
    c.interc_F0              = getf(v, 'interc_F0', 0.15);
    c.canopy_water_scheme    = char(getf(v, 'canopy_water_scheme', 'single'));
    % canopy-water diurnal reservoir split weights (moved from runtime.canopy_water 2026-07-05).
    % Ps_max via env PS_MAX (io.getEnvDouble: COMPASS_PS_MAX > STEMMUS_PS_MAX > 0.7); Ps_min fixed 0.
    c.Ps_max_diurnal         = io.getEnvDouble('PS_MAX', 0.7);   % [-] sunlit reservoir split weight
    c.Ps_min_diurnal         = 0.0;                              % [-] shaded reservoir split floor
    % canopy-layer / canopy-air feedback switches (vegetation-canopy options)
    c.canopy_layer_profiles  = toFlag(getf(v, 'canopy_layer_profiles', false));
    c.canopy_air_profile_eb  = toFlag(getf(v, 'canopy_air_profile_eb', false));
    c.canopy_air_ta_feedback = toFlag(getf(v, 'canopy_air_ta_feedback', true));
    c.canopy_air_ea_feedback = toFlag(getf(v, 'canopy_air_ea_feedback', false));
    c.canopy_air_co2_feedback = toFlag(getf(v, 'canopy_air_co2_feedback', false));
    % Leaf inclination distribution (Verhoef LIDF two-parameter form, |LIDFa|+|LIDFb|<=1).
    % RTMo/RTMt read these via veg.canopy to build lidf (equations.leafangles). Config-driven
    % so each site sets its own leaf-angle distribution; default = SCOPE spherical (-0.35/-0.15).
    c.LIDFa = getf(v, 'LIDFa', -0.35);
    c.LIDFb = getf(v, 'LIDFb', -0.15);
    c.clumping_index_source = lower(char(getf(v, 'clumping_index_source', 'pft_table')));
    pftClumpingIndex = vegetation.getPFTValue(s.PFT_type, 'clumping_index', 1.0);
    switch c.clumping_index_source
        case 'pft_table'
            c.clumping_index = pftClumpingIndex;
        case 'manual'
            c.clumping_index = getf(v, 'clumping_index', pftClumpingIndex);
        otherwise
            error('buildVeg:badClumpingIndexSource', ...
                'clumping_index_source "%s" must be pft_table or manual.', ...
                c.clumping_index_source);
    end
    if ~isscalar(c.clumping_index) || ~isfinite(c.clumping_index) || ...
            c.clumping_index <= 0 || c.clumping_index > 1
        error('buildVeg:badClumpingIndex', ...
            'clumping_index must be a finite scalar in (0,1], got %s.', ...
            mat2str(c.clumping_index));
    end
    veg.canopy = c;

    % --- phenology: GDD model params (model LAI) ---
    ph = struct();
    ph.T_BASE     = getf(v, 'T_BASE', 5);
    ph.T_OPT      = getf(v, 'T_OPT', 25);
    ph.PHU        = getf(v, 'PHU', 4555);
    ph.GMHU       = getf(v, 'GMHU', 300);
    ph.sowing_day = getf(v, 'sowing_day', 104);
    ph.FRGW1    = getf(v, 'FRGW1', 0.1);
    ph.LAIMX1   = getf(v, 'LAIMX1', 0.15);
    ph.FRGW2    = getf(v, 'FRGW2', 0.5);
    ph.LAIMX2   = getf(v, 'LAIMX2', 0.75);
    ph.DLAI     = getf(v, 'DLAI', 0.45);
    ph.ALAI_MIN = getf(v, 'ALAI_MIN', 0);
    ph.beta_decay = getf(v, 'beta_decay', 1.0);
    ph.yr_cur   = getf(v, 'yr_cur', 11);
    ph.yr_full  = getf(v, 'yr_full', 12);
    defaultHarvest = 'thermal_maturity';
    isWinterWheat = vegetation.getPFTValue(s.PFT_type,'winter_wheat_cn',0)==1;
    if isWinterWheat && strcmpi(c.LAI_mode,'data')
        defaultHarvest = 'observed_lai_zero';
    end
    ph.crop_harvest_rule = lower(char(getf(v,'crop_harvest_rule',defaultHarvest)));
    if ~ismember(ph.crop_harvest_rule,{'thermal_maturity','observed_lai_zero'})
        error('buildVeg:badCropHarvestRule','Unknown crop_harvest_rule "%s".',ph.crop_harvest_rule);
    end
    if strcmp(ph.crop_harvest_rule,'observed_lai_zero') ...
            && (~isWinterWheat || ~strcmpi(c.LAI_mode,'data'))
        error('buildVeg:observedHarvestRequiresWheatData', ...
            'observed_lai_zero requires Winter Wheat with LAI_mode=data.');
    end
    veg.phenology = ph;

    % --- hydraulic: plant hydraulic scheme + dynamic params + tree fallback ---
    %   friendly dynamic names -> legacy dyn_* (retention_slope->dyn_a_retention,
    %   water_capacity->dyn_W_capacity, path_*->dyn_dz_*, stem_capacitance_lag->dyn_tau_stem_s)
    hy = struct();
    hy.hydraulic_scheme   = vegetation.hydraulics.normalizeHydraulicScheme( ...
        getf(v, 'hydraulic_scheme', 'steady_spac'));
    hy.hydraulic_gpmax    = getf(v, 'hydraulic_gpmax', 8e-8);
    hy.hydraulic_p50      = getf(v, 'hydraulic_p50', -3.8);
    hy.hydraulic_aa       = getf(v, 'hydraulic_aa', 3);
    hy.hydraulic_psil_min = getf(v, 'hydraulic_psil_min', -6);
    % Validation (Init_Plant.m:93-96): psil_min is the lowest leaf water potential
    % in the vulnerability lookup; a non-negative or non-finite value is invalid.
    if ~isfinite(hy.hydraulic_psil_min) || hy.hydraulic_psil_min >= 0
        error('buildVeg:badPsilMin', ...
            'hydraulic_psil_min must be a finite negative value in MPa, got %g.', hy.hydraulic_psil_min);
    end
    hy.dyn_a_retention_root = getf(v, 'retention_slope_root', 0.002);
    hy.dyn_a_retention_stem = getf(v, 'retention_slope_stem', 0.002);
    hy.dyn_a_retention_leaf = getf(v, 'retention_slope_leaf', 0.002);
    hy.dyn_W_capacity_root  = getf(v, 'water_capacity_root', 0.003);
    hy.dyn_W_capacity_stem  = getf(v, 'water_capacity_stem', 0.005);
    hy.dyn_W_capacity_leaf  = getf(v, 'water_capacity_leaf', 0.001);
    hy.dyn_dz_rs       = getf(v, 'path_root_to_stem', 0.5);
    hy.dyn_dz_sl       = getf(v, 'path_stem_to_leaf', 0.5);
    hy.dyn_tau_stem_s  = getf(v, 'stem_capacitance_lag', 0);
    hy.tree_DBH   = getf(v, 'tree_DBH', 4.2);
    hy.tree_CAH   = getf(v, 'tree_CAH', 4.13);
    hy.tree_Ksmax = getf(v, 'tree_Ksmax', 5);
    veg.hydraulic = hy;

    % --- root: system + soil-root conductance params (scalars; RAI is V2) ---
    r = struct();
    r.root_dist_mode   = getf(v, 'root_dist_mode', 0);
    % Validation (Init_Plant.m:113-121): only 0 (Jackson) / 1 (measured CSV) are
    % defined; reject other values instead of silently routing them to measured.
    if ~ismember(r.root_dist_mode, [0 1])
        error('buildVeg:invalidRootDistMode', ...
            'root_dist_mode must be 0 (Jackson) or 1 (measured root CSV), got %g.', r.root_dist_mode);
    end
    r.root_measured_file = char(getf(v, 'root_measured_file', ''));
    r.root_type        = getf(v, 'root_type', 1);
    r.root_max_depth   = getf(v, 'root_max_depth', 200);
    r.root_rroot_mm    = getf(v, 'root_rroot_mm', 0.5);
    r.gsr_formula      = char(getf(v, 'gsr_formula', 'katul'));
    r.root_hydraulic_krmax     = getf(v, 'root_hydraulic_krmax', 0);
    r.root_side_mode           = char(getf(v, 'root_side_mode', 'none'));
    r.root_radial_impedance_Pr = getf(v, 'root_radial_impedance_Pr', 4e11);
    r.root_contact_enabled     = toFlag(getf(v, 'root_contact_enabled', false));
    r.root_contact_psi_start   = getf(v, 'root_contact_psi_start', -0.05);
    r.root_contact_psi_end     = getf(v, 'root_contact_psi_end', -1.5);
    % Validation (Init_Plant.m:40-48): when the root-contact reduction is enabled,
    % both thresholds must be finite negative MPa with psi_start > psi_end
    % (drier onset above wetter completion).
    if r.root_contact_enabled ~= 0
        if ~isfinite(r.root_contact_psi_start) || ~isfinite(r.root_contact_psi_end) || ...
                r.root_contact_psi_start >= 0 || r.root_contact_psi_end >= 0 || ...
                r.root_contact_psi_start <= r.root_contact_psi_end
            error('buildVeg:badRootContact', ...
                ['root_contact thresholds must be finite negative MPa with ', ...
                 'root_contact_psi_start > root_contact_psi_end.']);
        end
    end
    r.root_solver_tol      = getf(v, 'root_solver_tol', 1e-6);
    r.root_solver_max_iter = getf(v, 'root_solver_max_iter', 50);
    veg.root = r;

    % --- photosynthesis: leaf photosynthesis + stomatal (Medlyn) + LUNA ---
    p = struct();
    p.Photo_type           = char(getf(v, 'photo_type', 'C3'));
    p.Photo_Vcmax_flag     = vcmaxFlag(getf(v, 'vcmax_source', 'manual'));    % manual=2, pft_table=0
    p.Photo_Vcmax_manual   = getf(v, 'vcmax_manual', 105);
    p.lLUNA                = toFlag(getf(v, 'luna', false));                  % bool -> 0/1
    p.cab_from_luna        = toFlag(getf(v, 'cab_from_luna', false));
    p.Photo_LeafResp_Method = leafRespMethod(getf(v, 'leaf_resp_method', 'vcmax_proportional')); % 0/2
    p.photo_Rd_C3 = getf(v, 'photo_Rd_C3', 0.015);     % no legacy VegetationData prop (passthrough)
    p.photo_Rd_C4 = getf(v, 'photo_Rd_C4', 0.025);
    p.gs_g0       = getf(v, 'gs_g0', 0.0001);
    p.gs_g1       = getf(v, 'gs_g1', 4.86);
    p.gs_iter_init = getf(v, 'gs_iter_init', 0.4);
    % Optional heat-induced stomatal opening (additive g_heat after Medlyn).
    % Default OFF until site-calibrated values exist; see
    % docs/plant/stomata_heat_induced_gs_design.md.
    p.heat_gs_enable   = toFlag(getf(v, 'heat_gs_enable', false));
    p.heat_gs_T_on_C   = getf(v, 'heat_gs_T_on_C', 32.0);    % [degC]
    p.heat_gs_T_full_C = getf(v, 'heat_gs_T_full_C', 38.0);  % [degC]
    p.heat_gs_max      = getf(v, 'heat_gs_max', 0.08);       % [mol m-2 s-1]
    p.heat_gs_wt_off   = getf(v, 'heat_gs_wt_off', 0.4);     % [-]
    p.heat_gs_wt_on    = getf(v, 'heat_gs_wt_on', 0.8);      % [-]
    p.heat_gs_D_lo_kPa = getf(v, 'heat_gs_D_lo_kPa', 1.5);   % [kPa]
    p.heat_gs_D_hi_kPa = getf(v, 'heat_gs_D_hi_kPa', 3.5);   % [kPa]
    p.heat_gs_A_lo     = getf(v, 'heat_gs_A_lo', 1.0);       % [umol m-2 s-1]
    p.heat_gs_A_hi     = getf(v, 'heat_gs_A_hi', 5.0);       % [umol m-2 s-1]
    p.heat_gs_cap      = getf(v, 'heat_gs_cap', 0.6);        % [mol m-2 s-1]
    veg.photosynthesis = p;

    % --- V2: per-layer root distribution RAI (needs grid + measured root) ---
    if nargin >= 2 && ~isempty(built) && isfield(built, 'domain')
        [veg.root.RAI, veg.root.rootProfile] = buildRootDistribution(veg.root, built.domain, raw);
    end
end

% ---------------------------------------------------------------------------
function [RAI, prof] = buildRootDistribution(rootP, domain, raw)
% Per-layer root area index [m2/m2], bottom-up (layer 1 = deepest), matching the
% legacy convention. Two sources by root_dist_mode:
%   0 jackson  : Jackson 1996 beta profile via vegetation.hydraulics.calculateRootAreaIndex1
%   1 measured : overlap-weighted RLD from raw.root mapped to model layers, then
%                RAI = 2*pi*r_root*RLD*dz (cylinder lateral area per ground area).
% MEASURED is the HEAD-correct logic (RL=tot_depth, Lmax=root_max_depth); it is
% NOT gate-A'd element-for-element against the working-tree calculateRootDistribution
% because that file's measured branch currently references undefined RL/Lmax
% (mid-edit). Verified by physical checks instead; restore element-wise gate-A
% once calculateRootDistribution stabilises. Ref (HEAD): src/+vegetation/+root/
% calculateRootDistribution.m + src/+vegetation/+hydraulics/calculateRootAreaIndex1.m.
% Editor: Jizixing Zhu; Date: 2026-06-27.
    RL    = domain.tot_depth;        % [cm] total soil profile depth
    Lmax  = rootP.root_max_depth;    % [cm] maximum rooting depth
    NN    = domain.NN;
    NL    = domain.NL;
    DeltZ = domain.dz(:);            % [cm, NL x 1] bottom-up (dz(1)=deepest)

    prof = struct('mode', '', 'root_length_density_m_m3', [], 'total_root_length_m_m2', NaN);

    if rootP.root_dist_mode == 0
        % --- Jackson baseline (reuse the legacy helper; gate-A confirms wiring) ---
        prof.mode = 'jackson';
        RAI = vegetation.hydraulics.calculateRootAreaIndex1(rootP.root_type, Lmax, RL, NN, DeltZ);
        RAI = RAI(:);
        return;
    end

    % --- measured: overlap-weighted RLD from raw.root onto the model layers ---
    prof.mode = 'measured';
    if ~istable(raw.root) || isempty(raw.root)
        error('buildVeg:noMeasuredRoot', 'root_dist_mode=1 needs raw.root (root.csv).');
    end
    depth_cm = raw.root.depth_cm(:);
    rld      = raw.root.root_length_density_m_m3(:);
    ok = isfinite(depth_cm) & isfinite(rld) & depth_cm >= 0 & rld >= 0;
    depth_cm = depth_cm(ok);  rld = rld(ok);
    if isempty(depth_cm)
        error('buildVeg:emptyMeasuredRoot', 'root.csv has no valid non-negative rows.');
    end
    [depth_cm, order] = sort(depth_cm);  rld = rld(order);
    [depth_cm, ia]    = unique(depth_cm, 'stable');  rld = rld(ia);
    meas_top = [0; depth_cm(1:end-1)];   % measured interval [top, bottom] from surface
    meas_bot = depth_cm;

    % model-layer top/bottom from surface (bottom-up: ML=1 deepest)
    layer_top = zeros(NL, 1);  layer_bot = zeros(NL, 1);  layer_ctr = zeros(NL, 1);
    from_bottom = 0;
    for ML = 1:NL
        layer_top(ML) = RL - (from_bottom + DeltZ(ML));
        layer_bot(ML) = RL - from_bottom;
        layer_ctr(ML) = 0.5 * (layer_top(ML) + layer_bot(ML));
        from_bottom   = from_bottom + DeltZ(ML);
    end

    rld_layer = zeros(NL, 1);
    for ML = 1:NL
        lt = max(layer_top(ML), 0);
        lb = min(layer_bot(ML), Lmax);
        if lb <= lt; continue; end
        wsum = 0; osum = 0;
        for k = 1:numel(rld)
            ov = min(lb, meas_bot(k)) - max(lt, meas_top(k));
            if ov > 0; wsum = wsum + rld(k) * ov;  osum = osum + ov; end
        end
        if osum > 0
            rld_layer(ML) = wsum / osum;
        else
            rld_layer(ML) = interp1(depth_cm, rld, layer_ctr(ML), 'linear', 0);
        end
    end
    rld_layer = max(rld_layer, 0);

    r_root_m = 5e-4;
    if rootP.root_rroot_mm > 0; r_root_m = rootP.root_rroot_mm * 1e-3; end
    DeltZ_m = DeltZ / 100;                       % cm -> m
    RAI = 2 * pi * r_root_m * rld_layer .* DeltZ_m;     % [m2/m2]
    if sum(RAI) <= 0
        error('buildVeg:zeroMeasuredRAI', 'measured root profile gave zero RAI on the grid.');
    end
    prof.root_length_density_m_m3 = rld_layer;
    prof.total_root_length_m_m2   = sum(rld_layer .* DeltZ_m);
end

% ---------------------------------------------------------------------------
function val = getf(s, name, default)
% Field by name with a fallback default (TOML key may be absent).
% Editor: Jizixing Zhu; Date: 2026-06-27.
    if isfield(s, name); val = s.(name); else; val = default; end
end

% ---------------------------------------------------------------------------
function f = toFlag(x)
% Logical / numeric -> 0/1 double.
% Editor: Jizixing Zhu; Date: 2026-06-27.
    f = double(logical(x));
end

% ---------------------------------------------------------------------------
function flag = vcmaxFlag(src)
% Vcmax25 source string -> legacy Photo_Vcmax_flag (VegetationData: 0=PFT/2=manual).
% Editor: Jizixing Zhu; Date: 2026-06-27.
    src = lower(strtrim(char(src)));
    switch src
        case 'manual',    flag = 2;
        case 'pft_table', flag = 0;
        otherwise
            error('buildVeg:badVcmaxSource', 'vcmax_source "%s" not in {manual, pft_table}', src);
    end
end

% ---------------------------------------------------------------------------
function m = leafRespMethod(name)
% leaf Rd method string -> legacy Photo_LeafResp_Method (0=vcmax_proportional, 2=atkin2015).
% Editor: Jizixing Zhu; Date: 2026-06-27.
    name = lower(strtrim(char(name)));
    switch name
        case 'vcmax_proportional', m = 0;
        case 'atkin2015',          m = 2;
        otherwise
            error('buildVeg:badLeafResp', 'leaf_resp_method "%s" not in {vcmax_proportional, atkin2015}', name);
    end
end
