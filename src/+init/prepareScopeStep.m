function step = prepareScopeStep(model, kt)
%PREPARESCOPESTEP  Build SCOPE/ebal input structs for one time step.
%
%   step = init.prepareScopeStep(model, kt) slices the initialized model at time
%   index kt and returns the dynamic SCOPE input bundle formerly mixed into
%   io.select_input/io.loadTimeSeries. It does not call RTMo/RTMf/RTMt and does
%   not update prognostic state.
%
%   Inputs:
%     model [in] struct  output of runCOMPASS/init.prepareRuntime, with
%                       model.built (incl. built.p.scope_static static optics cache) and model.runtime.scope.
%     kt    [in] scalar  1-based model time index.
%
%   Returns step [struct]:
%     .Soil_bare  [struct] immutable current-step bare-soil template before
%                          any accepted/predicted snow surface is blended
%     .Soil       [struct] working soil/snow optical-aerodynamic state for
%                          RT/ebal; stageCanopyEnergy may rebuild it from
%                          Soil_bare with the current precipitation predictor
%     .leafbio    [struct] leaf optical + photosynthesis baseline fields
%     .canopy     [struct] canopy structure and current LAI/hc/SAI
%     .meteo      [struct] meteorological driver slice at kt
%     .atmo       [struct] MODTRAN atmosphere table plus current air temperature
%     .angles     [struct] solar/view geometry slice at kt [deg]
%     .xyt        [struct] SCOPE time vectors from runtime.scope
%     .rad/.thermal/.fluxes/.profiles  per-step output containers
%
%   Closure:
%     The returned structs are per-step working values. Forcing slices live in
%     step.meteo/step.angles, static optical caches stay in built.p.scope_static,
%     and accepted prognostic state remains in model.built.Y.  The initial Soil
%     view contains accepted snow only as a seed; it is not the authoritative
%     current-step snow boundary. This function does not mutate model.built.Y,
%     runtime.step, or diagnostics.
%
%   Editor: Jizixing Zhu; Date: 2026-06-28.
%   Change: 2026-07-05 Jizixing Zhu -- read the static optics cache from built.p.scope_static (was model.runtime.scope_static).
%   Change: 2026-07-05 Jizixing Zhu -- buildStepLeafBio LUNA view sourced from Y.plant_acclimation.luna_state + built.p.plant_acclimation.luna_enabled (was rt.luna).
%   Change: 2026-07-17 Jizixing Zhu -- blend accepted ClimaLand bulk-snow
%     cover, optics, temperature, and thermal inertia into the SCOPE surface.
%   Change: 2026-07-17 Jizixing Zhu -- retain an immutable bare-surface
%     template so current-step snowfall can replace the accepted-snow view
%     before RTMo without compounding optical blends.
%   Change: 2026-07-17 Jizixing Zhu -- document Soil_bare, the replaceable
%     accepted-snow seed, and current-step SCOPE surface ownership.

% Change: 2026-09-29 Wenjie Liu -- explicit BGC replay context preserves
% accepted skin-temperature memory at the forcing-cycle seam only.
    validateScopeModel(model, kt);

    built = model.built;
    scope = built.parameters.scope;
    veg = built.parameters.veg;
    rt = model.runtime;

    scopeStatic = built.p.scope_static;   % static optics cache used by RTMo/RTMt/ebal
    step = struct();
    carryAtCycleStart=isfield(rt,'biogeochemistry_cycle') ...
        &&rt.biogeochemistry_cycle.active;
    step.Soil_bare = buildStepSoil(built, scope, kt, scopeStatic,carryAtCycleStart);
    step.Soil    = step.Soil_bare;
    step.Soil    = surface.snow.applySurfaceProperties( ...
        step.Soil, built.Y.snow, built.parameters.snow, scopeStatic.spectral);
    step.leafbio = buildStepLeafBio(built, scope, veg);
    step.leafbio.emis = scopeStatic.emis_static;   % leaf thermal emissivity for RTMt/ebal (COMPASS.m:849)
    step.leafopt = scopeStatic.leafopt_s;           % static fluspect leaf optics for RTMo (COMPASS.m:850)
    step.canopy  = buildStepCanopy(built, veg, kt, scopeStatic);
    step.meteo   = buildStepMeteo(built, kt);
    step.atmo    = buildStepAtmo(scopeStatic, step.meteo);
    step.angles  = buildStepAngles(built, scope, kt);
    step.xyt     = rt.scope.xyt;

    step.rad      = rt.scope.output_template.rad;
    step.thermal  = rt.scope.output_template.thermal;
    step.fluxes   = rt.scope.output_template.fluxes;
    step.profiles = rt.scope.output_template.profiles;
    step.scope_static = scopeStatic;
end

% -------------------------------------------------------------------------
function validateScopeModel(model, kt)
% Guard the runtime contract before building per-step SCOPE inputs.
% Editor: Jizixing Zhu; Date: 2026-06-28.
% Change: 2026-07-16 Jizixing Zhu — pass the static foliage clumping index to all canopy radiative-transfer kernels.
% Change: 2026-07-05 Jizixing Zhu -- guard checks built.p.scope_static (was model.runtime.scope_static); runtime.scope still required.
    if nargin < 1 || ~isstruct(model) || ~isfield(model, 'built') || ~isfield(model, 'runtime')
        error('prepareScopeStep:badModel', ...
            'prepareScopeStep expects a model returned by init.prepareRuntime.');
    end
    if ~isfield(model.built, 'p') || ~isfield(model.built.p, 'scope_static') || ~isfield(model.runtime, 'scope')
        error('prepareScopeStep:missingScopeRuntime', ...
            'model.built.p.scope_static and model.runtime.scope are required.');
    end
    tS = model.built.control.time.tS;
    if kt < 1 || kt > tS || kt ~= round(kt)
        error('prepareScopeStep:badKt', 'kt must be an integer in [1,%d], got %g.', tS, kt);
    end
end

% -------------------------------------------------------------------------
function Soil = buildStepSoil(built, scope, kt, scopeStatic,carryAtCycleStart)
% Build the per-step SCOPE Soil struct from current state and soil baselines.
% Editor: Jizixing Zhu; Date: 2026-06-28.
% Change: 2026-07-03 Jizixing Zhu -- add Soil.refl for RTMo: static optical band plus per-step moisture-corrected thermal band; ports COMPASS.m:868-888.
% Change: 2026-07-08 Jizixing Zhu -- Soil.Ts is now the [2x1] sunlit/shaded soil skin-temperature carry from the previous accepted ebal step; kt=1 falls back to meteo.Ta.
    so = scope.soil_optical;
    options = scope.options;
    cs = getFieldOr(so, 'cs', NaN);
    rhos = getFieldOr(so, 'rhos', NaN);
    lambdas = getFieldOr(so, 'lambdas', NaN);
    Soil = struct();
    Soil.rss           = so.rss;
    Soil.rs_thermal    = so.rs_thermal;
    Soil.CSSOIL        = getFieldOr(so, 'CSSOIL', built.parameters.veg.structure.CSSOIL);
    Soil.rbs           = so.rbs;
    Soil.SMC           = surfaceLiquidWater(built.Y);
    Soil.BSMBrightness = so.BSMBrightness;
    Soil.BSMlat        = so.BSMlat;
    Soil.BSMlon        = so.BSMlon;
    % Soil.Ts [degC, 2x1]: sunlit/shaded soil skin temperature consumed by RTMt/ebal.
    % Legacy COMPASS.m:889-893 cold-starts kt=1 from meteo.Ta, then carries the
    % previous accepted ebal thermal.Ts into kt>1. RTMo does not read this field.
    Soil.Ts            = previousAcceptedSurfaceTs(built, kt,carryAtCycleStart);

    if options.soil_heat_method == 1
        Soil.GAM = dynamicSoilThermalInertia(built);
        if ~isfinite(Soil.GAM)
            theta_s_top = built.parameters.soil.theta_s(end);
            Soil.GAM = equations.Soil_Inertia1(Soil.SMC, theta_s_top);
        end
    elseif options.soil_heat_method == 0 && all(isfinite([cs, rhos, lambdas]))
        Soil.GAM = equations.Soil_Inertia0(cs, rhos, lambdas);
    else
        Soil.GAM = NaN;
    end

    % --- Soil.refl: soil reflectance spectrum RTMo reads as rs (COMPASS.m:868-888) ---
    % rs_thermal_dyn: moisture-corrected thermal-band soil reflectance (COMPASS.m:869-871),
    % interpolating Dry->Wet by the surface saturation. Surface saturation uses the same
    % start-of-step surface liquid water (Soil.SMC) and top saturated water content the
    % soil-inertia branch above uses, matching legacy Theta_LL(NL,2)/Theta_s(NL).
    theta_s_surface = built.parameters.soil.theta_s(end);
    Sat = max(0, min(1, Soil.SMC / max(theta_s_surface, 1e-6)));
    rs_thermal_dyn = so.Rho_Thermal_Soil_Dry + ...
        (so.Rho_Thermal_Soil_Wet - so.Rho_Thermal_Soil_Dry) * Sat;
    if options.soilspectrum == 0
        rs_cur = scopeStatic.rs_static;                          % static optical band (rs_static)
        if options.rt_thermal == 0
            rs_cur(scopeStatic.spectral.IwlT) = rs_thermal_dyn;  % overwrite thermal band per step (COMPASS.m:884-885)
        end
        Soil.refl = rs_cur;
    else
        % Legacy soilspectrum==1 rebuilds the optical band each step via BSM(Soil,optipar,
        % {SMC_capacity_percent=25,film=0.015}) and then fills the thermal band by rt_thermal.
        % The dynamic surface water content is Soil.SMC; SMC_capacity_percent is the
        % fixed SCOPE BSM water-film capacity parameter, not the current soil water.
        spectral = scopeStatic.spectral;
        rs_cur = zeros(numel(spectral.wlS), 1);
        soilemp = struct('SMC_capacity_percent', 25, 'film', 0.015);
        rs_cur(spectral.IwlP) = surface.radiation.BSM(Soil, scopeStatic.library.optipar, soilemp);
        rslast_dyn = rs_cur(numel(spectral.wlP));
        switch options.rt_thermal
            case 0
                rs_cur(spectral.IwlT) = rs_thermal_dyn;
            case 1
                rs_cur(spectral.IwlT) = rslast_dyn;
            otherwise
                error('prepareScopeStep:badRtThermal', ...
                    'options.rt_thermal must be 0 or 1, got %g.', options.rt_thermal);
        end
        Soil.refl = rs_cur;
    end
end

% -------------------------------------------------------------------------
function GAM = dynamicSoilThermalInertia(built)
% Return top-surface thermal inertia from the current soil thermal auxiliary cache.
% Editor: Jizixing Zhu; Date: 2026-07-08.
    GAM = NaN;
    if ~isfield(built, 'workspace') || ~isfield(built.workspace, 'derived')
        return
    end

    d = built.workspace.derived;
    lambda_top = d.Lambda_eff(end, 2);
    rhoC_top = d.EHCAP(end, 2);
    if ~isfinite(rhoC_top) || rhoC_top <= 0
        rhoC_top = d.c_unsat(end, 2);
    end
    if lambda_top <= 0 || rhoC_top <= 0 || ~isfinite(lambda_top) || ~isfinite(rhoC_top)
        return
    end

    % Same thermal-inertia identity as equations.Soil_Inertia0, using the
    % dynamic soil-thermal cache instead of fixed optical-parameter inputs.
    lambda_si = lambda_top * 100;   % [W/(cm K)] -> [W/(m K)]
    Cv_si = rhoC_top * 1e6;         % [J/(cm3 K)] -> [J/(m3 K)]
    GAM0 = sqrt(lambda_si * Cv_si); % [J m-2 K-1 s-1/2]
    if isfinite(GAM0) && GAM0 >= 300 && GAM0 <= 3000
        GAM = GAM0;
    end
end


% -------------------------------------------------------------------------
function Ts = previousAcceptedSurfaceTs(built, kt,carryAtCycleStart)
% Return the SCOPE/ebal soil skin-temperature carry for this step.
% Editor: Jizixing Zhu; Date: 2026-07-08.
    Ts = [];
    if (kt > 1 || carryAtCycleStart) ...
            && isfield(built.Y, 'surface') && isfield(built.Y.surface, 'Ts')
        Ts0 = built.Y.surface.Ts;
        if isnumeric(Ts0) && numel(Ts0) == 2 && all(isfinite(Ts0(:)))
            Ts = Ts0(:);
        end
    end
    if isempty(Ts)
        Ts = built.forcing.hourly.Ta(kt) * ones(2, 1);
    end
end

% -------------------------------------------------------------------------
function leafbio = buildStepLeafBio(built, scope, veg)
% Build leaf optical and photosynthesis inputs for the current SCOPE step.
%
%   Inputs:
%     built [in] struct  canonical model container; supplies accepted LUNA state
%     scope [in] struct  scope.leaf_optical static leaf optics and SCOPE flags
%     veg   [in] struct  vegetation photosynthesis and structure parameters
%
%   Output:
%     leafbio [out] struct consumed by fluspect, RTMf, and vegetation.energy.ebal.
%
%   Role:
%     This function does not compute LUNA. It only reads the accepted LUNA state
%     produced after the 10-day climate buffer and daily runDailyLuna update.
%     On cold start or before the first LUNA update, Vcmo comes from manual
%     configuration or the PFT table.
%
%   Editor: Jizixing Zhu; Date: 2026-07-08.
%   Change: 2026-07-20 Jizixing Zhu — pass optional heat_gs_* params from veg.photosynthesis to leafbio (default off).
% Change: 2026-09-30 Wenjie Liu -- consume compiled leaf physiology when
% PlantCN is active; final actual-green-N bounds are staged after final LAI.
% Change: 2026-09-30 Wenjie Liu -- roll back this round's leaf-N LUNA
%   extension to the frozen COMPASS implementation; historical records above
%   are retained. Winter-wheat C/N allocation is handled separately.
    lo = scope.leaf_optical;
    ps = veg.photosynthesis;
    leafbio = lo;
    if ps.Photo_Vcmax_flag == 2
        leafbio.Vcmo = ps.Photo_Vcmax_manual;
    else
        leafbio.Vcmo = vegetation.getPFTValue(veg.structure.PFT_type, 'Vcmax25');
    end
    % LUNA capacity is lagged: runDailyLuna updates Y.plant_acclimation.luna_state
    % after a daily boundary and a warm 10-day climate buffer. The next accepted
    % step reads that state here and passes Vcmax25_LUNA/Jmax25_LUNA to biochemical.
    luna_enabled = built.p.plant_acclimation.luna_enabled;
    luna_state = built.Y.plant_acclimation.luna_state;
    leafbio.use_luna = 0;
    if luna_enabled && ~luna_state.first_day && ~isempty(luna_state.Vcmax25_z)
        leafbio.use_luna     = 1;
        leafbio.Vcmax25_LUNA = luna_state.Vcmax25_z(1);
        leafbio.Jmax25_LUNA  = luna_state.Jmax25_z(1);
        leafbio.Vcmo         = luna_state.Vcmax25_z(1);
    end
    leafbio.medlyn_g0 = ps.gs_g0;
    leafbio.medlyn_g1 = ps.gs_g1;
    % Optional heat-induced gs (additive after Medlyn). Defaults off via buildVeg.
    % See docs/plant/stomata_heat_induced_gs_design.md.
    leafbio.heat_gs_enable   = getFieldOr(ps, 'heat_gs_enable', 0);
    leafbio.heat_gs_T_on_C   = getFieldOr(ps, 'heat_gs_T_on_C', 32.0);
    leafbio.heat_gs_T_full_C = getFieldOr(ps, 'heat_gs_T_full_C', 38.0);
    leafbio.heat_gs_max      = getFieldOr(ps, 'heat_gs_max', 0.08);
    leafbio.heat_gs_wt_off   = getFieldOr(ps, 'heat_gs_wt_off', 0.4);
    leafbio.heat_gs_wt_on    = getFieldOr(ps, 'heat_gs_wt_on', 0.8);
    leafbio.heat_gs_D_lo_kPa = getFieldOr(ps, 'heat_gs_D_lo_kPa', 1.5);
    leafbio.heat_gs_D_hi_kPa = getFieldOr(ps, 'heat_gs_D_hi_kPa', 3.5);
    leafbio.heat_gs_A_lo     = getFieldOr(ps, 'heat_gs_A_lo', 1.0);
    leafbio.heat_gs_A_hi     = getFieldOr(ps, 'heat_gs_A_hi', 5.0);
    leafbio.heat_gs_cap      = getFieldOr(ps, 'heat_gs_cap', 0.6);
    leafbio.Type = ps.Photo_type;
    % Mesophyll conductance g_m [mol m-2 s-1 (~bar-1)]. A finite g_m turns on the mesophyll
    % CO2 drawdown (Cc = Ci - A/g_m < Ci -> lower carboxylation-site CO2 -> lower
    % photosynthesis) in vegetation.biochemical; Inf disables it (standard Farquhar Cc=Ci). Per-PFT
    % gm25 values live in the PFT table (getPFTValue 'g_m', from Knauer 2019 GCB25 Table 3);
    % the Inf default guards any PFT without a g_m entry.
    leafbio.g_m = vegetation.getPFTValue(veg.structure.PFT_type, 'g_m', Inf);
    if strcmpi(ps.Photo_type, 'C4')
        leafbio.Rdparam = ps.photo_Rd_C4;
    else
        leafbio.Rdparam = ps.photo_Rd_C3;
    end
    if scope.options.Cca_function_of_Cab
        leafbio.Cca = 0.25 * leafbio.Cab;
    end
    if scope.options.calc_PSI == 1
        leafbio.fqe = [leafbio.fqe / 5, leafbio.fqe];
    end
end

% -------------------------------------------------------------------------
function canopy = buildStepCanopy(built, veg, kt, scopeStatic)
% Build canopy state, aerodynamic inputs, and RT discretization for one step.
%
% LAI, hc, and SAI are the step-entry vegetation values. runOneStep may later
% update LAI/hc from accepted phenology or LUNA state and then recompute hot and
% LAI_layer before RTMo/ebal. The angular RT discretization and LIDF are static
% inputs precomputed by buildScopeStatic.
%
% Editor: Jizixing Zhu; Date: 2026-06-28.
    s = veg.structure;
    cv = veg.canopy;
    ls = built.vegetation_state.lai_state;
    canopy = struct();
    canopy.LAI = currentSeriesValue(ls.LAI, kt);
    canopy.hc  = currentSeriesValue(ls.hc, kt);
    canopy.SAI = currentSeriesValue(ls.SAI, kt);
    canopy.Method = s.Method;
    canopy.SurfaceType = s.SurfaceType;
    canopy.leafwidth = s.w_leaf_max;
    canopy.Cd = s.Cd;
    canopy.CR = s.CR;
    canopy.CD1 = s.CD1;
    canopy.Psicor = s.Psicor;
    canopy.rwc = s.rwc;
    canopy.kV = getFieldOr(cv, 'kV', 0.6);
    canopy.zo = getFieldOr(cv, 'zo', 0.1 * canopy.hc);
    canopy.d  = getFieldOr(cv, 'd', 0.67 * canopy.hc);
    canopy.hot = canopy.leafwidth / max(canopy.hc, eps);
    canopy.clumping_index = getFieldOr(cv, 'clumping_index', 1.0); % [-] old checkpoints retain random-foliage behavior

    % Atkin-style leaf respiration branch inputs consumed by vegetation.biochemical
    % when Photo_LeafResp_Method == 2; otherwise Rd25 remains Vcmax-proportional.
    canopy.leafresp_atkin      = double(getFieldOr(veg.photosynthesis, 'Photo_LeafResp_Method', 0) == 2);
    canopy.lmr_intercept_atkin = 1.756;   % [umol m-2 s-1] CLM5-FruitTree citrus
    canopy.lnc                 = 1.81;     % [gN m-2 leaf]

    % Static RT discretization and leaf-inclination distribution.
    rtc = scopeStatic.canopy_rt;
    canopy.nlayers      = rtc.nlayers;
    canopy.nlincl       = rtc.nlincl;
    canopy.nlazi        = rtc.nlazi;
    canopy.xl           = rtc.xl;
    canopy.litab        = rtc.litab;
    canopy.lazitab      = rtc.lazitab;
    canopy.layer_weight = rtc.layer_weight;
    canopy.lidf         = scopeStatic.lidf_static;   % [13x1] leaf inclination distribution

    % ebal recomputes zo/d through equations.zo_and_d each call. The step-entry
    % values above are retained as fallbacks and ebal-side static diagnostics.
end

% -------------------------------------------------------------------------
function meteo = buildStepMeteo(built, kt)
% Slice hourly forcing and derived atmospheric quantities for SCOPE.
% Editor: Jizixing Zhu; Date: 2026-06-28.
    h = built.forcing.hourly;
    d = built.forcing.derived;
    cfg = built.forcing.config;
    meteo = struct();
    meteo.z   = cfg.Z_measure;
    meteo.Rin = h.Rin(kt);
    meteo.Ta  = h.Ta(kt);
    meteo.Ts  = h.Ts(kt);
    meteo.Rli = h.Rli(kt);
    meteo.p   = h.Pg(kt);
    meteo.ea  = d.ea(kt) / 100;   % [hPa], SCOPE convention
    meteo.RH  = h.RH(kt);         % [0-1] rain/snow phase adapter
    meteo.u   = h.U(kt);
    meteo.Ca  = h.Ca(kt);
    meteo.Oa  = h.Oa(kt);
end

% -------------------------------------------------------------------------
function atmo = buildStepAtmo(scopeStatic, meteo)
% Build the MODTRAN atmosphere input struct consumed by RTMo.
%
%   Inputs:
%     scopeStatic [in] struct  static SCOPE optical cache; supplies library.atmo.M
%     meteo       [in] struct  current forcing slice; supplies meteo.Ta [degC]
%
%   Output:
%     atmo [out] struct
%       .M  [nwlM x 6] aggregated MODTRAN functions t1/t3/t4/t5/t12/t16
%       .Ta [degC]     current-step air temperature used for Planck sky radiance
%
%   Role:
%     atmo.Ta duplicates meteo.Ta only to preserve the legacy SCOPE RTMo interface
%     (COMPASS.m:894-899). It is forcing, not prognostic state. Changwu uses one
%     atmosphere file, so atmo.M is constant across the run; per-step atmosphere
%     file reloads can be restored later by replacing scopeStatic.library.atmo.M.
%
% Editor: Jizixing Zhu; Date: 2026-07-03.
    atmo = struct();
    atmo.M  = scopeStatic.library.atmo.M;
    atmo.Ta = meteo.Ta;
end

% -------------------------------------------------------------------------
function angles = buildStepAngles(built, scope, kt)
% Build the SCOPE geometry input struct for the current step.
%
%   Inputs:
%     built [in] struct  canonical model container; supplies forcing.solar.tts
%     scope [in] struct  scope.geometry static viewing configuration
%     kt    [in] scalar  1-based model time index
%
%   Output:
%     angles [out] struct
%       .tts [deg] solar zenith angle, computed from site lat/lon/time
%       .tto [deg] observer zenith angle (0 = nadir), from site geometry config
%       .psi [deg] sun-view relative azimuth, from site geometry config
%
%   Role:
%     tts is a per-step forcing-derived geometry. tto and psi are observation
%     geometry settings, currently static for Changwu (`observer_zenith=0`,
%     `relative_azimuth=90` in site.toml). They are read by RTMo/RTMf/RTMz/RTMt
%     for directional scattering and radiance calculations.
%
% Editor: Jizixing Zhu; Date: 2026-06-28.
    angles = struct();
    angles.tts = built.forcing.solar.tts(kt);
    angles.tto = scope.geometry.tto;
    angles.psi = scope.geometry.psi;
end

% -------------------------------------------------------------------------
function value = currentSeriesValue(series, kt)
% Read a scalar time series at kt, clamping kt to the available length.
% Editor: Jizixing Zhu; Date: 2026-06-28.
    idx = min(max(kt, 1), numel(series));
    value = series(idx);
end

% -------------------------------------------------------------------------
function smc = surfaceLiquidWater(Y)
% Surface-layer liquid water content from bottom-up node/layer state.
% Editor: Jizixing Zhu; Date: 2026-06-28.
    if isfield(Y, 'Theta_L') && ~isempty(Y.Theta_L)
        smc = Y.Theta_L(end, 1);
    else
        smc = Y.theta(end);
    end
end

% -------------------------------------------------------------------------
function value = getFieldOr(s, name, fallback)
% Read a struct field with fallback for not-yet-parameterized SCOPE fields.
% Editor: Jizixing Zhu; Date: 2026-06-28.
    if isfield(s, name) && ~isempty(s.(name))
        value = s.(name);
    else
        value = fallback;
    end
end
