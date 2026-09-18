function built = buildState(raw)
%BUILDSTATE  Stage-2: turn raw site containers into grid-ready model state.
%
%   built = init.buildState(raw) takes the raw containers from io.loadSite
%   (already validated by io.checkSite) and performs grid construction,
%   interpolation, van-Genuchten inversion, name mapping and unit conversion.
%
%   Container layout (single-writer by builder):
%       built.domain      grid geometry/config + site location
%       built.control     switches / numerics / time
%       built.boundary    BC codes/values + thermal-bottom Tav/At
%       built.sources     management source events
%       built.forcing     meteorological drivers + prescribed boundary series
%       built.parameters  soil, soil_n, veg, scope, plant_cn parameters
%       built.Y           prognostic initial state
%       built.p           recomputable auxiliary cache
%       built.vegetation_state  phenology/LAI daily state
%
%   Transition aliases (legacy short paths; dropped once the main loop reads the
%   canonical containers): built.switches/numerics/time == built.control.*,
%   built.grid == built.domain.grid. raw.config stays raw (not in compute logic).
%
%   Environment ownership:
%       Tav/At                    -> built.boundary.thermal_bottom
%       Rho_Thermal_Soil_Dry/Wet  -> built.parameters.scope.soil_optical
%       Z_measure/Wind_*/Rin/Ca/Oa -> built.forcing.config
%
%   Execution order:
%       domain -> control -> boundary/sources -> forcing -> soil/soil_n ->
%       veg/vegetation_state -> scope -> Y -> plant_hydraulic/plant_cn ->
%       canopy_air/scope_soil (SCOPE/canopy prognostic state).
%
%   No file reading here (that is io.loadSite); buildState is pure transform.
%
%   Editor: Jizixing Zhu; Date: 2026-06-26.
%   Change: 2026-07-05 Jizixing Zhu — add built.p.scope_static (static SCOPE optics cache via init.buildScopeStatic, moved from runtime.scope_static; state-ownership migration increment 4).
%   Change: 2026-07-05 Jizixing Zhu -- add built.p.plant_acclimation (LUNA/leaf-resp derived config via init.buildPlantAcclimationConfig; state moved to Y.plant_acclimation); increment 6a.
%   Change: 2026-07-11 Jizixing Zhu -- store prescribed water-boundary driver data under built.forcing.water_boundary.
%   Change: 2026-07-11 Jizixing Zhu -- skip plant C/N parameter and pool assembly when plant_cn_cycle is false.
%   Change: 2026-07-12 Jizixing Zhu -- pass the site-selected soil.csv vertical mapping into stage B.
%   Change: 2026-07-17 Jizixing Zhu -- add the ClimaLand-derived conservative
%     bulk-snow parameter group from the site [snow] contract.
%   Change: 2026-07-18 Jizixing Zhu -- align snow solid-phase resolution with
%     the coupled internal-energy balance tolerance.
%   Change: 2026-08-29 Codex -- initialize native AOM/soil C/N, DOM and
%     MEMC/MEND states from the measured-SOC-constrained COMPASS contract.
% Change: 2026-09-14 Jizixing Zhu -- initialize explicitly selected source-hour
%   dynamic DOM and intracellular enzyme owners after the native pool builders.

    built = struct();

    % --- domain: grid config + discretized geometry + site location ---
    gridCfg = init.buildGrid(raw.config);
    built.domain = init.buildDomain(gridCfg);     % geometry: dz/NL/NN/z_node/sumdz/tot_depth
    built.domain.grid = gridCfg;                  % the grid config that produced the geometry
    site = raw.config.site;
    if isfield(site, 'elevation'); elev = site.elevation; else; elev = NaN; end
    built.domain.site = struct('latitude', site.latitude, 'longitude', site.longitude, ...
        'elevation', elev);

    % --- control: run-control metadata, not process parameters ---
    % switches: process on/off and enum choices; numerics: solver tolerances and
    % kernels; time: simulation window, forcing clock and UTC anchors.
    built.control = struct();
    built.control.switches = init.buildSwitches(raw.config);
    built.control.numerics = init.buildNumerics(raw.config);
    built.control.time     = init.buildTime(raw.config);

    % --- boundary conditions (top/bottom BC codes + values) ---
    built.boundary = init.buildBoundary(raw.config);

    % --- sources / sinks: management events (fertilization, ...) ---
    built.sources = init.buildSources(raw);

    % --- environment scalars split into their final physical owners ---
    envGroups = init.buildEnv(raw.config);
    built.boundary.thermal_bottom = envGroups.boundary.thermal_bottom;

    % --- forcing: meteorological drivers (raw -> time axis -> solar -> hourly -> derived) ---
    built.forcing = struct('config', envGroups.forcing);
    built.forcing = init.buildForcing(raw.forcing, built);
    rawWaterBoundary = [];
    if isfield(raw, 'water_boundary')
        rawWaterBoundary = raw.water_boundary;
    end
    built.forcing.water_boundary = init.buildWaterBoundaryForcing(rawWaterBoundary);

    % --- transition aliases (legacy short paths) ---
    built.switches = built.control.switches;
    built.numerics = built.control.numerics;
    built.time     = built.control.time;
    built.grid     = built.domain.grid;

    % --- parameters: soil horizons -> per-layer parameters (bottom-up) ---
    if isfield(gridCfg, 'NS') && height(raw.soil) ~= gridCfg.NS
        error('buildState:soilRowMismatch', ...
            'soil.csv row count (%d) must match grid.NS (%d).', height(raw.soil), gridCfg.NS);
    end
    soilConstants = soil.water.constants();
    soilParams = init.buildSoilParameters(raw.soil, built.domain, raw.config);

    % decide the bulk-density source once; everything below uses it
    bd = init.buildBulkDensity(raw.initial, soilParams, built.domain, raw.config);
    soilParams.bulk_density = bd.layer;
    soilParams.porosity     = bd.porosity;

    soilParams = init.buildMoistureChar(soilParams);   % characteristic water contents (field cap/wilt)
    soilParams = init.buildSoilThermal(soilParams);    % static thermal params (needs wilting + bd)
    soilParams.constants = soilConstants;              % dry-end hd/hm thresholds used by water auxiliaries

    built.parameters = struct();
    built.parameters.soil = soilParams;
    built.parameters.snow = init.buildSnowParameters(raw.config);
    built.parameters.snow.solid_phase_energy_resolution_J_cm2 = ...
        built.control.numerics.tol_energy_balance_J_cm2;

    % --- legacy soil mode support: fixed params live under soil; state lives in Y ---
    [built.parameters.soil, soilModeState] = ...
        init.buildSoilModes(built, raw.config);

    % --- D-soil: static soil-nitrogen solver params (RHOKG/DIF1111/rates/...) ---
    built.parameters.soil_n = init.buildNitrogenSoilParams(built);

    % --- parameters.veg (V1 params + V2 RAI) + vegetation_state (V3 driver + V4 LAI) ---
    built.parameters.veg   = init.buildVeg(raw, built);                        % needs built.domain
    built.vegetation_state = init.buildVegState(built, built.parameters.veg, raw);  % needs forcing + control.time

    % --- SCOPE (S1 config + S2 spectral/canopy_rt/optical_library; needs veg.canopy) ---
    built.parameters.scope = init.buildScope(raw, built);

    % --- p: auxiliary cache ---
    built.p = struct();
    % static SCOPE optics precompute -> the auxiliary cache (moved from runtime.scope_static 2026-07-05).
    % Built here (inputs ready: parameters.scope above + parameters.veg) so runtime stays a pure allocator.
    built.p.scope_static = init.buildScopeStatic(built);
    % static LUNA/leaf-resp acclimation config -> the auxiliary cache (derived, not checkpointed;
    % moved from runtime.luna/leaf_resp 2026-07-05 increment 6a). State lives in Y.plant_acclimation.
    built.p.plant_acclimation = init.buildPlantAcclimationConfig(built);

    % --- Y: soil initial conditions (flat); plant state is added to built.Y below ---
    built.Y = init.buildInitialState(raw.initial, built, raw.config, bd);
    if built.control.switches.lSoilOrganic == 1
        [built.Y.soil_organic, built.parameters.soil_organic, organicReport] = ...
            init.buildSoilCarbon(raw, built, bd);
        built.initialization.soil_organic = organicReport;
    else
        built.parameters.soil_organic = struct();
        built.initialization.soil_organic = struct();
    end
    if built.control.switches.lSoilDOM == 1
        [built.Y.soil_dom, built.parameters.soil_dom] = ...
            init.buildSoilDOM(raw, built, bd);
        built=init.reserveInitialDOM(built);
    else
        built.parameters.soil_dom = struct();
    end
    if built.control.switches.lExplicitSoilMicrobes == 1
        [built.Y.soil_microbes, built.parameters.soil_microbes] = ...
            init.buildSoilMicrobes(built, raw.config);
    else
        built.parameters.soil_microbes = struct();
    end
    built=init.buildMendHourlyState(built,raw.config);
    built.parameters.litter = init.buildLitterParameters( ...
        raw.config, built.control.switches.lLitterInput == 1);
    [built.Y, soilModeState] = init.applySoilModesToInitialState( ...
        built.Y, soilModeState, built, raw.config);
    built.Y.soil_modes = soilModeState;

    % --- D-plant: hydraulic state/cache (P1) -> C/N params (P3) -> C/N pools (P2) ---
    [built.Y.plant_hydraulic, built.p.plant_hydraulic] = init.buildPlantHydraulicState(raw, built);
    % static legacy-shaped plantHydro view for the ebal cluster (built once; = Init_Plant [F]).
    % Reads Y.plant_hydraulic here while it still holds the config init seeds (constant thereafter).
    built.p.plantHydro = init.assemblePlantHydro(built);
    if built.control.switches.lPlantCN == 1
        built.parameters.plant_cn = init.buildPlantCNParams(raw, built);
        built.Y.plant_cn = init.buildInitialPlantPools(raw, built);
    else
        built.parameters.plant_cn = struct();
        built.Y.plant_cn = struct();
    end

    % --- SCOPE/canopy prognostic state: canopy-air OLD layers + soil-heat history ---
    scopeState = init.buildScopeState(built);
    built.Y.canopy_air = scopeState.canopy_air;
    built.Y.scope_soil = scopeState.scope_soil;
end
