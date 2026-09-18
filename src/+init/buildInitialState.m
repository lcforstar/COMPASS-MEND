function Y = buildInitialState(rawInit, built, cfg, bd, plan)
%BUILDINITIALSTATE  Interpolate initial conditions onto the grid (stage C).
%
%   Y = init.buildInitialState(rawInit, built, cfg, bd) reproduces the legacy
%   init.InitialProfile.load_and_interpolate: it maps the measured per-depth
%   initial profile (raw.initial) onto the model nodes, computes matric head
%   from water content via van Genuchten, converts nitrogen to liquid-phase
%   concentration, and flips everything to the bottom-up node convention.
%
%   Inputs:
%     rawInit  table   raw.initial from io.loadSite (columns depth_cm, theta,
%                      T_C, nh4, no3, p_atm_pa, bulk_density_gcm3)
%     built    struct  built.domain (grid) + built.parameters.soil (VG per layer)
%     cfg      struct  raw.config (uses switches.solve_soil_heat as the thermal
%                      flag for the h->theta retention correction)
%     bd       struct  shared bulk density from init.buildBulkDensity: bd.node
%                      [g/cm3, NN] drives Y.bd and the nitrogen conversion
%
%   Returns Y (prognostic state), all [NN x 1] BOTTOM-UP (node 1 = deepest):
%     Y.theta   [m3/m3]   volumetric water content
%     Y.h       [cm]      matric head (negative unsaturated), from theta via VG
%     Y.T       [degC]    temperature
%     Y.p_gas   [g/cm/s2] gas pressure (= Pa * 10, legacy internal unit)
%     Y.nh4     [mg/L]    ammonium liquid concentration (numeric ug/cm3 water)
%     Y.no3     [mg/L]    nitrate liquid concentration (numeric ug/cm3 water)
%     Y.bd      [g/cm3]   bulk density at nodes
%     Y.Theta_L [NL x 2]  layer water content (node i, node i+1)
%     Y.surface_water .DSTOR/.DSTOR0/.RS [cm] surface depression storage + runoff
%     Y.soil_gas.C_O2 [mol/cm3] gas-phase O2 concentration at nodes
%
%   Conventions / unit conversions (match legacy InitialProfile):
%     - depths_target is top-down (surface=0) during interpolation, then flipped.
%     - p_atm_pa * 10 : Pa -> g cm^-1 s^-2 (InitialProfile.m:91).
%     - nitrogen: read as total [mg/kg] or dissolved [mg/L], then stored as
%       liquid-phase concentration [mg/L] for the runtime N state.
%     - matric head recomputed from theta with per-node VG params (theta_to_h).
% Editor: Jizixing Zhu; Date: 2026-06-27.
% Change: 2026-06-28 Jizixing Zhu — add prognostic Y.surface_water (DSTOR/DSTOR0/RS).
% Change: 2026-07-05 Jizixing Zhu — add prognostic Y.snow.SWE (moved from runtime.snow; state-ownership migration template).
% Change: 2026-07-05 Jizixing Zhu — add prognostic Y.surface.Tss=NaN (moved from runtime.thermal; state-ownership migration increment 2).
% Change: 2026-07-05 Jizixing Zhu — add prognostic Y.canopy_water.{state_h2ocan,_sun,_sha}=zeros(1,nCell) (moved from runtime.canopy_water; state-ownership migration increment 3 stage 2).
% Change: 2026-07-05 Jizixing Zhu — add prognostic Y.vegetation.{phenology(whole struct),lai_state} (moved from runtime.phenology/lai_state; state-ownership migration increment 5).
% Change: 2026-07-05 Jizixing Zhu -- add prognostic Y.plant_acclimation.{leaf_resp.Ta_history, luna_state(+last_* NaN)} (moved from runtime.leaf_resp/luna.state); climate_buffer stays runtime until 6b; increment 6a.
% Change: 2026-07-14 Jizixing Zhu -- add prognostic Y.U (nodal internal energy).
% Change: 2026-07-06 Jizixing Zhu -- add Y.plant_acclimation.luna_climate_buffer = Climate10dBuffer(steps_per_day,10).toStruct() (LUNA buffer handle serialized into Y; runtime.luna deleted); increment 6b.
% Change: 2026-07-08 Jizixing Zhu -- add Y.surface.Ts [2x1] as the accepted SCOPE/ebal soil skin-temperature carry; Tss remains the scalar guarded surface temperature.
% Change: 2026-07-11 Jizixing Zhu -- skip vegetation and acclimation state assembly for the laboratory component plan.
% Change: 2026-07-16 Jizixing Zhu -- document the descriptive soil initialization filename.
% Change: 2026-07-17 Jizixing Zhu -- replace scalar SWE with the conservative
%   bulk-snow prognostic triplet SWE/liquid_water/energy.
% Change: 2026-07-20 Jizixing Zhu -- allocate the detached accepted soil-N memory schema.

    domain = built.domain;
    soil   = built.parameters.soil;
    NL = domain.NL;
    NN = domain.NN;

    % top-down node depths (surface=0 ... bottom=tot_depth), == legacy depths_target
    depths_target = flip(domain.z_node(:));

    % --- inputs (column vectors) ---
    depths_input   = rawInit.depth_cm(:);
    nh4_input      = getTableColumn(rawInit, 'nh4',  zeros(size(depths_input)));
    no3_input      = getTableColumn(rawInit, 'no3',  zeros(size(depths_input)));
    moisture_input = rawInit.theta(:);
    temp_input     = rawInit.T_C(:);
    pressure_input = getTableColumn(rawInit, 'p_atm_pa', depths_input) * 10;   % Pa -> g/cm/s2

    % node bulk density comes pre-decided from init.buildBulkDensity (shared bd);
    % bd.node is bottom-up, the interpolation frame here is top-down -> flip.
    bd_node_td = flip(bd.node(:));

    % --- per-layer VG params in TOP-DOWN order (legacy uses flip(SoilProperties.*)) ---
    theta_s_vec = flip(soil.theta_s(:));
    theta_r_vec = flip(soil.theta_r(:));
    alpha_vec   = flip(soil.alpha(:));
    n_vec       = flip(soil.n(:));

    % --- C.1 temperature: linear interpolation + extrapolation ---
    T_target = init.calcSoilTemp(depths_input, temp_input, depths_target);

    % --- C.2 water content ---
    % Restart profiles already on the grid (theta in [0,1] on identical depths):
    % use theta directly to conserve storage. Otherwise interpolate via matric head.
    is_theta    = all(isfinite(moisture_input)) && all(moisture_input >= -1e-9) ...
                  && all(moisture_input <= 1 + 1e-9);
    same_grid   = numel(depths_input) == numel(depths_target) ...
                  && max(abs(depths_input - depths_target)) < 1e-9;
    if is_theta && same_grid
        theta_target = moisture_input;
    else
        enableThermal = getSwitch(cfg, 'solve_soil_heat', true);
        theta_target = moistureViaMatricHead(depths_input, moisture_input, depths_target, ...
            T_target, theta_s_vec, theta_r_vec, alpha_vec, n_vec, enableThermal);
    end

    % --- C.3 gas pressure ---
    p_gas_target = init.interpLinearWithExtrapolation(depths_input, pressure_input, depths_target);

    % --- C.5 nitrogen: interpolate, convert to LIQUID-phase concentration ---
    % The model state C_NH4/C_NO3 is the liquid-phase concentration [mg/L]; total
    % storage is (theta + RHOKG)*C (NH4) and theta*C (NO3) (calculateStorage.m),
    % where RHOKG = bulk_density * Kd is the (dimensionless) NH4 adsorption buffer.
    % switches.nitrogen_input_phase selects how soil_initial_profile.csv columns are read:
    %   'dissolved' : columns are already liquid concentration [mg/L] -> use directly
    %   'total'     : columns are total mineral N per dry soil [mg/kg] (e.g. KCl-
    %                 extractable) -> per bulk volume (x bulk_density), then divide
    %                 by the buffer to recover the liquid concentration.
    phase = getSwitchString(cfg, 'nitrogen_input_phase', 'total');
    switch lower(phase)
        case 'dissolved'
            nh4_target = max(init.interpLinearWithExtrapolation(depths_input, nh4_input, depths_target), 0);
            no3_target = max(init.interpLinearWithExtrapolation(depths_input, no3_input, depths_target), 0);
        case 'total'
            % per bulk volume [ugN/cm3] = interp(mg/kg) * bulk_density  (legacy conversion)
            [nh4_bulk, no3_bulk] = init.calcNitrogenConcentration(depths_input, nh4_input, ...
                no3_input, depths_target, bd_node_td, struct());
            % NH4 adsorption buffer RHOKG = bulk_density * Kd at nodes (Kd = soil kg_cm3g).
            % Kd is per layer (bottom-up): map to nodes with the thickness-weighted
            % element->node helper, then flip to the top-down working frame.
            if isfield(soil, 'kg')
                Kd_node = flip(equations.interpolateElementToNode(soil.kg(:), domain.dz(:)));
            else
                Kd_node = zeros(NN, 1);   % no adsorption data -> NH4 treated like NO3
            end
            rhokg_node = bd_node_td .* Kd_node;        % [-]
            nh4_target = nh4_bulk ./ max(theta_target + rhokg_node, 1e-12);   % [mg/L] liquid
            no3_target = no3_bulk ./ max(theta_target, 1e-12);               % [mg/L] liquid
        otherwise
            error('buildInitialState:badPhase', ...
                'unknown switches.nitrogen_input_phase "%s" (expected total|dissolved)', phase);
    end

    % --- C.6 matric head from theta (per-node VG), then physical-range check ---
    h_target = recomputeHeadFromTheta(theta_target, theta_s_vec, theta_r_vec, alpha_vec, n_vec);
    [theta_target, h_target] = validatePhysical(theta_target, h_target, ...
        theta_s_vec, theta_r_vec, alpha_vec, n_vec);

    % --- C.7 flip everything to bottom-up (node 1 = deepest) ---
    Y = struct();
    Y.theta = flip(theta_target(:));
    Y.h     = flip(h_target(:));
    Y.T     = flip(T_target(:));
    Y.p_gas = flip(p_gas_target(:));
    Y.U     = nan(NN, 1);   % [J/cm³] nodal internal energy (computed after SFCC init)
    Y.nh4   = flip(nh4_target(:));
    Y.no3   = flip(no3_target(:));
    Y.urea  = zeros(NN, 1);   % [ug/cm3, NN x 1] urea pool (0 at IC; enters via fertilizer, hydrolyses to NH4)
    Y.bd    = bd.node(:);

    % layer water content NL x 2 (node i, node i+1), bottom-up
    Y.Theta_L = [Y.theta(1:NL), Y.theta(2:NL+1)];

    % --- prognostic N pools that carry cross-step memory (zero at IC) ---
    % These are accepted state (granule fertilizer buffer + freeze-entrapped N),
    % so they belong in Y, NOT in the per-step workspace. Legacy creates them as
    % zeros in the StateVariables constructor (StateVariables.m:253-255); the Y
    % single-writer reproduces that here so init.allocateWorkspace never has to
    % touch Y. Both are 0 at t=0 (no fertilizer applied, no ice yet).
    Y.C_NH4_granule = zeros(NN, 1);   % [ug/cm3, NN x 1]  undissolved NH4 granule pool
    Y.C_NH4_ice     = zeros(NL, 2);   % [mg/L eq, NL x 2] NH4 entrapped in ice
    Y.C_NO3_ice     = zeros(NL, 2);   % [mg/L eq, NL x 2] NO3 entrapped in ice

    % --- prognostic surface-water store (legacy setBoundaryCondition.m:109-111) ---
    % Depression storage and surface runoff carry cross-step memory, so they are
    % accepted state. Zero at t=0 (no ponding, no runoff yet); DSTOR0 holds the
    % previous-step DSTOR for the storage-change term.
    Y.surface_water = struct('DSTOR', 0, 'DSTOR0', 0, 'RS', 0);   % DSTOR/DSTOR0 [cm], RS [cm/s]
    Y.snow = struct( ...
        'SWE', 0, ...                 % [cm] total snow water equivalent
        'liquid_water', 0, ...        % [cm] liquid water retained in snow
        'energy', 0);                 % [J cm-2] bulk snow internal energy
    Y.surface = struct( ...
        'Ts',  NaN(2, 1), ...   % [degC, 2x1] last accepted ebal sunlit/shaded soil skin temperatures; feeds next step Soil.Ts
        'Tss', NaN);            % [degC] last accepted guarded scalar surface temperature; NaN -> kt=1 guard fallback uses Ta
    includeEcosystem = true;
    if nargin >= 5 && isfield(plan, 'is_laboratory') && plan.is_laboratory
        includeEcosystem = false;
    end
    if includeEcosystem
    nCell = 1;
    if isfield(domain, 'nCell') && ~isempty(domain.nCell), nCell = domain.nCell; end   % single-column -> 1
    Y.canopy_water = struct( ...   % [mm] canopy interception reservoir (prognostic, moved from runtime.canopy_water 2026-07-05)
        'state_h2ocan', zeros(1, nCell), 'state_h2ocan_sun', zeros(1, nCell), 'state_h2ocan_sha', zeros(1, nCell));

    % vegetation running state (moved from runtime.phenology + runtime.lai_state 2026-07-05 increment 5).
    % phenology = the WHOLE preparePhenology() equivalent {state,year_mtemp_max,year,now};
    % lai_state = ALL built.vegetation_state.lai_state fields, including wt_day_min from buildLaiState.
    % The guard below is compatibility for older lai_state structs/checkpoints that predate wt_day_min.
    laiState0 = built.vegetation_state.lai_state;
    if ~isfield(laiState0, 'wt_day_min'); laiState0.wt_day_min = 1.0; end   % [-] daily min(wt_sun,wt_sha), reset to 1 each day
    Y.vegetation = struct( ...
        'phenology', struct('state', vegetation.integrated.phenology.PhenologyState(), ...
                            'year_mtemp_max', -999, 'year', 0, 'now', 0), ...
        'lai_state', laiState0);

    % plant-acclimation running state (moved from runtime.leaf_resp + runtime.luna 2026-07-05,
    % increment 6a state/config + 6b buffer). leaf_resp.Ta_history = per-step growth-temp memory (was
    % runtime.leaf_resp.Ta_history); luna_state = the whole runtime.luna.state {Vcmax25_z/Jmax25_z/
    % PNlc_z/enzs_z/first_day} PLUS the 6 last_* daily caches promoted to a FIXED field set (init NaN
    % = no LUNA daily update yet; ~isfinite(last_Cab_nlc) is the gate). luna_climate_buffer = the LUNA
    % 10-day Climate10dBuffer serialized to a plain struct (was the runtime.luna.climate_buffer HANDLE;
    % increment 6b) -- reconstructed transiently on the accept path (fromStruct/push/toStruct), never a handle in Y.
    lunaBuf0 = vegetation.integrated.LUNA.Climate10dBuffer(built.p.plant_acclimation.steps_per_day, 10).toStruct();
    Y.plant_acclimation = struct( ...
        'leaf_resp', struct('Ta_history', NaN(built.control.time.tS, 1)), ...
        'luna_state', struct( ...
            'Vcmax25_z', [], 'Jmax25_z', [], 'PNlc_z', [], 'enzs_z', [], 'first_day', true, ...
            'last_Vcmax25', NaN, 'last_Jmax25', NaN, 'last_FNCa', NaN, ...
            'last_PNlc', NaN, 'last_PNcb', NaN, 'last_Cab_nlc', NaN), ...
        'luna_climate_buffer', lunaBuf0);
    else
        Y.canopy_water = struct();
        Y.vegetation = struct();
        Y.plant_acclimation = struct();
    end

    % Soil O2 is a prognostic soil-gas state. Direct measurements are rarely
    % available, so the default initial guess is atmospheric equilibrium at the
    % initial gas pressure and temperature profile. Soil O2 then evolves through
    % soil.gas.solveO2Transport; this initial value is not reused as a boundary.
    Y.soil_gas = struct();
    Y.soil_gas.C_O2 = initialSoilOxygenProfile(Y.p_gas, Y.T, built);

    % Soil nitrogen owns accepted cross-step chemistry/transport memory.  The
    % arrays are materialized after the first soil-auxiliary refresh, when the
    % effective liquid/ice endpoint is available.  Keeping the container in Y
    % from construction makes checkpoint ownership explicit without copying
    % static kinetic parameters into prognostic state.
    % This function retains a local variable named `soil`; use FEVAL so that
    % MATLAB resolves the +soil package rather than that local parameter view.
    Y.soil_nitrogen = feval('soil.nitrogen.emptyAcceptedState');
end

% ---------------------------------------------------------------------------
function C_O2 = initialSoilOxygenProfile(p_gas, T, built)
%INITIALSOILOXYGENPROFILE  Atmospheric-equilibrium initial soil O2 profile.
%
% Inputs:
%   p_gas  [g cm-1 s-2, NN x 1]  initial gas pressure (= Pa*10)
%   T      [degC, NN x 1]        initial soil temperature
%   built  struct                forcing.gas / forcing.hourly / forcing.derived
%
% Returns:
%   C_O2   [mol/cm3, NN x 1]     initial gas-phase O2 concentration
%
% Closure:
%   Use the same moist-air partial-pressure convention as
%   init.buildAtmosphericGasDrivers and soil.gas.calcAtmO2. The atmospheric O2
%   mole fraction is a dry-air mixing ratio; water vapour displaces dry gases.
% Editor: Jizixing Zhu; Date: 2026-06-27.
    phys = equations.physicalConstants();
    Oa = phys.Oa_mmolmol_default;     % [mmol/mol]
    ea = 0;                           % [Pa]

    if isfield(built, 'forcing') && isfield(built.forcing, 'hourly') ...
            && isfield(built.forcing.hourly, 'Oa') && ~isempty(built.forcing.hourly.Oa)
        Oa0 = built.forcing.hourly.Oa(1);
        if isfinite(Oa0) && Oa0 > 0
            Oa = Oa0;
        end
    end
    if isfield(built, 'forcing') && isfield(built.forcing, 'derived') ...
            && isfield(built.forcing.derived, 'ea') && ~isempty(built.forcing.derived.ea)
        ea0 = built.forcing.derived.ea(1);
        if isfinite(ea0) && ea0 >= 0
            ea = ea0;
        end
    end

    C_O2 = soil.gas.calcAtmO2(p_gas(:) / 10, T(:), Oa, ea);
end

% ---------------------------------------------------------------------------
function h = recomputeHeadFromTheta(theta, theta_s_vec, theta_r_vec, alpha_vec, n_vec)
% Per-node matric head from water content (legacy InitialProfile.m:196-201).
% Node i uses layer VG params at index min(max(i-1,1), NL).
% Editor: Jizixing Zhu; Date: 2026-06-27.
    NN = numel(theta);
    NL = numel(theta_s_vec);
    h = zeros(NN, 1);
    for i = 1:NN
        k = min(max(i - 1, 1), NL);
        h(i) = equations.van_genuchten_theta_to_h(theta(i), ...
            theta_s_vec(k), theta_r_vec(k), alpha_vec(k), n_vec(k));
    end
end

% ---------------------------------------------------------------------------
function [theta, h] = validatePhysical(theta, h, theta_s_vec, theta_r_vec, alpha_vec, n_vec)
% Clamp theta into [theta_r, theta_s] and recompute h (legacy validatePhysicalRange).
% Field-capacity fallback for over-saturation = 0.75*theta_s (legacy fallback when
% SoilProperties.Field_mc is empty).
% Editor: Jizixing Zhu; Date: 2026-06-27.
    NN = numel(theta);
    NL = numel(theta_s_vec);
    for i = 1:NN
        k = min(max(i - 1, 1), NL);
        ts = theta_s_vec(k); tr = theta_r_vec(k);
        changed = false;
        if theta(i) > ts
            theta(i) = 0.75 * ts;   % field-capacity fallback
            changed = true;
        end
        if theta(i) < tr
            theta(i) = tr + 1e-6;
            changed = true;
        end
        if changed
            h(i) = equations.van_genuchten_theta_to_h(theta(i), ts, tr, alpha_vec(k), n_vec(k));
        end
    end
end

% ---------------------------------------------------------------------------
function theta_interp = moistureViaMatricHead(depths_in, theta_in, depths_target, ...
        T_target, theta_s_vec, theta_r_vec, alpha_vec, n_vec, enableThermal)
% Measured-profile path (legacy calcMoistureMatricHead): theta->h at measurement
% points, log-space interpolation of h, h->theta at target nodes.
% Editor: Jizixing Zhu; Date: 2026-06-27.
    h_measured = zeros(numel(theta_in), 1);
    for i = 1:numel(depths_in)
        k = findNearestNode(depths_in(i), depths_target);
        h_measured(i) = equations.van_genuchten_theta_to_h(theta_in(i), ...
            theta_s_vec(k), theta_r_vec(k), alpha_vec(k), n_vec(k));
    end
    log_h        = log10(abs(h_measured) + 1e-10);
    log_h_interp = init.interpLinearWithExtrapolation(depths_in, log_h, depths_target);
    h_interp     = -10 .^ log_h_interp;
    theta_interp = zeros(numel(h_interp), 1);
    for i = 1:numel(depths_target)
        k = findNearestNode(depths_target(i), depths_target);
        theta_interp(i) = equations.van_genuchten_h_to_theta(h_interp(i), ...
            theta_s_vec(k), theta_r_vec(k), alpha_vec(k), n_vec(k), T_target(i), enableThermal);
    end
end
