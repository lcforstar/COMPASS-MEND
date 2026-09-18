function out = ebal(in)
%EBAL  One time-step canopy–soil energy balance (two-leaf: sunlit + shaded).
%
%   Easy picture:
%     Sunlight hits leaves and soil. This function finds temperatures so that
%     energy nearly balances:  Rn ≈ H + lE + G
%     and plant water supply can limit stomata / transpiration.
%
%   Call chain (you usually do not call this directly):
%     assembleEbalInput -> solveCanopySoilEnergy -> ebal(in) -> applyEbalOutput
%
%   Usage:
%     out = vegetation.energy.ebal(in)
%
%   Returns (out):
%     .iteration         How the solve went (counts, residual [W m-2])
%     .fluxes            Energy [W m-2], CO2, water diagnostics
%     .radiation         Radiation after thermal update
%     .gap               Gap fractions (from input, passed through)
%     .thermal           Temperatures [C], resistances
%     .profiles          Optional layer diagnostics
%     .soil_surface      Updated soil surface (SCOPE fields)
%     .root_water_uptake Root uptake per soil node (later -> Srt)
%     .frac              Component fractions [-]
%     .hydraulic         Leaf/root potentials [MPa], stress wt, supply ratio
%
%   Inputs (in) — packed by assembleEbalInput:
%     .time              kt [-], dt_s [s]
%     .soil.state         h [cmH2O], theta_liquid [m3 m-3], root_sink, G Fourier [W m-2]
%     .soil.params       theta_s/r, field_capacity, thermal soil reflectance
%     .soil.conductivity Soil hydraulic conductivity for RWU
%     .domain            nl, nn, dz_cm [cm], lat/lon [deg]
%     .control           tol_eb_W_m2, anderson, freeze couple, canopy-air flags
%     .canopy_params     lai_max, w_leaf_max_m, crop_code
%     .iteration         max_iter, tol_W_m2, relax_weight
%     .scope_static      options, spectral
%     .scope_step        rad, gap, meteo, canopy, soil_surface, interception, ...
%     .plant_hydro       Plant hydraulic parameters / init
%
%   Method (code blocks below match these numbers):
%     [A] Unpack in
%     [B] Setup: wet-leaf fraction, first-guess temperatures, soil surface props
%     [C] Iterate until energy residual is small (or max steps)
%         C1 Rn  C2 resistances  C3 photosynthesis  C4 H/lE + hydraulics
%         C5 stability L  C6 residual  C7 update T  (+ soil G inside loop)
%     [D] Warnings if not converged
%     [E] Optional vertical profiles
%     [F] Canopy totals, fluxes/thermal, root water uptake
%     [G] Pack out
%
%   Sign (Rn = H + lE + G) [W m-2]:
%     Rn>0 into surface; H>0 heats air; lE>0 ET; G>0 into soil; lE<0 condensation
%     With snow, Ts is the bare-soil skin; prognostic snow temperature is
%     diagnosed from snow enthalpy and shares only the canopy RT solve.
%
%   Reference:
%     Wallace and Verhoef (2000) two-source resistance
%     Monin and Obukhov (1954) stability
%     Oleson et al. (2013) CLM4.5 aerodynamic resistance
%     van der Tol et al. (2009) SCOPE energy balance
%     docs/plant_hydraulic_method.md
%     plant/+energy/EBAL_IO_CONTRACT.md
%
%   Original authors: C. van der Tol, J. Timmermans (2007-2013).
%   CLM4.5 enhancements: Jizixing Zhu, 2026-04-18 to 2026-04-22.
% Editor: Jizixing Zhu; Date: 2026-07-03.
% Change: 2026-07-13 Jizixing Zhu — retry a converged but nonphysical angle-resolved sunlit solution with the conservative layer-mean closure.
% Change: 2026-07-12 Jizixing Zhu — use angles.tts for the low-sun energy-balance guard.
% Change: 2026-07-03 Jizixing Zhu — created as a verbatim copy of ebal_clm45.m (renamed ebal).
% Change: 2026-07-03 Jizixing Zhu — narrow signature (drop dead Rl; SoilVariables/VanGenuchten megastructs -> new-chain soil_state{hh,Theta_LL,Srt} + vg{Theta_s,Theta_r,Theta_f}) and modernize old-style dependencies (equations.physicalConstants() struct); body physics unchanged.
% Change: 2026-07-18 Jizixing Zhu -- remove the retired Fourier-G blend from ground heat flux closure.
% Change: 2026-07-03 Jizixing Zhu — replace legacy cfg megastruct with flat new-chain ebcfg (13 scalar fields); drop isobject/isprop handle-class dual-paths (new chain is plain struct); the one cfg->assembleRWU cascade now passes a tiny struct('Grid',...) shell (sub-function unchanged); physics unchanged.
% Change: 2026-07-14 Jizixing Zhu — use conditioned Anderson least squares, defer the first-pass non-convergence warning when runOneStep owns retry, and raise the extreme-temperature fallback to 90 degC / 60 degC above air.
% Change: 2026-07-16 Jizixing Zhu — rewrite English doc header to coding-standard template; clarify kernel role under Phase-A assemble/solve/apply adapter.
% Change: 2026-07-16 Jizixing Zhu — Phase C: rename body locals to target names (soil_state.h, soil_params, domain, control, plant_hydro, soil_surface).
% Change: 2026-07-16 Jizixing Zhu — Phase B: single-struct entry out=ebal(in); unpack target names at top, pack out at end; body physics unchanged.
% Change: 2026-07-16 Jizixing Zhu — rewrite body comments into simple block sections (A-G); no physics change.
% Change: 2026-07-16 Jizixing Zhu — convert all body comments to English only.
% Change: 2026-07-16 Jizixing Zhu — drop packLegacySoilBridge; pass soil_state/soil_params/domain directly.
% Change: 2026-07-16 Jizixing Zhu — rename ebal to ebal; drop Clm45 from name.
% Change: 2026-07-20 Jizixing Zhu — apply_phase3 passes optional heat_gs_* to biochemical (default off).
% Change: 2026-07-16 Jizixing Zhu — remove migration-era unused local aliases.

% Change: 2026-07-17 Jizixing Zhu -- split bare-soil and snow vapor
% paths while retaining the shared SCOPE ground temperature.
% Change: 2026-07-18 Jizixing Zhu -- solve the bare-soil skin balance
%   while routing independent snow radiation, H, LE, and vapor exchange.
% Change: 2026-07-18 Jizixing Zhu -- expose bare/snow grid contributions
%   of the ground heat flux without changing the force-restore closure.
% Change: 2026-07-19 Jizixing Zhu -- expose the bare/snow external-energy
%   fluxes and exact local temperature tangents; empirical force-restore G
%   remains an SCOPE predictor/diagnostic and is no longer the soil BC target.
% Change: 2026-07-20 Jizixing Zhu -- materialize wet/dry and sunlit/shaded canopy
%   latent-heat components from the same accepted leaf fluxes as lEctot.
    %% ========================================================================
    %%  [A] Unpack input struct 'in'
    %%  Split the big input into local variables (names only; values unchanged)
    %% ========================================================================
    assert(isstruct(in) && isfield(in, 'soil') && isfield(in, 'scope_step'), ...
        'ebal:badInput', ...
        'in must be a target-name struct from assembleEbalInput (EBAL_IO_CONTRACT §5.2).');

    ss = in.scope_step;
    st = in.scope_static;

    iteration = struct( ...
        'maxit', in.iteration.max_iter, ...
        'maxEBer', in.iteration.tol_W_m2, ...
        'Wc', in.iteration.relax_weight, ...
        'counter', 0);
    options  = st.options;
    spectral = st.spectral;
    rad      = ss.rad;
    gap      = ss.gap;
    leafopt  = ss.leafopt;
    angles   = ss.angles;
    meteo    = ss.meteo;
    soil_surface = ss.soil_surface;
    canopy   = ss.canopy;
    leafbio  = ss.leafbio;
    xyt      = ss.xyt;
    kt       = in.time.kt;
    k        = kt;
    profiles = ss.profiles;
    dt_s     = in.time.dt_s;

    soil_state = in.soil.state;          % .h / .theta_liquid / .root_sink
    soil_params = in.soil.params;      % .theta_s / .theta_r / .field_capacity / rho_thermal_*
    domain = in.domain;                % .nl / .nn / .dz_cm / lat/lon
    control = in.control;              % anderson / tol / canopy-air flags
    canopy_params = in.canopy_params;  % .lai_max / .w_leaf_max_m / .crop_code
    plant_hydro = in.plant_hydro;
    soil_conductivity = in.soil.conductivity;
    if isfield(ss, 'interception') && ~isempty(ss.interception)
        interception = ss.interception;
    else
        interception = [];
    end

    %% ========================================================================
    %%  [B] Pre-iteration setup
    %%  Grid size, tolerances, wet-leaf fraction, first-guess T, air density,
    %%  soil surface resistance / roughness
    %% ========================================================================
    NL        = domain.nl;
    DeltZ     = domain.dz_cm;


    counter         = 0;              % outer EB iteration count
    maxit           = iteration.maxit;  % max outer iterations
    maxit_base      = maxit;
    tolEBer         = get_ebal_tolerance(control, iteration.maxEBer);  % [W m-2] residual tolerance
    Wc              = iteration.Wc;   % temperature under-relaxation (0-1)
    % --- Wet-leaf fraction fwet (0 = dry, 1 = fully wet) ---
    % Intercepted water evaporates from the leaf surface: not through stomata,
    % and does not request root water. Single-store: one fwet for sun/shade.
    % Four-quadrant: sun and shade may each have their own fwet.
    if isstruct(interception)
        if isfield(interception, 'fwet_sun')
            fwet_sun_canopy = max(0, min(1, interception.fwet_sun));
        elseif isfield(interception, 'fwet')
            fwet_sun_canopy = max(0, min(1, interception.fwet));
        else
            fwet_sun_canopy = 0;
        end
        if isfield(interception, 'fwet_sha')
            fwet_sha_canopy = max(0, min(1, interception.fwet_sha));
        elseif isfield(interception, 'fwet')
            fwet_sha_canopy = max(0, min(1, interception.fwet));
        else
            fwet_sha_canopy = 0;
        end
    else
        fwet_sun_canopy = 0;
        fwet_sha_canopy = 0;
    end
    use_ebal_anderson = logical(control.use_anderson);

    CONT            = 1;              % 1 = keep iterating; 0 = done (balanced or maxit)
    Ta              = meteo.Ta;
    ea              = meteo.ea;
    Ca              = meteo.Ca;
    Ts              = soil_surface.Ts;
    p               = meteo.p;
    if options.soil_heat_method < 2 && options.simulation == 1
        Deltat          = dt_s;           %           Duration of the time interval (s)
        x       = [1:12; 1:12]' .* Deltat;
        Tsold = soil_surface.Tsold;
    end

    nl = canopy.nlayers;

    Rnuc  = rad.Rnuc;
    GAM   = soil_surface.GAM;
    Tch   = (Ta + .1) * ones(nl, 1);       % [C] shaded leaf T first guess
    Tcu   = (Ta + .3) * ones(size(Rnuc)); % [C] sunlit leaf T first guess
    ech   = ea * ones(nl, 1);            % [hPa] shaded leaf vapour first guess
    ecu   = ea * ones(size(Rnuc));      % [hPa] sunlit leaf vapour first guess
    Cch   = Ca * ones(nl, 1);            % shaded leaf CO2 first guess
    Ccu   = Ca * ones(size(Rnuc));      % sunlit leaf CO2 first guess
    L     = -1;                       % [m] Monin-Obukhov length (stability)
    % --- Physical constants: MW water/air, cp, g, Stefan-Boltzmann, ... ---
    phys = equations.physicalConstants();
    Oa    = get_struct_field(meteo, 'Oa', phys.Oa_mmolmol_default);
    O2a_Pa = Oa .* 1e-3 .* max(p - ea .* 100, 1);

    MH2O  = phys.MW_wv;
    Mair  = phys.MW_da;
    rhoa_std = phys.Pstd ./ (phys.R_da_SI .* (phys.Tr + phys.C_to_K));
    cp    = phys.Cp_air_SI;
    g     = phys.g_SI;
    kappa = phys.k_VK;
    sigmaSB = phys.sigma;
    Ps    = gap.Ps;
    nl    = canopy.nlayers;

    SoilHeatMethod = options.soil_heat_method;
    if ~(options.simulation == 1)
        SoilHeatMethod = 2;
    end

    kV   = canopy.kV;
    xl   = canopy.xl;

    % --- Geometry / unit helpers: e->q, air density, sun/shade weights ---
    e_to_q          = 100 * MH2O / Mair ./ p;       %           Conversion of vapour pressure [hPa] to absolute humidity [kg kg-1]
    qa_air          = max(0, ea .* e_to_q);
    Tv_air          = (Ta + phys.C_to_K) .* (1 + 0.61 .* qa_air);
    rhoa            = p ./ (phys.R_da_SI .* Tv_air); % CLM-style forcing-state air density [kg m-3]
    rhoa(~isfinite(rhoa) | rhoa <= 0) = rhoa_std;
    Fs              = [1 - Ps(end), Ps(end)];      %           Matrix containing values for 1-Ps and Ps of soil
    Fc              = (1 - Ps(1:end - 1))' / nl;      %           Matrix containing values for Ps of canopy



    if isfield(canopy, 'layer_weight') && ~isempty(canopy.layer_weight)
        cum_lai_top = [0; cumsum(canopy.layer_weight(1:end-1))];  % [nl x 1] cumulative LAI above each layer top
        fVh = exp(-kV * cum_lai_top);
    else
        fVh = exp(kV * xl(1:end - 1));   % no LAI profile: equal-layer decay
    end
    fVu             = ones(13, 36, nl);

    for i = 1:nl
        fVu(:, :, i) = fVh(i);
    end

    LAI = canopy.LAI;

    hydraulic = struct( ...
        'Tp_sun', 0, 'Tp_sha', 0, 'Tp_t', 0, ...
        'Tp_sun_demand', 0, 'Tp_sha_demand', 0, 'Tp_demand', 0, ...
        'supply_ratio', 1, ...
        'psil_sun', plant_hydro.init.psil_sun, ...
        'psil_sha', plant_hydro.init.psil_sha, ...
        'psil', 0.5 * (plant_hydro.init.psil_sun + plant_hydro.init.psil_sha), ...
        'psir', plant_hydro.init.psir, ...
        'psis', plant_hydro.init.psir, ...
        'wt_sun', 1, 'wt_sha', 1, ...
        'wt_sun_layer', ones(nl, 1), 'wt_sha_layer', ones(nl, 1), ...
        'gsr_total', 0, ...
        'Kplant_diff_sun', 0, ...
        'Kplant_diff_sha', 0, ...
        'Kplant_diff_flux_weighted', 0, ...
        'Kplant_app_sun', NaN, ...
        'Kplant_app_sha', NaN, ...
        'Kplant_app_total', NaN, ...
        'Kwhole_app_total', NaN, ...
        'Ksoilroot_app_total', NaN, ...
        'rwu_diag', struct());
    sfactor_sha = hydraulic.wt_sha;
    sfactor_sun = hydraulic.wt_sun;
    PSIs = soil_state.h;
    % Surface soil water potential: cm water column -> J/kg (psi = g*h, cm to m)
    PSIss = PSIs(NL, 1) * 9.81 / 100;
    if isfield(canopy, 'LAI_layer') && ~isempty(canopy.LAI_layer)
        lai_layer_for_sun = canopy.LAI_layer(:);
    else
        lai_layer_for_sun = ones(nl, 1) * LAI / nl;
    end
    % Solar altitude from prepareScopeStep (angles.tts = zenith angle)
    if isfield(angles, 'tts') && ~isempty(angles.tts)
        sun_altitude_deg = 90 - angles.tts;
    else
        sun_altitude_deg = 90;
    end
    sunlit_LAI_eff = sum(lai_layer_for_sun .* max(0, Ps(1:nl)));
    low_sun_skip = (sun_altitude_deg < 5) || (sunlit_LAI_eff < 0.05);
    skip_sunlit_ebal = low_sun_skip;
    sunlit_angle_fallback = false;

    % --- Soil surface resistance rss: drier surface -> larger resistance ---
    Theta_surf_now = soil_state.theta_liquid(NL, 2);
    soil_surface.rss = equations.Calculate_Soil_Resistance( ...
        Theta_surf_now, soil_params.theta_r(NL), soil_params.theta_s(NL), soil_params.field_capacity(NL));

    % Resistance model uses RTMo direct-beam k (same radiative geometry)
    if isfield(gap, 'k') && ~isempty(gap.k) && isfinite(gap.k)
        kb_direct = gap.k;
        canopy.Kb_direct = kb_direct;                                             % SCOPE RTMo direct-beam k
        canopy.Kb_source = 'SCOPE_RTMo_gap_k';
    end

    % --- Soil thermal IR: interpolate dry/wet reflectance by surface saturation ---
    rho_thermal_soil_dry = get_struct_field(soil_params, 'rho_thermal_dry', 0.086);
    rho_thermal_soil_wet = get_struct_field(soil_params, 'rho_thermal_wet', 0.020);
    Sat_surf = max(0, min(1, Theta_surf_now ./ max(soil_params.theta_s(NL), 1e-6)));
    soil_surface.rs_thermal_static = soil_surface.rs_thermal;
    bareRsThermal = rho_thermal_soil_dry + ...
        (rho_thermal_soil_wet - rho_thermal_soil_dry) .* Sat_surf;
    soil_surface.rs_thermal_bare = bareRsThermal;
    snowFraction = 0;
    snowRsThermal = bareRsThermal;
    if isfield(soil_surface, 'snow') && isstruct(soil_surface.snow)
        snowFraction = min(max(get_struct_field( ...
            soil_surface.snow, 'cover_fraction', 0), 0), 1);
        snowRsThermal = 1 - get_struct_field( ...
            soil_surface.snow, 'emissivity', 1 - bareRsThermal);
    end
    soil_surface.rs_thermal_bare = bareRsThermal;
    soil_surface.rs_thermal_snow = snowRsThermal;
    if ~isfield(soil_surface, 'refl_bare')
        soil_surface.refl_bare = soil_surface.refl;
    end
    if ~isfield(soil_surface, 'refl_snow')
        soil_surface.refl_snow = soil_surface.refl;
    end
    rtSurface = soil_surface;
    thermalIndex = find(spectral.wlS == 10000, 1);
    if ~isempty(thermalIndex)
        rtSurface.refl(thermalIndex) = (1 - snowFraction) .* bareRsThermal ...
            + snowFraction .* snowRsThermal;
    end
    soil_surface.rs_thermal = (1 - snowFraction) .* bareRsThermal ...
        + snowFraction .* snowRsThermal;
    soil_surface.epsilon_s = 1 - soil_surface.rs_thermal;
    soil_surface.Sat_surf = Sat_surf;
    bareEmissivity = 1 - bareRsThermal;
    snowEmissivity = 1 - snowRsThermal;
    snowTemperature = Fs * Ts;
    if isfield(soil_surface, 'snow') && isstruct(soil_surface.snow) ...
            && isfield(soil_surface.snow, 'temperature_C') ...
            && ~isempty(soil_surface.snow.temperature_C)
        snowTemperature = soil_surface.snow.temperature_C;
    end

    % Soil boundary resistance rbs:
    %   switch off: keep input; switch on: update each step from Ts-Ta and near-soil wind
    soil_surface.rbs_static = soil_surface.rbs;
    soil_surface.rbs_init = soil_surface.rbs;
    soil_surface.rbs_scheme = 'input';

    % --- Roughness zo and displacement d: update with LAI/SAI ---
    if ~isfield(canopy, 'SAI') || isempty(canopy.SAI)
        canopy.SAI = 0.5;        % deciduous orchard default (matches config_vegetation.txt)
    end
    canopy.SAI = max(canopy.SAI, 0);
    canopy.L_plus_S = max(canopy.LAI + canopy.SAI, 1e-6);
    canopy.zo_static = canopy.zo;                                            % Save initial zo from loadTimeSeries.
    canopy.d_static  = canopy.d;                                             % Save initial displacement height.
    [canopy.zo, canopy.d] = equations.zo_and_d(soil_surface, canopy);

    %% ========================================================================
    %%  [C] Energy-balance main loop
    %%  Repeat: radiation -> resistances -> photosynthesis -> heat fluxes
    %%  -> plant hydraulics -> residual -> update T
    %%  until Rn ≈ H + lE (+G) or max iterations
    %% ========================================================================

    aaEbal = [];
    aaEbalAccept = 0;
    aaEbalReject = 0;
    aaEbalOpts = struct( ...
        'enabled', use_ebal_anderson && ~skip_sunlit_ebal, ...
        'memory', 3, ...
        'damping', 0.4, ...
        'regularization', 1e-8, ...
        'minSingularRatio', 1e-8, ...
        'maxGammaNorm', 10, ...
        'maxStepNorm', 8, ...
        'minHistory', 3);

    while CONT                          % continue while not balanced and under maxit

        % ------------------------------------------------------------------
        %  [C1] Net radiation (update thermal IR and Rn from current leaf/soil T)
        % ------------------------------------------------------------------
        effectiveGroundT = surface.radiation.effectiveGroundTemperature( ...
            Ts, snowTemperature, snowFraction, ...
            bareEmissivity, snowEmissivity);
        rad = surface.radiation.RTMt_sb(spectral, rad, rtSurface, leafopt, ...
            canopy, gap, angles, Tcu, Tch, ...
            effectiveGroundT(2), effectiveGroundT(1), 1);

        Rnhct = rad.Rnhct;
        Rnuct = rad.Rnuct;
        Rnhst = rad.Rnhst;
        Rnust = rad.Rnust;

        Rnhc = rad.Rnhc;
        Rnuc = rad.Rnuc;
        Rnhs = rad.Rnhs;
        Rnus = rad.Rnus;

        Rnch        = Rnhc + Rnhct;             % [W m-2] shaded-leaf net radiation
        Rncu        = Rnuc + Rnuct;             % [W m-2] sunlit-leaf net radiation
        tileRadiation = surface.radiation.partitionGroundRadiation( ...
            spectral, rad, soil_surface, Ts, snowTemperature);
        RnsBare = tileRadiation.bare_net_W_m2;
        RnsSnow = tileRadiation.snow_net_W_m2;
        Rns = (1 - snowFraction) .* RnsBare ...
            + snowFraction .* RnsSnow;

        % ------------------------------------------------------------------
        %  [C2] Aerodynamic resistances (wind, stability, roughness -> rah/raw/...)
        % ------------------------------------------------------------------

        resist_in.u   = max(meteo.u, .2);%0.2
        resist_in.L   = L;
        resist_in.LAI = canopy.LAI;
        resist_in.rbs = soil_surface.rbs;
        resist_in.rss = soil_surface.rss;
        resist_in.rwc = canopy.rwc;
        resist_in.zo  = canopy.zo;
        resist_in.d   = canopy.d;
        resist_in.z   = meteo.z;
        resist_in.hc  = canopy.hc;
        resist_in.w   = canopy.leafwidth;
        resist_in.Cd  = canopy.Cd;
        if isfield(canopy_params, 'lai_max') && ~isempty(canopy_params.lai_max)
            resist_in.LAI_max = canopy_params.lai_max;
        end
        if isfield(canopy_params, 'w_leaf_max_m') && ~isempty(canopy_params.w_leaf_max_m)
            resist_in.w_leaf_max = canopy_params.w_leaf_max_m;
        end
        if isfield(canopy_params, 'crop_code') && ~isempty(canopy_params.crop_code)
            resist_in.CropCode = canopy_params.crop_code;
        end
        % Tao = RTMo direct-beam k; CropCode==1 uses user leaf width
        if isfield(canopy, 'Kb_direct') && ~isempty(canopy.Kb_direct)
            resist_in.Tao = canopy.Kb_direct;
            resist_in.Kb_direct = canopy.Kb_direct;
        end

        resist_out  = vegetation.energy.resistancesClm45(resist_in);
        if isfield(options, 'calc_rss_rbs') && options.calc_rss_rbs
            % Optional: Kustas-Norman soil boundary resistance
            soil_surface.rbs = equations.calcSoilBoundaryResistanceKustas(Ts, Ta, resist_out.uz0);
            soil_surface.rbs_scheme = 'kustas_norman_1999';
            resist_in.rbs = soil_surface.rbs;
            resist_out = vegetation.energy.resistancesClm45(resist_in);
        end

        ustar = resist_out.ustar;
        raa   = resist_out.raa;
        rah   = resist_out.rah;
        raw   = resist_out.raw;
        rawc  = resist_out.rawc;
        raws  = resist_out.raws;

        % ------------------------------------------------------------------
        %  [C3] Photosynthesis + stomata (sun and shade -> stomatal resistance rcw)
        % ------------------------------------------------------------------
        biochem_in.Fluorescence_model = options.Fluorescence_model;
        biochem_in.Type         = leafbio.Type;
        biochem_in.p            = p;
        biochem_in.Oa           = Oa;
        biochem_in.O2a_Pa       = O2a_Pa;
        biochem_in.Rdparam      = leafbio.Rdparam;
        biochem_in.tempcor      = options.apply_T_corr;

        % --- Shaded-leaf photosynthesis ---
        biochem_in = apply_phase3_biochem_fields(biochem_in, leafbio, canopy, domain, xyt, k);
        biochem_in.T    = Tch;
        biochem_in.eb   = ech;
        biochem_in.Vcmo = fVh .* leafbio.Vcmo;
        biochem_in.Cs   = Cch;
        biochem_in.Q    = rad.Pnh_Cab * 1E6;
        biochem_out = vegetation.biochemical(biochem_in, sfactor_sha, []);
        biochem_diag_sha = biochem_out;

        Ah   = biochem_out.A;
        Ahh  = biochem_out.Ag;
        Cih  = biochem_out.Ci;
        Fh   = biochem_out.eta;
        rcwh = biochem_out.rcw;
        qEh  = biochem_out.qE;
        Knh  = biochem_out.Kn;

        % --- Sunlit-leaf photosynthesis ---
        biochem_in = apply_phase3_biochem_fields(biochem_in, leafbio, canopy, domain, xyt, k);
        biochem_in.T    = Tcu;
        biochem_in.eb   = ecu;
        biochem_in.Vcmo = fVu .* leafbio.Vcmo;
        biochem_in.Cs   = Ccu;
        biochem_in.Q    = rad.Pnu_Cab * 1E6;
        biochem_out = vegetation.biochemical(biochem_in, sfactor_sun, []);
        biochem_diag_sun = biochem_out;

        Au                  = biochem_out.A;
        Auu                 = biochem_out.Ag;   % gross photosynthesis (for GPP)
        Ciu                 = biochem_out.Ci;
        Fu                  = biochem_out.eta;
        rcwu                = biochem_out.rcw;
        qEu                 = biochem_out.qE;
        Knu                 = biochem_out.Kn;

        Pinh                = rad.Pnh;
        Pinu                = rad.Pnu;
        Pinh_Cab            = rad.Pnh_Cab;
        Pinu_Cab            = rad.Pnu_Cab;
        Rnh_PAR             = rad.Rnh_PAR;
        Rnu_PAR             = rad.Rnu_PAR;

        % ------------------------------------------------------------------
        %  [C4] Sensible/latent heat + plant hydraulics
        %  Compute H, lE; limit dry-leaf demand by water supply; update stomata
        % ------------------------------------------------------------------
        rss  = soil_surface.rss;
        rac_h   = (LAI + 1) * (rah + rawc);
        rac_w   = (LAI + 1) * (raw + rawc);
        ras_h   = (LAI + 1) * (rah + raws);
        ras_w   = (LAI + 1) * (raw + raws);
        rac     = rac_h;
        % This path: fluxes use above-canopy forcing; canopy-air profiles are diagnostic
        ra_leaf_h = rac_h;
        ra_leaf_w = rac_w;
        ra_soil_h = ras_h;
        ra_soil_w = ras_w;
        hydraulic_inner_maxit = 6;
        hydraulic_wt_tol = 1e-3;
        for hyd_iter = 1:hydraulic_inner_maxit
            air_source = struct( ...
                'H_layer_ground',  zeros(nl, 1), ...
                'LE_layer_ground', zeros(nl, 1), ...
                'A_layer_ground',  zeros(nl, 1), ...
                'H_soil_ground',   0, ...
                'LE_soil_ground',  0, ...
                'e_to_q', e_to_q, ...
                'lambda', (2.501 - 0.002361 * Ta) * 1E6);
            if exist('Hch', 'var') && exist('Hcu', 'var') && exist('Hs', 'var') ...
                    && exist('lEch', 'var') && exist('lEcu', 'var') && exist('lEs', 'var')
                Hcu1d_air_prev  = equations.meanleaf(canopy, Hcu,  'angles');
                lEcu1d_air_prev = equations.meanleaf(canopy, lEcu, 'angles');
                if isfield(canopy, 'LAI_layer') && ~isempty(canopy.LAI_layer)
                    LAI_layer_air = canopy.LAI_layer(:);
                else
                    LAI_layer_air = ones(nl, 1) * LAI / nl;
                end
                air_source.H_layer_ground  = LAI_layer_air .* ((1 - Ps(1:nl)) .* Hch(:)  + Ps(1:nl) .* Hcu1d_air_prev(:));
                air_source.LE_layer_ground = LAI_layer_air .* ((1 - Ps(1:nl)) .* lEch(:) + Ps(1:nl) .* lEcu1d_air_prev(:));
                air_source.H_soil_ground   = Fs * Hs;
                air_source.LE_soil_ground  = Fs * lEs;
                if exist('Ahh', 'var') && exist('Auu', 'var')
                    Auu1d_air_prev = equations.meanleaf(canopy, Auu, 'angles');
                    air_source.A_layer_ground = LAI_layer_air .* ...
                        ((1 - Ps(1:nl)) .* max(0, Ahh(:)) + Ps(1:nl) .* max(0, Auu1d_air_prev(:)));
                end
            end

            % Heat fluxes use classic forcing; diagnose canopy-air after EB converges
            for air_iter = 1:1
                Ta_ebal_layer = Ta * ones(nl, 1);
                ea_ebal_layer = ea * ones(nl, 1);
                Ta_ebal_layer_3d = Ta * ones(size(Tcu));
                ea_ebal_layer_3d = ea * ones(size(Tcu));

                % Fluxes use above-canopy forcing; leaf psi in Kelvin term (MPa->J/kg)
                % Dry leaf = stomatal T (can go to RWU); wet leaf = interception E (no RWU)
                Ca_leaf_sha = Ca;
                Ca_leaf_sun = Ca;
                PSI = hydraulic.psil * 1000;
                [lEch, Hch, ech, Cch, lambdah, sh, lEch_dry, lEch_wet] = vegetation.energy.heatFluxesClm45( ...
                    ra_leaf_h, rcwh, Tch, ea_ebal_layer, Ta_ebal_layer, e_to_q, PSI, Ca_leaf_sha, Cih, ra_leaf_w, rhoa, fwet_sha_canopy);
                [lEcu, Hcu, ecu, Ccu, lambdau, su, lEcu_dry, lEcu_wet] = vegetation.energy.heatFluxesClm45( ...
                    ra_leaf_h, rcwu, Tcu, ea_ebal_layer_3d, Ta_ebal_layer_3d, e_to_q, PSI, Ca_leaf_sun, Ciu, ra_leaf_w, rhoa, fwet_sun_canopy);
                groundSurface = struct( ...
                    'temperature_C', Ts, ...
                    'bare_temperature_C', Ts, ...
                    'snow_temperature_C', snowTemperature, ...
                    'soil_water_potential_J_kg', PSIss, ...
                    'soil_resistance_s_m', rss);
                if isfield(soil_surface, 'snow')
                    groundSurface.snow = soil_surface.snow;
                end
                groundAtmosphere = struct( ...
                    'ra_heat_s_m', ra_soil_h, ...
                    'ra_water_s_m', ra_soil_w, ...
                    'vapor_pressure_hPa', ea, ...
                    'air_temperature_C', Ta, ...
                    'vapor_to_humidity_hPa_inverse', e_to_q, ...
                    'air_density_kg_m3', rhoa);
                groundFlux = vegetation.energy.groundSurfaceFluxes( ...
                    groundSurface, groundAtmosphere);
                lEs = groundFlux.ground_upward_latent_heat_flux_W_m2;
                Hs = groundFlux.ground_upward_sensible_heat_flux_W_m2;

                Hcu1d_air  = equations.meanleaf(canopy, Hcu,  'angles');
                lEcu1d_air = equations.meanleaf(canopy, lEcu, 'angles');
                if isfield(canopy, 'LAI_layer') && ~isempty(canopy.LAI_layer)
                    LAI_layer_air = canopy.LAI_layer(:);
                else
                    LAI_layer_air = ones(nl, 1) * LAI / nl;
                end
                air_source.H_layer_ground  = LAI_layer_air .* ((1 - Ps(1:nl)) .* Hch(:)  + Ps(1:nl) .* Hcu1d_air(:));
                air_source.LE_layer_ground = LAI_layer_air .* ((1 - Ps(1:nl)) .* lEch(:) + Ps(1:nl) .* lEcu1d_air(:));
                air_source.H_soil_ground   = Fs * Hs;
                air_source.LE_soil_ground  = Fs * lEs;
            end

            % Canopy-integrated H (same LAI weights as final output; keeps L consistent)
            Hstot       = Fs * Hs;
            Hcu1d_loop  = equations.meanleaf(canopy, Hcu, 'angles');
            if isfield(canopy, 'layer_weight') && ~isempty(canopy.layer_weight)
                Hctot = sum(canopy.LAI_layer .* ((1 - Ps(1:nl)) .* Hch + Ps(1:nl) .* Hcu1d_loop));
            else
                Hctot = LAI * (Fc * Hch + sum(Ps(1:nl) .* Hcu1d_loop) / nl);
            end
            Htot        = Hstot + Hctot;
            lambda_ebal = (2.501 - 0.002361 * Ta) * 1E6;
            if isfinite(lambda_ebal) && lambda_ebal > 0
                % Hydraulics use dry-leaf LE only -> transpiration demand
                % Wet-leaf interception evaporation does not request root water
                lEcu1d_hyd     = equations.meanleaf(canopy, lEcu_dry, 'angles');
                if isfield(canopy, 'layer_weight') && ~isempty(canopy.layer_weight)
                    LAI_layer_hyd = canopy.LAI_layer(:);
                    layer_weight_hyd = canopy.layer_weight(:);
                    LE_sha_layer_hyd = (1 - Ps(1:nl)) .* lEch_dry(:);
                    LE_sun_layer_hyd = Ps(1:nl) .* lEcu1d_hyd(:);
                    LE_sha_for_hyd = sum(layer_weight_hyd .* LE_sha_layer_hyd);
                    LE_sun_for_hyd = sum(layer_weight_hyd .* LE_sun_layer_hyd);
                    Tp_sha_layer_cm_s = LAI_layer_hyd .* max(LE_sha_layer_hyd, 0) / lambda_ebal / 1000 * 100;
                    Tp_sun_layer_cm_s = LAI_layer_hyd .* max(LE_sun_layer_hyd, 0) / lambda_ebal / 1000 * 100;
                else
                    LE_sha_for_hyd = Fc * lEch_dry;
                    LE_sun_for_hyd = sum(Ps(1:nl) .* lEcu1d_hyd) / nl;
                    Tp_sha_layer_cm_s = ones(nl, 1) * LAI / nl .* max((1 - Ps(1:nl)) .* lEch_dry(:), 0) / lambda_ebal / 1000 * 100;
                    Tp_sun_layer_cm_s = ones(nl, 1) * LAI / nl .* max(Ps(1:nl) .* lEcu1d_hyd(:), 0) / lambda_ebal / 1000 * 100;
                    layer_weight_hyd = ones(nl, 1) / nl;
                end
                if ~isreal(LE_sha_for_hyd) || ~isreal(LE_sun_for_hyd)
                    error('PlantHydraulic:ComplexLE', ...
                        ['LE is complex before hydraulic transpiration conversion. ', ...
                         'isreal(lEch)=%d, isreal(lEcu)=%d, isreal(Tch)=%d, isreal(Tcu)=%d, ', ...
                         'isreal(rcwh)=%d, isreal(rcwu)=%d, isreal(Cih)=%d, isreal(Ciu)=%d, ', ...
                         'isreal(rac)=%d, isreal(raa)=%d, isreal(rawc)=%d, isreal(e_to_q)=%d, isreal(ea)=%d.'], ...
                        isreal(lEch), isreal(lEcu), isreal(Tch), isreal(Tcu), ...
                        isreal(rcwh), isreal(rcwu), isreal(Cih), isreal(Ciu), ...
                        isreal(rac), isreal(raa), isreal(rawc), isreal(e_to_q), isreal(ea));
                end
                % Hydraulic topology: shared stem + sun/shade leaves; optional layered stress
                hyd_mode = 'aggregate';
                if isfield(options, 'use_layer_hydraulic_stress') && options.use_layer_hydraulic_stress
                    hyd_mode = 'LAI_proportional';
                end
                hydraulic_iter = vegetation.hydraulics.solveEbalHydraulicsLayered( ...
                    plant_hydro, soil_state, soil_conductivity, kt, canopy.hc, ...
                    Tp_sun_layer_cm_s, Tp_sha_layer_cm_s, layer_weight_hyd, hyd_mode, soil_params);
                hydraulic = hydraulic_iter;
                if hydraulic.supply_ratio < 0.05
                    maxit = max(maxit, 900);
                elseif hydraulic.supply_ratio < 0.2
                    maxit = max(maxit, 300);
                else
                    maxit = maxit_base;
                end
                sfactor_sha_old = sfactor_sha;
                sfactor_sun_old = sfactor_sun;
                sfactor_sha = hydraulic.wt_sha;
                sfactor_sun = hydraulic.wt_sun;

                wt_delta = max([abs(sfactor_sha(:) - sfactor_sha_old(:)); abs(sfactor_sun(:) - sfactor_sun_old(:))]);
                if wt_delta <= hydraulic_wt_tol || hyd_iter == hydraulic_inner_maxit
                    break;
                end

                % Hydraulic inner loop: update stress wt then stomata/photosynthesis
                % --- Shaded leaves ---
                biochem_in = apply_phase3_biochem_fields(biochem_in, leafbio, canopy, domain, xyt, k);
                biochem_in.T    = Tch;
                biochem_in.eb   = ech;
                biochem_in.Vcmo = fVh .* leafbio.Vcmo;
                biochem_in.Cs   = Cch;
                biochem_in.Q    = rad.Pnh_Cab * 1E6;
                sfactor_sha_bio = hydraulic.wt_sha_layer;
                biochem_out = vegetation.biochemical(biochem_in, sfactor_sha_bio, []);
                biochem_diag_sha = biochem_out;

                Ah   = biochem_out.A;
                Ahh  = biochem_out.Ag;
                Cih  = biochem_out.Ci;
                Fh   = biochem_out.eta;
                rcwh = biochem_out.rcw;
                qEh  = biochem_out.qE;
                Knh  = biochem_out.Kn;

                % --- Sunlit leaves ---
                biochem_in = apply_phase3_biochem_fields(biochem_in, leafbio, canopy, domain, xyt, k);
                biochem_in.T    = Tcu;
                biochem_in.eb   = ecu;
                biochem_in.Vcmo = fVu .* leafbio.Vcmo;
                biochem_in.Cs   = Ccu;
                biochem_in.Q    = rad.Pnu_Cab * 1E6;
                sfactor_sun_bio = reshape(hydraulic.wt_sun_layer, 1, 1, nl);
                biochem_out = vegetation.biochemical(biochem_in, sfactor_sun_bio, []);
                biochem_diag_sun = biochem_out;

                Au                  = biochem_out.A;
                Auu                 = biochem_out.Ag;
                Ciu                 = biochem_out.Ci;
                Fu                  = biochem_out.eta;
                rcwu                = biochem_out.rcw;
                qEu                 = biochem_out.qE;
                Knu                 = biochem_out.Kn;
                if wt_delta > hydraulic_wt_tol && hyd_iter < hydraulic_inner_maxit
                    continue;
                end
            end
                break  % Converged, or lambda_ebal is invalid
        end

        % --- Soil heat flux G ---
        if SoilHeatMethod == 2
            G = 0.30 * RnsBare;       % legacy option, per bare-soil tile
            dG_dTs = 0;
        else
            % force-restore: use soil surface temperature history Tsold
            G = GAM / sqrt(pi) * 2 * sum(([Ts'; Tsold(1:end - 1, :)] - Tsold) ./ Deltat .* (sqrt(x) - sqrt(x - Deltat)));
            G = G';
            dG_dTs = 2 * GAM / sqrt(pi * Deltat);
        end

        % ------------------------------------------------------------------
        %  [C5] Monin-Obukhov length L (from total H; next pass uses it in resistances)
        % ------------------------------------------------------------------
        L           = -rhoa * cp * ustar.^3 .* (Ta + 273.15) ./ (kappa * g * Htot);           % [1]
        L(L < -1E3)   = -1E3;                                                     % [1]
        L(L > 1E2)    =  1E2;                                                     % [1]
        L(isnan(L)) = -1;                                                       % [1]

        % ------------------------------------------------------------------
        %  [C6] Residuals: is Rn - H - lE (-G) small enough on each surface?
        % ------------------------------------------------------------------
        if skip_sunlit_ebal
            Rnu1d_skip = equations.meanleaf(canopy, Rncu, 'angles');
            if isfield(canopy, 'LAI_layer') && ~isempty(canopy.LAI_layer)
                LW_skip = canopy.LAI_layer(:);
            else
                LW_skip = ones(nl, 1) * LAI / nl;
            end
            sunlit_Rn_extra = sum(LW_skip .* Ps(1:nl) .* (Rnu1d_skip(:) - Rnch(:))) / max(sum(LW_skip), eps);
            Rnch_eff = Rnch(:) + sunlit_Rn_extra;
            Rnch = Rnch_eff;
            Tcu = repmat(reshape(Tch(:), 1, 1, nl), size(Tcu, 1), size(Tcu, 2), 1);
            Rncu = repmat(reshape(Rnch_eff(:), 1, 1, nl), size(Tcu, 1), size(Tcu, 2), 1);
            lEcu = repmat(reshape(lEch(:), 1, 1, nl), size(Tcu, 1), size(Tcu, 2), 1);
            Hcu = repmat(reshape(Hch(:), 1, 1, nl), size(Tcu, 1), size(Tcu, 2), 1);
        end
        EBerch      = Rnch - lEch - Hch;
        EBercu      = Rncu - lEcu - Hcu;
        EBers = RnsBare ...
            - groundFlux.bare_ground_upward_latent_heat_flux_tile_W_m2 ...
            - groundFlux.bare_ground_upward_sensible_heat_flux_tile_W_m2 - G;
        if snowFraction >= 1 - 1e-12
            EBers(:) = 0;
        end

        counter     = counter + 1;
        % Sunlit residual: angle-weighted per layer, then max over layers
        % maxEBercu_perangle is diagnostic only
        EBercu_layer        = equations.meanleaf(canopy, EBercu, 'angles');
        maxEBercu           = max(abs(EBercu_layer));
        maxEBercu_perangle  = max(abs(EBercu(:)));
        maxEBerch   = max(abs(EBerch));
        maxEBers    = max(abs(EBers));

        CONT        = (maxEBercu >   tolEBer    | ...
                       maxEBerch >   tolEBer    | ...
                       maxEBers  >   tolEBer)    & ...
                        counter   <   maxit + 1;

        % Guard: layer residual OK but angle T is nonphysical -> layer-mean sunlit T
        % (not a simple temperature clamp)
        detailed_ebal_closed = max([maxEBercu, maxEBerch, maxEBers]) <= tolEBer;
        sunlit_temperature_invalid = any(~isfinite(Tcu(:))) || ...
            any(Tcu(:) < -60) || any(Tcu(:) > 90) || any(abs(Tcu(:) - Ta) > 60);
        if ~CONT && detailed_ebal_closed && ~skip_sunlit_ebal && sunlit_temperature_invalid
            skip_sunlit_ebal = true;
            sunlit_angle_fallback = true;
            aaEbal = [];
            aaEbalOpts.enabled = false;
            Tcu = repmat(reshape(Tch(:), 1, 1, nl), size(Tcu, 1), size(Tcu, 2), 1);
            CONT = true;
            continue
        end
        if counter == 5
            Wc = 0.6;
        end
        if counter == 10
            Wc = 0.4;
        end
        if counter == 30
            Wc = 0.1;
        end

        % ------------------------------------------------------------------
        %  [C7] Update temperatures (Picard under-relaxation + optional Anderson)
        %  Only if another iteration will run (keep returned T consistent with fluxes)
        % ------------------------------------------------------------------
        if CONT
            Tch_size = size(Tch);
            Tcu_size = size(Tcu);
            Ts_size = size(Ts);
            nTch = numel(Tch);
            nTcu = numel(Tcu);
            nTs = numel(Ts);
            xEbalOld = [Tch(:); Tcu(:); Ts(:)];

            Tch_picard = Tch + Wc * EBerch ./ ((rhoa * cp) ./ ra_leaf_h + rhoa * lambdah * e_to_q .* sh ./ (ra_leaf_w + rcwh) + 4 * leafbio.emis * sigmaSB * (Tch + 273.15).^3);
            Tcu_picard = Tcu + Wc * EBercu ./ ((rhoa * cp) ./ ra_leaf_h + rhoa * lambdau * e_to_q .* su ./ (ra_leaf_w + rcwu) + 4 * leafbio.emis * sigmaSB * (Tcu + 273.15).^3);
            Ts_picard  = Ts  + Wc * EBers  ./ (rhoa * cp ./ ra_soil_h ...
                + groundFlux.dlatent_bare_tile_dT_W_m2_K ...
                + 4 * bareEmissivity * sigmaSB ...
                * (Ts + 273.15).^3 + dG_dTs);
            if snowFraction >= 1 - 1e-12
                Ts_picard = Ts;
            end
            Tch_picard(abs(Tch_picard) > 100) = Ta;
            Tcu_picard(abs(Tcu_picard) > 100) = Ta;
            Ts_picard(abs(Ts_picard) > 100) = Ta;

            xEbalPicard = [Tch_picard(:); Tcu_picard(:); Ts_picard(:)];
            [xEbalNext, aaEbal] = solver.andersonMixUpdate(xEbalOld, xEbalPicard, aaEbal, aaEbalOpts);

            useCandidate = isfield(aaEbal, 'lastAccepted') && aaEbal.lastAccepted && ...
                all(isfinite(xEbalNext)) && all(abs(xEbalNext) <= 100);
            if useCandidate
                Tch = reshape(xEbalNext(1:nTch), Tch_size);
                Tcu = reshape(xEbalNext(nTch + (1:nTcu)), Tcu_size);
                Ts = reshape(xEbalNext(nTch + nTcu + (1:nTs)), Ts_size);
                aaEbalAccept = aaEbalAccept + 1;
            else
                Tch = Tch_picard;
                Tcu = Tcu_picard;
                Ts = Ts_picard;
                aaEbalReject = aaEbalReject + 1;
            end
        end
        if any(isnan(Tch)) || any(isnan(Tcu(:)))
            warning('Canopy temperature gives NaNs');
        end
        if any(isnan(Ts))
            warning('Soil temperature gives NaNs');
        end

    end

    iteration.counter = counter;
    iteration.ebalAA_accept = aaEbalAccept;
    iteration.ebalAA_reject = aaEbalReject;
    iteration.lowSunSkip = low_sun_skip;
    iteration.sunlitAngleFallback = sunlit_angle_fallback;
    iteration.sunlitEbalCollapsed = skip_sunlit_ebal;
    iteration.sunlit_LAI_eff = sunlit_LAI_eff;
    iteration.sun_altitude_deg = sun_altitude_deg;
    iteration.maxEBercu = maxEBercu;
    iteration.maxEBercu_perangle = maxEBercu_perangle;
    iteration.maxEBerch = maxEBerch;
    iteration.maxEBers = maxEBers;
    iteration.maxEBer = max([maxEBercu, maxEBerch, maxEBers]);
    ebal_converged = iteration.maxEBer <= tolEBer;
    profiles.etah = Fh;
    profiles.etau = Fu;

    if SoilHeatMethod < 2 && ebal_converged
        Tsold(2:end, :) = soil_surface.Tsold(1:end - 1, :);
        Tsold(1, :)  = Ts(:);
        if isnan(Ts)
            Tsold(1, :) = Tsold(2, :);
        end
        if isreal(Ts)
            soil_surface.Tsold = Tsold;
        else
            Tsold(1, :) = Tsold(2, :);
            soil_surface.Tsold = Tsold;
        end
    end

    Tbr         = (rad.Eoutte / sigmaSB)^0.25;
    Lot_        = equations.Planck(spectral.wlS', Tbr);
    rad.LotBB_  = Lot_;           % Note that this is the blackbody radiance!

    %% ========================================================================
    %%  [D] Non-convergence warning (may be deferred by runStepEbal)
    %% ========================================================================
    deferNonConvergenceWarning = logical(control.defer_nonconvergence_warning);
    if ~deferNonConvergenceWarning && counter >= maxit && iteration.maxEBer > tolEBer
        warning(['\n warning: maximum number of iteratations exceeded', ...
                 '\n Energy balance error sunlit vegetation = %4.2f W m-2 ', ...
                 '\n Energy balance error shaded vegetation = %4.2f W m-2 ', ...
                 '\n Energy balance error soil              = %4.2f W m-2 '], maxEBercu, maxEBerch, maxEBers);

    end

    %% ========================================================================
    %%  [E] Optional vertical profiles (layer T/fluxes/hydraulics/canopy air)
    %%  Canopy-air fields here are diagnostic only (not fed back into EB)
    %% ========================================================================
    if io.canopyLayerProfilesEnabled(options)
        [Hcu1d]           = equations.meanleaf(canopy, Hcu,          'angles');   % [nl] mean sensible heat, sunlit leaves
        [lEcu1d]          = equations.meanleaf(canopy, lEcu,         'angles');   % [nl] mean latent heat, sunlit leaves
        [Au1d]            = equations.meanleaf(canopy, Au,           'angles');   % [nl] mean photosynthesis, sunlit leaves
        [Fu_Pn1d]         = equations.meanleaf(canopy, Fu .* Pinu_Cab, 'angles'); % [nl] mean fluorescence-weighted APAR, sunlit leaves
        [qEuL]            = equations.meanleaf(canopy, qEu,          'angles');   % [nl] mean fluorescence yield, sunlit leaves
        [Rnu1d]           = equations.meanleaf(canopy, Rncu,         'angles');   % [nl] mean net radiation, sunlit leaves
        [Tcu1d]           = equations.meanleaf(canopy, Tcu,          'angles');   % [nl] mean temperature, sunlit leaves
        [ecu1d]           = equations.meanleaf(canopy, ecu,          'angles');   % [nl] mean vapour pressure, sunlit leaves
        [Ccu1d]           = equations.meanleaf(canopy, Ccu,          'angles');   % [nl] mean CO2, sunlit leaves

        % Layer transpiration-demand diagnostics (optional for two-leaf hydraulics)
        if isfield(canopy, 'LAI_layer') && ~isempty(canopy.LAI_layer)
            LAI_layer_profile = canopy.LAI_layer(:);
        else
            LAI_layer_profile = ones(nl, 1) * LAI / nl;
        end
        if isfield(canopy, 'layer_weight') && ~isempty(canopy.layer_weight)
            layer_weight_profile = canopy.layer_weight(:);
        else
            layer_weight_profile = ones(nl, 1) / nl;
        end
        lambda_profile = (2.501 - 0.002361 * Ta) * 1E6;  % [J kg-1]
        Ps_profile = Ps(1:nl);
        lEch_profile = lEch(:);
        lEcu_profile = lEcu1d(:);
        lEsha1d = (1 - Ps_profile) .* lEch_profile;
        lEsun1d = Ps_profile .* lEcu_profile;
        profiles.layer_weight = layer_weight_profile;              % [-]      layer LAI fraction
        profiles.LAI_layer    = LAI_layer_profile;                 % [m2 m-2] layer LAI, top to bottom
        profiles.Ps1d         = Ps_profile;                        % [-]      sunlit fraction per layer
        profiles.lEch1d       = lEch_profile;                      % [W m-2 leaf] shaded-leaf LE before Ps weighting
        profiles.lEcu1d       = lEcu_profile;                      % [W m-2 leaf] sunlit-leaf LE before Ps weighting
        profiles.lEsha1d      = lEsha1d;                           % [W m-2 leaf] shaded contribution per layer
        profiles.lEsun1d      = lEsun1d;                           % [W m-2 leaf] sunlit contribution per layer
        profiles.Tp_sha_layer_demand = LAI_layer_profile .* max(lEsha1d, 0) / lambda_profile / 1000 * 100; % [cm s-1]
        profiles.Tp_sun_layer_demand = LAI_layer_profile .* max(lEsun1d, 0) / lambda_profile / 1000 * 100; % [cm s-1]
        profiles.Tp_layer_demand     = profiles.Tp_sha_layer_demand + profiles.Tp_sun_layer_demand;       % [cm s-1]
        profiles.Tp_sha_layer        = hydraulic.Tp_sha_layer;       % [cm s-1] hydraulic-limited shaded transpiration
        profiles.Tp_sun_layer        = hydraulic.Tp_sun_layer;       % [cm s-1] hydraulic-limited sunlit transpiration
        profiles.Tp_layer            = hydraulic.Tp_layer;           % [cm s-1] hydraulic-limited total layer transpiration
        profiles.wt_sha_layer        = hydraulic.wt_sha_layer;       % [-]      shaded hydraulic stress per layer
        profiles.wt_sun_layer        = hydraulic.wt_sun_layer;       % [-]      sunlit hydraulic stress per layer
        profiles.psil_sha_layer      = hydraulic.psil_sha_layer;     % [MPa]    shaded leaf water potential per layer
        profiles.psil_sun_layer      = hydraulic.psil_sun_layer;     % [MPa]    sunlit leaf water potential per layer
        profiles.layer_hydraulic_mode = hydraulic.layer_hydraulic_mode; % [-] layered hydraulic topology flag

        % Canopy-air profiles after convergence (inspection / future coupling)
        Auu1d_profile = equations.meanleaf(canopy, Auu, 'angles');
        air_source_profile = struct( ...
            'H_layer_ground',  LAI_layer_profile .* ((1 - Ps_profile) .* Hch(:)  + Ps_profile .* Hcu1d(:)), ...
            'LE_layer_ground', LAI_layer_profile .* ((1 - Ps_profile) .* lEch(:) + Ps_profile .* lEcu1d(:)), ...
            'A_layer_ground',  LAI_layer_profile .* ((1 - Ps_profile) .* max(0, Ahh(:)) + Ps_profile .* max(0, Auu1d_profile(:))), ...
            'H_soil_ground',   Fs * Hs, ...
            'LE_soil_ground',  Fs * lEs, ...
            'e_to_q', e_to_q, ...
            'lambda', lambda_profile);
        air_profile = equations.calc_canopy_air_profile( ...
            canopy, gap, resist_out, Ta, ea, Ca, Tch, Tcu1d, ech, ecu1d, Cch, Ccu1d, air_source_profile);
        profiles.Ta_layer       = air_profile.Ta_layer;       % [degC] layer canopy-air temperature
        profiles.ea_layer       = air_profile.ea_layer;       % [hPa]  layer canopy-air vapour pressure
        profiles.CO2_layer      = air_profile.CO2_layer;      % [umol mol-1] layer canopy-air CO2
        profiles.Kh_layer       = air_profile.Kh_layer;       % [m2 s-1] layer turbulent diffusivity
        profiles.z_layer_m      = air_profile.z_layer_m;      % [m] layer center height above ground
        profiles.depth_geom     = air_profile.depth_geom;     % [-] 0=canopy top, 1=bottom
        profiles.depth_lai      = air_profile.depth_lai;      % [-] cumulative LAI-weighted depth
        profiles.air_mix_source = air_profile.mix_air_source; % [-] diagnostic source relaxation strength
        profiles.Tleaf_layer    = air_profile.Tleaf_layer;    % [degC] weighted two-leaf temperature source
        profiles.eleaf_layer    = air_profile.eleaf_layer;    % [hPa] weighted two-leaf vapour source
        profiles.Cleaf_layer    = air_profile.Cleaf_layer;    % [umol mol-1] weighted two-leaf CO2 source
        profiles.canopy_air_profile_method = air_profile.method;

        profiles.Tchave     = mean(Tch);                                           % [1]               mean temp shaded leaves
        profiles.Tch        = Tch;                                                 % [nl]
        profiles.Tcu1d      = Tcu1d;                                               % [nl]
        profiles.ecu1d      = ecu1d;                                               % [nl]
        profiles.Ccu1d      = Ccu1d;                                               % [nl]
        profiles.Tc1d       = (1 - Ps(1:nl)) .* Tch       + Ps(1:nl) .* (Tcu1d);         % [nl]              mean temp leaves, per layer
        profiles.Hc1d       = (1 - Ps(1:nl)) .* Hch       + Ps(1:nl) .* (Hcu1d);         % [nl]              mean sens heat leaves, per layer
        profiles.lEc1d      = (1 - Ps(1:nl)) .* lEch      + Ps(1:nl) .* (lEcu1d);        % [nl]              mean latent heat leaves, per layer
        profiles.A1d        = (1 - Ps(1:nl)) .* Ah        + Ps(1:nl) .* (Au1d);          % [nl]              mean photos leaves, per layer
        profiles.F_Pn1d     = ((1 - Ps(1:nl)) .* Fh .* Pinh_Cab + Ps(1:nl) .* (Fu_Pn1d));  % [nl]           mean fluor leaves, per layer
        profiles.qE         = ((1 - Ps(1:nl)) .* qEh      + Ps(1:nl) .* (qEuL));         % [nl]           mean fluor leaves, per layer
        profiles.Rn1d       = ((1 - Ps(1:nl)) .* Rnch     + Ps(1:nl) .* (Rnu1d));        % [nl]
    end

    %% ========================================================================
    %%  [F] Integrate canopy totals, write fluxes/thermal, assemble RWU
    %% ========================================================================
    % --- Canopy integral: weight by LAI_layer if present, else uniform 1/nl ---

    Ps_nl = Ps(1:nl);   % [nl x 1] sunlit fraction per layer
    if isfield(canopy, 'layer_weight') && ~isempty(canopy.layer_weight)
        LW   = canopy.LAI_layer;    % [nl x 1] LAI per layer, sums to LAI
        WW   = canopy.layer_weight; % [nl x 1] normalised weight, sums to 1
        Rnu1d_f   = equations.meanleaf(canopy, Rncu,                      'angles');
        lEcu1d_f  = equations.meanleaf(canopy, lEcu,                      'angles');
        Hcu1d_f   = equations.meanleaf(canopy, Hcu,                       'angles');
        Auu1d_f   = equations.meanleaf(canopy, Auu,                       'angles');
        Tcu1d_f   = equations.meanleaf(canopy, Tcu,                       'angles');
        Pinu1d_f  = equations.meanleaf(canopy, Pinu,                      'angles');
        PinuC1d_f = equations.meanleaf(canopy, Pinu_Cab,                  'angles');
        RnuP1d_f  = equations.meanleaf(canopy, Rnu_PAR,                   'angles');
        etRnP1d_f = equations.meanleaf(canopy, profiles.etau .* Rnu_PAR,  'angles');
        % Wet-leaf interception evaporation accumulated separately (lE_wet_canopy)
        lEcu_wet1d_f = equations.meanleaf(canopy, lEcu_wet,               'angles');

        Rnctot       = sum(LW .* ((1-Ps_nl) .* Rnch         + Ps_nl .* Rnu1d_f));
        lEctot       = sum(LW .* ((1-Ps_nl) .* lEch         + Ps_nl .* lEcu1d_f));
        lE_wet_sha = sum(LW .* (1 - Ps_nl) .* lEch_wet);
        lE_wet_sun = sum(LW .* Ps_nl .* lEcu_wet1d_f);
        lE_dry_sha = sum(LW .* (1 - Ps_nl) .* lEch) - lE_wet_sha;
        lE_dry_sun = sum(LW .* Ps_nl .* lEcu1d_f) - lE_wet_sun;
        lE_wet_canopy_total = lE_wet_sha + lE_wet_sun;
        Hctot        = sum(LW .* ((1-Ps_nl) .* Hch          + Ps_nl .* Hcu1d_f));
        Actot        = sum(LW .* ((1-Ps_nl) .* Ahh          + Ps_nl .* Auu1d_f));
        % Leaf dark respiration Rdtot = Ag - A (layered branch)
        Au1d_f       = equations.meanleaf(canopy, Au,  'angles');
        Rdtot        = sum(LW .* ((1-Ps_nl) .* (Ahh - Ah)   + Ps_nl .* (Auu1d_f - Au1d_f)));
        Tcave        = sum(WW .* ((1-Ps_nl) .* Tch          + Ps_nl .* Tcu1d_f));
        Pntot        = sum(LW .* ((1-Ps_nl) .* Pinh         + Ps_nl .* Pinu1d_f));
        Pntot_Cab    = sum(LW .* ((1-Ps_nl) .* Pinh_Cab     + Ps_nl .* PinuC1d_f));
        Rntot_PAR    = sum(LW .* ((1-Ps_nl) .* Rnh_PAR      + Ps_nl .* RnuP1d_f));
        aPAR_Cab_eta = sum(LW .* ((1-Ps_nl) .* profiles.etah .* Rnh_PAR + Ps_nl .* etRnP1d_f));
    else
        Rnctot       = LAI * (Fc * Rnch      + equations.meanleaf(canopy, Rncu,                     'angles_and_layers', Ps));
        lEctot       = LAI * (Fc * lEch      + equations.meanleaf(canopy, lEcu,                     'angles_and_layers', Ps));
        % Wet-leaf interception evaporation (uniform-canopy branch)
        lE_wet_sha = LAI * Fc * lEch_wet;
        lE_wet_sun = LAI * equations.meanleaf( ...
            canopy, lEcu_wet, 'angles_and_layers', Ps);
        lE_dry_sha = LAI * Fc * lEch - lE_wet_sha;
        lE_dry_sun = LAI * equations.meanleaf( ...
            canopy, lEcu, 'angles_and_layers', Ps) - lE_wet_sun;
        lE_wet_canopy_total = lE_wet_sha + lE_wet_sun;
        Hctot        = LAI * (Fc * Hch       + equations.meanleaf(canopy, Hcu,                      'angles_and_layers', Ps));
        Actot        = LAI * (Fc * Ahh       + equations.meanleaf(canopy, Auu,                      'angles_and_layers', Ps));
        % Leaf dark respiration Rdtot (uniform branch); used by NPP downstream
        Rdtot        = LAI * (Fc * (Ahh - Ah) + equations.meanleaf(canopy, Auu - Au,                 'angles_and_layers', Ps));
        Tcave        =       (Fc * Tch       + equations.meanleaf(canopy, Tcu,                      'angles_and_layers', Ps));
        Pntot        = LAI * (Fc * Pinh      + equations.meanleaf(canopy, Pinu,                     'angles_and_layers', Ps));
        Pntot_Cab    = LAI * (Fc * Pinh_Cab  + equations.meanleaf(canopy, Pinu_Cab,                 'angles_and_layers', Ps));
        Rntot_PAR    = LAI * (Fc * Rnh_PAR   + equations.meanleaf(canopy, Rnu_PAR,                  'angles_and_layers', Ps));
        aPAR_Cab_eta = LAI * (Fc * (profiles.etah .* Rnh_PAR) + equations.meanleaf(canopy, profiles.etau .* Rnu_PAR, 'angles_and_layers', Ps));
    end
    Rnstot          = Fs * Rns;           % [W m-2] soil net radiation
    Rnstot_bare = (1 - snowFraction) * (Fs * RnsBare);
    Rnstot_snow = snowFraction * (Fs * RnsSnow);
    lEstot          = Fs * lEs;           % [W m-2] soil latent heat
    lEstot_bare = Fs * groundFlux. ...
        bare_ground_upward_latent_heat_flux_contribution_W_m2;
    lEstot_snow = Fs * groundFlux. ...
        snow_upward_latent_heat_flux_contribution_W_m2;
    qv_bare_net_kg_m2_s = Fs * groundFlux. ...
        bare_ground_upward_water_vapor_flux_contribution_kg_m2_s;
    qv_snow_net_kg_m2_s = Fs * groundFlux. ...
        snow_upward_water_vapor_flux_contribution_kg_m2_s;
    evaporation_bare_kg_m2_s = max(qv_bare_net_kg_m2_s, 0);
    bare_dew_kg_m2_s = max(-qv_bare_net_kg_m2_s, 0);
    sublimation_snow_kg_m2_s = max(qv_snow_net_kg_m2_s, 0);
    snow_deposition_kg_m2_s = max(-qv_snow_net_kg_m2_s, 0);
    Hstot_bare = Fs * groundFlux. ...
        bare_ground_upward_sensible_heat_flux_contribution_W_m2;
    Hstot_snow = Fs * groundFlux. ...
        snow_upward_sensible_heat_flux_contribution_W_m2;
    bareExternalEnergyFlux = Fs * (RnsBare ...
        - groundFlux.bare_ground_upward_sensible_heat_flux_tile_W_m2 ...
        - groundFlux.bare_ground_upward_latent_heat_flux_tile_W_m2);
    dBareExternalEnergyFlux_dT = Fs * ( ...
        tileRadiation.bare_dnet_dtemperature_W_m2_K ...
        - groundFlux.dsensible_bare_tile_dT_W_m2_K ...
        - groundFlux.dlatent_bare_tile_dT_W_m2_K);
    snowExternalEnergyFlux = Fs * (RnsSnow ...
        - groundFlux.snow_upward_sensible_heat_flux_tile_W_m2 ...
        - groundFlux.snow_upward_latent_heat_flux_tile_W_m2);
    dSnowExternalEnergyFlux_dT = Fs * ( ...
        tileRadiation.snow_dnet_dtemperature_W_m2_K ...
        - groundFlux.dsensible_snow_tile_dT_W_m2_K ...
        - groundFlux.dlatent_snow_tile_dT_W_m2_K);
    Gtot_bare_surface = Fs * G;
    Gtot_snow_surface = snowExternalEnergyFlux;
    Gtot_bare = (1 - snowFraction) * Gtot_bare_surface;
    Gtot_snow = snowFraction * Gtot_snow_surface;
    Gtot = Gtot_bare + Gtot_snow;
    TbareSave = Fs * Ts;
    Tsave = (1 - snowFraction) * TbareSave ...
        + snowFraction * snowTemperature;
    % Soil respiration: not computed here (needs main-loop precip window); Resp=0 standalone
    Resp            = 0;                  %             Soil respiration (set by COMPASS.m main loop)

    % --- Canopy + soil totals ---
    Atot            = Actot;            %                   GPP
    Rntot           = Rnctot + Rnstot;  %                   Net radiation
    lEtot           = lEctot + lEstot;  %                   Latent heat
    Htot            = Hctot  + Hstot;   %                   Sensible heat

    fluxes.Rntot    = Rntot;  % [W m-2]             total net radiation (canopy + soil)
    fluxes.lEtot    = lEtot;  % [W m-2]             total latent heat flux (canopy + soil)
    fluxes.Htot     = Htot;   % [W m-2]             total sensible heat flux (canopy + soil)
    fluxes.Atot     = Atot;   % [umol m-2 s-1]      total net CO2 uptake (canopy + soil)
    fluxes.Rnctot   = Rnctot; % [W m-2]             canopy net radiation
    fluxes.lEctot   = lEctot; % [W m-2]             canopy latent heat flux (dry-leaf transpiration + wet-leaf interception evaporation)
    fluxes.lE_wet_canopy = lE_wet_canopy_total;  % [W m-2] (v29) independent wet-leaf fwet interception evaporation component
    fluxes.lE_dry_canopy = lEctot - lE_wet_canopy_total;  % [W m-2] dry-leaf stomatal transpiration component
    % Exact diagnostic partition: wet+dry and sunlit+shaded both reconstruct
    % lEctot. These are pre-hydraulic ebal components; Tp_t is below.
    fluxes.lE_wet_sun = lE_wet_sun;  % [W m-2] sunlit wet-leaf component
    fluxes.lE_wet_sha = lE_wet_sha;  % [W m-2] shaded wet-leaf component
    fluxes.lE_dry_sun = lE_dry_sun;  % [W m-2] sunlit dry-leaf component
    fluxes.lE_dry_sha = lE_dry_sha;  % [W m-2] shaded dry-leaf component
    fluxes.Hctot    = Hctot;  % [W m-2]             canopy sensible heat flux
    fluxes.Actot    = Actot;  % [umol m-2 s-1]      canopy net CO2 uptake
    fluxes.Rnstot   = Rnstot; % [W m-2]             soil net radiation
    fluxes.Rnstot_bare = Rnstot_bare; % [W m-2] bare-area contribution
    fluxes.Rnstot_snow = Rnstot_snow; % [W m-2] snow-area contribution
    fluxes.lEstot   = lEstot; % [W m-2]             soil latent heat flux
    fluxes.lEstot_bare = lEstot_bare; % [W m-2] bare-ground contribution
    fluxes.lEstot_snow = lEstot_snow; % [W m-2] snow-covered contribution
    fluxes.evaporation_bare_kg_m2_s = evaporation_bare_kg_m2_s;
    fluxes.bare_dew_kg_m2_s = bare_dew_kg_m2_s;
    fluxes.qv_bare_net_kg_m2_s = qv_bare_net_kg_m2_s;
    fluxes.sublimation_snow_kg_m2_s = sublimation_snow_kg_m2_s;
    fluxes.snow_deposition_kg_m2_s = snow_deposition_kg_m2_s;
    fluxes.qv_snow_net_kg_m2_s = qv_snow_net_kg_m2_s;
    fluxes.Hstot    = Hstot;  % [W m-2]             soil sensible heat flux
    fluxes.Hstot_bare = Hstot_bare; % [W m-2] bare-area contribution
    fluxes.Hstot_snow = Hstot_snow; % [W m-2] snow-area contribution
    fluxes.Gtot     = Gtot;   % [W m-2]             soil heat flux
    fluxes.Gtot_bare_surface = Gtot_bare_surface; % [W m-2 per bare area]
    fluxes.Gtot_snow_surface = Gtot_snow_surface; % [W m-2 per snow area]
    fluxes.Gtot_bare = Gtot_bare; % [W m-2] bare-area grid contribution
    fluxes.Gtot_snow = Gtot_snow; % [W m-2] snow-area grid contribution
    fluxes.bare_ground_external_energy_flux_W_m2 = ...
        bareExternalEnergyFlux;
    fluxes.bare_ground_external_energy_flux_temperature_derivative_W_m2_K = ...
        dBareExternalEnergyFlux_dT;
    fluxes.bare_ground_upward_latent_heat_flux_temperature_derivative_W_m2_K = ...
        Fs * groundFlux.dlatent_bare_contribution_dT_W_m2_K;
    fluxes.snow_external_energy_flux_W_m2 = snowExternalEnergyFlux;
    fluxes.snow_external_energy_flux_temperature_derivative_W_m2_K = ...
        dSnowExternalEnergyFlux_dT;
    fluxes.scope_empirical_bare_ground_heat_flux_W_m2 = Gtot_bare_surface;
    % Canonical physical interface names. Legacy SCOPE names above remain
    % compatibility aliases only; new coupling and diagnostics use these
    % unit/direction/area-explicit fields.
    fluxes.ground_net_radiation_W_m2 = Rnstot;
    fluxes.bare_ground_net_radiation_contribution_W_m2 = Rnstot_bare;
    fluxes.snow_net_radiation_contribution_W_m2 = Rnstot_snow;
    fluxes.ground_upward_latent_heat_flux_W_m2 = lEstot;
    fluxes.bare_ground_upward_latent_heat_flux_contribution_W_m2 = lEstot_bare;
    fluxes.snow_upward_latent_heat_flux_contribution_W_m2 = lEstot_snow;
    fluxes.ground_upward_sensible_heat_flux_W_m2 = Hstot;
    fluxes.bare_ground_upward_sensible_heat_flux_contribution_W_m2 = Hstot_bare;
    fluxes.snow_upward_sensible_heat_flux_contribution_W_m2 = Hstot_snow;
    fluxes.soil_top_downward_energy_flux_W_m2 = Gtot;
    fluxes.bare_ground_downward_soil_heat_flux_contribution_W_m2 = Gtot_bare;
    fluxes.snow_downward_soil_energy_flux_contribution_W_m2 = Gtot_snow;
    fluxes.bare_ground_upward_water_vapor_flux_kg_m2_s = qv_bare_net_kg_m2_s;
    fluxes.snow_upward_water_vapor_flux_kg_m2_s = qv_snow_net_kg_m2_s;
    fluxes.Tp_t      = hydraulic.Tp_t;      % [cm s-1]  actual transpiration (hydraulic-limited)
    fluxes.Tp_demand = hydraulic.Tp_demand; % [cm s-1]  potential transpiration demand from atmosphere
    fluxes.supply_ratio = hydraulic.supply_ratio; % [-] actual / demand transpiration
    fluxes.wt_sun    = hydraulic.wt_sun;    % [-]        hydraulic stress factor sunlit leaves
    fluxes.wt_sha    = hydraulic.wt_sha;    % [-]        hydraulic stress factor shaded leaves
    fluxes.wt_sun_layer = hydraulic.wt_sun_layer; % [-]  layer-resolved sunlit hydraulic stress
    fluxes.wt_sha_layer = hydraulic.wt_sha_layer; % [-]  layer-resolved shaded hydraulic stress
    % Stomatal conductance diag: gs ≈ 1/rcw * 41.6 [mol m-2 s-1] (leaf mean, not canopy-integrated)
    gsu_arr = 1 ./ rcwu(isfinite(rcwu) & rcwu > 0);
    gsh_arr = 1 ./ rcwh(isfinite(rcwh) & rcwh > 0);
    if isempty(gsu_arr), fluxes.gs_sun = 0; else, fluxes.gs_sun = 41.6 * mean(gsu_arr(:)); end  % [mol m-2 s-1] mean sunlit leaf gs
    if isempty(gsh_arr), fluxes.gs_sha = 0; else, fluxes.gs_sha = 41.6 * mean(gsh_arr(:)); end  % [mol m-2 s-1] mean shaded leaf gs
    % [DIAG] Temporary online Medlyn/thermal validation (remove after archive)
    fluxes.online_gs_diag = assemble_online_gs_diag(biochem_diag_sun, biochem_diag_sha);
    if isnan(fluxes.Rntot)
        fluxes.lEtot = Rntot;
        fluxes.Htot = Rntot;
        fluxes.Atot = Rntot;
        fluxes.Rnctot   = Rntot; % [W m-2]             canopy net radiation
        fluxes.lEctot   = Rntot; % [W m-2]             canopy latent heat flux
        fluxes.Hctot    = Rntot;  % [W m-2]             canopy sensible heat flux
        fluxes.Actot    = Rntot;  % [umol m-2 s-1]      canopy net CO2 uptake
        fluxes.Rnstot   = Rntot; % [W m-2]             soil net radiation
        fluxes.lEstot   = Rntot; % [W m-2]             soil latent heat flux
        fluxes.Hstot    = Rntot;  % [W m-2]             soil sensible heat flux
        fluxes.Gtot     = Rntot;   % [W m-2]             soil heat flux
        fluxes.lE_wet_canopy = Rntot;  % v29 fwet interception evaporation component
        fluxes.lE_dry_canopy = Rntot;  % v29 dry-leaf transpiration component
    end
    fluxes.Resp     = Resp;   % [umol m-2 s-1]      soil respiration
    fluxes.aPAR     = Pntot;  % [mol m-2 s-1]       absorbed PAR (diagnostic writer converts to umol)
    fluxes.aPAR_Cab = Pntot_Cab; % [mol m-2 s-1]    chlorophyll-absorbed PAR
    fluxes.aPAR_Wm2 = Rntot_PAR; % [W m-2]      absorbed PAR
    fluxes.aPAR_Cab_eta = aPAR_Cab_eta;
    fluxes.GPP    = Actot * 12 / 1000000000;  % [kg C m-2 s-1] canopy GPP (Actot = Ag; biochem Ag = A + Rd)
    fluxes.Rd_leaf = Rdtot * 12 / 1000000000; % [kg C m-2 s-1] canopy leaf dark respiration (for NPP downstream)
    % NEE is not computed in ebal; only GPP / Rd_leaf / Resp placeholders here

    thermal.Ta    = Ta;       % [oC]                air temperature (as in input)
    thermal.rhoa  = rhoa;     % [kg m-3]            forcing-state air density for H/LE
    thermal.Ts    = Ts;       % [oC] legacy alias: bare-soil skin [shade; sun]
    thermal.Tbare = Ts;       % [oC] bare-soil skin [shade; sun]
    thermal.Tsnow = snowTemperature; % [oC] prognostic snow temperature
    thermal.Trad = Fs * effectiveGroundT; % [oC] equivalent RT ground temperature
    thermal.bare_ground_skin_temperature_C = Ts;
    thermal.bare_ground_effective_temperature_C = TbareSave;
    if snowFraction > 0
        thermal.snow_surface_temperature_C = snowTemperature;
    else
        thermal.snow_surface_temperature_C = NaN;
    end
    thermal.ground_radiative_temperature_C = Fs * effectiveGroundT;
    thermal.Tcave = Tcave;    % [oC]                weighted average canopy temperature
    thermal.Tsave = Tsave;    % [oC]                weighted average soil temperature
    thermal.raa   = raa;      % [s m-1]             original Wallace-Verhoef above-canopy resistance
    thermal.rah   = rah;      % [s m-1]             CLM4.5-style heat aerodynamic resistance
    thermal.raw   = raw;      % [s m-1]             CLM4.5-style water-vapor aerodynamic resistance
    thermal.z0h   = resist_out.z0h;
    thermal.z0w   = resist_out.z0w;
    thermal.zo    = canopy.zo;
    thermal.d     = canopy.d;
    thermal.zo_static = canopy.zo_static;
    thermal.d_static = canopy.d_static;
    thermal.SAI   = canopy.SAI;
    thermal.L_plus_S = canopy.L_plus_S;
    if isfield(canopy, 'Kb_direct')
        thermal.Kb_direct = canopy.Kb_direct; % [-] direct-beam extinction coefficient used as Tao
    end
    if isfield(canopy, 'Kb_source')
        thermal.Kb_source = canopy.Kb_source;
    end
    thermal.Tao   = resist_out.Tao; % [-] extinction coefficient used in resistances_clm45
    thermal.Va    = resist_out.Va;
    thermal.ustar_clm = resist_out.ustar_clm;
    thermal.Kh    = resist_out.Kh;  % [m2 s-1] canopy-top/roughness-layer heat diffusivity
    thermal.uh    = resist_out.uh;  % [m s-1]  wind speed at canopy top
    thermal.n_canopy = resist_out.n_canopy; % [-] W&V canopy wind/K extinction coefficient
    if isfield(profiles, 'Ta_layer')
        thermal.Ta_layer  = profiles.Ta_layer;  % [degC] diagnostic layer canopy-air temperature
        thermal.ea_layer  = profiles.ea_layer;  % [hPa]  diagnostic layer canopy-air vapour pressure
        thermal.CO2_layer = profiles.CO2_layer; % [umol mol-1] diagnostic layer canopy-air CO2
        thermal.Kh_layer  = profiles.Kh_layer;  % [m2 s-1] diagnostic layer turbulent diffusivity
    end
    thermal.rawc  = rawc;     % [s m-1]             aerodynamic resistance below canopy for canopy
    thermal.raws  = raws;     % [s m-1]             aerodynamic resistance below canopy for soil
    thermal.rbc   = resist_out.rbc;     % [s m-1]   canopy boundary-layer resistance
    thermal.rbc_wv = resist_out.rbc_wv; % [s m-1]   previous Wallace-Verhoef empirical rbc
    thermal.rss   = soil_surface.rss;            % [s m-1]   dynamic soil surface resistance
    thermal.rbs   = soil_surface.rbs;            % [s m-1]   dynamic soil boundary resistance
    thermal.rbs_init = soil_surface.rbs_init;    % [s m-1]   input rbs at ebal entry
    thermal.rbs_static = soil_surface.rbs_static; % [s m-1]  input rbs before dynamic update
    thermal.Sat_surf = soil_surface.Sat_surf;    % [-]       surface soil saturation
    thermal.epsilon_s = soil_surface.epsilon_s;  % [-]       dynamic soil thermal emissivity
    thermal.rs_thermal = soil_surface.rs_thermal; % [-]      dynamic soil thermal reflectance
    thermal.rs_thermal_static = soil_surface.rs_thermal_static; % [-] original input thermal reflectance
    thermal.snow_cover_fraction = groundFlux.cover_fraction; % [-]
    thermal.ustar = ustar;    % [m s-1]             friction velocity
    thermal.Tcu   = Tcu;
    thermal.Tch   = Tch;

    fluxes.Au     = Au;
    fluxes.Ah     = Ah;

    % --- Root water uptake RWU: split actual T onto soil nodes (may include HR) ---
    % Positive = soil->root, negative = root->soil; enters Richards as sink
    [RWU, hydraulic] = vegetation.root.assembleRWU( ...
        hydraulic, soil_state, plant_hydro, domain, DeltZ);
    frac = zeros(size(PSIs));

    profiles.Knu     = Knu;
    profiles.Knh     = Knh;

    %% ========================================================================
    %%  [G] Pack output 'out' (applyEbalOutput writes it to scope_step)
    %% ========================================================================
    out = struct( ...
        'iteration', iteration, ...
        'fluxes', fluxes, ...
        'radiation', rad, ...
        'gap', gap, ...
        'thermal', thermal, ...
        'profiles', profiles, ...
        'soil_surface', soil_surface, ...
        'root_water_uptake', RWU, ...
        'frac', frac, ...
        'hydraulic', hydraulic);
end

function diag = assemble_online_gs_diag(bio_sun, bio_sha)
%ASSEMBLE_ONLINE_GS_DIAG  Scalar online Medlyn and thermal diagnostics.
%
% Inputs:
%   bio_sun  [in]  struct from vegetation.biochemical for sunlit leaves
%   bio_sha  [in]  struct from vegetation.biochemical for shaded leaves
%
% Output:
%   diag     [out] scalar fields written only when STEMMUS_ONLINE_GS_DIAG=1
%
% [DIAG 2026-06-19] Temporary validation support. Delete after the annual
% Medlyn/SOX/thermal-acclimation validation is archived.
    diag = struct();
    diag.medlyn_gs_model_sun    = mean_field(bio_sun, 'gs');
    diag.medlyn_gs_model_sha    = mean_field(bio_sha, 'gs');
    diag.medlyn_gs_from_rcw_sun = mean_field(bio_sun, 'gs_from_rcw');
    diag.medlyn_gs_from_rcw_sha = mean_field(bio_sha, 'gs_from_rcw');
    diag.medlyn_gs_residual_sun = mean_field(bio_sun, 'gs_medlyn_residual');
    diag.medlyn_gs_residual_sha = mean_field(bio_sha, 'gs_medlyn_residual');
    diag.medlyn_a_sun           = mean_field(bio_sun, 'A');
    diag.medlyn_a_sha           = mean_field(bio_sha, 'A');
    diag.medlyn_ag_sun          = mean_field(bio_sun, 'Ag');
    diag.medlyn_ag_sha          = mean_field(bio_sha, 'Ag');
    diag.medlyn_cs_sun          = mean_field(bio_sun, 'Cs');
    diag.medlyn_cs_sha          = mean_field(bio_sha, 'Cs');
    diag.medlyn_ci_sun          = mean_field(bio_sun, 'Ci');
    diag.medlyn_ci_sha          = mean_field(bio_sha, 'Ci');
    diag.medlyn_vpd_sun         = mean_field(bio_sun, 'VPD');
    diag.medlyn_vpd_sha         = mean_field(bio_sha, 'VPD');
    diag.medlyn_tleaf_sun       = mean_field(bio_sun, 'Tleaf');
    diag.medlyn_tleaf_sha       = mean_field(bio_sha, 'Tleaf');
    diag.medlyn_q_sun           = mean_field(bio_sun, 'Q');
    diag.medlyn_q_sha           = mean_field(bio_sha, 'Q');
    diag.medlyn_sfactor_sun     = mean_field(bio_sun, 'sfactor');
    diag.medlyn_sfactor_sha     = mean_field(bio_sha, 'sfactor');
    diag.medlyn_g0              = mean_pair(bio_sun, bio_sha, 'g0_medlyn');
    diag.medlyn_g1              = mean_pair(bio_sun, bio_sha, 'g1_medlyn');
    % Optional heat-induced gs (0 when heat_gs_enable=false)
    diag.g_heat_sun             = mean_field(bio_sun, 'g_heat');
    diag.g_heat_sha             = mean_field(bio_sha, 'g_heat');
    diag.gs0_medlyn_sun         = mean_field(bio_sun, 'gs0_medlyn');
    diag.gs0_medlyn_sha         = mean_field(bio_sha, 'gs0_medlyn');
    diag.thermal_tgrowth        = mean_pair(bio_sun, bio_sha, 'Tgrowth');
    diag.thermal_tgrowth_clamped = mean_pair(bio_sun, bio_sha, 'Tgrowth_clamped');
    diag.thermal_delta_s_vcmax  = mean_pair(bio_sun, bio_sha, 'DeltaS_Vcmax');
    diag.thermal_delta_s_jmax   = mean_pair(bio_sun, bio_sha, 'DeltaS_Jmax');
    diag.thermal_topt_vcmax     = mean_pair(bio_sun, bio_sha, 'Topt_Vcmax');
    diag.thermal_topt_jmax      = mean_pair(bio_sun, bio_sha, 'Topt_Jmax');
    diag.thermal_vcmax25_sun    = mean_field(bio_sun, 'Vcmax25_base');
    diag.thermal_vcmax25_sha    = mean_field(bio_sha, 'Vcmax25_base');
    diag.thermal_jmax25_sun     = mean_field(bio_sun, 'Jmax25');
    diag.thermal_jmax25_sha     = mean_field(bio_sha, 'Jmax25');
    diag.thermal_rd25_sun       = mean_field(bio_sun, 'Rd25');
    diag.thermal_rd25_sha       = mean_field(bio_sha, 'Rd25');
    diag.thermal_vcmax_sun      = mean_field(bio_sun, 'Vcmax');
    diag.thermal_vcmax_sha      = mean_field(bio_sha, 'Vcmax');
    diag.thermal_jmax_sun       = mean_field(bio_sun, 'Jmax');
    diag.thermal_jmax_sha       = mean_field(bio_sha, 'Jmax');
    diag.thermal_rd_sun         = mean_field(bio_sun, 'Rd');
    diag.thermal_rd_sha         = mean_field(bio_sha, 'Rd');
end

function value = mean_pair(s1, s2, field_name)
    value = finite_mean([mean_field(s1, field_name), mean_field(s2, field_name)]);
end

function value = mean_field(s, field_name)
    if ~isstruct(s) || ~isfield(s, field_name)
        value = NaN;
        return
    end
    x = s.(field_name);
    if isempty(x) || ~(isnumeric(x) || islogical(x)) || ~isreal(x)
        value = NaN;
        return
    end
    value = finite_mean(double(x(:)));
end

function value = finite_mean(x)
    x = x(isfinite(x));
    if isempty(x)
        value = NaN;
    else
        value = mean(x);
    end
end

function tol = get_ebal_tolerance(control, fallback)
    tol = fallback;
    if isstruct(control) && isfield(control, 'tol_eb_W_m2') && ...
            isfinite(control.tol_eb_W_m2) && control.tol_eb_W_m2 > 0
        tol = control.tol_eb_W_m2;
    end
end
function biochem_in = apply_phase3_biochem_fields(biochem_in, leafbio, canopy, domain, xyt, k)
    biochem_in.Nefc = get_struct_field(leafbio, 'photo_Nefc', 0);
    biochem_in.FN = get_struct_field(leafbio, 'photo_FN', []);
    biochem_in.Tgrowth = get_struct_field(canopy, 'Tgrowth', []);
    % [Atkin2015 leaf Rd 2026-06-14] growth-T acclimated leaf Rd25 toggle + params (biochemical.m)
    biochem_in.leafresp_atkin      = get_struct_field(canopy, 'leafresp_atkin', 0);
    biochem_in.lmr_intercept_atkin = get_struct_field(canopy, 'lmr_intercept_atkin', []);
    biochem_in.lnc                 = get_struct_field(canopy, 'lnc', []);
    biochem_in.latitude = domain.latitude_deg;

    % Solar declination via surface.radiation.calculatesunpos_ultimate_azi_atan22
    idx = max(1, min(k, numel(xyt.t)));
    t_doy = xyt.t(idx);
    if isfield(xyt, 'year') && ~isempty(xyt.year)
        yr_idx = min(idx, numel(xyt.year));
        yr = xyt.year(yr_idx);
        if yr < 1900, yr = 2021; end
    else
        yr = 2021;
    end
    doy_int = floor(t_doy) + 1;
    hour_frac = (t_doy - floor(t_doy)) * 24;
    dt0 = datetime(yr, 1, 1) + days(doy_int - 1);
    [~, ~, ~, ~, ~, ~, ~, ~, ~, delta_deg] = ...
        surface.radiation.calculatesunpos_ultimate_azi_atan22( ...
            yr, month(dt0), day(dt0), hour_frac, ...
            domain.latitude_deg, domain.longitude_deg, 1);
    biochem_in.decl = deg2rad(delta_deg);

    biochem_in.gs_model_type = get_struct_field(leafbio, 'gs_model_type', 1);
    biochem_in.g0_medlyn = get_struct_field(leafbio, 'medlyn_g0', 1e-4);
    biochem_in.g1_medlyn = get_struct_field(leafbio, 'medlyn_g1', 4.86);
    biochem_in.g_m = get_struct_field(leafbio, 'g_m', Inf);  % [mol m-2 s-1]; Inf disables mesophyll limitation
    % Optional heat-induced gs after Medlyn (default off). leafbio fields from prepareScopeStep.
    biochem_in.heat_gs_enable   = get_struct_field(leafbio, 'heat_gs_enable', 0);
    biochem_in.heat_gs_T_on_C   = get_struct_field(leafbio, 'heat_gs_T_on_C', 32.0);
    biochem_in.heat_gs_T_full_C = get_struct_field(leafbio, 'heat_gs_T_full_C', 38.0);
    biochem_in.heat_gs_max      = get_struct_field(leafbio, 'heat_gs_max', 0.08);
    biochem_in.heat_gs_wt_off   = get_struct_field(leafbio, 'heat_gs_wt_off', 0.4);
    biochem_in.heat_gs_wt_on    = get_struct_field(leafbio, 'heat_gs_wt_on', 0.8);
    biochem_in.heat_gs_D_lo_kPa = get_struct_field(leafbio, 'heat_gs_D_lo_kPa', 1.5);
    biochem_in.heat_gs_D_hi_kPa = get_struct_field(leafbio, 'heat_gs_D_hi_kPa', 3.5);
    biochem_in.heat_gs_A_lo     = get_struct_field(leafbio, 'heat_gs_A_lo', 1.0);
    biochem_in.heat_gs_A_hi     = get_struct_field(leafbio, 'heat_gs_A_hi', 5.0);
    biochem_in.heat_gs_cap      = get_struct_field(leafbio, 'heat_gs_cap', 0.6);

    % [v5 2026-05-27] LUNA strict path (CLM5 PhotosynthesisMod.F90:1515-1517).
    % Without these 3 fields, biochemical.m path A is never true (falls back to Kattge R_JV).
    % With use_luna=1 and non-empty Jmax25_LUNA: path A (LUNA independent Jmax).
    biochem_in.use_luna     = get_struct_field(leafbio, 'use_luna', 0);
    biochem_in.Vcmax25_LUNA = get_struct_field(leafbio, 'Vcmax25_LUNA', []);
    biochem_in.Jmax25_LUNA  = get_struct_field(leafbio, 'Jmax25_LUNA', []);
end


function value = get_struct_field(s, field_name, default_value)
    if isstruct(s) && isfield(s, field_name) && ~isempty(s.(field_name))
        value = s.(field_name);
    else
        value = default_value;
    end
end
