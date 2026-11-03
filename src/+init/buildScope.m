function scope = buildScope(raw, built)
% Editor: Wenjie Liu; Date: 2026-10-17.
% Change: 2026-10-17 Wenjie Liu -- materialize the sampled surface-resistance multiplier.
%BUILDSCOPE  Build the SCOPE radiative-transfer parameters (S1 + S2).
%
%   scope = init.buildScope(raw, built) builds the SCOPE parameter container.
%   S1 (always): the four config sections [scope]/[leaf_optical]/[soil_optical]/
%   [geometry] -> scope.{options,leaf_optical,soil_optical,geometry}, mapped to
%   legacy consumer names with bool->0/1. S2 (when built supplied): the hard-coded
%   spectral grid (io.define_bands), the RT discretization (nlayers/litab/lazitab
%   + static layer_weight from the veg LAI profile) and the external optical
%   NetCDF resource paths ([model_library]) -> scope.{spectral,canopy_rt,
%   optical_library}.
%   Pure config/resource assembly -- no physics.
%
%   Inputs:
%     raw    struct  io.loadSite output (uses .config sections + .siteDir)
%     built  struct  needs parameters.veg.canopy (LAI_profile_* + LAI_constant)
%                    for the static layer_weight. Optional: omit for S1 only.
%
%   Returns scope with sub-structs:
%     scope.options/leaf_optical/soil_optical/geometry   (S1)
%     scope.spectral         wlS/wlP/wlT/wlF/wlPAR/... (io.define_bands)     (S2)
%     scope.canopy_rt        nlayers/nlincl/nlazi/litab/lazitab/x/xl +
%                            layer_weight/shape/z (static; LAI/LAI_layer stay veg) (S2)
%     scope.optical_library  resolved [model_library] NetCDF file paths      (S2)
%
%   BASELINES ONLY: Cab/Cca (LUNA may override per step), rss/rbs/rs_thermal (ebal
%   recomputes from soil moisture), tts (per-hour) and LAI_layer (= layer_weight x
%   LAI, runtime) are NOT frozen here. Reference: select_input / io.define_bands /
%   vegetation.calc_lai_profile / run_stemmus_initialize canopy section.
%
%   Editor: Jizixing Zhu; Date: 2026-06-26.
%   Change: 2026-07-16 Jizixing Zhu -- resolve model_library.directional_angles_file
%     (.nc under data/model_library_nc/) for BRDF when calc_directional.
%   Change: 2026-07-16 Jizixing Zhu -- carry the site canopy-profile switch
%     into the SCOPE options consumed by the energy-balance kernels.
%   Change: 2026-07-03 Jizixing Zhu — options field Fluorescence_model (capital F, SCOPE convention read by ebal/biochemical) was mistakenly lowercased; nothing read the lowercase key.

    rawConfig = raw.config;
    scope = struct();

    % --- options: SCOPE calculation switches (booleans -> 0/1; enums kept) ---
    s = rawConfig.scope;
    o = struct();
    o.simulation           = s.simulation;
    o.calc_ebal            = toFlag(s.calc_ebal);
    o.calc_fluor           = toFlag(s.calc_fluor);
    o.calc_PSI             = toFlag(s.calc_PSI);
    o.calc_directional     = toFlag(s.calc_directional);
    o.calc_planck          = toFlag(s.calc_planck);
    o.calc_xanthophyll_abs = toFlag(s.calc_xanthophyll_abs);
    o.calc_rss_rbs         = toFlag(s.calc_rss_rbs);
    o.Fluorescence_model   = s.fluorescence_model;
    o.apply_T_corr         = toFlag(s.apply_T_corr);
    o.soilspectrum         = s.soilspectrum;
    o.rt_thermal           = s.rt_thermal;
    o.Cca_function_of_Cab  = toFlag(s.cca_function_of_cab);
    o.soil_heat_method     = s.soil_heat_method;
    o.canopy_layer_profiles = toFlag( ...
        rawConfig.vegetation.canopy_layer_profiles);
    scope.options = o;

    % --- leaf_optical: PROSPECT/fluspect baselines (legacy leafbio names) ---
    l = rawConfig.leaf_optical;
    lo = struct();
    lo.Cab         = l.Cab;
    lo.Cca         = l.Cca;
    lo.Cdm         = l.Cdm;
    lo.Cw          = l.Cw;
    lo.Cs          = l.Cs;
    lo.Cant        = l.Cant;
    lo.N           = l.N_mesophyll;       % friendly N_mesophyll -> leafbio.N
    lo.rho_thermal = l.leaf_rho_thermal;  % -> leafbio.rho_thermal
    lo.tau_thermal = l.leaf_tau_thermal;  % -> leafbio.tau_thermal
    lo.fqe         = l.fqe;
    lo.V2Z         = l.V2Z;
    lo.beta        = l.fluor_beta;        % fluor NPQ tuning -> leafbio.*
    lo.kNPQs       = l.fluor_kNPQs;
    lo.qLs         = l.fluor_qLs;
    lo.stressfactor = l.fluor_stressfactor;
    scope.leaf_optical = lo;

    % --- soil_optical: reflectance + evaporation resistance baselines ---
    so = rawConfig.soil_optical;
    sm = struct();
    sm.rss_scale=1.0;
    if isfield(rawConfig.scope,'soil_optical') && isfield(rawConfig.scope.soil_optical,'rss_scale')
        sm.rss_scale=rawConfig.scope.soil_optical.rss_scale;
    end
    validateattributes(sm.rss_scale,{'numeric'},{'real','scalar','positive','finite'});
    sm.spectrum      = so.soil_spectrum_index;   % column in soil_spectrum_file
    sm.rss           = so.rss;                    % baseline (ebal recomputes per KT)
    sm.rbs           = so.rbs;                    % baseline
    sm.rs_thermal    = so.rs_thermal;             % baseline (recomputed from moisture)
    sm.Rho_Thermal_Soil_Dry = rawConfig.environment.rho_thermal_soil_dry;
    sm.Rho_Thermal_Soil_Wet = rawConfig.environment.rho_thermal_soil_wet;
    sm.BSMBrightness = so.BSM_brightness;
    sm.BSMlat        = so.BSM_lat;                % BSM spectral-shape param (NOT site lat)
    sm.BSMlon        = so.BSM_lon;
    scope.soil_optical = sm;

    % --- geometry: viewing angles (solar zenith tts auto-computed elsewhere) ---
    g = rawConfig.geometry;
    scope.geometry = struct( ...
        'tto', g.observer_zenith, ...   % viewing zenith
        'psi', g.relative_azimuth);     % sun-viewing relative azimuth

    % --- S2: spectral grid + RT discretization + external optical resources ---
    if nargin >= 2 && ~isempty(built)
        scope.spectral        = io.define_bands();                          % hard-coded wlS/wlP/...
        scope.canopy_rt       = buildCanopyRT(built.parameters.veg.canopy); % discretization + layer_weight
        scope.optical_library = resolveOpticalLibrary(raw);                 % [model_library] paths
    end
end

% ---------------------------------------------------------------------------
function rt = buildCanopyRT(canopy)
% RT discretization (run_stemmus_initialize canopy section) + static layer_weight.
% layer_weight/shape from vegetation.calc_lai_profile at a REFERENCE LAI (LAI_constant);
% it is LAI-independent, so the runtime LAI_layer = layer_weight * canopy.LAI.
% Editor: Jizixing Zhu; Date: 2026-06-27.
    rt = struct();
    rt.nlayers = 30;
    rt.nlincl  = 13;
    rt.nlazi   = 36;
    nl = rt.nlayers;
    rt.x  = (-1/nl : -1/nl : -1)';
    rt.xl = [0; rt.x];
    rt.litab   = [5:10:75 81:2:89]';
    rt.lazitab = (5:10:355)';

    [~, rt.layer_weight, rt.LAI_profile_shape, rt.LAI_profile_z] = vegetation.calc_lai_profile( ...
        canopy.LAI_constant, nl, canopy.LAI_profile_mode, canopy.LAI_profile_peak, ...
        canopy.LAI_profile_sigma, canopy.LAI_profile_skew, canopy.LAI_profile_k);
    rt.LAI_profile_mode         = canopy.LAI_profile_mode;
    rt.LAI_profile_reference_LAI = canopy.LAI_constant;     % reference only; runtime LAI lives in veg
end

% ---------------------------------------------------------------------------
function ol = resolveOpticalLibrary(raw)
% Resolve shared SCOPE reference-data file paths from [model_library].
% data/model_library_nc lives alongside data/sites (two levels up from the site dir).
% Editor: Jizixing Zhu; Date: 2026-06-27.
    ml    = raw.config.model_library;
    mldir = fullfile(fileparts(fileparts(raw.siteDir)), 'model_library_nc');
    format = 'netcdf';
    if isfield(ml, 'format')
        format = lower(strtrim(char(ml.format)));
    end
    if ~strcmp(format, 'netcdf')
        error('buildScope:unsupportedLibraryFormat', ...
            'Only model_library.format = "netcdf" is supported; got "%s".', format);
    end

    ol = struct();
    ol.format = format;
    ol.model_library_dir = mldir;
    ol.leaf_optics_file   = resolveLibraryFile(mldir, ml.leaf_optics_file, 'leaf_optics_file');
    ol.atmosphere_file    = resolveLibraryFile(mldir, ml.atmosphere_file, 'atmosphere_file');
    ol.soil_spectrum_file = resolveLibraryFile(mldir, ml.soil_spectrum_file, 'soil_spectrum_file');
    % OPTIONAL: measured leaf inclination distribution file (file-based LIDF). Empty when the
    % site uses the analytic LIDFa/LIDFb parameters instead (legacy COMPASS.m:226-231 branch).
    ol.leaf_angles_file   = resolveOptionalLibraryFile(mldir, ml, 'leaf_angles_file');
    % BRDF observer angles (NetCDF only). Default relative path when site omits the key.
    % Loaded only if scope.options.calc_directional (init.buildScopeStatic).
    if isfield(ml, 'directional_angles_file') ...
            && ~isempty(ml.directional_angles_file) ...
            && ~isempty(strtrim(char(ml.directional_angles_file)))
        ol.directional_angles_file = resolveLibraryFile( ...
            mldir, ml.directional_angles_file, 'directional_angles_file');
    else
        ol.directional_angles_file = resolveLibraryFile( ...
            mldir, 'directional/brdf_angles2.nc', 'directional_angles_file');
    end
end

% ---------------------------------------------------------------------------
function filePath = resolveOptionalLibraryFile(rootDir, ml, key)
% Resolve an OPTIONAL model_library file path: '' when the key is absent/blank, else a
% required .nc path (same rules as resolveLibraryFile). Used for leaf_angles_file, which
% only file-based-LIDF sites set.
% Editor: Jizixing Zhu; Date: 2026-07-03.
    if ~isfield(ml, key) || isempty(ml.(key)) || isempty(strtrim(char(ml.(key))))
        filePath = '';
        return;
    end
    filePath = resolveLibraryFile(rootDir, ml.(key), key);
end

% ---------------------------------------------------------------------------
function filePath = resolveLibraryFile(rootDir, relPath, key)
% Resolve one data/model_library_nc relative path and require a NetCDF file.
% Editor: Jizixing Zhu; Date: 2026-06-27.
    relPath = char(relPath);
    if startsWith(relPath, filesep) || ~isempty(regexp(relPath, '^[A-Za-z]:', 'once'))
        error('buildScope:absoluteLibraryPath', ...
            'model_library.%s must be relative to data/model_library_nc, got "%s".', key, relPath);
    end

    parts = regexp(relPath, '[\\/]', 'split');
    filePath = fullfile(rootDir, parts{:});
    if ~endsWith(lower(filePath), '.nc')
        error('buildScope:nonNetcdfLibraryFile', ...
            'model_library.%s must point to a .nc file, got "%s".', key, relPath);
    end
    if ~isfile(filePath)
        error('buildScope:missingLibraryFile', ...
            'model_library.%s does not exist: %s', key, filePath);
    end
end

% ---------------------------------------------------------------------------
function f = toFlag(x)
% Logical / numeric switch -> 0/1 double.
% Editor: Jizixing Zhu; Date: 2026-06-27.
    f = double(logical(x));
end
