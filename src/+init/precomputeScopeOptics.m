function optics = precomputeScopeOptics(scope, veg, library, spectral)
%PRECOMPUTESCOPEOPTICS  Precompute the once-per-run static SCOPE optical quantities.
%
%   optics = init.precomputeScopeOptics(scope, veg, library, spectral) computes
%   the leaf/canopy/soil optical resources that legacy COMPASS.m built ONCE before
%   the time loop (COMPASS.m:218-269) and that RTMo/RTMt then consume each step.
%   These depend only on constant config (leaf optical params, LIDF, soil spectrum)
%   so they are computed here and cached in built.p.scope_static (via init.buildScopeStatic),
%   not rebuilt per time step.
%
%   Inputs:
%     scope    [in] struct  built.parameters.scope: .leaf_optical (Cab/Cw/.../
%                          rho_thermal/tau_thermal), .options (calc_PSI/rt_thermal/
%                          soilspectrum).
%     veg      [in] struct  built.parameters.veg: .canopy.LIDFa/LIDFb.
%     library  [in] struct  init.loadScopeStaticLibrary output: .optipar,
%                          .soil.rs (soil reflectance column, soilspectrum==0).
%     spectral [in] struct  scope.spectral incl. .IwlP/.IwlT band index sets.
%
%   Returns optics [struct]:
%     .lidf_static  [13x1]  leaf inclination distribution over litab; from a measured LIDF
%                          file (library.leaf_angles) if configured, else equations.leafangles
%                          on LIDFa/LIDFb; static (COMPASS.m:226-231)
%     .emis_static  [-]     leaf thermal emissivity 1-rho_thermal-tau_thermal
%                          (COMPASS.m:232)
%     .leafopt_s    [struct] static fluspect leaf optics: .refl/.tran [nwl x1],
%                          .kChlrel [nwlP x1], .reflZ/.tranZ [nwl x1] (COMPASS.m:246-268)
%     .rs_static    [nwl x1] soil reflectance; optical band from the soil-spectrum
%                          file, thermal band a per-step-overwritten placeholder
%                          (present only when soilspectrum==0; COMPASS.m:250-269)
%
%   Method / reference:
%     Faithful port of the "[7] Static SCOPE spectral pre-computation" block of
%     legacy src/COMPASS.m:218-269. No state is advanced; this is immutable
%     reference data for the RT kernels.
%
%   Editor: Jizixing Zhu; Date: 2026-07-03.
%   Change: 2026-07-03 Jizixing Zhu — lidf_static now supports both a measured LIDF file (library.leaf_angles) and the analytic LIDFa/LIDFb, mirroring the legacy FileName.LIDF_file if/else.
%   Change: 2026-07-05 Jizixing Zhu — doc: output cached in built.p.scope_static (was runtime.scope_static); increment 4.

    optics = struct();

    % --- lidf_static: leaf inclination distribution over litab (COMPASS.m:226-231) ---
    % Two ways, mirroring the legacy if/else: a MEASURED file (library.leaf_angles, present
    % when the site sets model_library.leaf_angles_file) OR the analytic Verhoef two-parameter
    % form on LIDFa/LIDFb. LIDFa/LIDFb source + defaults match buildStepCanopy (prepareScopeStep.m).
    if isfield(library, 'leaf_angles') && ~isempty(library.leaf_angles)
        optics.lidf_static = library.leaf_angles(:);           % file method (measured LIDF, .nc)
    else
        lidfa = getFieldOr(veg.canopy, 'LIDFa', -0.35);
        lidfb = getFieldOr(veg.canopy, 'LIDFb', -0.15);
        optics.lidf_static = equations.leafangles(lidfa, lidfb);   % parameter method (Verhoef)
    end

    % --- emis_static: leaf thermal-band emissivity (COMPASS.m:232) ---
    optics.emis_static = 1 - scope.leaf_optical.rho_thermal - scope.leaf_optical.tau_thermal;

    % --- leafopt_s: static fluspect leaf optics on the full band grid (COMPASS.m:221-268) ---
    % Fluspect runs on the raw config leaf_optical: Cca is the config value, NOT the
    % Cca_function_of_Cab per-step override (legacy applies that only in select_input /
    % the cab_dyn feedback, i.e. AFTER this static point, so the static optics use raw Cca).
    % V2Z=0 gives the base optics; V2Z=1 (leafoptZ) supplies the xanthophyll-transition
    % bins 1:300 of reflZ/tranZ (consumed by RTMz, not RTMo).
    options = scope.options;
    if getFieldOr(options, 'calc_PSI', 0)
        fversion = @surface.radiation.fluspect_B_CX;
    else
        fversion = @surface.radiation.fluspect_B_CX_PSI_PSII_combined;
    end
    lb = scope.leaf_optical;
    lb.V2Z = 0; leafopt_s  = fversion(spectral, lb, library.optipar);
    lb.V2Z = 1; leafoptZ_s = fversion(spectral, lb, library.optipar);

    IwlP = spectral.IwlP;  IwlT = spectral.IwlT;
    nwlP = numel(spectral.wlP);  nwlT = numel(spectral.wlT);  nwl = numel(spectral.wlS);
    rho = zeros(nwl, 1);  tau = zeros(nwl, 1);
    rho(IwlP) = leafopt_s.refl;  tau(IwlP) = leafopt_s.tran;   % optical band from fluspect (COMPASS.m:248)
    rt_thermal = getFieldOr(options, 'rt_thermal', 0);
    switch rt_thermal
        case 0   % constant leaf thermal reflectance/transmittance (COMPASS.m:256-257)
            rho(IwlT) = ones(nwlT, 1) * scope.leaf_optical.rho_thermal;
            tau(IwlT) = ones(nwlT, 1) * scope.leaf_optical.tau_thermal;
        case 1   % extend the last optical band into the thermal range (COMPASS.m:260-261)
            rho(IwlT) = ones(nwlT, 1) * rho(nwlP);
            tau(IwlT) = ones(nwlT, 1) * tau(nwlP);
    end
    leafopt_s.refl = rho;  leafopt_s.tran = tau;               % full-grid (kChlrel stays optical-length; COMPASS.m:264)
    reflZ = rho;  tranZ = tau;                                 % Z-version bins 1:300 from V2Z=1 (COMPASS.m:265-268)
    reflZ(1:300) = leafoptZ_s.refl(1:300);
    tranZ(1:300) = leafoptZ_s.tran(1:300);
    leafopt_s.reflZ = reflZ;  leafopt_s.tranZ = tranZ;
    optics.leafopt_s = leafopt_s;

    % --- rs_static: soil reflectance on the full band grid (COMPASS.m:250-269) ---
    % File-based soil spectrum only (soilspectrum==0; BSM==1 is built per step). The
    % THERMAL band is a static placeholder: when rt_thermal==0 buildStepSoil overwrites
    % it per step with the moisture-corrected rs_thermal_dyn; the optical band persists.
    if getFieldOr(options, 'soilspectrum', 0) == 0
        rs = zeros(nwl, 1);
        rs(IwlP) = library.soil.rs;                            % optical band = rsfile(:,spectrum+1) (COMPASS.m:251)
        switch rt_thermal
            case 0, rs(IwlT) = ones(nwlT, 1) * scope.soil_optical.rs_thermal;   % (COMPASS.m:258)
            case 1, rs(IwlT) = ones(nwlT, 1) * rs(nwlP);                        % (COMPASS.m:262)
        end
        optics.rs_static = rs;
    end
end

% -------------------------------------------------------------------------
function v = getFieldOr(s, name, default)
% Return s.(name) when present and non-empty, else the default.
% Editor: Jizixing Zhu; Date: 2026-07-03.
    if isstruct(s) && isfield(s, name) && ~isempty(s.(name))
        v = s.(name);
    else
        v = default;
    end
end
