function s = buildScopeStatic(built)
%BUILDSCOPESTATIC  Precompute the immutable SCOPE static optical resources into the init cache.
%
%   s = init.buildScopeStatic(built) loads the SCOPE static optical library and precomputes the
%   once-per-run leaf/canopy/soil optics (lidf_static / emis_static / leafopt_s / rs_static) that
%   RTMo/RTMt then read every step. Static: derived from built.parameters.scope + built.parameters.veg,
%   never mutated per step. Stored in built.p.scope_static (the recomputable auxiliary cache; runtime
%   no longer owns it, 2026-07-05 increment 4).
%
%   Inputs:
%     built  [in] struct  needs built.parameters.scope (.options/.optical_library/.spectral/.canopy_rt)
%                         + built.parameters.veg (.canopy.LIDFa/LIDFb [-])
%
%   Returns:
%     s      [struct]  scope_static cache: .status/.options/.optical_library/.spectral/.canopy_rt/.library
%                      + precompute outputs .lidf_static [13x1] / .emis_static [-] / .leafopt_s [struct]
%                      / .rs_static [nwl x1] (present only when soilspectrum==0)
%
%   Method / reference:
%     Verbatim move of the former prepareRuntime.prepareScopeStatic local (no behaviour change), which
%     ports the "[7] Static SCOPE spectral pre-computation" block of legacy src/COMPASS.m:218-269.
%     Static reference data for the RT kernels; no state is advanced.
%
%   Editor: Jizixing Zhu; Date: 2026-07-05.
%   Change: 2026-07-16 Jizixing Zhu — load directional BRDF angle template when
%     options.calc_directional (for surface.radiation.calc_brdf / runStepSIF).
    scope = built.parameters.scope;
    s = struct();
    s.status = "library_loaded";
    s.options = scope.options;
    s.optical_library = scope.optical_library;
    s.spectral = scope.spectral;
    s.canopy_rt = scope.canopy_rt;
    s.library = init.loadScopeStaticLibrary(scope);

    % Static optical precompute (leaf/canopy/soil), computed once per run and read
    % by RTMo/RTMt each step (legacy COMPASS.m:218-269). Merge all returned fields
    % so future additions (leafopt_s/rs_static) flow through without editing here.
    optics = init.precomputeScopeOptics(scope, built.parameters.veg, s.library, s.spectral);
    fn = fieldnames(optics);
    for i = 1:numel(fn)
        s.(fn{i}) = optics.(fn{i});
    end

    % Multi-angle BRDF template (tto/psi) when calc_directional — consumed by
    % surface.radiation.runStepSIF -> surface.radiation.calc_brdf (COMPASS.m / run_stemmus_initialize).
    s.directional = loadDirectionalTemplate(scope, s.optical_library);
end

function directional = loadDirectionalTemplate(scope, optical_library)
%LOADDIRECTIONALTEMPLATE  Load observer angles for calc_directional BRDF.
%
%   Returns empty .noa=0 when calc_directional is off; otherwise reads the
%   site-resolved NetCDF path optical_library.directional_angles_file
%   (from site.toml [model_library], default directional/brdf_angles2.nc
%   under data/model_library_nc/). Variable "data": N x 2 [tto, psi] deg.
% Editor: Jizixing Zhu; Date: 2026-07-16.
% Change: 2026-07-16 Jizixing Zhu — site [model_library] NetCDF only (no .dat).
    directional = struct('tto', [], 'psi', [], 'noa', 0, 'enabled', false);
    options = scope.options;
    if ~(isfield(options, 'calc_directional') && logical(options.calc_directional))
        return;
    end
    if ~isstruct(optical_library) ...
            || ~isfield(optical_library, 'directional_angles_file') ...
            || isempty(optical_library.directional_angles_file)
        error('buildScopeStatic:directionalAnglesUnresolved', ...
            ['calc_directional=true requires model_library.directional_angles_file ', ...
             'resolved under data/model_library_nc/ (format=netcdf).']);
    end
    angle_file = char(optical_library.directional_angles_file);
    if exist(angle_file, 'file') ~= 2
        error('buildScopeStatic:brdfAnglesMissing', ...
            'Directional angles NetCDF not found: %s', angle_file);
    end
    [~, ~, ext] = fileparts(angle_file);
    if ~strcmpi(ext, '.nc')
        error('buildScopeStatic:brdfAnglesNotNetcdf', ...
            'model_library.directional_angles_file must be .nc (format=netcdf), got %s', ...
            angle_file);
    end

    % Same layout as soil_spectrum / leaf_angles NetCDF: variable "data".
    angles_file = double(ncread(angle_file, 'data'));
    if size(angles_file, 1) == 2 && size(angles_file, 2) > 2
        angles_file = angles_file';
    end
    if size(angles_file, 2) < 2 || size(angles_file, 1) < 1
        error('buildScopeStatic:brdfAnglesShape', ...
            'Directional angles must be N x 2 [tto psi], got %s from %s.', ...
            mat2str(size(angles_file)), angle_file);
    end
    directional.tto = double(angles_file(:, 1));
    directional.psi = double(angles_file(:, 2));
    directional.noa = size(angles_file, 1);
    directional.enabled = true;
    directional.source_file = angle_file;
end
