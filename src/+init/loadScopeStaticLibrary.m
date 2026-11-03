function library = loadScopeStaticLibrary(scope)
%LOADSCOPESTATICLIBRARY  Read SCOPE static optical resources from NetCDF files.
%
%   library = init.loadScopeStaticLibrary(scope) loads the shared model-library
%   NetCDF resources resolved by init.buildScope:
%       leaf optics  -> library.optipar  (PROSPECT/Fluspect constants)
%       atmosphere   -> library.atmo.M   (SCOPE-band aggregated MODTRAN columns)
%       soil spectrum -> library.soil.rsfile / rs when scope.options.soilspectrum=0
%
%   Inputs:
%     scope [struct]  built.parameters.scope with fields:
%       optical_library.leaf_optics_file      [.nc]
%       optical_library.atmosphere_file       [.nc]
%       optical_library.soil_spectrum_file    [.nc]
%       spectral.SCOPEspec                    SCOPE aggregation bands
%       options.soilspectrum                  0=file, 1=BSM model
%       soil_optical.spectrum                 soil reflectance column index
%
%   Returns:
%     library [struct] runtime static optical data. This is read-only input for
%     RT/SCOPE kernels; physical code should not mutate it.
%
%   Closure:
%     No state is advanced. This function only materializes immutable shared
%     reference data that older code loaded with load()/importdata().
%
%   Editor: Jizixing Zhu; Date: 2026-07-03.
%   Change: 2026-07-03 Jizixing Zhu — read optional library.leaf_angles from ol.leaf_angles_file (file-based LIDF); empty when the site uses analytic LIDFa/LIDFb.

    requireScopeFields(scope);
    ol = scope.optical_library;

    library = struct();
    library.status = "loaded";
    library.files = ol;
    library.optipar = readLeafOptics(ol.leaf_optics_file);
    library.atmo = struct('M', readAtmosphere(ol.atmosphere_file, scope.spectral.SCOPEspec));

    if scope.options.soilspectrum == 0
        rsfile = readMatrixData(ol.soil_spectrum_file);
        idx = scope.soil_optical.spectrum + 1;  % legacy: column 1 is wavelength
        if idx < 1 || idx > size(rsfile, 2)
            error('loadScopeStaticLibrary:badSoilSpectrumIndex', ...
                'soil spectrum column %d is outside rsfile width %d.', idx, size(rsfile, 2));
        end
        library.soil = struct('loaded', true, 'rsfile', rsfile, ...
            'spectrum_index', scope.soil_optical.spectrum, 'rs', rsfile(:, idx));
    else
        library.soil = struct('loaded', false, 'reason', "soilspectrum != 0");
    end

    % Optional measured leaf inclination distribution (file-based LIDF). Same 'data' NetCDF
    % variable as the soil file; returns the [nangles x 1] lidf table (13 bins). Empty when the
    % site uses the analytic LIDFa/LIDFb parameters (precomputeScopeOptics picks the branch).
    if isfield(ol, 'leaf_angles_file') && ~isempty(ol.leaf_angles_file)
        library.leaf_angles = readMatrixData(ol.leaf_angles_file);
    else
        library.leaf_angles = [];
    end
end

% -------------------------------------------------------------------------
function requireScopeFields(scope)
% Check the minimal scope contract needed to load static SCOPE resources.
    req = {'optical_library','spectral','options','soil_optical'};
    for i = 1:numel(req)
        if ~isfield(scope, req{i})
            error('loadScopeStaticLibrary:missingScopeField', ...
                'scope.%s is required.', req{i});
        end
    end
    if ~isfield(scope.spectral, 'SCOPEspec')
        error('loadScopeStaticLibrary:missingSpectralSpec', ...
            'scope.spectral.SCOPEspec is required for atmosphere aggregation.');
    end
end

% -------------------------------------------------------------------------
function optipar = readLeafOptics(ncPath)
% Read all variables from a PROSPECT/Fluspect NetCDF into an optipar struct.
    info = ncinfo(ncPath);
    optipar = struct();
    for i = 1:numel(info.Variables)
        name = info.Variables(i).Name;
        value = double(ncread(ncPath, name));
        value = orientOpticalValue(value);
        optipar.(name) = value;
    end
end

% -------------------------------------------------------------------------
function value = orientOpticalValue(value)
% Match legacy MAT orientation: wavelength vectors are column vectors; GSV is
% [nw x 3], not [3 x nw].
    if isvector(value)
        value = value(:);
    elseif size(value, 1) < size(value, 2)
        value = value.';
    end
end

% -------------------------------------------------------------------------
function data = readMatrixData(ncPath)
% Read a generic NetCDF matrix stored as variable "data" and orient rows as
% records/wavelengths, matching legacy loadtxt/importdata output.
    data = double(ncread(ncPath, 'data'));
    if size(data, 1) < size(data, 2)
        data = data.';
    end
end

% -------------------------------------------------------------------------
function M = readAtmosphere(ncPath, SCOPEspec)
% Read FLEX-S3/MODTRAN NetCDF and aggregate to SCOPE wavelengths.
    atm = readMatrixData(ncPath);      % legacy .atm numeric table, [nrow x 20]
    if size(atm, 2) < 20
        error('loadScopeStaticLibrary:badAtmosphereShape', ...
            'Atmosphere table must have at least 20 columns; got %d.', size(atm, 2));
    end

    wlM = atm(:, 2);
    T = atm(:, 3:20);
    U = [T(:, 1), T(:, 3), T(:, 4), T(:, 5), T(:, 12), T(:, 16)];

    nreg = SCOPEspec.nreg;
    streg = SCOPEspec.start;
    enreg = SCOPEspec.end;
    width = SCOPEspec.res;
    nwreg = int32((enreg - streg) ./ width) + 1;
    off = int32(zeros(nreg, 1));
    for i = 2:nreg
        off(i) = off(i - 1) + nwreg(i - 1);
    end

    nwS = sum(nwreg);
    n = zeros(nwS, 1);
    S = zeros(nwS, 6);
    j = int32(zeros(nreg, 1));
    for iwl = 1:numel(wlM)
        w = wlM(iwl);
        for r = 1:nreg
            j(r) = int32(round(w - streg(r)) ./ width(r)) + 1;
            if j(r) > 0 && j(r) <= nwreg(r)
                k = j(r) + off(r);
                S(k, :) = S(k, :) + U(iwl, :);
                n(k) = n(k) + 1;
            end
        end
    end

    M = zeros(size(S, 1), 6);
    for i = 1:6
        M(:, i) = S(:, i) ./ n;
    end
end
