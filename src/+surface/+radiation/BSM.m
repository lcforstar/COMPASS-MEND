function rwet = BSM(soilpar, spec, emp)
%BSM  Compute wet soil reflectance with the SCOPE brightness-shape model.
%
%   rwet = surface.radiation.BSM(soilpar, spec, emp) returns the optical soil
%   reflectance spectrum used by SCOPE when options.soilspectrum == 1.
%
%   Inputs:
%     soilpar [in] struct with:
%       .BSMBrightness [-]        Soil brightness scaling.
%       .BSMlat        [deg]      BSM spectral-shape latitude.
%       .BSMlon        [deg]      BSM spectral-shape longitude.
%       .SMC           [m3/m3]    Surface volumetric soil liquid water.
%     spec [in] struct with:
%       .GSV           [-, nwl x 3] Global Soil Vectors spectra.
%       .Kw            [-, nwl x 1] Water absorption spectrum.
%       .nw            [-, nwl x 1] Water refraction index spectrum.
%     emp [in] struct with:
%       .SMC_capacity_percent [%] SCOPE BSM water-film capacity parameter.
%                                Legacy field .SMC is also accepted.
%       .film          [-]        Optical thickness of one water film.
%
%   Output:
%     rwet             [-, nwl x 1] Wet soil reflectance spectrum.
%
%   Method / reference:
%     Port of SCOPE-2.0 src/RTMs/BSM.m from the local reference checkout
%     D:/OneDrive/.../Model/SCOPE-2.0. The soil-water film calculation follows
%     the SCOPE BSM implementation by Wout Verhoef (version 1.0, September
%     2012) with the Peiqi Yang January 2020 truncation at k = 0:6 water films.
%     This function is a pure optical auxiliary: it does not update model state.
%
%   Editor: Jizixing Zhu; Date: 2026-07-08.

    requireFields(soilpar, {'BSMBrightness','BSMlat','BSMlon','SMC'}, 'soilpar');
    requireFields(spec, {'GSV','Kw','nw'}, 'spec');
    requireEmpFields(emp);

    GSV = double(spec.GSV);
    kw = double(spec.Kw(:));
    nw = double(spec.nw(:));
    if size(GSV, 2) < 3 || size(GSV, 1) ~= numel(kw) || numel(kw) ~= numel(nw)
        error('Radiation:BSM:badSpectralShape', ...
            'spec.GSV must be [nwl x 3+] and match spec.Kw/spec.nw length.');
    end

    B = double(soilpar.BSMBrightness);
    lat = double(soilpar.BSMlat);
    lon = double(soilpar.BSMlon);
    SMp = double(soilpar.SMC) * 1e2;   % [m3/m3] -> volume percent, as in SCOPE BSM

    SMC_capacity_percent = double(getEmpCapacityPercent(emp));
    film = double(emp.film);
    if ~isscalar(SMp) || ~isfinite(SMp) || ~isscalar(SMC_capacity_percent) ...
            || SMC_capacity_percent <= 0 || ~isfinite(SMC_capacity_percent)
        error('Radiation:BSM:badMoistureInput', ...
            'soilpar.SMC must be finite scalar and emp.SMC_capacity_percent must be positive scalar.');
    end
    if ~isscalar(film) || ~isfinite(film) || film <= 0
        error('Radiation:BSM:badFilmInput', ...
            'emp.film must be a positive finite scalar.');
    end

    f1 = B * sind(lat);
    f2 = B * cosd(lat) * sind(lon);
    f3 = B * cosd(lat) * cosd(lon);
    rdry = f1 * GSV(:, 1) + f2 * GSV(:, 2) + f3 * GSV(:, 3);

    rwet = soilwat(rdry, nw, kw, SMp, SMC_capacity_percent, film);
end

% -------------------------------------------------------------------------
function requireEmpFields(emp)
% Check the empirical BSM parameter contract.
    if ~isstruct(emp)
        error('Radiation:BSM:badInput', 'emp must be a struct.');
    end
    if ~isfield(emp, 'SMC_capacity_percent') && ~isfield(emp, 'SMC')
        error('Radiation:BSM:missingField', ...
            'emp.SMC_capacity_percent is required; legacy emp.SMC is also accepted.');
    end
    if ~isfield(emp, 'film')
        error('Radiation:BSM:missingField', 'emp.film is required.');
    end
end

% -------------------------------------------------------------------------
function value = getEmpCapacityPercent(emp)
% Return the SCOPE BSM water-film capacity parameter in percent units.
    if isfield(emp, 'SMC_capacity_percent')
        value = emp.SMC_capacity_percent;
    else
        value = emp.SMC;
    end
end

% -------------------------------------------------------------------------
function rwet = soilwat(rdry, nw, kw, SMp, SMC_capacity_percent, deleff)
% Apply the SCOPE BSM water-film correction to dry soil reflectance.

    k = 0:6;                        % number of water films; 0 is dry soil
    mu = (SMp - 5) / SMC_capacity_percent;  % Poisson mean used by SCOPE BSM
    if mu <= 0
        rwet = rdry;
        return
    end

    % SCOPE-2.0 src/RTMs/BSM.m soilwat: background reflectance and water-film
    % Fresnel terms. The local poissonWeights call replaces poisspdf to avoid a
    % Statistics Toolbox dependency while preserving the k = 0:6 truncation.
    rbac = 1 - (1 - rdry) .* (rdry .* tav(90, 2.0 ./ nw) / tav(90, 2.0) + 1 - rdry);
    p = 1 - tav(90, nw) ./ nw.^2;
    Rw = 1 - tav(40, nw);

    fmul = poissonWeights(k, mu);
    tw = exp(-2 * kw * deleff .* k);
    Rwet_k = Rw + (1 - Rw) .* (1 - p) .* tw .* rbac ./ (1 - p .* tw .* rbac);
    rwet = rdry * fmul(1) + Rwet_k(:, 2:end) * fmul(2:end);
end

% -------------------------------------------------------------------------
function w = poissonWeights(k, mu)
% Return Poisson probabilities for integer k and scalar mu.
    w = exp(-mu) * (mu .^ k) ./ factorial(k);
    w = w(:);
end

% -------------------------------------------------------------------------
function Tav = tav(alfa, nr)
% Diffuse transmittance through an interface, copied from the SCOPE BSM code.

    n2 = nr.^2;
    np = n2 + 1;
    nm = n2 - 1;

    a = ((nr + 1).^2) / 2;
    k = -((n2 - 1).^2) / 4;
    sin_a = sind(alfa);

    if alfa ~= 0
        B2 = sin_a^2 - np / 2;
        B1 = (alfa ~= 90) * sqrt(B2.^2 + k);

        b = B1 - B2;
        b3 = b.^3;
        a3 = a.^3;

        ts = (k.^2 ./ (6 * b3) + k ./ b - b ./ 2) - ...
            (k.^2 ./ (6 * a3) + k ./ a - a ./ 2);

        tp1 = -2 * n2 .* (b - a) ./ (np.^2);
        tp2 = -2 * n2 .* np .* log(b ./ a) ./ (nm.^2);
        tp3 = n2 .* (1 ./ b - 1 ./ a) ./ 2;
        tp4 = 16 * n2.^2 .* (n2.^2 + 1) .* ...
            log((2 * np .* b - nm.^2) ./ (2 * np .* a - nm.^2)) ./ (np.^3 .* nm.^2);
        tp5 = 16 * n2.^2 .* n2 .* ...
            (1 ./ (2 * np .* b - nm.^2) - 1 ./ (2 * np .* a - nm.^2)) ./ (np.^3);
        tp = tp1 + tp2 + tp3 + tp4 + tp5;
        Tav = (ts + tp) ./ (2 * sin_a.^2);
    else
        Tav = 4 * nr ./ ((nr + 1) .* (nr + 1));
    end
end

% -------------------------------------------------------------------------
function requireFields(s, names, owner)
% Check required struct fields.
    if ~isstruct(s)
        error('Radiation:BSM:badInput', '%s must be a struct.', owner);
    end
    for i = 1:numel(names)
        if ~isfield(s, names{i})
            error('Radiation:BSM:missingField', '%s.%s is required.', owner, names{i});
        end
    end
end
