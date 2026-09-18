function writeScopeStep(registry, built, kt)
%WRITESCOPESTEP  Export one accepted SCOPE radiation and energy-balance step.
%
%   The canopy solver has already evaluated the four SCOPE output structs:
%   fluxes, thermal, profiles, and rad. This function only maps those accepted
%   values to the diagnostics registry; it does not rerun radiative transfer,
%   photosynthesis, or the canopy energy balance.
%
%   Surface-energy closure follows the same equation as legacy COMPASS:
%
%       residual = Rn - LE - H - G.                         [W m^-2]
%
%   Red-edge arrays use the physical 590--850 nm wavelength interval rather
%   than hard-coded array positions. SIF is zero below its 640 nm support.
%   Editor: Jizixing Zhu; Date: 2026-07-16.
%   Change: 2026-07-18 Jizixing Zhu -- pass normalized process switches so
%     unavailable ecosystem respiration remains a missing diagnostic.
    if ~isfield(built.workspace, 'scope_step') ...
            || isempty(built.workspace.scope_step)
        return;
    end
    scope = built.workspace.scope_step;
    fluxes = readStruct(scope, 'fluxes');
    thermal = readStruct(scope, 'thermal');
    profiles = readStruct(scope, 'profiles');
    rad = readStruct(scope, 'rad');

    diagnostics.writeScopeStructs( ...
        registry, fluxes, thermal, profiles, rad, kt, ...
        built.control.switches);
    writeEnergyClosure(registry, fluxes, kt);
    writeRedEdgeSpectra(registry, scope, rad, kt);

    registeredNames = registry.list();
    onlineGsRegistered = any(strcmp( ...
        registeredNames, 'medlyn_gs_model_sun'));
    if onlineGsRegistered && isfield(fluxes, 'online_gs_diag') ...
            && isstruct(fluxes.online_gs_diag)
        diagnostics.writeOnlineGsValidation( ...
            registry, fluxes.online_gs_diag, kt);
    end
end

function writeEnergyClosure(registry, fluxes, kt)
% Write Rn - LE - H - G only when all four accepted fluxes are finite.
    names = {'Rntot', 'lEtot', 'Htot', 'Gtot'};
    values = nan(1, numel(names));
    for i = 1:numel(names)
        if ~isfield(fluxes, names{i}) || ~isscalar(fluxes.(names{i}))
            return;
        end
        values(i) = fluxes.(names{i});
    end
    if all(isfinite(values))
        residual = values(1) - values(2) - values(3) - values(4);
        registry.write('energy_balance_residual', residual, kt);
    end
end

function writeRedEdgeSpectra(registry, scope, rad, kt)
% Map full SCOPE spectra to the registered 590--850 nm output grid.
    if ~isfield(scope, 'scope_static') ...
            || ~isfield(scope.scope_static, 'spectral')
        return;
    end
    spectral = scope.scope_static.spectral;
    if ~isfield(spectral, 'wlS')
        return;
    end
    wavelength = spectral.wlS(:);
    redIndex = find(wavelength >= 590 & wavelength <= 850);
    if numel(redIndex) ~= 261
        return;
    end
    tryWriteBand(registry, 'toc_radiance_re', rad, 'Lo_', redIndex, kt);
    tryWriteBand(registry, 'esun_re', rad, 'Esun_', redIndex, kt);
    tryWriteBand(registry, 'esky_re', rad, 'Esky_', redIndex, kt);

    if ~isfield(rad, 'LoF_') || isempty(rad.LoF_) ...
            || ~isfield(spectral, 'wlF')
        return;
    end
    fluorescenceWavelength = spectral.wlF(:);
    fluorescence = rad.LoF_(:);
    if numel(fluorescenceWavelength) ~= numel(fluorescence)
        return;
    end
    sif = zeros(numel(redIndex), 1);
    [matchedWavelength, redLocation, fluorescenceLocation] = intersect( ...
        wavelength(redIndex), fluorescenceWavelength, 'stable'); %#ok<ASGLU>
    sif(redLocation) = fluorescence(fluorescenceLocation);
    registry.write('sif_radiance_re', sif, kt);
end

function tryWriteBand(registry, outputName, source, fieldName, index, kt)
% Write one red-edge band when its full-spectrum source is available.
    if ~isfield(source, fieldName) || isempty(source.(fieldName))
        return;
    end
    value = source.(fieldName);
    value = value(:);
    if numel(value) >= max(index)
        registry.write(outputName, value(index), kt);
    end
end

function value = readStruct(container, name)
% Return a struct field or an empty struct for optional SCOPE branches.
    value = struct();
    if isfield(container, name) && isstruct(container.(name))
        value = container.(name);
    end
end
