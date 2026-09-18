function [spec, legacy] = normalizeSpec(raw, useHydrostaticBottom)
%NORMALIZESPEC  Build the typed static dry-air boundary policy.
%
%   Public pressure values use Pa. Prescribed flux values are dry-air mass
%   fluxes [g cm^-2 s^-1], positive into the soil; they are never pressure
%   rates. Integer NBCP codes exist only in the returned compatibility view.
%
%   Editor: Jizixing Zhu; Date: 2026-07-16.

    removed = {'air_top','air_bottom','P_top','P_bottom'};
    if any(isfield(raw, removed))
        error('normalizeDryAirBoundarySpec:legacyInputRemoved', ...
            ['Dry-air boundary integers were removed. Use dry_air_top_kind, ', ...
             'dry_air_bottom_kind, dry_air_top_value, and dry_air_bottom_value.']);
    end
    required = {'dry_air_top_kind','dry_air_bottom_kind', ...
        'dry_air_top_value','dry_air_bottom_value'};
    assert(isstruct(raw) && all(isfield(raw, required)), ...
        'normalizeDryAirBoundarySpec:missingField', ...
        'Dry-air boundary config requires %s.', strjoin(required, ', '));

    top = endpoint('top', raw.dry_air_top_kind, raw.dry_air_top_value);
    bottom = endpoint('bottom', raw.dry_air_bottom_kind, ...
        raw.dry_air_bottom_value);
    bottom.use_hydrostatic_atmosphere = logical(useHydrostaticBottom);
    spec = struct('schema_version', 1, 'provenance', 'flat_kind_v1', ...
        'top', top, 'bottom', bottom);
    legacy = struct('NBCP_top', top.legacy_code, ...
        'NBCP_bottom', bottom.legacy_code, ...
        'P_top', top.configured_value, ...
        'P_bottom', bottom.configured_value);
end

function out = endpoint(side, rawKind, value)
% Map one named dry-air policy to its retained matrix-adapter code.
% Editor: Jizixing Zhu; Date: 2026-07-16.
    kind = char(string(rawKind));
    assert(isnumeric(value) && isscalar(value) && isfinite(value), ...
        'normalizeDryAirBoundarySpec:badValue', ...
        'Dry-air %s boundary value must be a finite scalar.', side);
    if strcmp(side, 'top')
        kinds = {'fixed_pressure','specified_mass_flux','atmospheric_pressure'};
    else
        kinds = {'hydrostatic_pressure','specified_mass_flux','zero_mass_flux'};
    end
    code = find(strcmp(kind, kinds), 1);
    assert(~isempty(code), 'normalizeDryAirBoundarySpec:badKindForPosition', ...
        'Unsupported dry-air %s kind "%s".', side, kind);
    out = struct('kind', kind, 'configured_value', double(value), ...
        'legacy_code', code);
end
