function [spec, legacy] = normalizeSpec(raw)
%NORMALIZESPEC  Build the typed static heat-boundary policy.
%
%   Public configuration names the physical condition directly:
%
%     top    = atmospheric | surface_temperature | specified_heat_flux
%     bottom = seasonal_bottom_temperature | fixed_temperature
%              | specified_heat_flux | free_drainage_internal_energy
%
%   Temperatures use degC and prescribed total-energy fluxes use
%   J cm^-2 s^-1, positive into the soil. Integer NBCT codes are generated
%   only in the returned compatibility adapter; they are not public policy.
%   surface_energy_response is an experimental field-surface policy. It
%   requires the affine radiative/turbulent/snow response assembled by the
%   surface domain and never treats SCOPE's empirical G as a soil boundary.
%   surface_temperature is an explicit diagnostic/legacy operator-split
%   policy and remains a strict frozen Dirichlet condition even when a
%   response is available.
%   atmospheric is the natural surface policy: the surface-energy tangent is
%   matched upstream to the dynamic soil response and materialized as one
%   fixed inward heat flux for the coupled soil step.
%
%   Editor: Jizixing Zhu; Date: 2026-07-16.
%   Change: 2026-07-17 Jizixing Zhu -- document the snow-aware dynamic
%           realization of the exposed-soil surface-temperature policy.
%   Change: 2026-07-17 Jizixing Zhu -- add the ClimaLand-style liquid internal-
%           energy outflow paired with a free-drainage water boundary.
%   Change: 2026-07-18 Jizixing Zhu -- make surface_temperature the strict
%           snow-aware policy: exposed bare soil uses Tss and any resolved
%           snow cover uses the explicit two-tile soil-interface energy flux.
%   Change: 2026-07-19 Jizixing Zhu -- add an accurately named production
%           surface-energy-response policy without changing legacy options.
%   Change: 2026-07-19 Jizixing Zhu -- restore strict surface-temperature
%           semantics as the production operator-split boundary contract.
%   Change: 2026-07-20 Jizixing Zhu -- name the natural coupled heat policy
%           atmospheric, matching the public water-boundary vocabulary.

    removed = {'temp_top','temp_bottom','T_top','T_bottom'};
    if any(isfield(raw, removed))
        error('normalizeHeatBoundarySpec:legacyInputRemoved', ...
            ['Heat boundary integers were removed. Use heat_top_kind, ', ...
             'heat_bottom_kind, heat_top_value, and heat_bottom_value.']);
    end
    required = {'heat_top_kind','heat_bottom_kind', ...
        'heat_top_value','heat_bottom_value'};
    assert(isstruct(raw) && all(isfield(raw, required)), ...
        'normalizeHeatBoundarySpec:missingField', ...
        'Heat boundary config requires %s.', strjoin(required, ', '));

    top = endpoint('top', raw.heat_top_kind, raw.heat_top_value);
    bottom = endpoint('bottom', raw.heat_bottom_kind, raw.heat_bottom_value);
    spec = struct('schema_version', 1, 'provenance', 'flat_kind_v1', ...
        'top', top, 'bottom', bottom);
    legacy = struct('NBCT_top', top.legacy_code, ...
        'NBCT_bottom', bottom.legacy_code, ...
        'T_top', top.configured_value, ...
        'T_bottom', bottom.configured_value);
end

function out = endpoint(side, rawKind, value)
% Map one named heat policy to its retained matrix-adapter code.
% Editor: Jizixing Zhu; Date: 2026-07-16.
    kind = char(string(rawKind));
    assert(isnumeric(value) && isscalar(value) && isfinite(value), ...
        'normalizeHeatBoundarySpec:badValue', ...
        'Heat %s boundary value must be a finite scalar.', side);
    if strcmp(side, 'top')
        kinds = {'atmospheric','surface_temperature', ...
            'specified_heat_flux','surface_energy_response'};
    else
        kinds = {'seasonal_bottom_temperature','fixed_temperature', ...
            'specified_heat_flux','free_drainage_internal_energy'};
    end
    kindIndex = find(strcmp(kind, kinds), 1);
    assert(~isempty(kindIndex), 'normalizeHeatBoundarySpec:badKindForPosition', ...
        'Unsupported heat %s kind "%s".', side, kind);
    code = kindIndex;
    if strcmp(kind, 'atmospheric')
        code = 2;
    elseif strcmp(kind, 'surface_energy_response')
        code = 3;
    end
    out = struct('kind', kind, 'configured_value', double(value), ...
        'legacy_code', code);
end
