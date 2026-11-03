function [spec, legacy] = normalizeSpec(rawBoundary)
%NORMALIZESPEC  Validate the single typed water-boundary configuration.
%
%   [spec, legacy] = soil.water.boundary.normalizeSpec(rawBoundary) creates
%   one canonical water-boundary definition while generating the flat aliases
%   required internally by the transitional Richards adapter. Current site
%   files use one mode per endpoint: constant_head, variable_head,
%   constant_flux, variable_flux, atmospheric/zero_flux, or free_drainage.
%   All prescribed values come from the fixed top_value/bottom_value CSV
%   columns. Constant modes require a constant column; variable head uses
%   linear interpolation and variable flux uses previous interpolation.
%   Historical integer fields, the
%   earlier split kind/source schema, and the nested typed draft are rejected.
%   This function performs configuration validation and unit/sign translation
%   only; it does not evaluate forcing, update state, or solve a boundary flux.
%
%   Inputs:
%     rawBoundary [in] scalar struct from raw.config.boundary. Required fields
%                       are water_top_kind and water_bottom_kind. Prescribed
%                       values are external forcing in water_boundary.csv.
%
%   Returns:
%     spec   canonical plain-data struct. top/bottom each contain kind,
%            source, legacy_code [-], head_cm [cm],
%            flux_into_domain_cm_s [cm/s], series_column and interpolation.
%            Positive public flux always means into the soil domain. A
%            negative top flux therefore prescribes upward removal from the
%            column, while a negative bottom flux prescribes downward outflow.
%     legacy transitional struct with NBCh_top/NBCh_bottom [-] and
%            h_top/h_bottom in the legacy solver sign/value convention.
%
%   Solver sign adapter after CSV materialization:
%     top legacy BCh    = -flux_into_domain_cm_s;
%     bottom legacy BCh =  flux_into_domain_cm_s.
%   The conversion is confined to this internal compatibility boundary.
%
%   Architecture reference:
%     ClimaLand.jl src/standalone/Soil/boundary_conditions.jl:25-98
%     (separate water-BC policies) and src/shared_utilities/drivers.jl:44-83
%     (drivers are separate from boundary-condition policy).
%
%   Editor: Jizixing Zhu; Date: 2026-07-11.
%   Change: 2026-07-11 Jizixing Zhu -- expose one flat typed [boundary] schema only; reject historical integer and nested draft inputs while retaining generated legacy solver aliases.
%   Change: 2026-07-11 Jizixing Zhu -- allow inactive prepared values to remain beside the active selector and ignore the configured CSV path during policy normalization.
%   Change: 2026-07-11 Jizixing Zhu -- collapse the public schema to one HYDRUS-style mode per endpoint and derive series columns/interpolation automatically.
%   Change: 2026-07-11 Jizixing Zhu -- make top_value/bottom_value CSV columns the only prescribed-value source, including constant modes.
%   Change: 2026-07-11 Jizixing Zhu -- document signed prescribed outflow under the single into-domain-positive convention.

    if ~isstruct(rawBoundary) || ~isscalar(rawBoundary)
        error('normalizeWaterBoundarySpec:badInput', ...
            'rawConfig.boundary must be a scalar struct.');
    end

    legacyFields = {'water_top','water_bottom','h_top','h_bottom', ...
        'q_top_cm_s','q_bottom_cm_s'};
    hasLegacy = any(isfield(rawBoundary, legacyFields));
    if hasLegacy
        error('normalizeWaterBoundarySpec:legacyInputRemoved', ...
            ['Historical flat water fields are no longer accepted in site.toml. ' ...
             'Migrate to water_top_* and water_bottom_* fields in [boundary].']);
    end
    if isfield(rawBoundary, 'water')
        error('normalizeWaterBoundarySpec:nestedInputRemoved', ...
            ['Nested [boundary.water.top/bottom] input is no longer accepted. ' ...
             'Keep all water_top_* and water_bottom_* fields in [boundary].']);
    end
    splitFields = {'water_top_source','water_top_series_column', ...
        'water_top_interpolation','water_bottom_source', ...
        'water_bottom_series_column','water_bottom_interpolation'};
    presentSplit = splitFields(isfield(rawBoundary, splitFields));
    if ~isempty(presentSplit)
        error('normalizeWaterBoundarySpec:splitInputRemoved', ...
            ['The split source/column/interpolation input has been removed. ' ...
             'Select constant_head, variable_head, constant_flux, or ' ...
             'variable_flux with water_top_kind/water_bottom_kind. ' ...
             'Remove field(s): %s.'], strjoin(presentSplit, ', '));
    end
    required = {'water_top_kind','water_bottom_kind'};
    if any(~isfield(rawBoundary, required))
        error('normalizeWaterBoundarySpec:missingSchema', ...
            '[boundary] requires water_top_kind and water_bottom_kind.');
    end
    checkUnknownWaterFields(rawBoundary);
    water = struct( ...
        'top', collectEndpoint(rawBoundary, 'top'), ...
        'bottom', collectEndpoint(rawBoundary, 'bottom'));
    [spec, legacy] = normalizeTyped(water);
end

% -------------------------------------------------------------------------
function endpoint = collectEndpoint(boundary, position)
% Convert one public endpoint mode into the canonical internal policy fields.
    modeField = ['water_' position '_kind'];
    mode = readToken(boundary.(modeField), modeField);
    switch mode
        case 'constant_head'
            endpoint = struct('kind','fixed_head','source','series', ...
                'series_column',[position '_value'],'interpolation','constant');
        case 'variable_head'
            endpoint = struct('kind','fixed_head','source','series', ...
                'series_column',[position '_value'],'interpolation','linear');
        case 'constant_flux'
            endpoint = struct('kind','fixed_flux','source','series', ...
                'series_column',[position '_value'],'interpolation','constant');
        case 'variable_flux'
            endpoint = struct('kind','fixed_flux','source','series', ...
                'series_column',[position '_value'],'interpolation','previous');
        case 'zero_flux'
            endpoint = struct('kind','zero_flux','source','none');
        case 'atmospheric'
            if ~strcmp(position, 'top')
                invalidPublicMode(mode, position);
            end
            endpoint = struct('kind','atmospheric','source','natural_surface');
        case 'free_drainage'
            if ~strcmp(position, 'bottom')
                invalidPublicMode(mode, position);
            end
            endpoint = struct('kind','free_drainage','source','none');
        otherwise
            invalidPublicMode(mode, position);
    end
end

% -------------------------------------------------------------------------
function invalidPublicMode(mode, position)
% Report the complete public mode set for one endpoint position.
    if strcmp(position, 'top')
        allowed = ['atmospheric, constant_head, variable_head, ' ...
            'constant_flux, variable_flux, zero_flux'];
    else
        allowed = ['free_drainage, constant_head, variable_head, ' ...
            'constant_flux, variable_flux, zero_flux'];
    end
    error('normalizeWaterBoundarySpec:badKindForPosition', ...
        'Water boundary mode "%s" is invalid at the %s boundary; expected {%s}.', ...
        mode, position, allowed);
end

% -------------------------------------------------------------------------
function checkUnknownWaterFields(boundary)
% Fail on misspelled water_* fields while leaving other process BCs untouched.
    allowed = { ...
        'water_boundary_file', ...
        'water_top_kind','water_bottom_kind'};
    names = fieldnames(boundary);
    waterNames = names(startsWith(names, 'water_'));
    unknown = setdiff(waterNames, allowed, 'stable');
    if ~isempty(unknown)
        error('normalizeWaterBoundarySpec:unknownField', ...
            'Unknown [boundary] water field(s): %s.', strjoin(unknown, ', '));
    end
end

% -------------------------------------------------------------------------
function [spec, legacy] = normalizeTyped(water)
% Validate constants, explicit series metadata, natural atmospheric, and no-flux policies.
    if ~isstruct(water) || ~isscalar(water)
        error('normalizeWaterBoundarySpec:badTypedSchema', ...
            'boundary.water must be a scalar struct.');
    end
    checkUnknownFields(water, {'top','bottom'}, 'boundary.water');
    requireFields(water, {'top','bottom'}, 'boundary.water');

    top = typedEndpoint(water.top, 'top');
    bottom = typedEndpoint(water.bottom, 'bottom');
    spec = makeSpec('flat_mode_v2', top, bottom);
    legacy = makeLegacy(top.legacy_code, bottom.legacy_code, ...
        endpointLegacyValue(top, 'top'), endpointLegacyValue(bottom, 'bottom'));
end

% -------------------------------------------------------------------------
function endpoint = typedEndpoint(raw, position)
% Validate one tagged endpoint; unsupported sources fail instead of falling back.
    if ~isstruct(raw) || ~isscalar(raw)
        error('normalizeWaterBoundarySpec:badEndpoint', ...
            'boundary.water.%s must be a scalar struct.', position);
    end
    allowedFields = {'kind','source','head_cm','flux_into_domain_cm_s', ...
        'series_column','interpolation'};
    checkUnknownFields(raw, allowedFields, ['boundary.water.' position]);
    requireFields(raw, {'kind','source'}, ['boundary.water.' position]);

    kind = readToken(raw.kind, [position '.kind']);
    source = readToken(raw.source, [position '.source']);
    validateKindPosition(kind, position);
    validatePreparedFields(raw, position);

    switch kind
        case 'fixed_head'
            if strcmp(source, 'constant')
                endpoint = makeEndpoint(kind, source, kindCode(kind, position), ...
                    readRequiredFinite(raw, 'head_cm'), NaN);
            elseif strcmp(source, 'series')
                [column, interpolation] = readSeriesMeta(raw, position, kind, ...
                    {'constant','previous','linear'});
                endpoint = makeEndpoint(kind, source, kindCode(kind, position), ...
                    NaN, NaN, column, interpolation);
            else
                unsupportedSource(source, position, kind, {'constant','series'});
            end
        case 'fixed_flux'
            if strcmp(source, 'constant')
                endpoint = makeEndpoint(kind, source, kindCode(kind, position), NaN, ...
                    readRequiredFinite(raw, 'flux_into_domain_cm_s'));
            elseif strcmp(source, 'series')
                [column, interpolation] = readSeriesMeta(raw, position, kind, ...
                    {'constant','previous'});
                endpoint = makeEndpoint(kind, source, kindCode(kind, position), ...
                    NaN, NaN, column, interpolation);
            else
                unsupportedSource(source, position, kind, {'constant','series'});
            end
        case 'atmospheric'
            requireSource(source, 'natural_surface', position, kind);
            endpoint = makeEndpoint(kind, source, kindCode(kind, position), NaN, NaN);
        case {'zero_flux','free_drainage'}
            requireSource(source, 'none', position, kind);
            flux = NaN;
            if strcmp(kind, 'zero_flux'), flux = 0; end
            endpoint = makeEndpoint(kind, source, kindCode(kind, position), NaN, flux);
        otherwise
            error('normalizeWaterBoundarySpec:badKind', ...
                'Unsupported water-boundary kind "%s".', kind);
    end
end

% -------------------------------------------------------------------------
function spec = makeSpec(provenance, top, bottom)
% Assemble immutable plain data; no function handles are stored in config/state.
    spec = struct('schema_version', 1, ...
        'provenance', provenance, ...
        'sign_convention', 'flux_into_domain_positive', ...
        'top', top, 'bottom', bottom);
end

% -------------------------------------------------------------------------
function endpoint = makeEndpoint(kind, source, legacyCode, head, flux, seriesColumn, interpolation)
% Use a fixed field layout so MATLAB structs map directly to future Julia data.
    if nargin < 6, seriesColumn = ''; end
    if nargin < 7, interpolation = 'none'; end
    endpoint = struct('kind', kind, 'source', source, ...
        'legacy_code', legacyCode, ...
        'head_cm', head, ...
        'flux_into_domain_cm_s', flux, ...
        'series_column', seriesColumn, ...
        'interpolation', interpolation);
end

% -------------------------------------------------------------------------
function legacy = makeLegacy(topCode, bottomCode, topValue, bottomValue)
% Transitional aliases consumed by the existing solver adapter.
    legacy = struct('NBCh_top', topCode, 'NBCh_bottom', bottomCode, ...
        'h_top', topValue, 'h_bottom', bottomValue);
end

% -------------------------------------------------------------------------
function value = endpointLegacyValue(endpoint, position)
% Convert canonical value to the current positional BCh convention.
    switch endpoint.kind
        case 'fixed_head'
            value = endpoint.head_cm;
        case 'fixed_flux'
            value = toLegacyFlux(endpoint.flux_into_domain_cm_s, position);
        otherwise
            value = 0;
    end
end

% -------------------------------------------------------------------------
function value = toLegacyFlux(fluxInto, position)
% Top legacy flux is upward-positive; bottom legacy flux is into-domain-positive.
    if strcmp(position, 'top')
        value = -fluxInto;
    else
        value = fluxInto;
    end
end

% -------------------------------------------------------------------------
function code = kindCode(kind, position)
% Map canonical policy names to transitional integer codes.
    switch kind
        case 'fixed_head', code = 1;
        case 'fixed_flux', code = 2;
        case {'atmospheric','free_drainage'}, code = 3;
        case 'zero_flux', code = 4;
        otherwise
            error('normalizeWaterBoundarySpec:badKind', ...
                'Unsupported %s water-boundary kind "%s".', position, kind);
    end
end

% -------------------------------------------------------------------------
function validateKindPosition(kind, position)
% Top and bottom endpoints intentionally support different policy sets.
    if strcmp(position, 'top')
        allowed = {'fixed_head','fixed_flux','zero_flux','atmospheric'};
    else
        allowed = {'fixed_head','fixed_flux','zero_flux','free_drainage'};
    end
    if ~ismember(kind, allowed)
        error('normalizeWaterBoundarySpec:badKindForPosition', ...
            'Water-boundary kind "%s" is invalid at the %s boundary.', ...
            kind, position);
    end
end

% -------------------------------------------------------------------------
function requireSource(actual, expected, position, kind)
% Require the one source policy supported by this endpoint kind.
    if ~strcmp(actual, expected)
        error('normalizeWaterBoundarySpec:unsupportedSource', ...
            'boundary.water.%s kind "%s" requires source "%s", got "%s".', ...
            position, kind, expected, actual);
    end
end

% -------------------------------------------------------------------------
function unsupportedSource(actual, position, kind, allowed)
% Report the full supported source set for a typed endpoint.
    error('normalizeWaterBoundarySpec:unsupportedSource', ...
        'boundary.water.%s kind "%s" source must be one of {%s}, got "%s".', ...
        position, kind, strjoin(allowed, ', '), actual);
end

% -------------------------------------------------------------------------
function [column, interpolation] = readSeriesMeta(raw, position, kind, allowedInterpolation)
% Validate explicit series ownership and interpolation policy.
    requireFields(raw, {'series_column','interpolation'}, ...
        ['boundary.water.' position]);
    column = readName(raw.series_column, [position '.series_column']);
    interpolation = readToken(raw.interpolation, [position '.interpolation']);
    if ~ismember(interpolation, allowedInterpolation)
        error('normalizeWaterBoundarySpec:badInterpolation', ...
            ['boundary.water.%s kind "%s" interpolation must be one of ' ...
             '{%s}, got "%s".'], position, kind, ...
            strjoin(allowedInterpolation, ', '), interpolation);
    end
end

% -------------------------------------------------------------------------
function validatePreparedFields(raw, position)
% Validate stored alternatives without making them part of the active endpoint.
    numericFields = {'head_cm','flux_into_domain_cm_s'};
    for i = 1:numel(numericFields)
        field = numericFields{i};
        if isfield(raw, field) && ~isempty(raw.(field))
            value = raw.(field);
            if ~isnumeric(value) || ~isscalar(value) || ~isfinite(value)
                error('normalizeWaterBoundarySpec:badValue', ...
                    'water_%s_%s must be a finite numeric scalar.', position, field);
            end
        end
    end
    if isfield(raw, 'series_column') && ~isempty(raw.series_column)
        readName(raw.series_column, [position '.series_column']);
    end
    if isfield(raw, 'interpolation') && ~isempty(raw.interpolation)
        interpolation = readToken(raw.interpolation, [position '.interpolation']);
        if ~ismember(interpolation, {'constant','previous','linear'})
            error('normalizeWaterBoundarySpec:badInterpolation', ...
                'water_%s_interpolation must be constant, previous, or linear.', position);
        end
    end
end

% -------------------------------------------------------------------------
function value = readRequiredFinite(s, field)
% Read one required finite numeric scalar.
    if ~isfield(s, field) || isempty(s.(field)) || ~isnumeric(s.(field)) ...
            || ~isscalar(s.(field)) || ~isfinite(s.(field))
        error('normalizeWaterBoundarySpec:badValue', ...
            '%s must be a finite numeric scalar.', field);
    end
    value = double(s.(field));
end

% -------------------------------------------------------------------------
function token = readToken(value, label)
% Normalize enum-like strings to lower-case char tokens.
    if ~(ischar(value) || (isstring(value) && isscalar(value)))
        error('normalizeWaterBoundarySpec:badToken', ...
            '%s must be a character vector or scalar string.', label);
    end
    token = lower(strtrim(char(value)));
    if isempty(token)
        error('normalizeWaterBoundarySpec:badToken', ...
            '%s must not be empty.', label);
    end
end

% -------------------------------------------------------------------------
function name = readName(value, label)
% Preserve case in external column names while rejecting empty/non-scalar names.
    if ~(ischar(value) || (isstring(value) && isscalar(value)))
        error('normalizeWaterBoundarySpec:badToken', ...
            '%s must be a character vector or scalar string.', label);
    end
    name = strtrim(char(value));
    if isempty(name)
        error('normalizeWaterBoundarySpec:badToken', ...
            '%s must not be empty.', label);
    end
end

% -------------------------------------------------------------------------
function requireFields(s, names, label)
% Fail once with the complete missing-field list.
    missing = names(~isfield(s, names));
    if ~isempty(missing)
        error('normalizeWaterBoundarySpec:missingField', ...
            '%s is missing field(s): %s.', label, strjoin(missing, ', '));
    end
end

% -------------------------------------------------------------------------
function checkUnknownFields(s, allowed, label)
% Strict typed schemas prevent misspellings from becoming silent no-ops.
    unknown = setdiff(fieldnames(s), allowed);
    if ~isempty(unknown)
        error('normalizeWaterBoundarySpec:unknownField', ...
            '%s has unsupported field(s): %s.', label, strjoin(unknown, ', '));
    end
end
