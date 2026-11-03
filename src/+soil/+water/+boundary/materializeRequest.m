function request = materializeRequest(spec, drivers, tStart_s, dt_s)
%MATERIALIZEREQUEST  Build one immutable water-boundary request for a trial.
%
%   request = soil.water.boundary.materializeRequest(spec, drivers,
%   tStart_s, dt_s) evaluates canonical top and bottom boundary definitions
%   over exactly one trial interval. The function is pure: it does not update
%   driver cursors, model clocks, prognostic state, or solver scratch.
%
%   Inputs:
%     spec      [in] canonical struct from normalizeSpec. Public fluxes are
%                    [cm/s], positive into the soil domain.
%     drivers   [in] plain-data struct. Optional fields:
%                    .natural_surface.liquid_supply_cm_s [cm/s, scalar]
%                    .natural_surface.evap_demand_cm_s   [cm/s, scalar]
%                    .boundary_series with time_s [s] and named value columns.
%                    boundary_series may be a scalar struct or table.
%     tStart_s  [in] [s] trial start relative to the simulation start.
%     dt_s      [in] [s] positive trial duration.
%
%   Returns:
%     request   plain-data struct with t_start_s/t_end_s/dt_s, sign convention,
%               and top/bottom endpoint requests. A variable-flux series is
%               returned as the exact interval-average piecewise-constant rate,
%               so rate*dt is invariant to retry/substep tiling. Constant modes
%               require an unchanged CSV column. A variable-head series is
%               evaluated at the trial end with linear interpolation.
%               Atmospheric supply and evaporation demand remain
%               separate; requested_flux_into_domain_cm_s is their net request.
%
%   Series convention:
%     value(i) applies on [time_s(i), time_s(i+1)) for fixed fluxes. The last
%     time is an explicit coverage endpoint. No extrapolation is allowed.
%
%   Architecture reference:
%     ClimaLand.jl src/shared_utilities/drivers.jl:44-83,175-187,1383-1391
%     and src/standalone/Soil/boundary_conditions.jl:164-211.
%
%   Editor: Jizixing Zhu; Date: 2026-07-11.
%   Change: 2026-07-11 Jizixing Zhu -- materialize constant head/flux from CSV and reject any change in a constant-mode column.

    validateTrialTime(tStart_s, dt_s);
    validateCanonicalSpec(spec);
    if nargin < 2 || isempty(drivers), drivers = struct(); end
    if ~isstruct(drivers) || ~isscalar(drivers)
        error('materializeWaterBoundaryRequest:badDrivers', ...
            'drivers must be a scalar struct.');
    end

    tEnd_s = tStart_s + dt_s;
    request = struct( ...
        'schema_version', 1, ...
        't_start_s', double(tStart_s), ...
        't_end_s', double(tEnd_s), ...
        'dt_s', double(dt_s), ...
        'sign_convention', 'flux_into_domain_positive', ...
        'top', materializeEndpoint(spec.top, 'top', drivers, tStart_s, tEnd_s), ...
        'bottom', materializeEndpoint(spec.bottom, 'bottom', drivers, tStart_s, tEnd_s));
end

% -------------------------------------------------------------------------
function endpoint = materializeEndpoint(policy, position, drivers, tStart_s, tEnd_s)
% Materialize one policy without changing the policy or driver data.
    head = NaN;
    requestedFlux = NaN;
    liquidSupply = NaN;
    evapDemand = NaN;

    switch policy.kind
        case 'fixed_head'
            if strcmp(policy.source, 'constant')
                head = policy.head_cm;
            elseif strcmp(policy.source, 'series')
                series = requireBoundarySeries(drivers);
                [time_s, values] = readSeries(series, policy.series_column);
                head = evaluateHead(time_s, values, tEnd_s, policy.interpolation);
            else
                unsupportedPolicy(policy, position);
            end
        case 'fixed_flux'
            if strcmp(policy.source, 'constant')
                requestedFlux = policy.flux_into_domain_cm_s;
            elseif strcmp(policy.source, 'series')
                series = requireBoundarySeries(drivers);
                [time_s, values] = readSeries(series, policy.series_column);
                if strcmp(policy.interpolation, 'constant')
                    requestedFlux = constantSeriesValue( ...
                        time_s, values, tStart_s, tEnd_s, position);
                else
                    requestedFlux = averagePreviousFlux( ...
                        time_s, values, tStart_s, tEnd_s);
                end
            else
                unsupportedPolicy(policy, position);
            end
        case 'zero_flux'
            requestedFlux = 0;
        case 'atmospheric'
            if ~strcmp(position, 'top') || ~strcmp(policy.source, 'natural_surface')
                unsupportedPolicy(policy, position);
            end
            [liquidSupply, evapDemand] = readNaturalSurface(drivers);
            requestedFlux = liquidSupply - evapDemand;
        case 'free_drainage'
            if ~strcmp(position, 'bottom') || ~strcmp(policy.source, 'none')
                unsupportedPolicy(policy, position);
            end
        otherwise
            unsupportedPolicy(policy, position);
    end

    endpoint = struct( ...
        'kind', policy.kind, ...
        'source', policy.source, ...
        'head_cm', head, ...
        'requested_flux_into_domain_cm_s', requestedFlux, ...
        'liquid_supply_cm_s', liquidSupply, ...
        'evap_demand_cm_s', evapDemand);
end

% -------------------------------------------------------------------------
function [liquidSupply, evapDemand] = readNaturalSurface(drivers)
% Natural input is a request from upstream components, not an observed h_sur.
    if ~isfield(drivers, 'natural_surface') ...
            || ~isstruct(drivers.natural_surface) || ~isscalar(drivers.natural_surface)
        error('materializeWaterBoundaryRequest:missingNaturalSurface', ...
            'Atmospheric top BC requires drivers.natural_surface.');
    end
    natural = drivers.natural_surface;
    liquidSupply = readFiniteScalar(natural, 'liquid_supply_cm_s', ...
        'materializeWaterBoundaryRequest:badNaturalSurface');
    evapDemand = readFiniteScalar(natural, 'evap_demand_cm_s', ...
        'materializeWaterBoundaryRequest:badNaturalSurface');
    if liquidSupply < 0 || evapDemand < 0
        error('materializeWaterBoundaryRequest:badNaturalSurface', ...
            'Natural liquid supply and evaporation demand must be nonnegative [cm/s].');
    end
end

% -------------------------------------------------------------------------
function series = requireBoundarySeries(drivers)
% Series data are external forcing, separate from the immutable policy.
    if ~isfield(drivers, 'boundary_series') || isempty(drivers.boundary_series)
        error('materializeWaterBoundaryRequest:missingSeries', ...
            'Series boundary source requires drivers.boundary_series.');
    end
    series = drivers.boundary_series;
    if ~(isstruct(series) && isscalar(series)) && ~istable(series)
        error('materializeWaterBoundaryRequest:badSeries', ...
            'drivers.boundary_series must be a scalar struct or table.');
    end
end

% -------------------------------------------------------------------------
function [time_s, values] = readSeries(series, column)
% Extract and validate one named series with its shared monotonic time axis.
    time_s = readSeriesColumn(series, 'time_s');
    values = readSeriesColumn(series, column);
    time_s = double(time_s(:));
    values = double(values(:));
    if numel(time_s) < 2 || numel(values) ~= numel(time_s) ...
            || any(~isfinite(time_s)) || any(~isfinite(values)) ...
            || any(diff(time_s) <= 0)
        error('materializeWaterBoundaryRequest:badSeries', ...
            ['Boundary series time_s and "%s" must be finite, equal-length ' ...
             'vectors with at least two rows and strictly increasing time.'], column);
    end
end

% -------------------------------------------------------------------------
function values = readSeriesColumn(series, name)
% Read a case-sensitive field/variable name from a struct or table.
    if istable(series)
        if ~ismember(name, series.Properties.VariableNames)
            error('materializeWaterBoundaryRequest:missingSeriesColumn', ...
                'Boundary series is missing column "%s".', name);
        end
        values = series.(name);
    else
        if ~isfield(series, name)
            error('materializeWaterBoundaryRequest:missingSeriesColumn', ...
                'Boundary series is missing field "%s".', name);
        end
        values = series.(name);
    end
    if ~isnumeric(values) || ~isvector(values)
        error('materializeWaterBoundaryRequest:badSeries', ...
            'Boundary series "%s" must be a numeric vector.', name);
    end
end

% -------------------------------------------------------------------------
function averageFlux = averagePreviousFlux(time_s, values, tStart_s, tEnd_s)
% Integrate piecewise-constant rates exactly over the requested interval.
    assertSeriesCoverage(time_s, tStart_s, tEnd_s);
    integral = 0;
    covered = 0;
    for i = 1:(numel(time_s) - 1)
        overlapStart = max(tStart_s, time_s(i));
        overlapEnd = min(tEnd_s, time_s(i + 1));
        overlap = max(0, overlapEnd - overlapStart);
        if overlap > 0
            integral = integral + values(i) * overlap;
            covered = covered + overlap;
        end
    end
    dt_s = tEnd_s - tStart_s;
    tol_s = 1e-9 * max(1, dt_s);
    if abs(covered - dt_s) > tol_s
        error('materializeWaterBoundaryRequest:seriesOutOfRange', ...
            'Flux series does not cover the full trial interval.');
    end
    averageFlux = integral / dt_s;
end

% -------------------------------------------------------------------------
function head = evaluateHead(time_s, values, tEnd_s, interpolation)
% Evaluate prescribed head at the implicit trial endpoint without extrapolation.
    assertSeriesCoverage(time_s, tEnd_s, tEnd_s);
    switch interpolation
        case 'constant'
            head = constantSeriesValue(time_s, values, tEnd_s, tEnd_s, 'head');
        case 'previous'
            idx = find(time_s <= tEnd_s, 1, 'last');
            head = values(idx);
        case 'linear'
            head = interp1(time_s, values, tEnd_s, 'linear');
        otherwise
            error('materializeWaterBoundaryRequest:badInterpolation', ...
                'Unsupported head interpolation "%s".', interpolation);
    end
    if ~isscalar(head) || ~isfinite(head)
        error('materializeWaterBoundaryRequest:badSeries', ...
            'Head series evaluation must return one finite value.');
    end
end

% -------------------------------------------------------------------------
function value = constantSeriesValue(time_s, values, tStart_s, tEnd_s, label)
% Require one unchanged CSV value while retaining the common time-axis contract.
    assertSeriesCoverage(time_s, tStart_s, tEnd_s);
    reference = values(1);
    tol = 8 * eps(max(1, abs(reference)));
    if any(abs(values - reference) > tol)
        error('materializeWaterBoundaryRequest:nonConstantSeries', ...
            ['The %s constant boundary requires every value in its CSV ' ...
             'column to equal %.15g.'], label, reference);
    end
    value = reference;
end

% -------------------------------------------------------------------------
function assertSeriesCoverage(time_s, tStart_s, tEnd_s)
% Require explicit coverage; never extrapolate a laboratory boundary protocol.
    tol_s = 1e-9 * max([1; abs(time_s(:)); abs(tStart_s); abs(tEnd_s)]);
    if tStart_s < time_s(1) - tol_s || tEnd_s > time_s(end) + tol_s
        error('materializeWaterBoundaryRequest:seriesOutOfRange', ...
            'Trial interval [%.15g, %.15g] s is outside series coverage [%.15g, %.15g] s.', ...
            tStart_s, tEnd_s, time_s(1), time_s(end));
    end
end

% -------------------------------------------------------------------------
function value = readFiniteScalar(s, field, errorId)
% Read one required finite scalar with a caller-specific contract error.
    if ~isfield(s, field) || ~isnumeric(s.(field)) || ~isscalar(s.(field)) ...
            || ~isfinite(s.(field))
        error(errorId, '%s must be a finite numeric scalar.', field);
    end
    value = double(s.(field));
end

% -------------------------------------------------------------------------
function validateTrialTime(tStart_s, dt_s)
% Trial time is explicit and independent of any mutable provider cursor.
    if ~isnumeric(tStart_s) || ~isscalar(tStart_s) || ~isfinite(tStart_s) ...
            || tStart_s < 0 || ~isnumeric(dt_s) || ~isscalar(dt_s) ...
            || ~isfinite(dt_s) || dt_s <= 0
        error('materializeWaterBoundaryRequest:badTime', ...
            'tStart_s must be finite and >=0; dt_s must be finite and >0.');
    end
end

% -------------------------------------------------------------------------
function validateCanonicalSpec(spec)
% Accept only the canonical version produced by normalizeSpec.
    if ~isstruct(spec) || ~isscalar(spec) ...
            || ~isfield(spec, 'schema_version') || spec.schema_version ~= 1 ...
            || ~isfield(spec, 'sign_convention') ...
            || ~strcmp(spec.sign_convention, 'flux_into_domain_positive') ...
            || ~isfield(spec, 'top') || ~isfield(spec, 'bottom')
        error('materializeWaterBoundaryRequest:badSpec', ...
            'spec must be canonical water-boundary schema version 1.');
    end
end

% -------------------------------------------------------------------------
function unsupportedPolicy(policy, position)
% Unsupported combinations fail instead of falling back to another BC.
    error('materializeWaterBoundaryRequest:unsupportedPolicy', ...
        'Unsupported %s policy kind="%s", source="%s".', ...
        position, policy.kind, policy.source);
end
