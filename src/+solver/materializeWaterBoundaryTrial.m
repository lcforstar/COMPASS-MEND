function request = materializeWaterBoundaryTrial(built, currentFlux, tStart_s, dt_s)
%MATERIALIZEWATERBOUNDARYTRIAL  Materialize one water BC request for a trial.
%
%   request = solver.materializeWaterBoundaryTrial(built, currentFlux,
%   tStart_s, dt_s) is the scheduler-side adapter between model data and the
%   pure boundary-policy materializer. It reads the canonical policy from
%   built.boundary.water, the current outer-step natural surface drivers from
%   currentFlux, and any prescribed laboratory series from
%   built.forcing.water_boundary. It does not mutate any input, provider
%   cursor, model clock, prognostic state, or solver scratch.
%
%   Natural-surface components remain separate: Throughfall is liquid supply
%   and EVAP is evaporation demand, both [cm/s]. The downstream boundary
%   package owns the public sign convention (flux positive into the soil).
%
%   Editor: Jizixing Zhu; Date: 2026-07-11.

    assert(isstruct(built) && isfield(built, 'boundary') ...
        && isfield(built.boundary, 'water'), ...
        'materializeWaterBoundaryTrial:missingSpec', ...
        'built.boundary.water must contain the canonical water-boundary policy.');
    assert(isstruct(currentFlux) && isfield(currentFlux, 'Throughfall') ...
        && isfield(currentFlux, 'EVAP'), ...
        'materializeWaterBoundaryTrial:missingNaturalFlux', ...
        'currentFlux must contain this-step Throughfall and EVAP.');

    drivers = struct();
    drivers.natural_surface = struct( ...
        'liquid_supply_cm_s', requireScalar(currentFlux.Throughfall, 'Throughfall'), ...
        'evap_demand_cm_s', requireScalar(currentFlux.EVAP, 'EVAP'));
    if isfield(built, 'forcing') && isfield(built.forcing, 'water_boundary') ...
            && ~isempty(built.forcing.water_boundary)
        drivers.boundary_series = built.forcing.water_boundary;
    end

    request = soil.water.boundary.materializeRequest( ...
        built.boundary.water, drivers, tStart_s, dt_s);
end

% -------------------------------------------------------------------------
function value = requireScalar(raw, name)
% Current water-only scheduler supports exactly one soil column.
    assert(isnumeric(raw) && isscalar(raw) && isfinite(raw), ...
        'materializeWaterBoundaryTrial:nonScalarFlux', ...
        '%s must be one finite scalar [cm/s] in the water-only scheduler.', name);
    value = double(raw);
end
