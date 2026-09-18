function result = buildResult(request, outcome)
%BUILDRESULT  Build one water-heat-dry-air boundary result for a solver trial.
%
%   Requested values and solved fluxes remain separate. All public actual
%   fluxes are positive into the soil domain. This object is audit data only;
%   accepted accumulation and rollback remain scheduler responsibilities.
%
%   Editor: Jizixing Zhu; Date: 2026-07-15.
%   Change: 2026-07-16 Jizixing Zhu -- carry applied heat and dry-air
%           boundary kinds into the unified audit result.

    result = struct('schema_version', 1, ...
        'sign_convention', 'flux_into_domain_positive', ...
        't_start_s', request.t_start_s, 't_end_s', request.t_end_s, ...
        'dt_s', request.dt_s, 'accepted', logical(outcome.converged), ...
        'water', struct(), 'heat', struct(), 'dry_air', struct());

    if isfield(outcome, 'waterFlux') && ~isempty(outcome.waterFlux)
        applied = struct();
        if isfield(outcome, 'waterBoundaryApplied')
            applied = outcome.waterBoundaryApplied;
        end
        result.water = soil.water.boundary.buildResult( ...
            request.water, outcome.waterFlux, outcome.converged, applied);
    end
    if request.heat.enabled && isfield(outcome, 'heatFlux') ...
            && isstruct(outcome.heatFlux) ...
            && ~isempty(fieldnames(outcome.heatFlux))
        applied = struct();
        if isfield(outcome, 'heatBoundaryApplied')
            applied = outcome.heatBoundaryApplied;
        end
        result.heat = soil.thermal.boundary.buildResult( ...
            request.heat, outcome.heatFlux, outcome.converged, applied);
    end
    if request.dry_air.enabled && isfield(outcome, 'airFlux') ...
            && isstruct(outcome.airFlux) ...
            && ~isempty(fieldnames(outcome.airFlux))
        applied = struct();
        if isfield(outcome, 'airBoundaryApplied')
            applied = outcome.airBoundaryApplied;
        end
        result.dry_air = soil.gas.boundary.buildResult( ...
            request.dry_air, outcome.airFlux, outcome.converged, applied);
    end
end
