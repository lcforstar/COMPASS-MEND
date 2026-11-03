function result = buildResult(request, heatFlux, accepted, applied)
%BUILDRESULT  Package actual heat boundary fluxes with one public sign.
%
%   Legacy QET is outward-positive at the top and QEB is inward-positive at
%   the bottom. Public actual fluxes are positive into the soil at both ends.
%   Editor: Jizixing Zhu; Date: 2026-07-15.
%   Change: 2026-07-16 Jizixing Zhu -- distinguish requested heat policy,
%           applied mathematical BC, and solved reaction flux.
%   Change: 2026-07-17 Jizixing Zhu -- report the named free-drainage internal-
%           energy Robin condition without calling it a fixed flux.
%   Change: 2026-07-19 Jizixing Zhu -- retain affine surface-energy response
%           temperature and tangent metadata beside the solved reaction.

    if nargin < 4 || isempty(applied) ...
            || (isstruct(applied) && isempty(fieldnames(applied)))
        applied = inferApplied(request);
    end

    assert(isstruct(heatFlux) && all(isfield(heatFlux, {'QET','QEB'})), ...
        'buildHeatBoundaryResult:badFlux', 'heatFlux needs QET and QEB.');
    result = struct('schema_version', 1, ...
        'sign_convention', 'flux_into_domain_positive', ...
        'accepted', logical(accepted), ...
        'top', endpoint(request.top, applied.top_kind, -double(heatFlux.QET)), ...
        'bottom', endpoint(request.bottom, applied.bottom_kind, double(heatFlux.QEB)));
end

function out = endpoint(request, appliedKind, actual)
%ENDPOINT  Pair one requested heat boundary with its solved flux.
% Editor: Jizixing Zhu; Date: 2026-07-15.
% Change: 2026-07-16 Jizixing Zhu -- include the applied boundary kind.
% Change: 2026-07-16 Jizixing Zhu -- retain the frozen value source and
%         temperature-fallback audit fields.
% Change: 2026-07-17 Jizixing Zhu -- retain the distinct SCOPE skin
%         temperature and snow fraction behind a dynamic snow-insulated BC.
    out = struct('requested_kind', request.kind, ...
        'applied_kind', appliedKind, ...
        'requested_value_source', request.value_source, ...
        'request_fallback_used', logical(request.fallback_used), ...
        'requested_temperature_C', request.temperature_C, ...
        'surface_skin_temperature_C', requestField( ...
            request, 'surface_skin_temperature_C', NaN), ...
        'snow_cover_fraction', requestField( ...
            request, 'snow_cover_fraction', 0), ...
        'requested_flux_into_domain_J_cm2_s', ...
            request.requested_flux_into_domain_J_cm2_s, ...
        'reference_temperature_C', requestField( ...
            request, 'reference_temperature_C', NaN), ...
        'flux_temperature_derivative_J_cm2_s_K', requestField( ...
            request, 'flux_temperature_derivative_J_cm2_s_K', NaN), ...
        'actual_flux_into_domain_J_cm2_s', actual);
end

function value = requestField(request, name, fallback)
% Read one optional audit field from older or current typed requests.
% Editor: Jizixing Zhu; Date: 2026-07-17.
    value = fallback;
    if isfield(request, name) && ~isempty(request.(name))
        value = request.(name);
    end
end

function applied = inferApplied(request)
% Infer the deterministic applied kind for direct compatibility callers.
% Editor: Jizixing Zhu; Date: 2026-07-16.
    applied = struct( ...
        'top_kind', requestKindToApplied(request.top.kind), ...
        'bottom_kind', requestKindToApplied(request.bottom.kind));
end

function kind = requestKindToApplied(requestKind)
% Translate a physical heat policy to its Dirichlet/Neumann row type.
% Editor: Jizixing Zhu; Date: 2026-07-16.
    if ismember(requestKind, { ...
            'surface_temperature','seasonal_bottom_temperature','fixed_temperature'})
        kind = 'fixed_temperature';
    elseif strcmp(requestKind, 'free_drainage_internal_energy')
        kind = 'free_drainage_internal_energy';
    elseif strcmp(requestKind, 'surface_energy_response')
        kind = 'linearized_surface_energy_flux';
    else
        kind = 'fixed_total_energy_flux';
    end
end
