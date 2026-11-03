function request = materializeRequest(spec, drivers, tStart_s, dt_s)
%MATERIALIZEREQUEST  Freeze heat BC values for one solver trial.
%
%   Editor: Jizixing Zhu; Date: 2026-07-15.
%   Change: 2026-07-15 Jizixing Zhu -- freeze SCOPE ground heat flux for
%           surface-energy-balance mode instead of rebuilding that balance
%           inside every Picard iteration.
%   Change: 2026-07-16 Jizixing Zhu -- record the exact value source and
%           surface-temperature fallback used by the immutable request.
%   Change: 2026-07-17 Jizixing Zhu -- keep SCOPE Tss as the exposed-soil
%           Dirichlet condition, but switch prognostic snow cover to the
%           conservative snow-insulated heat flux entering the soil.
%   Change: 2026-07-17 Jizixing Zhu -- require snow cold-content/latent energy to
%           exceed the accepted energy tolerance before changing BC type.
%   Change: 2026-07-17 Jizixing Zhu -- materialize the ClimaLand-style bottom
%           free-drainage internal-energy policy without freezing its flux.
%   Change: 2026-07-18 Jizixing Zhu -- after the explicit two-tile surface
%           became available, use the aggregate soil-interface energy flux
%           for every resolved partial/full snow cover while preserving the
%           original bare-soil Tss boundary exactly when snow is absent.
%   Change: 2026-07-19 Jizixing Zhu -- consume one affine surface-energy
%           response for bare and snow-covered ground; retain the old
%           temperature/fixed-flux branches only for direct legacy callers.
%   Change: 2026-07-19 Jizixing Zhu -- require that response explicitly for
%           the production surface_energy_response configuration.
%   Change: 2026-07-19 Jizixing Zhu -- preserve explicit temperature and
%           frozen-flux policies when an optional response is also present.
%   Change: 2026-07-20 Jizixing Zhu -- rename the natural matched-flux policy
%           atmospheric to align heat and water public configuration.

    assert(isstruct(spec) && isfield(spec, 'top') && isfield(spec, 'bottom'), ...
        'materializeHeatBoundaryRequest:badSpec', 'Heat spec needs top and bottom endpoints.');
    request = struct('schema_version', 1, 'enabled', true, ...
        't_start_s', tStart_s, 't_end_s', tStart_s + dt_s, 'dt_s', dt_s, ...
        'top', materializeTop(spec.top, drivers), ...
        'bottom', materializeBottom(spec.bottom, drivers));
end

function out = materializeTop(spec, d)
%MATERIALIZETOP  Resolve one top heat policy against frozen trial drivers.
% Editor: Jizixing Zhu; Date: 2026-07-15.
    out = emptyEndpoint(spec.kind);
    out.surface_skin_temperature_C = scalarField( ...
        d, 'surface_temperature_C', NaN);
    out.snow_cover_fraction = scalarField(d, 'snow_cover_fraction', 0);
    out.snow_ice_water_equivalent_cm = scalarField( ...
        d, 'snow_ice_water_equivalent_cm', 0);
    out.snow_thermal_energy_J_cm2 = scalarField( ...
        d, 'snow_thermal_energy_J_cm2', NaN);
    out.snow_thermal_tolerance_J_cm2 = scalarField( ...
        d, 'snow_thermal_tolerance_J_cm2', 0);
    out.snow_thermal_active = hasResolvedSnowThermalState(out);
    switch spec.kind
        case 'surface_energy_response'
            assert(hasTopHeatResponse(d), ...
                'materializeHeatBoundaryRequest:missingSurfaceResponse', ...
                ['surface_energy_response requires the typed affine response ', ...
                 'assembled by the surface domain.']);
            out = attachTopHeatResponse(out, d.top_heat_response);
        case 'surface_temperature'
            value = out.surface_skin_temperature_C;
            if isfinite(value)
                out.value_source = 'scope_surface_temperature';
            else
                value = d.air_temperature_C;
                out.value_source = 'air_temperature_fallback';
                out.fallback_used = true;
            end
            out.temperature_C = value;
        case 'specified_heat_flux'
            out.requested_flux_into_domain_J_cm2_s = spec.configured_value;
            out.value_source = 'configured_total_energy_flux';
        case 'atmospheric'
            out.temperature_C = out.surface_skin_temperature_C;
            out.requested_flux_into_domain_J_cm2_s = ...
                requiredGroundHeatFlux(d) * 1e-4;
            out.value_source = 'matched_atmospheric_soil_heat_flux';
        otherwise
            error('materializeHeatBoundaryRequest:badTopKind', ...
                'Unsupported heat top kind %s.', spec.kind);
    end
end

function out = materializeBottom(spec, d)
%MATERIALIZEBOTTOM  Resolve one bottom heat policy for the trial.
% Editor: Jizixing Zhu; Date: 2026-07-15.
    out = emptyEndpoint(spec.kind);
    switch spec.kind
        case 'seasonal_bottom_temperature'
            out.temperature_C = d.bottom_temperature_C;
            out.value_source = 'seasonal_bottom_temperature';
        case 'fixed_temperature'
            out.temperature_C = spec.configured_value;
            out.value_source = 'configured_temperature';
        case 'specified_heat_flux'
            out.requested_flux_into_domain_J_cm2_s = spec.configured_value;
            out.value_source = 'configured_total_energy_flux';
        case 'free_drainage_internal_energy'
            % The flux depends on the current Picard QMB and bottom T. It is
            % assembled later with the water and heat rows, not frozen here.
            out.value_source = 'coupled_bottom_water_flux_and_temperature';
        otherwise
            error('materializeHeatBoundaryRequest:badBottomKind', ...
                'Unsupported heat bottom kind %s.', spec.kind);
    end
end

function out = emptyEndpoint(kind)
%EMPTYENDPOINT  Allocate the fixed-shape heat request endpoint.
% Editor: Jizixing Zhu; Date: 2026-07-15.
    out = struct('kind', kind, 'temperature_C', NaN, ...
        'requested_flux_into_domain_J_cm2_s', NaN, ...
        'reference_temperature_C', NaN, ...
        'flux_temperature_derivative_J_cm2_s_K', NaN, ...
        'surface_skin_temperature_C', NaN, ...
        'snow_cover_fraction', 0, ...
        'snow_ice_water_equivalent_cm', 0, ...
        'snow_thermal_energy_J_cm2', NaN, ...
        'snow_thermal_tolerance_J_cm2', 0, ...
        'snow_thermal_active', false, ...
        'soil_evaporation_demand_cm_s', NaN, ...
        'soil_evaporation_latent_heat_flux_W_m2', NaN, ...
        'soil_evaporation_latent_heat_flux_temperature_derivative_W_m2_K', NaN, ...
        'value_source', '', 'fallback_used', false);
end

function tf = hasTopHeatResponse(drivers)
% Detect the typed production response without guessing from snow state.
% Editor: Jizixing Zhu; Date: 2026-07-19.
    tf = isfield(drivers, 'top_heat_response') ...
        && isstruct(drivers.top_heat_response) ...
        && ~isempty(fieldnames(drivers.top_heat_response));
end

function out = attachTopHeatResponse(out, response)
% Convert W m-2 response fields to the heat kernel's J cm-2 s-1 units.
% Editor: Jizixing Zhu; Date: 2026-07-19.
    required = {'reference_temperature_C', ...
        'reference_flux_into_soil_W_m2', ...
        'dflux_dtemperature_W_m2_K'};
    for i = 1:numel(required)
        name = required{i};
        assert(isfield(response, name) && isscalar(response.(name)) ...
                && isfinite(response.(name)), ...
            'materializeHeatBoundaryRequest:badSurfaceResponse', ...
            'top_heat_response.%s must be a finite scalar.', name);
    end
    assert(response.dflux_dtemperature_W_m2_K <= 1e-10, ...
        'materializeHeatBoundaryRequest:badSurfaceResponseDerivative', ...
        'Surface inward energy flux must not increase with temperature.');
    out.kind = 'surface_energy_response';
    out.reference_temperature_C = response.reference_temperature_C;
    out.temperature_C = response.reference_temperature_C;
    out.requested_flux_into_domain_J_cm2_s = ...
        response.reference_flux_into_soil_W_m2 * 1e-4;
    out.flux_temperature_derivative_J_cm2_s_K = ...
        response.dflux_dtemperature_W_m2_K * 1e-4;
    out.soil_evaporation_demand_cm_s = optionalResponseScalar( ...
        response, 'soil_evaporation_demand_cm_s', NaN);
    out.soil_evaporation_latent_heat_flux_W_m2 = optionalResponseScalar( ...
        response, 'soil_evaporation_latent_heat_flux_W_m2', NaN);
    out.soil_evaporation_latent_heat_flux_temperature_derivative_W_m2_K = ...
        optionalResponseScalar(response, ...
        'soil_evaporation_latent_heat_flux_temperature_derivative_W_m2_K', NaN);
    out.value_source = 'coupled_surface_energy_temperature_response';
end

function value = optionalResponseScalar(response, name, fallback)
% Read optional response metadata used by the coupled dry-surface branch.
% Editor: Jizixing Zhu; Date: 2026-07-19.
    value = fallback;
    if isfield(response, name) && isscalar(response.(name)) ...
            && (isfinite(response.(name)) || isnan(response.(name)))
        value = double(response.(name));
    end
end

function active = hasResolvedSnowThermalState(endpoint)
% Switch the shared column for any physically resolved snow-covered fraction.
%
% The explicit bare/snow surface tiles provide the area-integrated energy
% flux for partial cover. With no resolved snow the exact legacy bare-soil
% Dirichlet Tss path is retained; snow therefore cannot perturb snow-free
% hours merely because the snow process is enabled.
% Editor: Jizixing Zhu; Date: 2026-07-17.
% Change: 2026-07-18 Jizixing Zhu -- activate partial snow now that the
%   two-tile aggregate is explicit, while retaining a strict zero-snow limit.
    coverTolerance = 64 * eps(max(1, abs(endpoint.snow_cover_fraction)));
    iceTolerance = 64 * eps(max(1, ...
        abs(endpoint.snow_ice_water_equivalent_cm)));
    active = endpoint.snow_cover_fraction > coverTolerance ...
        && endpoint.snow_ice_water_equivalent_cm > iceTolerance;
    if isfinite(endpoint.snow_thermal_energy_J_cm2)
        active = active && abs(endpoint.snow_thermal_energy_J_cm2) ...
            > endpoint.snow_thermal_tolerance_J_cm2;
    else
        active = false;
    end
end

function value = requiredGroundHeatFlux(d)
% Read the finite SCOPE/snow flux that actually reaches the soil surface.
% Editor: Jizixing Zhu; Date: 2026-07-17.
    assert(isfield(d, 'ground_heat_flux_into_soil_W_m2') ...
            && isscalar(d.ground_heat_flux_into_soil_W_m2) ...
            && isfinite(d.ground_heat_flux_into_soil_W_m2), ...
        'materializeHeatBoundaryRequest:missingGroundHeatFlux', ...
        ['A flux-driven heat boundary requires the finite SCOPE/snow ', ...
         'ground heat flux entering soil for this outer step.']);
    value = double(d.ground_heat_flux_into_soil_W_m2);
end

function value = scalarField(s, name, fallback)
% Read one optional finite scalar driver, retaining NaN only for skin T.
% Editor: Jizixing Zhu; Date: 2026-07-17.
    value = fallback;
    if isfield(s, name) && ~isempty(s.(name)), value = double(s.(name)); end
    assert(isscalar(value) && (isfinite(value) || isnan(value)), ...
        'materializeHeatBoundaryRequest:badDriver', ...
        'Heat driver %s must be scalar.', name);
    if strcmp(name, 'snow_cover_fraction')
        assert(isfinite(value) && value >= 0 && value <= 1, ...
            'materializeHeatBoundaryRequest:badSnowCover', ...
            'snow_cover_fraction must lie in [0,1].');
    elseif strcmp(name, 'snow_ice_water_equivalent_cm')
        assert(isfinite(value) && value >= 0, ...
            'materializeHeatBoundaryRequest:badSnowIce', ...
            'snow_ice_water_equivalent_cm must be finite and nonnegative.');
    elseif strcmp(name, 'snow_thermal_tolerance_J_cm2')
        assert(isfinite(value) && value >= 0, ...
            'materializeHeatBoundaryRequest:badSnowThermalTolerance', ...
            'snow_thermal_tolerance_J_cm2 must be finite and nonnegative.');
    end
end
