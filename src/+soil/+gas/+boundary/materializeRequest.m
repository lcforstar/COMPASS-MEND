function request = materializeRequest(spec, drivers, tStart_s, dt_s)
%MATERIALIZEREQUEST  Freeze dry-air BC values for one solver trial.
%
%   Editor: Jizixing Zhu; Date: 2026-07-15.
%   Change: 2026-07-16 Jizixing Zhu -- record pressure provenance,
%           hydrostatic application, and mass-flux quantity semantics.

    request = struct('schema_version', 1, 'enabled', true, ...
        't_start_s', tStart_s, 't_end_s', tStart_s + dt_s, 'dt_s', dt_s, ...
        'top', materializeTop(spec.top, drivers), ...
        'bottom', materializeBottom(spec.bottom, drivers));
end

function out = materializeTop(spec, d)
%MATERIALIZETOP  Resolve one top dry-air policy for the trial.
% Editor: Jizixing Zhu; Date: 2026-07-15.
    out = emptyEndpoint(spec.kind);
    switch spec.kind
        case 'fixed_pressure'
            out.pressure_Pa = spec.configured_value;
            out.value_source = 'configured_pressure';
        case 'specified_mass_flux'
            out.requested_flux_into_domain_g_cm2_s = spec.configured_value;
            out.value_source = 'configured_mass_flux';
            out.flux_quantity = 'dry_air_mass';
        case 'atmospheric_pressure'
            out.pressure_Pa = d.atmospheric_pressure_Pa;
            out.value_source = 'atmospheric_pressure_forcing';
        otherwise
            error('materializeDryAirBoundaryRequest:badTopKind', ...
                'Unsupported dry-air top kind %s.', spec.kind);
    end
end

function out = materializeBottom(spec, d)
%MATERIALIZEBOTTOM  Resolve one bottom dry-air policy for the trial.
% Editor: Jizixing Zhu; Date: 2026-07-15.
    out = emptyEndpoint(spec.kind);
    switch spec.kind
        case 'hydrostatic_pressure'
            if spec.use_hydrostatic_atmosphere
                rhoAir_g_cm3 = 1.2e-3;
                pressureRisePa = rhoAir_g_cm3 * equations.getPhysicalConstant('g') ...
                    * d.column_depth_cm / 10;
                out.pressure_Pa = d.atmospheric_pressure_Pa + pressureRisePa;
                out.value_source = 'hydrostatic_atmospheric_pressure';
                out.hydrostatic_enabled = true;
            else
                out.pressure_Pa = d.atmospheric_pressure_Pa;
                out.value_source = 'atmospheric_pressure_no_hydrostatic';
            end
        case 'specified_mass_flux'
            out.requested_flux_into_domain_g_cm2_s = spec.configured_value;
            out.value_source = 'configured_mass_flux';
            out.flux_quantity = 'dry_air_mass';
        case 'zero_mass_flux'
            out.requested_flux_into_domain_g_cm2_s = 0;
            out.value_source = 'zero_mass_flux';
            out.flux_quantity = 'dry_air_mass';
        otherwise
            error('materializeDryAirBoundaryRequest:badBottomKind', ...
                'Unsupported dry-air bottom kind %s.', spec.kind);
    end
end

function out = emptyEndpoint(kind)
%EMPTYENDPOINT  Allocate the fixed-shape dry-air request endpoint.
% Editor: Jizixing Zhu; Date: 2026-07-15.
    out = struct('kind', kind, 'pressure_Pa', NaN, ...
        'requested_flux_into_domain_g_cm2_s', NaN, ...
        'value_source', '', 'hydrostatic_enabled', false, ...
        'flux_quantity', '');
end
