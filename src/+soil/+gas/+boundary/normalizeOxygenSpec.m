function spec = normalizeOxygenSpec(carrierSpec)
%NORMALIZEOXYGENSPEC Build the static O2 boundary policy.
%
%   The carrier pressure or mass-flux policy does not determine species
%   composition. This function records the current production O2 behavior as
%   an explicit compatibility policy so the nitrogen solver no longer infers
%   species boundaries from dry-air fields. A later physical change may
%   replace these endpoint kinds without changing the carrier contract.
%
%   Input:
%     carrierSpec [struct] normalized dry-air boundary policy.
%
%   Return:
%     spec [struct] O2 policy. Composition is a dry mole fraction [-]; applied
%       concentration and molar flux use [mol cm^-3 gas] and
%       [mol cm^-2 s^-1], respectively.
%
%   Reference: pre-refactor active behavior in
%   src/+soil/+nitrogen/solveNitrogen.m; architecture decision in
%   docs/soil/multigas_carrier_species_coupling_decision.md sections 11-12.
%
%   Editor: Jizixing Zhu; Date: 2026-08-20.
%   Change: 2026-08-20 Jizixing Zhu -- map a zero-free-gas carrier boundary
%   to dissolved O2 free-drainage outflow rather than zero total O2 flux.

    assert(isstruct(carrierSpec) && isfield(carrierSpec, 'bottom') ...
            && isfield(carrierSpec.bottom, 'kind'), ...
        'normalizeOxygenBoundarySpec:missingCarrier', ...
        'The O2 compatibility policy requires a normalized carrier boundary.');

    top = endpoint('fixed_atmospheric_dry_mole_fraction', ...
        'atmospheric_o2_forcing', false);
    carrierBottomKind = char(string(carrierSpec.bottom.kind));
    switch carrierBottomKind
        case 'hydrostatic_pressure'
            bottom = endpoint('fixed_atmospheric_dry_mole_fraction', ...
                'compatibility_mapping_from_carrier_pressure', true);
        case 'zero_mass_flux'
            bottom = endpoint('zero_total_flux', ...
                'compatibility_mapping_from_closed_carrier', false);
        case 'zero_gas_phase_flux'
            bottom = endpoint('zero_gas_phase_flux_with_liquid_outflow', ...
                'mapping_from_zero_free_gas_carrier', true);
        case 'specified_mass_flux'
            value = carrierSpec.bottom.configured_value;
            if value == 0
                bottom = endpoint('zero_total_flux', ...
                    'compatibility_mapping_from_zero_carrier_flux', false);
            else
                bottom = endpoint('external_composition_required', ...
                    'carrier_inflow_composition_not_configured', false);
            end
        otherwise
            error('normalizeOxygenBoundarySpec:unsupportedCarrierBottom', ...
                'Unsupported carrier bottom kind for O2: %s.', carrierBottomKind);
    end

    spec = struct('schema_version', 1, 'species', 'O2', ...
        'quantity', 'molar', 'provenance', ...
        'explicit_current_behavior_v1', 'top', top, 'bottom', bottom);
end

function out = endpoint(kind, source, candidateDependent)
% Allocate one fixed-shape species endpoint policy.
% Editor: Jizixing Zhu; Date: 2026-08-20.
    out = struct('kind', kind, 'composition_source', source, ...
        'candidate_state_dependent', logical(candidateDependent));
end

% Change: 2026-10-23 Wenjie Liu -- restore the COMPASS-COMPASS_MATLAB
% oxygen transport/source boundary implementation; retain current carrier physics.
