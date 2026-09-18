function ledger = buildSurfaceSoilEnergyLedger(interface_closure)
%BUILDSURFACESOILENERGYLEDGER Publish one accepted surface--soil heat ledger.
%
%   Inputs:
%     interface_closure [struct] Accepted closure audit produced by the soil
%                       scheduler for one outer interval.
%
%   Returns:
%     ledger [struct] Actual common flux and mismatch where the dynamic heat
%            interface applies; otherwise an explicit unavailable record.
%
%   The ledger is a diagnostic record of the flux actually passed to and
%   returned by the formal soil solve.  It must not be used to overwrite an
%   ebal diagnostic G or to reconstruct a boundary condition after acceptance.
%
%   Editor: Jizixing Zhu; Date: 2026-07-20.
    ledger = struct( ...
        'schema_version', 1, ...
        'sign_convention', 'flux_into_soil_positive', ...
        'available', false, ...
        'kind', 'not_available', ...
        'flux_into_soil_W_m2', NaN, ...
        'surface_flux_W_m2', NaN, ...
        'residual_W_m2', NaN, ...
        'interface_temperature_C', NaN);
    if ~isstruct(interface_closure) || isempty(fieldnames(interface_closure))
        return;
    end
    required = {'soil_flux_W_m2','surface_flux_W_m2', ...
        'final_residual_W_m2','final_temperature_C'};
    for i = 1:numel(required)
        if ~isfield(interface_closure, required{i}) ...
                || ~isscalar(interface_closure.(required{i})) ...
                || ~isfinite(interface_closure.(required{i}))
            return;
        end
    end
    ledger.available = true;
    ledger.kind = readReason(interface_closure);
    ledger.flux_into_soil_W_m2 = interface_closure.soil_flux_W_m2;
    ledger.surface_flux_W_m2 = interface_closure.surface_flux_W_m2;
    ledger.residual_W_m2 = interface_closure.final_residual_W_m2;
    ledger.interface_temperature_C = interface_closure.final_temperature_C;
end

function reason = readReason(closure)
% Read the diagnostic reason without permitting an empty dynamic type.
% Editor: Jizixing Zhu; Date: 2026-07-20.
    reason = 'legacy_or_disabled';
    if isfield(closure, 'reason') && (ischar(closure.reason) ...
            || isstring(closure.reason)) && strlength(string(closure.reason)) > 0
        reason = char(string(closure.reason));
    end
end
