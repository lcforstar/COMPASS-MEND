function waterFlux = resolveCandidateBoundaryFlux( ...
        recoveredFlux, request, applied)
%RESOLVECANDIDATEBOUNDARYFLUX  Select prescribed fluxes or solved reactions.
%
%   waterFlux = soil.water.resolveCandidateBoundaryFlux(
%       recoveredFlux, request, applied)
%   preserves the coupled h/T/P_g/storage reaction at state-dependent
%   boundaries. At a prescribed Neumann boundary, however, the physical
%   flux is input data rather than an unknown reaction, so the exact request
%   replaces the recovered value. The subsequent column mass residual then
%   tests whether candidate storage is consistent with that prescribed flux.
%
%   Legacy signs are retained: QMT is positive outward at the top and QMB is
%   positive inward at the bottom, both [cm s^-1]. Atmospheric
%   supply_limited applies the exact net inward flux P_liq-E.
%
%   Editor: Jizixing Zhu; Date: 2026-08-09.

    assert(isstruct(recoveredFlux) ...
            && isfield(recoveredFlux, 'QMT') ...
            && isfield(recoveredFlux, 'QMB') ...
            && isscalar(recoveredFlux.QMT) ...
            && isscalar(recoveredFlux.QMB) ...
            && isfinite(recoveredFlux.QMT) ...
            && isfinite(recoveredFlux.QMB), ...
        'resolveCandidateBoundaryFlux:badRecoveredFlux', ...
        'recoveredFlux must contain finite scalar QMT and QMB values.');
    assert(isstruct(request) && isfield(request, 'top') ...
            && isfield(request, 'bottom') ...
            && isstruct(applied) && isfield(applied, 'top_kind') ...
            && isfield(applied, 'bottom_kind'), ...
        'resolveCandidateBoundaryFlux:badBoundaryContract', ...
        'request and applied boundary records are required.');

    waterFlux = recoveredFlux;
    topKind = char(applied.top_kind);
    switch topKind
        case 'fixed_flux'
            waterFlux.QMT = prescribedTopOutwardFlux(request.top);
        case 'zero_flux'
            waterFlux.QMT = 0;
    end

    bottomKind = char(applied.bottom_kind);
    switch bottomKind
        case 'fixed_flux'
            waterFlux.QMB = prescribedBottomInwardFlux(request.bottom);
        case 'zero_flux'
            waterFlux.QMB = 0;
    end
end

function qTopOut = prescribedTopOutwardFlux(top)
% Translate one prescribed top request to the legacy outward-positive sign.
    switch char(top.kind)
        case 'fixed_flux'
            requireFiniteField(top, 'requested_flux_into_domain_cm_s');
            qTopOut = -double(top.requested_flux_into_domain_cm_s);
        case 'atmospheric'
            requireFiniteField(top, 'liquid_supply_cm_s');
            requireFiniteField(top, 'evap_demand_cm_s');
            qTopOut = double(top.evap_demand_cm_s) ...
                - double(top.liquid_supply_cm_s);
        otherwise
            error('resolveCandidateBoundaryFlux:badTopFixedFlux', ...
                ['An applied fixed_flux top requires a fixed_flux or ', ...
                 'atmospheric request, not "%s".'], char(top.kind));
    end
end

function qBottomIn = prescribedBottomInwardFlux(bottom)
% Translate one prescribed bottom request to the legacy inward-positive sign.
    if ~strcmp(char(bottom.kind), 'fixed_flux')
        error('resolveCandidateBoundaryFlux:badBottomFixedFlux', ...
            ['An applied fixed_flux bottom requires a fixed_flux request, ', ...
             'not "%s".'], char(bottom.kind));
    end
    requireFiniteField(bottom, 'requested_flux_into_domain_cm_s');
    qBottomIn = double(bottom.requested_flux_into_domain_cm_s);
end

function requireFiniteField(record, name)
% Fail before a malformed request can enter a conservation ledger.
    assert(isfield(record, name) && isnumeric(record.(name)) ...
            && isscalar(record.(name)) && isfinite(record.(name)), ...
        'resolveCandidateBoundaryFlux:badPrescribedFlux', ...
        'Boundary field %s must be a finite scalar [cm s^-1].', name);
end
