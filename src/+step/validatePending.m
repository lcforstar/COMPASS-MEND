function validatePending(pending)
%VALIDATEPENDING  Fail fast when an accepted step lacks required candidates.
%
%   Input:
%     pending [in] Destination-keyed transaction from step.emptyPending.
%
%   Required candidates are produced on every enabled or disabled physics
%   path. An empty required value therefore indicates a staging wiring gap,
%   not a physical zero. Optional gated candidates are not checked here.
%   Editor: Jizixing Zhu; Date: 2026-07-16.
%   Change: 2026-07-18 Jizixing Zhu -- require the finalized SCOPE soil-history
%           candidate on every accepted canopy/surface path.
%   Change: 2026-07-20 Jizixing Zhu -- validate typed inter-domain requests instead
%     of duplicated ownerless workspace flux fields.
    requireValue(pending.candidate.vegetation.luna_climate_entry, ...
        'candidate.vegetation.luna_climate_entry');
    requireValue(pending.candidate.surface.snow, 'candidate.surface.snow');
    requireValue(pending.candidate.surface.surface, ...
        'candidate.surface.surface');
    requireValue(pending.candidate.surface.surface_water, ...
        'candidate.surface.surface_water');
    requireValue(pending.candidate.surface.force_restore, ...
        'candidate.surface.force_restore');
    requireValue(pending.candidate.vegetation.canopy_water, ...
        'candidate.vegetation.canopy_water');
    requireValue(pending.candidate.vegetation.vegetation, ...
        'candidate.vegetation.vegetation');
    requireValue(pending.candidate.vegetation.leaf_resp_Ta, ...
        'candidate.vegetation.leaf_resp_Ta');

    requireValue(pending.request.surface_to_soil.water, ...
        'request.surface_to_soil.water');
    requireValue(pending.request.vegetation_to_soil.water, ...
        'request.vegetation_to_soil.water');
    requireValue(pending.request.vegetation_to_soil.nitrogen, ...
        'request.vegetation_to_soil.nitrogen');
    requireValue(pending.diagnostics.root_demand, ...
        'diagnostics.root_demand');
    requireValue(pending.management.fertilizer_applied, ...
        'management.fertilizer_applied');
end

function requireValue(value, path)
% Require one fixed-schema pending field to be populated.
    if isempty(value)
        error('runOneStep:pendingMissing', ...
            'step.commitPending: required pending.%s is empty (staging wiring gap).', ...
            path);
    end
end
