function built = commitCandidate(built, candidate)
%COMMITCANDIDATE Install the accepted surface-domain state candidate.
%
%   Inputs:
%     built     [struct] Legacy state shell before outer-step acceptance.
%     candidate [struct] Surface-owned snow, surface, and force-restore state.
%
%   Returns:
%     built [struct] Shell with the accepted surface state installed.
%
%   Physical contract:
%     Snow S/S_l/U, surface thermal state, and force-restore history are one
%     surface-domain state package.  They are installed only after the coupled
%     soil endpoint has accepted; this routine computes neither a flux nor a
%     snow transition.
%
%   Editor: Jizixing Zhu; Date: 2026-07-20.
%   Change: 2026-07-20 Jizixing Zhu -- commit depression storage/runoff through the
%     surface candidate instead of inside the soil scheduler.
    names = {'snow','surface','surface_water','force_restore'};
    for i = 1:numel(names)
        assert(isstruct(candidate) && isfield(candidate, names{i}) ...
                && ~isempty(candidate.(names{i})), ...
            'surface:commitCandidate:missingCandidate', ...
            'Surface candidate.%s is required at accepted commit.', names{i});
    end
    built.Y.snow = candidate.snow;
    built.Y.surface = candidate.surface;
    built.Y.surface_water = candidate.surface_water;
    built.Y.scope_soil = candidate.force_restore;
end
