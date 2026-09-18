function candidate = materializeCandidate(post_ebal)
%MATERIALIZECANDIDATE Build the surface-owned accepted-state candidate.
%
%   Input:
%     post_ebal [struct] Cross-zone energy-stage result after the final
%               ordinary or bounded snow-correction pass.
%
%   Returns:
%     candidate [struct] Surface state: snow, tile thermal state, and the
%               force-restore history required by the next energy solve.
%
%   Physical contract:
%     The atmospheric-surface energy balance is solved by +coupling because
%     it consumes atmospheric forcing and vegetation resistance.  This
%     routine does not solve it again.  It selects only the resulting states
%     whose unique owner is the surface domain, so neither coupling nor
%     vegetation can silently install a surface state.
%
%   Editor: Jizixing Zhu; Date: 2026-07-20.
    required = {'snow', 'surface', 'scope_soil', 'surface_water'};
    for i = 1:numel(required)
        name = required{i};
        assert(isstruct(post_ebal) && isfield(post_ebal, name) ...
                && ~isempty(post_ebal.(name)), ...
            'surface:materializeCandidate:missingPostEbalField', ...
            'Final surface-energy result.%s is required.', name);
    end
    candidate = struct( ...
        'snow', post_ebal.snow, ...
        'surface', post_ebal.surface, ...
        'surface_water', post_ebal.surface_water, ...
        'force_restore', post_ebal.scope_soil);
end
