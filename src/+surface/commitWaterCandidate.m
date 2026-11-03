function built = commitWaterCandidate(built, candidate)
%COMMITWATERCANDIDATE Install an accepted standalone surface-water result.
%
%   The production multi-domain path commits surface water through
%   surface.commitCandidate. This narrow installer exists for the standalone
%   water experiment, which has no snow/vegetation transaction.
%
%   Editor: Jizixing Zhu; Date: 2026-07-20.
    required = {'DSTOR','DSTOR0','RS'};
    for i = 1:numel(required)
        name = required{i};
        assert(isstruct(candidate) && isfield(candidate, name) ...
                && isscalar(candidate.(name)) && isfinite(candidate.(name)), ...
            'surface:commitWaterCandidate:invalidCandidate', ...
            'Accepted surface-water candidate.%s must be finite.', name);
    end
    built.Y.surface_water = candidate;
end
