function built = refreshSoilAuxiliary(built, varargin)
%REFRESHSOILAUXILIARY  Compatibility wrapper for canonical soil-water refresh.
%
%   built = init.refreshSoilAuxiliary(built) delegates to
%   soil.water.updateWaterContentAuxiliary(built). New runtime code should call
%   soil.water.updateWaterContentAuxiliary directly; this wrapper is kept only
%   for existing initialization tests and transitional callers.
%
%   Optional name-value:
%     "apply_initial_freezing" [logical] default false. When true, first call
%     soil.water.initializeFreezingHead(built), then refresh the auxiliary cache.
%
%   Editor: Jizixing Zhu; Date: 2026-06-28.
%   Change: Reduced +init refresh to a thin compatibility wrapper around the
%   canonical +soil/+water calculation-layer functions.

    applyInitialFreezing = parseApplyInitialFreezing(varargin{:});
    if applyInitialFreezing
        built = soil.water.initializeFreezingHead(built);
    end
    built = soil.water.updateWaterContentAuxiliary(built);
end

function value = parseApplyInitialFreezing(varargin)
% Parse the one transitional option supported by this wrapper.
    value = false;
    if mod(numel(varargin), 2) ~= 0
        error('refreshSoilAuxiliary:badOptions', 'Options must be name-value pairs.');
    end
    for k = 1:2:numel(varargin)
        name = lower(char(varargin{k}));
        switch name
            case {'apply_initial_freezing','applyinitialfreezing'}
                value = logical(varargin{k + 1});
            otherwise
                error('refreshSoilAuxiliary:badOption', 'Unknown option "%s".', name);
        end
    end
end
