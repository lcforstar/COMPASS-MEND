function [x, meta] = packSoilPicardState(SoilVariables, opts)
%PACKSOILPICARDSTATE Pack primary soil Picard unknowns into one vector.
%
%   [x, meta] = solver.packSoilPicardState(SoilVariables, opts)
%
%   Only primary unknowns are packed. Derived fields such as Theta_LL, KfL_h,
%   RHOV, Se, and hh_frez must be recomputed by the normal physics update
%   after the state is unpacked.
%
%   opts fields:
%       includePressure - include P_gg in the state, default false
%       includeThermal  - include TT or UU, default true
%       includeEnergy   - pack UU instead of TT, default false
%       scale_hh        - scale for hh, default 100
%       scale_TT        - scale for TT, default 10
%       scale_UU        - scalar or nodal U scale; generic fallback is 1.
%                         The production h/U/P_g caller supplies accepted
%                         nodal U_T times a 10 degC equivalent increment.
%       scale_P_gg      - scale for P_gg, default 1e5
%
%   Editor: Jizixing Zhu; Date: 2026-07-14.
%   Change: 2026-07-14 Jizixing Zhu -- support prognostic U packing with a
%           frozen scalar or nodal energy scale for the h/U/P_g path.
%   Change: 2026-07-14 Jizixing Zhu -- omit the inactive thermal coordinate
%           from water and water-dry-air Anderson vectors.
%   Change: 2026-07-15 Jizixing Zhu -- distinguish the generic scale_UU
%           fallback from the capacity-based production scaling contract.

if nargin < 2 || isempty(opts)
    opts = struct();
end

includePressure = getOption(opts, 'includePressure', false);
includeThermal = getOption(opts, 'includeThermal', true);
includeEnergy = getOption(opts, 'includeEnergy', false);
scale_hh = getOption(opts, 'scale_hh', 100);
scale_TT = getOption(opts, 'scale_TT', 10);
scale_UU = getOption(opts, 'scale_UU', 1);
scale_P_gg = getOption(opts, 'scale_P_gg', 1e5);

hh = SoilVariables.hh(:);
n = numel(hh);
validateScale(scale_hh, n, 'scale_hh');

assert(~includeEnergy || includeThermal, ...
    'solver:packSoilPicardState:InactiveEnergy', ...
    'includeEnergy requires includeThermal=true.');

xParts = {hh ./ expandScale(scale_hh, n)};
names = {'hh'};
lengths = n;
scaleVectors = {expandScale(scale_hh, n)};
scales = scaleVectors{1}(1);

if includeThermal
    if includeEnergy
        if ~isfield(SoilVariables, 'UU')
            error('solver:packSoilPicardState:MissingEnergy', ...
                'SoilVariables.UU is required when includeEnergy is true.');
        end
        thermalState = SoilVariables.UU(:);
        thermalName = 'UU';
        thermalScale = scale_UU;
    else
        thermalState = SoilVariables.TT(:);
        thermalName = 'TT';
        thermalScale = scale_TT;
    end
    if numel(thermalState) ~= n
        error('solver:packSoilPicardState:SizeMismatch', ...
            'The hydraulic and thermal primary states must have the same length.');
    end
    validateScale(thermalScale, n, ['scale_' thermalName]);
    thermalScaleVector = expandScale(thermalScale, n);
    xParts{end+1} = thermalState ./ thermalScaleVector;
    names{end+1} = thermalName;
    lengths(end+1) = n;
    scaleVectors{end+1} = thermalScaleVector;
    scales(end+1) = thermalScaleVector(1);
end

if includePressure
    if ~isfield(SoilVariables, 'P_gg')
        error('solver:packSoilPicardState:MissingPressure', ...
            'SoilVariables.P_gg is required when includePressure is true.');
    end
    P_gg = SoilVariables.P_gg(:);
    if numel(P_gg) ~= n
        error('solver:packSoilPicardState:SizeMismatch', ...
            'SoilVariables.P_gg must have the same length as hh.');
    end
    validateScale(scale_P_gg, n, 'scale_P_gg');
    pressureScale = expandScale(scale_P_gg, n);
    xParts{end+1} = P_gg ./ pressureScale;
    names{end+1} = 'P_gg';
    lengths(end+1) = n;
    scales(end+1) = pressureScale(1);
    scaleVectors{end+1} = pressureScale;
end

x = vertcat(xParts{:});

meta = struct();
meta.names = names;
meta.lengths = lengths;
meta.scales = scales;
meta.n = n;
meta.includePressure = includePressure;
meta.includeThermal = includeThermal;
meta.includeEnergy = includeEnergy;
meta.scaleVectors = scaleVectors;
end

function value = getOption(opts, name, defaultValue)
% Read one optional pack control without mutating the caller's struct.
% Editor: Jizixing Zhu; Date: 2026-07-14.
if isfield(opts, name) && ~isempty(opts.(name))
    value = opts.(name);
else
    value = defaultValue;
end
end

function validateScale(scale, n, name)
% Require a positive scalar or one positive value per packed node.
% Editor: Jizixing Zhu; Date: 2026-07-14.
    ok = isvector(scale) && (isscalar(scale) || numel(scale) == n) ...
        && all(isfinite(scale(:))) && all(scale(:) > 0);
    assert(ok, 'solver:packSoilPicardState:BadScale', ...
        '%s must be a finite positive scalar or n-vector.', name);
end

function expanded = expandScale(scale, n)
% Expand scalar scaling while preserving a nodal scale vector.
% Editor: Jizixing Zhu; Date: 2026-07-14.
    if isscalar(scale)
        expanded = repmat(scale, n, 1);
    else
        expanded = scale(:);
    end
end
