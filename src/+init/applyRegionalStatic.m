function built = applyRegionalStatic(built, regional)
%APPLYREGIONALSTATIC  Apply NetCDF regional overrides to an expanded built state.
%
%   built = init.applyRegionalStatic(built, regional) applies static-file
%   overrides produced by init.readRegionalStatic after
%   init.expandBuiltToRegional has created the regional container shapes. This
%   function is initialization-only: it validates and writes initial arrays, but
%   does not recompute physics or advance the time solver.
%
%   Inputs:
%     built     [inout] struct  regional built container with domain.nCell
%     regional  [in]    struct  .overrides from init.readRegionalStatic
%
%   Returns:
%     built with selected Y/parameters/forcing fields replaced by NetCDF values.
% Editor: Jizixing Zhu; Date: 2026-06-27.

    if nargin < 2 || ~isstruct(regional) || ~isfield(regional, 'overrides')
        return;
    end

    NN = built.domain.NN;
    NL = built.domain.NL;
    nD = built.domain.grid.nD;
    nCell = built.domain.nCell;
    [nAllCell, activeIndex] = getPackingInfo(built.domain, nCell);
    overrides = regional.overrides;

    if isfield(overrides, 'soil')
        built.parameters.soil = applyLayerOverrides( ...
            built.parameters.soil, overrides.soil, NL, nCell, nAllCell, activeIndex, 'soil');
    end
    if isfield(overrides, 'soil_n')
        built.parameters.soil_n = applyLayerOverrides( ...
            built.parameters.soil_n, overrides.soil_n, NL, nCell, nAllCell, activeIndex, 'soil_n');
    end
    if isfield(overrides, 'veg')
        built = applyVegOverrides(built, overrides.veg, NL, nCell, nAllCell, activeIndex);
    end
    if isfield(overrides, 'Y')
        built.Y = applyStateOverrides(built.Y, overrides.Y, NN, NL, nD, nCell, nAllCell, activeIndex);
    end
    if isfield(overrides, 'forcing')
        built.forcing = applyForcingOverrides(built.forcing, overrides.forcing, nCell, nAllCell, activeIndex);
    end

    built.domain.regional_static_path = regional.path;
end

% -------------------------------------------------------------------------
function target = applyLayerOverrides(target, source, NL, nCell, nAllCell, activeIndex, groupName)
% Apply [NL x nCell] layer-parameter overrides.
% Editor: Jizixing Zhu; Date: 2026-06-27.
    names = fieldnames(source);
    for i = 1:numel(names)
        name = names{i};
        value = source.(name);
        if isstruct(value)
            if ~isfield(target, name)
                target.(name) = struct();
            end
            target.(name) = applyLayerOverrides(target.(name), value, NL, nCell, ...
                nAllCell, activeIndex, [groupName '.' name]);
        else
            target.(name) = orientLayerCell(value, NL, nCell, nAllCell, activeIndex, [groupName '.' name]);
        end
    end
end

% -------------------------------------------------------------------------
function built = applyVegOverrides(built, veg, NL, nCell, nAllCell, activeIndex)
% Apply vegetation static overrides.
% Editor: Jizixing Zhu; Date: 2026-06-27.
    if isfield(veg, 'root') && isfield(veg.root, 'RAI')
        built.parameters.veg.root.RAI = orientLayerCell( ...
            veg.root.RAI, NL, nCell, nAllCell, activeIndex, 'veg.root.RAI');
    end
    if isfield(veg, 'structure') && isfield(veg.structure, 'PFT_type')
        pft = orientCellVector(veg.structure.PFT_type, nCell, nAllCell, activeIndex, 'veg.PFT_type');
        if numel(pft) ~= nCell
            error('applyRegionalStatic:badPFT', ...
                'veg.PFT_type must have nCell=%d entries.', nCell);
        end
        built.parameters.veg.structure.PFT_type = pft;
    end
end

% -------------------------------------------------------------------------
function Y = applyStateOverrides(Y, source, NN, NL, nD, nCell, nAllCell, activeIndex)
% Apply prognostic initial-state overrides.
% Editor: Jizixing Zhu; Date: 2026-06-27.
    names = fieldnames(source);
    for i = 1:numel(names)
        name = names{i};
        value = source.(name);
        if ismember(name, {'h','T','p_gas','nh4','no3','theta','bd','C_NH4_granule'})
            Y.(name) = orientNodeCell(value, NN, nCell, nAllCell, activeIndex, ['Y.' name]);
        elseif ismember(name, {'Theta_L','Theta_U','Theta_I','C_NH4_ice','C_NO3_ice'})
            Y.(name) = orientLayerState(value, NL, nD, nCell, nAllCell, activeIndex, ['Y.' name]);
        elseif strcmp(name, 'soil_gas') && isstruct(value)
            if isfield(value, 'C_O2')
                Y.soil_gas.C_O2 = orientNodeCell(value.C_O2, NN, nCell, nAllCell, activeIndex, ...
                    'Y.soil_gas.C_O2');
            end
        end
    end
end

% -------------------------------------------------------------------------
function forcing = applyForcingOverrides(forcing, source, nCell, nAllCell, activeIndex)
% Apply forcing time-series overrides.
% Editor: Jizixing Zhu; Date: 2026-06-27.
    names = fieldnames(source);
    for i = 1:numel(names)
        name = names{i};
        value = source.(name);
        if isstruct(value)
            if ~isfield(forcing, name)
                forcing.(name) = struct();
            end
            forcing.(name) = applyForcingOverrides(forcing.(name), value, nCell, nAllCell, activeIndex);
        else
            forcing.(name) = orientTimeCell(value, nCell, nAllCell, activeIndex, ['forcing.' name]);
        end
    end
end

% -------------------------------------------------------------------------
function A = orientLayerCell(A, NL, nCell, nAllCell, activeIndex, name)
% Orient a layer-by-cell field to [NL x nCell].
% Editor: Jizixing Zhu; Date: 2026-06-27.
    A = double(A);
    if isvector(A) && numel(A) == NL
        A = repmat(A(:), 1, nCell);
    elseif isequal(size(A), [NL, nCell])
        % already oriented
    elseif isequal(size(A), [nCell, NL])
        A = A.';
    elseif isequal(size(A), [NL, nAllCell])
        A = A(:, activeIndex);
    elseif isequal(size(A), [nAllCell, NL])
        A = A(activeIndex, :).';
    elseif ndims(A) == 3 && size(A, 1) == NL && prodDims(A, [2, 3]) == nAllCell
        A = reshape(A, NL, nAllCell);
        A = A(:, activeIndex);
    elseif ndims(A) == 3 && size(A, 3) == NL && prodDims(A, [1, 2]) == nAllCell
        A = reshape(A, nAllCell, NL);
        A = A(activeIndex, :).';
    else
        error('applyRegionalStatic:badLayerCell', ...
            '%s must be layer-by-cell or layer-by-grid.', name);
    end
end

% -------------------------------------------------------------------------
function A = orientNodeCell(A, NN, nCell, nAllCell, activeIndex, name)
% Orient a node-by-cell field to [NN x nCell].
% Editor: Jizixing Zhu; Date: 2026-06-27.
    A = double(A);
    if isvector(A) && numel(A) == NN
        A = repmat(A(:), 1, nCell);
    elseif isequal(size(A), [NN, nCell])
        % already oriented
    elseif isequal(size(A), [nCell, NN])
        A = A.';
    elseif isequal(size(A), [NN, nAllCell])
        A = A(:, activeIndex);
    elseif isequal(size(A), [nAllCell, NN])
        A = A(activeIndex, :).';
    elseif ndims(A) == 3 && size(A, 1) == NN && prodDims(A, [2, 3]) == nAllCell
        A = reshape(A, NN, nAllCell);
        A = A(:, activeIndex);
    elseif ndims(A) == 3 && size(A, 3) == NN && prodDims(A, [1, 2]) == nAllCell
        A = reshape(A, nAllCell, NN);
        A = A(activeIndex, :).';
    else
        error('applyRegionalStatic:badNodeCell', ...
            '%s must be node-by-cell or node-by-grid.', name);
    end
end

% -------------------------------------------------------------------------
function A = orientLayerState(A, NL, nD, nCell, nAllCell, activeIndex, name)
% Orient a layer-state field to [NL x nD x nCell].
% Editor: Jizixing Zhu; Date: 2026-06-27.
    A = double(A);
    if isequal(size(A), [NL, nD])
        A = repmat(A, 1, 1, nCell);
    elseif isequal(size(A), [NL, nD, nCell])
        % already oriented
    elseif isequal(size(A), [nCell, NL, nD])
        A = permute(A, [2, 3, 1]);
    elseif isequal(size(A), [NL, nD, nAllCell])
        A = A(:, :, activeIndex);
    elseif ndims(A) == 4 && size(A, 1) == NL && size(A, 2) == nD ...
            && prodDims(A, [3, 4]) == nAllCell
        A = reshape(A, NL, nD, nAllCell);
        A = A(:, :, activeIndex);
    else
        error('applyRegionalStatic:badLayerState', ...
            '%s must be layer-state-by-cell or layer-state-by-grid.', name);
    end
end

% -------------------------------------------------------------------------
function A = orientTimeCell(A, nCell, nAllCell, activeIndex, name)
% Orient a time-by-cell forcing field to [nTime x nCell].
% Editor: Jizixing Zhu; Date: 2026-06-27.
    A = double(A);
    if isvector(A)
        A = repmat(A(:), 1, nCell);
    elseif size(A, 2) == nCell
        % already oriented
    elseif size(A, 1) == nCell
        A = A.';
    elseif size(A, 2) == nAllCell
        A = A(:, activeIndex);
    elseif size(A, 1) == nAllCell
        A = A(activeIndex, :).';
    elseif ndims(A) == 3 && prodDims(A, [2, 3]) == nAllCell
        A = reshape(A, size(A, 1), nAllCell);
        A = A(:, activeIndex);
    elseif ndims(A) == 3 && prodDims(A, [1, 2]) == nAllCell
        nTime = size(A, 3);
        A = reshape(A, nAllCell, nTime);
        A = A(activeIndex, :).';
    else
        error('applyRegionalStatic:badTimeCell', ...
            '%s must be time-by-cell or time-by-grid.', name);
    end
end

% -------------------------------------------------------------------------
function x = orientCellVector(x, nCell, nAllCell, activeIndex, name)
% Orient a cell vector or grid to [nCell x 1].
% Editor: Jizixing Zhu; Date: 2026-06-27.
    if isscalar(x)
        x = repmat(x, nCell, 1);
    elseif isvector(x) && numel(x) == nCell
        x = x(:);
    elseif isvector(x) && numel(x) == nAllCell
        x = x(:);
        x = x(activeIndex);
    elseif numel(x) == nAllCell
        x = x(:);
        x = x(activeIndex);
    else
        error('applyRegionalStatic:badCellVector', ...
            '%s must be scalar, nCell vector, or full-grid vector.', name);
    end
end

% -------------------------------------------------------------------------
function [nAllCell, activeIndex] = getPackingInfo(domain, nCell)
% Return full-grid packing metadata when available.
% Editor: Jizixing Zhu; Date: 2026-06-27.
    nAllCell = nCell;
    activeIndex = (1:nCell)';
    if isfield(domain, 'space')
        if isfield(domain.space, 'nAllCell') && ~isempty(domain.space.nAllCell)
            nAllCell = domain.space.nAllCell;
        end
        if isfield(domain.space, 'active_index') && ~isempty(domain.space.active_index)
            activeIndex = domain.space.active_index(:);
        end
    end
end

% -------------------------------------------------------------------------
function n = prodDims(A, dims)
% Product of selected array dimensions without relying on vector-dim size().
% Editor: Jizixing Zhu; Date: 2026-06-27.
    sz = size(A);
    n = 1;
    for i = 1:numel(dims)
        dim = dims(i);
        if dim <= numel(sz)
            n = n * sz(dim);
        end
    end
end
