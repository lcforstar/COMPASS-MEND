function gridData = scatterRegionalColumns(columnData, domain, fillValue)
%SCATTERREGIONALCOLUMNS  Scatter land-only columns back to the original grid.
%
%   gridData = init.scatterRegionalColumns(columnData, domain, fillValue)
%   maps an array whose last dimension is the active land-cell column dimension
%   back to the original regional/global grid described by
%   domain.space.grid_shape and domain.space.active_index.
%
%   Inputs:
%     columnData  [numeric]  land-only data. Supported shapes include:
%                            [nCell x 1], [N x nCell],
%                            [N x M x nCell], [N x M x S x nCell].
%     domain      [struct]   built.domain with .nCell and .space.active_index
%     fillValue   [numeric]  value used for inactive/non-land cells (default NaN)
%
%   Returns:
%     gridData  [numeric]  same leading dimensions as columnData, but with the
%                          final land dimension replaced by domain.space.grid_shape.
%
%   Closure:
%     This is an initialization/output-shape utility only. It does not update
%     prognostic state, recompute parameters, or advance the solver.

    if nargin < 3
        fillValue = NaN;
    end
    if nargin < 2 || ~isstruct(domain) || ~isfield(domain, 'space')
        error('scatterRegionalColumns:badDomain', ...
            'domain.space with grid_shape and active_index is required.');
    end
    if ~isfield(domain.space, 'grid_shape') || ~isfield(domain.space, 'active_index')
        error('scatterRegionalColumns:missingPacking', ...
            'domain.space.grid_shape and domain.space.active_index are required.');
    end

    nCell = getDomainDim(domain, 'nCell', numel(domain.space.active_index));
    activeIndex = domain.space.active_index(:);
    gridShape = double(domain.space.grid_shape(:)).';
    nAllCell = prod(gridShape);
    if numel(activeIndex) ~= nCell
        error('scatterRegionalColumns:badActiveIndex', ...
            'active_index length (%d) must match nCell (%d).', numel(activeIndex), nCell);
    end
    if any(activeIndex < 1) || any(activeIndex > nAllCell)
        error('scatterRegionalColumns:activeIndexRange', ...
            'active_index entries must be within 1:nAllCell.');
    end

    [flatData, leadingShape] = orientColumns(columnData, nCell);
    flatGrid = fillLike(columnData, [size(flatData, 1), nAllCell], fillValue);
    flatGrid(:, activeIndex) = flatData;

    if isempty(leadingShape)
        gridData = reshape(flatGrid, gridShape);
    else
        gridData = reshape(flatGrid, [leadingShape, gridShape]);
    end
end

% -------------------------------------------------------------------------
function [flatData, leadingShape] = orientColumns(columnData, nCell)
% Move the land-cell dimension to the last dimension and flatten leading dims.
    sz = size(columnData);
    if isvector(columnData) && numel(columnData) == nCell
        flatData = reshape(columnData, 1, nCell);
        leadingShape = [];
    elseif sz(end) == nCell
        leadingShape = sz(1:end-1);
        flatData = reshape(columnData, prod(leadingShape), nCell);
    elseif sz(1) == nCell
        leadingShape = sz(2:end);
        perm = [2:ndims(columnData), 1];
        moved = permute(columnData, perm);
        flatData = reshape(moved, prod(leadingShape), nCell);
    else
        error('scatterRegionalColumns:badColumnShape', ...
            'columnData must have nCell as the first or last dimension.');
    end
end

% -------------------------------------------------------------------------
function A = fillLike(reference, sz, fillValue)
% Allocate a numeric fill array matching the reference class.
    A = zeros(sz, 'like', reference);
    A(:) = fillValue;
end

% -------------------------------------------------------------------------
function value = getDomainDim(domain, name, fallback)
% Read one domain dimension with fallback.
    value = fallback;
    if isfield(domain, name) && ~isempty(domain.(name))
        value = domain.(name);
    end
end
