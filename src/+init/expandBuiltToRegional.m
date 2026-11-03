function built = expandBuiltToRegional(built, space)
%EXPANDBUILTTOREGIONAL  Broadcast a column built state onto horizontal cells.
%
%   built = init.expandBuiltToRegional(built, space) converts an init.buildState
%   column template into a regional/global initialization container. No physics
%   is changed: values are copied to each cell so the result is a valid
%   multi-cell initial condition. Future regional NetCDF loaders should replace
%   these broadcasted arrays with spatial fields before the solver runs.
%
%   Inputs:
%     built  [inout] struct  init.buildState output before allocateWorkspace
%     space  [struct]        lat/lon/mask/cell_id descriptor (see
%                            init.initializeRegional)
%
%   Returns built with:
%     domain.kind/nCell/space set, Y node fields [NN x nCell], layer fields
%     [NL x nD x nCell], soil fields [NL x nCell], forcing fields
%     [nTime x nCell], and no stale workspace/diagnostics.
% Editor: Jizixing Zhu; Date: 2026-06-27.

    space = normalizeSpace(space);
    nCell = numel(space.cell_id);
    NN = built.domain.NN;
    NL = built.domain.NL;
    nD = built.domain.grid.nD;
    nS = built.domain.grid.nS;

    built.domain.kind = 'regional';
    built.domain.nCell = nCell;
    built.domain.space = space;
    built.domain.reference_site = built.domain.site;

    built.Y = broadcastState(built.Y, NN, NL, nD, nCell);
    if isfield(built.parameters, 'soil')
        built.parameters.soil = broadcastLayerStruct(built.parameters.soil, NN, NL, nD, nS, nCell);
    end
    if isfield(built.parameters, 'soil_n')
        built.parameters.soil_n = broadcastLayerStruct(built.parameters.soil_n, NN, NL, nD, nS, nCell);
    end
    if isfield(built.parameters, 'veg') && isfield(built.parameters.veg, 'root') ...
            && isfield(built.parameters.veg.root, 'RAI')
        built.parameters.veg.root.RAI = repmat(built.parameters.veg.root.RAI(:), 1, nCell);
    end
    if isfield(built.parameters, 'veg') && isfield(built.parameters.veg, 'structure') ...
            && isfield(built.parameters.veg.structure, 'PFT_type')
        built.parameters.veg.structure.PFT_type = repmat(built.parameters.veg.structure.PFT_type, nCell, 1);
    end
    built.forcing = broadcastForcing(built.forcing, nCell);

    if isfield(built, 'workspace')
        built = rmfield(built, 'workspace');
    end
    if isfield(built, 'diagnostics')
        built = rmfield(built, 'diagnostics');
    end
end

% ---------------------------------------------------------------------------
function space = normalizeSpace(space)
% Validate and normalize the horizontal space descriptor.
% Editor: Jizixing Zhu; Date: 2026-06-27.
    if ~isfield(space, 'lat') || ~isfield(space, 'lon')
        error('expandBuiltToRegional:missingLatLon', 'space.lat and space.lon are required.');
    end
    lat = space.lat(:);
    lon = space.lon(:);
    if numel(lat) ~= numel(lon)
        error('expandBuiltToRegional:latLonSize', 'space.lat and space.lon must have the same length.');
    end
    nCell = numel(lat);
    if nCell < 1
        error('expandBuiltToRegional:emptySpace', 'space must contain at least one cell.');
    end
    if isfield(space, 'mask')
        mask = logical(space.mask(:));
    else
        mask = true(nCell, 1);
    end
    if numel(mask) ~= nCell
        error('expandBuiltToRegional:maskSize', 'space.mask must have nCell entries.');
    end
    if isfield(space, 'cell_id')
        cellId = space.cell_id(:);
    else
        cellId = (1:nCell)';
    end
    if numel(cellId) ~= nCell
        error('expandBuiltToRegional:cellIdSize', 'space.cell_id must have nCell entries.');
    end
    normalized = struct('lat', lat, 'lon', lon, 'mask', mask, 'cell_id', cellId);
    optionalFields = {'grid_shape','nAllCell','active_index','full_active_index','pack','is_subset'};
    for i = 1:numel(optionalFields)
        name = optionalFields{i};
        if isfield(space, name)
            normalized.(name) = space.(name);
        end
    end
    space = normalized;
end

% ---------------------------------------------------------------------------
function Y = broadcastState(Y, NN, NL, nD, nCell)
% Broadcast prognostic state fields to regional shapes.
% Editor: Jizixing Zhu; Date: 2026-06-27.
% Change: 2026-07-05 Jizixing Zhu -- broadcast the state-ownership migration scalar containers snow/surface/canopy_water to [1 x nCell] (increment 3; also closes the latent increments 1-2 regional gap).
% Change: 2026-07-14 Jizixing Zhu -- broadcast prognostic U to regional [NN x nCell].
    nodeFields = {'h','T','p_gas','U','nh4','no3','urea','theta','bd','C_NH4_granule'};
    for i = 1:numel(nodeFields)
        name = nodeFields{i};
        if isfield(Y, name)
            Y.(name) = repmat(Y.(name)(:), 1, nCell);
        end
    end

    elemFields = {'Theta_L','Theta_U','Theta_I','C_NH4_ice','C_NO3_ice'};
    for i = 1:numel(elemFields)
        name = elemFields{i};
        if isfield(Y, name)
            Y.(name) = reshapeLayerField(Y.(name), NL, nD, nCell, name);
        end
    end

    if isfield(Y, 'soil_gas') && isfield(Y.soil_gas, 'C_O2')
        Y.soil_gas.C_O2 = repmat(Y.soil_gas.C_O2(:), 1, nCell);
    end

    % canopy-air OLD profile: per-layer [nlayers x 1] -> [nlayers x nCell]
    if isfield(Y, 'canopy_air')
        canopyLayerFields = {'Ta_layer','ea_layer','CO2_layer'};
        for i = 1:numel(canopyLayerFields)
            name = canopyLayerFields{i};
            if isfield(Y.canopy_air, name)
                Y.canopy_air.(name) = repmat(Y.canopy_air.(name)(:), 1, nCell);
            end
        end
    end

    % SCOPE soil-heat history: [12 x 2] -> [12 x 2 x nCell]
    if isfield(Y, 'scope_soil') && isfield(Y.scope_soil, 'Tsold')
        Y.scope_soil.Tsold = repmat(Y.scope_soil.Tsold, 1, 1, nCell);
    end

    if isfield(Y, 'soil_modes')
        Y.soil_modes = broadcastSoilModeState(Y.soil_modes, NN, NL, nD, nCell);
    end

    % nested scalar-leaf state containers -> [1 x nCell] per cell (plant pools + the state-ownership
    % migration scalars snow/surface/canopy_water, 2026-07-05).
    nested = {'plant_hydraulic','plant_cn','snow','surface','canopy_water'};
    for i = 1:numel(nested)
        name = nested{i};
        if isfield(Y, name)
            Y.(name) = broadcastNumericScalars(Y.(name), nCell);
        end
    end

    validateNodeField(Y.h, NN, nCell, 'Y.h');
end

% ---------------------------------------------------------------------------
function s = broadcastSoilModeState(s, NN, NL, nD, nCell)
% Broadcast nested soil-mode state arrays to regional shapes.
% Editor: Jizixing Zhu; Date: 2026-06-28.
    names = fieldnames(s);
    for i = 1:numel(names)
        name = names{i};
        value = s.(name);
        if isstruct(value)
            s.(name) = broadcastSoilModeState(value, NN, NL, nD, nCell);
        elseif isnumeric(value) && ~isempty(value)
            s.(name) = broadcastSoilModeNumeric(value, NN, NL, nD, nCell, name);
        end
    end
end

% ---------------------------------------------------------------------------
function value = broadcastSoilModeNumeric(value, NN, NL, nD, nCell, name)
% Broadcast vectors/matrices used by hysteresis and dual-porosity mode memory.
% Editor: Jizixing Zhu; Date: 2026-06-28.
    sz = size(value);
    if isvector(value) && numel(value) == NN
        value = repmat(value(:), 1, nCell);
    elseif isvector(value) && numel(value) == NL
        value = repmat(value(:), 1, nCell);
    elseif ismatrix(value) && size(value, 1) == NL && size(value, 2) >= nD
        value = repmat(value, 1, 1, nCell);
    elseif isscalar(value)
        value = repmat(value, 1, nCell);
    else
        error('expandBuiltToRegional:badSoilModeState', ...
            'Y.soil_modes.%s has unsupported shape [%s].', name, strjoin(string(sz), ' x '));
    end
end

% ---------------------------------------------------------------------------
function s = broadcastLayerStruct(s, NN, NL, nD, nS, nCell)
% Broadcast common layer/node parameter arrays while leaving true scalars alone.
% Editor: Jizixing Zhu; Date: 2026-06-27.
    names = fieldnames(s);
    for i = 1:numel(names)
        name = names{i};
        value = s.(name);
        if isstruct(value)
            s.(name) = broadcastLayerStruct(value, NN, NL, nD, nS, nCell);
        elseif isnumeric(value) && ~isempty(value)
            s.(name) = broadcastLayerNumeric(value, NN, NL, nD, nS, nCell);
        end
    end
end

% ---------------------------------------------------------------------------
function value = broadcastLayerNumeric(value, NN, NL, nD, nS, nCell)
% Broadcast arrays that clearly live on node/layer/solute dimensions.
% Editor: Jizixing Zhu; Date: 2026-06-27.
    sz = size(value);
    if isvector(value) && numel(value) == NN
        value = repmat(value(:), 1, nCell);
    elseif isvector(value) && numel(value) == NL
        value = repmat(value(:), 1, nCell);
    elseif isequal(sz, [NL, nD])
        value = repmat(value, 1, 1, nCell);
    elseif ndims(value) == 3 && isequal(sz, [1, NL, nS])
        value = repmat(value, 1, 1, 1, nCell);
    end
end

% ---------------------------------------------------------------------------
function forcing = broadcastForcing(forcing, nCell)
% Broadcast forcing time-series vectors to [nTime x nCell] arrays.
% Editor: Jizixing Zhu; Date: 2026-06-27.
    names = fieldnames(forcing);
    for i = 1:numel(names)
        name = names{i};
        value = forcing.(name);
        if isstruct(value)
            forcing.(name) = broadcastForcing(value, nCell);
        elseif isnumeric(value) && isvector(value) && numel(value) > 1
            forcing.(name) = repmat(value(:), 1, nCell);
        end
    end
end

% ---------------------------------------------------------------------------
function s = broadcastNumericScalars(s, nCell)
% Broadcast numeric scalar leaves in nested plant state containers.
% Editor: Jizixing Zhu; Date: 2026-06-27.
    names = fieldnames(s);
    for i = 1:numel(names)
        name = names{i};
        value = s.(name);
        if isstruct(value)
            s.(name) = broadcastNumericScalars(value, nCell);
        elseif isnumeric(value) && isscalar(value)
            s.(name) = repmat(value, 1, nCell);
        end
    end
end

% ---------------------------------------------------------------------------
function A = reshapeLayerField(A, NL, nD, nCell, name)
% Broadcast one NL x nD layer-state field to NL x nD x nCell.
% Editor: Jizixing Zhu; Date: 2026-06-27.
    if ~isequal(size(A), [NL, nD])
        error('expandBuiltToRegional:badLayerState', ...
            '%s must be [%d x %d] before regional expansion.', name, NL, nD);
    end
    A = repmat(A, 1, 1, nCell);
end

% ---------------------------------------------------------------------------
function validateNodeField(A, NN, nCell, name)
% Validate one regional node-state field.
% Editor: Jizixing Zhu; Date: 2026-06-27.
    if ~isequal(size(A), [NN, nCell])
        error('expandBuiltToRegional:badNodeState', ...
            '%s must be [%d x %d] after regional expansion.', name, NN, nCell);
    end
end
