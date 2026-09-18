function regional = readRegionalStatic(ncPath, options)
%READREGIONALSTATIC  Read regional/grid static initialization fields from NetCDF.
%
%   regional = init.readRegionalStatic(ncPath, options) reads the horizontal
%   land-cell descriptor and optional initialization overrides from a NetCDF
%   static file. This is an initialization I/O adapter only: it does not build
%   parameters, update prognostic state, or advance the solver.
%
%   Required variables:
%     lat      [degree_north, any shape]   cell latitude
%     lon      [degree_east,  same shape]  cell longitude
%
%   Optional variables:
%     mask     [logical/int, same shape] land mask, default true
%     cell_id  [-, same shape]           cell identifiers, default 1:nAllCell
%
%   Optional options struct:
%     active_index [nBlock x 1]  subset of full-grid linear land indices to read
%     max_cells    [scalar]      first max_cells active land cells after masking
%
%   Optional override naming convention:
%     soil_<field>             -> built.parameters.soil.<field>
%     soil_n_<field>           -> built.parameters.soil_n.<field>
%     veg_root_RAI             -> built.parameters.veg.root.RAI
%     veg_PFT_type             -> built.parameters.veg.structure.PFT_type
%     Y_<field>                -> built.Y.<field>
%     Y_soil_gas_C_O2          -> built.Y.soil_gas.C_O2
%     forcing_hourly_<field>   -> built.forcing.hourly.<field>
%
%   The returned .space is land-packed when mask is present: only selected
%   mask==true cells become model columns. The original grid shape and active
%   indices are stored in .grid and .space so output writers can scatter columns
%   back to a 2-D/3-D map later. Arrays are oriented later by
%   init.applyRegionalStatic against the already built domain sizes. This reader
%   deliberately keeps raw NetCDF numeric values untouched so units remain the
%   responsibility of the static-file producer.
% Editor: Jizixing Zhu; Date: 2026-06-27.

    if nargin < 1 || ~(ischar(ncPath) || isstring(ncPath))
        error('readRegionalStatic:badPath', 'ncPath must be a NetCDF file path.');
    end
    if nargin < 2 || isempty(options)
        options = struct();
    end
    ncPath = char(ncPath);
    if ~isfile(ncPath)
        error('readRegionalStatic:fileNotFound', 'NetCDF static file not found: %s', ncPath);
    end

    vars = ncinfo(ncPath);
    varNames = string({vars.Variables.Name});
    requireVar(varNames, "lat", ncPath);
    requireVar(varNames, "lon", ncPath);

    latRaw = double(ncread(ncPath, 'lat'));
    lonRaw = double(ncread(ncPath, 'lon'));
    if ~isequal(size(latRaw), size(lonRaw))
        error('readRegionalStatic:latLonSize', ...
            'lat and lon must have the same shape in %s.', ncPath);
    end
    gridShape = size(latRaw);
    nAllCell = numel(latRaw);

    maskRaw = true(gridShape);
    if any(varNames == "mask")
        maskRaw = logical(ncread(ncPath, 'mask'));
    end
    if ~isequal(size(maskRaw), gridShape)
        error('readRegionalStatic:maskSize', ...
            'mask must have the same shape as lat/lon in %s.', ncPath);
    end
    fullActiveIndex = find(maskRaw(:));
    activeIndex = selectActiveIndex(fullActiveIndex, options, nAllCell);
    if isempty(activeIndex)
        error('readRegionalStatic:emptyMask', 'selected mask contains no active land cells in %s.', ncPath);
    end
    nCell = numel(activeIndex);

    cellIdRaw = reshape(1:nAllCell, gridShape);
    if any(varNames == "cell_id")
        cellIdRaw = ncread(ncPath, 'cell_id');
    end
    if ~isequal(size(cellIdRaw), gridShape)
        error('readRegionalStatic:cellIdSize', ...
            'cell_id must have the same shape as lat/lon in %s.', ncPath);
    end

    regional = struct();
    regional.path = ncPath;
    regional.grid = struct();
    regional.grid.shape = gridShape;
    regional.grid.nAllCell = nAllCell;
    regional.grid.mask = maskRaw;
    regional.grid.active_index = activeIndex;
    regional.grid.full_active_index = fullActiveIndex;
    regional.grid.pack = 'land_only';
    regional.grid.is_subset = numel(activeIndex) ~= numel(fullActiveIndex);
    regional.space = struct('lat', latRaw(activeIndex), 'lon', lonRaw(activeIndex), ...
        'mask', true(nCell, 1), 'cell_id', cellIdRaw(activeIndex), ...
        'grid_shape', gridShape, 'nAllCell', nAllCell, ...
        'active_index', activeIndex, 'full_active_index', fullActiveIndex, ...
        'pack', 'land_only', 'is_subset', regional.grid.is_subset);
    regional.overrides = struct('soil', struct(), 'soil_n', struct(), ...
        'veg', struct(), 'Y', struct(), 'forcing', struct());

    for i = 1:numel(varNames)
        name = char(varNames(i));
        if ismember(name, {'lat','lon','mask','cell_id'})
            continue;
        end
        value = ncread(ncPath, name);
        regional.overrides = addOverride(regional.overrides, name, value);
    end
end

% -------------------------------------------------------------------------
function activeIndex = selectActiveIndex(fullActiveIndex, options, nAllCell)
% Select an active land-cell subset for block/chunk initialization.
% Editor: Jizixing Zhu; Date: 2026-06-27.
    activeIndex = fullActiveIndex(:);
    if isfield(options, 'active_index') && ~isempty(options.active_index)
        requested = options.active_index(:);
        if any(requested < 1) || any(requested > nAllCell)
            error('readRegionalStatic:activeIndexRange', ...
                'options.active_index entries must be within 1:nAllCell.');
        end
        if ~all(ismember(requested, fullActiveIndex))
            error('readRegionalStatic:activeIndexMask', ...
                'options.active_index must refer only to mask==true land cells.');
        end
        activeIndex = requested;
    end
    if isfield(options, 'max_cells') && ~isempty(options.max_cells)
        maxCells = round(options.max_cells);
        if maxCells < 1
            error('readRegionalStatic:maxCells', 'options.max_cells must be positive.');
        end
        activeIndex = activeIndex(1:min(maxCells, numel(activeIndex)));
    end
end

% -------------------------------------------------------------------------
function requireVar(varNames, name, ncPath)
% Require one NetCDF variable by name.
% Editor: Jizixing Zhu; Date: 2026-06-27.
    if ~any(varNames == name)
        error('readRegionalStatic:missingVariable', ...
            'Required variable "%s" is missing in %s.', name, ncPath);
    end
end

% -------------------------------------------------------------------------
function overrides = addOverride(overrides, name, value)
% Route one NetCDF variable into the override container.
% Editor: Jizixing Zhu; Date: 2026-06-27.
    if startsWith(name, 'soil_n_')
        field = name(8:end);
        overrides.soil_n.(field) = value;
    elseif startsWith(name, 'soil_')
        field = name(6:end);
        overrides.soil.(field) = value;
    elseif strcmp(name, 'veg_root_RAI')
        if ~isfield(overrides.veg, 'root')
            overrides.veg.root = struct();
        end
        overrides.veg.root.RAI = value;
    elseif strcmp(name, 'veg_PFT_type')
        if ~isfield(overrides.veg, 'structure')
            overrides.veg.structure = struct();
        end
        overrides.veg.structure.PFT_type = value;
    elseif strcmp(name, 'Y_soil_gas_C_O2')
        if ~isfield(overrides.Y, 'soil_gas')
            overrides.Y.soil_gas = struct();
        end
        overrides.Y.soil_gas.C_O2 = value;
    elseif startsWith(name, 'Y_')
        field = name(3:end);
        overrides.Y.(field) = value;
    elseif startsWith(name, 'forcing_hourly_')
        field = name(16:end);
        if ~isfield(overrides.forcing, 'hourly')
            overrides.forcing.hourly = struct();
        end
        overrides.forcing.hourly.(field) = value;
    end
end
