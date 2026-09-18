classdef Registry < handle
%REGISTRY  COMPASS diagnostic output registry (ClimaLand-style register/write/flush).
%
% Replaces the 103-place scattered `Sim_output.X(out_idx, ...) = ...` pattern with
% a centralised registry, following ClimaLand's `define_diagnostics.jl` design
% and the project coding standard §7 (diagnostics container convention).
%
% Lifecycle (3 phases):
%
%   1. REGISTER (once, before main loop):
%        reg = diagnostics.Registry(tS, NL, NN);
%        diagnostics.registerAllDiagnostics(reg, cfg);
%
%   2. WRITE (every KT inside main loop):
%        reg.write('soil_moisture', SoilVariables.Theta_LL, KT);
%        reg.write('nitrification_flux', NitrogenFlux.cvSnit(KT), KT);
%
%   3. FLUSH (once, after main loop):
%        diagnostics.flushToNetCDF(reg, output_dir, run_metadata);
%        diagnostics.flushToCsv(reg, output_dir);    % compatibility for old plots
%
% Output: data/output_<tag>/<group>.nc per physical group + metadata.yaml.
%   File layout follows the 8 design decisions in
%   docs/netcdf_output_design_decisions.md.
%
% Author: Zhu Jizixing, 2026-05-15
%
% [PERF 2026-06-15] Hot-path storage rewrite. The original design held each
% time-buffer inside a `containers.Map` (name → array). Every write() did
%     buf = obj.buffer(name);   % Map get → copies the FULL (tS × dims) array
%     buf(kt, :) = value;
%     obj.buffer(name) = buf;   % Map set → copies the FULL array back
% so for node-dimensioned variables (tS × NN ≈ 8760 × 257 ≈ 18 MB) each of the
% ~185 writes/step deep-copied ~36 MB. Over 8760 steps this dominated runtime
% (~2000 s self-time, profiler 2026-06-15). Buffers now live in a cell array
% indexed by an integer id, so `obj.buf_data{id}(kt,:) = value` updates the row
% in place without copying the whole array. Public API (register/write/get/
% getMeta/list/listGroup/groupNames/registerDim) is unchanged.

    properties (SetAccess = private)
        tS              % [-] total time-step count (length of time axis)
        NL              % [-] soil layer count
        NN              % [-] soil node count (= NL + 1)
        entries         % containers.Map: name → metadata struct (cold path: register/getMeta/flush)
        groups          % containers.Map: group_name → cell array of variable names
        custom_dims     % containers.Map: dim_name → integer size (for canopy_layer etc.)

        % --- hot-path storage (id-indexed, in-place writes; see header note) ---
        name2id         % containers.Map: name → integer index into the cell arrays below
        buf_data        % cell: id → numeric array, size = [tS, dims_after_transform]
        buf_transform   % cell: id → transform char (cached so write() skips the entries lookup)
        buf_ndim        % [1×n] cached ndims of each buffer (2 or 3)
        buf_ncol        % [1×n] cached size(buf,2) for 2-D buffers (shape validation in write)
        n_vars          % [-] number of registered variables (= numel(buf_data))
        warn_unknown    % [logical] report writes omitted by a sparse output profile
    end

    methods
        function obj = Registry(tS, NL, NN, warnUnknown)
            %REGISTRY  Construct the registry with run dimensions.
            %
            % Inputs:
            %   tS  [in] [-] total outer time steps
            %   NL  [in] [-] number of soil layers
            %   NN  [in] [-] number of soil nodes (= NL + 1 typically)

            if nargin < 3
                error('diagnostics:Registry:Args', 'Registry(tS, NL, NN) needs 3 args');
            end
            if nargin < 4
                warnUnknown = true;
            end
            obj.tS = tS;
            obj.NL = NL;
            obj.NN = NN;
            obj.entries     = containers.Map();
            obj.groups      = containers.Map();
            obj.custom_dims = containers.Map();

            obj.name2id       = containers.Map();
            obj.buf_data      = {};
            obj.buf_transform = {};
            obj.buf_ndim      = [];
            obj.buf_ncol      = [];
            obj.n_vars        = 0;
            obj.warn_unknown  = logical(warnUnknown);
        end

        function registerDim(obj, dim_name, dim_size)
            %REGISTERDIM  Declare a custom dimension before any variable uses it.
            %   Use for canopy_layer (=30 in SCOPE), wavelength bins, etc.
            %
            %   reg.registerDim('canopy_layer', 30);
            obj.custom_dims(dim_name) = dim_size;
        end

        function register(obj, entry)
            %REGISTER  Add a diagnostic variable definition.
            %
            % `entry` is a struct with fields (CF-Convention style):
            %   name          [char]  short snake_case identifier (also NetCDF var name)
            %   long_name     [char]  human-readable name
            %   standard_name [char]  CF standard name (or '' if no match)
            %   units         [char]  CF UDUNITS-2 string (e.g. 'm3 m-3', 'cm s-1')
            %   dims          [cell]  NetCDF dimension names {'time','node'} / {'time'} / ...
            %   group         [char]  physical group, → <group>.nc file
            %   source        [char]  in-memory source (e.g. 'SoilVariables.Theta_LL')
            %   source_code   [char]  file:line where it's set
            %   comment       [char]  free text (optional)
            %   fill_value    [num]   missing-value marker (default -9999)
            %   transform     [char]  pre-write op: 'none' | 'edge_to_node' | 'flip_node' (default 'none')
            %
            % After register, buf_data{id} is allocated as nan(tS, ...).

            if ~isfield(entry, 'name') || isempty(entry.name)
                error('diagnostics:Registry:NoName', 'entry must have a .name');
            end
            if obj.entries.isKey(entry.name)
                warning('diagnostics:Registry:Duplicate', ...
                    'Variable "%s" registered twice; second registration ignored', entry.name);
                return;
            end

            % --- fill defaults ---
            if ~isfield(entry, 'long_name'),     entry.long_name = entry.name; end
            if ~isfield(entry, 'standard_name'), entry.standard_name = ''; end
            if ~isfield(entry, 'units'),         entry.units = ''; end
            if ~isfield(entry, 'dims'),          entry.dims = {'time'}; end
            if ~isfield(entry, 'group'),         entry.group = 'misc'; end
            if ~isfield(entry, 'source'),        entry.source = ''; end
            if ~isfield(entry, 'source_code'),   entry.source_code = ''; end
            if ~isfield(entry, 'comment'),       entry.comment = ''; end
            if ~isfield(entry, 'fill_value'),    entry.fill_value = -9999; end
            if ~isfield(entry, 'transform'),     entry.transform = 'none'; end

            % --- compute buffer shape from dims ---
            shape = obj.computeShape(entry.dims);
            entry.shape = shape;

            % --- allocate buffer + register hot-path metadata ---
            obj.entries(entry.name) = entry;
            id = obj.n_vars + 1;
            obj.n_vars              = id;
            obj.name2id(entry.name) = id;
            obj.buf_data{id}        = nan(shape);   % NaN sentinel; flush replaces with fill_value
            obj.buf_transform{id}   = entry.transform;
            obj.buf_ndim(id)        = numel(shape); % computeShape pads to ≥2, so 2 or 3
            obj.buf_ncol(id)        = shape(2);

            % --- index by group ---
            if obj.groups.isKey(entry.group)
                grp = obj.groups(entry.group);
                grp{end+1} = entry.name;
                obj.groups(entry.group) = grp;
            else
                obj.groups(entry.group) = {entry.name};
            end
        end

        function write(obj, name, value, kt)
            %WRITE  Store one time-step value for variable `name`.
            %
            %   reg.write('soil_moisture', SoilVariables.Theta_LL, KT)
            %
            % `value` shape must match (entry.shape without time dim).
            % Transforms (e.g. edge_to_node for NL×2 → NN) are applied here.
            %
            % Hot path: one scalar Map lookup (name→id), then an in-place row
            % assignment into buf_data{id} — no full-buffer copy (see header note).

            if ~obj.name2id.isKey(name)
                % silent ignore to avoid breaking main loop on typos —
                % but log it (once per name) so the user sees it.
                if obj.warn_unknown
                    obj.warnUnknown(name);
                end
                return;
            end
            id = obj.name2id(name);

            % apply transform only when needed ('none' is the common case)
            tr = obj.buf_transform{id};
            if ~strcmp(tr, 'none')
                value = diagnostics.transforms.apply(value, tr, obj.NL, obj.NN);
            end

            % slot value into row kt — in-place, no whole-array copy
            if obj.buf_ndim(id) == 2
                ncol = obj.buf_ncol(id);
                if isscalar(value) && ncol == 1
                    obj.buf_data{id}(kt, 1) = value;
                elseif numel(value) == ncol
                    obj.buf_data{id}(kt, :) = value(:).';
                else
                    error('diagnostics:Registry:Shape', ...
                        'write(%s): value shape mismatch (expected col %d, got %s)', ...
                        name, ncol, mat2str(size(value)));
                end
            else  % buf_ndim(id) == 3
                obj.buf_data{id}(kt, :, :) = value;
            end
        end

        function val = get(obj, name)
            %GET  Retrieve the full accumulated buffer for one variable.
            if ~obj.name2id.isKey(name)
                error('diagnostics:Registry:Unknown', 'unknown variable %s', name);
            end
            val = obj.buf_data{obj.name2id(name)};
        end

        function md = getMeta(obj, name)
            %GETMETA  Get the metadata struct for one variable.
            md = obj.entries(name);
        end

        function names = list(obj)
            %LIST  All registered variable names.
            names = obj.entries.keys;
        end

        function names = listGroup(obj, group_name)
            %LISTGROUP  Variable names in one group, or {} if group not registered.
            if obj.groups.isKey(group_name)
                names = obj.groups(group_name);
            else
                names = {};
            end
        end

        function gnames = groupNames(obj)
            %GROUPNAMES  All registered group names.
            gnames = obj.groups.keys;
        end
    end

    methods (Access = private)
        function shape = computeShape(obj, dims)
            % Map dim names → integer sizes
            shape = zeros(1, numel(dims));
            for k = 1:numel(dims)
                d = dims{k};
                switch d
                    case 'time'
                        shape(k) = obj.tS;
                    case 'node'
                        shape(k) = obj.NN;
                    case 'layer'
                        shape(k) = obj.NL;
                    case 'nD'
                        shape(k) = 2;
                    otherwise
                        if obj.custom_dims.isKey(d)
                            shape(k) = obj.custom_dims(d);
                        else
                            error('diagnostics:Registry:Dim', ...
                                'unknown dimension %s — register it first via reg.registerDim()', d);
                        end
                end
            end
            if numel(shape) == 1
                % scalar-time series needs to be (tS, 1) not (tS,)
                shape(2) = 1;
            end
        end

        function warnUnknown(obj, name)
            persistent warned
            if isempty(warned), warned = containers.Map(); end
            if ~warned.isKey(name)
                warning('diagnostics:Registry:Unknown', ...
                    'write() called for unregistered variable "%s" (silenced after first hit)', name);
                warned(name) = true;  %#ok<NASGU>
            end
        end
    end
end
