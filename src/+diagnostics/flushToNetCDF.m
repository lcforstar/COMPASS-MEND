function flushToNetCDF(reg, output_dir, run_metadata, maxKT)
%FLUSHTONETCDF  Write all registered diagnostics into per-group NetCDF files.
%
% Layout (per design doc docs/netcdf_output_design_decisions.md):
%   data/output_<tag>/
%     <group_name>.nc        one file per physical group
%     metadata.yaml          human-readable summary across all groups
%
% Conventions:
%   - CF-1.10 compliant
%   - dimensions: time (tS), node (NN), layer (NL), nD (2)
%   - time coord: float datetime, units = "hours since <start>"
%   - depth_node coord: float [NN], cm, positive="down"
%   - data type: single (float32) by default; conservation diagnostics use
%     double so cumulative closure residuals retain their acceptance precision
%
% Inputs:
%   reg          [in]  diagnostics.Registry
%   output_dir   [in]  char, e.g. 'data/output_v33_yearfull_logsor_tol50'
%   run_metadata [in]  struct with fields: run_tag, start_time, end_time,
%                      delt_t_seconds, cfg_lNefc, cfg_lFreezing, cfg_lSoilairefc
%   maxKT        [in]  (optional) [-] write only time rows 1:maxKT instead of the
%                      full tS buffer. Periodic partial snapshots pass this so they
%                      do not serialise the trailing all-NaN future steps. Default
%                      = reg.tS (full buffer == end-of-run flush behaviour).
% Change: 2026-08-30 Codex -- persist conservation-group diagnostics as
%   float64 to prevent false cumulative-closure failures from float32 output.

    if ~exist(output_dir, 'dir'); mkdir(output_dir); end

    % --- depth_node coord [cm], positive=down (top to bottom) ---
    %     COMPASS internal memory varies, but we always export top-to-bottom.
    %     If actual depth array isn't passed, fall back to placeholder (replaced
    %     by caller setting run_metadata.depth_node).
    if isfield(run_metadata, 'depth_node') && ~isempty(run_metadata.depth_node)
        depth_node = run_metadata.depth_node(:);
    else
        depth_node = linspace(0, 1000, reg.NN)';   % fallback uniform
    end

    % --- time-axis length: full buffer, or trimmed to maxKT for partial snapshots ---
    nT = reg.tS;
    if nargin >= 4 && ~isempty(maxKT)
        nT = max(1, min(round(maxKT), reg.tS));
    end

    % --- time coord [hours since start] ---
    start_time = getf(run_metadata, 'start_time', '2021-01-01T00:00:00');
    time_units = sprintf('hours since %s', start_time);
    time_vals  = (0:nT-1)';

    % --- write one .nc per group (按 group 子目录, 2026-05-16) ---
    % 布局: output_<tag>/nc/<group>.nc  (不再堆 root, 出图脚本可直接 listdir)
    nc_root = fullfile(output_dir, 'nc');
    if ~exist(nc_root, 'dir'); mkdir(nc_root); end
    group_names = reg.groupNames();
    for gi = 1:numel(group_names)
        gname = group_names{gi};
        var_names = reg.listGroup(gname);

        % [2026-05-16] Skip group if ALL fields are全 NaN (e.g. deep_learning never
        % populated by a non-DL run; soil_oxygen when O2 module idle).
        % 跟 flushToCsv.m 行为对齐, 避免 root 里有 20 个 .nc 但 csv/ 只有 18 个子目录的不对称.
        any_with_data = false;
        for vi = 1:numel(var_names)
            buf = reg.get(var_names{vi});
            if ~all(isnan(buf(:)))
                any_with_data = true;
                break;
            end
        end
        if ~any_with_data
            fprintf('[diag] skipped %s.nc (all fields NaN — group not active this run)\n', gname);
            continue;
        end

        nc_path = fullfile(nc_root, sprintf('%s.nc', gname));
        % nuke if exists (re-runs)
        if exist(nc_path, 'file'); delete(nc_path); end

        % [P1.5-a 2026-06-16] Accumulate all variable attributes here and write
        % them (plus the globals) in ONE netcdf define-mode session at the end
        % (writeNcAttrs), instead of ~1400 per-file high-level ncwriteatt calls
        % (each opens/syncs the file). Data + _FillValue still use nccreate/ncwrite.
        varAttrs = struct('name', {}, 'attrs', {});

        % --- write time coord ---
        nccreate(nc_path, 'time', ...
            'Dimensions', {'time', nT}, ...
            'Datatype', 'double');
        ncwrite(nc_path, 'time', time_vals);
        varAttrs = addVarAttrs(varAttrs, 'time', { ...
            'units',         time_units; ...
            'long_name',     'time'; ...
            'standard_name', 'time'; ...
            'calendar',      'standard'});

        % --- write depth_node coord (only if any var uses 'node' dim) ---
        any_node = false;
        for vi = 1:numel(var_names)
            md = reg.getMeta(var_names{vi});
            if any(strcmp(md.dims, 'node')); any_node = true; break; end
        end
        if any_node
            nccreate(nc_path, 'depth_node', ...
                'Dimensions', {'node', reg.NN}, ...
                'Datatype', 'single');
            ncwrite(nc_path, 'depth_node', single(depth_node));
            varAttrs = addVarAttrs(varAttrs, 'depth_node', { ...
                'units',     'cm'; ...
                'long_name', 'Depth of node from surface'; ...
                'positive',  'down'; ...
                'axis',      'Z'});
        end

        % --- write variables ---
        for vi = 1:numel(var_names)
            name = var_names{vi};
            md   = reg.getMeta(name);
            buf  = reg.get(name);

            % trim trailing (un-computed) steps for partial snapshots
            if nT < reg.tS
                if ndims(buf) <= 2
                    buf = buf(1:nT, :);
                else
                    buf = buf(1:nT, :, :);
                end
            end

            % NaN → fill_value
            buf(isnan(buf)) = md.fill_value;

            % Conservation residuals are often 1e-10--1e-6 differences of
            % large cumulative ledgers. Float32 quantization previously made
            % a closed NH4 ledger appear to fail by about 1e-3 ug cm^-2.
            dataType='single'; fillValue=single(md.fill_value); writeBuf=single(buf);
            if strcmp(md.group,'conservation')
                dataType='double'; fillValue=double(md.fill_value); writeBuf=double(buf);
            end
            % create with deflate=5 for compression
            nccreate(nc_path, name, ...
                'Dimensions', dimToMatlab(md.dims, reg, nT), ...
                'Datatype', dataType, ...
                'DeflateLevel', 5, ...
                'FillValue', fillValue);
            ncwrite(nc_path, name, writeBuf);

            % CF attributes (collected; written in one session below)
            attrs = {'long_name', md.long_name};
            if ~isempty(md.standard_name)
                attrs(end+1, :) = {'standard_name', md.standard_name}; %#ok<AGROW>
            end
            attrs(end+1, :) = {'units', md.units}; %#ok<AGROW>
            if ~isempty(md.comment)
                attrs(end+1, :) = {'comment', md.comment}; %#ok<AGROW>
            end
            attrs(end+1, :) = {'source', md.source}; %#ok<AGROW>
            if ~isempty(md.source_code)
                attrs(end+1, :) = {'source_code', md.source_code}; %#ok<AGROW>
            end
            varAttrs = addVarAttrs(varAttrs, name, attrs);
        end

        % --- global attributes (CF-1.10), collected ---
        globalAttrs = { ...
            'Conventions',     'CF-1.10'; ...
            'title',           sprintf('COMPASS output: %s', gname); ...
            'source',          'COMPASS v33 (Zhu Jizixing, Northwest A&F University)'; ...
            'run_tag',         getf(run_metadata, 'run_tag', 'unknown'); ...
            'start_time',      getf(run_metadata, 'start_time', ''); ...
            'end_time',        getf(run_metadata, 'end_time', ''); ...
            'delt_t_seconds',  double(getf(run_metadata, 'delt_t_seconds', 3600)); ...
            'cfg_lFreezing',   double(getf(run_metadata, 'cfg_lFreezing', 0)); ...
            'cfg_lNefc',       double(getf(run_metadata, 'cfg_lNefc', 0)); ...
            'cfg_lSoilairefc', double(getf(run_metadata, 'cfg_lSoilairefc', 0)); ...
            'institution',     'Northwest A&F University'; ...
            'contact',         'Zhu Jizixing <zjzx3100362274@gmail.com>, <3100362274@nwafu.edu.cn>'; ...
            'history',         sprintf('Created %s by diagnostics.flushToNetCDF', ...
                char(datetime('now', 'Format', 'yyyy-MM-dd HH:mm:ss')))};

        % --- write ALL attributes for this file in one define-mode session ---
        writeNcAttrs(nc_path, varAttrs, globalAttrs);

        fprintf('[diag] wrote %s (%d vars, %.2f MB)\n', nc_path, numel(var_names), ...
            mb(dir(nc_path)));
    end

    % --- metadata.yaml (human-readable cross-group summary) ---
    writeMetadataYaml(reg, output_dir, run_metadata, depth_node);
end


function dim_pairs = dimToMatlab(dim_names, reg, nT)
% Convert {'time', 'node'} → {'time', nT, 'node', NN} for nccreate.
% nT is the (possibly trimmed) time-axis length; defaults to reg.tS.
    if nargin < 3 || isempty(nT); nT = reg.tS; end
    dim_pairs = {};
    for k = 1:numel(dim_names)
        d = dim_names{k};
        switch d
            case 'time';  sz = nT;
            case 'node';  sz = reg.NN;
            case 'layer'; sz = reg.NL;
            case 'nD';    sz = 2;
            otherwise
                if reg.custom_dims.isKey(d)
                    sz = reg.custom_dims(d);
                else
                    error('diagnostics:flush:Dim', 'unknown dim %s', d);
                end
        end
        dim_pairs{end+1} = d;  %#ok<AGROW>
        dim_pairs{end+1} = sz; %#ok<AGROW>
    end
end


function v = getf(s, name, default)
    if isfield(s, name) && ~isempty(s.(name)); v = s.(name); else; v = default; end
end


function varAttrs = addVarAttrs(varAttrs, name, attrs)
% Append one variable's attribute set to the accumulator for writeNcAttrs.
%   name   [in]  char, NetCDF variable name (already created via nccreate)
%   attrs  [in]  cell {attName, attValue; ...}, one attribute per row
    varAttrs(end+1).name = name; %#ok<AGROW>
    varAttrs(end).attrs  = attrs;
end


function writeNcAttrs(nc_path, varAttrs, globalAttrs)
%WRITENCATTRS  Write all variable + global attributes in one define-mode session.
%
% [P1.5-a 2026-06-16] Replaces ~1400 per-file high-level ncwriteatt calls (each
% of which opens/syncs the file) with a single netcdf.open/reDef/putAtt/endDef/
% close session, cutting flushToNetCDF attribute-write cost ~10x. Attribute
% names/values/types are identical to the ncwriteatt path. _FillValue is NOT
% written here — it is set at creation via nccreate 'FillValue'.
%
%   nc_path      [in]  char, path to the already-created .nc file
%   varAttrs     [in]  struct array, fields .name (char) + .attrs (cell {n,v;...})
%   globalAttrs  [in]  cell {attName, attValue; ...} for the NC_GLOBAL group
    ncid = netcdf.open(nc_path, 'WRITE');
    closer = onCleanup(@() netcdf.close(ncid)); %#ok<NASGU>
    netcdf.reDef(ncid);
    for i = 1:numel(varAttrs)
        varid = netcdf.inqVarID(ncid, varAttrs(i).name);
        a = varAttrs(i).attrs;
        for k = 1:size(a, 1)
            netcdf.putAtt(ncid, varid, a{k, 1}, a{k, 2});
        end
    end
    gid = netcdf.getConstant('NC_GLOBAL');
    for k = 1:size(globalAttrs, 1)
        netcdf.putAtt(ncid, gid, globalAttrs{k, 1}, globalAttrs{k, 2});
    end
    netcdf.endDef(ncid);
end


function s = mb(dir_struct)
    if isempty(dir_struct); s = 0; else; s = dir_struct(1).bytes / 1024 / 1024; end
end


function writeMetadataYaml(reg, output_dir, run_metadata, depth_node)
    yaml_path = fullfile(output_dir, 'metadata.yaml');
    fid = fopen(yaml_path, 'w');
    if fid < 0; return; end
    cleanup = onCleanup(@() fclose(fid));

    fprintf(fid, '# COMPASS output metadata (generated by diagnostics.flushToNetCDF)\n');
    fprintf(fid, 'model: COMPASS\n');
    fprintf(fid, 'version: v33\n');
    fprintf(fid, 'run_tag: %s\n', getf(run_metadata, 'run_tag', 'unknown'));
    fprintf(fid, 'start_time: %s\n', getf(run_metadata, 'start_time', ''));
    fprintf(fid, 'end_time:   %s\n', getf(run_metadata, 'end_time', ''));
    fprintf(fid, 'delt_t_seconds: %d\n', double(getf(run_metadata, 'delt_t_seconds', 3600)));
    fprintf(fid, 'grid:\n');
    fprintf(fid, '  NL: %d\n', reg.NL);
    fprintf(fid, '  NN: %d\n', reg.NN);
    fprintf(fid, '  depth_top_cm: %g\n', depth_node(1));
    fprintf(fid, '  depth_bottom_cm: %g\n', depth_node(end));
    fprintf(fid, 'cfg:\n');
    fprintf(fid, '  lFreezing:    %d\n', double(getf(run_metadata, 'cfg_lFreezing', 0)));
    fprintf(fid, '  lNefc:        %d\n', double(getf(run_metadata, 'cfg_lNefc', 0)));
    fprintf(fid, '  lSoilairefc:  %d\n', double(getf(run_metadata, 'cfg_lSoilairefc', 0)));
    fprintf(fid, 'groups:\n');
    gnames = reg.groupNames();
    for gi = 1:numel(gnames)
        g = gnames{gi};
        vars = reg.listGroup(g);
        fprintf(fid, '  %s:\n', g);
        fprintf(fid, '    file: %s.nc\n', g);
        fprintf(fid, '    n_variables: %d\n', numel(vars));
        fprintf(fid, '    variables:\n');
        for vi = 1:numel(vars)
            md = reg.getMeta(vars{vi});
            fprintf(fid, '      %s:\n', md.name);
            fprintf(fid, '        long_name: %s\n', escape(md.long_name));
            if ~isempty(md.standard_name)
                fprintf(fid, '        standard_name: %s\n', md.standard_name);
            end
            fprintf(fid, '        units: %s\n', md.units);
            fprintf(fid, '        dimensions: [%s]\n', strjoin(md.dims, ', '));
            if ~isempty(md.source); fprintf(fid, '        source: %s\n', md.source); end
            if ~isempty(md.source_code); fprintf(fid, '        source_code: %s\n', md.source_code); end
            if ~isempty(md.comment); fprintf(fid, '        comment: |\n          %s\n', escape(md.comment)); end
        end
    end
end


function s = escape(s)
    s = strrep(s, '\n', ' ');
    s = strrep(s, '"', '''');
end
