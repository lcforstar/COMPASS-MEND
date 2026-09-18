function reg = loadDLMetricsToRegistry(metricsJsonPath, outputDir, runMeta)
%LOADDLMETRICSTOREGISTRY  Pump Python train_transolver.py output into NetCDF.
%
% Bridges the asymmetric "DL group" (epoch-indexed, written by Python) to
% the COMPASS diagnostics.Registry / flushToNetCDF pipeline so deep_learning.nc
% lands in the same data/output_<tag>/ folder as physical group outputs.
%
% Workflow:
%   Python train_transolver.py → metrics.json (history: list of per-epoch dict)
%                              → loadDLMetricsToRegistry (this file)
%                              → diagnostics.Registry buffer (epoch axis)
%                              → diagnostics.flushToNetCDF
%                              → data/output_<run>/deep_learning.nc
%
% Inputs:
%   metricsJsonPath  [in]  char    Absolute path to Python-emitted metrics.json
%   outputDir        [in]  char    Where deep_learning.nc + metadata.yaml land
%   runMeta          [in]  struct  Run-level metadata; passed to flushToNetCDF.
%                                  Should include: run_tag, start_time, end_time.
%                                  See COMPASS.m end-of-run flush for the canonical form.
%
% Returns:
%   reg  [out]  diagnostics.Registry  Populated registry (epoch dim sized to history)
%
% Side effects:
%   - Writes <outputDir>/deep_learning.nc (CF-1.10) via flushToNetCDF.
%   - Writes <outputDir>/metadata.yaml.
%   - Does NOT write CSVs (flushToCsv); DL training run is not part of the
%     1-month legacy CSV compat period.
%
% Example:
%   meta = struct('run_tag', 'transolver_freezethaw_full_year', ...
%                 'start_time', '2026-05-15T00:00', ...
%                 'end_time',   '2026-05-15T12:00', ...
%                 'delt_t_seconds', 0);  % epoch-indexed; no physical Δt
%   diagnostics.loadDLMetricsToRegistry( ...
%       'data/warmstart_transolver_freezethaw_extreme_full_year/metrics.json', ...
%       'data/warmstart_transolver_freezethaw_extreme_full_year', ...
%       meta);
%
% Author: Zhu Jizixing, 2026-05-15

    if nargin < 3
        runMeta = struct();
    end

    if ~exist(metricsJsonPath, 'file')
        error('loadDLMetricsToRegistry:NotFound', ...
            'metrics.json not found: %s', metricsJsonPath);
    end

    raw = fileread(metricsJsonPath);
    metrics = jsondecode(raw);
    if ~isfield(metrics, 'history') || isempty(metrics.history)
        error('loadDLMetricsToRegistry:NoHistory', ...
            'metrics.json missing history field or empty');
    end

    history = metrics.history;
    if iscell(history)
        % older jsondecode emits cell-of-struct when rows have mixed schemas
        n_epochs = numel(history);
        rowGetter = @(i) history{i};
    else
        n_epochs = numel(history);
        rowGetter = @(i) history(i);
    end
    if n_epochs == 0
        error('loadDLMetricsToRegistry:EmptyHistory', 'history has 0 rows');
    end

    % Build a minimal Registry just for the DL group.
    %   tS = n_epochs (placeholder; flushToNetCDF reads reg.tS but DL uses 'epoch' dim)
    %   NL / NN are unused by epoch-only fields but Registry constructor needs them
    reg = diagnostics.Registry(n_epochs, 1, 2);
    reg.registerDim('epoch', n_epochs);
    diagnostics.groups.registerDeepLearning(reg, []);

    % Python row dict → MATLAB registry field name mapping.
    % Most names match directly; these 5 differ because Python uses pre-mapping
    % legacy names while MATLAB registry uses unified k_/loss_*/train_loss_* convention.
    %   Python (train_transolver.py:1500-1551)    MATLAB (registerDeepLearning.m)
    alias = containers.Map( ...
        { 'train_mse_scaled', ...
          'train_data_mse_scaled', ...
          'train_pino_residual_loss', ...
          'val_noharm_pred_rmse_h_cm', ...
          'test_noharm_pred_rmse_h_cm' }, ...
        { 'train_loss_total', ...
          'train_loss_supervised', ...
          'train_loss_pino', ...
          'val_noharm_rmse_h_cm', ...
          'test_noharm_rmse_h_cm' } );

    % Forward each Python row → reg.write per registered DL field.
    dlNames = reg.listGroup('deep_learning');
    dlSet = containers.Map(dlNames, num2cell(true(size(dlNames))));
    skipped_unmapped = {};
    for i = 1:n_epochs
        row = rowGetter(i);
        epoch_kt = i;   % registry "kt" axis = epoch index
        rowFields = fieldnames(row);
        for j = 1:numel(rowFields)
            pyName = rowFields{j};
            % Resolve alias (Python → MATLAB)
            if alias.isKey(pyName)
                regName = alias(pyName);
            else
                regName = pyName;
            end
            % Skip if registry doesn't know this field (e.g. 'epoch', 'loss_mode')
            if ~dlSet.isKey(regName)
                continue;
            end
            value = row.(pyName);
            if isnumeric(value) && isscalar(value) && isfinite(value)
                reg.write(regName, value, epoch_kt);
            elseif isnumeric(value) && isscalar(value)
                % NaN/Inf — buffer stays NaN (already default)
                continue;
            else
                % non-scalar or non-numeric (e.g. loss_mode string) — skip
                skipped_unmapped{end+1} = pyName; %#ok<AGROW>
            end
        end
    end

    % Top-level metrics → derive scalar epoch-indexed fields that aren't in row.
    if isfield(metrics, 'best_checkpoint') && ischar(metrics.best_checkpoint)
        % best_checkpoint path like ".../best_epoch_53.pt"; parse the index out.
        tok = regexp(metrics.best_checkpoint, 'epoch[_-]?(\d+)', 'tokens', 'once');
        if ~isempty(tok)
            best_ep = str2double(tok{1});
            if isfinite(best_ep)
                % Write same value across every epoch (constant time series; the
                % "story" of when it crystallised is in selection_score min).
                for i = 1:n_epochs
                    reg.write('best_epoch', best_ep, i);
                end
            end
        end
    end

    if ~isempty(skipped_unmapped)
        u = unique(skipped_unmapped);
        fprintf('[dl-load] %d row fields skipped (non-scalar/string): %s\n', ...
            numel(u), strjoin(u, ', '));
    end

    if ~exist(outputDir, 'dir'); mkdir(outputDir); end

    % Default minimal runMeta (start_time used by flushToNetCDF for time units)
    if ~isfield(runMeta, 'run_tag')
        [~, runTag] = fileparts(outputDir);
        runMeta.run_tag = runTag;
    end
    if ~isfield(runMeta, 'start_time')
        runMeta.start_time = char(datetime('now', 'Format', 'yyyy-MM-dd''T''HH:mm'));
    end

    diagnostics.flushToNetCDF(reg, outputDir, runMeta);
    fprintf('[dl-load] wrote deep_learning.nc (n_epochs=%d, fields_with_data=%d) to %s\n', ...
        n_epochs, sum_with_data(reg, dlNames), outputDir);
end


function n = sum_with_data(reg, names)
% Count how many DL fields actually got non-NaN data from this JSON.
    n = 0;
    for k = 1:numel(names)
        buf = reg.get(names{k});
        if any(~isnan(buf(:)))
            n = n + 1;
        end
    end
end
