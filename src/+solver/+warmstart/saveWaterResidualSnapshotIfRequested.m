function saveWaterResidualSnapshotIfRequested(WaterMatrix, RHS, SoilVariables, CHKWATER, cfg, KT, Delt_t, stage)
%SAVEWATERRESIDUALSNAPSHOTIFREQUESTED Export assembled water residual data.
% Part of the +deep_learning/+warmstart package. Called optionally from
% +soil/+water/solveSoilMoistureBalance to export the assembled tridiagonal
% system (C4/C4_a/RHS) plus current hh, for Python-side PINO training and
% MATLAB <-> Python residual alignment.
%
% Call:
%   solver.warmstart.saveWaterResidualSnapshotIfRequested( ...
%       WaterMatrix, RHS, SoilVariables, CHKWATER, cfg, KT, Delt_t, stage)
%
% Disabled by default. Enable with:
%   STEMMUS_WATER_RESIDUAL_SNAPSHOT_DIR=F:\COMPASSV1.0\data\water_residual_snapshots
%
% Optional filters:
%   STEMMUS_WATER_RESIDUAL_SNAPSHOT_KT=123      % match a single KT
%   STEMMUS_WATER_RESIDUAL_SNAPSHOT_STAGE=after_solve   % "before_solve" / "after_solve"
%   STEMMUS_WATER_RESIDUAL_SNAPSHOT_MAX=20      % cap total files written (default 200)
%
% The snapshot is for Python/MATLAB residual alignment. It saves the exact
% post-boundary C4/C4_a/RHS system and the current hh profile. Each KT/stage
% combination may produce multiple files because the Picard loop and adaptive
% sub-step can re-enter solveSoilMoistureBalance several times per outer KT.
% The "call" counter in the file name reflects that physical re-entry order.

out_dir = getenv('STEMMUS_WATER_RESIDUAL_SNAPSHOT_DIR');
if isempty(out_dir)
    return;
end

kt_filter = getenv('STEMMUS_WATER_RESIDUAL_SNAPSHOT_KT');
if ~isempty(kt_filter) && str2double(kt_filter) ~= KT
    return;
end

stage_filter = getenv('STEMMUS_WATER_RESIDUAL_SNAPSHOT_STAGE');
if ~isempty(stage_filter) && ~strcmp(stage_filter, stage)
    return;
end

if ~exist(out_dir, 'dir')
    mkdir(out_dir);
end

% Cap total files written to avoid runaway disk usage during full-year runs
max_files_env = getenv('STEMMUS_WATER_RESIDUAL_SNAPSHOT_MAX');
if isempty(max_files_env)
    max_files = 200;
else
    max_files = str2double(max_files_env);
    if ~isfinite(max_files) || max_files <= 0
        max_files = 200;
    end
end

persistent call_count
if isempty(call_count)
    call_count = 0;
end
call_count = call_count + 1;
if call_count > max_files
    return;
end

hh = SoilVariables.hh(:);
linear_residual = solver.warmstart.calculateLinearWaterResidual(WaterMatrix, RHS, hh, cfg);
C4 = WaterMatrix.C4;
C4_a = WaterMatrix.C4_a;
RHS_snapshot = RHS(:);
indxBotmLayer = cfg.Grid.indxBotmLayer;
NN = cfg.Grid.NN;
NL = cfg.Grid.NL;

path = fullfile(out_dir, sprintf('water_residual_%s_KT%05d_call%04d.mat', stage, KT, call_count));
save(path, ...
    'C4', 'C4_a', 'RHS_snapshot', 'hh', 'linear_residual', ...
    'CHKWATER', 'Delt_t', 'KT', 'stage', 'indxBotmLayer', 'NN', 'NL');
end
