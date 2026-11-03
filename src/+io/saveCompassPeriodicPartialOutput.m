function saveCompassPeriodicPartialOutput(reg, cfg, ~, lastCompletedKT, nextKT)
%SAVECOMPASSPERIODICPARTIALOUTPUT  Write a tagged partial NetCDF snapshot (crash-recovery).
%
% Inputs:
%   reg               [in]    diagnostics.Registry (already populated up to lastCompletedKT)
%   cfg               [in]    config struct (uses cfg.OutputPath)
%   ~                 [in]    char, checkpoint directory (reserved for future; ignored here)
%   lastCompletedKT   [in]    last KT whose outputs are valid in the registry
%   nextKT            [in]    KT about to start (used only in tag/log)
%
% Side effects:
%   Creates a tagged subdirectory under cfg.OutputPath with NetCDF only (no CSV).
%   [PERF 2026-06-15] The per-variable legacy CSV dump (diagnostics.flushToCsv)
%   was dropped from periodic snapshots: it re-wrote the whole tS-length buffer
%   of every variable every ~30 days via writematrix (~750 s/run, profiler
%   2026-06-15) and nothing downstream reads the periodic-partial CSVs. The
%   end-of-run flush in COMPASS.m still writes CSV for back-compat.
%
% Author: Zhu Jizixing, 2026-05-15 (replaces evalin-based legacy version)

    try
        outputFolder = io.resolveOutputFolder(cfg);
        tag = sprintf('periodic_partial_cutoff_KT%05d_before_KT%05d', ...
            lastCompletedKT, nextKT);
        tagFolder = fullfile(outputFolder, tag);
        if ~exist(tagFolder, 'dir'); mkdir(tagFolder); end

        runMeta = struct( ...
            'tag',                tag, ...
            'run_tag',            tag, ...
            'last_completed_index', lastCompletedKT, ...
            'start_time',         char(cfg.Time.start_time), ...
            'end_time',           char(cfg.Time.end_time), ...
            'delt_t_seconds',     cfg.Time.Delt_t, ...
            'cfg_lFreezing',      double(helpers.isfield_or_prop_true(cfg, 'lFreezing')), ...
            'cfg_lNefc',          double(helpers.isfield_or_prop_true(cfg, 'lNefc')), ...
            'cfg_lSoilairefc',    double(helpers.isfield_or_prop_true(cfg, 'lSoilairefc')));
        % [PERF 2026-06-15] Trim to 1:lastCompletedKT so the snapshot does not
        % serialise the trailing all-NaN future steps of the full-year buffer.
        diagnostics.flushToNetCDF(reg, tagFolder, runMeta, lastCompletedKT);
        % No flushToCsv here — NetCDF is enough for a recovery snapshot; the
        % legacy per-variable CSV stays only in the end-of-run flush.

        fprintf('[PARTIAL-OUTPUT] Saved %s (last completed KT=%d).\n', ...
            tag, lastCompletedKT);
    catch partialErr
        warning('COMPASS:PartialOutputFailed', ...
            'Periodic partial output failed before KT=%d: %s', nextKT, partialErr.message);
    end
end
