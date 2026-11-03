function saveCompassResumeCheckpoint(reg, cfg, ckptDir, lastCompletedKT, failedKT, ResumeState)
%SAVECOMPASSRESUMECHECKPOINT  Persist outputs + resume state on outer-step rollback.
%
% Called by COMPASS.m on OUTER-STEP-ROLLBACK (Picard divergence). Saves
% (a) the diagnostic registry up to lastCompletedKT (NetCDF + legacy CSV) and
% (b) a resume .mat snapshot that lets a follow-up run continue from failedKT.
%
% Inputs:
%   reg               [in]      diagnostics.Registry (populated up to lastCompletedKT)
%   cfg               [in]      legacy ModelSettings (uses cfg.OutputPath)
%   ckptDir           [in] char Checkpoint directory (resume*.mat go here)
%   lastCompletedKT   [in] [-]  Last KT whose outputs are valid in the registry
%   failedKT          [in] [-]  KT that did NOT converge (will be retried)
%   ResumeState       [in] struct Snapshot caller assembled (SoilVariables,
%                                  PicardSummary, h2ocan_history, current_date,
%                                  current_day, current_hour, current_year,
%                                  SUMTIME, NoTime, DayHour, StepFailed, ...).
%                                  Must include .failed_KT, .failed_TIME,
%                                  .failed_reason fields.
%
% Side effects:
%   1. data/output_<run>/partial_before_nonconv_KT<failed>_cutoff_KT<last>/
%        ├── *.nc          (one per group, NetCDF)
%        └── Sim_output_*.csv  (legacy compat)
%   2. <ckptDir>/resume_before_KT<failedKT>.mat
%   3. <ckptDir>/resume_latest.mat (always overwritten — used by auto-resume)
%
% Author: Zhu Jizixing, 2026-05-15 (replaces evalin-based legacy version)
%   Previously used 40+ evalin('caller', ...) calls — violated coding standard §3.3.

    try
        % ---- (1) Output flush: registry → NetCDF + CSV ----
        outputFolder = io.resolveOutputFolder(cfg);
        tag = sprintf('partial_before_nonconv_KT%05d_cutoff_KT%05d', ...
            failedKT, lastCompletedKT);
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
        diagnostics.flushToNetCDF(reg, tagFolder, runMeta);
        diagnostics.flushToCsv(reg, tagFolder);
        fprintf('[PARTIAL-OUTPUT] Saved %s (last completed KT=%d).\n', ...
            tag, lastCompletedKT);

        % ---- (2) Resume state snapshot ----
        if ~exist(ckptDir, 'dir'); mkdir(ckptDir); end
        resumeFile = fullfile(ckptDir, sprintf('resume_before_KT%05d.mat', failedKT));
        save(resumeFile, 'ResumeState', '-v7.3');
        save(fullfile(ckptDir, 'resume_latest.mat'), 'ResumeState', '-v7.3');
        fprintf('[RESUME-CKPT] Saved resume state before KT=%d: %s\n', ...
            failedKT, resumeFile);
    catch err
        warning('COMPASS:ResumeCheckpointFailed', ...
            'Resume checkpoint save failed at KT=%d: %s', failedKT, err.message);
    end
end
