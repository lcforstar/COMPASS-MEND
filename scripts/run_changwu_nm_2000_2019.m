%RUN_CHANGWU_NM_2000_2019  Full NM treatment simulation and C/N acceptance run.
%
% The configured interval is 2000-06-15 00:00 through 2019-06-15 00:00.
% Annual urea, cattle-manure, and straw events are read exclusively from the
% Changwu site.toml. The terminal 2019-06-15 10:00 management event lies after
% the observation endpoint and is intentionally not applied.
%
% Editor: Codex; Date: 2026-09-03.

root = fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(root, 'src'));

siteDir = fullfile(root, 'data', 'sites', 'changwu');
startTime = datetime(2000, 6, 15, 0, 0, 0);
endTime = datetime(2019, 6, 15, 0, 0, 0);
nSteps = hours(endTime - startTime);
assert(nSteps == 166536, 'Unexpected 2000--2019 hourly step count.');

outputDir = getenv('COMPASS_NM_2000_2019_OUTPUT');
if isempty(outputDir)
    runStamp = char(datetime('now', 'Format', 'yyyyMMdd_HHmmss'));
    outputDir = fullfile(root, 'data', 'output', ...
        ['changwu_nm_manure_2000_2019_', runStamp]);
end

result = runCOMPASSBiogeochemistrySpinup(siteDir, ...
    cycles=1, n_steps=nSteps, output_dir=outputDir, ...
    progress_interval=24 * 30, auto_stop=false);

save(fullfile(outputDir, 'full_period_acceptance_result.mat'), ...
    'result', 'startTime', 'endTime', 'nSteps', '-v7.3');
fprintf('[NM 2000--2019 COMPLETE] output=%s\n', outputDir);
