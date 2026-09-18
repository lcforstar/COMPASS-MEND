%RUN_WINTER_WHEAT_FULL_YEAR  Background acceptance run for Changwu wheat.
% The output path can be overridden with COMPASS_WHEAT_ANNUAL_OUTPUT.
% Editor: Jizixing Zhu; Date: 2026-08-24.
root = fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(root, 'src'));
outputDir = getenv('COMPASS_WHEAT_ANNUAL_OUTPUT');
if isempty(outputDir)
    outputDir = fullfile(root, 'data', 'output', 'changwu_winter_wheat_full_year');
end
runCOMPASSBiogeochemistrySpinup(fullfile(root, 'data', 'sites', 'changwu'), ...
    cycles=1, n_steps=8760, output_dir=outputDir, progress_interval=24);
