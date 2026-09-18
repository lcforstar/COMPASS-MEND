function outputFolder = resolveOutputFolder(cfg)
%RESOLVEOUTPUTFOLDER  Decide where COMPASS writes its run outputs.
%
% Inputs:
%   cfg     [in]    config struct (uses cfg.OutputPath as fallback)
%
% Resolution order:
%   1. io.getEnv('OUTPUT')  = COMPASS_OUTPUT, then STEMMUS_OUTPUT fallback
%   2. cfg.OutputPath
%
% Creates the directory if it does not exist.
%
% Author: Zhu Jizixing, 2026-05-15
%   Extracted from io.writeCompassOutputs to be reusable by
%   diagnostics.flushToNetCDF / flushToCsv path.
% Editor: Jizixing Zhu; Date: 2026-07-01.
% Change: 2026-07-01 Jizixing Zhu — env override via io.getEnv('OUTPUT') (COMPASS_ first, STEMMUS_ fallback) instead of raw STEMMUS_OUTPUT, so partial-output/flush honour the SAME dir as runtime.output/checkpoint (fixes the COMPASS_/STEMMUS_ split).

    outputEnv = io.getEnv('OUTPUT');   % COMPASS_OUTPUT, then STEMMUS_OUTPUT fallback
    if ~isempty(outputEnv)
        outputFolder = outputEnv;
    else
        outputFolder = cfg.OutputPath;
    end

    if ~exist(outputFolder, 'dir')
        mkdir(outputFolder);
    end
end
