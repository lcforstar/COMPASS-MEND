function root_sink = resolveCurrentSrt(built)
%RESOLVECURRENTSRT  Return the current staged root-water sink for soil solve.
%
%   Input:
%     built [in] Built model after vegetation.root.stageWaterUptake.
%
%   Return:
%     root_sink [1/s, NLx2] Current layer-endpoint root-water sink.
%
%   Prior-step committed flux and PlantCN are deliberately not fallback
%   sources. Missing current Srt is a staging contract violation.
%   Editor: Jizixing Zhu; Date: 2026-07-16.
    assert(isfield(built.workspace, 'scope_step') ...
        && isfield(built.workspace.scope_step, 'Srt') ...
        && ~isempty(built.workspace.scope_step.Srt), ...
        'runOneStep:missingCurrentSrt', ...
        'vegetation.root.stageWaterUptake must produce current Srt before solve.');
    root_sink = built.workspace.scope_step.Srt;
end

