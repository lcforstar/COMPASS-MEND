function [xGuess, aa] = andersonApplyGuess(xCurrent, aa, opts)
%ANDERSONAPPLYGUESS Apply an Anderson extrapolated initial guess.
%
%   [xGuess, aa] = solver.andersonApplyGuess(xCurrent, aa, opts)
%
%   This is the aa-master style "apply" phase. It is called at the start of
%   a Picard iteration, before evaluating the model map. The map history is
%   supplied by solver.andersonCommitMap after previous Picard passes.
%
%   Editor: Jizixing Zhu; Date: 2026-07-14.
%   Change: 2026-07-14 Jizixing Zhu -- replace normal equations with a
%           conditioned least-squares solve and Picard fallback.

if nargin < 2 || isempty(aa)
    aa = struct();
end
if nargin < 3 || isempty(opts)
    opts = struct();
end

xCurrent = xCurrent(:);
enabled = getOption(opts, 'enabled', true);
memory = max(1, round(getOption(opts, 'memory', 3)));
minHistory = max(2, round(getOption(opts, 'minHistory', memory)));
damping = min(max(getOption(opts, 'damping', 0.5), 0), 1);
regularization = getOption(opts, 'regularization', 0);
minSingularRatio = getOption(opts, 'minSingularRatio', 1e-8);
maxGammaNorm = getOption(opts, 'maxGammaNorm', Inf);
maxStepNorm = getOption(opts, 'maxStepNorm', Inf);

aa = ensureState(aa);
aa.lastApplied = false;
aa.lastApplyRejectReason = '';

if ~enabled || damping == 0 || size(aa.output, 2) < minHistory
    xGuess = xCurrent;
    aa.lastApplyRejectReason = 'not_ready';
    return;
end

if numel(xCurrent) ~= size(aa.output, 1)
    xGuess = xCurrent;
    aa = resetHistory(aa);
    aa.lastApplyRejectReason = 'size_changed';
    return;
end

% Keep only the requested memory window before solving the small LS problem.
maxStored = memory + 1;
if size(aa.output, 2) > maxStored
    aa.input = aa.input(:, end-maxStored+1:end);
    aa.output = aa.output(:, end-maxStored+1:end);
    aa.residual = aa.residual(:, end-maxStored+1:end);
end

dResidual = diff(aa.residual, 1, 2);
dOutput = diff(aa.output, 1, 2);
residual = aa.residual(:, end);
base = aa.output(:, end);

if isempty(dResidual) || any(~isfinite(dResidual(:))) || ...
        any(~isfinite(dOutput(:))) || any(~isfinite(residual(:)))
    xGuess = xCurrent;
    aa.lastApplyRejectReason = 'nonfinite_history';
    aa.stats.nRejectNonfinite = aa.stats.nRejectNonfinite + 1;
    return;
end

% Reject dependent history before solving; normal equations would square the condition number.
[gamma, lsInfo] = solver.solveAndersonLeastSquares( ...
    dResidual, residual, regularization, minSingularRatio);
aa.lastSingularRatio = lsInfo.singularRatio;
if ~lsInfo.accepted
    xGuess = xCurrent;
    aa.stats.nRejectIllConditioned = aa.stats.nRejectIllConditioned + 1;
    aa = resetHistory(aa);
    aa.lastApplyRejectReason = 'ill_conditioned_history';
    return;
end

gammaNorm = norm(gamma);
aa.lastGammaNorm = gammaNorm;

if any(~isfinite(gamma))
    xGuess = xCurrent;
    aa.lastApplyRejectReason = 'nonfinite_gamma';
    aa.stats.nRejectNonfinite = aa.stats.nRejectNonfinite + 1;
    return;
end

if gammaNorm > maxGammaNorm
    xGuess = xCurrent;
    aa.lastApplyRejectReason = 'gamma_norm';
    aa.stats.nRejectGammaNorm = aa.stats.nRejectGammaNorm + 1;
    return;
end

xAnderson = base - dOutput * gamma;
xCandidate = base + damping * (xAnderson - base);
correction = xCandidate - xCurrent;
correctionNorm = norm(correction);

if correctionNorm > maxStepNorm
    correction = correction * (maxStepNorm / correctionNorm);
    xCandidate = xCurrent + correction;
    aa.stats.nStepClipped = aa.stats.nStepClipped + 1;
end

if any(~isfinite(xCandidate))
    xGuess = xCurrent;
    aa.lastApplyRejectReason = 'nonfinite_candidate';
    aa.stats.nRejectNonfinite = aa.stats.nRejectNonfinite + 1;
    return;
end

xGuess = xCandidate;
aa.lastApplied = true;
aa.lastInputBeforeApply = xCurrent;
aa.lastGuess = xGuess;
aa.lastReferenceResidualNorm = aa.lastResidualNorm;
aa.lastCorrectionNorm = norm(xGuess - xCurrent);
aa.stats.nApply = aa.stats.nApply + 1;
end

function aa = ensureState(aa)
if ~isfield(aa, 'input'), aa.input = []; end
if ~isfield(aa, 'output'), aa.output = []; end
if ~isfield(aa, 'residual'), aa.residual = []; end
if ~isfield(aa, 'stats'), aa.stats = initializeStats(); end
aa.stats = ensureStatsFields(aa.stats);
if ~isfield(aa, 'lastResidualNorm'), aa.lastResidualNorm = Inf; end
end

function aa = resetHistory(aa)
aa.input = [];
aa.output = [];
aa.residual = [];
aa.lastResidualNorm = Inf;
aa.lastApplied = false;
if ~isfield(aa, 'stats'), aa.stats = initializeStats(); end
aa.stats = ensureStatsFields(aa.stats);
aa.stats.nReset = aa.stats.nReset + 1;
end

function stats = initializeStats()
stats.nApply = 0;
stats.nCommit = 0;
stats.nRejectNonfinite = 0;
stats.nRejectIllConditioned = 0;
stats.nRejectGammaNorm = 0;
stats.nRejectSafeguard = 0;
stats.nStepClipped = 0;
stats.nReset = 0;
end

function stats = ensureStatsFields(stats)
defaults = initializeStats();
names = fieldnames(defaults);
for i = 1:numel(names)
    if ~isfield(stats, names{i})
        stats.(names{i}) = defaults.(names{i});
    end
end
end

function value = getOption(opts, name, defaultValue)
if isfield(opts, name) && ~isempty(opts.(name))
    value = opts.(name);
else
    value = defaultValue;
end
end
