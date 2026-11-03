function [xNext, aa] = andersonMixUpdate(xOld, xPicard, aa, opts)
%ANDERSONMIXUPDATE One Anderson mixing update for a fixed-point iteration.
%
%   [xNext, aa] = solver.andersonMixUpdate(xOld, xPicard, aa, opts)
%
%   The caller provides one ordinary Picard step:
%
%       xPicard = G(xOld)
%
%   This function stores recent Picard outputs and residuals
%
%       r = G(x) - x
%
%   and returns an Anderson-corrected next iterate. It does not call G
%   itself, so it can be inserted after an existing block Picard pass without
%   wrapping the whole model in a function handle.
%
%   Required inputs:
%       xOld    - vector state before the Picard pass
%       xPicard - vector state after the Picard pass
%
%   Optional state:
%       aa      - history struct returned by the previous call. Use [] to
%                 reset history.
%
%   Optional opts fields:
%       memory         - number of history differences to use, default 3
%       damping        - fraction of the Anderson correction, default 0.5
%       regularization - Tikhonov regularization for the least-squares
%                        problem, default 0
%       maxGammaNorm   - reject AA when norm(gamma) is too large, default Inf
%       maxStepNorm    - maximum norm of xNext - xPicard, default Inf
%       minHistory     - minimum stored states before acceleration, default 2
%       enabled        - false returns the Picard iterate and updates history
%       reset          - true clears AA history before using this update
%       safeguardFactor       - reject AA if predicted residual norm grows
%                               beyond this factor, default Inf
%       previousResidualNorm  - residual norm from the previous accepted
%                               state, used by safeguardFactor
%
%   Notes:
%       Scale the state vector before calling this routine when variables
%       have different units, for example [hh/100; TT/10].
%
%   Editor: Jizixing Zhu; Date: 2026-07-14.
%   Change: 2026-07-14 Jizixing Zhu -- replace normal equations with a
%           conditioned least-squares solve and Picard fallback.

if nargin < 3 || isempty(aa)
    aa = struct();
end
if nargin < 4 || isempty(opts)
    opts = struct();
end

memory = getOption(opts, 'memory', 3);
damping = getOption(opts, 'damping', 0.5);
regularization = getOption(opts, 'regularization', 0);
minSingularRatio = getOption(opts, 'minSingularRatio', 1e-8);
maxGammaNorm = getOption(opts, 'maxGammaNorm', Inf);
maxStepNorm = getOption(opts, 'maxStepNorm', Inf);
minHistory = getOption(opts, 'minHistory', 2);
enabled = getOption(opts, 'enabled', true);
reset = getOption(opts, 'reset', false);
safeguardFactor = getOption(opts, 'safeguardFactor', Inf);
previousResidualNorm = getOption(opts, 'previousResidualNorm', []);

memory = max(1, round(memory));
damping = min(max(damping, 0), 1);
minHistory = max(2, round(minHistory));

xOld = xOld(:);
xPicard = xPicard(:);

if numel(xOld) ~= numel(xPicard)
    error('solver:andersonMixUpdate:SizeMismatch', ...
        'xOld and xPicard must have the same number of elements.');
end

if reset
    aa = resetHistory(aa);
end

residual = xPicard - xOld;

if ~isfield(aa, 'picard') || ~isfield(aa, 'residual')
    aa.picard = zeros(numel(xOld), 0);
    aa.residual = zeros(numel(xOld), 0);
end
if ~isfield(aa, 'stats')
    aa.stats = initializeStats();
end
aa.stats = ensureStatsFields(aa.stats);

aa.picard = [aa.picard, xPicard];
aa.residual = [aa.residual, residual];

maxStored = memory + 1;
if size(aa.picard, 2) > maxStored
    aa.picard = aa.picard(:, end-maxStored+1:end);
    aa.residual = aa.residual(:, end-maxStored+1:end);
end

aa.lastResidualNorm = norm(residual);
aa.lastAccepted = false;
aa.lastGamma = [];
aa.lastRejectReason = '';

if ~enabled || size(aa.picard, 2) < minHistory || damping == 0
    xNext = xPicard;
    aa.lastRejectReason = 'not_ready';
    return;
end

dResidual = diff(aa.residual, 1, 2);
dPicard = diff(aa.picard, 1, 2);

if isempty(dResidual) || any(~isfinite(dResidual(:))) || ...
        any(~isfinite(dPicard(:))) || any(~isfinite(residual(:)))
    xNext = xPicard;
    aa.lastRejectReason = 'nonfinite_history';
    aa.stats.nRejectNonfinite = aa.stats.nRejectNonfinite + 1;
    return;
end

% Reject dependent history before solving; normal equations would square the condition number.
[gamma, lsInfo] = solver.solveAndersonLeastSquares( ...
    dResidual, residual, regularization, minSingularRatio);
aa.lastSingularRatio = lsInfo.singularRatio;
if ~lsInfo.accepted
    xNext = xPicard;
    aa.stats.nRejectIllConditioned = aa.stats.nRejectIllConditioned + 1;
    aa = resetHistory(aa);
    aa.lastRejectReason = 'ill_conditioned_history';
    return;
end

gammaNorm = norm(gamma);
aa.lastGammaNorm = gammaNorm;

if any(~isfinite(gamma))
    xNext = xPicard;
    aa.lastRejectReason = 'nonfinite_gamma';
    aa.stats.nRejectNonfinite = aa.stats.nRejectNonfinite + 1;
    return;
end

if gammaNorm > maxGammaNorm
    xNext = xPicard;
    aa.lastRejectReason = 'gamma_norm';
    aa.stats.nRejectGammaNorm = aa.stats.nRejectGammaNorm + 1;
    return;
end

xAnderson = xPicard - dPicard * gamma;
xCandidate = xPicard + damping * (xAnderson - xPicard);
correction = xCandidate - xPicard;
correctionNorm = norm(correction);

if correctionNorm > maxStepNorm
    correction = correction * (maxStepNorm / correctionNorm);
    xCandidate = xPicard + correction;
    aa.stats.nStepClipped = aa.stats.nStepClipped + 1;
end

if any(~isfinite(xCandidate))
    xNext = xPicard;
    aa.lastRejectReason = 'nonfinite_candidate';
    aa.stats.nRejectNonfinite = aa.stats.nRejectNonfinite + 1;
    return;
end

candidateResidualNorm = norm(xCandidate - xOld);
aa.lastCandidateResidualNorm = candidateResidualNorm;
if ~isempty(previousResidualNorm) && isfinite(safeguardFactor) && ...
        candidateResidualNorm > safeguardFactor * previousResidualNorm
    xNext = xPicard;
    aa.lastRejectReason = 'safeguard';
    aa.stats.nRejectSafeguard = aa.stats.nRejectSafeguard + 1;
    return;
end

xNext = xCandidate;
aa.lastAccepted = true;
aa.lastGamma = gamma;
aa.lastCorrectionNorm = norm(xNext - xPicard);
aa.stats.nAccept = aa.stats.nAccept + 1;
end

function value = getOption(opts, name, defaultValue)
if isfield(opts, name) && ~isempty(opts.(name))
    value = opts.(name);
else
    value = defaultValue;
end
end

function aa = resetHistory(aa)
aa.picard = [];
aa.residual = [];
aa.lastAccepted = false;
aa.lastGamma = [];
aa.lastRejectReason = 'reset';
if ~isfield(aa, 'stats')
    aa.stats = initializeStats();
end
aa.stats = ensureStatsFields(aa.stats);
aa.stats.nReset = aa.stats.nReset + 1;
end

function stats = initializeStats()
stats.nAccept = 0;
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
