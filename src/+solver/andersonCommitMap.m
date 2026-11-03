function aa = andersonCommitMap(xInput, xMap, aa, opts)
%ANDERSONCOMMITMAP Store one Picard map evaluation for Anderson history.
%
%   aa = solver.andersonCommitMap(xInput, xMap, aa, opts)
%
%   xInput is the state sent into the current Picard map, and xMap is the
%   state returned by the ordinary water-air-heat Picard pass.

if nargin < 3 || isempty(aa)
    aa = struct();
end
if nargin < 4 || isempty(opts)
    opts = struct();
end

xInput = xInput(:);
xMap = xMap(:);

enabled = getOption(opts, 'enabled', true);
memory = max(1, round(getOption(opts, 'memory', 3)));
reset = getOption(opts, 'reset', false);
safeguardFactor = getOption(opts, 'safeguardFactor', Inf);

aa = ensureState(aa);
aa.lastCommitRejectReason = '';

if reset
    aa = resetHistory(aa);
end

if numel(xInput) ~= numel(xMap)
    error('solver:andersonCommitMap:SizeMismatch', ...
        'xInput and xMap must have the same number of elements.');
end

if ~enabled
    return;
end

residual = xMap - xInput;
residualNorm = norm(residual);

if any(~isfinite(xInput)) || any(~isfinite(xMap)) || ~isfinite(residualNorm)
    aa = resetHistory(aa);
    aa.lastCommitRejectReason = 'nonfinite_map';
    aa.stats.nRejectNonfinite = aa.stats.nRejectNonfinite + 1;
    return;
end

if isfield(aa, 'lastApplied') && aa.lastApplied && ...
        isfinite(safeguardFactor) && isfinite(aa.lastReferenceResidualNorm) && ...
        residualNorm > safeguardFactor * aa.lastReferenceResidualNorm
    aa = resetHistory(aa);
    aa.lastCommitRejectReason = 'safeguard';
    aa.stats.nRejectSafeguard = aa.stats.nRejectSafeguard + 1;
    return;
end

aa.input = [aa.input, xInput];
aa.output = [aa.output, xMap];
aa.residual = [aa.residual, residual];

maxStored = memory + 1;
if size(aa.output, 2) > maxStored
    aa.input = aa.input(:, end-maxStored+1:end);
    aa.output = aa.output(:, end-maxStored+1:end);
    aa.residual = aa.residual(:, end-maxStored+1:end);
end

aa.lastResidualNorm = residualNorm;
aa.lastApplied = false;
aa.stats.nCommit = aa.stats.nCommit + 1;
end

function aa = ensureState(aa)
if ~isfield(aa, 'input'), aa.input = []; end
if ~isfield(aa, 'output'), aa.output = []; end
if ~isfield(aa, 'residual'), aa.residual = []; end
if ~isfield(aa, 'stats'), aa.stats = initializeStats(); end
if ~isfield(aa, 'lastResidualNorm'), aa.lastResidualNorm = Inf; end
if ~isfield(aa, 'lastApplied'), aa.lastApplied = false; end
if ~isfield(aa, 'lastReferenceResidualNorm'), aa.lastReferenceResidualNorm = Inf; end
end

function aa = resetHistory(aa)
aa.input = [];
aa.output = [];
aa.residual = [];
aa.lastResidualNorm = Inf;
aa.lastApplied = false;
aa.lastReferenceResidualNorm = Inf;
if ~isfield(aa, 'stats')
    aa.stats = initializeStats();
end
aa.stats.nReset = aa.stats.nReset + 1;
end

function stats = initializeStats()
stats.nApply = 0;
stats.nCommit = 0;
stats.nRejectNonfinite = 0;
stats.nRejectGammaNorm = 0;
stats.nRejectSafeguard = 0;
stats.nStepClipped = 0;
stats.nReset = 0;
end

function value = getOption(opts, name, defaultValue)
if isfield(opts, name) && ~isempty(opts.(name))
    value = opts.(name);
else
    value = defaultValue;
end
end
