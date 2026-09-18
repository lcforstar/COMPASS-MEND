function [correction, active, reaction] = solveLowerBoundCorrection( ...
        bank, rhs, current, initialActive, allowRelease)
%SOLVELOWERBOUNDCORRECTION Solve a tridiagonal correction with x>=0.
%
%   The unconstrained coordinate is x_new=current+correction. Nodes whose
%   unconstrained value is negative are placed on x_new=0. Their known
%   correction is eliminated from neighboring free rows, and the remaining
%   tridiagonal system is solved again with the repository Thomas kernel.
%   The active set grows monotonically within this frozen linearization.
%
%   This is the primal active-set step used when the dry-air component
%   disappears. It changes only the numerical coordinate solve; constitutive
%   dry-air storage and flux laws remain unchanged.
%
%   Editor: Jizixing Zhu; Date: 2026-08-12.
%   Change: 2026-08-13 Jizixing Zhu -- optionally hold an inherited active
%     set fixed until the enclosing nonlinear branch has converged.

    rhs = rhs(:);
    current = current(:);
    n = numel(rhs);
    validateBank(bank, n);
    assert(numel(current) == n && all(isfinite(rhs)) ...
            && all(isfinite(current)) && all(current >= 0), ...
        'solveLowerBoundCorrection:badInput', ...
        'rhs/current must be finite n-vectors and current must be nonnegative.');

    if nargin < 4 || isempty(initialActive)
        initialActive = false(n,1);
    else
        initialActive = logical(initialActive(:));
        assert(numel(initialActive) == n, ...
            'solveLowerBoundCorrection:badInitialActive', ...
            'initialActive must be empty or an n-vector logical mask.');
    end
    if nargin < 5 || isempty(allowRelease)
        allowRelease = true;
    end
    assert(isscalar(allowRelease), ...
        'solveLowerBoundCorrection:badAllowRelease', ...
        'allowRelease must be one scalar logical value.');
    allowRelease = logical(allowRelease);
    correction = solver.solveTridiagonal( ...
        bank.lower, bank.matrix(:,1), bank.matrix(:,2), rhs, [1,n]);
    scale = max(max(current), 1);
    active = initialActive | current + correction < -64 * eps(scale);
    reaction = zeros(n,1);
    for iteration = 1:2*n+2
        if ~any(active)
            break;
        end
        [reduced, reducedRhs] = imposeActiveValues( ...
            bank, rhs, active, -current);
        correction = solver.solveTridiagonal( ...
            reduced.lower, reduced.matrix(:,1), ...
            reduced.matrix(:,2), reducedRhs, [1,n]);
        reaction = solver.applyTridiagonalBank(bank, correction) - rhs;
        reactionScale = max(max(abs(reaction)), 1);
        release = allowRelease & active ...
            & reaction < -128 * eps(reactionScale);
        newlyActive = ~active ...
            & current + correction < -64 * eps(scale);
        if ~any(release) && ~any(newlyActive)
            break;
        end
        active = (active & ~release) | newlyActive;
        if ~any(active)
            correction = solver.solveTridiagonal( ...
                bank.lower, bank.matrix(:,1), ...
                bank.matrix(:,2), rhs, [1,n]);
            active = current + correction < -64 * eps(scale);
            if ~any(active)
                reaction = zeros(n,1);
                break;
            end
        end
    end
    correction(active) = -current(active);
    candidate = current + correction;
    assert(all(candidate >= -128 * eps(scale)), ...
        'solveLowerBoundCorrection:activeSetFailure', ...
        'The frozen active-set solve did not satisfy the nonnegative bound.');
    correction(candidate < 0) = -current(candidate < 0);
    reaction = solver.applyTridiagonalBank(bank, correction) - rhs;
end

function [reduced, reducedRhs] = imposeActiveValues( ...
        bank, rhs, active, fixedCorrection)
% Eliminate active columns from free rows and replace active rows by identity.
    reduced = bank;
    reducedRhs = rhs(:);
    n = numel(rhs);
    activeIndex = find(active(:)).';
    for i = activeIndex
        value = fixedCorrection(i);
        if i > 1 && ~active(i-1)
            reducedRhs(i-1) = reducedRhs(i-1) ...
                - bank.matrix(i-1,2) * value;
            reduced.matrix(i-1,2) = 0;
        end
        if i < n && ~active(i+1)
            reducedRhs(i+1) = reducedRhs(i+1) ...
                - bank.lower(i) * value;
            reduced.lower(i) = 0;
        end
    end
    for i = activeIndex
        reduced.matrix(i,1) = 1;
        reduced.matrix(i,2) = 0;
        reducedRhs(i) = fixedCorrection(i);
        if i > 1
            reduced.lower(i-1) = 0;
        end
    end
end

function validateBank(bank, n)
% Validate the shared compact tridiagonal contract.
    assert(isstruct(bank) && isfield(bank, 'matrix') ...
            && isfield(bank, 'lower') ...
            && isequal(size(bank.matrix), [n,2]) ...
            && numel(bank.lower) == max(n-1,0) ...
            && all(isfinite(bank.matrix(:))) ...
            && all(isfinite(bank.lower(:))), ...
        'solveLowerBoundCorrection:badBank', ...
        'bank must be a finite compact tridiagonal n-by-2/lower structure.');
end
