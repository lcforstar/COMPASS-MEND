function [x, diagonalEliminated] = solveTridiagonal(lower, diagonal, upper, rhs, activeRange)
%SOLVETRIDIAGONAL  Solve a scalar tridiagonal system with a guarded Thomas sweep.
%
%   [x, diagonalEliminated] = solver.solveTridiagonal(lower, diagonal,
%       upper, rhs, activeRange) solves A*x=rhs over the inclusive node range
%   [first,last]. Values outside that range are returned unchanged.
%
%   Inputs:
%     lower      [coefficient, N-1 x 1] sub-diagonal A(i+1,i)
%     diagonal   [coefficient, N x 1] main diagonal A(i,i)
%     upper      [coefficient, N x 1] super-diagonal A(i,i+1); upper(N) unused
%     rhs        [equation unit, N x 1] right-hand side
%     activeRange[-, 1 x 2] inclusive [first,last] solve range
%
%   Returns:
%     x                  [unknown unit, N x 1] solution
%     diagonalEliminated [coefficient, N x 1] forward-eliminated diagonal
%
%   Method:
%     Pure scalar Thomas elimination with a sign-preserving 1e-20 pivot
%     guard. It has no knowledge of water head, temperature, gas pressure,
%     boundary conditions, damping, or accepted-state transactions.
%
%   Reference:
%     soil.water.solveTridiagonalMatrixEquations
%     soil.thermal.solveTridiagonalMatrixEquations
%     soil.gas.solveTridiagonalMatrixEquations
%   Editor: Jizixing Zhu; Date: 2026-07-13.

    diagonalEliminated = diagonal(:);
    rhsWork = rhs(:);
    upper = upper(:);
    lower = lower(:);
    N = numel(diagonalEliminated);
    assert(numel(rhsWork) == N && numel(upper) >= N - 1 && numel(lower) >= N - 1, ...
        'solveTridiagonal:badShape', 'Tridiagonal arrays have incompatible lengths.');
    if nargin < 5 || isempty(activeRange)
        activeRange = [1, N];
    end
    assert(numel(activeRange) == 2 && activeRange(1) >= 1 ...
        && activeRange(2) <= N && activeRange(1) <= activeRange(2), ...
        'solveTridiagonal:badActiveRange', 'activeRange must be [first,last] within 1:N.');

    first = activeRange(1);
    last = activeRange(2);
    pivotFloor = 1e-20;
    diagonalEliminated(first) = guardPivot(diagonalEliminated(first), pivotFloor);
    rhsWork(first) = rhsWork(first) / diagonalEliminated(first);

    for node = first + 1:last
        diagonalEliminated(node) = diagonalEliminated(node) ...
            - lower(node - 1) * upper(node - 1) / diagonalEliminated(node - 1);
        diagonalEliminated(node) = guardPivot(diagonalEliminated(node), pivotFloor);
        rhsWork(node) = (rhsWork(node) - lower(node - 1) * rhsWork(node - 1)) ...
            / diagonalEliminated(node);
    end
    for node = last - 1:-1:first
        rhsWork(node) = rhsWork(node) ...
            - upper(node) * rhsWork(node + 1) / diagonalEliminated(node);
    end
    x = rhsWork;
end

function value = guardPivot(value, floorValue)
% Preserve pivot sign while preventing division by numerical zero.
    pivotSign = sign(value);
    if pivotSign == 0
        pivotSign = 1;
    end
    value = max(abs(value), floorValue) * pivotSign;
end
