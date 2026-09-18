function y = applyTridiagonalBank(bank, x)
%APPLYTRIDIAGONALBANK  Multiply one COMPASS tridiagonal bank by a vector.
%
%   y = solver.applyTridiagonalBank(bank, x) evaluates y=A*x without
%   expanding A.  The bank contract is:
%     bank.matrix(:,1)  main diagonal, row i / column i
%     bank.matrix(:,2)  upper diagonal, row i / column i+1
%     bank.lower(i)     lower diagonal, row i+1 / column i
%
%   Inputs and output retain the caller's physical units.  This routine is
%   algebra only; it applies no boundary condition and changes no state.
%
%   Editor: Jizixing Zhu; Date: 2026-07-14.

    x = x(:);
    n = numel(x);
    validateBank(bank, n);

    y = bank.matrix(:, 1) .* x;
    if n > 1
        y(1:n-1) = y(1:n-1) ...
            + bank.matrix(1:n-1, 2) .* x(2:n);
        y(2:n) = y(2:n) + bank.lower(:) .* x(1:n-1);
    end
end

function validateBank(bank, n)
% Validate the shared main/upper/lower storage contract.
% Editor: Jizixing Zhu; Date: 2026-07-14.
    assert(isstruct(bank) && isfield(bank, 'matrix') && isfield(bank, 'lower') ...
        && isequal(size(bank.matrix), [n, 2]) ...
        && numel(bank.lower) == max(n - 1, 0) ...
        && all(isfinite(bank.matrix(:))) && all(isfinite(bank.lower(:))), ...
        'applyTridiagonalBank:badBank', ...
        'bank.matrix must be finite n-by-2 and bank.lower finite (n-1)-by-1.');
end
