function scaled = scaleTridiagonalColumns(bank, columnScale)
%SCALETRIDIAGONALCOLUMNS  Right-scale a tridiagonal matrix by a diagonal.
%
%   scaled = solver.scaleTridiagonalColumns(bank, s) returns the compact
%   bank for A*diag(s).  Entries are scaled by their COLUMN node:
%     lower(i)      <- lower(i)      * s(i)
%     diagonal(i)   <- diagonal(i)   * s(i)
%     upper(i)      <- upper(i)      * s(i+1)
%
%   This is the exact compact operation required by the matrix-level T-to-U
%   coordinate transform in docs/compass_soil_water_heat_physics.md Sec. 5.9.
%   It preserves the tridiagonal sparsity pattern and performs no state or
%   boundary mutation.
%
%   Editor: Jizixing Zhu; Date: 2026-07-14.

    columnScale = columnScale(:);
    n = numel(columnScale);
    assert(isstruct(bank) && isfield(bank, 'matrix') && isfield(bank, 'lower') ...
        && isequal(size(bank.matrix), [n, 2]) ...
        && numel(bank.lower) == max(n - 1, 0) ...
        && all(isfinite(bank.matrix(:))) && all(isfinite(bank.lower(:))) ...
        && all(isfinite(columnScale)), ...
        'scaleTridiagonalColumns:badInput', ...
        'The bank and columnScale must be finite and dimensionally consistent.');

    scaled = bank;
    scaled.matrix(:, 1) = bank.matrix(:, 1) .* columnScale;
    if n > 1
        scaled.matrix(1:n-1, 2) = bank.matrix(1:n-1, 2) .* columnScale(2:n);
        scaled.matrix(n, 2) = 0;
        scaled.lower = bank.lower(:) .* columnScale(1:n-1);
    else
        scaled.matrix(1, 2) = 0;
        scaled.lower = zeros(0, 1);
    end
end
