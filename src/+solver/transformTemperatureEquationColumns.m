function transformed = transformTemperatureEquationColumns(row, U_h, U_T, U_P)
%TRANSFORMTEMPERATUREEQUATIONCOLUMNS  Change one FEM row from T to U coordinates.
%
%   transformed = solver.transformTemperatureEquationColumns(
%       row, U_h, U_T, U_P)
%   applies the nodal relation
%       dT = -(U_h/U_T) dh + (1/U_T) dU - (U_P/U_T) dP_g
%   to one already assembled correction equation
%       J_h dh + J_T dT + J_P dP_g = rhs.
%
%   row.h, row.T, and row.P are COMPASS tridiagonal banks.  The result has
%   row-coordinate banks transformed.h, transformed.U, transformed.P:
%       Jhat_h = J_h - J_T*diag(U_h/U_T)
%       Jhat_U = J_T*diag(1/U_T)
%       Jhat_P = J_P - J_T*diag(U_P/U_T).
%
%   U_h [J cm^-3 cm^-1], U_T [J cm^-3 K^-1], and
%   U_P [J cm^-3 pressure^-1] are nodal derivatives of the SAME internal
%   energy state used by the residual.  Non-positive U_T is a physical
%   closure failure and is never replaced by a numerical floor.
%
%   Reference: docs/compass_soil_water_heat_physics.md, Eqs. 5.28-5.30.
%   Editor: Jizixing Zhu; Date: 2026-07-14.

    U_h = U_h(:);
    U_T = U_T(:);
    U_P = U_P(:);
    n = numel(U_T);
    assert(numel(U_h) == n && numel(U_P) == n ...
        && all(isfinite(U_h)) && all(isfinite(U_T)) && all(isfinite(U_P)) ...
        && all(U_T > 0), ...
        'transformTemperatureEquationColumns:badEnergyDerivative', ...
        'U_h/U_T/U_P must be finite nodal vectors and U_T must be positive.');
    validateRow(row, n);

    a_h = U_h ./ U_T;
    a_U = 1 ./ U_T;
    a_P = U_P ./ U_T;

    T_h = solver.scaleTridiagonalColumns(row.T, a_h);
    T_U = solver.scaleTridiagonalColumns(row.T, a_U);
    T_P = solver.scaleTridiagonalColumns(row.T, a_P);

    transformed = struct();
    transformed.h = subtractBanks(row.h, T_h);
    transformed.U = T_U;
    transformed.P = subtractBanks(row.P, T_P);
    transformed.rhs = row.rhs(:);
    transformed.transform = struct('a_h', a_h, 'a_U', a_U, 'a_P', a_P);
end

function validateRow(row, n)
% Validate all three coordinate banks and the correction right-hand side.
% Editor: Jizixing Zhu; Date: 2026-07-14.
    assert(isstruct(row) && all(isfield(row, {'h','T','P','rhs'})) ...
        && numel(row.rhs) == n && all(isfinite(row.rhs(:))), ...
        'transformTemperatureEquationColumns:badRow', ...
        'row must contain finite h/T/P banks and an n-vector rhs.');
    fields = {'h','T','P'};
    for k = 1:numel(fields)
        bank = row.(fields{k});
        assert(isstruct(bank) && isfield(bank, 'matrix') && isfield(bank, 'lower') ...
            && isequal(size(bank.matrix), [n,2]) ...
            && numel(bank.lower) == max(n-1,0) ...
            && all(isfinite(bank.matrix(:))) && all(isfinite(bank.lower(:))), ...
            'transformTemperatureEquationColumns:badBank', ...
            'Every row bank must use the finite n-by-2/(n-1)-vector contract.');
    end
end

function out = subtractBanks(a, b)
% Subtract two banks with an identical sparsity contract.
% Editor: Jizixing Zhu; Date: 2026-07-14.
    out = struct('matrix', a.matrix - b.matrix, ...
        'lower', a.lower(:) - b.lower(:));
end
