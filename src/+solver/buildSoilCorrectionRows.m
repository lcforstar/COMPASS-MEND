function rows = buildSoilCorrectionRows(matrices, systems, state, dt_s, heatEnabled, airEnabled)
%BUILDSOILCORRECTIONROWS  Form complete h/T/P_g correction equations.
%
%   rows = solver.buildSoilCorrectionRows(matrices, systems, state, dt_s,
%       heatEnabled, airEnabled)
%   converts the three post-boundary absolute scalar systems into residual
%   equations at the current Picard iterate:
%       J_h dh + J_T dT + J_P dP_g = b - A_self*x_self.
%
%   The nonlinear storage residual already present in each absolute RHS is
%   retained exactly.  The complete local storage tangent is then restored
%   to every active correction column:
%     water:  C_hh, C_hT
%     dry air: C_ah, C_aT, C_aa
%     energy: C_Th, C_TT, C_Ta.
%   This is not a second storage increment: the exact accepted-to-current
%   storage difference defines the residual, while C/dt linearizes only the
%   next correction.
%
%   Dirichlet rows are detected from the already applied scalar system.
%   Cross-coordinate banks are zeroed on water-head and pressure Dirichlet
%   rows.  On temperature Dirichlet rows the T bank remains the identity;
%   after the T-to-U transform it becomes
%       -a_h dh + a_U dU - a_P dP_g = T_bc - T_current.
%
%   Inputs use the existing COMPASS cm-g-s equation units.  This function
%   performs algebra only and mutates no model state.
%
%   Reference: docs/compass_soil_water_heat_physics.md Sec. 5.9.1-5.9.5.
%   Editor: Jizixing Zhu; Date: 2026-07-14.
%   Change: 2026-07-17 Jizixing Zhu -- retain the strict-lower gas-pressure
%           gradient coefficient in the water correction row.
%   Change: 2026-08-12 Jizixing Zhu -- add Vvh/VvT only to correction
%           Jacobian columns while evaluating the physical residual from
%           the conservative Kha pressure-flux bank exactly once.

    assert(isscalar(dt_s) && isfinite(dt_s) && dt_s > 0, ...
        'buildSoilCorrectionRows:badTimeStep', ...
        'dt_s must be a finite positive scalar.');
    h = state.h(:);
    T = state.T(:);
    P = state.P(:);
    n = numel(h);
    assert(numel(T) == n && numel(P) == n ...
        && all(isfinite(h)) && all(isfinite(T)) && all(isfinite(P)), ...
        'buildSoilCorrectionRows:badState', ...
        'h, T, and P must be finite nodal vectors of equal length.');

    rows = struct();

    % Water row: the post-boundary self bank already contains C1/dt.
    waterResidualH = bankFromSystem(systems.water.matrix, 'C4', 'C4_a', n);
    waterVaporH = optionalBankFromMatrix( ...
        matrices.water, 'C4_vapor_jacobian', ...
        'C4_vapor_jacobian_a', n);
    waterT = addLumpedStorage( ...
        bankFromMatrix(matrices.water, 'C5', 'C5_a', n), ...
        matrices.water.C2, dt_s);
    waterT = addBanks(waterT, optionalBankFromMatrix( ...
        matrices.water, 'C5_vapor_jacobian', ...
        'C5_vapor_jacobian_a', n));
    waterP = bankFromSymmetricGradient( ...
        matrices.water, 'C6', n);
    waterDirichlet = boundaryDirichletRows(waterResidualH);
    waterVaporH = zeroRows(waterVaporH, waterDirichlet);
    waterH = addBanks(waterResidualH, waterVaporH);
    waterT = zeroRows(waterT, waterDirichlet);
    waterP = zeroRows(waterP, waterDirichlet);
    rows.water = struct('h', waterH, 'T', waterT, 'P', waterP, ...
        'rhs', systems.water.rhs(:) ...
            - solver.applyTridiagonalBank(waterResidualH, h), ...
        'dirichletRows', waterDirichlet);

    if airEnabled
        airSelf = bankFromSystem(systems.air.matrix, 'C6', 'C6_a', n);
        airH = addLumpedStorage( ...
            bankFromMatrix(matrices.dry_air, 'C4', 'C4_a', n), ...
            matrices.dry_air.C1, dt_s);
        airT = addLumpedStorage( ...
            bankFromMatrix(matrices.dry_air, 'C5', 'C5_a', n), ...
            matrices.dry_air.C2, dt_s);
        airDirichlet = boundaryDirichletRows(airSelf);
        airH = zeroRows(airH, airDirichlet);
        airT = zeroRows(airT, airDirichlet);
        rows.air = struct('h', airH, 'T', airT, 'P', airSelf, ...
            'rhs', systems.air.rhs(:) ...
                - solver.applyTridiagonalBank(airSelf, P), ...
            'dirichletRows', airDirichlet);
    else
        rows.air = struct();
    end

    if heatEnabled
        heatSelf = bankFromSystem(systems.heat.matrix, 'C5', 'C5_a', n);
        heatH = addLumpedStorage( ...
            bankFromMatrix(matrices.heat, 'C4', 'C4_a', n), ...
            matrices.heat.C1, dt_s);
        heatP = addLumpedStorage( ...
            bankFromMatrix(matrices.heat, 'C6', 'C6_a', n), ...
            matrices.heat.C3, dt_s);
        heatDirichlet = boundaryDirichletRows(heatSelf);
        heatH = zeroRows(heatH, heatDirichlet);
        heatP = zeroRows(heatP, heatDirichlet);
        rows.heat = struct('h', heatH, 'T', heatSelf, 'P', heatP, ...
            'rhs', systems.heat.rhs(:) ...
                - solver.applyTridiagonalBank(heatSelf, T), ...
            'dirichletRows', heatDirichlet);
    else
        rows.heat = struct();
    end
end

function bank = optionalBankFromMatrix(matrixStruct, matrixField, lowerField, n)
% Read an optional correction-only bank, defaulting to an exact zero bank.
% Editor: Jizixing Zhu; Date: 2026-08-12.
    if ~isfield(matrixStruct, matrixField)
        bank = struct('matrix', zeros(n, 2), ...
            'lower', zeros(max(n - 1, 0), 1));
        return;
    end
    bank = bankFromMatrix(matrixStruct, matrixField, lowerField, n);
end

function out = addBanks(a, b)
% Add two tridiagonal banks with identical nodal contracts.
% Editor: Jizixing Zhu; Date: 2026-08-12.
    assert(isequal(size(a.matrix), size(b.matrix)) ...
            && isequal(size(a.lower), size(b.lower)), ...
        'buildSoilCorrectionRows:bankShapeMismatch', ...
        'Tridiagonal banks must have identical shapes before addition.');
    out = struct('matrix', a.matrix + b.matrix, ...
        'lower', a.lower + b.lower);
end

function bank = bankFromSymmetricGradient(matrixStruct, matrixField, n)
% Convert a symmetric element-gradient block whose lower equals prior upper.
% Editor: Jizixing Zhu; Date: 2026-07-17.
    bank = bankFromMatrix(matrixStruct, matrixField, '', n);
    if n > 1
        % Water C6 is assembled from Kha*dP_g/dz. Its element matrix is
        % symmetric, so row i+1's strict-lower entry is row i's stored upper.
        bank.lower = matrixStruct.(matrixField)(1:n-1, 2);
    end
end

function bank = bankFromSystem(matrixStruct, matrixField, lowerField, n)
% Read a post-boundary self-coordinate bank.
% Editor: Jizixing Zhu; Date: 2026-07-14.
    bank = bankFromMatrix(matrixStruct, matrixField, lowerField, n);
end

function bank = bankFromMatrix(matrixStruct, matrixField, lowerField, n)
% Convert the equation matrix naming contract to the shared bank contract.
% Editor: Jizixing Zhu; Date: 2026-07-14.
    assert(isfield(matrixStruct, matrixField), ...
        'buildSoilCorrectionRows:missingMatrix', ...
        'Missing matrix field %s.', matrixField);
    matrix = matrixStruct.(matrixField);
    assert(isequal(size(matrix), [n,2]) && all(isfinite(matrix(:))), ...
        'buildSoilCorrectionRows:badMatrix', ...
        '%s must be a finite n-by-2 matrix.', matrixField);
    if isempty(lowerField)
        lower = zeros(max(n-1,0), 1);
    else
        assert(isfield(matrixStruct, lowerField), ...
            'buildSoilCorrectionRows:missingLower', ...
            'Missing lower-diagonal field %s.', lowerField);
        rawLower = matrixStruct.(lowerField);
        assert(numel(rawLower) >= max(n-1,0) ...
            && all(isfinite(rawLower(1:max(n-1,0)))), ...
            'buildSoilCorrectionRows:badLower', ...
            '%s must contain at least n-1 finite entries.', lowerField);
        lower = rawLower(1:max(n-1,0));
    end
    matrix(n,2) = 0;
    bank = struct('matrix', matrix, 'lower', lower(:));
end

function bank = addLumpedStorage(bank, storage, dt_s)
% Add a lumped endpoint storage matrix to one spatial bank.
% Editor: Jizixing Zhu; Date: 2026-07-14.
    n = size(bank.matrix, 1);
    assert(isequal(size(storage), [n,2]) && all(isfinite(storage(:))) ...
        && all(storage(:,2) == 0), ...
        'buildSoilCorrectionRows:notLumped', ...
        'Correction-coordinate storage matrices must be finite and lumped.');
    bank.matrix(:,1) = bank.matrix(:,1) + storage(:,1) / dt_s;
end

function mask = boundaryDirichletRows(bank)
% Identify identity rows after scalar Dirichlet elimination.
%
% The ordinary soil column has only endpoint Dirichlet conditions. An
% external MODFLOW water table is an interior bottom-up node, while rows
% beneath it are not Richards unknowns. Detecting every isolated identity
% row keeps correction coordinates consistent with the scalar system.
% Editor: Jizixing Zhu; Date: 2026-07-14.
% Change: 2026-07-20 Jizixing Zhu -- recognize the interior MODFLOW fixed-head row.
    n = size(bank.matrix, 1);
    mask = false(n,1);
    scale = max(1, max(abs(bank.matrix(:))));
    tol = 100 * eps(scale);
    for row = 1:n
        upperIsZero = row == n || abs(bank.matrix(row,2)) <= tol;
        lowerIsZero = row == 1 || abs(bank.lower(row-1)) <= tol;
        mask(row) = abs(bank.matrix(row,1)-1) <= tol ...
            && upperIsZero && lowerIsZero;
    end
end

function bank = zeroRows(bank, mask)
% Remove cross-coordinate contributions from prescribed self rows.
% Editor: Jizixing Zhu; Date: 2026-07-14.
    n = numel(mask);
    ids = find(mask(:));
    bank.matrix(ids,1) = 0;
    bank.matrix(ids,2) = 0;
    lowerRows = ids(ids >= 2) - 1;
    lowerRows = lowerRows(lowerRows >= 1 & lowerRows <= n-1);
    bank.lower(lowerRows) = 0;
end
