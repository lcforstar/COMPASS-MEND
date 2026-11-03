function [deltaH, deltaU, deltaM, active, reaction] = ...
        solveCoupledDryAirMassCorrection( ...
            waterRow, heatRow, airRow, rhsWater, rhsHeat, rhsAir, ...
            currentMass, initialActive)
%SOLVECOUPLEDDRYAIRMASSCORRECTION Solve the frozen h/U/Ma Picard block system.
%
%   The three transformed conservative rows are solved as one block-
%   tridiagonal system. Each nodal diagonal, lower, and upper block is 3-by-3
%   with equation order (water, heat, dry air) and unknown order (h, U, M_a).
%   A dry-air lower-bound node replaces only its air row by
%   deltaM=-currentMass. The primal-dual active set is then updated from the
%   original conservative air row.
%
%   This is an algebraic solve of the already assembled COMPASS Jacobian. It
%   changes neither equations nor signs, and it avoids the non-contractive
%   scalar h->M_a->U->h Gauss-Seidel cycle at component disappearance.
%
%   Units are inherited from the row contracts: h [cm], U [J cm^-2], and
%   nodal M_a [g cm^-2]. The coupled matrix is equilibrated in equation and
%   coordinate units before partial-pivoting solve; no accepted state is
%   mutated here.
%
%   Editor: Jizixing Zhu; Date: 2026-08-13.

    rhsWater = rhsWater(:);
    rhsHeat = rhsHeat(:);
    rhsAir = rhsAir(:);
    currentMass = currentMass(:);
    n = numel(currentMass);
    assert(numel(rhsWater) == n && numel(rhsHeat) == n ...
            && numel(rhsAir) == n && all(isfinite(currentMass)) ...
            && all(currentMass >= 0), ...
        'solveCoupledDryAirMassCorrection:badInput', ...
        'All right-hand sides and currentMass must be compatible n-vectors.');
    validateRow(waterRow, n);
    validateRow(heatRow, n);
    validateRow(airRow, n);
    if nargin < 8 || isempty(initialActive)
        active = false(n,1);
    else
        active = logical(initialActive(:));
        assert(numel(active) == n, ...
            'solveCoupledDryAirMassCorrection:badInitialActive', ...
            'initialActive must contain n logical values.');
    end

    original = buildBlocks(waterRow, heatRow, airRow, ...
        rhsWater, rhsHeat, rhsAir);
    scale = max(max(currentMass), 1);
    reaction = zeros(n,1);
    for iteration = 1:2*n+2
        system = imposeMassBound(original, active, currentMass);
        correction = solveEquilibrated(system);
        deltaH = correction(1,:).';
        deltaU = correction(2,:).';
        deltaM = correction(3,:).';
        reaction = applyAirRow(airRow, deltaH, deltaU, deltaM) - rhsAir;
        reactionScale = max(max(abs(reaction)), 1);
        release = active & reaction < -128 * eps(reactionScale);
        newlyActive = ~active ...
            & currentMass + deltaM < -64 * eps(scale);
        if ~any(release) && ~any(newlyActive)
            break;
        end
        active = (active & ~release) | newlyActive;
    end
    deltaM(active) = -currentMass(active);
    assert(all(currentMass + deltaM >= -128 * eps(scale)), ...
        'solveCoupledDryAirMassCorrection:activeSetFailure', ...
        'The block active-set solve did not preserve nonnegative dry-air mass.');
end

function system = buildBlocks(water, heat, air, rhsW, rhsH, rhsA)
    n = numel(rhsW);
    system.diagonal = zeros(3,3,n);
    system.upper = zeros(3,3,max(n-1,0));
    system.lower = zeros(3,3,max(n-1,0));
    system.rhs = [rhsW.'; rhsH.'; rhsA.'];
    rows = {water, heat, air};
    names = {'h','U','M'};
    for equation = 1:3
        for variable = 1:3
            bank = rows{equation}.(names{variable});
            system.diagonal(equation,variable,:) = bank.matrix(:,1);
            if n > 1
                system.upper(equation,variable,:) = bank.matrix(1:n-1,2);
                system.lower(equation,variable,:) = bank.lower(:);
            end
        end
    end
end

function system = imposeMassBound(original, active, currentMass)
    system = original;
    n = numel(active);
    for node = find(active(:)).'
        system.diagonal(3,:,node) = 0;
        system.diagonal(3,3,node) = 1;
        system.rhs(3,node) = -currentMass(node);
        if node > 1
            system.lower(3,:,node-1) = 0;
        end
        if node < n
            system.upper(3,:,node) = 0;
        end
    end
end

function solution = solveEquilibrated(system)
    n = size(system.diagonal,3);
    matrix = spalloc(3*n, 3*n, 27*n);
    for node = 1:n
        index = (3*node-2):(3*node);
        matrix(index,index) = system.diagonal(:,:,node);
        if node < n
            next = index + 3;
            matrix(index,next) = system.upper(:,:,node);
            matrix(next,index) = system.lower(:,:,node);
        end
    end
    rhs = system.rhs(:);
    columnNorm = full(max(abs(matrix), [], 1)).';
    assert(all(isfinite(columnNorm)) && all(columnNorm > 0), ...
        'solveCoupledDryAirMassCorrection:zeroColumn', ...
        'Every coupled correction coordinate must have a finite nonzero column.');
    coordinateScale = 1 ./ columnNorm;
    scaled = matrix * spdiags(coordinateScale, 0, 3*n, 3*n);
    rowNorm = full(max(abs(scaled), [], 2));
    assert(all(isfinite(rowNorm)) && all(rowNorm > 0), ...
        'solveCoupledDryAirMassCorrection:zeroRow', ...
        'Every coupled correction equation must have a finite nonzero row.');
    scaled = spdiags(1 ./ rowNorm, 0, 3*n, 3*n) * scaled;
    scaledRhs = rhs ./ rowNorm;
    reciprocalCondition = rcond(full(scaled));
    assert(isfinite(reciprocalCondition) && reciprocalCondition > 1e-14, ...
        'solveCoupledDryAirMassCorrection:singularSystem', ...
        'The equilibrated coupled correction matrix is singular (rcond %.3e).', ...
        reciprocalCondition);
    raw = coordinateScale .* (scaled \ scaledRhs);
    solution = reshape(raw, 3, n);
end

function value = applyAirRow(row, deltaH, deltaU, deltaM)
    value = solver.applyTridiagonalBank(row.h, deltaH) ...
        + solver.applyTridiagonalBank(row.U, deltaU) ...
        + solver.applyTridiagonalBank(row.M, deltaM);
end

function validateRow(row, n)
    assert(isstruct(row) && all(isfield(row, {'h','U','M','rhs'})), ...
        'solveCoupledDryAirMassCorrection:badRow', ...
        'Every equation row must contain h/U/M banks and rhs.');
    for name = {'h','U','M'}
        bank = row.(name{1});
        assert(isstruct(bank) && isfield(bank, 'matrix') ...
                && isfield(bank, 'lower') ...
                && isequal(size(bank.matrix), [n,2]) ...
                && numel(bank.lower) == max(n-1,0) ...
                && all(isfinite(bank.matrix(:))) ...
                && all(isfinite(bank.lower(:))), ...
            'solveCoupledDryAirMassCorrection:badBank', ...
            'Every block bank must use the finite compact tridiagonal contract.');
    end
end
