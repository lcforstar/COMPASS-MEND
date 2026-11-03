function [correction, info] = solveLaggedPicardMassBlock( ...
        rows, currentMass, initialActive)
%SOLVELAGGEDPICARDMASSBLOCK Solve one frozen h/U/M_a Picard active set.
%
% The constitutive and transport coefficients are already frozen by the
% enclosing Modified-Picard epoch. This function only solves that linear
% block while enforcing M_a >= 0. Updating the bound-active set is an
% algebraic complementarity solve; it does not refresh coefficients, form a
% nonlinear Jacobian, or perform a Newton iteration.
% ROWS use water/energy/dry_air equation order and h [cm], U [J cm-3],
% conservative nodal M_a [g cm-2] coordinates. CURRENTMASS and the returned
% M correction share [g cm-2]. Singular systems or an infeasible bound-active
% set fail the trial; this function never installs or commits a model state.
%
% Editor: Jizixing Zhu; Date: 2026-08-15.
    currentMass = currentMass(:);
    n = numel(currentMass);
    assert(isstruct(rows) && all(isfield(rows, ...
            {'water','energy','dry_air'})) ...
            && all(isfinite(currentMass)) && all(currentMass >= 0), ...
        'solveLaggedPicardMassBlock:badInput', ...
        'Complete h/U/M rows and nonnegative current mass are required.');
    if nargin < 3 || isempty(initialActive)
        active = false(n,1);
    else
        active = logical(initialActive(:));
        assert(numel(active) == n, ...
            'solveLaggedPicardMassBlock:badActiveSet', ...
            'The initial active set must contain one flag per node.');
    end
    original = buildSystem(rows, n);
    reaction = zeros(n,1);
    minimumRcond = Inf;
    for iteration = 1:(2*n + 2)
        constrained = imposeMassBound(original, active, currentMass);
        [solution, condition] = solveEquilibrated(constrained);
        minimumRcond = min(minimumRcond, condition);
        deltaH = solution(1,:).';
        deltaU = solution(2,:).';
        deltaM = solution(3,:).';
        reaction = applyAirRow(rows.dry_air, ...
            deltaH, deltaU, deltaM) - rows.dry_air.rhs(:);
        reactionScale = max(max(abs(reaction)), 1);
        massScale = max(max(currentMass), 1);
        release = active & reaction < -128 * eps(reactionScale);
        newlyActive = ~active ...
            & currentMass + deltaM < -64 * eps(massScale);
        if ~any(release) && ~any(newlyActive)
            break;
        end
        active = (active & ~release) | newlyActive;
    end
    deltaM(active) = -currentMass(active);
    massScale = max(max(currentMass), 1);
    assert(all(currentMass + deltaM >= -128 * eps(massScale)), ...
        'solveLaggedPicardMassBlock:activeSetFailure', ...
        'The frozen Picard active set violated M_a >= 0.');
    residual = applyFullSystem(original, [deltaH.'; deltaU.'; deltaM.']) ...
        - original.rhs;
    residual(3,active) = reaction(active);
    correction = struct('h', deltaH, 'U', deltaU, 'M', deltaM, ...
        'active', active, 'reaction', reaction);
    info = struct('minimum_block_rcond', minimumRcond, ...
        'max_abs_linear_residual', max(abs(residual(~massResidualMask(active)))), ...
        'relative_linear_residual', ...
            norm(residual(~massResidualMask(active)),2) ...
            / max(norm(original.rhs(:),2), realmin), ...
        'active_node_count', nnz(active), ...
        'newton_used', false, ...
        'method', 'equilibrated_block_picard_active_set');
end

function system = buildSystem(rows, n)
    system.diagonal = zeros(3,3,n);
    system.upper = zeros(3,3,max(n-1,0));
    system.lower = zeros(3,3,max(n-1,0));
    system.rhs = [rows.water.rhs(:).'; rows.energy.rhs(:).'; ...
        rows.dry_air.rhs(:).'];
    equations = {'water','energy','dry_air'};
    coordinates = {'h','U','M'};
    for equation = 1:3
        row = rows.(equations{equation});
        for coordinate = 1:3
            bank = row.(coordinates{coordinate});
            system.diagonal(equation,coordinate,:) = bank.matrix(:,1);
            if n > 1
                system.upper(equation,coordinate,:) = ...
                    bank.matrix(1:n-1,2);
                system.lower(equation,coordinate,:) = bank.lower(:);
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

function [solution, reciprocalCondition] = solveEquilibrated(system)
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
        'solveLaggedPicardMassBlock:zeroColumn', ...
        'Every h/U/M_a coordinate must have a finite nonzero column.');
    coordinateScale = 1 ./ columnNorm;
    scaled = matrix * spdiags(coordinateScale, 0, 3*n, 3*n);
    rowNorm = full(max(abs(scaled), [], 2));
    assert(all(isfinite(rowNorm)) && all(rowNorm > 0), ...
        'solveLaggedPicardMassBlock:zeroRow', ...
        'Every frozen Picard equation must have a finite nonzero row.');
    scaled = spdiags(1 ./ rowNorm, 0, 3*n, 3*n) * scaled;
    scaledRhs = rhs ./ rowNorm;
    reciprocalCondition = rcond(full(scaled));
    assert(isfinite(reciprocalCondition) && reciprocalCondition > 1e-14, ...
        'solveLaggedPicardMassBlock:singularSystem', ...
        'The equilibrated h/U/M_a block is singular (rcond %.3e).', ...
        reciprocalCondition);
    raw = coordinateScale .* (scaled \ scaledRhs);
    solution = reshape(raw, 3, n);
end

function value = applyAirRow(row, deltaH, deltaU, deltaM)
    value = solver.applyTridiagonalBank(row.h, deltaH) ...
        + solver.applyTridiagonalBank(row.U, deltaU) ...
        + solver.applyTridiagonalBank(row.M, deltaM);
end

function value = applyFullSystem(system, x)
    n = size(x,2);
    value = zeros(3,n);
    for node = 1:n
        value(:,node) = system.diagonal(:,:,node) * x(:,node);
        if node < n
            value(:,node) = value(:,node) ...
                + system.upper(:,:,node) * x(:,node+1);
        end
        if node > 1
            value(:,node) = value(:,node) ...
                + system.lower(:,:,node-1) * x(:,node-1);
        end
    end
end

function mask = massResidualMask(active)
    mask = false(3,numel(active));
    mask(3,active) = true;
end
