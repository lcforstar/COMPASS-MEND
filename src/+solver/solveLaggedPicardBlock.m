function [correction, info] = solveLaggedPicardBlock(linearization)
%SOLVELAGGEDPICARDBLOCK Solve one frozen h/U/P_g block-tridiagonal system.
%
%   This is a block Thomas/Schur elimination of the already assembled
%   Modified-Picard rows. Constitutive coefficients remain frozen; no
%   derivative is refreshed inside the solve and no Newton iteration occurs.
%   The input rows use equation order water/energy/dry_air and coordinate
%   order h [cm], U [J cm-3], and P_g [internal pressure]. The returned
%   correction has the same nodal coordinates; INFO reports pivot conditioning
%   and the linear residual. Singular or ill-conditioned pivots fail the trial.
%
%   Editor: Jizixing Zhu; Date: 2026-08-14.
    assert(isstruct(linearization) && linearization.available, ...
        'solveLaggedPicardBlock:missingLinearization', ...
        'A complete frozen h/U/P_g linearization is required.');
    rows = linearization.rows;
    equationNames = {'water','energy','dry_air'};
    coordinateNames = {'h','U','P'};
    n = numel(rows.water.rhs);
    diagonal = zeros(3,3,n);
    upper = zeros(3,3,max(n-1,0));
    lower = zeros(3,3,max(n-1,0));
    rhs = zeros(3,n);
    for equation = 1:3
        row = rows.(equationNames{equation});
        rhs(equation,:) = row.rhs(:).';
        for coordinate = 1:3
            bank = row.(coordinateNames{coordinate});
            diagonal(equation,coordinate,:) = bank.matrix(:,1);
            if n > 1
                upper(equation,coordinate,:) = bank.matrix(1:n-1,2);
                lower(equation,coordinate,:) = bank.lower(:);
            end
        end
    end

    reducedDiagonal = diagonal;
    reducedRhs = rhs;
    minimumRcond = Inf;
    for node = 2:n
        prior = reducedDiagonal(:,:,node-1);
        condition = rcond(prior);
        minimumRcond = min(minimumRcond, condition);
        assert(isfinite(condition) && condition > 1e-14, ...
            'solveLaggedPicardBlock:singularBlockPivot', ...
            'Lagged Picard block pivot %d is singular or ill-conditioned.', ...
            node - 1);
        multiplier = lower(:,:,node-1) / prior;
        reducedDiagonal(:,:,node) = diagonal(:,:,node) ...
            - multiplier * upper(:,:,node-1);
        reducedRhs(:,node) = rhs(:,node) ...
            - multiplier * reducedRhs(:,node-1);
    end
    lastCondition = rcond(reducedDiagonal(:,:,n));
    minimumRcond = min(minimumRcond, lastCondition);
    assert(isfinite(lastCondition) && lastCondition > 1e-14, ...
        'solveLaggedPicardBlock:singularBlockPivot', ...
        'Final lagged Picard block pivot is singular or ill-conditioned.');

    solution = zeros(3,n);
    solution(:,n) = reducedDiagonal(:,:,n) \ reducedRhs(:,n);
    for node = n-1:-1:1
        solution(:,node) = reducedDiagonal(:,:,node) ...
            \ (reducedRhs(:,node) - upper(:,:,node) ...
                * solution(:,node+1));
    end
    linearResidual = applyBlock(diagonal, upper, lower, solution) - rhs;
    correction = struct('h', solution(1,:).', ...
        'U', solution(2,:).', 'P', solution(3,:).');
    info = struct('minimum_block_rcond', minimumRcond, ...
        'max_abs_linear_residual', max(abs(linearResidual(:))), ...
        'relative_linear_residual', norm(linearResidual(:),2) ...
            / max(norm(rhs(:),2), realmin), ...
        'newton_used', false, 'method', 'block_thomas_schur');
end

function value = applyBlock(diagonal, upper, lower, x)
    n = size(x,2);
    value = zeros(3,n);
    for node = 1:n
        value(:,node) = diagonal(:,:,node) * x(:,node);
        if node < n
            value(:,node) = value(:,node) ...
                + upper(:,:,node) * x(:,node+1);
        end
        if node > 1
            value(:,node) = value(:,node) ...
                + lower(:,:,node-1) * x(:,node-1);
        end
    end
end
