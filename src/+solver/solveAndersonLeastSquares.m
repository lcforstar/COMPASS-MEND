function [gamma, info] = solveAndersonLeastSquares(dResidual, residual, regularization, minSingularRatio)
%SOLVEANDERSONLEASTSQUARES  Solve the small Anderson problem without normal equations.
%
%   [gamma, info] = solver.solveAndersonLeastSquares(dResidual, residual,
%   regularization, minSingularRatio) checks the singular spectrum of the
%   residual-difference history, then solves the Tikhonov problem through an
%   augmented least-squares system. It never forms dResidual' * dResidual.
%
%   Inputs:
%     dResidual       [in]  [-] scaled residual-difference matrix, nState x nHistory
%     residual        [in]  [-] current scaled residual, nState x 1
%     regularization  [in]  [-] nonnegative Tikhonov coefficient
%     minSingularRatio[in]  [-] minimum sigma_min / sigma_max for accepting history
%
%   Returns:
%     gamma [out] [-] Anderson coefficients; empty when history is rejected
%     info  [out] struct with accepted, reason, singularRatio, sigmaMin, and sigmaMax
%
%   Closure:
%     Accepted history solves min ||D*gamma-r||^2 + lambda||gamma||^2. Rank-deficient or nearly
%     dependent history is rejected before division so the caller can clear history and use the
%     ordinary Picard iterate without a near-singular-matrix warning.
%
%   Editor: Jizixing Zhu; Date: 2026-07-14.
    if nargin < 3 || isempty(regularization), regularization = 0; end
    if nargin < 4 || isempty(minSingularRatio), minSingularRatio = 1e-8; end

    gamma = [];
    info = struct('accepted', false, 'reason', 'not_evaluated', ...
        'singularRatio', 0, 'sigmaMin', 0, 'sigmaMax', 0);

    residual = residual(:);
    nHistory = size(dResidual, 2);
    if isempty(dResidual) || size(dResidual, 1) ~= numel(residual) || ...
            nHistory == 0 || any(~isfinite(dResidual(:))) || any(~isfinite(residual))
        info.reason = 'invalid_input';
        return;
    end
    if regularization < 0 || ~isfinite(regularization) || ...
            minSingularRatio < 0 || ~isfinite(minSingularRatio)
        info.reason = 'invalid_option';
        return;
    end

    singularValues = svd(dResidual, 'econ');
    if numel(singularValues) < nHistory || isempty(singularValues)
        info.reason = 'rank_deficient';
        return;
    end
    info.sigmaMax = singularValues(1);
    info.sigmaMin = singularValues(end);
    if info.sigmaMax > 0
        info.singularRatio = info.sigmaMin / info.sigmaMax;
    end
    if info.sigmaMax <= eps(class(info.sigmaMax)) || info.singularRatio < minSingularRatio
        info.reason = 'ill_conditioned';
        return;
    end

    if regularization > 0
        augmentedMatrix = [dResidual; sqrt(regularization) * eye(nHistory)];
        augmentedResidual = [residual; zeros(nHistory, 1)];
        gamma = augmentedMatrix \ augmentedResidual;
    else
        gamma = dResidual \ residual;
    end

    if any(~isfinite(gamma))
        gamma = [];
        info.reason = 'nonfinite_solution';
        return;
    end
    info.accepted = true;
    info.reason = 'accepted';
end
