function [T, result] = invertTemperatureFromInternalEnergy( ...
        U_target, h, P_g, T_guess, refreshState, dz, L_f, opts)
%INVERTTEMPERATUREFROMINTERNALENERGY  Vectorized active-set Newton T ← U.
%
%   [T, result] = soil.thermal.invertTemperatureFromInternalEnergy(...
%       U_target, h, P_g, T_guess, refreshState, dz, L_f)
%   inverts the nodal internal-energy state function
%       U_target(i) = U_nodal(h(i), T(i), P_g(i))
%   for T at every node simultaneously, using protected Newton with
%   bisection fallback.  The function refreshes the full-column SFCC and
%   thermodynamic state ONCE per iteration (NOT per node).
%
%   Inputs:
%     U_target    [J/cm³]    target nodal internal energy (NN×1)
%     h           [cm]       current matric head at nodes (NN×1)
%     P_g         [g·cm⁻¹·s⁻²]  current gas pressure at nodes (NN×1)
%                            (solver internal coordinate; NOT Pa)
%     T_guess     [°C]       initial temperature (NN×1)
%     refreshState  function_handle  refreshState(h, T, P_g, L_f)
%                    refreshes the full-column SFCC + gas/thermal state
%                    and returns structs (estate, ederiv) ready for
%                    evaluateInternalEnergyEndpoints.
%     dz          [cm]       element thicknesses, bottom-up (NL×1)
%     L_f         [J/kg]     per-step latent heat of fusion
%     opts        struct     optional: .T_min [-60], .T_max [90],
%                            .maxIter [50], .tol_T [1e-6],
%                            .use_local_fast_path [false], and
%                            .fast_max_iter [8].  The local fast path skips
%                            global endpoint evaluations only while every
%                            active node accepts safeguarded Newton updates.
%
%   Returns:
%     T            [°C]     inverted temperature (NN×1)
%     result       struct   .allConverged (logical scalar)
%                           .convergedMask (NN×1 logical)
%                           .iterations   (NN×1 double)
%                           .energyResidual (NN×1, J/cm³)
%                           .tempEquivalentResidual (NN×1, °C)
%                           .failureReason (NN×1 cellstr)
%
%   Physics guard: U_T ≤ 0 or non-finite → physics closure failure for
%   that node (no silent replacement with a positive floor).
%
%   Reference: docs/compass_soil_water_heat_physics.md §5.9.6.
%   Editor: Jizixing Zhu; Date: 2026-07-14.
%   Change: 2026-07-14 Jizixing Zhu -- keep non-monotonic energy brackets
%           outside the Newton active set and preserve their failure reason.
%   Change: 2026-07-17 Jizixing Zhu -- add an opt-in local Newton fast path
%           with exact fallback to the original global-bracket solver.
%
%   See also: soil.thermal.evaluateInternalEnergyEndpoints,
%             soil.thermal.assembleNodalInternalEnergy

    arguments
        U_target      (:,1) double
        h             (:,1) double
        P_g           (:,1) double
        T_guess       (:,1) double
        refreshState  function_handle
        dz            (:,1) double
        L_f           (1,1) double
        opts          struct = struct()
    end

    % --- inversion parameters ---
    T_MIN   = -60;    % [°C]
    T_MAX   =  90;    % [°C]
    MAX_ITER = 50;
    USE_LOCAL_FAST_PATH = false;
    FAST_MAX_ITER = 8;
    TOL_T   = 1e-6;   % [°C]

    if isfield(opts, 'T_min'),   T_MIN   = opts.T_min;   end
    if isfield(opts, 'T_max'),   T_MAX   = opts.T_max;   end
    if isfield(opts, 'maxIter'), MAX_ITER = opts.maxIter; end
    if isfield(opts, 'tol_T'),   TOL_T   = opts.tol_T;   end
    if isfield(opts, 'use_local_fast_path')
        USE_LOCAL_FAST_PATH = logical(opts.use_local_fast_path);
    end
    if isfield(opts, 'fast_max_iter')
        FAST_MAX_ITER = opts.fast_max_iter;
    end

    NL = numel(dz);
    NN = NL + 1;
    assert(numel(U_target) == NN && numel(h) == NN ...
        && numel(P_g) == NN && numel(T_guess) == NN, ...
        'invertTemperatureFromInternalEnergy:badDimensions', ...
        'U_target, h, P_g, T_guess must all be NN×1.');

    T = T_guess(:);

    % --- input guards ---
    assert(T_MIN < T_MAX, ...
        'invertTemperatureFromInternalEnergy:badBounds', ...
        'T_MIN (%g) must be < T_MAX (%g).', T_MIN, T_MAX);
    assert(L_f >= 0 && isfinite(L_f), ...
        'invertTemperatureFromInternalEnergy:badLf', ...
        'L_f must be finite nonnegative; got %g.', L_f);
    assert(all(T >= T_MIN & T <= T_MAX), ...
        'invertTemperatureFromInternalEnergy:badTguess', ...
        'T_guess must be in [%g, %g].', T_MIN, T_MAX);
    assert(isscalar(USE_LOCAL_FAST_PATH), ...
        'invertTemperatureFromInternalEnergy:badFastPathSwitch', ...
        'use_local_fast_path must be scalar logical-compatible.');
    assert(isscalar(FAST_MAX_ITER) && isfinite(FAST_MAX_ITER) ...
            && FAST_MAX_ITER >= 1 && FAST_MAX_ITER == floor(FAST_MAX_ITER), ...
        'invertTemperatureFromInternalEnergy:badFastMaxIter', ...
        'fast_max_iter must be a positive integer.');

    ctl = struct('L_f', L_f, 'validate', false);  % kernel validation off

    % --- initialise per-node diagnostics ---
    iterations       = zeros(NN, 1);
    energyResidual   = zeros(NN, 1);
    tempEquivResid   = zeros(NN, 1);
    failureReason    = repmat({''}, NN, 1);
    convergedMask    = false(NN, 1);

    % Production h/U/P_g iterations start from the previous Picard
    % temperature, so the safeguarded Newton step normally converges without
    % evaluating U at both global temperature limits. Return from this path
    % only when every node converges without bisection. Otherwise discard its
    % diagnostics and restart the unchanged protected algorithm below from
    % the caller's original temperature guess.
    if USE_LOCAL_FAST_PATH
        [T_fast, fastResult, fastSucceeded] = tryLocalFastNewton( ...
            U_target, h, P_g, T_guess, refreshState, dz, L_f, ctl, ...
            T_MIN, T_MAX, FAST_MAX_ITER, TOL_T);
        if fastSucceeded
            T = T_fast;
            result = fastResult;
            return;
        end
        T = T_guess(:);
    end

    % --- bracket: evaluate U at T_MIN and T_MAX for the whole column ---
    [estate_lo, ~] = refreshState(h, T_MIN * ones(NN,1), P_g, L_f);
    ederiv_lo = emptyDeriv(NL);
    [U_ep_lo] = soil.thermal.evaluateInternalEnergyEndpoints( ...
        estate_lo, ederiv_lo, ctl);
    [U_lo] = soil.thermal.assembleNodalInternalEnergy( ...
        U_ep_lo, zeros(NL,2), zeros(NL,2), zeros(NL,2), dz);

    [estate_hi, ~] = refreshState(h, T_MAX * ones(NN,1), P_g, L_f);
    [U_ep_hi] = soil.thermal.evaluateInternalEnergyEndpoints( ...
        estate_hi, ederiv_lo, ctl);
    [U_hi] = soil.thermal.assembleNodalInternalEnergy( ...
        U_ep_hi, zeros(NL,2), zeros(NL,2), zeros(NL,2), dz);

    % Bracket sanity: internal energy must be monotonic in T.
    nonMonotonic = U_lo >= U_hi;
    if any(nonMonotonic)
        failureReason(nonMonotonic) = { ...
            'U(T_min) >= U(T_max) — thermodynamics closure failure (non-monotonic)'};
        convergedMask(nonMonotonic) = false;
    end

    % Nodes with a valid monotonic bracket but an unreachable target indicate
    % an input error or failed trial.  Non-monotonic nodes retain the distinct
    % thermodynamic-closure failure assigned above.
    outOfRange = ~nonMonotonic & ...
        (U_target < U_lo | U_target > U_hi);
    if any(outOfRange)
        failureReason(outOfRange) = {'target U outside [U(T_min), U(T_max)]'};
        convergedMask(outOfRange) = false;
    end

    % --- active set: nodes still being solved ---
    active = ~(nonMonotonic | outOfRange);
    if ~any(active)
        result = packResult();
        return;
    end

    % Per-node bisection bracket (only for active nodes).
    % a/b track temperature bounds; U values are obtained from the
    % single refreshState evaluation at the top of each iteration.
    a = T_MIN * ones(NN, 1);
    b = T_MAX * ones(NN, 1);

    % ---- vectorized Newton loop (one refreshState per iteration) ----
    % At the TOP of each iteration, the current T (set by the previous
    % iteration's Newton/bisection step, or T_guess on iter=1) is
    % evaluated.  The result serves dual purpose: convergence check AND
    % bracket update.  Then a new T is computed for the next iteration.
    % This guarantees exactly ONE full-column SFCC refresh per iteration.
    for iter = 1:MAX_ITER
        if ~any(active), break; end

        % ---- single evaluation at current T ----
        [estate, ederiv] = refreshState(h, T, P_g, L_f);
        [U_ep, CTh, CTT, CTa] = soil.thermal.evaluateInternalEnergyEndpoints( ...
            estate, ederiv, ctl);
        [U_cur, U_T, ~, ~] = soil.thermal.assembleNodalInternalEnergy( ...
            U_ep, CTh, CTT, CTa, dz);

        % ---- update bracket from this evaluation ----
        shrinkLo = active & U_cur < U_target;
        shrinkHi = active & U_cur >= U_target;
        a(shrinkLo) = T(shrinkLo);
        b(shrinkHi) = T(shrinkHi);

        % ---- residual and diagnostics ----
        R = U_cur - U_target;  % NN×1
        idx = find(active);
        iterations(idx) = iter;
        energyResidual(idx) = abs(R(idx));
        tempEquivResid(idx) = abs(R(idx)) ./ max(abs(U_T(idx)), eps);

        % ---- physics guard: U_T must be strictly positive ----
        badUT = active & (U_T <= 0 | ~isfinite(U_T));
        if any(badUT)
            failureReason(badUT) = {'U_T <= 0 or non-finite (physics closure failure)'};
            convergedMask(badUT) = false;
            active(badUT) = false;
        end
        if ~any(active), break; end

        % ---- convergence check ----
        newConv = active & (abs(R) < TOL_T .* U_T);
        if any(newConv)
            convergedMask(newConv) = true;
            energyResidual(newConv)  = abs(R(newConv));
            tempEquivResid(newConv)  = abs(R(newConv)) ./ max(U_T(newConv), eps);
            iterations(newConv) = iter;
            active(newConv) = false;
        end
        if ~any(active), break; end

        % ---- Newton step for all active nodes ----
        dT = -R ./ U_T;  % NN×1; only active entries matter

        % Accept Newton step where it stays within the bracket and is
        % not too large relative to the bracket width.
        inBracket = T + dT > a & T + dT < b;
        stepOk    = abs(dT) < 0.5 * (b - a);
        useNewton = active & inBracket & stepOk;

        % Bisection fallback for the rest.
        useBisect = active & ~useNewton;

        if any(useNewton)
            T(useNewton) = T(useNewton) + dT(useNewton);
        end
        if any(useBisect)
            T(useBisect) = (a(useBisect) + b(useBisect)) / 2;
        end
        % T is now set for the next iteration's evaluation at the top.
    end

    % Nodes still active after MAX_ITER → did not converge.
    if any(active)
        failureReason(active) = {'exceeded max iterations'};
        convergedMask(active) = false;
    end

    result = packResult();

    % ---- nested helpers --------------------------------------------------
    function r = packResult()
        r = struct( ...
            'allConverged', all(convergedMask), ...
            'convergedMask', convergedMask, ...
            'iterations', iterations, ...
            'energyResidual', energyResidual, ...
            'tempEquivalentResidual', tempEquivResid, ...
            'failureReason', {failureReason});
    end
end

% -------------------------------------------------------------------------
function [T, result, succeeded] = tryLocalFastNewton( ...
        U_target, h, P_g, T_guess, refreshState, dz, L_f, ctl, ...
        T_MIN, T_MAX, maxIter, tolT)
%TRYLOCALFASTNEWTON  Invert U(T) without eager global endpoint evaluations.
%
%   This numerical fast path uses the same local residual, derivative,
%   convergence test, and Newton admissibility test as the protected solver.
%   It succeeds only if every active node remains on Newton updates.  A need
%   for bisection, a nonpositive U_T, or exhaustion returns succeeded=false
%   so the caller can restart the complete global-bracket algorithm.
%
%   Inputs and outputs retain the units documented by
%   invertTemperatureFromInternalEnergy; no physical state is committed.
%   Editor: Jizixing Zhu; Date: 2026-07-17.

    NL = numel(dz);
    NN = NL + 1;
    T = T_guess(:);
    a = T_MIN * ones(NN, 1);       % [degC] local lower guard
    b = T_MAX * ones(NN, 1);       % [degC] local upper guard
    active = true(NN, 1);
    iterations = zeros(NN, 1);
    energyResidual = zeros(NN, 1); % [J cm^-2]
    tempEquivResid = zeros(NN, 1); % [degC]
    failureReason = repmat({''}, NN, 1);
    convergedMask = false(NN, 1);
    succeeded = false;

    for iter = 1:maxIter
        [estate, ederiv] = refreshState(h, T, P_g, L_f);
        [U_ep, CTh, CTT, CTa] = ...
            soil.thermal.evaluateInternalEnergyEndpoints(estate, ederiv, ctl);
        [U_cur, U_T] = soil.thermal.assembleNodalInternalEnergy( ...
            U_ep, CTh, CTT, CTa, dz);

        R = U_cur - U_target;       % [J cm^-2]
        idx = find(active);
        iterations(idx) = iter;
        energyResidual(idx) = abs(R(idx));
        tempEquivResid(idx) = abs(R(idx)) ./ max(abs(U_T(idx)), eps);

        if any(active & (U_T <= 0 | ~isfinite(U_T) | ~isfinite(R)))
            result = packFastResult();
            return;
        end

        shrinkLo = active & U_cur < U_target;
        shrinkHi = active & U_cur >= U_target;
        a(shrinkLo) = T(shrinkLo);
        b(shrinkHi) = T(shrinkHi);

        newConv = active & abs(R) < tolT .* U_T;
        if any(newConv)
            convergedMask(newConv) = true;
            active(newConv) = false;
        end
        if ~any(active)
            succeeded = true;
            result = packFastResult();
            return;
        end

        dT = -R ./ U_T;             % [degC]
        candidateT = T + dT;        % [degC]
        useNewton = active & isfinite(dT) ...
            & candidateT > a & candidateT < b ...
            & abs(dT) < 0.5 * (b - a);
        if any(active & ~useNewton)
            result = packFastResult();
            return;
        end
        T(useNewton) = candidateT(useNewton);
    end

    result = packFastResult();

    function r = packFastResult()
        r = struct( ...
            'allConverged', all(convergedMask), ...
            'convergedMask', convergedMask, ...
            'iterations', iterations, ...
            'energyResidual', energyResidual, ...
            'tempEquivalentResidual', tempEquivResid, ...
            'failureReason', {failureReason});
    end
end

% -------------------------------------------------------------------------
function ederiv = emptyDeriv(NL)
% Return a zero-filled ederiv struct for calls that only need U (not CTh/CTT/CTa).
% Editor: Jizixing Zhu; Date: 2026-07-14.
    Z2  = zeros(NL, 2);
    ZN  = zeros(NL+1, 1);
    ederiv = struct( ...
        'DThetaLiquidDh', Z2, 'DThetaLiquidDT', Z2, ...
        'DThetaIceDh',    Z2, 'DThetaIceDT',    Z2, ...
        'DThetaGasDh',    Z2, 'DThetaGasDT',    Z2, ...
        'DRHOVh', ZN, 'DRHOVT', ZN, ...
        'Xah',    ZN, 'XaT',    ZN, 'Xaa',     ZN, ...
        'dW_dh',  Z2, 'dW_dTheta', Z2);
end
