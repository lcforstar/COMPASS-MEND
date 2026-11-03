function plan = winterWheatGrainNitrogen(C_prev, N_prev, enter_grain_fill, params)
%WINTERWHEATGRAINNITROGEN  One-time, conservative crop N remobilization plan.
%
% Inputs (read-only):
%   C_prev [scalar struct] Accepted beginning-of-step organ carbon [g C/m2]:
%       C_leaf, C_livestem, C_root_fine. All must be finite and nonnegative.
%   N_prev [scalar struct] Corresponding accepted organ nitrogen [g N/m2]:
%       N_leaf, N_livestem, N_root_fine. All must be finite and nonnegative.
%   enter_grain_fill [logical or numeric scalar 0/1] One accepted-to-pending
%       transition into live winter-wheat grain filling. The caller owns this
%       edge; a phase flag held true for a whole day is NOT an event edge.
%   params [scalar struct] Finite positive residual-organ C:N [g C/g N]:
%       CN_leaf_grainfill, CN_stem_grainfill, CN_root_grainfill. Required:
%       there are no hidden defaults or additional time-scale parameters.
%       Pools and C:N values must be real floating-point scalars; integer
%       storage is rejected to prevent rounded integer division/subtraction.
%
% Output:
%   plan [scalar struct] Nonnegative step-integrated transfers [g N/m2]:
%       leaf, stem, root are donor debits; to_labile is their exact sum.
%       No carbon moves, no external N is added, and no state is mutated.
%       The caller applies each donor debit and the matching mobile-N credit
%       in the SAME accepted transaction. The values are amounts, not rates.
%
% Mathematics, on the single grain-fill event:
%   debit_i = max(0, N_i_actual - C_i / CN_i_grainfill).
%   delta_N_donor_i = -debit_i; delta_N_labile = sum(debit_i).
% This preserves each donor's residual-N floor when that floor is initially
% affordable, never draws more N than is present, and leaves an already
% N-poor donor unchanged. On other steps every output amount is zero.
%
% Source and deliberate COMPASS adaptation:
%   CoLM commit 626347a9156ab97db4ca846f47cd1b7f2728f705,
%   main/BGC/MOD_BGC_Veg_NutrientCompetition.F90:418-443, 506-532:
%   https://github.com/CoLM-SYSU/CoLM/blob/626347a9156ab97db4ca846f47cd1b7f2728f705/main/BGC/MOD_BGC_Veg_NutrientCompetition.F90#L418
%   CoLM transfers C_i*(1/CN_base_i-1/CN_grainfill_i) once, expressed as a
%   rate by division by deltim. Its old retransn pool supplies current demand;
%   its state update later adds the newly mobilized N. This kernel does not
%   prescribe the caller's allocation timing or change its soil-uptake flux.
%   CoLM main/MOD_Const_PFT.F90:1482-1525 gives winter-wheat (index 21)
%   grain-fill C:N = 65, 100, 40 for leaf, stem, fine root respectively:
%   https://github.com/CoLM-SYSU/CoLM/blob/626347a9156ab97db4ca846f47cd1b7f2728f705/main/MOD_Const_PFT.F90#L1482
%   COMPASS FlexibleCN does not guarantee N_actual=C/CN_base. Therefore this
%   kernel uses the actual donor N and allows only donor-to-mobile transfer.
%   In particular, CoLM's base fine-root C:N=42 and grain-fill value=40 would
%   give a negative original flux; this kernel does NOT manufacture root N
%   through that reverse transfer. This is an explicit conservative adaptation,
%   not a claim of complete CoLM parity. Equality to the original integrated
%   flux holds when N_actual=C/CN_base and that original flux is nonnegative.
%
% Errors: malformed structs, absent fields, invalid event flags, nonfinite
% or negative pools, nonpositive C:N parameters, and nonfinite total transfer
% are rejected. Event=false does not bypass validation of scientific inputs.
%
% Editor: Wenjie Liu; Date: 2026-09-30.

    assert(isstruct(C_prev) && isscalar(C_prev), ...
        'winterWheatGrainNitrogen:invalidCarbonState', ...
        'C_prev must be a scalar struct of accepted organ carbon pools.');
    assert(isstruct(N_prev) && isscalar(N_prev), ...
        'winterWheatGrainNitrogen:invalidNitrogenState', ...
        'N_prev must be a scalar struct of accepted organ nitrogen pools.');
    assert(isstruct(params) && isscalar(params), ...
        'winterWheatGrainNitrogen:invalidParameters', ...
        'params must be a scalar struct with explicit grain-fill C:N values.');
    assert((islogical(enter_grain_fill) || isnumeric(enter_grain_fill)) ...
        && isreal(enter_grain_fill) && isscalar(enter_grain_fill) ...
        && isfinite(enter_grain_fill) ...
        && (enter_grain_fill == 0 || enter_grain_fill == 1), ...
        'winterWheatGrainNitrogen:invalidEvent', ...
        'enter_grain_fill must be a scalar logical or numeric 0/1 event edge.');

    carbon = [requiredScalar(C_prev, 'C_leaf', false), ...
        requiredScalar(C_prev, 'C_livestem', false), ...
        requiredScalar(C_prev, 'C_root_fine', false)];
    nitrogen = [requiredScalar(N_prev, 'N_leaf', false), ...
        requiredScalar(N_prev, 'N_livestem', false), ...
        requiredScalar(N_prev, 'N_root_fine', false)];
    cnResidual = [requiredScalar(params, 'CN_leaf_grainfill', true), ...
        requiredScalar(params, 'CN_stem_grainfill', true), ...
        requiredScalar(params, 'CN_root_grainfill', true)];

    debit = zeros(1, 3);
    if enter_grain_fill
        debit = max(0, nitrogen - carbon ./ cnResidual);
    end
    plan = struct('leaf', debit(1), 'stem', debit(2), ...
        'root', debit(3), 'to_labile', sum(debit));
    assert(all(isfinite(debit)) && isfinite(plan.to_labile) ...
        && all(debit >= 0) && all(debit <= nitrogen), ...
        'winterWheatGrainNitrogen:invalidTransfer', ...
        'Grain-fill N transfers must be finite and bounded by actual donor N.');
end

function value = requiredScalar(container, name, strictlyPositive)
% Validate one scalar pool or C:N parameter without defaults or conversions.
    assert(isfield(container, name), ...
        'winterWheatGrainNitrogen:missingField', ...
        'Required scientific input %s is absent.', name);
    value = container.(name);
    assert(isfloat(value) && isreal(value) && isscalar(value) ...
        && isfinite(value), ...
        'winterWheatGrainNitrogen:invalidScalar', ...
        '%s must be a finite real floating-point scalar.', name);
    if strictlyPositive
        assert(value > 0, 'winterWheatGrainNitrogen:invalidCN', ...
            '%s must be a positive C:N ratio.', name);
    else
        assert(value >= 0, 'winterWheatGrainNitrogen:negativePool', ...
            '%s must be a nonnegative organ pool.', name);
    end
end
