function [limited,info]=limitNitrificationByNodeBudget(rate,nodeBudget_ug_cm2_s,dz)
%LIMITNITRIFICATIONBYNODEBUDGET Limit consistent FE loads before O2/N ledgers.
% RATE is a finite nonnegative bottom-up NL-by-2 endpoint NH4 oxidation rate
% [ug N/cm3 soil/s]. NODEBUDGET_UG_CM2_S is a finite nonnegative NN-vector of
% shared-node available reaction loads [ug N/cm2/s], frozen for the caller's
% reaction interval. DZ is a finite positive bottom-up NL-vector [cm].
%
% Assemble with the ADE's consistent linear-element quadrature. Each element
% receives the minimum request/budget factor at its two incident nodes and
% both endpoint rates are scaled together. This preserves stoichiometry
% when the caller derives O2, NH4, NO3, N2O and tracer ledgers from LIMITED.
% A zero node allowance can suppress its high-concentration neighbor in the
% same element; this is the conservative consistent-load closure, not an
% independent kinetic threshold. This helper does not change physical
% stocks, root sinks, boundary exchanges or transport coefficients.
%
% Feasible input RATE is returned directly. INFO records requested, applied
% and reduced node/column loads [ug N/cm2/s], endpoint rate reductions
% [ug N/cm3/s], and dimensionless NL-vector element_scale. Invalid shapes,
% negative rates/budgets or nonfinite data raise the named input error.
% Source: limitNitrificationAtNodes.m, identical ADE source quadrature.
% Editor: Wenjie Liu; Date: 2026-11-01.
    dz=dz(:);nodeBudget_ug_cm2_s=nodeBudget_ug_cm2_s(:);nl=numel(dz);
    assert(isnumeric(rate)&&isreal(rate)&&isnumeric(dz)&&isreal(dz) ...
        &&isnumeric(nodeBudget_ug_cm2_s)&&isreal(nodeBudget_ug_cm2_s) ...
        &&isequal(size(rate),[nl,2])&&numel(nodeBudget_ug_cm2_s)==nl+1 ...
        &&all(isfinite([rate(:);nodeBudget_ug_cm2_s;dz])) ...
        &&all(rate>=0,'all')&&all(nodeBudget_ug_cm2_s>=0)&&all(dz>0), ...
        'soil:nitrogen:nodalNitrificationBudgetInputs', ...
        'Expected nonnegative NL-by-2 rates, NN node budgets and positive NL thicknesses.');
    load=@(v)[(2*v(:,1)+v(:,2)).*dz/6;0] ...
        +[0;(v(:,1)+2*v(:,2)).*dz/6];
    requested=load(rate);
    factor=ones(nl+1,1);
    active=requested>nodeBudget_ug_cm2_s;
    if any(active)
        factor(active)=nodeBudget_ug_cm2_s(active)./requested(active)*(1-8*eps);
        element=min(factor(1:end-1),factor(2:end));
        limited=rate.*element;
    else
        element=ones(nl,1);
        limited=rate;
    end
    applied=load(limited);
    reduced=requested-applied;
    info=struct('node_available_ug_cm2_s',nodeBudget_ug_cm2_s, ...
        'node_requested_ug_cm2_s',requested,'node_applied_ug_cm2_s',applied, ...
        'node_reduced_ug_cm2_s',reduced,'element_scale',element, ...
        'column_requested_ug_cm2_s',sum(requested), ...
        'column_applied_ug_cm2_s',sum(applied), ...
        'column_reduced_ug_cm2_s',sum(reduced), ...
        'rate_reduction_ug_cm3_s',rate-limited, ...
        'source_quadrature','consistent linear finite element');
end
