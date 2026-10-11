function [summary,detail]=compareSpinupCycles(startVector,endVector,policy)
%COMPARESPINUPCYCLES Compare independent accepted stocks at equal cycle phase.
% Inputs are tables from soil.organic.spinupStateVector. POLICY must supply
% absolute_C_g_cm2, absolute_N_g_cm2 (nonnegative g element/cm2 per location),
% relative_tolerance in [0,1), and nonempty provenance explaining thresholds.
% No absolute floor or scientific tolerance is selected implicitly. Each row
% passes when |end-start| <= absolute_element + relative*max(|start|,|end|).
% The symmetric scale handles depletion/appearance without zero division;
% location-wise checks cannot be hidden by offsetting column changes.
% SUMMARY.state_gate_pass reports this gate only. Its explicit
% periodic_equilibrium_certified=false prevents a stock match from asserting
% full annual convergence: caller must check full-duration same-phase forcing,
% management, water/heat, flux stability and conservation independently.
% Rejects missing/invalid policy, mismatched geometry/layout/units, or NaN/Inf.
% Does not modify state or the supplied tables. Editor: Wenjie Liu; Date: 2026-09-30.
    required={'absolute_C_g_cm2','absolute_N_g_cm2','relative_tolerance','provenance'};
    require(isstruct(policy)&&isscalar(policy)&&all(isfield(policy,required)), ...
        'policy','Explicit unit-aware tolerances and provenance are required.');
    numbers=[policy.absolute_C_g_cm2,policy.absolute_N_g_cm2,policy.relative_tolerance];
    require(isnumeric(numbers)&&isreal(numbers)&&numel(numbers)==3 ...
        &&all(isfinite(numbers))&&all(numbers>=0)&&numbers(3)<1, ...
        'policy','Invalid absolute or relative tolerances.');
    provenance=string(policy.provenance);
    require(isscalar(provenance)&&~ismissing(provenance)&&strlength(strtrim(provenance))>0, ...
        'policy','Tolerance provenance must be nonempty.');
    keys={'state_path','element','location','element_index','endpoint_index', ...
        'source_index','native_unit','canonical_unit'};
    for item={startVector,endVector}
        v=item{1};
        require(istable(v)&&height(v)>0&&all(ismember([keys {'mass_g_cm2'}], ...
            v.Properties.VariableNames)),'vector','Input must be a spinupStateVector table.');
        require(isequal(v.Properties.UserData.schema_version,1) ...
            &&all(ismember(v.element,["C","N"])) ...
            &&isequal(v.canonical_unit,"g_"+v.element+"/cm2") ...
            &&isa(v.mass_g_cm2,'double')&&isreal(v.mass_g_cm2) ...
            &&all(isfinite(v.mass_g_cm2)),'vector','Invalid canonical stock units or values.');
    end
    require(isequal(startVector(:,keys),endVector(:,keys)) ...
        &&isequaln(startVector.Properties.UserData,endVector.Properties.UserData), ...
        'layout','Cycle geometry, owner, phase layout or unit contract differs.');
    a=startVector.mass_g_cm2;b=endVector.mass_g_cm2;scale=max(abs(a),abs(b));
    absolute=repmat(double(policy.absolute_C_g_cm2),height(startVector),1);
    absolute(startVector.element=="N")=double(policy.absolute_N_g_cm2);
    delta=b-a;change=abs(delta);allowed=absolute+double(policy.relative_tolerance)*scale;
    require(all(isfinite(allowed)),'policy','Tolerance arithmetic overflow.');
    passed=change<=allowed;ratio=zeros(size(change));active=allowed>0;
    ratio(active)=change(active)./allowed(active);ratio(~active&change>0)=Inf;
    relative=zeros(size(change));active=scale>0;relative(active)=change(active)./scale(active);
    detail=startVector(:,keys);detail.start_g_cm2=a;detail.end_g_cm2=b;
    detail.delta_g_cm2=delta;detail.absolute_change_g_cm2=change;
    detail.absolute_tolerance_g_cm2=absolute;detail.allowed_change_g_cm2=allowed;
    detail.symmetric_relative_change=relative;detail.tolerance_ratio=ratio;detail.passed=passed;
    [worst,index]=max(ratio);
    summary=struct('schema_version',1,'state_gate_pass',all(passed), ...
        'periodic_equilibrium_certified',false,'n_compared_states',height(detail), ...
        'n_failed_states',sum(~passed),'max_tolerance_ratio',worst, ...
        'worst_state_path',char(detail.state_path(index)), ...
        'worst_element',char(detail.element(index)), ...
        'worst_element_index',detail.element_index(index), ...
        'worst_endpoint_index',detail.endpoint_index(index), ...
        'worst_absolute_change_g_cm2',change(index), ...
        'max_symmetric_relative_change',max(relative),'policy',policy, ...
        'criterion','abs(end-start) <= abs_element_g_cm2 + rtol*max(abs(start),abs(end))', ...
        'requires_external_gates','full_cycle_same_phase_water_heat_management_flux_conservation');
    for el=["C","N"]
        match=detail.element==el;
        summary.([char(el) '_net_change_g_cm2'])=sum(delta(match));
        summary.([char(el) '_sum_absolute_change_g_cm2'])=sum(change(match));
    end
end

function require(condition,id,message)
    assert(condition,['soil:organic:compareSpinup:' id],'%s',message);
end
