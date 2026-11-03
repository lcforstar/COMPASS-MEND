function validateOrganicDOMObservations(built,raw)
%VALIDATEORGANICDOMOBSERVATIONS Check optional DOM against resolved CSV priors.
% Final accepted t0 liquid capacity is used for free-D or aqueous inputs.
% Detailed priors own mass; optional profile columns are consistency checks,
% never an additional source or a reason to repartition MOC/MON. Cold format
% keeps the existing DOM input semantics. Finite nonnegative values required.
% Editor: Wenjie Liu; Date: 2026-10-17.
    if ~isfield(built.parameters,'soil_organic') ...
            ||~isfield(built.parameters.soil_organic,'initial_dom'),return;end
    d=built.Y.soil_dom;rho=built.parameters.soil.bulk_density(:);
    theta=mean(built.workspace.derived.Theta_LL,2);
    z=.5*(built.domain.z_node(1:end-1)+built.domain.z_node(2:end));
    for kind={'doc','don'}
        k=kind{1};p=soil.dom.phaseState(d,k,theta,rho,built.parameters.soil_dom);
        names={[k '_mgkg'],[k '_mgL']};values={d.([k '_total_g_cm3'])./rho*1e6,p.aqueous_g_cm3_water*1e6};
        if strcmp(k,'doc'),names{end+1}='doc_free_mgkg';values{end+1}=p.free_g_cm3_bulk./rho*1e6;end
        for j=1:numel(names)
            name=names{j};tbl=raw.initial;
            if ~ismember(name,tbl.Properties.VariableNames)||all(ismissing(tbl.(name))),continue;end
            v=tbl.(name);assert(isnumeric(v)&&all(isfinite(v))&&all(v>=0), ...
                'organicInitial:DOMObservation','Optional %s must be wholly blank or finite and nonnegative.',name);
            target=init.interpLinearWithExtrapolation(tbl.depth_cm,v,z);
            assert(all(abs(target-values{j})<=1e-12.*max(abs(target),1)), ...
                'organicInitial:DOMConflict','Optional %s conflicts with resolved organic priors.',name);
        end
    end
end
