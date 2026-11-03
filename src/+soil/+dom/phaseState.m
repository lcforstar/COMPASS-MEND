function phase=phaseState(state,kind,theta,bulkDensity,parameters,rewetting)
%PHASESTATE Read the sole DOM owner's free and sorbed phase inventories.
% State totals and optional dynamic sorbed stores are bottom-up NL columns
% [g element/cm3 bulk soil]. theta [-], bulkDensity [g/cm3] and Kd [cm3/g]
% have the existing DOM contract. This reader does not repartition dynamic
% stores: free=total-sorbed, Cwater=free/theta, Csorbed=sorbed/rho_b.
% Equilibrium states use soil.dom.equilibrate unchanged. Dynamic states must
% explicitly carry sorption_model='mend_dynamic' and schema_version>=2.
% Zero liquid capacity with nonzero free mass is rejected; no water floor,
% instantaneous adsorption or hidden state migration is introduced.
% Editor: Jizixing Zhu; Date: 2026-09-14.
    assert(ismember(string(kind),["doc","don"]), ...
        'soil:dom:phaseKind','DOM kind must be doc or don.');
    kind=char(kind);total=state.([kind '_total_g_cm3']);
    theta=theta(:);bulkDensity=bulkDensity(:);
    assert(iscolumn(total)&&isequal(size(total),size(theta),size(bulkDensity)) ...
        &&all(isfinite([total;theta;bulkDensity]))&&all(total>=0) ...
        &&all(theta>=0)&&all(bulkDensity>0), ...
        'soil:dom:phaseShape','DOM phase inventory and drivers must align.');
    mode=string(parameters.sorption_model);
    tagged=isfield(state,'sorption_model');
    if tagged
        assert(isscalar(string(state.sorption_model))&&string(state.sorption_model)==mode, ...
            'soil:dom:phaseMode','DOM state and parameter sorption modes differ.');
    end
    switch mode
        case "linear_equilibrium"
            assert(~isfield(state,'doc_sorbed_g_cm3')&&~isfield(state,'don_sorbed_g_cm3'), ...
                'soil:dom:phaseMode','Dynamic phase stores cannot be silently re-equilibrated.');
            if nargin<6,rewetting=false(size(theta));end
            phase=soil.dom.equilibrate(total,theta,bulkDensity,parameters.(['kd_' kind '_cm3_g']),rewetting);
            phase.free_g_cm3_bulk=theta.*phase.aqueous_g_cm3_water;
            phase.sorbed_g_cm3_bulk=bulkDensity.*phase.sorbed_g_g_soil;
        case "mend_dynamic"
            assert(tagged&&isfield(state,'schema_version')&&state.schema_version>=2 ...
                &&isfield(state,[kind '_sorbed_g_cm3']), ...
                'soil:dom:dynamicSchema','Dynamic DOM requires explicit accepted phase stores.');
            sorbed=state.([kind '_sorbed_g_cm3']);
            assert(isequal(size(sorbed),size(total))&&all(isfinite(sorbed)) ...
                &&all(sorbed>=0)&&all(sorbed<=total), ...
                'soil:dom:dynamicInventory','Sorbed DOM must be nonnegative and not exceed total.');
            free=total-sorbed;wet=theta>0;
            assert(all(free(wet==0)==0),'soil:dom:dynamicDryFreeMass', ...
                'Nonzero free DOM in a zero-liquid element requires an explicit dry-phase policy.');
            aqueous=zeros(size(free));aqueous(wet)=free(wet)./theta(wet);
            phase=struct('aqueous_g_cm3_water',aqueous,'aqueous_mg_L',aqueous*1e6, ...
                'sorbed_g_g_soil',sorbed./bulkDensity,'total_g_cm3_bulk',total, ...
                'free_g_cm3_bulk',free,'sorbed_g_cm3_bulk',sorbed, ...
                'closure_residual_g_cm3',total-free-sorbed);
        otherwise
            error('soil:dom:phaseMode','Unsupported DOM sorption model %s.',mode);
    end
end
