function spec=normalizeColumnBoundary(b,prefixes,particle)
%NORMALIZECOLUMNBOUNDARY Organic-only typed boundary configuration adapter.
% Types follow the supported DIN interface: top 1 prescribed concentration,
% top 2 prescribed inflow concentration/nonvolatile evaporative boundary;
% bottom 1 prescribed concentration, 2 inflow concentration/outflow gradient
% zero, 3 outflow-only gradient zero. Top 3 is NOT advertised: DIN assembly
% and flux recovery disagree there. Type 4 has no implemented DIN equation.
% Site concentrations are mg element/L water (not mg/kg soil); convert once
% to g/cm3 water. Particle vectors are [MB BD ENZN], with explicit C AND N.
% Current water integration is physical free drainage/liquid_only, as before.
% Editor: Wenjie Liu; Date: 2026-09-29.
    if particle,owner='bioparticle';else,owner='dom';end
    assert(~isfield(b,[owner '_bottom_kind']),'soil:dom:ambiguousBoundary', ...
        'Do not mix numbered organic boundaries with the legacy %s_bottom_kind.',owner);
    assert(isfield(b,[owner '_bottom_water_phase']) ...
        &&strcmp(b.([owner '_bottom_water_phase']),'liquid_only') ...
        &&isfield(b,'water_bottom_kind')&&strcmp(b.water_bottom_kind,'free_drainage'), ...
        'soil:dom:boundaryWaterPolicy', ...
        'Numbered organic boundaries currently require free_drainage and explicit liquid_only.');
    spec=struct('schema_version',2,'bottom_kind','configured','top_kind','configured', ...
        'water_phase','liquid_only','water_bottom_driver','accepted_physical_QMB', ...
        'water_top_driver','accepted_top_element_QL_like_DIN_node_adapter', ...
        'concentration_units','g_element_cm3_water','ledger_version',2,'carriers',struct());
    allowed={[owner '_bottom_water_phase']};
    for j=1:numel(prefixes)
        prefix=prefixes{j};
        for edge={'top','bottom'}
            side=edge{1};key=[prefix '_' side];
            assert(isfield(b,key)&&isnumeric(b.(key))&&isscalar(b.(key)), ...
                'soil:dom:boundaryConfig','Explicit %s type required.',key);
            types=[1 2 3];if strcmp(side,'top'),types=[1 2];end
            assert(ismember(b.(key),types),'soil:dom:boundaryConfig', ...
                '%s supports only types %s.',key,mat2str(types));
            allowed{end+1}=key; %#ok<AGROW>
            if particle
                elements={'carbon','nitrogen'};groups={'MB','BD','ENZN'};
            else
                elements={'value'};groups={prefix};
            end
            for e=1:numel(elements)
                valueKey=[key '_' elements{e} '_mg_L'];
                assert(isfield(b,valueKey)&&isnumeric(b.(valueKey)) ...
                    &&numel(b.(valueKey))==numel(groups) ...
                    &&all(isfinite(b.(valueKey)))&&all(b.(valueKey)>=0), ...
                    'soil:dom:boundaryConcentration','Explicit nonnegative %s required.',valueKey);
                values=b.(valueKey);
                assert(b.(key)~=3||all(values==0),'soil:dom:unusedBoundaryValue', ...
                    'Type 3 uses no external concentration; its value must be zero.');
                allowed{end+1}=valueKey; %#ok<AGROW>
                for g=1:numel(groups)
                    name=groups{g};if particle,name=[name '_' elements{e}];end
                    spec.carriers.(name).(side)=struct('type',b.(key), ...
                        'concentration_g_cm3_water',double(values(g))*1e-6);
                end
            end
        end
    end
    if particle
        for side={'top','bottom'}
            c=b.(['bioparticle_' side{1} '_carbon_mg_L']);
            n=b.(['bioparticle_' side{1} '_nitrogen_mg_L']);
            assert(all((c==0)==(n==0)),'bioparticle:boundaryCN', ...
                'Every external particle cargo must specify paired positive C/N or paired zero.');
            assert(c(3)==0||c(1)+c(2)>0,'bioparticle:boundaryCargo', ...
                'External intracellular enzyme requires an external mobile cell carrier.');
        end
    end
    keys=fieldnames(b);selected=false(size(keys));
    for p=[prefixes,{owner}],selected=selected|startsWith(keys,[p{1} '_']);end
    assert(isempty(setdiff(keys(selected),allowed)),'soil:dom:boundaryConfig', ...
        'Unknown organic boundary key; refusing an unconsumed setting.');
end
