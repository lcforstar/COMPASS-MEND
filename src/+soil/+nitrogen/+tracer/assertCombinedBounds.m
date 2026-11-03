function assertCombinedBounds(tags,total)
%ASSERTCOMBINEDBOUNDS Bound mutually exclusive passive-source cohorts jointly.
% tags is a cell of at most three detached tracer states (mineral fertilizer,
% manure organic N, manure TAN). total is the identical accepted owner bundle
% used by assertBounds, including total.aom for manure. Nodal units are ug N
% cm^-3 water and element units are g N cm^-3 bulk; surface AOM is g N cm^-2.
% This validation never redistributes labels or changes physical states.
% Editor: Wenjie Liu; Date: 2026-10-30.
    assert(iscell(tags)&&~isempty(tags)&&numel(tags)<=3, ...
        'soil:nitrogen:tracer:cohortSet','One to three exclusive label states are required.');
    for i=1:numel(tags),soil.nitrogen.tracer.assertBounds(tags{i},total);end
    combined=tags{1};
    for name={'nh4','no3','urea','nh4_granule','nh4_ice','no3_ice'}
        key=name{1};combined.(key)=zeros(size(combined.(key)));
        for i=1:numel(tags),combined.(key)=combined.(key)+tags{i}.(key);end
    end
    for name=fieldnames(combined.organic_n)'
        key=name{1};combined.organic_n.(key)=zeros(size(combined.organic_n.(key)));
        for i=1:numel(tags)
            combined.organic_n.(key)=combined.organic_n.(key)+tags{i}.organic_n.(key);
        end
    end
    for name={'don_total_g_cm3','don_sorbed_g_cm3'}
        key=name{1};if ~isfield(combined.dom,key),continue;end
        combined.dom.(key)=zeros(size(combined.dom.(key)));
        for i=1:numel(tags),combined.dom.(key)=combined.dom.(key)+tags{i}.dom.(key);end
    end
    hasAom=cellfun(@(x)isfield(x,'aom'),tags);
    if any(hasAom)
        combined.aom=struct('manure',struct('surface_N',0, ...
            'element_N',zeros(size(total.aom.manure.element_N))));
        for i=reshape(find(hasAom),1,[])
            combined.aom.manure.surface_N=combined.aom.manure.surface_N+tags{i}.aom.manure.surface_N;
            combined.aom.manure.element_N=combined.aom.manure.element_N+tags{i}.aom.manure.element_N;
        end
    end
    soil.nitrogen.tracer.assertBounds(combined,total);
end
