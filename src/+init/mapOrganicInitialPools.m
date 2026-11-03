function prior=mapOrganicInitialPools(tbl,centersTop,socTop,mobileFraction)
%MAPORGANICINITIALPOOLS Materialize resolved CSV priors once on soil elements.
% Inputs: table at positive top-down centers [cm], centersTop NLx1 [cm],
% socTop NLx1 [g C/cm3 bulk soil], scalar default mobile fraction. Returns
% C [g/cm3] top-down, N [g/cm3] bottom-up, parent-relative phase fractions,
% AOM [g/cm3 element or g/cm2 surface]. Endpoint-held linear interpolation;
% exact matching centers are copied without interpolation roundoff.
% Invalid/mixed columns, closure, phases or extrapolated negativity fail.
% No MAT I/O, state overwrites, integration or cumulative flux restoration.
% Editor: Wenjie Liu; Date: 2026-10-17.
    spec=io.adapters.organicPoolSpec();names=tbl.Properties.VariableNames;
    detailed=any(ismember({'frac_POM1','frac_POM2','frac_DOC'},names));
    prior=struct('detailed',detailed);
    if ~detailed,return;end
    assert(all(ismember(spec.detailed,names)),'organicInitial:columns', ...
        'Detailed priors require all twelve carbon fraction columns.');
    assert(~ismember('frac_POC',names)&&~ismember('depth_top_cm',names), ...
        'organicInitial:mixedFormat','Detailed priors use layer-center rows, without cold POC or horizon columns.');
    x=tbl.depth_cm(:);z=centersTop(:);n=numel(z);
    assert(isnumeric(x)&&all(isfinite(x))&&all(x>=0)&&all(diff(x)>0), ...
        'organicInitial:depth','Depths must be finite, nonnegative and strictly increasing.');
    fractions=tbl{:,strcat('frac_',spec.carbon)};
    assert(isnumeric(fractions)&&all(isfinite(fractions),'all') ...
        &&all(fractions>=0&fractions<=1,'all') ...
        &&all(abs(sum(fractions,2)-1)<=16*eps), ...
        'organicInitial:carbonClosure','Twelve carbon fractions must sum to one to floating-point precision.');
    prior.soc_top_down=socTop(:);prior.carbon=struct();prior.nitrogen=struct();
    for j=1:numel(spec.carbon)
        key=spec.carbon{j};prior.carbon.(key)=socTop(:).*mapped(['frac_' key]);
    end
    for j=1:numel(spec.nitrogen)
        key=spec.nitrogen{j};a=['nfrac_' key];b=['nc_' key];
        assert(ismember(a,names)||ismember(b,names),'organicInitial:nitrogen', ...
            'Provide %s (N/SOC) or %s (paired N/C).',a,b);
        if ismember(a,names),value=socTop(:).*mapped(a);
        else,value=prior.carbon.(spec.paired_carbon{j}).*mapped(b);end
        if all(ismember({a,b},names))
            other=prior.carbon.(spec.paired_carbon{j}).*mapped(b);
            assert(all(abs(value-other)<=1e-13.*max(abs(value),realmin)), ...
                'organicInitial:ambiguousNitrogen','N/SOC and N/C disagree for %s.',key);
        end
        prior.nitrogen.(key)=flip(value);
    end
    for group={'MB','BD','ENZN'}
        g=group{1};parentC=prior.carbon.(g);parentN=flip(prior.nitrogen.([g '_N']));
        cf=phaseFractions(g,'');nf=phaseFractions(g,'_N',cf);
        for j=1:3
            key=[g '_' spec.phases{j}];
            prior.carbon.(key)=parentC.*cf(:,j);
            prior.nitrogen.([key '_N'])=flip(parentN.*nf(:,j));
        end
    end
    for material={'litter','straw','manure'}
        g=material{1};prior.aom.(g)=aom(['aom_' g]);
        if strcmp(g,'litter')
            componentNames={'metabolic','cellulose','lignin','cwd'};
            present=any(startsWith(names,'aom_litter_metabolic_')| ...
                startsWith(names,'aom_litter_cellulose_')|startsWith(names,'aom_litter_lignin_')| ...
                startsWith(names,'aom_litter_cwd_'));
            if present
                for c=componentNames
                    prior.aom.litter.components.(c{1})=aom(['aom_litter_' c{1}]);
                end
                prior.aom.litter=soil.litter.componentTotals(prior.aom.litter,true);
            end
        end
    end
    prior.carbon.POC=prior.carbon.POM1+prior.carbon.POM2;
    prior.carbon.EP=prior.carbon.ENZP1+prior.carbon.ENZP2;
    prior.nitrogen.PON=prior.nitrogen.PON1+prior.nitrogen.PON2;
    prior.nitrogen.EP_N=prior.nitrogen.ENZP1_N+prior.nitrogen.ENZP2_N;
    prior.nitrogen.NEC_C=zeros(n,1);prior.nitrogen.NEC_N=zeros(n,1);
    prior.nitrogen.IOM_C=flip(prior.carbon.IOM);
    prior.dom=struct('doc_total_g_cm3',flip(prior.carbon.DOC), ...
        'don_total_g_cm3',prior.nitrogen.DON);

    function value=mapped(name)
        v=tbl.(name);assert(isnumeric(v)&&all(isfinite(v))&&all(v>=0), ...
            'organicInitial:value','%s must be a finite nonnegative numeric column.',name);
        [found,index]=ismember(z,x);
        value=zeros(n,1);value(found)=v(index(found));
        if any(~found)
            if numel(x)==1,value(~found)=v(1);
            else,value(~found)=interp1(x,v,min(max(z(~found),x(1)),x(end)),'linear');end
        end
        assert(all(isfinite(value))&&all(value>=0),'organicInitial:interpolation','Invalid mapped %s.',name);
    end
    function f=phaseFractions(g,suffix,default)
        cols=strcat(['frac_' g '_'],spec.phases,suffix);has=ismember(cols,names);
        if ~any(has)
            if nargin<3,f=repmat([mobileFraction,1-mobileFraction,0],n,1);else,f=default;end
            return;
        end
        assert(all(has),'organicInitial:phaseColumns','Provide all three %s%s phase fractions.',g,suffix);
        f=cell2mat(cellfun(@mapped,cols,'UniformOutput',false));
        assert(all(f<=1,'all')&&all(abs(sum(f,2)-1)<=32*eps), ...
            'organicInitial:phaseClosure','Parent-relative phase fractions must sum to one.');
    end
    function state=aom(prefix)
        state=struct('surface_C',0,'surface_N',0,'element_C',zeros(n,1),'element_N',zeros(n,1));
        for element={'C','N'}
            e=element{1};density=[prefix '_' e '_g_cm3'];ratio=[prefix '_' e '_SOC'];
            assert(~all(ismember({density,ratio},names)),'organicInitial:aomUnits','Use AOM density or ratio, not both.');
            if ismember(density,names),state.(['element_' e])=flip(mapped(density));
            elseif ismember(ratio,names),state.(['element_' e])=flip(socTop(:).*mapped(ratio));end
            surface=[prefix '_surface_' e '_g_cm2'];
            if ismember(surface,names)
                v=tbl.(surface);assert(isnumeric(v)&&all(isfinite(v))&&all(v>=0)&&all(v(2:end)==0), ...
                    'organicInitial:aomSurface','Surface AOM is stored only in the first row.');
                state.(['surface_' e])=v(1);
            end
        end
    end
end
