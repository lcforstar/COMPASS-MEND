function vector=spinupStateVector(Y,domain,parameters)
%SPINUPSTATEVECTOR Map accepted independent C/N stores to areal mass.
% Y is a complete accepted single-column state; domain.dz [cm,NL-by-1]
% and parameters.soil.bulk_density [g/cm3,NL-by-1] are bottom-up.
% Returns one row per independent pool and spatial location, in g C or N
% per cm2 ground. Layer concentration is multiplied by dz; top-down MEND
% mg C/g dry soil is flipped and multiplied by rho*dz/1000. DIN uses the
% accepted effective-liquid/ice snapshots and FE half-element capacities;
% urea and granules are bulk-volume stores. Plant g/m2 is divided by 1e4.
% Dynamic DOM and microbial phases are disjoint; their aggregate aliases
% are checked but never counted again. Respiratory IC and cumulative flux
% ledgers are excluded. Plant xsmr is signed C debt, included separately.
% This read-only adapter does not certify periodic equilibrium: forcing
% phase, full-window duration, water/heat, events and fluxes need own gates.
% Rejects unsupported schema/units/layout, invalid pools and broken aliases.
% Sources: init.buildSoilCarbon, init.buildSoilMicrobes, init.buildSoilDOM,
% init.buildMendHourlyState, init.buildInitialPlantPools,
% soil.nitrogen.transport.calculateStorage, soil.nitrogen.observationProfile.
% Editor: Wenjie Liu; Date: 2026-09-30.
% Change: 2026-09-30 Wenjie Liu -- withdraw the extra annual PlantCN pool
% schema; retain the original soil stock mapping and unknown-owner checks.
    dz=domain.dz(:); n=numel(dz);
    require(n>0&&all(isfinite(dz))&&all(dz>0) ...
        &&domain.NL==n,'geometry','Invalid bottom-up layer geometry.');
    rho=parameters.soil.bulk_density(:);
    require(isequal(size(rho),[n,1])&&all(isfinite(rho))&&all(rho>0), ...
        'geometry','Invalid layer dry bulk density.');
    o=Y.soil_organic; m=Y.soil_microbes; d=Y.soil_dom;
    require(ismember(o.schema_version,[4 5])&&ismember(m.schema_version,[3 4]) ...
        &&ismember(d.schema_version,[1 2]),'schema','Unsupported C/N state schema.');
    require(strcmp(o.grid_order,'bottom_up')&&strcmp(m.grid_order,'top_down') ...
        &&strcmp(o.volume_units,'g_per_cm3_bulk_soil') ...
        &&strcmp(o.surface_units,'g_per_cm2_ground') ...
        &&strcmp(m.units,'mg_C_g_dry_soil'),'units','Unsupported C/N units or orientation.');
    require(all(pool(m,'QOC',[n 1])==0)&&all(pool(o,'QON',[n 1])==0), ...
        'duplicateDOM','QOC/QON must not duplicate the sole DOM owner.');
    state_path=strings(0,1); element=strings(0,1); location=strings(0,1);
    element_index=zeros(0,1); endpoint_index=zeros(0,1);
    source_index=zeros(0,1); native_value=zeros(0,1);
    native_unit=strings(0,1); mass_g_cm2=zeros(0,1);
    function append(path,el,loc,v,mult,unit,idx,endpoint,source)
        count=numel(v); v=v(:); mult=mult(:);
        require(numel(mult)==count&&all(isfinite(mult))&&all(mult>=0), ...
            'capacity','Invalid stock conversion capacity.');
        state_path=[state_path;repmat(string(path),count,1)];
        element=[element;repmat(string(el),count,1)];
        location=[location;repmat(string(loc),count,1)];
        element_index=[element_index;idx(:)]; endpoint_index=[endpoint_index;endpoint(:)];
        source_index=[source_index;source(:)]; native_value=[native_value;v];
        native_unit=[native_unit;repmat(string(unit),count,1)];
        mass_g_cm2=[mass_g_cm2;v.*mult];
    end
    function layer(path,el,v,unit,microbe)
        idx=(1:n)'; source=idx; factor=dz;
        if microbe,v=flip(v);source=flip(idx);factor=rho.*dz/1000;end
        append(path,el,'layer_bottom_up',v,factor,unit,idx,zeros(n,1),source);
    end
    function aom(path,s)
        for el=["C","N"]
            key=['surface_' char(el)];v=pool(s,key,[1 1]);
            append(path+"."+key,el,'surface',v,1,'g/cm2',0,0,1);
            key=['element_' char(el)];v=pool(s,key,[n 1]);
            layer(path+"."+key,el,v,'g/cm3_bulk',false);
        end
    end
    require(isequal(sort(fieldnames(o.aom)),sort({'litter';'straw';'manure'})), ...
        'schema','Unknown AOM material owner.');
    for material={'litter','straw','manure'}
        k=material{1}; s=o.aom.(k); p="soil_organic.aom."+k;
        if isfield(s,'components')
            keys={'metabolic','cellulose','lignin','cwd'};
            require(isequal(sort(fieldnames(s.components)),sort(keys(:))), ...
                'schema','Unknown litter component owner.');
            for field={'surface_C','surface_N','element_C','element_N'}
                f=field{1}; shape=[n 1];if startsWith(f,'surface'),shape=[1 1];end
                total=zeros(shape);
                for j=1:numel(keys),total=total+pool(s.components.(keys{j}),f,shape);end
                alias(pool(s,f,shape),total,p+"."+f);
            end
            for j=1:numel(keys),aom(p+".components."+keys{j},s.components.(keys{j}));end
        else
            aom(p,s);
        end
    end
    layer('soil_organic.NEC_C','C',pool(o,'NEC_C',[n 1]),'g/cm3_bulk',false);
    layer('soil_organic.IOM_C','C',pool(o,'IOM_C',[n 1]),'g/cm3_bulk',false);
    layer('soil_organic.NEC_N','N',pool(o,'NEC_N',[n 1]),'g/cm3_bulk',false);
    sorption=string(parameters.soil_dom.sorption_model);
    require(isscalar(sorption)&&ismember(sorption,["linear_equilibrium","mend_dynamic"]), ...
        'schema','Unsupported DOM sorption owner.');
    for kind={'doc','don'}
        k=kind{1}; el='C';if strcmp(k,'don'),el='N';end
        total=pool(d,[k '_total_g_cm3'],[n 1]);hasPhase=isfield(d,[k '_sorbed_g_cm3']);
        require(hasPhase==(sorption=="mend_dynamic"),'schema','DOM phase schema disagrees with parameters.');
        if hasPhase
            sorbed=pool(d,[k '_sorbed_g_cm3'],[n 1]);free=total-sorbed;
            require(all(free>=0),'pool','Sorbed DOM exceeds its total owner.');
            layer("soil_dom."+k+"_free_g_cm3",el,free,'g/cm3_bulk',false);
            layer("soil_dom."+k+"_sorbed_g_cm3",el,sorbed,'g/cm3_bulk',false);
        else
            layer("soil_dom."+k+"_total_g_cm3",el,total,'g/cm3_bulk',false);
        end
    end
    carbon={'POM1','POM2','MOC','ENZP1','ENZP2','EM'};
    nitrogen={'PON1','PON2','MON','ENZP1_N','ENZP2_N','EM_N'};
    for j=1:numel(carbon)
        layer("soil_microbes."+carbon{j},'C',pool(m,carbon{j},[n 1]),'mg_C/g_dry_soil',true);
        layer("soil_organic."+nitrogen{j},'N',pool(o,nitrogen{j},[n 1]),'g/cm3_bulk',false);
    end
    alias(pool(m,'POC',[n 1]),m.POM1+m.POM2,'soil_microbes.POC');
    alias(pool(m,'EP',[n 1]),m.ENZP1+m.ENZP2,'soil_microbes.EP');
    alias(pool(o,'PON',[n 1]),o.PON1+o.PON2,'soil_organic.PON');
    alias(pool(o,'EP_N',[n 1]),o.ENZP1_N+o.ENZP2_N,'soil_organic.EP_N');
    groups={'MB','BD'};
    require(isfield(m,'ENZN')==isfield(o,'ENZN_N'),'schema','Unpaired intracellular enzyme C/N.');
    require(isfield(m,'ENZN')==(m.schema_version==4)&&isfield(o,'ENZN_N')==(o.schema_version==5), ...
        'schema','Intracellular enzyme owner disagrees with schema.');
    if isfield(m,'ENZN'),groups{end+1}='ENZN';end
    for j=1:numel(groups)
        g=groups{j};csum=zeros(n,1);nsum=csum;
        for phase={'mobile','attached','retained'}
            ckey=[g '_' phase{1}];nkey=[ckey '_N'];
            c=pool(m,ckey,[n 1]);nn=pool(o,nkey,[n 1]);csum=csum+c;nsum=nsum+nn;
            layer("soil_microbes."+ckey,'C',c,'mg_C/g_dry_soil',true);
            layer("soil_organic."+nkey,'N',nn,'g/cm3_bulk',false);
        end
        alias(pool(m,g,[n 1]),csum,"soil_microbes."+g);
        alias(pool(o,[g '_N'],[n 1]),nsum,"soil_organic."+g+"_N");
    end
    [~,~,din]=soil.nitrogen.observationProfile(Y,parameters.soil_n,dz);
    require(isequal(parameters.soil_n.RHO_bulk(:),rho), ...
        'geometry','DIN and organic dry bulk-density contracts differ.');
    nn=n+1;ni=(1:nn)';iz=zeros(nn,1);
    append('nh4','N','node_bottom_up',pool(Y,'nh4',[nn 1]), ...
        din.nh4_capacity_cm*1e-6,'ug/cm3_water',ni,iz,ni);
    append('no3','N','node_bottom_up',pool(Y,'no3',[nn 1]), ...
        din.no3_capacity_cm*1e-6,'ug/cm3_water',ni,iz,ni);
    for key={'urea','C_NH4_granule'}
        append(key{1},'N','node_bottom_up',pool(Y,key{1},[nn 1]), ...
            din.node_volume_cm*1e-6,'ug/cm3_bulk',ni,iz,ni);
    end
    ice=pool(Y.soil_nitrogen,'theta_ice_prev',[n 2]);
    for key={'C_NH4_ice','C_NO3_ice'}
        append(key{1},'N','layer_endpoint_bottom_up',pool(Y,key{1},[n 2]), ...
            ice.*(dz/2)*1e-6,'ug/cm3_ice',repmat((1:n)',2,1), ...
            [ones(n,1);2*ones(n,1)],(1:2*n)');
    end
    organs={'labile';'leaf';'root_fine';'livestem';'deadstem';'livecroot'; ...
        'deadcroot';'fruit';'root_storage';'stem_storage'};
    for el=["C","N"]
        p=Y.plant_cn.(el);keys=organs;if el=="C",keys{end+1}='xsmr';end
        require(isequal(sort(fieldnames(p)),sort(keys)),'schema','Unknown plant C/N pool owner.');
        for j=1:numel(keys)
            k=keys{j};v=pool(p,k,[1 1],strcmp(k,'xsmr'));
            if strcmp(k,'xsmr'),require(v<=0,'pool','xsmr must be nonpositive carbon debt.');end
            append("plant_cn."+el+"."+k,el,'plant',v,1e-4,'g/m2',0,0,1);
        end
        p=Y.plant_cn.transfer.(el);keys={'leaf';'root';'livestem'};
        require(isequal(sort(fieldnames(p)),sort(keys)),'schema','Unknown plant transfer owner.');
        for j=1:numel(keys)
            k=keys{j};v=pool(p,k,[1 1]);
            append("plant_cn.transfer."+el+"."+k,el,'plant',v,1e-4,'g/m2',0,0,1);
        end
    end
    require(all(isfinite(mass_g_cm2)),'pool','Canonical mass overflow.');
    canonical_unit="g_"+element+"/cm2";
    vector=table(state_path,element,location,element_index,endpoint_index, ...
        source_index,native_value,native_unit,mass_g_cm2,canonical_unit);
    vector.Properties.UserData=struct('schema_version',1,'dz_cm',dz, ...
        'rho_g_cm3',rho,'nh4_sorption_buffer',parameters.soil_n.RHOKG(:), ...
        'dom_sorption_model',char(sorption),'alias_validation_ulps',64, ...
        'scope','soil_organic_DOM_microbes_DIN_urea_granules_plant', ...
        'periodic_equilibrium_certified',false);
end

function value=pool(owner,key,shape,signed)
    if nargin<4,signed=false;end
    require(isstruct(owner)&&isfield(owner,key),'missingPool',"Missing pool "+key);
    value=owner.(key);
    require(isa(value,'double')&&isreal(value)&&isequal(size(value),shape) ...
        &&all(isfinite(value),'all')&&(signed||all(value>=0,'all')), ...
        'pool',"Invalid double-precision pool "+key);
end

function alias(a,b,name)
% Pure arithmetic consistency check, not a spin-up tolerance. 64 ulps covers
% these short (<=9 term) positive sums without a unit-dependent absolute floor.
    bound=64*eps(max(abs(a),abs(b)));
    require(all(abs(a-b)<=bound,'all'),'alias',"Inconsistent aggregate "+name);
end

function require(condition,id,message)
    assert(condition,['soil:organic:spinupState:' id],'%s',message);
end
