function inventory=biogeochemistryInventory(built,Y)
%BIOGEOCHEMISTRYINVENTORY Disjoint accepted coupled C/N inventory [g cm^-2].
% built supplies immutable bottom-up dz [cm], bulk density [g cm^-3], and
% NH4 sorption capacity RHOKG [-]. Y is an accepted full state (defaults to
% built.Y). MEND carbon is top-down mg C/g; organic N/DOM bottom-up g/cm3;
% AOM includes surface g/cm2; PlantCN uses g/m2. Mineral N uses the exact
% accepted FE endpoint water memory, including ice, granules and urea.
% Sum physical pools once: compatibility aggregates, microbial phases and
% fertilizer-source labels are not extra stores. MEND IC is cumulative
% respired C, not stored C. Negative plant xsmr is debt, reported separately.
% Scope: one coupled field column, all organic/DOM/microbial/PlantCN owners
% enabled. Unknown owner schemas or inconsistent duplicate owners are errors.
% No state changes, time integration, component residuals or flux inversion.
% Editor: Wenjie Liu; Date: 2026-09-30.
    if nargin<2,Y=built.Y;end
    sw=built.control.switches;
    required={'lSoilOrganic','lSoilDOM','lExplicitSoilMicrobes','lPlantCN','lNefc'};
    for j=1:numel(required)
        assert(isfield(sw,required{j})&&sw.(required{j})==1, ...
            'biogeochemistryInventory:scope','Required coupled owner %s is off.',required{j});
    end
    dz=positiveColumn(built.domain.dz,built.domain.NL,'dz');n=numel(dz);
    rho=positiveColumn(built.parameters.soil.bulk_density,n,'rho');
    rhokg=pool(built.parameters.soil_n,'RHOKG',n);
    m=Y.soil_microbes;o=Y.soil_organic;d=Y.soil_dom;
    assert(m.schema_version>=4&&o.schema_version>=5 ...
        &&strcmp(m.grid_order,'top_down')&&strcmp(o.grid_order,'bottom_up'), ...
        'biogeochemistryInventory:schema','Source-hour MEND owner schemas are required.');
    assert(all(pool(m,'QOC',n)==0)&&all(pool(o,'QON',n)==0), ...
        'biogeochemistryInventory:duplicateDOM','DOM must have one physical owner.');
    c=struct();nn=struct();
    carbon={'POM1','POM2','MOC','MB','BD','ENZP1','ENZP2','EM','ENZN'};
    nitrogen={'PON1','PON2','MON','MB_N','BD_N','ENZP1_N','ENZP2_N','EM_N','ENZN_N'};
    for j=1:numel(carbon)
        c.(carbon{j})=sum(flip(pool(m,carbon{j},n)).*rho.*dz)/1000;
        nn.(nitrogen{j})=sum(pool(o,nitrogen{j},n).*dz);
    end
    for group={'MB','BD','ENZN'}
        key=group{1};phaseC=zeros(n,1);phaseN=phaseC;
        for phase={'mobile','attached','retained'}
            phaseC=phaseC+pool(m,[key '_' phase{1}],n);
            phaseN=phaseN+pool(o,[key '_' phase{1} '_N'],n);
        end
        same(phaseC,pool(m,key,n),[key ' C phases']);
        same(phaseN,pool(o,[key '_N'],n),[key ' N phases']);
    end
    same(pool(m,'POC',n),pool(m,'POM1',n)+pool(m,'POM2',n),'POC aggregate');
    same(pool(m,'EP',n),pool(m,'ENZP1',n)+pool(m,'ENZP2',n),'EP aggregate');
    same(pool(o,'PON',n),pool(o,'PON1',n)+pool(o,'PON2',n),'PON aggregate');
    same(pool(o,'EP_N',n),pool(o,'ENZP1_N',n)+pool(o,'ENZP2_N',n),'EP_N aggregate');
    c.IOM=sum(pool(o,'IOM_C',n).*dz);c.NEC=sum(pool(o,'NEC_C',n).*dz);
    nn.NEC=sum(pool(o,'NEC_N',n).*dz);
    c.DOM=sum(pool(d,'doc_total_g_cm3',n).*dz);
    nn.DOM=sum(pool(d,'don_total_g_cm3',n).*dz);
    for material={'litter','straw','manure'}
        key=material{1};a=o.aom.(key);
        if isfield(a,'components')
            a=soil.litter.componentTotals(a,true);
        end
        c.(['aom_' key '_surface'])=pool(a,'surface_C',1);
        c.(['aom_' key '_soil'])=sum(pool(a,'element_C',n).*dz);
        nn.(['aom_' key '_surface'])=pool(a,'surface_N',1);
        nn.(['aom_' key '_soil'])=sum(pool(a,'element_N',n).*dz);
    end
    theta=endpoint(Y.soil_nitrogen.theta_liquid_prev,n,'liquid');
    ice=endpoint(Y.soil_nitrogen.theta_ice_prev,n,'ice');
    nh4=pool(Y,'nh4',n+1);no3=pool(Y,'no3',n+1);
    nh4ice=endpoint(Y.C_NH4_ice,n,'NH4 ice');
    no3ice=endpoint(Y.C_NO3_ice,n,'NO3 ice');
    nn.nh4_dissolved_adsorbed=sum(dz.*((theta(:,1)+rhokg).*nh4(1:n) ...
        +(theta(:,2)+rhokg).*nh4(2:n+1)))/2e6;
    nn.no3_dissolved=sum(dz.*(theta(:,1).*no3(1:n)+theta(:,2).*no3(2:n+1)))/2e6;
    nn.nh4_ice=sum(dz.*sum(ice.*nh4ice,2))/2e6;
    nn.no3_ice=sum(dz.*sum(ice.*no3ice,2))/2e6;
    volume=([dz;0]+[0;dz])/2;
    nn.nh4_granule=sum(volume.*pool(Y,'C_NH4_granule',n+1))/1e6;
    nn.urea=sum(volume.*pool(Y,'urea',n+1))/1e6;
    organs={'labile','leaf','root_fine','livestem','deadstem','livecroot', ...
        'deadcroot','fruit','root_storage','stem_storage'};
    for j=1:numel(organs)
        key=organs{j};c.(['plant_' key])=pool(Y.plant_cn.C,key,1)/1e4;
        nn.(['plant_' key])=pool(Y.plant_cn.N,key,1)/1e4;
    end
    for organ={'leaf','root','livestem'}
        key=organ{1};c.(['transfer_' key])=pool(Y.plant_cn.transfer.C,key,1)/1e4;
        nn.(['transfer_' key])=pool(Y.plant_cn.transfer.N,key,1)/1e4;
    end
    debt=Y.plant_cn.C.xsmr;
    assert(isscalar(debt)&&isreal(debt)&&isfinite(debt), ...
        'biogeochemistryInventory:debt','Plant xsmr must be finite.');
    c.plant_xsmr_positive=max(debt,0)/1e4;
    inventory=struct('schema_version',1,'units','g element/cm2 ground', ...
        'carbon_pools_g_cm2',c,'nitrogen_pools_g_cm2',nn, ...
        'physical_carbon_g_cm2',sum(cell2mat(struct2cell(c))), ...
        'nitrogen_g_cm2',sum(cell2mat(struct2cell(nn))), ...
        'plant_carbon_debt_g_cm2',min(debt,0)/1e4, ...
        'mend_respired_carbon_g_cm2',sum(flip(pool(m,'IC',n)).*rho.*dz)/1000);
    inventory.accounting_carbon_g_cm2=inventory.physical_carbon_g_cm2 ...
        +inventory.plant_carbon_debt_g_cm2;
end

function value=pool(s,key,n)
    assert(isfield(s,key),'biogeochemistryInventory:missingPool','Missing owner %s.',key);
    value=s.(key);
    assert(isnumeric(value)&&isreal(value)&&numel(value)==n ...
        &&all(isfinite(value(:)))&&all(value(:)>=0), ...
        'biogeochemistryInventory:invalidPool','Invalid physical owner %s.',key);
    value=value(:);
end
function value=positiveColumn(value,n,key)
    value=pool(struct('value',value),'value',n);
    assert(all(value>0),'biogeochemistryInventory:grid','%s must be positive.',key);
end
function value=endpoint(value,n,key)
    assert(isequal(size(value),[n,2])&&isreal(value) ...
        &&all(isfinite(value(:)))&&all(value(:)>=0), ...
        'biogeochemistryInventory:endpoint','Invalid accepted %s endpoint memory.',key);
end
function same(a,b,key)
    assert(all(abs(a-b)<=128*eps(max(1,max(abs([a(:);b(:)]))))), ...
        'biogeochemistryInventory:aggregate','Inconsistent %s.',key);
end

% Change: 2026-09-30 Wenjie Liu -- withdraw this round's APSIM stem WSC
%   and standing-senescent-leaf extension; retain the original COMPASS litter
%   and CoLM crop C/N allocation/harvest contracts. Earlier Change records are historical.
