function [Y,audit]=applyAnchoredBiogeochemistryState(built,carryY,anchorY)
%APPLYANCHOREDBIOGEOCHEMISTRYSTATE 在轮间重置指定慢库，完整继承快库和水热。
% 输入为同一站点、网格、参数下的接受状态；built 提供固定容重、层厚与吸附容量。
% POM1/POM2/MOC 为自上向下 mg C/g；PON1/PON2/MON/IOM_C 为自下向上 g/cm3。
% 无机氮以初始接受态的有限元干土浓度锚定，按继承水分的容量重建水相浓度。
% 这是初始化投影，不推进时间，不是物理通量，不保证总 SOC 或总氮保持初始值。
% 原模型没有 IOM_N 独立库存；不得虚构该库。快库 C/N 及全部相分量逐值保留。
% Editor: Wenjie Liu; Date: 2026-10-23.
% Change: 2026-10-23 Wenjie Liu -- 慢库锚定、无机氮质量一致重建及非物理重置审计。
    sw=built.control.switches;
    assert(sw.lPlantCN==0&&sw.lLitterInput==0, ...
        'applyAnchoredBiogeochemistryState:scope','该锚定模式要求关闭整株碳氮及凋落输入。');
    assert(~isfield(sw,'lFertilizerSourceTracer')||sw.lFertilizerSourceTracer==0, ...
        'applyAnchoredBiogeochemistryState:tracer','锚定模式不允许未经同步投影的施肥来源标签。');
    assert(~isfield(carryY.soil_organic,'IOM_N')&&~isfield(anchorY.soil_organic,'IOM_N'), ...
        'applyAnchoredBiogeochemistryState:IOMN','发现新 IOM_N owner，必须先明确其锚定与守恒合同。');
    dz=built.domain.dz(:);rho=built.parameters.soil.bulk_density(:);n=numel(dz);
    assert(isequal(size(rho),[n,1])&&all(isfinite(rho))&&all(rho>0), ...
        'applyAnchoredBiogeochemistryState:geometry','容重与层厚不匹配。');
    before=diagnostics.biogeochemistryInventory(built,carryY);
    Y=carryY;reset=struct();
    for key={'POM1','POM2','MOC'}
        k=key{1};v=validPool(anchorY.soil_microbes,k,[n 1]);
        old=validPool(carryY.soil_microbes,k,[n 1]);Y.soil_microbes.(k)=v;
        reset.(k)=sum(flip(v-old).*rho.*dz)/1000;
    end
    for key={'PON1','PON2','MON','IOM_C'}
        k=key{1};v=validPool(anchorY.soil_organic,k,[n 1]);
        old=validPool(carryY.soil_organic,k,[n 1]);Y.soil_organic.(k)=v;
        reset.(k)=sum((v-old).*dz);
    end
    % 只重建兼容汇总，不改变任何快库或累计呼吸量 IC。
    m=Y.soil_microbes;m.POC=m.POM1+m.POM2;
    m.Tot=m.POC+m.MOC+m.QOC+m.MB+m.BD+m.ENZP1+m.ENZP2+m.EM;
    if isfield(m,'ENZN'),m.Tot=m.Tot+m.ENZN;end
    Y.soil_microbes=m;Y.soil_organic.PON=Y.soil_organic.PON1+Y.soil_organic.PON2;
    [targetNH4,targetNO3,targetCapacity]=soil.nitrogen.observationProfile( ...
        anchorY,built.parameters.soil_n,dz);
    [~,~,capacity]=soil.nitrogen.observationProfile(Y,built.parameters.soil_n,dz);
    assert(isequal(capacity.dry_mass_g_cm2,targetCapacity.dry_mass_g_cm2), ...
        'applyAnchoredBiogeochemistryState:dryMass','锚定与继承态必须具有相同干土质量。');
    Y.nh4=concentration(targetNH4.*capacity.dry_mass_g_cm2,capacity.nh4_capacity_cm,'NH4');
    Y.no3=concentration(targetNO3.*capacity.dry_mass_g_cm2,capacity.no3_capacity_cm,'NO3');
    for key={'urea','C_NH4_granule'}
        k=key{1};Y.(k)=validPool(anchorY,k,[n+1 1]);
        reset.(k)=sum((Y.(k)-carryY.(k)).*capacity.node_volume_cm)*1e-6;
    end
    % 冰中浓度也按初始质量重建；初始零冰氮保持零，不复制不相容浓度。
    ice=validPool(Y.soil_nitrogen,'theta_ice_prev',[n 2]);
    anchorIce=validPool(anchorY.soil_nitrogen,'theta_ice_prev',[n 2]);
    for key={'C_NH4_ice','C_NO3_ice'}
        k=key{1};mass=validPool(anchorY,k,[n 2]).*anchorIce;
        Y.(k)=concentration(mass,ice,k);
        reset.(k)=sum(sum((Y.(k).*ice-carryY.(k).*ice).*(dz/2)))*1e-6;
    end
    [nh4,no3]=soil.nitrogen.observationProfile(Y,built.parameters.soil_n,dz);
    roundtrip=[max(abs(nh4-targetNH4)),max(abs(no3-targetNO3))];
    assert(all(roundtrip<=128*eps(max(1,max([targetNH4;targetNO3])))), ...
        'applyAnchoredBiogeochemistryState:DINRoundtrip','干土基准无机氮锚定未达到机器精度。');
    [oldNH4,oldNO3]=dinStorage(built,carryY);[newNH4,newNO3]=dinStorage(built,Y);
    referenceShift=[newNH4-oldNH4,newNO3-oldNO3];
    reset.nh4_dissolved_adsorbed=referenceShift(1)*1e-6-reset.C_NH4_granule-reset.C_NH4_ice;
    reset.no3_dissolved=referenceShift(2)*1e-6-reset.C_NO3_ice;
    % cv_*_init 是库存参照，非物理累计通量。显式平移参照以隔离轮间重置。
    f=Y.soil_nitrogen.flux;oldFlux=f;
    for pair={{'nh4',newNH4,referenceShift(1)},{'no3',newNO3,referenceShift(2)}}
        p=pair{1};k=p{1};stock=p{2};delta=p{3};reference=['cv_' k '_init_mg_cm2'];
        if isfield(f,reference),f.(reference)=f.(reference)+delta;end
        f.(['cv' k])(:)=stock;f.(['cv' k 'OLD'])(:)=stock;
    end
    Y.soil_nitrogen.flux=f;
    for key=fieldnames(oldFlux)'
        k=key{1};if startsWith(k,'cum_')
            assert(isequaln(f.(k),oldFlux.(k)), ...
                'applyAnchoredBiogeochemistryState:flux','锚定不得修改累计物理通量 %s。',k);
        end
    end
    after=diagnostics.biogeochemistryInventory(built,Y);
    deltaC=reset.POM1+reset.POM2+reset.MOC+reset.IOM_C;
    deltaN=reset.PON1+reset.PON2+reset.MON+reset.nh4_dissolved_adsorbed+reset.no3_dissolved ...
        +reset.C_NH4_granule+reset.C_NH4_ice+reset.C_NO3_ice+reset.urea;
    assert(abs((after.physical_carbon_g_cm2-before.physical_carbon_g_cm2)-deltaC)<=1e-10 ...
        &&abs((after.nitrogen_g_cm2-before.nitrogen_g_cm2)-deltaN)<=1e-12, ...
        'applyAnchoredBiogeochemistryState:resetAudit','非物理重置差额与独立库存核算不一致。');
    % 逐字段核验只修改规定慢库、DIN 及它们的汇总和参照快照。
    permittedTop={'soil_microbes','soil_organic','soil_nitrogen','nh4','no3','urea', ...
        'C_NH4_granule','C_NH4_ice','C_NO3_ice'};
    unchangedFields(carryY,Y,permittedTop);
    unchangedFields(carryY.soil_microbes,Y.soil_microbes,{'POM1','POM2','MOC','POC','Tot'});
    unchangedFields(carryY.soil_organic,Y.soil_organic,{'PON1','PON2','MON','IOM_C','PON'});
    unchangedFields(carryY.soil_nitrogen,Y.soil_nitrogen,{'flux'});
    audit=struct('schema_version',1,'mode','anchored_slow_pools', ...
        'physical_flux_added',false,'cumulative_physical_fluxes_unchanged',true, ...
        'slow_pools_exact',true,'fast_pools_exact',true,'water_heat_exact',true, ...
        'DIN_anchor_basis','initial_accepted_FE_dry_soil_profile', ...
        'nh4_roundtrip_max_abs_mgkg',roundtrip(1),'no3_roundtrip_max_abs_mgkg',roundtrip(2), ...
        'nitrogen_reference_shift_ug_cm2',struct('nh4',referenceShift(1),'no3',referenceShift(2)), ...
        'carbon_reset_g_cm2',deltaC,'nitrogen_reset_g_cm2',deltaN, ...
        'pre_inventory',before,'post_inventory',after,'reset_by_pool_g_cm2',reset, ...
        'IOM_N_status','not_represented_in_current_model', ...
        'total_SOC_anchored',false,'total_N_anchored',false, ...
        'interpretation','指定慢库与初始矿物氮锚定；总SOC及总氮随快库变化，重置差额仅属初始化账目。');
end

function value=validPool(owner,key,shape)
    assert(isfield(owner,key),'applyAnchoredBiogeochemistryState:pool','缺少库存 %s。',key);
    value=owner.(key);
    assert(isnumeric(value)&&isreal(value)&&isequal(size(value),shape) ...
        &&all(isfinite(value(:)))&&all(value(:)>=0), ...
        'applyAnchoredBiogeochemistryState:pool','库存 %s 的形状或数值无效。',key);
end

function value=concentration(mass,capacity,label)
    assert(isequal(size(mass),size(capacity))&&all(isfinite(mass(:))) ...
        &&all(mass(:)>=0)&&all(isfinite(capacity(:)))&&all(capacity(:)>=0) ...
        &&all(mass(capacity==0)==0), ...
        'applyAnchoredBiogeochemistryState:dryCapacity','%s 的非零初始库存缺少可用相容量。',label);
    value=zeros(size(mass));wet=capacity>0;value(wet)=mass(wet)./capacity(wet);
end

function [nh4,no3]=dinStorage(built,Y)
    [nh4,no3]=soil.nitrogen.transport.calculateStorage( ...
        Y.soil_nitrogen.theta_liquid_prev,Y.nh4,Y.no3,built.parameters.soil_n.RHOKG, ...
        built.domain.dz,built.domain.NL,Y.soil_nitrogen.theta_ice_prev, ...
        Y.C_NH4_ice,Y.C_NO3_ice,Y.C_NH4_granule);
end

function unchangedFields(before,after,permitted)
    assert(isequal(fieldnames(before),fieldnames(after)), ...
        'applyAnchoredBiogeochemistryState:unexpectedOwner','锚定不得新增或删除接受态 owner。');
    for key=fieldnames(before)'
        k=key{1};if ~ismember(k,permitted)
            assert(isequaln(before.(k),after.(k)), ...
                'applyAnchoredBiogeochemistryState:unexpectedChange','锚定意外改变 %s。',k);
        end
    end
end
