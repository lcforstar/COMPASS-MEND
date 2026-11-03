function [candidate,organic,dom,mend,feedback]=advanceSourceCandidate( ...
        accepted,organic,dom,theta,bd,dp,base,dt_s,temperature,head,freeDIN,remaining,scale)
%ADVANCESOURCECANDIDATE Adapt sole COMPASS owners to a paired source-hour map.
% Accepted C is top-down mg/g; organic N and total/sorbed DOM bottom-up g/cm3.
% Free DIN and root-reserved inventory are top-down mg/g. One oxygen factor
% per layer scales the COMPLETE conservative map increment, not its rates.
% Includes all organic pools, overflow and intracellular enzyme C/N. Input
% management and DOM transport have already been applied once by the caller.
% COMPASS has no NO2 state: only resolved NH4/NO3 enter enzyme allocation.
% Returns detached candidates and actual DIN/CO2 ledgers; no state commits.
% Source and interface differences: MEND_SOURCE_PRODUCTION_CONTRACT.md.
% Editor: Jizixing Zhu; Date: 2026-09-14.
% Change: 2026-09-16 Jizixing Zhu -- linear-equilibrium mode derives both
% phases from total/Kd, disables source adsorption rates, and returns only
% total DOM. Uptake uses the start-of-reaction aqueous inventory in the
% retained explicit hourly map; all later phase reads re-equilibrate.
% Change: 2026-10-12 Wenjie Liu -- pass separate nodal allocation ceilings
% into both N uptake stages, preserving physical and kinetic DIN inventories.
% Change: 2026-10-12 Wenjie Liu -- retain applied intracellular-enzyme gross
% C transfers by cell phase for passive N tracing, without endpoint inference.
% Editor: Wenjie Liu; Date: 2026-10-17.
% Change: 2026-10-17 Wenjie Liu -- carry configured C:N through source map,
%   oxygen scaling, intracellular phase transfers and the passive N ledger.
% Change: 2026-10-27 Wenjie Liu -- pass the environmental layer dormancy exponent.
% Change: 2026-10-31 Wenjie Liu -- expand scalar/layer sorption references and
%   multiply Kba and Kdes elementwise, retaining the source temperature laws.
    bd=bd(:);bt=flip(bd);scale=scale(:);n=numel(bd);
    assert(isequal(size(scale),size(bt))&&all(isfinite(scale))&&all(scale>=0&scale<=1), ...
        'soil:microbe:sourceScale','Oxygen reaction factor must be NL-by-1 in [0,1].');
    assert(accepted.schema_version>=4&&organic.schema_version>=5 ...
        &&ismember(string(dp.sorption_model),["mend_dynamic","linear_equilibrium"]) ...
        &&all(accepted.QOC==0)&&all(organic.QON==0), ...
        'soil:microbe:sourceOwners','Source-hour owner schemas/sole DOM ownership are required.');
    cp=soil.dom.phaseState(dom,'doc',theta,bd,dp);np=soil.dom.phaseState(dom,'don',theta,bd,dp);
    source={'POMO','POMH','MOM','MBA','MBD','ENZP1','ENZP2','ENZM','ENZNmt'};
    carbon={'POM1','POM2','MOC','MB','BD','ENZP1','ENZP2','EM','ENZN'};
    nitrogen={'PON1','PON2','MON','MB_N','BD_N','ENZP1_N','ENZP2_N','EM_N','ENZN_N'};
    c=struct();nn=struct();
    for j=1:numel(source)
        c.([source{j} '_C'])=accepted.(carbon{j});
        nn.([source{j} '_N'])=flip(organic.(nitrogen{j}))./bt.*1000;
    end
    c.CO2=accepted.IC;c.DOMT_C=flip(cp.free_g_cm3_bulk)./bt.*1000;
    c.QOM_C=flip(cp.sorbed_g_cm3_bulk)./bt.*1000;
    nn.DOMT_N=flip(np.free_g_cm3_bulk)./bt.*1000;
    nn.QOM_N=flip(np.sorbed_g_cm3_bulk)./bt.*1000;
    e=soil.microbe.mendEnvironmentalParameters(base,flip(temperature(:)), ...
        flip(head(:)).*9.80665e-5,base.pH_top_down);
    ref=base.hourly;t=flip(temperature(:));r=8.31446261815324;
    qmax=layerSorption(ref.Qmax,n,'Qmax');
    kba=layerSorption(ref.Kba,n,'Kba');
    kdes=layerSorption(ref.Kdes_ref,n,'Kdes_ref');
    ar=@(ea)exp(ea*1000/r.*(1/293.15-1./(t+273.15)));
    p=struct('VP1',e.VdPOM_h(:,1),'VP2',e.VdPOM_h(:,2), ...
        'KP1',e.KsPOM_mg_g(:,1),'KP2',e.KsPOM_mg_g(:,2), ...
        'VM',e.VdMOM_h,'KM',e.KsMOM_mg_g,'Vg',e.Vg_h,'Vm',e.Vm_h, ...
        'Vm_dorm',e.VmD_h,'KD',e.KsDOM_mg_g,'Yg',e.Yg, ...
        'Kads',kba.*kdes.*ar(5),'Kdes',kdes.*ar(20),'Qmax',qmax, ...
        'rMORT',e.rMORT_h,'pEP',base.pENZP,'pEM',base.pENZM, ...
        'rEP1',base.rENZP_h(1),'rEP2',base.rENZP_h(2),'rEM',base.rENZM_h, ...
        'VmA2D',e.VmA2D_h,'VmD2A',e.VmD2A_h,'fD',base.frPOM2DOM,'gD',base.frMB2DOM, ...
        'VNup_MB',e.VNup_MB_h,'KsNH4_MB',e.KsNH4_MB_mg_g,'KsNO3_MB',e.KsNO3_MB_mg_g, ...
        'VNup_VG',0,'wdorm',e.wdorm,'fpENZN',ref.fpENZN, ...
        'CN_ENZP',base.CN_ENZP(1),'CN_ENZM',base.CN_ENZM, ...
        'CN_MB_min',base.CN_MB_min,'CN_MB_max',base.CN_MB_max);
    equilibrium=strcmp(dp.sorption_model,'linear_equilibrium');
    if equilibrium
        % No second kinetic sorption process. QOM here is a temporary
        % algebraic partition, not another prognostic organic-carbon owner.
        p.Kads=zeros(n,1);p.Kdes=zeros(n,1);
    end
    z=zeros(n,1);addC=struct('POMinp1',z,'POMinp2',z,'DOMinp',z);
    addN=addC;addN.NH4_dep=z;addN.NO3_dep=z;
    din=struct('NH4',freeDIN(:,1),'NO3',freeDIN(:,2),'NO2',z);
    inventory=struct('NH4',remaining.nh4_mg_g(:),'NO3',remaining.no3_mg_g(:));
    if isfield(remaining,'nh4_accessible_mg_g')||isfield(remaining,'no3_accessible_mg_g')
        assert(all(isfield(remaining,{'nh4_accessible_mg_g','no3_accessible_mg_g'})), ...
            'soil:microbe:sourceAccessibleDIN','Require both accessible DIN species.');
        inventory.accessible_NH4=remaining.nh4_accessible_mg_g(:);
        inventory.accessible_NO3=remaining.no3_accessible_mg_g(:);
    end
    [cc,newN,~,l]=soil.microbe.mendSourceHour(c,nn,din,p,addC,addN,ref.iHR,dt_s,inventory);
    for key=fieldnames(c)',k=key{1};cc.(k)=c.(k)+scale.*(cc.(k)-c.(k));end
    for key=fieldnames(nn)',k=key{1};newN.(k)=nn.(k)+scale.*(newN.(k)-nn.(k));end
    % All ledger values other than the dimensionless N response are amounts.
    phi=l.phi_n;
    for key=fieldnames(l)'
        k=key{1};
        if isstruct(l.(k))
            for name=fieldnames(l.(k))',f=name{1};l.(k).(f)=scale.*l.(k).(f);end
        else,l.(k)=scale.*l.(k);end
    end
    l.phi_n=phi;f=l.carbon_flux;nf=l.nitrogen_flux;
    s=accepted;
    for j=1:numel(source)
        s.(carbon{j})=cc.([source{j} '_C']);
        organic.(nitrogen{j})=flip(newN.([source{j} '_N']).*bt./1000);
    end
    s.IC=cc.CO2;s.POC=s.POM1+s.POM2;s.EP=s.ENZP1+s.ENZP2;
    s.Tot=s.POC+s.MOC+s.MB+s.BD+s.EP+s.EM+s.ENZN;
    organic.PON=organic.PON1+organic.PON2;organic.EP_N=organic.ENZP1_N+organic.ENZP2_N;
    if ~equilibrium
        dom.doc_sorbed_g_cm3=flip(cc.QOM_C.*bt./1000);
        dom.don_sorbed_g_cm3=flip(newN.QOM_N.*bt./1000);
    end
    dom.doc_total_g_cm3=flip((cc.DOMT_C+cc.QOM_C).*bt./1000);
    dom.don_total_g_cm3=flip((newN.DOMT_N+newN.QOM_N).*bt./1000);
    beforeEnzymeC=cc.MBA_C+l.enzyme_production_C-l.enzyme_turnover_C;
    beforeEnzymeN=newN.MBA_N+(l.enzyme_production_C-l.enzyme_turnover_C)./p.CN_ENZM;
    phaseC=phaseLedger(c.MBA_C,c.MBD_C,beforeEnzymeC,cc.MBD_C,f.MBA2MBD,f.MBD2MBA,f.CO2_maintn_MBD);
    phaseN=phaseLedger(nn.MBA_N,nn.MBD_N,beforeEnzymeN,newN.MBD_N,nf.MBA2MBD,nf.MBD2MBA,z);
    s=soil.microbe.carryReactionPhases(s,{phaseC},'C',bd);
    organic=soil.microbe.carryReactionPhases(organic,{phaseN},'N',bd);
    % Apply enzyme synthesis/turnover to the SAME biomass phase in C and N.
    % Aggregate scaling would move returned intracellular matter between cell
    % phases whenever enzyme and active-biomass phase shares differ.
    % These ledger amounts are already oxygen-scaled: NL-by-3 top-down
    % mg C/g dry soil, columns mobile/attached/retained, for this candidate.
    l.enzyme_production_by_phase_C=zeros(n,3);
    l.enzyme_turnover_by_phase_C=zeros(n,3);
    phases={'mobile','attached','retained'};
    for phaseIndex=1:numel(phases)
        k=phases{phaseIndex};share=zeros(n,1);active=beforeEnzymeC>0;
        share(active)=s.(['MB_' k])(active)./beforeEnzymeC(active);
        fallback=~active&accepted.MB>0;
        share(fallback)=accepted.(['MB_' k])(fallback)./accepted.MB(fallback);
        if strcmp(k,'attached'),share(~active&~fallback)=1;end
        production=l.enzyme_production_C.*share;
        turnover=accepted.(['ENZN_' k]).*scale.*base.rENZM_h;
        l.enzyme_production_by_phase_C(:,phaseIndex)=production;
        l.enzyme_turnover_by_phase_C(:,phaseIndex)=turnover;
        s.(['MB_' k])=s.(['MB_' k])-production+turnover;
        organic.(['MB_' k '_N'])=organic.(['MB_' k '_N']) ...
            +flip((turnover-production).*bt./(1000.*p.CN_ENZM));
        s.(['ENZN_' k])=accepted.(['ENZN_' k])+production-turnover;
        organic.(['ENZN_' k '_N'])=flip(s.(['ENZN_' k]).*bt./(1000.*p.CN_ENZM));
        assert(all(s.(['MB_' k])>=0)&&all(organic.(['MB_' k '_N'])>=0), ...
            'soil:microbe:sourceEnzymePhaseOverdraw','Enzyme C/N allocation overdraws its cell phase.');
    end
    assert(all(abs(s.ENZN-s.ENZN_mobile-s.ENZN_attached-s.ENZN_retained)<1e-11), ...
        'soil:microbe:sourceEnzymePhases','Intracellular enzyme phase C does not close.');
    % Dimensionless metadata is added AFTER all oxygen amount scaling.
    l.enzyme_cn=p.CN_ENZM;
    secondRequest=zeros(n,2);hasSecond=l.second_uptake>0;
    secondRequest(hasSecond,1)=l.second_deficit(hasSecond).*l.second_NH4(hasSecond)./l.second_uptake(hasSecond);
    secondRequest(hasSecond,2)=l.second_deficit(hasSecond).*l.second_NO3(hasSecond)./l.second_uptake(hasSecond);
    first=allocationStage([l.raw_first_NH4,l.raw_first_NO3],l.first_NH4,l.first_NO3);
    second=allocationStage(secondRequest,l.second_NH4,l.second_NO3);
    release=l.Nmn_MBA+l.Nmn_MBD;actual=l.Nim_NH4+l.Nim_NO3;
    unmet=max(l.raw_first_NH4-l.first_NH4,0)+max(l.raw_first_NO3-l.first_NO3,0) ...
        +max(l.second_deficit-l.second_uptake,0);
    ma=struct('first',first,'second',second,'nh4_actual_mg_g',l.Nim_NH4, ...
        'no3_actual_mg_g',l.Nim_NO3,'actual_mg_g',actual,'requested_mg_g',actual+unmet, ...
        'unmet_mg_g',unmet,'nh4_release_mg_g',release,'net_nh4_sink_mg_g',l.Nim_NH4-release, ...
        'net_no3_sink_mg_g',l.Nim_NO3,'policy','root_first_source_hour_free_kinetics_total_inventory');
    assert(all(ma.net_nh4_sink_mg_g<=inventory.NH4+1e-12) ...
        &&all(ma.net_no3_sink_mg_g<=inventory.NO3+1e-12), ...
        'soil:microbe:sourceInventory','Source candidate exceeded post-root DIN.');
    if isfield(inventory,'accessible_NH4')
        assert(all(ma.net_nh4_sink_mg_g<=inventory.accessible_NH4+1e-12) ...
            &&all(ma.net_no3_sink_mg_g<=inventory.accessible_NO3+1e-12), ...
            'soil:microbe:sourceAccessibleDIN','Source candidate exceeded nodal DIN allowance.');
    end
    mk=struct('phi_n_limitation',phi,'stoichiometric_demand_mg_g',l.second_deficit, ...
        'raw_kinetic_request_mg_g',[l.raw_first_NH4,l.raw_first_NO3], ...
        'requested_mg_g',[l.raw_first_NH4,l.raw_first_NO3]);
    mend=struct('din_allocation',ma,'din_kinetic_request',mk,'dt_s',dt_s, ...
        'nitrogen_residual_mg_g',l.N_residual,'source_hour_ledger',l, ...
        'substep_coupling_policy','paired_source_hour_then_common_oxygen_increment_factor');
    overflow=l.CO2_ovflow_MBA+l.CO2_ovflow_MBD;
    feedback=struct('overflow_co2_mg_g',overflow,'overflow_mba_mg_g',l.CO2_ovflow_MBA, ...
        'overflow_mbd_mg_g',l.CO2_ovflow_MBD);
    sourceFlux={'POMdec1','POMdec2','POMdec2DOM1','POMdec2DOM2','POMdec2MOM1','POMdec2MOM2', ...
        'MOM2DOM','DOM2MBA','CO2_growth','CO2_maintn_MBA','CO2_maintn_MBD','MBA_mortality', ...
        'MBA2DOM','MBA2POM1','MBA2POM2','MBA2EP1','MBA2EP2','MBA2EM', ...
        'EP2DOM1','EP2DOM2','EM2DOM','MBA2MBD','MBD2MBA'};
    targetFlux={'POM1_decomp','POM2_decomp','POM1_to_DOM','POM2_to_DOM','POM1_to_MOM','POM2_to_MOM', ...
        'MOM_to_DOM','DOM_to_MBA','CO2_growth','CO2_maintenance_MBA','CO2_maintenance_MBD','MBA_mortality', ...
        'MBA_to_DOM','MBA_to_POM1','MBA_to_POM2','MBA_to_ENZP1','MBA_to_ENZP2','MBA_to_ENZM', ...
        'ENZP1_to_DOM','ENZP2_to_DOM','ENZM_to_DOM','MBA_to_MBD','MBD_to_MBA'};
    flux=struct();for j=1:numel(targetFlux),flux.(targetFlux{j})=f.(sourceFlux{j});end
    % Preserve every accepted source-hour carbon route needed to diagnose
    % dissolved-carbon partitioning, intracellular enzyme turnover, and the
    % complete microbial respiration budget.
    % Change: 2026-09-18 Wenjie Liu -- expose complete MEND C flux ledger.
    flux.DOM_to_QOM=f.DOM2QOM;
    flux.QOM_to_DOM=f.QOM2DOM;
    flux.DOM_sorption_net=f.DOM2QOMnet;
    flux.ENZN_production=l.enzyme_production_C;
    flux.ENZN_to_MBA=l.enzyme_turnover_C;
    flux.N_overflow_CO2=overflow;
    rc=zeros(n,1);rn=rc;
    for key=fieldnames(c)',k=key{1};rc=rc+cc.(k)-c.(k);end
    for key=fieldnames(nn)',k=key{1};rn=rn+newN.(k)-nn.(k);end
    rn=rn-ma.net_nh4_sink_mg_g-ma.net_no3_sink_mg_g;
    assert(all(abs(rc)<1e-10)&&all(abs(rn)<1e-10), ...
        'soil:microbe:sourceCoupledClosure','Coupled source C/N increments failed closure.');
    mend.nitrogen_residual_mg_g=rn;
    if equilibrium
        % Validate the new algebraic phases after reaction/O2 scaling and
        % expose their redistribution separately from source kinetic fluxes.
        afterC=soil.dom.phaseState(dom,'doc',theta,bd,dp);
        afterN=soil.dom.phaseState(dom,'don',theta,bd,dp);
        mend.equilibrium_sorption=struct('policy','linear_equilibrium_stage_split', ...
            'doc_sorbed_change_g_cm3',afterC.sorbed_g_cm3_bulk-cp.sorbed_g_cm3_bulk, ...
            'don_sorbed_change_g_cm3',afterN.sorbed_g_cm3_bulk-np.sorbed_g_cm3_bulk);
        % Linear-equilibrium sorption is an algebraic repartition rather than
        % a pair of kinetic gross rates. Report its accepted net redistribution
        % while retaining zero gross adsorption/desorption diagnostics.
        flux.DOM_to_QOM(:)=0;
        flux.QOM_to_DOM(:)=0;
        flux.DOM_sorption_net=flip(afterC.sorbed_g_cm3_bulk- ...
            cp.sorbed_g_cm3_bulk)./bt.*1000;
    end
    candidate=struct('state',s,'doc_total_g_cm3',dom.doc_total_g_cm3, ...
        'co2_c_rate_ug_cm3_s',flip((s.IC-accepted.IC).*bt.*1000./dt_s), ...
        'carbon_residual_mg_g',rc,'integrated_flux_mg_g',flux, ...
        'substep_count',1,'rejected_substeps',0,'min_substep_s',dt_s, ...
        'rate_time_basis','hour','source_model','MEND','source_commit',base.source_commit, ...
        'source_hourly',true,'oxygen_factor_top_down',scale);
end

function value=layerSorption(value,n,name)
% Read-only references use source top-down mg/g and hourly units.
% Editor: Wenjie Liu; Date: 2026-10-31.
    assert(isnumeric(value)&&isreal(value) ...
        &&(isscalar(value)||isequal(size(value),[n,1])) ...
        &&all(isfinite(value),'all'), ...
        'soil:microbe:sourceSorptionShape','%s must be finite scalar or NL-by-1.',name);
    if isscalar(value),value=repmat(value,n,1);end
end

function item=phaseLedger(a,b,aa,bb,toDormant,toActive,loss)
% Source's signed MBA compensation becomes a reverse phase transfer.
% Editor: Jizixing Zhu; Date: 2026-09-14.
    item=struct('before',struct('MB',a,'BD',b),'after',struct('MB',aa,'BD',bb), ...
        'a_to_d',max(toDormant,0),'d_to_a',toActive+max(-toDormant,0),'dormant_loss',loss);
end

function s=allocationStage(request,nh4,no3)
% Match existing accepted diagnostic species ledgers, in top-down mg N/g.
% Editor: Jizixing Zhu; Date: 2026-09-14.
    s=struct('requested_mg_g',request,'nh4_actual_mg_g',nh4,'no3_actual_mg_g',no3, ...
        'actual_mg_g',nh4+no3,'unmet_by_species_mg_g',max(request-[nh4,no3],0));
end
