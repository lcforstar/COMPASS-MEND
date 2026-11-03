function report=verifyMendDynamicSorption(outputBase)
%VERIFYMENDDYNAMICSORPTION Exercise actual source-hour D/Q adapter contracts.
% Initialize the configured site once, then call the production adapter on
% detached one-hour probes. D/Q/DON/QON are separate donor inventories;
% site inputs, DIN chemistry, accepted state and model clock are untouched.
% Oracle: MEND-main MOD_MEND.F90 1189-1197,1229-1238,2243-2258 and
% MEND_mult MOD_MEND.F90 796-808,828-837,1190-1205 (local source archives).
% Tests target gross exchange, order, units, grid direction, Qmax response,
% temperature, dry/wet phase reads, paired C/N and caller immutability.
% Full independent Fortran trajectory parity is verifyMendSourceMatlab;
% this focused test does not relabel its own oracle as that trajectory test.
% Editor: Wenjie Liu; Date: 2026-10-23.
% Change: 2026-10-31 Wenjie Liu -- make sorption references scalar/profile
% safe and test heterogeneous profiles and single-layer response locality.
    root=fileparts(fileparts(fileparts(mfilename('fullpath'))));
    if nargin<1||isempty(outputBase)
        outputBase=fullfile(root,'data','verification','mend_dq_restoration_20261023_01');
    end
    oldPath=path;oldDirectory=pwd;
    cleanup=onCleanup(@()restoreProcess(oldPath,oldDirectory)); %#ok<NASGU>
    cd(root);addpath(root,fullfile(root,'src'));
    if ~isfolder(outputBase),mkdir(outputBase);end
    target=fullfile(outputBase,'mend_dynamic_sorption.json');
    report=struct('status','RUNNING','site_inputs_modified',false, ...
        'DIN_chemistry_modified',false,'model_time_advanced',false, ...
        'independent_Fortran_trajectory_test','run separately: verifyMendSourceMatlab');
    writeReport(target,report);
    try
        model=init.initialize(fullfile(root,'data','sites','changwu'));
        b=model.built;before=b;nl=b.domain.NL;p=b.parameters.soil_microbes;
        assert(strcmp(p.reaction_method,'source_hourly') ...
            &&strcmp(b.parameters.soil_dom.sorption_model,'mend_dynamic'), ...
            'verifyMendDynamic:mode','Expected configured dynamic source-hour site.');
        assert(b.control.switches.lPlantCN==0, ...
            'verifyMendDynamic:plantCN','PlantCN must remain disabled.');
        context=sourceContext(b);
        actual=runAdapter(context,p);
        report.site_state=checkResult(context,p,actual,'actual initialized site');

        % Alternating D-rich/Q-poor and D-poor/Q-rich layers ensure both
        % net directions occur and a reversed grid cannot pass silently.
        d=linspace(.12,.20,nl)';q=linspace(.02,.08,nl)';
        d(1:2:end)=.001;q(1:2:end)=linspace(.35,.45,numel(q(1:2:end)))';
        probe=setDOM(context,d,q,d./10,q./12);
        probe.temperature=20*ones(nl,1);
        baseline=runAdapter(probe,p);
        report.mixed_donors=checkResult(probe,p,baseline,'alternating D/Q donors');
        f=baseline.microbes.integrated_flux_mg_g;
        assert(all(f.DOM_to_QOM>0)&&all(f.QOM_to_DOM>0) ...
            &&any(f.DOM_sorption_net>0)&&any(f.DOM_sorption_net<0), ...
            'verifyMendDynamic:exchange','Probes must exercise both gross and net directions.');

        % Dynamic D/Q is independent of an equilibrium Kd. Holding actual
        % D/Q, temperature and water potential fixed also holds dry-soil
        % source kinetics fixed when liquid volume changes.
        changed=probe;
        changed.dom_parameters.kd_doc_cm3_g=37*changed.dom_parameters.kd_doc_cm3_g;
        changed.dom_parameters.kd_don_cm3_g=19*changed.dom_parameters.kd_don_cm3_g;
        assert(isequaln(baseline,runAdapter(changed,p)), ...
            'verifyMendDynamic:KdRepartition','Kd must not repartition a dynamic reaction state.');
        changed=probe;changed.theta=.4*probe.theta;
        assert(isequaln(baseline,runAdapter(changed,p)), ...
            'verifyMendDynamic:thetaRepartition','Water-volume change must not reinitialize D/Q.');
        report.Kd_and_water_volume_do_not_repartition_DQ=true;

        larger=p;larger.hourly.Qmax=2.*p.hourly.Qmax;
        roomy=runAdapter(probe,larger);
        report.double_Qmax=checkResult(probe,larger,roomy,'double source capacity');
        g=roomy.microbes.integrated_flux_mg_g;
        assert(isequal(f.DOM_to_MBA,g.DOM_to_MBA) ...
            &&all(g.DOM_to_QOM>f.DOM_to_QOM)&&all(g.QOM_to_DOM<f.QOM_to_DOM), ...
            'verifyMendDynamic:capacity','Qmax must affect source sorption, not direct D uptake.');
        qmax=layerReference(p.hourly.Qmax,nl,'Qmax');
        full=setDOM(probe,.15*ones(nl,1),qmax, ...
            .015*ones(nl,1),qmax./12);
        saturated=runAdapter(full,p);
        report.at_Qmax=checkResult(full,p,saturated,'Q at source capacity');
        assert(max(abs(saturated.microbes.integrated_flux_mg_g.DOM_to_QOM))<1e-14 ...
            &&all(saturated.microbes.integrated_flux_mg_g.QOM_to_DOM>0), ...
            'verifyMendDynamic:capacityLimit','At Qmax adsorption vanishes but desorption persists.');

        % Legacy scalars and an exactly uniform NL-column have the same
        % source map, even when the configured site uses layered inputs.
        refs={'Qmax','Kba','Kdes_ref'};scalar=p;uniform=p;
        for j=1:numel(refs)
            key=refs{j};value=layerReference(p.hourly.(key),nl,key);
            scalar.hourly.(key)=value(1);
            uniform.hourly.(key)=repmat(value(1),nl,1);
        end
        assert(isequaln(runAdapter(probe,scalar),runAdapter(probe,uniform)), ...
            'verifyMendDynamic:scalarProfileParity', ...
            'Scalar and constant-column sorption references must produce identical outputs.');
        report.scalar_uniform_profile_bitwise_parity=true;

        % All three profiles vary independently in source top-down order.
        % The original Fortran algebra remains the independent oracle.
        varied=p;
        varied.hourly.Qmax=qmax.*linspace(.8,1.3,nl)';
        varied.hourly.Kba=layerReference(p.hourly.Kba,nl,'Kba').*linspace(.6,1.4,nl)';
        varied.hourly.Kdes_ref=layerReference(p.hourly.Kdes_ref,nl,'Kdes_ref') ...
            .*linspace(.7,1.3,nl)';
        heterogeneous=runAdapter(probe,varied);
        report.heterogeneous_sorption=checkResult(probe,varied,heterogeneous, ...
            'heterogeneous Qmax/Kba/Kdes profiles');
        hf=heterogeneous.microbes.integrated_flux_mg_g;
        selected=max(1,floor(nl/2));other=true(nl,1);other(selected)=false;
        report.single_layer_response=struct();
        for j=1:numel(refs)
            key=refs{j};changed=varied;
            changed.hourly.(key)(selected)=1.5.*changed.hourly.(key)(selected);
            local=runAdapter(probe,changed);
            report.single_layer_response.(key)=checkResult(probe,changed,local, ...
                ['single-layer ' key ' response']);
            lf=local.microbes.integrated_flux_mg_g;
            assert(isequal(hf.DOM_to_MBA,lf.DOM_to_MBA) ...
                &&isequal(hf.DOM_to_QOM(other),lf.DOM_to_QOM(other)) ...
                &&isequal(hf.QOM_to_DOM(other),lf.QOM_to_DOM(other)) ...
                &&(~isequal(hf.DOM_to_QOM(selected),lf.DOM_to_QOM(selected)) ...
                    ||~isequal(hf.QOM_to_DOM(selected),lf.QOM_to_DOM(selected))), ...
                'verifyMendDynamic:layerLocality', ...
                'Changing %s at one source layer must affect sorption only at that layer.',key);
        end
        report.single_layer_response.source_top_down_layer=selected;

        thermal=cell(2,1);temperatures=[10 30];
        for j=1:2
            changed=probe;changed.temperature=temperatures(j)*ones(nl,1);
            thermal{j}=runAdapter(changed,p);
            report.(['temperature_' num2str(temperatures(j)) '_C'])= ...
                checkResult(changed,p,thermal{j},'source Arrhenius sorption');
        end
        assert(all(thermal{2}.microbes.integrated_flux_mg_g.QOM_to_DOM ...
                   >thermal{1}.microbes.integrated_flux_mg_g.QOM_to_DOM), ...
            'verifyMendDynamic:temperature','Source desorption must respond to temperature.');

        % A full Q donor cannot feed uptake before the original explicit
        % hour reaches desorption. It may supply D for the following hour.
        emptyD=setDOM(probe,zeros(nl,1),.3*ones(nl,1),zeros(nl,1),.025*ones(nl,1));
        delayed=runAdapter(emptyD,p);
        report.zero_D=checkResult(emptyD,p,delayed,'zero initial D with nonzero Q');
        fd=delayed.microbes.integrated_flux_mg_g;
        assert(all(fd.DOM_to_MBA==0)&&all(fd.DOM_to_QOM==0)&&all(fd.QOM_to_DOM>0), ...
            'verifyMendDynamic:hourOrder','Desorbed/new D must not enter earlier source uptake.');
        afterPhase=soil.dom.phaseState(delayed.dom,'doc',emptyD.theta, ...
            emptyD.bulk_density,emptyD.dom_parameters);
        assert(all(afterPhase.free_g_cm3_bulk>0), ...
            'verifyMendDynamic:nextHourSubstrate','Desorption must create next-hour free D.');
        assert(isequaln(b,before)&&isequaln(model.built,before), ...
            'verifyMendDynamic:mutation','Detached probes mutated the initialized model.');
        report.status='PASS';report.layers=nl;report.PlantCN_enabled=false;
        report.production_kernel_changed=false;report.production_adapter_changed=false;
        report.D_and_Q_independent=true;report.original_hour_order_preserved=true;
        report.gross_exchange_and_capacity_response=true;report.temperature_response=true;
        report.accepted_inputs_unchanged=true;
        writeReport(target,report);
        fprintf('Dynamic MEND D/Q: source adapter, gross exchange, order and C/N PASS.\n%s\n',target);
    catch err
        report.status='FAIL';report.failure=struct('identifier',err.identifier, ...
            'message',err.message,'stack',err.stack);writeReport(target,report);rethrow(err);
    end
end

function context=sourceContext(b)
    bd=b.parameters.soil.bulk_density(:);rhoKd=b.parameters.soil_n.RHOKG;
    theta=b.Y.soil_nitrogen.theta_liquid_prev;
    [nh4,no3]=soil.organic_n.availableDinMgG(b.Y.nh4,b.Y.no3,theta,rhoKd,bd);
    remaining=struct('nh4_mg_g',nh4,'no3_mg_g',no3);
    reservation=struct('nh4_mg_g',zeros(size(nh4)),'no3_mg_g',zeros(size(no3)));
    remaining=soil.nitrogen.limitMendDinInventoryAtNodes(remaining,struct('enabled',false), ...
        b.Y.nh4,b.Y.no3,theta,rhoKd,bd,b.domain.dz(:),3600);
    free=soil.nitrogen.mendKineticDinAfterRoot(remaining,reservation, ...
        b.Y.nh4,b.Y.no3,theta,rhoKd,bd);
    context=struct('Y',b.Y,'theta',mean(b.workspace.derived.Theta_LL,2), ...
        'bulk_density',bd,'dom_parameters',b.parameters.soil_dom,'dt_s',3600, ...
        'temperature',.5*(b.Y.T(1:end-1)+b.Y.T(2:end)), ...
        'head',.5*(b.Y.h(1:end-1)+b.Y.h(2:end)), ...
        'free_DIN',free,'remaining',remaining,'scale',ones(b.domain.NL,1));
end

function c=setDOM(c,d,q,nD,nQ)
% Explicit dry-soil mg/g probes -> independent bottom-up bulk inventories.
    rho=flip(c.bulk_density);
    c.Y.soil_dom.doc_total_g_cm3=flip((d+q).*rho/1000);
    c.Y.soil_dom.doc_sorbed_g_cm3=flip(q.*rho/1000);
    c.Y.soil_dom.don_total_g_cm3=flip((nD+nQ).*rho/1000);
    c.Y.soil_dom.don_sorbed_g_cm3=flip(nQ.*rho/1000);
end

function r=runAdapter(c,p)
    before={c,p};
    [microbes,organic,dom,mend,feedback]=soil.microbe.advanceSourceCandidate( ...
        c.Y.soil_microbes,c.Y.soil_organic,c.Y.soil_dom,c.theta,c.bulk_density, ...
        c.dom_parameters,p,c.dt_s,c.temperature,c.head,c.free_DIN,c.remaining,c.scale);
    assert(isequaln(before,{c,p}),'verifyMendDynamic:mutation','Adapter changed inputs.');
    r=struct('microbes',microbes,'organic',organic,'dom',dom,'mend',mend,'feedback',feedback);
end

function summary=checkResult(c,p,r,label)
% Independent algebra from original Fortran: uptake first, then fAds, then
% paired donor transfer. It checks the complete adapter's units/grid/owner,
% not merely a direct invocation of its internal kernel.
    rho=flip(c.bulk_density);d=flip(c.Y.soil_dom.doc_total_g_cm3 ...
        -c.Y.soil_dom.doc_sorbed_g_cm3)./rho*1000;
    q=flip(c.Y.soil_dom.doc_sorbed_g_cm3)./rho*1000;
    nD=flip(c.Y.soil_dom.don_total_g_cm3-c.Y.soil_dom.don_sorbed_g_cm3)./rho*1000;
    nQ=flip(c.Y.soil_dom.don_sorbed_g_cm3)./rho*1000;
    e=soil.microbe.mendEnvironmentalParameters(p,flip(c.temperature), ...
        flip(c.head)*9.80665e-5,p.pH_top_down);
    uptake=min((e.Vg_h+e.Vm_h)./e.Yg.*d.*c.Y.soil_microbes.MB ...
        ./(e.KsDOM_mg_g+d),d);
    rgas=8.31446261815324;t=flip(c.temperature);
    qmax=layerReference(p.hourly.Qmax,numel(d),'Qmax');
    kba=layerReference(p.hourly.Kba,numel(d),'Kba');
    kdes=layerReference(p.hourly.Kdes_ref,numel(d),'Kdes_ref');
    ads=kba.*kdes ...
        .*exp(5000/rgas*(1/293.15-1./(t+273.15))) ...
        .*(d-uptake).*(1-q./qmax);
    des=kdes.*exp(20000/rgas*(1/293.15-1./(t+273.15))) ...
        .*q./qmax;
    boundedDes=des>q+ads;des(boundedDes)=q(boundedDes)+ads(boundedDes);
    boundedAds=~boundedDes&ads>d-uptake+des;
    ads(boundedAds)=d(boundedAds)-uptake(boundedAds)+des(boundedAds);
    f=r.microbes.integrated_flux_mg_g;
    errors=[maximum(f.DOM_to_MBA-uptake),maximum(f.DOM_to_QOM-ads), ...
        maximum(f.QOM_to_DOM-des),maximum(f.DOM_sorption_net-(ads-des)), ...
        maximum(flip(r.dom.doc_sorbed_g_cm3)./rho*1000-(q+ads-des))];
    assert(max(errors)<1e-12,'verifyMendDynamic:sourceAlgebra', ...
        '%s differs from original Fortran uptake/fAds/order.',label);
    ratioD=zeros(size(d));ratioQ=ratioD;
    ratioD(d>0)=nD(d>0)./d(d>0);ratioQ(q>0)=nQ(q>0)./q(q>0);
    qNexpected=nQ+ads.*ratioD-des.*ratioQ;
    pairedNerror=maximum(flip(r.dom.don_sorbed_g_cm3)./rho*1000-qNexpected);
    assert(pairedNerror<1e-12,'verifyMendDynamic:pairedDonorN', ...
        'Original D and Q donor C:N must accompany sorption transfer.');
    [cb,nb]=reactiveInventory(c.Y.soil_microbes,c.Y.soil_organic,c.Y.soil_dom,c.bulk_density);
    [ca,na]=reactiveInventory(r.microbes.state,r.organic,r.dom,c.bulk_density);
    respiration=r.microbes.co2_c_rate_ug_cm3_s*c.dt_s*1e-6;
    din=r.mend.din_allocation;
    incomingN=flip(din.net_nh4_sink_mg_g+din.net_no3_sink_mg_g).*c.bulk_density/1000;
    balanceC=maximum(ca-cb+respiration);balanceN=maximum(na-nb-incomingN);
    assert(balanceC<1e-12&&balanceN<1e-12,'verifyMendDynamic:closure', ...
        'Adapter paired C/N mass balance failed.');
    assert(~isfield(r.mend,'equilibrium_sorption') ...
        &&all(r.dom.doc_sorbed_g_cm3>=0) ...
        &&all(r.dom.doc_sorbed_g_cm3<=r.dom.doc_total_g_cm3) ...
        &&all(r.dom.don_sorbed_g_cm3>=0) ...
        &&all(r.dom.don_sorbed_g_cm3<=r.dom.don_total_g_cm3), ...
        'verifyMendDynamic:owners','Dynamic phases must persist and close without equilibrium overwrite.');
    summary=struct('status','PASS','max_source_algebra_error_mg_g',max(errors), ...
        'max_paired_sorbed_N_error_mg_g',pairedNerror, ...
        'max_C_balance_g_cm3',balanceC,'max_N_balance_g_cm3',balanceN, ...
        'gross_adsorption_range_mg_g',[min(ads),max(ads)], ...
        'gross_desorption_range_mg_g',[min(des),max(des)], ...
        'net_sorption_range_mg_g',[min(ads-des),max(ads-des)]);
end

function value=layerReference(value,n,name)
% Expand only a finite real scalar or a top-down NL-by-1 reference.
% Editor: Wenjie Liu; Date: 2026-10-31.
    assert(isnumeric(value)&&isreal(value) ...
        &&(isscalar(value)||isequal(size(value),[n,1])) ...
        &&all(isfinite(value),'all'), ...
        'verifyMendDynamic:referenceShape','%s must be scalar or NL-by-1.',name);
    if isscalar(value),value=repmat(value,n,1);end
end

function [carbon,nitrogen]=reactiveInventory(m,n,d,rho)
    carbon=zeros(size(rho));nitrogen=zeros(size(rho));
    cfields={'POM1','POM2','MOC','MB','BD','ENZP1','ENZP2','EM','ENZN'};
    nfields={'PON1','PON2','MON','MB_N','BD_N','ENZP1_N','ENZP2_N','EM_N','ENZN_N'};
    for j=1:numel(cfields)
        carbon=carbon+flip(m.(cfields{j})).*rho/1000;
        nitrogen=nitrogen+n.(nfields{j});
    end
    carbon=carbon+d.doc_total_g_cm3;nitrogen=nitrogen+d.don_total_g_cm3;
end

function x=maximum(x)
    x=max(abs(x),[],'all');
end

function writeReport(file,report)
    fid=fopen(file,'w');assert(fid>=0,'verifyMendDynamic:output','Cannot write %s.',file);
    closer=onCleanup(@()fclose(fid)); %#ok<NASGU>
    fprintf(fid,'%s\n',jsonencode(report,PrettyPrint=true));
end

function restoreProcess(oldPath,oldDirectory)
    path(oldPath);cd(oldDirectory);
end
