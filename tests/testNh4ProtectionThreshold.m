function report=testNh4ProtectionThreshold(outputRoot)
%TESTNH4PROTECTIONTHRESHOLD Focused independent NH4 protection contracts.
% Initializes the real site without advancing its clock. Synthetic probes
% test units/grid direction, exact default/zero arithmetic, the actual ADE
% source assembly at a protected neighbor, signed deficits, gross recycling,
% invalid configurations, checkpoint policy and the existing snapshot restore.
% Coupled windows/resume are verified separately by the scripts harness.
% No production source or site input is changed and no MATLAB is spawned.
% Editor: Wenjie Liu; Date: 2026-11-01.
    root=fileparts(fileparts(mfilename('fullpath')));
    if nargin<1,outputRoot=fullfile(root,'data','output','changwu','nh4_threshold_unit');end
    if ~isfolder(outputRoot),mkdir(outputRoot);end
    oldPath=path;oldDirectory=pwd;
    cleanup=onCleanup(@()restoreProcess(oldPath,oldDirectory)); %#ok<NASGU>
    cd(root);addpath(root,fullfile(root,'src'));
    report=struct('status','RUNNING','model_clock_advanced',false,'site_inputs_modified',false);
    target=fullfile(outputRoot,'unit_receipt.json');writeReport(target,report);
    try
        model=init.initialize(fullfile(root,'data','sites','changwu'));
        b=model.built;cfg=model.raw.config;
        report.units_and_projection=checkUnits(b,cfg);
        report.zero_parity_and_phase=checkZeroAndPhase();
        report.node_bounds=checkNodeBounds();
        report.gross_recycling=checkGrossRecycling();
        report.configuration=checkConfiguration(model);
        report.checkpoints=checkCheckpoints(b,cfg);
        report.snapshot_restore=checkRestore(b);
        report.status='PASS';writeReport(target,report);
        fprintf('NH4 protection focused contracts PASS.\n%s\n',target);
    catch err
        report.status='FAIL';report.failure=struct('identifier',err.identifier, ...
            'message',err.message,'stack',err.stack);writeReport(target,report);rethrow(err);
    end
end

function result=checkUnits(b,cfg)
    s=b.parameters.soil;
    fields={'bulk_density','kg','disp','dif_nh4','dif_no3','k_mineralization', ...
        'vmax_nit','km_nit','vmax_denit','km_denit','k_nh4_dissolve','kurea', ...
        'ph','msoc_volume','field_capacity','clay'};
    for j=1:numel(fields),key=fields{j};s.(key)=s.(key)(1:2);end
    s.bulk_density=[1.1;1.6];s.clay=[.301;.1985];b.parameters.soil=s;
    b.domain.NL=2;b.domain.NN=3;b.domain.dz=[7;11];
    cfg.nitrogen=struct('nh4_protected_per_clay_mg_kg',13.95);
    p=init.buildNitrogenSoilParams(b,cfg);
    expectedS=[4.19895;2.769075];expectedM=[4.618845;4.43052];
    assert(max(abs(p.S_min_nh4_mg_kg-expectedS))<2e-14, ...
        'testNh4Protection:clayUnits','Mapped clay must not be divided by 100 twice.');
    assert(max(abs(p.M_min_nit_ug_cm3-expectedM))<2e-14, ...
        'testNh4Protection:volumeUnits','mg/kg to ug/cm3 uses final rho without a factor 1000.');
    theta=.2*ones(2,2);nh4=20*ones(3,1);no3=zeros(3,1);dz=[7;11];
    [a4,a3]=soil.organic_n.availableDinMgG(nh4,no3,theta,p.RHOKG,s.bulk_density);
    remaining=struct('nh4_mg_g',a4,'no3_mg_g',a3);
    [out,audit]=soil.nitrogen.limitMendDinInventoryAtNodes(remaining,struct('enabled',false), ...
        nh4,no3,theta,p.RHOKG,s.bulk_density,dz,3600,p.S_min_nh4_mg_kg);
    expectedNode=[4.618845*7/2;4.618845*7/2+4.43052*11/2;4.43052*11/2];
    assert(max(abs(audit.nh4.protected_node_mass_ug_cm2-expectedNode))<5e-13);
    assert(isequaln(out.nh4_mg_g,remaining.nh4_mg_g)&&isequaln(out.no3_mg_g,remaining.no3_mg_g));
    assert(all(out.nh4_accessible_mg_g<=max(a4-1e-3*flip(expectedS),0)+1e-15));
    result=struct('layer_mg_kg',p.S_min_nh4_mg_kg,'layer_ug_cm3',p.M_min_nit_ug_cm3, ...
        'node_protected_ug_cm2',expectedNode,'physical_inventory_unchanged',true);
end

function result=checkZeroAndPhase()
    bd=[1.2;1.5];kd=bd*3.5;theta=[.13 .13;.27 .27];dz=[10;20];
    nh4=[3;8;12];no3=[1;2;4];phys=equations.physicalConstants();
    args={20*ones(3,1),nh4,2,theta,[.3;.3],[.5;.5],kd,[1e-5;2e-5], ...
        [10;20],zeros(2,2),phys.R,phys.Tr,ones(2,2),0,[7;7],zeros(2,2), ...
        1,dz,ones(2,2),false,-100*ones(3,1)};
    a=cell(1,8);z=a;v=a;
    [a{:}]=soil.nitrogen.chemistry.nitrification(args{:});
    [z{:}]=soil.nitrogen.chemistry.nitrification(args{:},0);
    [v{:}]=soil.nitrogen.chemistry.nitrification(args{:},zeros(2,1));
    assert(isequaln(a,z)&&isequaln(a,v),'testNh4Protection:nitZero','All nitrification outputs must be bitwise identical.');
    [v{:}]=soil.nitrogen.chemistry.nitrification(args{:},100*ones(2,1));
    assert(all(v{1}==0,'all')&&all(v{4}==0,'all')&&all(v{7}==0,'all'));
    assert(isequaln(a{5},v{5}),'testNh4Protection:totalDiagnostic','Protection must not redefine M_total.');
    [t4,t3]=soil.organic_n.availableDinMgG(nh4,no3,theta,kd,bd);
    remaining=struct('nh4_mg_g',t4,'no3_mg_g',t3);
    reservation=struct('nh4_mg_g',zeros(2,1),'no3_mg_g',zeros(2,1));
    root=struct('enabled',false);
    [r,audit]=soil.nitrogen.limitMendDinInventoryAtNodes(remaining,root,nh4,no3,theta,kd,bd,dz,3600);
    [rz,az]=soil.nitrogen.limitMendDinInventoryAtNodes(remaining,root,nh4,no3,theta,kd,bd,dz,3600,0);
    [rv,av]=soil.nitrogen.limitMendDinInventoryAtNodes(remaining,root,nh4,no3,theta,kd,bd,dz,3600,zeros(2,1));
    assert(isequaln({r,audit},{rz,az})&&isequaln({r,audit},{rv,av}));
    free=soil.nitrogen.mendKineticDinAfterRoot(r,reservation,nh4,no3,theta,kd,bd);
    free0=soil.nitrogen.mendKineticDinAfterRoot(r,reservation,nh4,no3,theta,kd,bd,0);
    freev=soil.nitrogen.mendKineticDinAfterRoot(r,reservation,nh4,no3,theta,kd,bd,zeros(2,1));
    assert(isequaln(free,free0)&&isequaln(free,freev));
    % Independent dry-soil arithmetic, including heterogeneous layer direction.
    expectedTotal=[5.52*10/1.5;4.33*5.5/1.2]*1e-3;
    expectedFree=[.27*10/1.5;.13*5.5/1.2]*1e-3;
    assert(max(abs(t4-expectedTotal))<1e-15&&max(abs(free(:,1)-expectedFree))<1e-15);
    threshold=[4;6];before=r;
    protected=soil.nitrogen.mendKineticDinAfterRoot(r,reservation,nh4,no3,theta,kd,bd,threshold);
    expected=max(expectedTotal-1e-3*flip(threshold),0).*(expectedFree./expectedTotal);
    assert(max(abs(protected(:,1)-expected))<1e-15&&isequaln(r,before));
    assert(isequaln(protected(:,2),free(:,2)));
    result=struct('nitrification_all_eight_outputs_bitwise',true,'MEND_helper_outputs_bitwise',true, ...
        'liquid_fraction_retained',true,'physical_inventory_unchanged',true);
end

function result=checkNodeBounds()
    theta=[.25 .25];rhoKd=1.75;dz=1;dt=1;c=[.5;10];floorMass=[.5;.5];
    rate=[0 .3];zero=[0 0];
    old=soil.nitrogen.limitNitrificationAtNodes(rate,zero,c,theta,rhoKd,dz,dt);
    before=transportReaction(c,theta,rhoKd,dz,old,zero,dt);
    assert(before(1)>0&&before(1)<.5,'testNh4Protection:probe','The unprotected FE counterexample must remain positive but cross its floor.');
    [limited,ledger]=soil.nitrogen.limitNitrificationAtNodes(rate,zero,c,theta,rhoKd,dz,dt,floorMass);
    after=transportReaction(c,theta,rhoKd,dz,limited,zero,dt);
    assert(all(limited==0,'all')&&all(after>=floorMass-1e-14));
    [byBudget,~]=soil.nitrogen.limitNitrificationByNodeBudget(rate,[0;9.5],dz);
    assert(isequaln(limited,byBudget));
    % Positive production first fills an initial protection deficit. Clamping
    % old excess before adding this source would incorrectly permit uptake.
    cLow=[.25;.25];rate=ones(1,2);smallSource=.3*ones(1,2);
    low=soil.nitrogen.limitNitrificationAtNodes(rate,smallSource,cLow,theta,rhoKd,dz,dt,floorMass);
    assert(all(low==0,'all'));
    smallEnd=transportReaction(cLow,theta,rhoKd,dz,low,smallSource,dt);
    assert(max(abs(smallEnd-[.4;.4]))<1e-14);
    largeSource=.8*ones(1,2);
    capped=soil.nitrogen.limitNitrificationAtNodes(rate,largeSource,cLow,theta,rhoKd,dz,dt,floorMass);
    largeEnd=transportReaction(cLow,theta,rhoKd,dz,capped,largeSource,dt);
    assert(all(capped>0,'all')&&all(largeEnd>=floorMass-1e-14));
    [base,info]=soil.nitrogen.limitNitrificationAtNodes(rate,zero,c,theta,rhoKd,dz,dt);
    [explicit,info0]=soil.nitrogen.limitNitrificationAtNodes(rate,zero,c,theta,rhoKd,dz,dt,0);
    assert(isequaln({base,info},{explicit,info0}));
    result=struct('old_neighbor_still_positive_but_below_floor',before(1), ...
        'protected_neighbor',after(1),'small_source_endpoint',smallEnd, ...
        'large_source_endpoint',largeEnd,'node_budget',ledger.node_available_ug_cm2_s);
end

function after=transportReaction(c,theta,rhoKd,dz,rate,other,dt)
% Exercise the real lumped-storage / consistent-source ADE assembly.
    nl=numel(dz);nn=nl+1;z=zeros(nl,2);r=1+rhoKd./theta;
    [c1,c1o,c2,c2o,c3,c4,c3a,c4a,c5,c6,m,k,~,~,~,f]= ...
        soil.nitrogen.transport.assembleMatrixCoefficients(r,r,z,z,z,other-rate, ...
        ones(nl,2),ones(nl,2),z,z,z,z,dz,nl,nn,zeros(nl,1),zeros(nl,1),theta,theta);
    [rhs,rhs3,c3,c4]=soil.nitrogen.transport.assembleEquation( ...
        c1,c1o,c2,c2o,c3,c4,c5,c6,nl,nn,dt,c,zeros(nn,1),m,k,f,zeros(3,3,5));
    [after,~]=soil.nitrogen.transport.solveSoluteLinearSystem(c4,c3,c,zeros(nn,1),nn,nl,c4a,c3a,rhs,rhs3);
end

function result=checkGrossRecycling()
    pools={'POMO','POMH','MOM','QOM','DOMT','MBA','MBD','ENZP1','ENZP2','ENZM','ENZNmt'};
    c=struct();n=struct();
    for j=1:numel(pools),c.([pools{j} '_C'])=0;n.([pools{j} '_N'])=0;end
    c.CO2=0;c.POMO_C=1;c.POMH_C=1;c.MOM_C=1;c.MBA_C=1;c.MBD_C=1;
    n.POMO_N=.1;n.POMH_N=.1;n.MOM_N=.1;n.MBA_N=.02;n.MBD_N=.5;
    positive={'KP1','KP2','KM','KD','Qmax','KsNH4_MB','KsNO3_MB'};
    nonnegative={'VP1','VP2','VM','Vg','Vm','Vm_dorm','Kads','Kdes','rMORT', ...
        'pEP','pEM','rEP1','rEP2','rEM','VmA2D','VmD2A','VNup_MB','VNup_VG','fpENZN'};
    p=struct();for j=1:numel(positive),p.(positive{j})=1;end
    for j=1:numel(nonnegative),p.(nonnegative{j})=0;end
    p.Yg=.5;p.fD=.75;p.gD=.5;p.wdorm=1;
    p.CN_MB_min=5.5;p.CN_MB_max=14;p.CN_ENZP=3.2;p.CN_ENZM=3.2;
    din=struct('NH4',0,'NO3',0,'NO2',0);
    inventory=struct('NH4',.004,'NO3',0,'accessible_NH4',0,'accessible_NO3',0);
    addC=struct('POMinp1',0,'POMinp2',0,'DOMinp',0);
    addN=addC;addN.NH4_dep=0;addN.NO3_dep=0;
    [~,newN,newDin,l]=soil.microbe.mendSourceHour(c,n,din,p,addC,addN,1,3600,inventory);
    released=l.Nmn_MBA+l.Nmn_MBD;
    assert(l.Nim_NH4>0&&released>0&&l.first_NH4==0&&l.second_NH4>0);
    assert(l.Nim_NH4-released<=inventory.accessible_NH4+1e-14);
    independent=(sum(cell2mat(struct2cell(newN)))-sum(cell2mat(struct2cell(n)))) ...
        +newDin.NH4-inventory.NH4+newDin.NO3-inventory.NO3;
    assert(abs(independent)<1e-14&&abs(l.N_residual)<1e-14&&abs(l.C_residual)<1e-14);
    n.MBD_N=1/p.CN_MB_min;
    [~,~,~,pure]=soil.microbe.mendSourceHour(c,n,din,p,addC,addN,1,3600,inventory);
    assert(pure.Nim_NH4==0&&pure.Nmn_MBA==0&&pure.Nmn_MBD==0);
    result=struct('gross_NH4_uptake_mg_g',l.Nim_NH4,'same_hour_release_mg_g',released, ...
        'net_old_inventory_consumption_mg_g',l.Nim_NH4-released,'independent_N_residual',independent, ...
        'no_positive_source_gross_uptake_zero',true);
end

function result=checkConfiguration(model)
    b=model.built;raw=model.raw;cfg=raw.config;
    before=init.buildNitrogenSoilParams(b);
    zeroCfg=cfg;zeroCfg.nitrogen=struct('nh4_protected_per_clay_mg_kg',0);
    assert(isequaln(before,init.buildNitrogenSoilParams(b,zeroCfg)));
    zeroCfg.nitrogen.nh4_protection_scope='nitrification_only';
    assert(isequaln(before,init.buildNitrogenSoilParams(b,zeroCfg)));
    io.checkSite(raw,raw.input_plan);
    invalid={-1,NaN,Inf,[0 1],true,false,'0',"0",1+1i,[]};
    for j=1:numel(invalid)
        probe=cfg;probe.nitrogen=struct('nh4_protected_per_clay_mg_kg',invalid{j});
        expectError(@()init.buildNitrogenSoilParams(b,probe),'buildNitrogenSoilParams:nh4ProtectionValue');
        bad=raw;bad.config=probe;expectError(@()io.checkSite(bad,bad.input_plan),'checkSite:invalid');
    end
    probe=cfg;probe.nitrogen=struct('nh4_protected_per_clay_mg_kk',13.95);
    expectError(@()init.buildNitrogenSoilParams(b,probe),'buildNitrogenSoilParams:nitrogenKeys');
    bad=raw;bad.config=probe;expectError(@()io.checkSite(bad,bad.input_plan),'checkSite:invalid');
    probe=cfg;probe.nitrogen=struct('nh4_protected_per_clay_mg_kg',13.95);
    for invalidScope={true,12,'nit',["nitrification_only","nitrification_and_mend"]}
        probe.nitrogen.nh4_protection_scope=invalidScope{1};
        expectError(@()init.buildNitrogenSoilParams(b,probe),'buildNitrogenSoilParams:nh4ProtectionScope');
        bad=raw;bad.config=probe;expectError(@()io.checkSite(bad,bad.input_plan),'checkSite:invalid');
    end
    probe.nitrogen.nh4_protection_scope='nitrification_only';
    good=raw;good.config=probe;io.checkSite(good,good.input_plan);
    probe.microbe.reaction_method='rk4_equilibrium';
    expectError(@()init.buildNitrogenSoilParams(b,probe),'buildNitrogenSoilParams:nh4ProtectionMode');
    bad=raw;bad.config=probe;expectError(@()io.checkSite(bad,bad.input_plan),'checkSite:invalid');
    result=struct('invalid_scalar_cases',numel(invalid),'unknown_key_rejected',true, ...
        'unsupported_microbe_method_rejected_by_both_preflight_layers',true,'default_zero_parameters_bitwise',true);
end

function result=checkCheckpoints(b,cfg)
    base=b;base.parameters.soil_n=init.buildNitrogenSoilParams(b);
    signature0=io.checkpointCompatibilitySignature(base);
    assert(~isfield(signature0,'nh4_protection_threshold'));
    old=struct('schema_version',2,'compatibility',signature0);
    io.validateCheckpointCompatibility(base,old);
    cfg.nitrogen=struct('nh4_protected_per_clay_mg_kg',13.95);
    protected=b;protected.parameters.soil_n=init.buildNitrogenSoilParams(b,cfg);
    signature=io.checkpointCompatibilitySignature(protected);
    checkpoint=struct('schema_version',2,'compatibility',signature);
    io.validateCheckpointCompatibility(protected,checkpoint);
    assert(~isfield(signature.nh4_protection_threshold,'scope'));
    cfg.nitrogen.nh4_protection_scope='nitrification_and_mend';
    namedShared=b;namedShared.parameters.soil_n=init.buildNitrogenSoilParams(b,cfg);
    io.validateCheckpointCompatibility(namedShared,checkpoint);
    cfg.nitrogen.nh4_protection_scope='nitrification_only';
    nitOnly=b;nitOnly.parameters.soil_n=init.buildNitrogenSoilParams(b,cfg);
    nitCkpt=struct('schema_version',2,'compatibility',io.checkpointCompatibilitySignature(nitOnly));
    assert(strcmp(nitCkpt.compatibility.nh4_protection_threshold.scope,'nitrification_only'));
    io.validateCheckpointCompatibility(nitOnly,nitCkpt);
    expectError(@()io.validateCheckpointCompatibility(nitOnly,checkpoint),'validateCheckpointCompatibility:mismatch');
    expectError(@()io.validateCheckpointCompatibility(protected,nitCkpt),'validateCheckpointCompatibility:mismatch');
    expectError(@()io.validateCheckpointCompatibility(protected,old),'validateCheckpointCompatibility:mismatch');
    expectError(@()io.validateCheckpointCompatibility(base,checkpoint),'validateCheckpointCompatibility:mismatch');
    changed=protected;changed.parameters.soil_n.S_min_nh4_mg_kg(1)=1.01*changed.parameters.soil_n.S_min_nh4_mg_kg(1);
    expectError(@()io.validateCheckpointCompatibility(changed,checkpoint),'validateCheckpointCompatibility:mismatch');
    expectError(@()io.validateCheckpointCompatibility(protected,struct('schema_version',1)), ...
        'validateCheckpointCompatibility:legacyBgcCheckpoint');
    result=struct('old_schema2_to_zero_accepted',true,'same_nonzero_accepted',true, ...
        'zero_to_nonzero_rejected',true,'nonzero_to_zero_rejected',true, ...
        'changed_derived_profile_rejected',true,'schema1_nonzero_rejected',true);
end

function result=checkRestore(b)
    snapshot=solver.snapshotSubstepState(b);trial=b;
    trial.Y.nh4(:)=0;trial.Y.soil_nitrogen.flux.cum_cvSnit=123;
    trial.workspace.step.trial.nh4(:)=99;
    restored=solver.restoreSubstepState(trial,snapshot);
    assert(isequaln(restored.Y,b.Y)&&isequaln(restored.workspace.step,b.workspace.step));
    result=struct('full_accepted_Y_and_trial_restored_bitwise',true, ...
        'scope','existing snapshot/restore primitive; coupled reject path not injected here');
end

function expectError(action,identifier)
    caught=false;
    try,action();catch err
        assert(strcmp(err.identifier,identifier),'testNh4Protection:wrongError', ...
            'Expected %s, got %s: %s',identifier,err.identifier,err.message);caught=true;
    end
    assert(caught,'testNh4Protection:missingError','Expected named error %s.',identifier);
end

function writeReport(file,value)
    fid=fopen(file,'w','n','UTF-8');assert(fid>=0);cleanup=onCleanup(@()fclose(fid)); %#ok<NASGU>
    fprintf(fid,'%s\n',jsonencode(value,PrettyPrint=true));
end

function restoreProcess(oldPath,oldDirectory)
    path(oldPath);cd(oldDirectory);
end
