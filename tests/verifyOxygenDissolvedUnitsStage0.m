function report = verifyOxygenDissolvedUnitsStage0(outputRoot, baselineKernelPath)
%VERIFYOXYGENDISSOLVEDUNITSSTAGE0 Discriminate the dissolved-O2 unit repair.
% Analytic, no-model-clock fixtures exercise the production dilute-species
% kernel and O2 adapter. Henry pressure is evaluated independently in atm;
% the retained 82.06 cm3-atm gas constant's rounding is explicitly preserved.
% The frozen pre-repair kernel is renamed only in a new verification output
% directory. For identical fixed inputs, its y and full ledger must remain
% bitwise identical; this is not a claim of coupled annual-trajectory parity.
% No source/site parameter or production tolerance is modified. MATLAB is
% not spawned. Use a new output directory for every invocation.
% Editor: Wenjie Liu; Date: 2026-11-01.
    root=fileparts(fileparts(mfilename('fullpath')));
    if nargin<1||isempty(outputRoot)
        outputRoot=fullfile(root,'data','output','changwu', ...
            ['oxygen_units_stage0_unit_' char(datetime('now','Format','yyyyMMdd_HHmmss'))]);
    end
    if nargin<2||isempty(baselineKernelPath)
        baselineKernelPath=fullfile(root,'data','verification', ...
            'mend_o2_units_stage0_20261101_01','before','src', ...
            '+soil','+gas','+species','solveDryMoleFractionTransport.m');
    end
    assert(isfile(baselineKernelPath),'verifyOxygenStage0:missingBaseline', ...
        'The actual frozen pre-repair kernel is required; no parity check is skipped.');
    assert(~isfolder(outputRoot),'verifyOxygenStage0:noOverwrite', ...
        'Choose a new Stage0 unit output directory.');
    mkdir(outputRoot);
    oldPath=path;cleanup=onCleanup(@()restorePath(oldPath)); %#ok<NASGU>
    addpath(fullfile(root,'src'),'-begin');
    kernel=which('soil.gas.species.solveDryMoleFractionTransport');
    expectedKernel=fullfile(root,'src','+soil','+gas','+species','solveDryMoleFractionTransport.m');
    assert(strcmpi(strrep(kernel,'\','/'),strrep(expectedKernel,'\','/')), ...
        'verifyOxygenStage0:wrongSource','The kernel must resolve in this release.');
    p=soil.gas.species.oxygenProperties();
    ratio=p.gas_constant_cm3_atm_mol_K/(8.31446261815324*1e6/101325);
    report=struct('status','RUNNING','model_clock_advanced',false, ...
        'production_or_site_modified',false,'production_kernel',kernel, ...
        'baseline_kernel',baselineKernelPath,'baseline_sha256',fileHash(baselineKernelPath), ...
        'gas_constant_rounding_ratio',ratio,'stage','setup');
    receipt=fullfile(outputRoot,'oxygen_units_stage0_receipt.json');
    writeReceipt(receipt,report);
    try
        legacy=prepareLegacy(baselineKernelPath,outputRoot);
        report.stage='analytic_Henry_and_oxygen_limits';
        levels=[0,p.half_saturation_mol_L/(p.henry_reference_mol_L_atm*ratio),.2095];
        values=zeros(numel(levels),3);
        for j=1:numel(levels)
            phase=uniformPhase(25,101325,0,.20,.30);
            f=fixture(phase,phase,[levels(j);levels(j)],zeros(1,1),10,3600,levels(j),0);
            [y,c,a,b,ledger]=soil.gas.solveO2Transport(f.args{1:9});
            [~,endpointC]=inventory(phase,y,f.dz,p);
            checkClose(c,mean(endpointC,2),'henry', ...
                'Henry mol/L: current/expected must be 1, not 1000.');
            checkClose(y,f.args{1},'uniformState','A uniform no-reaction state must remain uniform.');
            checkClose(a,repmat(c./(c+p.half_saturation_mol_L),1,2), ...
                'aerobic','The O2 adapter must use the corrected mol/L concentration.');
            checkClose(a+b,ones(1,2),'complement','Oxygen factors must sum to one.');
            assert(all(a(:)>=0&a(:)<=1)&&all(b(:)>=0&b(:)<=1));
            assert(abs(ledger.balance_residual_mol_cm2)<1e-14,'verifyOxygenStage0:balance');
            values(j,:)=[c,a(1),b(1)];
        end
        assert(values(1,1)==0&&values(1,2)==0&&values(1,3)==1);
        checkClose(values(2,2),.5,'halfSaturation','C=20e-6 mol/L must give f_aerobic=0.5.');
        assert(all(diff(values(:,2))>0)&&all(diff(values(:,3))<0));
        report.air_and_limits=struct('y_levels',levels,'C_aerobic_anaerobic',values);

        report.stage='temperature_pressure_vapor_Henry';
        conditions=[0,90000,500;35,110000,4000];
        henryValues=zeros(2,1);
        for j=1:2
            phase=uniformPhase(conditions(j,1),conditions(j,2),conditions(j,3),.40,.05);
            f=fixture(phase,phase,[.2095;.2095],0,10,600,.2095,0);
            [y,c,ledger]=soil.gas.species.solveDryMoleFractionTransport(f.args{:},p);
            [expectedStorage,endpointC]=inventory(phase,y,f.dz,p);
            checkClose(c,mean(endpointC,2),'henry','Henry must use dry partial pressure and local temperature.');
            checkClose(ledger.storage_after_mol_cm2,expectedStorage, ...
                'phaseStorage','Gas plus dissolved inventory must use mol/cm3 soil.');
            henryValues(j)=c;
        end
        report.temperature_pressure_vapor=struct('conditions_C_Pa_PvPa',conditions, ...
            'C_mol_L',henryValues);

        report.stage='heterogeneous_phase_inventory_and_fixed_sink_parity';
        old=struct('theta_g',[.26 .25;.25 .23],'theta_l',[.24 .25;.25 .27], ...
            'temperature_C',[8;15;25],'porosity',[.5;.5], ...
            'gas_pressure_Pa',[99000;100500;101325],'vapor_pressure_Pa',[850;1300;2500]);
        new=struct('theta_g',[.22 .20;.20 .21],'theta_l',[.28 .30;.30 .29], ...
            'temperature_C',[12;18;24],'porosity',[.5;.5], ...
            'gas_pressure_Pa',[98000;100000;101000],'vapor_pressure_Pa',[900;1700;2400]);
        f=fixture(old,new,[.08;.13;.2095],[2e-11;3e-11],[7;13],1200,.2095,-2e-5);
        f.args{2}=[.10;.15;.2095];
        f.args{5}.free_dry_air_mass_flux_g_cm2_s=[3e-9;-1e-9];
        f.args{5}.gas_dispersion_cm2_s=[2e-5;5e-5];
        [y,c,ledger]=soil.gas.species.solveDryMoleFractionTransport(f.args{:},p);
        [yOld,cOld,ledgerOld]=legacy(f.args{:},p);
        assert(isequaln(y,yOld)&&isequaln(ledger,ledgerOld), ...
            'verifyOxygenStage0:fixedSinkParity', ...
            'Identical fixed inputs must preserve y and every discrete ledger field bitwise.');
        checkClose(c,cOld/1000,'outputRatio','Only the returned mol/L output must lose the extra factor.');
        [before,~,oldGas,oldLiquid]=inventory(old,f.args{1},f.dz,p);
        [after,endpointC,newGas,newLiquid]=inventory(new,y,f.dz,p);
        checkClose(c,mean(endpointC,2),'henry','Average endpoint Henry concentrations, not averaged drivers.');
        checkClose(ledger.storage_before_mol_cm2,before,'phaseStorage','Old phase inventory mismatch.');
        checkClose(ledger.storage_after_mol_cm2,after,'phaseStorage','New phase inventory mismatch.');
        bottomExpected=f.bottomLiquid*endpointC(1,1)/1000*f.dt;
        checkClose(ledger.bottom_input_mol_cm2,bottomExpected, ...
            'liquidOutflow','Liquid outflow carries C_liq/1000 mol/cm3 water.');
        referenceMean=mean([f.args{2}(1:2),f.args{2}(2:3)],2);
        newMean=mean([y(1:2),y(2:3)],2);
        reactionExpected=sum(f.dz.*f.args{6}.*newMean./referenceMean)*f.dt;
        checkClose(ledger.reaction_loss_mol_cm2,reactionExpected, ...
            'reactionLedger','Fixed-reference implicit reaction must match the applied layer removal.');
        assert(abs(ledger.balance_residual_mol_cm2)<1e-14,'verifyOxygenStage0:balance');
        assert(ledger.positivity_correction_mol_cm2==0&&min(y)>0&&max(y)<1);
        assert(ledger.reaction_loss_mol_cm2>0&&ledger.bottom_input_mol_cm2<0);
        report.fixed_sink=struct('bitwise_y_and_full_ledger',true, ...
            'gas_before_after_mol_cm2',[oldGas,newGas], ...
            'liquid_before_after_mol_cm2',[oldLiquid,newLiquid], ...
            'ledger',ledger,'new_C_mol_L',c,'legacy_C_mol_L',cOld);
        report.stage='complete';report.status='PASS';writeReceipt(receipt,report);
        fprintf('Stage0 dissolved-O2 analytic/fixed-input contracts PASS.\n%s\n',receipt);
    catch err
        report.status='FAIL';report.failure=struct('identifier',err.identifier, ...
            'message',err.message,'stack',err.stack);
        writeReceipt(receipt,report);rethrow(err);
    end
end

function phase=uniformPhase(T,Pg,Pv,gas,liquid)
    phase=struct('theta_g',[gas gas],'theta_l',[liquid liquid], ...
        'temperature_C',[T;T],'porosity',.5, ...
        'gas_pressure_Pa',[Pg;Pg],'vapor_pressure_Pa',[Pv;Pv]);
end

function f=fixture(old,new,y,sink,dz,dt,top,bottomLiquid)
    dz=dz(:);n=numel(dz);
    transport=struct('free_dry_air_mass_flux_g_cm2_s',zeros(n,1), ...
        'gas_dispersion_cm2_s',zeros(n,1),'liquid_flux_cm_s',bottomLiquid*ones(n,1));
    if bottomLiquid==0
        bottom=struct('applied_kind','fixed_total_molar_flux','flux_into_domain_mol_cm2_s',0);
    else
        bottom=struct('applied_kind','zero_gas_phase_flux_with_liquid_outflow', ...
            'liquid_flux_into_domain_cm_s',bottomLiquid);
    end
    boundary=struct('species','O2','top',struct('applied_kind','fixed_dry_mole_fraction', ...
        'dry_mole_fraction',top),'bottom',bottom);
    f=struct('args',{{y(:),y(:),old,new,transport,sink(:),dz,dt,boundary}}, ...
        'dz',dz,'dt',dt,'bottomLiquid',bottomLiquid);
end

function [total,Cend,gas,liquid]=inventory(phase,y,dz,p)
% Independent pressure-form Henry calculation, with the native R rounding.
    n=numel(dz);T=phase.temperature_C(:)+273.15;
    dry=phase.gas_pressure_Pa(:)-phase.vapor_pressure_Pa(:);
    KH=p.henry_reference_mol_L_atm.*exp(p.henry_enthalpy_over_R_K ...
        .*(1./T-1/p.henry_reference_temperature_K));
    partialAtm=dry.*y(:)/101325;
    Rrounding=p.gas_constant_cm3_atm_mol_K/(8.31446261815324*1e6/101325);
    Cnode=KH.*partialAtm.*Rrounding; % mol/L water, no soil-volume factor
    Cend=[Cnode(1:n),Cnode(2:n+1)];
    gasNode=dry.*y(:)./(8.31446261815324*T)*1e-6; % mol/cm3 gas
    gasEndpoint=[gasNode(1:n),gasNode(2:n+1)];
    gas=sum(dz(:).*mean(phase.theta_g.*gasEndpoint,2));
    liquid=sum(dz(:).*mean(phase.theta_l.*Cend/1000,2));
    total=gas+liquid;
end

function legacy=prepareLegacy(source,outputRoot)
    code=fileread(source);
    assert(contains(code,'yEndpoint * 1000;'),'verifyOxygenStage0:baselineNotOld', ...
        'The reference must contain the actual pre-repair dissolved-output bug.');
    pattern='(?m)^(\s*)solveDryMoleFractionTransport\(';
    assert(numel(regexp(code,pattern,'match'))==1,'verifyOxygenStage0:baselineSignature');
    code=regexprep(code,pattern,'$1stage0LegacyOxygenKernel(','once');
    folder=fullfile(outputRoot,'reference_kernel');mkdir(folder);
    target=fullfile(folder,'stage0LegacyOxygenKernel.m');
    fid=fopen(target,'w','n','UTF-8');assert(fid>=0);
    close=onCleanup(@()fclose(fid));fprintf(fid,'%s',code);clear close;
    clear('stage0LegacyOxygenKernel');addpath(folder,'-begin');rehash;
    legacy=str2func('stage0LegacyOxygenKernel');
end

function checkClose(actual,expected,key,message)
    scale=max([realmin;abs(actual(:));abs(expected(:))]);
    roundoff=256*eps(scale);
    delta=max(abs(actual(:)-expected(:)));
    assert(isequal(size(actual),size(expected))&&all(isfinite(actual(:))) ...
        &&delta<=roundoff,['verifyOxygenStage0:' key], ...
        '%s actual(1)=%.17g expected(1)=%.17g max_diff=%.17g roundoff=%.17g', ...
        message,actual(1),expected(1),delta,roundoff);
end

function h=fileHash(source)
    fid=fopen(source,'rb');assert(fid>=0);close=onCleanup(@()fclose(fid)); %#ok<NASGU>
    md=java.security.MessageDigest.getInstance('SHA-256');md.update(fread(fid,Inf,'*uint8'));
    h=lower(reshape(dec2hex(typecast(md.digest(),'uint8'),2).',1,[]));
end

function writeReceipt(target,report)
    fid=fopen(target,'w','n','UTF-8');assert(fid>=0);close=onCleanup(@()fclose(fid)); %#ok<NASGU>
    fprintf(fid,'%s\n',jsonencode(report,PrettyPrint=true));
end

function restorePath(oldPath)
    clear('stage0LegacyOxygenKernel');path(oldPath);
end
