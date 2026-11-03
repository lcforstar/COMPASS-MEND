function receipt=verifyStandaloneRelease(targetRoot,outputRoot,nSteps)
%VERIFYSTANDALONERELEASE Exercise only the selected standalone release.
% This portable harness uses the release root and src paths only. It calls
% the public run_site launcher with the normal full diagnostic registry and
% an explicit short-run MAX_KT override. No verification archive, calibration
% tool, server helper, old workspace path, altered site or relaxed physical
% tolerance is needed. It is not an annual or full-test-suite validation.
% Editor: Wenjie Liu; Date: 2026-10-31.
    if nargin<3,nSteps=4;end
    targetRoot=char(java.io.File(targetRoot).getCanonicalPath());
    assert(isfolder(targetRoot)&&isfile(fullfile(targetRoot,'run_site.m')));
    if nargin<2||isempty(outputRoot)
        outputRoot=fullfile(targetRoot,'data','output', ...
            ['release_smoke_' char(datetime('now','Format','yyyyMMdd_HHmmss'))]);
    end
    outputRoot=char(java.io.File(outputRoot).getCanonicalPath());
    assert(~isfolder(outputRoot),'standaloneRelease:noOverwrite','Choose a new smoke output directory.');
    allowed={fullfile(targetRoot,'data','output'),fullfile(targetRoot,'output')};
    assert(any(cellfun(@(p)within(outputRoot,p),allowed)), ...
        'standaloneRelease:outputBoundary','Smoke outputs must stay within a release output directory.');
    assert(isscalar(nSteps)&&isfinite(nSteps)&&nSteps>=1&&nSteps==fix(nSteps));

    originalPath=path;originalDirectory=pwd;
    suffixes={'OUTPUT','MAX_KT','DIAGNOSTIC_PROFILE','RESUME_FILE','WORKER_ROOT','THREEWAY_COMPARE'};
    envNames=cell(0,1);envValues=cell(0,1);
    for prefix={'COMPASS_','STEMMUS_'}
        for k=1:numel(suffixes)
            envNames{end+1,1}=[prefix{1} suffixes{k}]; %#ok<AGROW>
            envValues{end+1,1}=getenv(envNames{end}); %#ok<AGROW>
            setenv(envNames{end},'');
        end
    end
    cleanup=onCleanup(@()restoreProcess(originalPath,originalDirectory,envNames,envValues)); %#ok<NASGU>
    % In a fresh -batch process this removes inherited source-package paths.
    % The loaded harness may live anywhere; every production function must
    % resolve within the requested standalone target, checked below.
    restoredefaultpath;
    cd(targetRoot);addpath(targetRoot,fullfile(targetRoot,'src'),'-begin');
    assert(usejava('jvm'),'standaloneRelease:JVM','Normal MEND provenance/checkpoint hashing requires MATLAB JVM.');
    entries={'run_site','runCOMPASS','init.initialize','io.checkSite', ...
        'diagnostics.acceptedOrganicAudit','diagnostics.writeStep', ...
        'soil.microbe.advanceSourceCandidate','soil.nitrogen.tracer.advanceCandidate'};
    resolved=cell(numel(entries),1);
    for k=1:numel(entries)
        location=which(entries{k});
        assert(~isempty(location)&&within(location,targetRoot), ...
            'standaloneRelease:sourceShadow','%s resolves outside the release: %s',entries{k},location);
        resolved{k}=struct('entry',entries{k},'file',location,'sha256',fileHash(location));
    end
    assert(strcmp(which('run_site'),fullfile(targetRoot,'run_site.m')));
    setenv('COMPASS_OUTPUT',outputRoot);
    setenv('COMPASS_MAX_KT',num2str(nSteps));
    setenv('COMPASS_DIAGNOSTIC_PROFILE','full');
    setenv('COMPASS_THREEWAY_COMPARE','0');

    wall=tic;
    model=run_site('changwu');
    assert(model.check_report.ok);
    assert(strcmp(model.status,'run_completed'));
    assert(strcmp(model.main_loop.status,'completed'));
    assert(model.main_loop.last_completed_kt==nSteps);
    assert(model.runtime.run_control.last_completed_kt==nSteps);
    assert(~model.runtime.resume.requested);
    assert(strcmp(model.runtime.diagnostics.profile,'full'));
    assert(strcmp(model.built.parameters.soil_microbes.reaction_method,'source_hourly'));
    assert(strcmp(model.built.parameters.soil_dom.sorption_model,'mend_dynamic'));
    reports=model.main_loop.step_reports;
    assert(all([reports.accepted])&&~any([reports.step_failed]));
    assert(isequal([reports.kt],1:nSteps));
    for k=1:nSteps
        assert(reports(k).clock_time==model.built.control.time.start_clock ...
            +seconds(k*model.built.control.time.Delt_t));
    end
    conservation=model.main_loop.conservation;
    assert(strcmp(conservation.status,'PASS')&&conservation.last_kt==nSteps);
    assert(conservation.full_budget.passed&&conservation.full_budget.coverage_complete);
    assert(conservation.full_budget.carbon_tolerance_g_cm2==1e-10);
    assert(conservation.full_budget.nitrogen_tolerance_g_cm2==1e-12);
    assert(all(conservation.max_hourly_CN_residual_g_cm2<=[1e-10,1e-12]));
    assert(all(conservation.max_interval_DIN_residual_ug_cm2<1e-6));
    assert(all(conservation.max_cumulative_DIN_residual_ug_cm2<1e-6));
    assert(conservation.max_O2_residual_mol_cm2<1e-14);
    assert(abs(conservation.cumulative_O2_residual_mol_cm2)<1e-14);

    fields={'POM1_decomp','POM2_decomp','POM1_to_DOM','POM2_to_DOM', ...
        'POM1_to_MOM','POM2_to_MOM','MOM_to_DOM','DOM_to_MBA', ...
        'CO2_growth','CO2_maintenance_MBA','CO2_maintenance_MBD', ...
        'MBA_mortality','MBA_to_DOM','MBA_to_POM1','MBA_to_POM2', ...
        'MBA_to_ENZP1','MBA_to_ENZP2','MBA_to_ENZM', ...
        'ENZP1_to_DOM','ENZP2_to_DOM','ENZM_to_DOM', ...
        'MBA_to_MBD','MBD_to_MBA','DOM_to_QOM','QOM_to_DOM', ...
        'DOM_sorption_net','ENZN_production','ENZN_to_MBA','N_overflow_CO2'};
    registry=model.runtime.diagnostics.registry;
    fluxChecks=cell(numel(fields),1);
    for k=1:numel(fields)
        profile=['mend_' lower(fields{k}) '_profile'];
        rate=['mend_' lower(fields{k}) '_flux'];
        buffer=registry.get(profile);column=registry.get(rate);
        assert(all(isfinite(buffer(1:nSteps,:)),'all')&&all(isfinite(column(1:nSteps,:)),'all'));
        meta=registry.getMeta(profile);
        file=fullfile(outputRoot,'nc',[meta.group '.nc']);
        written=ncread(file,profile);
        assert(~isempty(written)&&all(isfinite(written),'all'));
        rateMeta=registry.getMeta(rate);
        rateFile=fullfile(outputRoot,'nc',[rateMeta.group '.nc']);
        writtenRate=ncread(rateFile,rate);
        assert(~isempty(writtenRate)&&all(isfinite(writtenRate),'all'));
        fluxChecks{k}=struct('field',fields{k},'registry_profile',profile, ...
            'registry_column_rate',rate,'native_nc_file',file, ...
            'native_rate_nc_file',rateFile,'finite_written',true);
    end
    for key={'soc_observation_profile','doc_carbon_profile','qoc_carbon_profile'}
        buffer=registry.get(key{1});assert(all(isfinite(buffer(1:nSteps,:)),'all'));
    end
    sourcePrefixes={'fertilizer_tracer_','manure_organic_tracer_','manure_tan_tracer_'};
    sourceState={model.built.Y.soil_fertilizer_n_tracer, ...
        model.built.Y.soil_manure_n_tracer.organic,model.built.Y.soil_manure_n_tracer.tan};
    sourceChecks=cell(3,1);
    for k=1:3
        name=[sourcePrefixes{k} 'cumulative_closure_residual'];
        residual=registry.get(name);residual=residual(1:nSteps,:);
        cumulative=sourceState{k}.cumulative;
        tolerance=2e-12*max(1,cumulative.fertilizer_input_g_cm2);
        assert(all(isfinite(residual),'all')&&all(abs(residual)<=tolerance,'all'));
        assert(cumulative.numerical_correction_g_cm2==0);
        meta=registry.getMeta(name);
        file=fullfile(outputRoot,'nc',[meta.group '.nc']);
        written=ncread(file,name);
        assert(~isempty(written)&&all(isfinite(written),'all') ...
            &&all(abs(written)<=tolerance,'all'));
        sourceChecks{k}=struct('prefix',sourcePrefixes{k}, ...
            'max_abs_column_closure_g_cm2',max(abs(residual),[],'all'), ...
            'original_tolerance_g_cm2',tolerance,'numerical_correction_g_cm2',0, ...
            'registry_residual',name,'native_nc_file',file,'finite_written',true);
    end
    f=model.built.forcing;n=numel(f.timeAxis.clockTime);
    assert(f.timeAxis.clockTime(1)==model.built.control.time.start_clock);
    assert(numel(f.solar.tts)==n&&f.solar.n_hours==n&&f.hourly.n_hours==n);

    receipt=struct('status','PASS','target_root',targetRoot,'output_root',outputRoot, ...
        'matlab_version',version,'matlab_release',version('-release'),'jvm_enabled',usejava('jvm'), ...
        'resolved_production_sources',{vertcat(resolved{:})}, ...
        'requested_steps',nSteps,'accepted_steps',nSteps,'elapsed_s',toc(wall), ...
        'public_entrypoint','run_site(''changwu'') with MAX_KT override', ...
        'normal_full_diagnostics',true,'mend_flux_contract_count',numel(fields), ...
        'mend_flux_checks',{vertcat(fluxChecks{:})}, ...
        'three_source_checks',{vertcat(sourceChecks{:})}, ...
        'native_conservation',conservation,'full_suite_run',false, ...
        'annual_run_performed',false,'site_or_source_modified',false, ...
        'coverage_limit','First June hours; this is an independent release smoke, not new September-event or annual validation');
    receiptFile=fullfile(outputRoot,'standalone_release_smoke_receipt.json');
    fid=fopen(receiptFile,'w','n','UTF-8');assert(fid>=0);close=onCleanup(@()fclose(fid)); %#ok<NASGU>
    fprintf(fid,'%s\n',jsonencode(receipt,PrettyPrint=true));
    disp(jsonencode(receipt));
end

function yes=within(candidate,root)
    candidate=strrep(char(candidate),'\','/');root=strrep(char(root),'\','/');
    if ispc,candidate=lower(candidate);root=lower(root);end
    yes=strcmp(candidate,root)||startsWith(candidate,[root '/']);
end

function h=fileHash(path)
    fid=fopen(path,'rb');assert(fid>=0);close=onCleanup(@()fclose(fid)); %#ok<NASGU>
    digest=java.security.MessageDigest.getInstance('SHA-256');
    while ~feof(fid),bytes=fread(fid,1048576,'*uint8');if ~isempty(bytes),digest.update(bytes);end,end
    h=lower(reshape(dec2hex(typecast(digest.digest(),'uint8'),2).',1,[]));
end

function restoreProcess(originalPath,originalDirectory,names,values)
    path(originalPath);cd(originalDirectory);
    for k=1:numel(names),setenv(names{k},values{k});end
end
