function receipt=verify_nh4_protection_windows(codeRoot,outputRoot,mode,variant,nSteps,baselineFile,resumeFile)
%VERIFY_NH4_PROTECTION_WINDOWS Capture one isolated coupled regression window.
% Run each CODE ROOT / MODE / VARIANT in a fresh MATLAB process. MODE is
% 'mend' or 'empirical'; VARIANT is 'default', 'zero', or 'protected' (13.95).
% Only an output-owned copy of the selected code root's Changwu site changes.
% BASELINEFILE optionally names another window_capture.mat: all numeric Y
% leaves and every baseline diagnostic series must be isequaln. Added numeric
% fields are reported separately, never substituted for baseline quantities.
% RESUMEFILE optionally resumes an accepted checkpoint; NSTEPS is the target
% absolute KT, not the number of additional steps. This bounded check does
% not certify the annual response or the scientific validity of the threshold.
% Editor: Wenjie Liu; Date: 2026-11-01.
    if nargin<3||isempty(mode),mode='mend';end
    if nargin<4||isempty(variant),variant='default';end
    if nargin<5||isempty(nSteps),nSteps=24;end
    if nargin<6,baselineFile='';end
    if nargin<7,resumeFile='';end
    codeRoot=char(java.io.File(codeRoot).getCanonicalPath());
    outputRoot=char(java.io.File(outputRoot).getCanonicalPath());
    mode=char(string(mode));variant=char(string(variant));
    assert(ismember(string(mode),["mend","empirical"]));
    assert(ismember(string(variant),["default","zero","protected"]));
    assert(isscalar(nSteps)&&isfinite(nSteps)&&nSteps>=1&&nSteps==fix(nSteps));
    assert(isfile(fullfile(codeRoot,'run_site.m')));
    assert(~isfolder(outputRoot),'verifyNh4Window:noOverwrite','Choose a fresh output folder.');
    oldPath=path;oldDirectory=pwd;
    suffixes={'OUTPUT','MAX_KT','DIAGNOSTIC_PROFILE','RESUME_FILE','WORKER_ROOT','THREEWAY_COMPARE'};
    names={};values={};
    for prefix={'COMPASS_','STEMMUS_'}
        for k=1:numel(suffixes)
            names{end+1}=[prefix{1} suffixes{k}]; %#ok<AGROW>
            values{end+1}=getenv(names{end}); %#ok<AGROW>
            setenv(names{end},'');
        end
    end
    cleanup=onCleanup(@()restoreProcess(oldPath,oldDirectory,names,values)); %#ok<NASGU>
    restoredefaultpath;cd(codeRoot);addpath(codeRoot,fullfile(codeRoot,'src'),'-begin');
    assert(startsWith(which('init.initialize'),[fullfile(codeRoot,'src') filesep]), ...
        'verifyNh4Window:sourceShadow','Production initialization resolves outside requested root.');
    mkdir(outputRoot);mkdir(fullfile(outputRoot,'data','sites'));
    siteDir=fullfile(outputRoot,'data','sites','changwu');
    copyfile(fullfile(codeRoot,'data','sites','changwu'),siteDir);
    % buildScope resolves its shared library two levels above the site.
    copyfile(fullfile(codeRoot,'data','model_library_nc'), ...
        fullfile(outputRoot,'data','model_library_nc'));
    configFile=fullfile(siteDir,'site.toml');config=fileread(configFile);
    if strcmp(mode,'empirical')
        config=regexprep(config,'(?m)^(\s*soil_carbon_model\s*=\s*)true','$1false');
        config=regexprep(config,'(?m)^(\s*manure_source_tracer\s*=\s*)true','$1false');
        config=regexprep(config,'(?m)^(\s*fertilizer_source_tracer\s*=\s*)true','$1false');
        config=regexprep(config,'(?m)^(\s*mend_manure_fast_input_enabled\s*=\s*)true','$1false');
    end
    assert(isempty(regexp(config,'(?m)^\s*\[nitrogen\]\s*$','once')), ...
        'verifyNh4Window:unexpectedNitrogenSection','Fixture source already has a nitrogen section.');
    if ~strcmp(variant,'default')
        threshold=0;if strcmp(variant,'protected'),threshold=13.95;end
        config=sprintf('%s\n[nitrogen]\nnh4_protected_per_clay_mg_kg = %.17g\n',config,threshold);
    end
    fid=fopen(configFile,'w','n','UTF-8');assert(fid>=0);fprintf(fid,'%s',config);fclose(fid);
    setenv('COMPASS_OUTPUT',fullfile(outputRoot,'run'));
    setenv('COMPASS_MAX_KT',num2str(nSteps));
    setenv('COMPASS_DIAGNOSTIC_PROFILE','full');setenv('COMPASS_THREEWAY_COMPARE','0');
    if ~isempty(resumeFile),setenv('COMPASS_RESUME_FILE',char(resumeFile));end
    receipt=struct('status','RUNNING','code_root',codeRoot,'mode',mode,'variant',variant, ...
        'target_kt',nSteps,'production_site_modified',false,'baseline_file',char(baselineFile), ...
        'resume_file',char(resumeFile));
    writeJson(fullfile(outputRoot,'window_receipt.json'),receipt);
    wall=tic;
    try
        model=run_site(siteDir,struct('stop_after_initialization',true));
        model=beginRun(model);startKt=model.runtime.run_control.next_kt;
        assert(startKt<=nSteps,'verifyNh4Window:emptyWindow','Requested target precedes the next accepted step.');
        reports=repmat(solver.emptyStepReport(),nSteps,1);
        initialY=model.built.Y;previousDIN=[0 0];auditRows=nan(nSteps,10);
        oxygen=nan(nSteps,1);auditOn=strcmp(mode,'mend');cumulative=struct();
        if startKt>1
            f=initialY.soil_nitrogen.flux;
            previousNH4=(f.cvnh4(end)-f.cv_nh4_init_mg_cm2) ...
                -(f.cum_nh4_in_mg_cm2-f.cum_nh4_out_mg_cm2);
            previousNO3=(f.cvno3(end)-f.cv_no3_init_mg_cm2) ...
                -(f.cum_no3_in_mg_cm2-f.cum_no3_out_mg_cm2);
            previousDIN=[previousNH4 previousNO3];
        end
        for kt=startKt:nSteps
            previousY=model.built.Y;
            [model,reports(kt)]=runOneStep(model,kt);
            assert(reports(kt).accepted&&~reports(kt).step_failed, ...
                'verifyNh4Window:rejected','Coupled window rejected KT=%d.',kt);
            nd=model.built.workspace.nitrogen;
            oxygen(kt)=nd.oxygen.balance_residual_mol_cm2;
            assert(isfinite(oxygen(kt))&&abs(oxygen(kt))<1e-14, ...
                'verifyNh4Window:oxygen','Accepted O2 residual failed at KT=%d.',kt);
            assert(abs(sum(oxygen(startKt:kt)))<1e-14, ...
                'verifyNh4Window:oxygenCumulative','Accepted cumulative O2 residual failed.');
            assert(all(abs([nd.cum_residual_nh4,nd.cum_residual_no3])<1e-6), ...
                'verifyNh4Window:DIN','Native cumulative DIN gate failed.');
            if auditOn
                [auditRows(kt,:),cumulative]=diagnostics.acceptedOrganicAudit( ...
                    model.built,previousY,initialY,previousDIN,kt,startKt==1,kt==startKt);
                previousDIN=auditRows(kt,6:7);
            end
            fprintf('NH4 window %s/%s KT=%d/%d accepted\n',mode,variant,kt,nSteps);
        end
        model.main_loop.step_reports=reports;model.main_loop.status='completed';
        model.main_loop.last_completed_kt=nSteps;model.main_loop.next_kt=nSteps+1;
        assert(model.runtime.run_control.last_completed_kt==nSteps ...
            &&strcmp(model.main_loop.status,'completed'),'verifyNh4Window:incomplete','Window did not complete.');
        reports=model.main_loop.step_reports;
        completed=[reports.kt]>0;
        assert(all([reports(completed).accepted])&&~any([reports(completed).step_failed]));
        assert(logical(model.built.control.switches.lSoilOrganic)==strcmp(mode,'mend'));
        capture=struct('schema_version',1,'mode',mode,'variant',variant, ...
            'target_kt',nSteps,'first_kt',startKt,'Y_numeric',numericLeaves(model.built.Y,'Y'), ...
            'nitrogen_numeric',numericLeaves(model.built.workspace.nitrogen,'nitrogen'), ...
            'diagnostics',struct(),'diagnostic_time_dimensions',struct());
        capture.oxygen_hourly_mol_cm2=oxygen;
        if auditOn,capture.organic_audit_rows=auditRows;end
        receipt.physical_window_completed=true;
        receipt.checkpoint_file=io.saveCompassCheckpoint(model,nSteps);
        % Preserve the accepted physical result before ancillary export work.
        acceptedY=model.built.Y; %#ok<NASGU>
        save(fullfile(outputRoot,'accepted_window_state.mat'),'acceptedY','capture','-v7.3');
        writeJson(fullfile(outputRoot,'window_receipt.json'),receipt);
        registry=model.runtime.diagnostics.registry;keys=sort(registry.list());
        for k=1:numel(keys)
            receipt.extracting_diagnostic=keys{k};
            value=registry.get(keys{k});
            if isnumeric(value)||islogical(value)
                meta=registry.getMeta(keys{k});timeDimension=find(strcmp(meta.dims,'time'));
                capture.diagnostic_time_dimensions.(keys{k})=timeDimension;
                if isempty(timeDimension)
                    capture.diagnostics.(keys{k})=value;
                else
                    assert(isscalar(timeDimension)&&size(value,timeDimension)>=nSteps, ...
                        'verifyNh4Window:diagnosticShape', ...
                        'Diagnostic %s has declared time dimension %s but size %s at KT=%d.', ...
                        keys{k},mat2str(timeDimension),mat2str(size(value)),nSteps);
                    indices=repmat({':'},1,max([2,ndims(value),timeDimension]));
                    indices{timeDimension}=1:nSteps;
                    capture.diagnostics.(keys{k})=value(indices{:});
                end
            end
        end
        % Native accepted-state O2 gates also apply to the empirical path,
        % which does not enter runCOMPASS's organic-budget audit branch.
        capture.oxygen_hourly_mol_cm2=oxygen;
        receipt.conservation=struct('status','PASS','first_kt',startKt,'last_kt',nSteps, ...
            'max_O2_residual_mol_cm2',max(abs(oxygen(startKt:nSteps))), ...
            'cumulative_O2_residual_mol_cm2',sum(oxygen(startKt:nSteps)));
        if auditOn
            capture.organic_audit_rows=auditRows;
            receipt.conservation.full_budget=cumulative;
            receipt.conservation.max_hourly_CN_residual_g_cm2=max(abs(auditRows(startKt:nSteps,2:3)),[],1);
            receipt.conservation.max_cumulative_DIN_residual_ug_cm2=max(abs(auditRows(startKt:nSteps,6:7)),[],1);
        end
        captureFile=fullfile(outputRoot,'window_capture.mat');
        save(captureFile,'capture','-v7.3');
        receipt.capture_file=captureFile;
        if ~isempty(baselineFile)
            reference=load(baselineFile,'capture');
            assert(strcmp(reference.capture.mode,mode)&&reference.capture.target_kt==nSteps);
            receipt.bitwise=compareCapture(reference.capture,capture);
        end
        model=finalizeRun(model);model.status='run_completed';
        receipt=rmfield(receipt,'extracting_diagnostic');
        receipt.elapsed_seconds=toc(wall);receipt.status='PASS';
        writeJson(fullfile(outputRoot,'window_receipt.json'),receipt);
        fprintf('NH4 window %s/%s PASS; KT=%d; wall=%.3fs\n%s\n',mode,variant,nSteps,receipt.elapsed_seconds,captureFile);
    catch err
        receipt.status='FAIL';receipt.elapsed_seconds=toc(wall);
        receipt.failure=struct('identifier',err.identifier,'message',err.message,'stack',err.stack);
        writeJson(fullfile(outputRoot,'window_receipt.json'),receipt);rethrow(err);
    end
end

function leaves=numericLeaves(value,prefix)
    leaves=struct('path',{},'value',{});
    if isnumeric(value)||islogical(value)
        leaves=struct('path',prefix,'value',value);
    elseif isstruct(value)
        fields=sort(fieldnames(value));
        for j=1:numel(value)
            for k=1:numel(fields)
                label=[prefix '.' fields{k}];
                if numel(value)>1,label=sprintf('%s(%d).%s',prefix,j,fields{k});end
                leaves=[leaves;numericLeaves(value(j).(fields{k}),label)]; %#ok<AGROW>
            end
        end
    elseif iscell(value)
        for j=1:numel(value)
            leaves=[leaves;numericLeaves(value{j},sprintf('%s{%d}',prefix,j))]; %#ok<AGROW>
        end
    end
    leaves=leaves(:);
end

function result=compareCapture(reference,current)
    result=struct('passed',true,'compared_numeric_leaves',0,'added_numeric_leaves',{{}});
    first=1;
    if isfield(current,'first_kt'),first=current.first_kt;end
    if isfield(reference,'first_kt'),first=max(first,reference.first_kt);end
    result.compared_first_kt=first;result.compared_last_kt=current.target_kt;
    for group={'Y_numeric','nitrogen_numeric'}
        name=group{1};a=reference.(name);b=current.(name);paths={b.path};
        for k=1:numel(a)
            where=find(strcmp(paths,a(k).path));
            assert(isscalar(where)&&isequaln(a(k).value,b(where).value), ...
                'verifyNh4Window:bitwise','Numeric baseline mismatch: %s.',a(k).path);
        end
        result.compared_numeric_leaves=result.compared_numeric_leaves+numel(a);
        result.added_numeric_leaves=[result.added_numeric_leaves,setdiff(paths,{a.path})];
    end
    names=fieldnames(reference.diagnostics);
    for k=1:numel(names)
        name=names{k};dimension=[];
        assert(isfield(current.diagnostics,name),'verifyNh4Window:bitwise', ...
            'Baseline diagnostic is missing: %s.',name);
        if isfield(current,'diagnostic_time_dimensions')
            dimension=current.diagnostic_time_dimensions.(name);
        elseif isfield(reference,'diagnostic_time_dimensions')
            dimension=reference.diagnostic_time_dimensions.(name);
        end
        a=reference.diagnostics.(name);b=current.diagnostics.(name);
        if ~isempty(dimension)&&first>1
            index=repmat({':'},1,max([2,ndims(a),dimension]));
            index{dimension}=first:current.target_kt;a=a(index{:});b=b(index{:});
        end
        assert(isfield(current.diagnostics,name)&&isequaln(a,b), ...
            'verifyNh4Window:bitwise','Diagnostic baseline mismatch: %s.',name);
    end
    result.compared_diagnostic_series=numel(names);
    result.added_diagnostic_series=setdiff(fieldnames(current.diagnostics),names);
    for key={'oxygen_hourly_mol_cm2','organic_audit_rows'}
        name=key{1};
        if isfield(reference,name)
            a=reference.(name);b=current.(name);
            if first>1,a=a(first:end,:);b=b(first:end,:);end
            assert(isfield(current,name)&&isequaln(a,b), ...
                'verifyNh4Window:bitwise','Accepted audit series mismatch: %s.',name);
        end
    end
end

function writeJson(file,value)
    fid=fopen(file,'w','n','UTF-8');assert(fid>=0);cleanup=onCleanup(@()fclose(fid)); %#ok<NASGU>
    fprintf(fid,'%s\n',jsonencode(value,PrettyPrint=true));
end

function restoreProcess(oldPath,oldDirectory,names,values)
    path(oldPath);cd(oldDirectory);
    for k=1:numel(names),setenv(names{k},values{k});end
end
