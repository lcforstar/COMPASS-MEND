function receipt = compass_batch_entry()
%COMPASS_BATCH_ENTRY Portable batch adapter around the unchanged run_site.
% Writes one launch-specific JSON next to its log, using feature('getpid').
% Physical stepping, diagnostics and conservation remain owned by src.
% Based on the COMPASS entry flow by Jizixing Zhu.
% Editor: Wenjie Liu; Date: 2026-10-31.
    projectRoot=fileparts(fileparts(mfilename('fullpath')));
    siteSelector=getenv('COMPASS_LAUNCH_SITE');
    if isempty(siteSelector),siteSelector='changwu';end
    receiptPath=getenv('COMPASS_RUN_RECEIPT');
    if isempty(receiptPath)
        error('compass:launcher:missingReceipt', ...
            'Use run.bat/run.sh or provide COMPASS_RUN_RECEIPT.');
    end
    receipt=struct('schema_version',1,'run_id',getenv('COMPASS_RUN_ID'), ...
        'project_root',projectRoot,'site_dir',siteSelector, ...
        'matlab_pid',feature('getpid'),'status','starting', ...
        'last_completed_kt',0,'last_kt_source','not_started', ...
        'target_steps',[],'configured_steps',[],'full_window_completed',false, ...
        'output_dir',getenv('COMPASS_OUTPUT'), ...
        'output_source',getenv('COMPASS_LAUNCH_OUTPUT_SOURCE'), ...
        'log_file',getenv('COMPASS_LAUNCH_LOG'),'started_utc',utcNow(), ...
        'updated_utc',utcNow(),'completed_utc','', ...
        'error',struct('identifier','','message','','report',''));
    writeReceipt(receiptPath,receipt);
    model=[];
    try
        cd(projectRoot);
        addpath(projectRoot,'-begin');
        addpath(fullfile(projectRoot,'src'),'-begin');
        assert(strcmp(which('run_site'),fullfile(projectRoot,'run_site.m')), ...
            'compass:launcher:foreignEntry','run_site did not resolve to this copy.');
        if ~isfolder(siteSelector)
            siteSelector=fullfile(projectRoot,'data','sites',siteSelector);
        end
        [ok,info]=fileattrib(siteSelector);
        assert(ok,'compass:launcher:missingSite','Selected site does not exist.');
        receipt.site_dir=info.Name;
        config=io.loadSiteConfig(receipt.site_dir);
        receipt.configured_steps=configuredSteps(config);
        receipt.target_steps=receipt.configured_steps;
        limit=str2double(io.getEnv('MAX_KT'));
        if isfinite(limit)&&limit>0&&limit<receipt.target_steps
            receipt.target_steps=limit;
        end
        receipt.diagnostic_profile=io.getEnv('DIAGNOSTIC_PROFILE');
        if isempty(receipt.diagnostic_profile),receipt.diagnostic_profile='full';end
        receipt.status='running';receipt.updated_utc=utcNow();
        writeReceipt(receiptPath,receipt);
        fprintf('[COMPASS-LAUNCH] run_id=%s MATLAB_PID=%d receipt=%s\n', ...
            receipt.run_id,receipt.matlab_pid,receiptPath);
        model=run_site(receipt.site_dir);
        assert(isfield(model,'runtime')&&isfield(model.runtime,'run_control'), ...
            'compass:launcher:noRunControl','run_site returned no run-control record.');
        control=model.runtime.run_control;
        receipt.target_steps=control.n_steps;
        receipt.last_completed_kt=control.last_completed_kt;
        receipt.last_kt_source='model.runtime.run_control';
        receipt.output_dir=model.runtime.output.output_dir;
        receipt.model_status=char(string(model.status));
        receipt.main_loop_status=char(string(model.main_loop.status));
        assert(strcmp(receipt.model_status,'run_completed') ...
            &&strcmp(receipt.main_loop_status,'completed') ...
            &&receipt.last_completed_kt==receipt.target_steps, ...
            'compass:launcher:incompleteRun', ...
            'run_site stopped before the requested accepted-step endpoint.');
        if isfield(model.runtime,'finalization')&&receipt.target_steps>0
            assert(model.runtime.finalization.flushed, ...
                'compass:launcher:unflushedDiagnostics','Final diagnostics were not flushed.');
        end
        if isfield(model.main_loop,'conservation')
            assert(strcmp(model.main_loop.conservation.status,'PASS'), ...
                'compass:launcher:nativeConservation','Native conservation report did not pass.');
            if isfield(model.main_loop.conservation,'full_budget')
                assert(model.main_loop.conservation.full_budget.passed, ...
                    'compass:launcher:nativeBudget','Native full budget did not pass.');
            end
        end
        receipt.full_window_completed=receipt.last_completed_kt==receipt.configured_steps;
        receipt.status='completed';receipt.updated_utc=utcNow();
        receipt.completed_utc=receipt.updated_utc;
        writeReceipt(receiptPath,receipt);
        fprintf('[COMPASS-LAUNCH] completed KT=%d/%d full_window=%d\n', ...
            receipt.last_completed_kt,receipt.target_steps,receipt.full_window_completed);
    catch exception
        if isstruct(model)&&isfield(model,'runtime')&&isfield(model.runtime,'run_control')
            receipt.last_completed_kt=model.runtime.run_control.last_completed_kt;
            receipt.last_kt_source='model.runtime.run_control';
        else
            receipt.last_completed_kt=lastAcceptedLogKt(receipt.log_file);
            receipt.last_kt_source='this_run_log_accepted_markers';
        end
        receipt.status='failed';receipt.updated_utc=utcNow();
        receipt.completed_utc=receipt.updated_utc;
        receipt.error=struct('identifier',exception.identifier, ...
            'message',exception.message,'report',getReport(exception,'extended','hyperlinks','off'));
        writeReceipt(receiptPath,receipt);
        rethrow(exception);
    end
end

function steps=configuredSteps(config)
    a=datetime(char(config.time.start),'InputFormat','yyyy-MM-dd''T''HH:mm:ss');
    b=datetime(char(config.time.end),'InputFormat','yyyy-MM-dd''T''HH:mm:ss');
    steps=seconds(b-a)/double(config.time.dt);
    assert(isfinite(steps)&&steps>=1&&steps==fix(steps), ...
        'compass:launcher:invalidWindow','The configured window must contain whole steps.');
end

function kt=lastAcceptedLogKt(path)
    kt=0;
    if isempty(path)||~isfile(path),return;end
    text=fileread(path);
    matches=regexp(text,'\[COMPASS\] KT=(\d+)/\d+[^\r\n]*accepted=1(?:\D|$)','tokens');
    if ~isempty(matches),kt=str2double(matches{end}{1});end
end

function value=utcNow()
    value=char(datetime('now','TimeZone','UTC', ...
        'Format','yyyy-MM-dd''T''HH:mm:ss.SSS''Z'''));
end

function writeReceipt(path,value)
    parent=fileparts(path);
    if ~isempty(parent)&&~isfolder(parent),mkdir(parent);end
    temporary=[path '.writing_' num2str(feature('getpid')) '.json'];
    fid=fopen(temporary,'w','n','UTF-8');
    assert(fid>=0,'compass:launcher:receiptOpen','Unable to open run receipt.');
    try
        fprintf(fid,'%s\n',jsonencode(value,'PrettyPrint',true));
        fclose(fid);
    catch exception
        fclose(fid);rethrow(exception);
    end
    [ok,message]=movefile(temporary,path,'f');
    assert(ok,'compass:launcher:receiptWrite','Unable to install run receipt: %s',message);
end
