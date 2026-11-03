function receipt=verifyDinTransportNumerics(targetRoot,outputRoot)
%VERIFYDINTRANSPORTNUMERICS Public-entry 00/10/01/11/missing switch regression.
% Exercises two accepted hours with original physical gates, checkpoint identity,
% true/false transport outcomes and the signed small-Peclet analytic limit.
% Each site is an isolated copy of the target-root Changwu example. A fresh
% output directory is required. No edits to user input or acceptance thresholds.
% Editor: Wenjie Liu; Date: 2026-11-03.
    if nargin<1,targetRoot=fileparts(fileparts(mfilename('fullpath')));end
    targetRoot=char(java.io.File(targetRoot).getCanonicalPath());
    if nargin<2,outputRoot=tempname;end
    assert(~isfolder(outputRoot),'testDinFlags:exists','Choose a fresh test directory.');
    origPath=path;oldDir=pwd;names={'COMPASS_OUTPUT','COMPASS_MAX_KT','COMPASS_RESUME_FILE','COMPASS_DIAGNOSTIC_PROFILE'};
    old=cellfun(@getenv,names,'UniformOutput',false);restore=onCleanup(@()resetTest(origPath,oldDir,names,old));
    cd(targetRoot);addpath(targetRoot,'-begin');addpath(fullfile(targetRoot,'src'),'-begin');
    assert(strcmp(which('run_site'),fullfile(targetRoot,'run_site.m')));
    mkdir(outputRoot);mkdir(fullfile(outputRoot,'sites'));
    copyfile(fullfile(targetRoot,'data','model_library_nc'),fullfile(outputRoot,'model_library_nc'));
    setenv('COMPASS_MAX_KT','2');setenv('COMPASS_RESUME_FILE','');setenv('COMPASS_DIAGNOSTIC_PROFILE','driver');
    cfgs={'00','10','01','11','missing'};results=struct;
    for i=1:numel(cfgs)
        id=cfgs{i};site=fullfile(outputRoot,'sites',['flags_' id]);
        copyfile(fullfile(targetRoot,'data','sites','changwu'),site);
        p=fullfile(site,'site.toml');buffer=fileread(p);
        for key={'upwind','artificial_disp'}
            flag=key{1};pattern=['(?m)^(' flag '[\t ]*=[\t ]*)(true|false)([^\n]*)$'];
            assert(numel(regexp(buffer,pattern,'start'))==1);
            if strcmp(id,'missing'),replacement='';
            else
                if strcmp(flag,'upwind'),bit=id(1);else,bit=id(2);end
                if bit=='1',value='true';else,value='false';end
                replacement=['$1' value '$3'];
            end
            buffer=regexprep(buffer,pattern,replacement);
        end
        h=fopen(p,'w','n','UTF-8');fprintf(h,'%s',buffer);fclose(h);
        setenv('COMPASS_OUTPUT',fullfile(outputRoot,['output_' id]));
        model=beginRun(run_site(site,struct('stop_after_initialization',true)));
        useUp=double(ismember(id,{'10','11'}));useArt=double(ismember(id,{'01','11'}));
        flags=soil.nitrogen.transport.numericalOptions(model.built.control);
        assert(isequal(flags,struct('lUpW',useUp,'lArtD',useArt)));
        sig=io.checkpointCompatibilitySignature(model.built);
        if useUp||useArt
            assert(sig.din_transport_numerics.lUpW==useUp&&sig.din_transport_numerics.lArtD==useArt);
            assert(isfield(sig.din_transport_numerics,'upwind_weight')==logical(useUp));
        else,assert(~isfield(sig,'din_transport_numerics'));end
        initial=model.built.Y;prev=[0 0];rows=zeros(2,10);
        for n=1:2
            oldY=model.built.Y;[model,r]=runOneStep(model,n);assert(r.accepted&&~r.step_failed);
            [rows(n,:),~]=diagnostics.acceptedOrganicAudit(model.built,oldY,initial,prev,n,true,n==1);prev=rows(n,6:7);
        end
        assert(abs(sum(rows(:,8)))<1e-14);
        results.(['mode' id])=struct('Dsh',model.built.Y.soil_nitrogen.variables.Dsh, ...
            'NH4',model.built.Y.nh4,'NO3',model.built.Y.no3,'rows',rows);
    end
    assert(isequal(results.mode00,results.modemissing));
    assert(~isequal(results.mode10.NH4,results.mode00.NH4) ...
        &&~isequal(results.mode10.Dsh,results.mode00.Dsh)&&~isequal(results.mode01.Dsh,results.mode00.Dsh));
    for imbalance=[-1e-12,0,1e-12]
        q=[1e-3;-1e-3*(1-imbalance)];x=mean(q)/.02;
        [~,weight]=soil.nitrogen.transport.pecletCourantNumber(1,[.01 .01],1,[.25 .25],q,1,3600,0,0,0,0,[1 1]);
        assert(abs(weight-x/3)<1e-20);
    end
    receipt=struct('status','PASS','site_modes',{cfgs},'accepted_hours_per_mode',2, ...
        'missing_zero_bitwise',true,'eligible_flags_change_ADE',true,'signed_small_Pe',true);
    h=fopen(fullfile(outputRoot,'din_flags_verified.json'),'w');fprintf(h,'%s',jsonencode(receipt));fclose(h);
end
function resetTest(route,folder,keys,values)
    path(route);cd(folder);for i=1:numel(keys),setenv(keys{i},values{i});end
end
