function result=runCOMPASSBiogeochemistrySpinup(siteDir,options)
%RUNCOMPASSBIOGEOCHEMISTRYSPINUP Spin up native AOM/DOM/MEMC/MEND C/N state.
%
%   RESULT = runCOMPASSBiogeochemistrySpinup(SITEDIR,OPTIONS) repeats a fixed
%   forcing window while carrying the complete accepted Y state between cycles.
%   Blank DOC/DON profiles are legal only here. Each cycle uses the ordinary
%   3600 s runOneStep transaction and writes a schema-2 full-state checkpoint.
%
%   OPTIONS: cycles (default 1), n_steps (default configured window), and
%   output_dir. progress_interval controls accepted-step telemetry (default
%   24 h). A checkpoint continues this cycle at its original accepted time;
%   it is NOT a measured-constrained production t0 initialization artifact.
%
%   Editor: Jizixing Zhu; Date: 2026-08-21.
%   Change: 2026-08-24 Jizixing Zhu -- add cycle-level DOM, MEMC, PlantCN,
%           litter, and mass-closure acceptance gates and finalize diagnostics.
%   Change: 2026-08-29 Codex -- document the native reaction path.
%   Change: 2026-09-02 Codex -- add pool-wise periodic convergence, optional
%   bounded slow-pool acceleration, automatic stop, and mandatory unaccelerated
%   re-equilibration. Production reaction equations remain unchanged.
%   Change: 2026-09-09 Jizixing Zhu -- optional cycle-only diagnostic observer
%   and extra convergence gate; no work added to the accepted-step loop.
%   Change: 2026-09-11 Jizixing Zhu -- resume an interrupted first cycle via
%   beginRun's accepted checkpoint contract; subsequent cycles start at KT=1.
%   resume_first_cycle_checkpoint [char] must belong to this first cycle and
%   the same site/forcing. Use a new output directory; old diagnostics are
%   not in the checkpoint and are not reconstructed before the resume KT.
%   Change: 2026-09-29 Wenjie Liu -- bind replay cycles to full immutable
%   identity, rebuild accepted-state caches/root-water memory, preserve canopy
%   seam history, isolate resume environment, and derive the default window.
    arguments
        siteDir (1,:) char
        options.cycles (1,1) double {mustBeInteger,mustBePositive}=1
        options.n_steps (1,1) double=NaN
        options.output_dir (1,:) char=''
        options.progress_interval (1,1) double {mustBeInteger,mustBePositive}=24
        options.auto_stop (1,1) logical=true
        options.relative_tolerance (1,1) double {mustBePositive}=1e-3
        options.absolute_floor (1,1) double {mustBePositive}=1e-12
        options.minimum_cycles (1,1) double {mustBeInteger,mustBePositive}=3
        options.consecutive_cycles (1,1) double {mustBeInteger,mustBePositive}=2
        options.slow_pool_acceleration (1,1) logical=false
        options.acceleration_start_cycle (1,1) double {mustBeInteger,mustBePositive}=3
        options.acceleration_max_factor (1,1) double {mustBeGreaterThanOrEqual(options.acceleration_max_factor,1)}=5
        options.final_rebalance_cycles (1,1) double {mustBeInteger,mustBeNonnegative}=2
        options.cycle_observer=[]
        options.resume_first_cycle_checkpoint (1,:) char=''
    end
    root=fileparts(fileparts(mfilename('fullpath')));
    if isempty(options.output_dir)
        [~,siteName]=fileparts(siteDir);
        outputRoot=fullfile(root,'data','output',[siteName,'_bgc_spinup']);
    else
        outputRoot=options.output_dir;
    end
    if ~isfolder(outputRoot), mkdir(outputRoot); end
    envNames={'COMPASS_OUTPUT','STEMMUS_OUTPUT','COMPASS_MAX_KT', ...
        'STEMMUS_MAX_KT','COMPASS_RESUME_FILE','STEMMUS_RESUME_FILE'};
    envValues=cellfun(@getenv,envNames,'UniformOutput',false);
    cleanup=onCleanup(@()restoreEnvironment(envNames,envValues)); %#ok<NASGU>
    for i=3:numel(envNames),setenv(envNames{i},'');end
    nSteps=options.n_steps; identity=[]; configuredSteps=NaN;
    carry=[]; checkpointFiles=strings(options.cycles,1);
    cycleAudits=cell(options.cycles,1);
    convergence=cell(options.cycles,1); acceleration=cell(options.cycles,1);
    history=cell(0,1); stableCount=0; lastAccelerationCycle=0;
    completedCycles=0; stopReason='maximum_cycles';
    for cycle=1:options.cycles
        cycleDir=fullfile(outputRoot,sprintf('cycle%02d',cycle));
        setenv('COMPASS_OUTPUT',cycleDir);
        model=init.initialize(siteDir);
        if cycle==1
            configuredSteps=model.built.control.time.tS;
            if isnan(nSteps),nSteps=configuredSteps;end
            assert(isfinite(nSteps)&&nSteps>=1&&nSteps==fix(nSteps) ...
                &&nSteps<=configuredSteps, ...
                'runCOMPASSBiogeochemistrySpinup:window', ...
                'n_steps must be an integer inside the configured forcing window.');
            setenv('COMPASS_MAX_KT',sprintf('%d',nSteps));
            identity=init.buildBiogeochemistryCycleIdentity(model);
        else
            assert(isequaln(identity,init.buildBiogeochemistryCycleIdentity(model)), ...
                'runCOMPASSBiogeochemistrySpinup:identity', ...
                'Model/source/forcing identity changed between cycles.');
        end
        fprintf('[BGC cycle-start] cycle=%d/%d steps=%d\n', ...
            cycle,options.cycles,nSteps);
        assert(model.built.control.time.Delt_t==3600, ...
            'runCOMPASSBiogeochemistrySpinup:outerStep', ...
            'Biogeochemistry spin-up requires a 3600 s outer step.');
        assert(model.built.control.switches.lSoilDOM==1, ...
            'runCOMPASSBiogeochemistrySpinup:domDisabled', ...
            'dissolved_organic_matter must be true.');
        model=init.prepareRuntime(model);
        model.runtime.biogeochemistry_spinup=true;
        model.runtime.biogeochemistry_spinup_identity=identity;
        model.runtime.biogeochemistry_spinup_cycle_steps=nSteps;
        if ~isempty(carry)
            model=init.applyBiogeochemistryCycleCarry(model,carry);
        end
        cycleInitial=model.built.Y;
        if ~isempty(options.resume_first_cycle_checkpoint)
            model.runtime.resume.requested=(cycle==1);
            model.runtime.resume.file=options.resume_first_cycle_checkpoint;
            if cycle==1
                saved=io.loadCompassCheckpoint(options.resume_first_cycle_checkpoint);
                assert(isfield(saved.runtime,'biogeochemistry_spinup_identity') ...
                    &&isequaln(saved.runtime.biogeochemistry_spinup_identity,identity) ...
                    &&isfield(saved.runtime,'biogeochemistry_spinup_cycle_steps') ...
                    &&saved.runtime.biogeochemistry_spinup_cycle_steps==nSteps, ...
                    'runCOMPASSBiogeochemistrySpinup:resumeIdentity', ...
                    'Resume requires an identity-bound checkpoint from this spin-up model.');
                if isfield(saved.runtime,'biogeochemistry_cycle')
                    error('runCOMPASSBiogeochemistrySpinup:resumeCycleIndex', ...
                        'resume_first_cycle_checkpoint cannot be a later replay cycle.');
                end
            end
        end
        model=beginRun(model);
        firstKt=model.runtime.run_control.next_kt;
        assert(firstKt>=1&&firstKt<=nSteps+1, ...
            'runCOMPASSBiogeochemistrySpinup:resumeRange','Checkpoint outside forcing window.');
        if cycle==1&&~isempty(options.resume_first_cycle_checkpoint)
            fprintf('[BGC resume] cycle=1 KT=%d/%d\n',firstKt-1,nSteps);
        end
        for kt=firstKt:nSteps
            [model,report]=runOneStep(model,kt);
            if report.step_failed
                error('runCOMPASSBiogeochemistrySpinup:stepFailed', ...
                    'Cycle %d rejected KT=%d: %s',cycle,kt, ...
                    char(report.termination_reason));
            end
            if mod(kt,options.progress_interval)==0 || kt==nSteps
                printProgress(model,cycle,kt,nSteps);
            end
        end
        cycleAudits{cycle}=validateBiogeochemistryCycle( ...
            model,cycleInitial,nSteps,cycle);
        cycleAudits{cycle}.resumed_at_kt=firstKt-1;
        cycleAudits{cycle}.full_configured_window=(nSteps==configuredSteps);
        model=finalizeRun(model);
        carriedY=model.built.Y;
        history{end+1}=carriedY; %#ok<AGROW>
        if cycle>=2
            convergence{cycle}=soil.organic.maxRelativePoolCarbonChange( ...
                history{end-1},history{end},options.absolute_floor);
            convergence{cycle}.tolerance=options.relative_tolerance;
            convergence{cycle}.converged= ...
                convergence{cycle}.max_relative_change<=options.relative_tolerance ...
                &&nSteps==configuredSteps;
        else
            convergence{cycle}=struct('max_relative_change',Inf, ...
                'state_path','','element_index',NaN,'C_start',NaN,'C_end',NaN, ...
                'absolute_floor_native_units',options.absolute_floor, ...
                'n_compared_states',0,'tolerance',options.relative_tolerance, ...
                'converged',false);
        end
        canAccelerate=options.slow_pool_acceleration&&cycle>=options.acceleration_start_cycle ...
            &&cycle<=options.cycles-options.final_rebalance_cycles&&numel(history)>=3;
        if canAccelerate
            [carriedY,acceleration{cycle}]=soil.organic.accelerateSlowPools( ...
                history{end-2},history{end-1},history{end}, ...
                options.acceleration_max_factor);
            if acceleration{cycle}.applied
                history{end}=carriedY; lastAccelerationCycle=cycle; stableCount=0;
            end
        else
            acceleration{cycle}=struct('applied',false,'accelerated_elements',0, ...
                'method','none','maximum_factor',1);
        end
        model.built.Y=carriedY;
        naturalCycles=cycle-lastAccelerationCycle;
        observerConverged=true;
        if ~isempty(options.cycle_observer)
            observerConverged=options.cycle_observer(model,cycleInitial, ...
                cycle,outputRoot,convergence{cycle},cycleAudits{cycle});
            assert(islogical(observerConverged)&&isscalar(observerConverged), ...
                'runCOMPASSBiogeochemistrySpinup:observerGate', ...
                'Cycle observer must return one logical convergence gate.');
        end
        if convergence{cycle}.converged&&observerConverged ...
                &&naturalCycles>=options.final_rebalance_cycles
            stableCount=stableCount+1;
        else
            stableCount=0;
        end
        carry=init.captureBiogeochemistryCycleCarry(model,cycle);
        assert(isequaln(carry.identity,identity), ...
            'runCOMPASSBiogeochemistrySpinup:identity', ...
            'Immutable model identity changed while integrating the cycle.');
        checkpointFiles(cycle)=string(io.saveCompassCheckpoint(model,nSteps));
        completedCycles=cycle;
        fprintf(['[BGC convergence] cycle=%d max_rel=%.8g path=%s ', ...
            'stable=%d/%d accelerated=%d rebalance=%d/%d\n'],cycle, ...
            convergence{cycle}.max_relative_change,convergence{cycle}.state_path, ...
            stableCount,options.consecutive_cycles,acceleration{cycle}.applied, ...
            naturalCycles,options.final_rebalance_cycles);
        if options.auto_stop&&cycle>=options.minimum_cycles ...
                &&stableCount>=options.consecutive_cycles
            stopReason='periodic_pool_convergence';
            break;
        end
    end
    checkpointFiles=checkpointFiles(1:completedCycles);
    cycleAudits=cycleAudits(1:completedCycles);
    convergence=convergence(1:completedCycles);
    acceleration=acceleration(1:completedCycles);
    result=struct('schema_version',3,'cycles',completedCycles, ...
        'maximum_cycles',options.cycles,'cycle_steps',nSteps, ...
        'configured_window_steps',configuredSteps, ...
        'periodic_convergence_eligible',nSteps==configuredSteps, ...
        'cycle_identity',identity, ...
        'checkpoint_files',checkpointFiles, ...
        'checkpoint_file',checkpointFiles(end),'cycle_audits',{cycleAudits}, ...
        'convergence',{convergence},'acceleration',{acceleration}, ...
        'stop_reason',stopReason,'final_state',carry.Y, ...
        'production_equations_modified',false, ...
        'acceleration_method','bounded_periodic_aitken_slow_pools', ...
        'colm_sasu_exact',false);
    save(fullfile(outputRoot,'biogeochemistry_spinup_result.mat'),'result','-v7.3');
end

function audit=validateBiogeochemistryCycle(model,initialY,nSteps,cycle)
% Hard scientific acceptance gates for a completed spin-up cycle.
% Editor: Jizixing Zhu; Date: 2026-08-24.
    Y=model.built.Y;
    audit=struct('cycle',cycle,'steps',nSteps);
    assert(allFiniteNonnegative(Y.soil_dom), ...
        'runCOMPASSBiogeochemistrySpinup:invalidDOM', ...
        'Cycle %d ended with non-finite or negative DOM state.',cycle);
    assert(allFiniteNonnegative(Y.soil_microbes), ...
        'runCOMPASSBiogeochemistrySpinup:invalidMicrobes', ...
        'Cycle %d ended with non-finite or negative microbial state.',cycle);

    pc=Y.plant_cn; a=pc.accum; a0=initialY.plant_cn.accum;
    physicalC=rmfield(pc.C,'xsmr');
    assert(allFiniteNonnegative(physicalC)&&allFiniteNonnegative(pc.N), ...
        'runCOMPASSBiogeochemistrySpinup:invalidPlantPools', ...
        'Cycle %d ended with non-finite or negative physical PlantCN pools.',cycle);
    assert(isfinite(pc.C.xsmr)&&pc.C.xsmr>=-1e-8, ...
        'runCOMPASSBiogeochemistrySpinup:plantCarbonDebt', ...
        'Cycle %d ended with unphysical C_xsmr=%g g C m^-2.',cycle,pc.C.xsmr);

    cNow=sumFields(pc.C)+sumFields(pc.transfer.C);
    nNow=sumFields(pc.N)+sumFields(pc.transfer.N);
    cLoss=a.cum_C_harvest+a.cum_C_leaf_litter+a.cum_C_root_litter ...
        +a.cum_C_prune_litter+a.cum_C_mortality;
    nLoss=a.cum_N_harvest+a.cum_N_leaf_litter+a.cum_N_root_litter ...
        +a.cum_N_prune_litter+a.cum_N_mortality;
    audit.plant_c_residual_g_m2=cNow+cLoss+a.cum_R_a ...
        -a.cum_C_xsmr_writeoff ...
        -pc.baseline.C_pool_init_const-a.cum_GPP;
    audit.plant_n_residual_g_m2=nNow+nLoss ...
        -pc.baseline.N_pool_init_const-a.cum_U_N;
    cTolerance=1e-8*max(1,pc.baseline.C_pool_init_const+a.cum_GPP);
    nTolerance=1e-8*max(1,pc.baseline.N_pool_init_const+a.cum_U_N);
    assert(abs(audit.plant_c_residual_g_m2)<=cTolerance, ...
        'runCOMPASSBiogeochemistrySpinup:plantCarbonClosure', ...
        'Cycle %d PlantCN carbon residual is %g g C m^-2.', ...
        cycle,audit.plant_c_residual_g_m2);
    assert(abs(audit.plant_n_residual_g_m2)<=nTolerance, ...
        'runCOMPASSBiogeochemistrySpinup:plantNitrogenClosure', ...
        'Cycle %d PlantCN nitrogen residual is %g g N m^-2.', ...
        cycle,audit.plant_n_residual_g_m2);

    audit.gpp_g_m2=accumDelta(a,a0,'cum_GPP');
    audit.n_uptake_g_m2=accumDelta(a,a0,'cum_U_N');
    audit.harvest_c_g_m2=accumDelta(a,a0,'cum_C_harvest');
    audit.harvest_n_g_m2=accumDelta(a,a0,'cum_N_harvest');
    audit.leaf_litter_c_g_m2=accumDelta(a,a0,'cum_C_leaf_litter');
    audit.root_litter_c_g_m2=accumDelta(a,a0,'cum_C_root_litter');
    audit.prune_litter_c_g_m2=accumDelta(a,a0,'cum_C_prune_litter');
    audit.mortality_litter_c_g_m2=accumDelta(a,a0,'cum_C_mortality');
    audit.leaf_litter_n_g_m2=accumDelta(a,a0,'cum_N_leaf_litter');
    audit.root_litter_n_g_m2=accumDelta(a,a0,'cum_N_root_litter');
    audit.prune_litter_n_g_m2=accumDelta(a,a0,'cum_N_prune_litter');
    audit.mortality_litter_n_g_m2=accumDelta(a,a0,'cum_N_mortality');
    audit.litter_c_g_m2=audit.leaf_litter_c_g_m2+audit.root_litter_c_g_m2 ...
        +audit.prune_litter_c_g_m2+audit.mortality_litter_c_g_m2;
    audit.litter_n_g_m2=audit.leaf_litter_n_g_m2+audit.root_litter_n_g_m2 ...
        +audit.prune_litter_n_g_m2+audit.mortality_litter_n_g_m2;
    if nSteps>=8760 && model.built.control.switches.lLitterInput==1
        assert(audit.litter_c_g_m2>0&&audit.litter_n_g_m2>0, ...
            'runCOMPASSBiogeochemistrySpinup:missingPlantLitter', ...
            ['Cycle %d completed a full year with litter_input enabled but ', ...
             'PlantCN returned no C/N litter.'],cycle);
    end
    pftType=char(string(model.built.parameters.veg.structure.PFT_type));
    if nSteps>=8760 && cycle>=2 && strcmpi(pftType,'Apple')
        assert(audit.n_uptake_g_m2>0, ...
            'runCOMPASSBiogeochemistrySpinup:missingAppleNitrogenUptake', ...
            'Cycle %d Apple PlantCN completed a full year with zero root N uptake.',cycle);
        assert(audit.leaf_litter_c_g_m2>0&&audit.root_litter_c_g_m2>0, ...
            'runCOMPASSBiogeochemistrySpinup:missingAppleTissueTurnover', ...
            ['Cycle %d Apple PlantCN completed a preconditioned full year ', ...
             'without both leaf and fine-root C litter.'],cycle);
    end
    if nSteps>=8760 && strcmpi(pftType,'Winter Wheat')
        assert(audit.n_uptake_g_m2>0, ...
            'runCOMPASSBiogeochemistrySpinup:missingWheatNitrogenUptake', ...
            'Cycle %d winter wheat completed a full crop year with zero root N uptake.',cycle);
        assert(audit.harvest_c_g_m2>0&&audit.harvest_n_g_m2>0, ...
            'runCOMPASSBiogeochemistrySpinup:missingWheatHarvest', ...
            'Cycle %d winter wheat completed a full crop year without C/N grain harvest.',cycle);
        assert(audit.leaf_litter_c_g_m2>=0&&audit.root_litter_c_g_m2>0, ...
            'runCOMPASSBiogeochemistrySpinup:missingWheatResidue', ...
            'Cycle %d winter wheat has negative shoot litter or no positive root residue.',cycle);
        ps=Y.vegetation.phenology.state;
        assert(~ps.crop_live&&ps.crop_phase==4, ...
            'runCOMPASSBiogeochemistrySpinup:wheatCalendarIncomplete', ...
            'Cycle %d ended before the winter-wheat harvest state was reached.',cycle);
        woodyC=pc.C.livestem+pc.C.deadstem+pc.C.livecroot+pc.C.deadcroot;
        storageC=pc.C.root_storage+pc.C.stem_storage;
        assert(abs(woodyC)<=cTolerance&&abs(storageC)<=cTolerance, ...
            'runCOMPASSBiogeochemistrySpinup:wheatPerennialPools', ...
            'Cycle %d winter wheat retained perennial wood/storage C.',cycle);
    end

    dz=model.built.domain.dz(:);
    audit.doc_inventory_g_cm2=sum(Y.soil_dom.doc_total_g_cm3(:).*dz);
    audit.don_inventory_g_cm2=sum(Y.soil_dom.don_total_g_cm3(:).*dz);
    audit.microbial_biomass_min_mg_g=min(Y.soil_microbes.MB(:));
    audit.microbial_biomass_max_mg_g=max(Y.soil_microbes.MB(:));
    mb0=sum(initialY.soil_microbes.MB(:));
    mb1=sum(Y.soil_microbes.MB(:));
    audit.microbial_biomass_relative=mb1/max(mb0,eps);
    if audit.microbial_biomass_relative<0.01
        warning('runCOMPASSBiogeochemistrySpinup:microbialCollapse', ...
            'Cycle %d microbial biomass fell below 1%% of its cycle-start sum.',cycle);
    end
    fprintf(['[BGC audit] cycle=%d litterC=%.8g gC m-2 litterN=%.8g gN m-2 ', ...
        'Cres=%+.3e Nres=%+.3e xsmr=%.8g DON=%.8g gN cm-2\n'], ...
        cycle,audit.litter_c_g_m2,audit.litter_n_g_m2, ...
        audit.plant_c_residual_g_m2,audit.plant_n_residual_g_m2, ...
        pc.C.xsmr,audit.don_inventory_g_cm2);
end

function tf=allFiniteNonnegative(s)
% True when every numeric state field is finite and nonnegative.
% Editor: Jizixing Zhu; Date: 2026-08-24.
    f=fieldnames(s); tf=true;
    for i=1:numel(f)
        value=s.(f{i});
        if isnumeric(value) && (~all(isfinite(value(:)))||any(value(:)<-1e-12))
            tf=false; return;
        end
    end
end

function value=sumFields(s)
% Sum numeric fields in a state struct.
% Editor: Jizixing Zhu; Date: 2026-08-24.
    f=fieldnames(s); value=0;
    for i=1:numel(f), value=value+sum(s.(f{i})(:)); end
end

function value=accumDelta(current,initial,name)
% Return one cycle-local cumulative PlantCN flux increment.
% Editor: Jizixing Zhu; Date: 2026-08-24.
    value=current.(name)-initial.(name);
    if abs(value)<1e-12, value=0; end
end

function printProgress(model,cycle,kt,nSteps)
% Print accepted-state biogeochemistry telemetry.
% Editor: Jizixing Zhu; Date: 2026-08-21.
% Accepted-state telemetry for long spin-up diagnostics; no state mutation.
    Y=model.built.Y; dz=model.built.domain.dz(:);
    docInventory=sum(Y.soil_dom.doc_total_g_cm3(:).*dz);
    donInventory=sum(Y.soil_dom.don_total_g_cm3(:).*dz);
    mb=Y.soil_microbes.MB(:);
    fprintf(['[BGC spin-up] cycle=%d KT=%d/%d DOC=%.8g gC cm-2 ', ...
        'DON=%.8g gN cm-2 MB=[%.8g,%.8g] mgC g-1\n'], ...
        cycle,kt,nSteps,docInventory,donInventory,min(mb),max(mb));
end

function restoreEnvironment(names,values)
% Restore process-local COMPASS run environment variables.
% Editor: Jizixing Zhu; Date: 2026-08-21.
    for i=1:numel(names),setenv(names{i},values{i});end
end

% Change: 2026-09-02 Jizixing Zhu -- include xsmr debt writeoff in the
% spin-up plant-C audit without classifying it as atmospheric respiration.
