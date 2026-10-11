function artifact=exportMeasuredBiogeochemistryInitialArtifact(fresh,donor,checkpointFile,artifactFile,evidence,options)
%EXPORTMEASUREDBIOGEOCHEMISTRYINITIALARTIFACT Bind a measured t0 prior to evidence.
% FRESH is an initialized/unprepared model, DONOR an accepted model carrying
% its pre-integration biogeochemistry_spinup_identity (or cycle identity).
% CHECKPOINTFILE must contain that exact accepted Y and identity. ARTIFACTFILE
% is a new MAT file, never overwritten. EVIDENCE contains optional JSON paths
% source_audit_file and convergence_file. Each release attestation requires
% status PASS, identity_signature, donor_checkpoint_sha256 and its named
% source_audit_verified/convergence_verified logical flag. This adapter binds
% and preserves attestations; it does not compute convergence or closure.
% Missing evidence produces an explicitly NOT_READY artifact. Periodic phase
% plus both verified attestations are required for scientific release.
% Organic states retain production units/orientation. Initialization C/N
% adjustments are g/cm2, outside all physical flux ledgers. Measured water,
% DIN and SOC remain the fresh targets; no time is advanced.
% Editor: Wenjie Liu; Date: 2026-09-30.
    arguments
        fresh (1,1) struct
        donor (1,1) struct
        checkpointFile (1,:) char
        artifactFile (1,:) char
        evidence (1,1) struct=struct()
        options.donor_kind (1,:) char {mustBeMember(options.donor_kind, ...
            {'periodic_state','diagnostic_state'})}='periodic_state'
    end
    prefix='exportMeasuredInitialArtifact';
    assert(~isfile(artifactFile),[prefix ':existingArtifact'], ...
        'Preserve existing artifacts; use a new filename.');
    assertMeasuredWater(fresh,prefix);
    identity=init.buildBiogeochemistryCycleIdentity(fresh);
    assert(isequaln(identity,init.buildBiogeochemistryCycleIdentity(donor)), ...
        [prefix ':identity'],'Donor and target immutable identities differ.');
    assert(isequaln(identity,capturedIdentity(donor.runtime,prefix)), ...
        [prefix ':identity'],'Donor identity changed after integration began.');
    checkpoint=io.loadCompassCheckpoint(checkpointFile);
    assert(isequaln(checkpoint.Y,donor.built.Y) ...
        &&checkpoint.kt==donor.runtime.run_control.last_completed_kt ...
        &&isequaln(checkpoint.runtime.step,donor.runtime.step) ...
        &&isequaln(identity,capturedIdentity(checkpoint.runtime,prefix)), ...
        [prefix ':checkpoint'],'Checkpoint is not this accepted donor and identity.');
    checkpointHash=fileHash(checkpointFile);
    [conditioned,adjustments]=init.materializeMeasuredBiogeochemistryInitialState( ...
        fresh,donor,donor_kind=options.donor_kind,donor_label=checkpointFile);
    owners={'soil_microbes','soil_organic','soil_dom','soil_fertilizer_n_tracer'};
    state=struct();
    for i=1:numel(owners)
        if isfield(conditioned.built.Y,owners{i})
            state.(owners{i})=conditioned.built.Y.(owners{i});
        end
    end
    source=attestation(evidence,'source_audit_file','source_audit_verified', ...
        identity.signature,checkpointHash,prefix);
    convergence=attestation(evidence,'convergence_file','convergence_verified', ...
        identity.signature,checkpointHash,prefix);
    periodic=strcmp(options.donor_kind,'periodic_state');
    ready=periodic&&source.verified&&convergence.verified;
    reasons=strings(0,1);
    if ~periodic,reasons(end+1)="donor_is_nonperiodic_diagnostic";end
    if ~source.verified,reasons(end+1)="verified_source_audit_missing";end
    if ~convergence.verified,reasons(end+1)="verified_convergence_evidence_missing";end
    policy=struct('method','measured_SOC_conditioned_organic_prior', ...
        'measured_water_DIN_preserved',true,'measured_SOC_preserved',true, ...
        'parameter_transfer','exact_identity_only_re_equilibrate_other_parameters', ...
        'AOM','donor_unreacted_AOM_outside_measured_SOC', ...
        'IOM','fresh_fixed_prior','clock','fresh_t0', ...
        'cumulative_ledgers','fresh_zero','fertilizer_tracer','fresh_period_only');
    artifact=struct('schema_version',1,'kind','measured_biogeochemistry_initial_artifact', ...
        'identity',identity,'measured_target_sha256',measuredHash(fresh), ...
        'donor_checkpoint',struct('path',checkpointFile,'sha256',checkpointHash, ...
            'accepted_kt',checkpoint.kt,'clock_time',char(checkpoint.runtime.step.clock_time)), ...
        'policy',policy,'donor_kind',options.donor_kind,'state',state, ...
        'adjustment_audit',adjustments,'source_audit',source,'convergence',convergence, ...
        'scientific_release',ready,'not_ready_reasons',reasons, ...
        'created_at',char(datetime('now',Format='yyyy-MM-dd HH:mm:ss')));
    artifact.signature=hashJSON(artifact);
    parent=fileparts(artifactFile);
    assert(isempty(parent)||isfolder(parent),[prefix ':directory'], ...
        'Create the explicit output directory before exporting.');
    save(artifactFile,'artifact','-v7.3');
end

function a=attestation(evidence,key,flag,identity,checkpointHash,prefix)
    a=struct('path','','sha256','','verified',false,'report',struct());
    if ~isfield(evidence,key)||strlength(string(evidence.(key)))==0,return;end
    file=char(evidence.(key));
    assert(isfile(file),[prefix ':evidenceMissing'],'Missing evidence file %s.',file);
    r=jsondecode(fileread(file));
    required={'status','identity_signature','donor_checkpoint_sha256',flag};
    assert(isstruct(r)&&isscalar(r)&&all(isfield(r,required)) ...
        &&islogical(r.(flag))&&isscalar(r.(flag)), ...
        [prefix ':evidenceSchema'],'Evidence %s does not meet its attestation contract.',file);
    assert(strcmp(r.identity_signature,identity) ...
        &&strcmp(r.donor_checkpoint_sha256,checkpointHash), ...
        [prefix ':evidenceIdentity'],'Evidence is from a different donor or model.');
    a=struct('path',file,'sha256',fileHash(file), ...
        'verified',strcmp(r.status,'PASS')&&r.(flag),'report',r);
end

function identity=capturedIdentity(runtime,prefix)
    if isfield(runtime,'biogeochemistry_spinup_identity')
        identity=runtime.biogeochemistry_spinup_identity;
    elseif isfield(runtime,'biogeochemistry_cycle')
        identity=runtime.biogeochemistry_cycle.identity;
    else
        error([prefix ':missingCapturedIdentity'], ...
            'Donor must carry identity captured before integration; cannot infer it retrospectively.');
    end
end

function assertMeasuredWater(model,prefix)
    assert(isfield(model,'status')&&string(model.status)=="initialized" ...
        &&~isfield(model,'runtime'),[prefix ':notFresh'],'Target must be fresh before runtime.');
    b=model.built;a=b.p.initialization.measured_soil;
    assert(a.applied&&a.measured_water_verified&&a.measured_DIN_verified, ...
        [prefix ':measuredProtocol'],'Measured water and total DIN conditioning is required.');
    dz=b.domain.dz(:);u=b.Y.Theta_U;
    theta=[u(1,1);(dz(1:end-1).*u(1:end-1,2)+dz(2:end).*u(2:end,1)) ...
        ./(dz(1:end-1)+dz(2:end));u(end,2)];
    target=init.interpLinearWithExtrapolation(model.raw.initial.depth_cm, ...
        model.raw.initial.theta,b.domain.z_node);
    assert(all(isfinite(theta))&&max(abs(theta-target))<1e-10 ...
        &&max(abs(b.Y.theta(:)-target(:)))<1e-10, ...
        [prefix ':measuredWater'],'Fresh water differs from measured targets.');
end

function value=measuredHash(model)
    input=model.raw.initial;target=struct();
    for key={'depth_cm','theta','nh4','no3','soc_total_mgkg'}
        target.(key{1})=input.(key{1});
    end
    target.nitrogen_input_phase=model.raw.config.switches.nitrogen_input_phase;
    value=hashJSON(target);
end
function value=hashJSON(value)
    value=digest(unicode2native(jsonencode(value,ConvertInfAndNaN=false),'UTF-8'));
end
function value=fileHash(file)
    fid=fopen(file,'rb');assert(fid>=0,'exportMeasuredInitialArtifact:read','Cannot read %s.',file);
    clean=onCleanup(@()fclose(fid));
    value=digest(fread(fid,Inf,'*uint8'));
end
function value=digest(bytes)
    engine=java.security.MessageDigest.getInstance('SHA-256');engine.update(uint8(bytes(:)));
    value=upper(reshape(dec2hex(typecast(engine.digest(),'uint8'),2).',1,[]));
end
