function [model,audit]=applyMeasuredBiogeochemistryInitialArtifact(fresh,artifactFile,options)
%APPLYMEASUREDBIOGEOCHEMISTRYINITIALARTIFACT Load a bound organic prior at new t0.
% FRESH is initialized/unprepared with measured water, total DIN and SOC.
% ARTIFACTFILE is created by exportMeasuredBiogeochemistryInitialArtifact.
% The exact source/parameter/grid/site/forcing/library identity, measured
% target hash, checkpoint hash and all supplied evidence hashes must match.
% No source-change or parameter-nearness bypass exists. Unreleased artifacts
% are rejected unless allow_diagnostic=true; this opt-in remains labelled
% scientific_release=false. It is intended for bounded technical checks.
% State owners retain production units/grid order; all other fresh Y,
% physical flux ledgers, runtime clock and management remain fresh. The SOC
% conditioner independently rechecks pool ownership and measured constraints.
% It does not certify convergence merely because this file can be applied.
% Editor: Wenjie Liu; Date: 2026-09-30.
    arguments
        fresh (1,1) struct
        artifactFile (1,:) char
        options.allow_diagnostic (1,1) logical=false
    end
    prefix='applyMeasuredInitialArtifact';
    assert(isfield(fresh,'status')&&string(fresh.status)=="initialized" ...
        &&~isfield(fresh,'runtime'),[prefix ':notFresh'],'Use fresh pre-runtime target.');
    data=load(artifactFile,'artifact');
    assert(isfield(data,'artifact'),[prefix ':schema'],'MAT file lacks artifact.');
    a=data.artifact;
    required={'schema_version','kind','identity','measured_target_sha256', ...
        'donor_checkpoint','policy','donor_kind','state','adjustment_audit', ...
        'source_audit','convergence','scientific_release','signature'};
    assert(isstruct(a)&&isscalar(a)&&all(isfield(a,required)) ...
        &&a.schema_version==1&&strcmp(a.kind,'measured_biogeochemistry_initial_artifact'), ...
        [prefix ':schema'],'Unknown measured initial artifact schema.');
    assert(strcmp(hashJSON(rmfield(a,'signature')),a.signature), ...
        [prefix ':integrity'],'Artifact metadata or state changed after export.');
    identity=init.buildBiogeochemistryCycleIdentity(fresh);
    assert(isequaln(identity,a.identity),[prefix ':identity'], ...
        'Source, parameters, site, forcing or library changed; regenerate and re-equilibrate.');
    assert(strcmp(measuredHash(fresh),a.measured_target_sha256), ...
        [prefix ':measuredTarget'],'Measured t0 target differs from exported target.');
    assert(strcmp(fileHash(a.donor_checkpoint.path),a.donor_checkpoint.sha256), ...
        [prefix ':checkpoint'],'The provenance checkpoint changed.');
    validateEvidence(a.source_audit,'source_audit_verified',a,prefix);
    validateEvidence(a.convergence,'convergence_verified',a,prefix);
    ready=strcmp(a.donor_kind,'periodic_state') ...
        &&a.source_audit.verified&&a.convergence.verified;
    assert(islogical(a.scientific_release)&&isscalar(a.scientific_release) ...
        &&a.scientific_release==ready,[prefix ':releaseFlag'],'Inconsistent release flag.');
    assert(ready||options.allow_diagnostic,[prefix ':notReady'], ...
        'Spinup source audit/convergence not verified; production application is refused.');
    assertMeasuredWater(fresh,prefix);
    before=fresh.built.Y;model=fresh;
    owners={'soil_microbes','soil_organic','soil_dom','soil_fertilizer_n_tracer'};
    assert(all(ismember(fieldnames(a.state),owners)) ...
        &&all(isfield(a.state,{'soil_microbes','soil_organic','soil_dom'})), ...
        [prefix ':stateOwners'],'Artifact may contain only declared organic/tracer owners.');
    for key=fieldnames(a.state)'
        assert(isfield(before,key{1}),[prefix ':stateOwners'],'Unknown target state owner.');
        model.built.Y.(key{1})=a.state.(key{1});
    end
    % A self-conditioning validation checks positive pools, aliases, measured
    % SOC/DIN and zero fresh ledgers. Retain exact exported bits afterwards.
    [~,check]=init.materializeMeasuredBiogeochemistryInitialState( ...
        model,model,donor_kind='initial_guess',donor_label=artifactFile);
    assert(~check.physical_flux_added,[prefix ':flux'],'Validation introduced a flux.');
    for key=fieldnames(before)'
        if ~ismember(key{1},owners)
            assert(isequaln(before.(key{1}),model.built.Y.(key{1})), ...
                [prefix ':protectedState'],'Protected fresh state changed.');
        end
    end
    y=model.built.Y;
    model.built=init.rebuildSoilAuxiliaryFromAccepted(model.built,0);
    assert(isequaln(model.built.Y,y),[prefix ':stateChanged'],'Cache rebuild changed state.');
    audit=struct('schema_version',1,'artifact_path',artifactFile, ...
        'artifact_file_sha256',fileHash(artifactFile),'artifact_signature',a.signature, ...
        'identity_signature',identity.signature,'measured_target_sha256',a.measured_target_sha256, ...
        'donor_checkpoint_sha256',a.donor_checkpoint.sha256, ...
        'source_audit_verified',a.source_audit.verified, ...
        'spinup_convergence_verified',a.convergence.verified, ...
        'scientific_release',ready,'diagnostic_opt_in',options.allow_diagnostic, ...
        'physical_flux_added',false,'clock','fresh_unprepared_t0', ...
        'adjustment_audit',a.adjustment_audit);
    model.built.initialization.measured_biogeochemistry=a.adjustment_audit;
    model.built.initialization.measured_initial_artifact=audit;
end

function validateEvidence(e,key,a,prefix)
    if isempty(e.path)
        assert(~e.verified,[prefix ':evidence'],'Missing evidence cannot be verified.');return;
    end
    assert(strcmp(fileHash(e.path),e.sha256),[prefix ':evidence'],'Evidence file changed.');
    r=jsondecode(fileread(e.path));
    assert(isequaln(r,e.report)&&isfield(r,key) ...
        &&strcmp(r.identity_signature,a.identity.signature) ...
        &&strcmp(r.donor_checkpoint_sha256,a.donor_checkpoint.sha256) ...
        &&e.verified==(strcmp(r.status,'PASS')&&r.(key)), ...
        [prefix ':evidence'],'Evidence no longer matches its donor and identity.');
end
function assertMeasuredWater(model,prefix)
    b=model.built;a=b.p.initialization.measured_soil;
    assert(a.applied&&a.measured_water_verified&&a.measured_DIN_verified, ...
        [prefix ':measuredProtocol'],'Measured initial water/DIN protocol required.');
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
    fid=fopen(file,'rb');assert(fid>=0,'applyMeasuredInitialArtifact:read','Cannot read %s.',file);
    clean=onCleanup(@()fclose(fid));
    value=digest(fread(fid,Inf,'*uint8'));
end
function value=digest(bytes)
    engine=java.security.MessageDigest.getInstance('SHA-256');engine.update(uint8(bytes(:)));
    value=upper(reshape(dec2hex(typecast(engine.digest(),'uint8'),2).',1,[]));
end
