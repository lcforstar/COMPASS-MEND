function [built,audit]=conditionMeasuredSoilInitialState(built,raw)
%CONDITIONMEASUREDSOILINITIALSTATE Zero-step warm field-profile conditioning.
% Called once after auxiliary and thermal-bottom initialization, before runtime.
% Raw theta is volumetric total water [-], NH4/NO3 total [mg/kg dry soil] or
% dissolved [mg/L] according to nitrogen_input_phase. Profiles are linearly
% interpolated/extrapolated to bottom-up nodes; this creates no new observations.
% Invert the actual endpoint hydraulic kernel (including temperature correction)
% to match half-element-volume projected water. Rebuild U/auxiliaries, then
% convert total DIN with the accepted FE effective-liquid capacities and dry
% mass from soil.nitrogen.observationProfile. No in-run reset or applied flux.
% SOC/DOM/microbial stores and all boundary/kinetic parameters remain unchanged.
% Scoped to one unfrozen field column with volumetric input; cold/laboratory/
% regional or nonvolumetric inputs retain existing behavior with explicit audit.
% Out-of-range/infeasible observations fail; no clipping to another initial value.
% Editor: Wenjie Liu; Date: 2026-09-29.
    audit=struct('schema_version',1,'applied',false,'reason','unsupported_profile', ...
        'measured_water_verified',false,'measured_DIN_verified',false);
    if size(built.Y.T,2)~=1||any(built.Y.T(:)<=0) ...
        ||(isfield(built.control,'simulation')&&built.control.simulation.is_laboratory)
        audit.reason='cold_laboratory_or_multicolumn_not_conditioned';return;
    end
    input=raw.initial;
    if ~all(isfinite(input.theta))||~all(input.theta>=0&input.theta<=1)
        audit.reason='nonvolumetric_water_input_not_conditioned';return;
    end
    assert(all(built.Y.soil_nitrogen.flux.cvnh4(:)==0) ...
        &&all(built.Y.soil_nitrogen.flux.cvno3(:)==0), ...
        'conditionMeasuredSoil:notFresh','Do not condition an integrated or resumed state.');
    dz=built.domain.dz(:); z=built.domain.z_node(:); n=numel(z);
    target=init.interpLinearWithExtrapolation(input.depth_cm,input.theta,z);
    assert(all(isfinite(target))&&all(target>=0&target<=1), ...
        'conditionMeasuredSoil:waterInput','Invalid interpolated volumetric water.');
    lower=-1e7*ones(n,1); upper=zeros(n,1);
    dry=waterAt(built,lower,dz); wet=waterAt(built,upper,dz);
    assert(all(target>=dry-1e-12&target<=wet+1e-12), ...
        'conditionMeasuredSoil:infeasibleWater', ...
        'Measured water is outside the active layered retention range.');
    for iteration=1:70
        middle=(lower+upper)/2;
        theta=waterAt(built,middle,dz);
        below=theta<target; lower(below)=middle(below); upper(~below)=middle(~below);
    end
    head=(lower+upper)/2; [theta,waterBuilt]=waterAt(built,head,dz);
    assert(max(abs(theta-target))<1e-10, ...
        'conditionMeasuredSoil:waterInversion','Water inversion did not converge.');
    y=built.Y; built=waterBuilt; built.Y=y;
    built.Y.h=head; built.Y.theta=target;
    built.Y.Theta_L=built.workspace.derived.Theta_LL;
    built.Y.Theta_U=built.workspace.derived.Theta_UU;
    built.Y.Theta_I=built.workspace.derived.Theta_II;
    built.workspace.step.trial.h=head; built.workspace.step.previous.h=head;
    % Re-seed only fresh N memory after the final water state has been rebuilt.
    built.Y.soil_nitrogen=soil.nitrogen.emptyAcceptedState();
    built=init.initializeSoilAux(built);
    built.Y.Theta_L=built.workspace.derived.Theta_LL;
    built.Y.Theta_U=built.workspace.derived.Theta_UU;
    built.Y.Theta_I=built.workspace.derived.Theta_II;
    assert(all(built.Y.Theta_I(:)==0), ...
        'conditionMeasuredSoil:frozen','Unfrozen initialization contract violated.');
    finalTheta=project(built.Y.Theta_U,dz);
    assert(max(abs(finalTheta-target))<1e-10, ...
        'conditionMeasuredSoil:waterRefresh','Water refresh changed the measured target.');
    [~,~,capacity]=soil.nitrogen.observationProfile( ...
        built.Y,built.parameters.soil_n,dz);
    phase=string(raw.config.switches.nitrogen_input_phase);
    expected=[init.interpLinearWithExtrapolation(input.depth_cm,input.nh4,z), ...
        init.interpLinearWithExtrapolation(input.depth_cm,input.no3,z)];
    assert(all(isfinite(expected),'all')&&all(expected>=0,'all'), ...
        'conditionMeasuredSoil:DINInput','Invalid interpolated measured DIN.');
    if phase=="total"
        assert(all(capacity.nh4_capacity_cm>0)&&all(capacity.no3_capacity_cm>0), ...
            'conditionMeasuredSoil:dryDIN','Positive capacities required for measured total DIN.');
        built.Y.nh4=expected(:,1).*capacity.dry_mass_g_cm2./capacity.nh4_capacity_cm;
        built.Y.no3=expected(:,2).*capacity.dry_mass_g_cm2./capacity.no3_capacity_cm;
        [nh4,no3]=soil.nitrogen.observationProfile(built.Y,built.parameters.soil_n,dz);
        dinError=max(abs([nh4,no3]-expected),[],1);
        assert(all(dinError<1e-10),'conditionMeasuredSoil:DINRoundtrip','DIN roundtrip failed.');
        audit.measured_DIN_verified=true;
        audit.DIN_roundtrip_max_abs_mgkg=dinError;
    elseif phase=="dissolved"
        built.Y.nh4=expected(:,1);built.Y.no3=expected(:,2);
    else
        error('conditionMeasuredSoil:DINPhase','Unsupported DIN input phase.');
    end
    % N concentration mirrors are disposable; rebuild them from the new Y.
    built=init.rebuildSoilAuxiliaryFromAccepted(built,0);
    audit.applied=true;audit.reason='warm_field_measured_profile_conditioned';
    audit.measured_water_verified=true;
    audit.water_max_abs_m3m3=max(abs(finalTheta-target));
    audit.water_target_bottom_up=target;audit.physical_flux_added=false;
    audit.DIN_observable=capacity.definition;
    built.p.initialization.measured_soil=audit;
end

function [theta,b]=waterAt(b,head,dz)
    b.workspace.step.trial.h=head;
    lf=0;
    if isfield(b.control.switches,'lFreezing')&&b.control.switches.lFreezing==1
        lf=equations.getPhysicalConstant('lambda_fus');
        override=str2double(io.getEnv('LFUSION_OVERRIDE'));
        if isfinite(override),lf=override;end
    end
    b=soil.water.updateWaterContentAuxiliary(b,0,lf);
    theta=project(b.workspace.derived.Theta_UU,dz);
end

function value=project(endpoint,dz)
    mass=endpoint.*(dz/2);
    volume=[dz/2;0]+[0;dz/2];
    value=([mass(:,1);0]+[0;mass(:,2)])./volume;
end
