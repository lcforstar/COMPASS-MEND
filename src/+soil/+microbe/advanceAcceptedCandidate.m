function candidate = advanceAcceptedCandidate( ...
        accepted, domAccepted, thetaBottomUp, bulkDensityBottomUp, ...
        domParameters, parameters, dt_s, temperatureBottomUp, headBottomUp)
%ADVANCEACCEPTEDCANDIDATE Advance the complete MEND carbon core with RK4.
%
% Solid pools are top-down [mg C g-1 dry soil]. Total DOM is owned by
% soil_dom in bottom-up [g C cm-3 bulk soil]. The internal RK DOC coordinate
% is total DOM [mg C g-1], while ledger stages expose aqueous DOC for MEND
% kinetics and the paired N solver. No accepted input is mutated.
%
% Editor: Jizixing Zhu; Date: 2026-08-21.
% Change: 2026-09-12 Jizixing Zhu -- replace MEMC plus detached dormancy by
%   the source MEND POM/MOM/DOM/MBA/MBD/enzyme RK system.
% Change: 2026-09-13 Codex -- integrate total DOC and repartition its aqueous
%   fraction at every RK stage under the COMPASS Kd equilibrium closure.
% Change: 2026-09-14 Codex -- transfer accepted MEND dormancy fluxes between
%   matching mobile/attached/retained phases before other phase scaling.
% Change: 2026-09-14 Jizixing Zhu -- reject dynamic DOM at the equilibrium
%   RK4 boundary instead of silently overwriting source phase semantics.
% Change: 2026-10-27 Wenjie Liu -- expose the complete diagnostic flux schema
%   with zero placeholders for processes absent from equilibrium RK4.
    assert((~isfield(domParameters,'sorption_model') ...
        ||strcmp(domParameters.sorption_model,'linear_equilibrium')) ...
        &&(~isfield(domAccepted,'sorption_model') ...
        ||strcmp(domAccepted.sorption_model,'linear_equilibrium')) ...
        &&~isfield(domAccepted,'doc_sorbed_g_cm3') ...
        &&~isfield(domAccepted,'don_sorbed_g_cm3'), ...
        'soil:microbe:dynamicDomRequiresSourceCoupling', ...
        'Dynamic DOM requires the paired source-hour coupling adapter, not equilibrium RK4.');
    bdBottom=bulkDensityBottomUp(:); thetaBottom=thetaBottomUp(:);
    n=numel(bdBottom);
    assert(numel(thetaBottom)==n&&all(isfinite([bdBottom;thetaBottom])) ...
        &&all(bdBottom>0)&&all(thetaBottom>=0), ...
        'soil:microbe:driverShape','MEND theta and bulk density must align.');
    defaulted=false;
    if nargin<8||isempty(temperatureBottomUp)
        temperatureBottomUp=repmat(20,n,1); defaulted=true;
    end
    if nargin<9||isempty(headBottomUp)
        headBottomUp=repmat(-0.033/9.80665e-5,n,1); defaulted=true;
    end
    temperatureBottomUp=temperatureBottomUp(:); headBottomUp=headBottomUp(:);
    assert(numel(temperatureBottomUp)==n&&numel(headBottomUp)==n ...
        &&all(isfinite([temperatureBottomUp;headBottomUp])), ...
        'soil:microbe:mendEnvironment','MEND T/head profiles are invalid.');
    assert(isfield(accepted,'schema_version')&&accepted.schema_version>=3, ...
        'soil:microbe:mendSchema', ...
        'MEND production requires soil_microbes schema version 3 or newer.');
    validateAggregates(accepted);
    assert(isfield(domAccepted,'doc_total_g_cm3') ...
        &&numel(domAccepted.doc_total_g_cm3)==n ...
        &&all(isfinite(domAccepted.doc_total_g_cm3(:))) ...
        &&all(domAccepted.doc_total_g_cm3(:)>=0), ...
        'soil:microbe:mendDomState', ...
        'MEND requires finite nonnegative total DOC on each layer.');
    phase=soil.dom.equilibrate(ones(n,1),thetaBottom, ...
        bdBottom,domParameters.kd_doc_cm3_g);
    aqueousFraction=flip(thetaBottom.*phase.aqueous_g_cm3_water);
    assert(all(abs(accepted.QOC)<=1e-14), ...
        'soil:microbe:duplicateQocOwner', ...
        'Accepted QOC must be zero because soil_dom owns total DOM.');
    bdTop=flip(bdBottom);
    state=accepted;
    state.DOC=flip(domAccepted.doc_total_g_cm3(:))./bdTop.*1000;
    state=canonical(state);
    p=soil.microbe.mendEnvironmentalParameters(parameters, ...
        flip(temperatureBottomUp),flip(headBottomUp).*9.80665e-5, ...
        parameters.pH_top_down);
    fields={'POM1','POM2','MOC','MB','BD','DOC','ENZP1','ENZP2','EM','IC'};
    carbonBefore=entityCarbon(state);
    integratedFlux=emptyFlux(state.MB);
    substeps={}; elapsed=0; count=0; rejected=0; minStep=Inf;
    while elapsed<dt_s-1e-12
        trial=min(parameters.reaction_substep_s,dt_s-elapsed);
        while true
            try
                [next,stepFlux,stages]=rk4( ...
                    state,p,trial/3600,fields,aqueousFraction);
                valid=validState(next,fields);
            catch exception
                if ~strcmp(exception.identifier,'soil:microbe:mendStageNegative')
                    rethrow(exception);
                end
                valid=false;
            end
            if valid, break; end
            trial=trial/2; rejected=rejected+1;
            if trial<1e-6
                error('soil:microbe:nonnegativeStep', ...
                    'MEND RK4 could not find a nonnegative carbon substep.');
            end
        end
        substeps{end+1}=struct('carbon_before',state, ... %#ok<AGROW>
            'carbon_after',next,'stages',stages,'integrated_flux',stepFlux, ...
            'dt_s',trial);
        state=next;
        integratedFlux=addFlux(integratedFlux,stepFlux);
        elapsed=elapsed+trial; count=count+1; minStep=min(minStep,trial);
    end
    state=canonical(state);
    carbonResidual=entityCarbon(state)-carbonBefore;
    assert(max(abs(carbonResidual))<=1e-10.*max(1,max(abs(carbonBefore))), ...
        'soil:microbe:carbonClosure','MEND physical-pool carbon identity failed.');
    docTotalBottom=flip(state.DOC.*bdTop./1000);
    icIncrease=state.IC-accepted.IC;
    co2RateBottom=flip(icIncrease.*bdTop.*1000./dt_s);
    acceptedState=rmfield(state,'DOC');
    acceptedState.QOC=zeros(size(acceptedState.QOC));
    phaseLedger=cell(size(substeps));
    for i=1:numel(substeps)
        item=substeps{i}; f=item.integrated_flux;
        phaseLedger{i}=struct( ...
            'before',struct('MB',item.carbon_before.MB,'BD',item.carbon_before.BD), ...
            'after',struct('MB',item.carbon_after.MB,'BD',item.carbon_after.BD), ...
            'a_to_d',f.MBA_to_MBD,'d_to_a',f.MBD_to_MBA, ...
            'dormant_loss',f.CO2_maintenance_MBD);
    end
    acceptedState=soil.microbe.carryReactionPhases( ...
        acceptedState,phaseLedger,'C',bulkDensityBottomUp);
    candidate=struct('state',acceptedState,'doc_total_g_cm3',docTotalBottom, ...
        'co2_c_rate_ug_cm3_s',co2RateBottom, ...
        'integrated_flux_mg_g',integratedFlux, ...
        'substep_ledger',{substeps},'carbon_residual_mg_g',carbonResidual, ...
        'substep_count',count,'rejected_substeps',rejected, ...
        'min_substep_s',minStep,'rate_time_basis','hour', ...
        'source_model','MEND','source_commit',parameters.source_commit, ...
        'environment_defaulted',defaulted, ...
        'dynamic_qom_disabled_for_compass_dom_owner',true);
end

function [out,integrated,ledger]=rk4(state,p,h,fields,aqueousFraction)
% One MEND RK4 substep; total DOC is the integrated conserved coordinate.
% Editor: Jizixing Zhu; Date: 2026-09-12.
% Change: 2026-09-13 Codex -- derive aqueous DOC from stage total before fluxes.
    s1=aqueousState(state,aqueousFraction);
    [k1,f1]=soil.microbe.mendDerivatives(s1,p);
    total2=addScaled(state,k1,h/2,fields); assertStage(total2,fields);
    s2=aqueousState(total2,aqueousFraction);
    [k2,f2]=soil.microbe.mendDerivatives(s2,p);
    total3=addScaled(state,k2,h/2,fields); assertStage(total3,fields);
    s3=aqueousState(total3,aqueousFraction);
    [k3,f3]=soil.microbe.mendDerivatives(s3,p);
    total4=addScaled(state,k3,h,fields); assertStage(total4,fields);
    s4=aqueousState(total4,aqueousFraction);
    [k4,f4]=soil.microbe.mendDerivatives(s4,p);
    out=state;
    for i=1:numel(fields)
        name=fields{i};
        out.(name)=state.(name)+h/6.*(k1.(name)+2*k2.(name) ...
            +2*k3.(name)+k4.(name));
    end
    integrated=weightedFlux(f1,f2,f3,f4,h);
    ledger=struct('state',{s1,s2,s3,s4},'flux',{f1,f2,f3,f4}, ...
        'weight_h',{h/6,h/3,h/3,h/6});
end

function stage=aqueousState(totalState,aqueousFraction)
% Convert the sole total-DOC reaction coordinate to the MEND aqueous donor.
% Editor: Jizixing Zhu; Date: 2026-09-13.
    stage=totalState;
    stage.DOC=aqueousFraction.*totalState.DOC;
end

function out=addScaled(state,k,scale,fields)
    out=state;
    for i=1:numel(fields)
        name=fields{i}; out.(name)=state.(name)+scale.*k.(name);
    end
end

function assertStage(state,fields)
    if ~validState(state,fields)
        error('soil:microbe:mendStageNegative', ...
            'MEND RK stage contains a negative or nonfinite physical pool.');
    end
end

function yes=validState(state,fields)
    yes=true;
    for i=1:numel(fields)
        v=state.(fields{i});
        yes=yes&&all(isfinite(v))&&all(v>=-1e-12);
    end
end

function state=canonical(state)
% Rebuild compatibility aggregates from authoritative MEND subpools.
% Editor: Jizixing Zhu; Date: 2026-09-12.
    state.POC=state.POM1+state.POM2;
    state.EP=state.ENZP1+state.ENZP2;
    state.Tot=state.POC+state.MOC+state.MB+state.BD+state.EP+state.EM;
end

function validateAggregates(state)
% Reject stale or ambiguous MEND aggregate ownership.
% Editor: Jizixing Zhu; Date: 2026-09-12.
    required={'POM1','POM2','POC','MOC','MB','BD','ENZP1','ENZP2','EP','EM','IC'};
    for i=1:numel(required)
        assert(isfield(state,required{i})&&iscolumn(state.(required{i})) ...
            &&all(isfinite(state.(required{i})))&&all(state.(required{i})>=0), ...
            'soil:microbe:mendState','Invalid accepted MEND pool %s.',required{i});
    end
    assert(all(abs(state.POC-state.POM1-state.POM2)<=1e-11.*max(1,state.POC)) ...
        &&all(abs(state.EP-state.ENZP1-state.ENZP2)<=1e-11.*max(1,state.EP)), ...
        'soil:microbe:mendAggregate','Accepted MEND aggregate is stale.');
end

function value=entityCarbon(state)
    value=state.POM1+state.POM2+state.MOC+state.MB+state.BD ...
        +state.DOC+state.ENZP1+state.ENZP2+state.EM+state.IC;
end

function out=weightedFlux(a,b,c,d,h)
    names=fieldnames(a); out=emptyFlux(a.POM1_decomp);
    for i=1:numel(names)
        name=names{i};
        if strcmp(name,'substrate_saturation'), continue; end
        out.(name)=h./6.*(a.(name)+2.*b.(name)+2.*c.(name)+d.(name));
    end
end

function out=emptyFlux(seed)
    names={'POM1_decomp','POM2_decomp','POM1_to_DOM','POM2_to_DOM', ...
        'POM1_to_MOM','POM2_to_MOM','MOM_to_DOM','DOM_to_MBA', ...
        'CO2_growth','CO2_maintenance_MBA','CO2_maintenance_MBD', ...
        'MBA_mortality','MBA_to_DOM','MBA_to_POM1','MBA_to_POM2', ...
        'MBA_to_ENZP1','MBA_to_ENZP2','MBA_to_ENZM', ...
        'ENZP1_to_DOM','ENZP2_to_DOM','ENZM_to_DOM', ...
        'MBA_to_MBD','MBD_to_MBA', ...
        'DOM_to_QOM','QOM_to_DOM','DOM_sorption_net', ...
        'ENZN_production','ENZN_to_MBA','N_overflow_CO2'};
    out=struct();
    for i=1:numel(names), out.(names{i})=zeros(size(seed)); end
end

function out=addFlux(out,value)
    names=fieldnames(out);
    for i=1:numel(names), out.(names{i})=out.(names{i})+value.(names{i}); end
end
