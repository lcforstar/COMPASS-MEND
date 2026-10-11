function [model,audit]=materializeMeasuredBiogeochemistryInitialState(fresh,donor,options)
%MATERIALIZEMEASUREDBIOGEOCHEMISTRYINITIALSTATE Condition a C/N prior at new t0.
%   [MODEL,AUDIT]=...(FRESH,DONOR) accepts an initialized, not-yet-prepared
%   production model and a same-parameter donor model. It scales independent
%   soil organic C/N pools by layer to FRESH's measured SOC target, retaining
%   FRESH IOM, DIN, water/heat and vegetation. Unreacted donor AOM is copied
%   separately, not counted in SOC. No integration, I/O or physical flux is
%   performed. Initialization adjustments are returned, never hidden in flux.
%
%   Supported contract: single-column source_hourly MEND, linear-equilibrium
%   total DOC/DON, explicit MB/BD/ENZN phases. Microbial C is top-down mg/g;
%   organic N/DOM are bottom-up g/cm3. AUDIT profiles are bottom-up g/cm3,
%   inventories and initialization adjustments g C/N per cm2 ground.
%   By default, measured total DOC/DON from FRESH remain bitwise unchanged.
%   options.dom_policy='donor_preheated' instead restores the donor DOC/DON
%   totals while keeping the fresh DOM schema and zero cumulative ledgers.
%   Each layer:
%     lambda=(SOC_target-IOM_fresh-DOC_selected) /
%            (SOC_donor-IOM_donor-DOC_donor).
%   Donor slow/biomass C and their paired N inventories share lambda. Alias
%   pools are rebuilt; cumulative respiration/exports and tracer are fresh.
%
%   donor_kind: periodic_state (default, complete same-calendar-phase window),
%   initial_guess (unstepped), or diagnostic_state (explicit nonperiodic test).
%   This adapter NEVER certifies spinup convergence or scientific release.
%   It preserves measured DIN using the accepted FE dry-mass observation
%   operator; the historical scorer discrepancy is retained as an audit only.
%   Errors reject incompatible identity/schema, invalid pools, nonfresh input,
%   wrong phase or infeasible SOC. No clipping, DIN borrowing or DOM re-reserve.
%   Sources: init.buildSoilCarbon/reserveInitialDOM/buildSoilMicrobes;
%   soil.microbe.advanceSourceCandidate; diagnostics.writeStep SOC operator.
%   Editor: Wenjie Liu; Date: 2026-09-29.
%   Change: 2026-09-29 Wenjie Liu -- verify corrected FE measured DIN mapping.
%   Change: 2026-10-10 Wenjie Liu -- preserve measured fresh DOC/DON at the
%     production t0 and condition only donor slow/biomass pools to measured SOC.
%   Change: 2026-10-11 Wenjie Liu -- add an explicit preheated-DOM policy so a
%     complete MEND prior can be conditioned to measured SOC without replacing
%     measured water or DIN and without carrying historical flux ledgers.
    arguments
        fresh (1,1) struct
        donor (1,1) struct
        options.donor_kind (1,:) char {mustBeMember(options.donor_kind, ...
            {'periodic_state','initial_guess','diagnostic_state'})}='periodic_state'
        options.donor_label (1,:) char='provided_prior'
        options.dom_policy (1,:) char {mustBeMember(options.dom_policy, ...
            {'fresh_measured','donor_preheated'})}='fresh_measured'
    end
    prefix='materializeMeasuredBiogeochemistryInitialState';
    assert(isfield(fresh,'status')&&string(fresh.status)=="initialized" ...
        &&~isfield(fresh,'runtime'),[prefix ':notFresh'], ...
        'Use an initialized model before prepareRuntime/beginRun, not a resume.');
    assert(isfield(donor,'built'),[prefix ':donor'],'Donor must contain built.');
    b=fresh.built; d=donor.built; y0=b.Y;
    assertCompatible(b,d,prefix);
    phase=donorPhase(fresh,donor,options.donor_kind,prefix);
    before=inventory(d,prefix); freshBefore=inventory(b,prefix);
    constraints=measuredConstraints(fresh,prefix);
    assertClose(freshBefore.soc_g_cm3,constraints.soc_g_cm3, ...
        [prefix ':freshSOC'],'Fresh SOC does not match its input mapping.');
    requireZero(y0.soil_organic.cumulative,[prefix ':dirtyFresh']);
    assert(all(y0.soil_microbes.IC==0),[prefix ':dirtyFresh'], ...
        'Fresh cumulative respiratory carbon must be zero.');
    for owner={'soil_dom','soil_microbes'}
        for field={'bottom_export_cumulative','boundary_exchange_cumulative'}
            if isfield(y0.(owner{1}),field{1})
                requireZero(y0.(owner{1}).(field{1}),[prefix ':dirtyFresh']);
            end
        end
    end

    if strcmp(options.dom_policy,'donor_preheated')
        selectedDOC=d.Y.soil_dom.doc_total_g_cm3;
        selectedDON=d.Y.soil_dom.don_total_g_cm3;
    else
        selectedDOC=y0.soil_dom.doc_total_g_cm3;
        selectedDON=y0.soil_dom.don_total_g_cm3;
    end
    donorDOC=d.Y.soil_dom.doc_total_g_cm3;
    available=constraints.soc_g_cm3-y0.soil_organic.IOM_C-selectedDOC;
    donorResponsive=before.soc_g_cm3-d.Y.soil_organic.IOM_C-donorDOC;
    assert(all(available>=0)&&all(donorResponsive>0), ...
        [prefix ':infeasibleSOC'],'SOC target or donor responsive C is infeasible.');
    scaleBottom=available./donorResponsive; scaleTop=flip(scaleBottom);
    assert(all(isfinite(scaleBottom)),[prefix ':infeasibleSOC'],'Invalid SOC scale.');
    y=y0; cfields=carbonFields(); nfields=nitrogenFields();
    for j=1:numel(cfields)
        y.soil_microbes.(cfields{j})=d.Y.soil_microbes.(cfields{j}).*scaleTop;
        y.soil_organic.(nfields{j})=d.Y.soil_organic.(nfields{j}).*scaleBottom;
    end
    for group={'MB','BD','ENZN'}
        for statePhase={'mobile','attached','retained'}
            key=[group{1} '_' statePhase{1}];
            y.soil_microbes.(key)=d.Y.soil_microbes.(key).*scaleTop;
            y.soil_organic.([key '_N'])=d.Y.soil_organic.([key '_N']).*scaleBottom;
        end
    end
    % Keep the current DOM schema and fresh zero ledgers. The selected total
    % inventories are either measured fresh values or the preheated donor prior.
    y.soil_dom=y0.soil_dom;
    y.soil_dom.doc_total_g_cm3=selectedDOC;
    y.soil_dom.don_total_g_cm3=selectedDON;
    y.soil_organic.NEC_C=d.Y.soil_organic.NEC_C.*scaleBottom;
    y.soil_organic.NEC_N=d.Y.soil_organic.NEC_N.*scaleBottom;
    y.soil_organic.aom=d.Y.soil_organic.aom;
    y.soil_microbes.POC=y.soil_microbes.POM1+y.soil_microbes.POM2;
    y.soil_microbes.EP=y.soil_microbes.ENZP1+y.soil_microbes.ENZP2;
    y.soil_microbes.Tot=zeros(b.domain.NL,1);
    for j=1:numel(cfields)
        y.soil_microbes.Tot=y.soil_microbes.Tot+y.soil_microbes.(cfields{j});
    end
    y.soil_organic.PON=y.soil_organic.PON1+y.soil_organic.PON2;
    y.soil_organic.EP_N=y.soil_organic.ENZP1_N+y.soil_organic.ENZP2_N;
    if b.control.switches.lFertilizerSourceTracer==1
        y.soil_fertilizer_n_tracer=soil.nitrogen.tracer.emptyState(b);
    end
    b.Y=y;
    after=inventory(b,prefix);
    assertClose(after.soc_g_cm3,constraints.soc_g_cm3, ...
        [prefix ':SOCClosure'],'Conditioned SOC must match the measured target.');
    if strcmp(options.dom_policy,'fresh_measured')
        assert(isequaln(y.soil_dom,y0.soil_dom), ...
            [prefix ':changedMeasuredDOM'],'Measured fresh DOC/DON must be preserved bitwise.');
    else
        assert(isequaln(y.soil_dom.doc_total_g_cm3,d.Y.soil_dom.doc_total_g_cm3) ...
            &&isequaln(y.soil_dom.don_total_g_cm3,d.Y.soil_dom.don_total_g_cm3), ...
            [prefix ':changedPreheatedDOM'],'Preheated donor DOC/DON must be preserved bitwise.');
    end
    allowed={'soil_microbes','soil_organic','soil_dom','soil_fertilizer_n_tracer'};
    for key=fieldnames(y0)'
        if ~ismember(key{1},allowed)
            assert(isequaln(y.(key{1}),y0.(key{1})),[prefix ':changedProtectedState'], ...
                'Protected fresh state %s changed.',key{1});
        end
    end
    assert(isequal(y.soil_organic.IOM_C,y0.soil_organic.IOM_C), ...
        [prefix ':changedIOM'],'IOM is a fixed fresh prior, not a fitted pool.');
    % Rebuild disposable caches from the final state, preserving all accepted Y.
    b=init.rebuildSoilAuxiliaryFromAccepted(b,0);
    assert(isequaln(b.Y,y),[prefix ':rebuildChangedState'],'Aux rebuild changed Y.');
    theta=mean(b.workspace.derived.Theta_LL,2);
    for kind={'doc','don'}
        p=soil.dom.phaseState(y.soil_dom,kind{1},theta, ...
            b.parameters.soil.bulk_density,b.parameters.soil_dom);
        assertClose(p.total_g_cm3_bulk,p.free_g_cm3_bulk+p.sorbed_g_cm3_bulk, ...
            [prefix ':DOMPhaseClosure'],'DOM phase reconstruction failed.');
    end
    adjustmentC=after.total_C_g_cm2-before.total_C_g_cm2;
    adjustmentN=after.total_N_g_cm2-before.total_N_g_cm2;
    audit=struct('schema_version',1,'method','measured_SOC_conditioned_organic_prior', ...
        'donor_kind',options.donor_kind,'donor_label',options.donor_label, ...
        'donor_phase',phase,'grid_order','bottom_up', ...
        'profile_units','g C_or_N/cm3_bulk_soil','inventory_units','g C_or_N/cm2_ground', ...
        'lambda_bottom_up',scaleBottom,'fresh_before',freshBefore, ...
        'donor_before',before,'conditioned_after',after, ...
        'initialization_adjustment_C_g_cm2',adjustmentC, ...
        'initialization_adjustment_N_g_cm2',adjustmentN, ...
        'organic_C_adjustment_profile_g_cm3',after.organic_C_g_cm3-before.organic_C_g_cm3, ...
        'organic_N_adjustment_profile_g_cm3',after.organic_N_g_cm3-before.organic_N_g_cm3, ...
        'physical_flux_added',false,'DIN_fresh_bitwise_preserved',true, ...
        'DOC_DON_fresh_bitwise_preserved',strcmp(options.dom_policy,'fresh_measured'), ...
        'DOC_DON_donor_bitwise_preserved',strcmp(options.dom_policy,'donor_preheated'), ...
        'DOM_initialization_policy',options.dom_policy, ...
        'water_heat_vegetation_fresh_bitwise_preserved',true, ...
        'tracer_policy','fresh_period_mineral_fertilizer_only_background_unlabelled', ...
        'parameter_grid_boundary_identity_verified',true, ...
        'campaign_source_forcing_hash_verified',false, ...
        'spinup_convergence_verified',false,'scientific_release',false, ...
        'measured_input_roundtrip',constraints.roundtrip, ...
        'legacy_scorer_t0_difference',constraints.legacy_scorer);
    model=fresh; model.built=b;
    model.built.initialization.measured_biogeochemistry=audit;
end

function assertCompatible(b,d,prefix)
    assert(isequaln(b.domain,d.domain)&&isequaln(b.boundary,d.boundary) ...
        &&isequaln(b.control.switches,d.control.switches), ...
        [prefix ':identity'],'Grid, boundary or active switches differ.');
    p=b.parameters; q=d.parameters;
    % This is historical cold-start metadata, not reaction/transport parameters.
    p.soil_organic=rmfield(p.soil_organic,'initial_carbon_top_down');
    q.soil_organic=rmfield(q.soil_organic,'initial_carbon_top_down');
    assert(isequaln(p,q),[prefix ':identity'], ...
        'Donor parameters differ; another parameter vector is not this prior.');
    for built={b,d}
        a=built{1};
        assert(a.parameters.soil_microbes.reaction_method=="source_hourly" ...
            &&a.parameters.soil_dom.sorption_model=="linear_equilibrium" ...
            &&a.Y.soil_microbes.schema_version==4 ...
            &&a.Y.soil_organic.schema_version==5 ...
            &&a.Y.soil_dom.schema_version==1, ...
            [prefix ':unsupportedSchema'],'Only the current source-hour linear DOM schema is supported.');
    end
end

function a=donorPhase(fresh,donor,kind,prefix)
    steps=0;
    if isfield(donor,'runtime')
        steps=donor.runtime.run_control.last_completed_kt;
    end
    assert(isscalar(steps)&&isfinite(steps)&&steps>=0&&steps==fix(steps), ...
        [prefix ':donorPhase'],'Invalid accepted donor step.');
    sourceTime=donor.built.control.time.start_time ...
        +seconds(steps*donor.built.control.time.Delt_t);
    targetTime=fresh.built.control.time.start_time;
    samePhase=month(sourceTime)==month(targetTime)&&day(sourceTime)==day(targetTime) ...
        &&timeofday(sourceTime)==timeofday(targetTime);
    if strcmp(kind,'periodic_state')
        assert(isfield(donor,'runtime')&&steps==donor.built.control.time.tS ...
            &&samePhase,[prefix ':donorPhase'], ...
            'Periodic donor must complete its configured window at the target calendar phase.');
    elseif strcmp(kind,'initial_guess')
        assert(~isfield(donor,'runtime')&&steps==0&&samePhase, ...
            [prefix ':donorPhase'],'Initial-guess donor must be unstepped at the target phase.');
    end
    a=struct('accepted_steps',steps,'source_time',char(sourceTime), ...
        'target_time',char(targetTime),'same_calendar_phase',samePhase);
end

function a=inventory(b,prefix)
    y=b.Y; n=b.domain.NL; rho=b.parameters.soil.bulk_density(:); dz=b.domain.dz(:);
    assert(~isfield(y.soil_dom,'doc_sorbed_g_cm3') ...
        &&~isfield(y.soil_dom,'don_sorbed_g_cm3'), ...
        [prefix ':unsupportedSchema'],'Linear DOM cannot carry dynamic sorbed owners.');
    for key={'nh4','no3','urea','C_NH4_granule','C_NH4_ice','C_NO3_ice'}
        v=y.(key{1});
        assert(all(isfinite(v(:)))&&all(v(:)>=0), ...
            [prefix ':invalidPool'],'DIN inventory must be finite and nonnegative.');
    end
    cfields=carbonFields(); nfields=nitrogenFields(); c=zeros(n,1); nn=c;
    assert(y.soil_microbes.grid_order=="top_down" ...
        &&y.soil_organic.grid_order=="bottom_up",[prefix ':gridOrder'],'Unsupported state orientation.');
    for j=1:numel(cfields)
        v=pool(y.soil_microbes,cfields{j},n,prefix);
        c=c+flip(v).*rho/1000;
        nn=nn+pool(y.soil_organic,nfields{j},n,prefix);
    end
    assert(all(pool(y.soil_microbes,'QOC',n,prefix)==0) ...
        &&all(pool(y.soil_organic,'QON',n,prefix)==0), ...
        [prefix ':duplicateDOM'],'QOC/QON must not duplicate the total DOM owner.');
    assertClose(y.soil_microbes.POC,y.soil_microbes.POM1+y.soil_microbes.POM2, ...
        [prefix ':aggregate'],'POC mismatch.');
    assertClose(y.soil_microbes.EP,y.soil_microbes.ENZP1+y.soil_microbes.ENZP2, ...
        [prefix ':aggregate'],'EP mismatch.');
    assertClose(y.soil_organic.PON,y.soil_organic.PON1+y.soil_organic.PON2, ...
        [prefix ':aggregate'],'PON mismatch.');
    assertClose(y.soil_organic.EP_N,y.soil_organic.ENZP1_N+y.soil_organic.ENZP2_N, ...
        [prefix ':aggregate'],'EP_N mismatch.');
    for group={'MB','BD','ENZN'}
        cc=zeros(n,1); pn=cc;
        for phase={'mobile','attached','retained'}
            key=[group{1} '_' phase{1}];
            cc=cc+pool(y.soil_microbes,key,n,prefix);
            pn=pn+pool(y.soil_organic,[key '_N'],n,prefix);
        end
        assertClose(cc,y.soil_microbes.(group{1}),[prefix ':phaseAggregate'],'C phases do not sum to total.');
        assertClose(pn,y.soil_organic.([group{1} '_N']),[prefix ':phaseAggregate'],'N phases do not sum to total.');
    end
    c=c+pool(y.soil_dom,'doc_total_g_cm3',n,prefix) ...
        +pool(y.soil_organic,'IOM_C',n,prefix)+pool(y.soil_organic,'NEC_C',n,prefix);
    nn=nn+pool(y.soil_dom,'don_total_g_cm3',n,prefix)+pool(y.soil_organic,'NEC_N',n,prefix);
    soc=c; aomC=zeros(n,1); aomN=aomC; surfaceC=0; surfaceN=0;
    assert(isequal(sort(fieldnames(y.soil_organic.aom)),sort({'litter';'straw';'manure'})), ...
        [prefix ':aomSchema'],'Unknown AOM owner; explicit inventory support is required.');
    for material={'litter','straw','manure'}
        s=y.soil_organic.aom.(material{1});
        if isfield(s,'components'),s=soil.litter.componentTotals(s,true);end
        aomC=aomC+pool(s,'element_C',n,prefix); aomN=aomN+pool(s,'element_N',n,prefix);
        surfaceC=surfaceC+pool(s,'surface_C',1,prefix);
        surfaceN=surfaceN+pool(s,'surface_N',1,prefix);
    end
    [nh4,no3]=soil.nitrogen.transport.calculateStorage(y.Theta_L,y.nh4,y.no3, ...
        b.parameters.soil_n.RHOKG,dz,n,y.Theta_I,y.C_NH4_ice,y.C_NO3_ice,y.C_NH4_granule);
    urea=sum(dz.*(y.urea(1:end-1)+y.urea(2:end))/2);
    a=struct('soc_g_cm3',soc,'organic_C_g_cm3',c+aomC, ...
        'organic_N_g_cm3',nn+aomN,'surface_AOM_C_g_cm2',surfaceC, ...
        'surface_AOM_N_g_cm2',surfaceN, ...
        'organic_C_g_cm2',sum((c+aomC).*dz)+surfaceC, ...
        'organic_N_g_cm2',sum((nn+aomN).*dz)+surfaceN, ...
        'DIN_plus_urea_g_cm2',(nh4+no3+urea)*1e-6);
    a.total_C_g_cm2=a.organic_C_g_cm2;
    a.total_N_g_cm2=a.organic_N_g_cm2+a.DIN_plus_urea_g_cm2;
    % Deliberately soil + surface AOM only: plant C/N is not this control volume.
end

function a=measuredConstraints(model,prefix)
    b=model.built; y=b.Y; raw=model.raw.initial; dz=b.domain.dz(:);
    assert(string(model.raw.config.switches.nitrogen_input_phase)=="total", ...
        [prefix ':inputPhase'],'Measured initial DIN must use the total dry-soil basis.');
    assert(all(y.Theta_I(:)==0),[prefix ':frozenInitialState'], ...
        'This measured-DIN roundtrip is scoped to the unfrozen Changwu t0.');
    z=(b.domain.z_node(1:end-1)+b.domain.z_node(2:end))/2;
    rho=b.parameters.soil.bulk_density(:);
    a.soc_g_cm3=init.interpLinearWithExtrapolation(raw.depth_cm,raw.soc_total_mgkg,z).*rho/1e6;
    expectedNH4=init.interpLinearWithExtrapolation(raw.depth_cm,raw.nh4,b.domain.z_node);
    expectedNO3=init.interpLinearWithExtrapolation(raw.depth_cm,raw.no3,b.domain.z_node);
    [actualNH4,actualNO3]=soil.nitrogen.observationProfile(y,b.parameters.soil_n,dz);
    assertClose(actualNH4,expectedNH4,[prefix ':freshDIN'],'Fresh NH4 differs from its measured input mapping.');
    assertClose(actualNO3,expectedNO3,[prefix ':freshDIN'],'Fresh NO3 differs from its measured input mapping.');
    a.roundtrip=struct('basis','accepted_FE_capacity_over_dry_mass', ...
        'nh4_max_abs_mgkg',max(abs(actualNH4-expectedNH4)), ...
        'no3_max_abs_mgkg',max(abs(actualNO3-expectedNO3)));
    % Record, do not silently repair, the independently discovered legacy
    % scorer mismatch (different node rho/RHOKG mapping and accepted theta).
    theta=[y.Theta_L(1,1); ...
        (dz(1:end-1).*y.Theta_L(1:end-1,2)+dz(2:end).*y.Theta_L(2:end,1)) ...
        ./(dz(1:end-1)+dz(2:end));y.Theta_L(end,2)];
    rhoN=[b.parameters.soil_n.RHO_bulk;b.parameters.soil_n.RHO_bulk(end)];
    buffer=[b.parameters.soil_n.RHOKG;b.parameters.soil_n.RHOKG(end)];
    a.legacy_scorer=struct('known_unresolved',false,'historical_formula_only',true, ...
        'nh4_max_abs_mgkg',max(abs(y.nh4.*(theta+buffer)./rhoN-expectedNH4)), ...
        'no3_max_abs_mgkg',max(abs(y.no3.*theta./rhoN-expectedNO3)));
end

function names=carbonFields()
    names={'POM1','POM2','MOC','MB','BD','ENZP1','ENZP2','EM','ENZN'};
end
function names=nitrogenFields()
    names={'PON1','PON2','MON','MB_N','BD_N','ENZP1_N','ENZP2_N','EM_N','ENZN_N'};
end
function v=pool(owner,key,n,prefix)
    assert(isfield(owner,key),[prefix ':missingPool'],'Missing pool %s.',key);
    v=owner.(key);
    assert(isnumeric(v)&&isequal(size(v),[n,1])&&all(isfinite(v))&&all(v>=0), ...
        [prefix ':invalidPool'],'Pool %s must be finite nonnegative %d-by-1.',key,n);
end
function assertClose(a,b,id,message)
    assert(isequal(size(a),size(b))&&all(isfinite(a(:)))&&all(isfinite(b(:))) ...
        &&all(abs(a(:)-b(:))<=1e-12*max(1,max(abs(b(:))))),id,message);
end
function requireZero(s,id)
    for field=fieldnames(s)'
        key=field{1}; if strcmp(key,'schema_version'),continue;end
        v=s.(key);
        if isstruct(v),requireZero(v,id);
        elseif isnumeric(v),assert(all(v(:)==0),id,'Fresh ledger %s must be zero.',key);
        end
    end
end
