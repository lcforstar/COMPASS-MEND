function [candidate, diagnostic, candidateApplied, microbeCandidate, domCandidate] = stageStep( ...
        accepted, schedule, applied, stepStart, dtSeconds, zplus, bulkDensityTop, ...
        microbeAccepted, domAccepted, microbeParameters, rootContext)
%STAGESTEP Stage manure/straw/root-residue C/N in AOM and TAN in DIN.
%
%   Event masses [g C/N cm^-2] retain their material identity. TAN is removed
%   before organic-N partition and returned as a top-down layer NH4 ledger.
%   Incorporated events use uniform-depth exact-overlap weights, independent
%   of bulk density. Surface additions remain surface AOM scalars [g/cm2];
%   incorporated additions enter bottom-up element AOM [g/cm3 bulk soil].
%   Microbes are unchanged. With optional manure_fast_input_enabled=false,
%   DOM is unchanged and AOM decomposition/leaching owns all organic release.
%   The opt-in incorporated-manure route sends manure_input_fraction(3) of
%   the SAME event C and organic N directly to total DOM; the remaining C/N
%   enters manure AOM. TAN is separate and is not split a second time.
%   Material-resolved added-density diagnostics are bottom-up NL-by-1
%   [g element cm^-3 bulk soil per step]; surface additions are [g cm^-2].
%   All returned states and APPLIED flags are detached candidates, committed
%   only after the complete outer transaction accepts. Non-unit legacy mixing
%   efficiencies fail explicitly.
%
% Editor: Jizixing Zhu; Date: 2026-08-28.
% Change: 2026-08-30 Codex -- stage manure TAN separately from organic AOM N.
% Change: 2026-10-10 Wenjie Liu -- direct manure C/N through the MEND input
%   adapter; uniform-depth placement for all organic events, no surface remainder.
% Editor: Wenjie Liu; Date: 2026-10-17.
% Change: 2026-10-17 Wenjie Liu -- restore manure and straw AOM buffering;
%   verify actual AOM increments plus TAN against event C/N at 1e-14.
% Change: 2026-10-30 Wenjie Liu -- stage external root C/N into existing
%   litter components using frozen bottom-up normalized Jackson weights.
%   Optional rootContext has weights_bottom_up [NL,1] and fractions [3,1]
%   (metabolic/cellulose/lignin). It is required only for root events.
% Change: 2026-10-30 Wenjie Liu -- opt-in paired manure fast DOC/DON input;
%   preserve the legacy arithmetic when disabled and expose disjoint source
%   additions for conservative passive manure-N tracing.

    candidate=accepted;
    if nargin<8, microbeAccepted=[]; end
    if nargin<9, domAccepted=[]; end
    if nargin<10, microbeParameters=[]; end
    if nargin<11, rootContext=[]; end
    microbeCandidate=microbeAccepted;
    domCandidate=domAccepted;
    candidateApplied=logical(applied(:).');
    diagnostic=emptyDiagnostic(numel(zplus));
    if schedule.n_events==0, return; end
    if isfield(schedule,'mixing_efficiency')
        efficiency=schedule.mixing_efficiency;
        assert(all(isfinite(efficiency(:)))&&all(efficiency(:)==1), ...
            'stageOrganicFertilizer:mixingEfficiency', ...
            ['Legacy partial mixing is not supported: mixing_efficiency must ', ...
             'equal 1. Use explicit surface or incorporated placement.']);
    end
    assert(isfield(candidate,'aom')&&candidate.schema_version>=2, ...
        'stageOrganicFertilizer:legacyState', ...
        'Organic management requires the native COMPASS AOM state.');
    if numel(candidateApplied)~=schedule.n_events
        error('stageOrganicFertilizer:flagSize', ...
            'Applied-flag count must match organic management events.');
    end
    stepEnd=stepStart+seconds(dtSeconds);
    events=find(~candidateApplied & schedule.time>=stepStart ...
        & schedule.time<stepEnd);
    if isempty(events), return; end

    dzTop=-diff([0;zplus(:)]);
    if any(~isfinite(dzTop))||any(dzTop<=0)
        error('stageOrganicFertilizer:badGrid', ...
            'zplus must define positive top-down COMPASS element thicknesses.');
    end
    if nargin<7||isempty(bulkDensityTop), bulkDensityTop=ones(size(dzTop)); end
    bulkDensityTop=bulkDensityTop(:);
    assert(numel(bulkDensityTop)==numel(dzTop) ...
        &&all(isfinite(bulkDensityTop))&&all(bulkDensityTop>0), ...
        'stageOrganicFertilizer:bulkDensity', ...
        'Top-down positive bulk density is required for MEND input conversion.');
    addedC=0; addedN=0; addedOrganicN=0; addedTanN=0;
    materials=cell(1,numel(events));
    materialC=struct('straw',0,'manure',0,'root_residue',0);
    materialN=struct('straw',0,'manure',0,'root_residue',0);
    layerC=zeros(numel(dzTop),1); layerN=zeros(numel(dzTop),1);
    layerTanN=zeros(numel(dzTop),1);
    directC=zeros(1,3); directN=zeros(1,3);
    nLayers=numel(dzTop);
    manureFastC=zeros(nLayers,1);manureFastN=zeros(nLayers,1);
    manureAomC=zeros(nLayers,1);manureAomN=zeros(nLayers,1);
    manureTanN=zeros(nLayers,1);surfaceManureC=0;surfaceManureN=0;
    fastEnabled=fastManureEnabled(microbeParameters);
    [beforeC,beforeN]=aomTotals(accepted.aom,flip(dzTop));
    for q=1:numel(events)
        event=events(q);
        material=eventMaterial(schedule,event);
        isRoot=strcmp(material,'root_residue');
        if ~isRoot&&(~ismember(material,{'straw','manure'})||~isfield(candidate.aom,material))
            error('stageOrganicFertilizer:material', ...
                'Unsupported organic material type "%s".',material);
        end
        carbon=schedule.carbon_g_cm2(event);
        nitrogen=schedule.nitrogen_g_cm2(event);
        organicNitrogen=nitrogen; tanNitrogen=0;
        if isfield(schedule,'organic_nitrogen_g_cm2')
            organicNitrogen=schedule.organic_nitrogen_g_cm2(event);
        end
        if isfield(schedule,'tan_nitrogen_g_cm2')
            tanNitrogen=schedule.tan_nitrogen_g_cm2(event);
        end
        assert(all(isfinite([carbon,nitrogen,organicNitrogen,tanNitrogen])) ...
            &&all([carbon,nitrogen,organicNitrogen,tanNitrogen]>=0), ...
            'stageOrganicFertilizer:negativeInput','Event C/N must be finite and nonnegative.');
        assert(abs(organicNitrogen+tanNitrogen-nitrogen)<= ...
            1e-14*max(1,nitrogen), 'stageOrganicFertilizer:tanClosure', ...
            'Organic N plus TAN must equal event total N.');
        placement=lower(string(schedule.placement{event}));
        assert(isscalar(placement)&&((isRoot&&placement=="root_profile") ...
            ||(~isRoot&&ismember(placement,["surface","incorporated"]))), ...
            'stageOrganicFertilizer:placement','Invalid organic-event placement.');
        isManure=strcmp(material,'manure');
        if fastEnabled&&isManure
            assert(placement=="incorporated", ...
                'stageOrganicFertilizer:fastManurePlacement', ...
                'Fast manure input supports incorporated manure only.');
        end
        if isRoot
            assert(tanNitrogen==0&&organicNitrogen==nitrogen, ...
                'stageOrganicFertilizer:rootTan','Root residue must contain organic N only, with zero TAN.');
            [weightsBottom,fractions]=rootProfile(rootContext,numel(dzTop));
            weights=flip(weightsBottom);
        elseif placement=="surface"
            weights=zeros(size(dzTop)); weights(1)=1;
        else
            weights=uniformDepthWeights(dzTop, ...
                schedule.depth_start_cm(event),schedule.depth_end_cm(event));
        end
        cTop=carbon.*weights./dzTop;
        nTop=organicNitrogen.*weights./dzTop;
        slowCTop=cTop;slowNTop=nTop;
        if fastEnabled&&isManure
            fastFraction=manureFastFraction(microbeParameters);
            if fastFraction>0
                validateDomCandidate(domCandidate,nLayers);
                fastC=fastFraction.*carbon;fastN=fastFraction.*organicNitrogen;
                fastCTop=fastC.*weights./dzTop;
                fastNTop=fastN.*weights./dzTop;
                % Remaining-first splitting preserves the supplied paired
                % event mass; no mineral N is created or borrowed.
                slowCTop=(carbon-fastC).*weights./dzTop;
                slowNTop=(organicNitrogen-fastN).*weights./dzTop;
                domCandidate.doc_total_g_cm3=domCandidate.doc_total_g_cm3+flip(fastCTop);
                domCandidate.don_total_g_cm3=domCandidate.don_total_g_cm3+flip(fastNTop);
                manureFastC=manureFastC+flip(fastCTop);
                manureFastN=manureFastN+flip(fastNTop);
                directC(3)=directC(3)+fastC;directN(3)=directN(3)+fastN;
            end
        end
        if isRoot
            assert(isfield(candidate.aom,'litter') ...
                &&isfield(candidate.aom.litter,'components'), ...
                'stageOrganicFertilizer:rootComponents','Root residue requires component litter AOM.');
            litter=soil.litter.componentTotals(candidate.aom.litter,true);
            names={'metabolic','cellulose','lignin'};
            for j=1:3
                name=names{j};
                litter.components.(name).element_C=litter.components.(name).element_C ...
                    +fractions(j).*flip(cTop);
                litter.components.(name).element_N=litter.components.(name).element_N ...
                    +fractions(j).*flip(nTop);
            end
            candidate.aom.litter=soil.litter.componentTotals(litter,false);
        elseif placement=="surface"
            candidate.aom.(material).surface_C=candidate.aom.(material).surface_C+carbon;
            candidate.aom.(material).surface_N=candidate.aom.(material).surface_N+organicNitrogen;
        else
            candidate.aom.(material).element_C=candidate.aom.(material).element_C+flip(slowCTop);
            candidate.aom.(material).element_N=candidate.aom.(material).element_N+flip(slowNTop);
        end
        if isManure
            if placement=="surface"
                surfaceManureC=surfaceManureC+carbon;
                surfaceManureN=surfaceManureN+organicNitrogen;
            else
                manureAomC=manureAomC+flip(slowCTop);
                manureAomN=manureAomN+flip(slowNTop);
            end
            manureTanN=manureTanN+flip(tanNitrogen.*weights./dzTop);
        end
        if isRoot||placement=="incorporated"
            layerC=layerC+cTop.*dzTop;
            layerN=layerN+nTop.*dzTop;
        end
        layerTanN=layerTanN+tanNitrogen.*weights;
        candidateApplied(event)=true;
        addedC=addedC+carbon; addedN=addedN+nitrogen;
        addedOrganicN=addedOrganicN+organicNitrogen;
        addedTanN=addedTanN+tanNitrogen;
        materialC.(material)=materialC.(material)+carbon;
        materialN.(material)=materialN.(material)+nitrogen;
        materials{q}=material;
    end
    expectedC=sum(schedule.carbon_g_cm2(events));
    expectedN=sum(schedule.nitrogen_g_cm2(events));
    [afterC,afterN]=aomTotals(candidate.aom,flip(dzTop));
    candidate.cumulative.input_C_g_cm2= ...
        candidate.cumulative.input_C_g_cm2+addedC;
    candidate.cumulative.input_N_g_cm2= ...
        candidate.cumulative.input_N_g_cm2+addedN;
    candidate.cumulative=addCumulative(candidate.cumulative, ...
        'input_TAN_N_g_cm2',addedTanN);
    for material={'straw','manure'}
        name=material{1};
        candidate.cumulative=addCumulative(candidate.cumulative, ...
            ['input_' name '_C_g_cm2'],materialC.(name));
        candidate.cumulative=addCumulative(candidate.cumulative, ...
            ['input_' name '_N_g_cm2'],materialN.(name));
    end
    if materialC.root_residue>0||isfield(candidate.cumulative,'input_root_residue_C_g_cm2')
        candidate.cumulative=addCumulative(candidate.cumulative, ...
            'input_root_residue_C_g_cm2',materialC.root_residue);
        candidate.cumulative=addCumulative(candidate.cumulative, ...
            'input_root_residue_N_g_cm2',materialN.root_residue);
    end
    % Retain exact legacy residual arithmetic when the new route is off.
    carbonResidual=afterC-beforeC-expectedC;
    nitrogenResidual=afterN-beforeN+addedTanN-expectedN;
    if fastEnabled
        carbonResidual=afterC-beforeC+directC(3)-expectedC;
        nitrogenResidual=afterN-beforeN+addedTanN+directN(3)-expectedN;
    end
    diagnostic=struct('triggered',true,'event_index',events(1), ...
        'event_count',numel(events),'event_indices',events, ...
        'material_type',{materials}, ...
        'carbon_added_g_cm2',addedC,'nitrogen_added_g_cm2',addedN, ...
        'organic_nitrogen_added_g_cm2',addedOrganicN, ...
        'tan_nitrogen_added_g_cm2',addedTanN, ...
        'straw_carbon_added_g_cm2',materialC.straw, ...
        'straw_nitrogen_added_g_cm2',materialN.straw, ...
        'manure_carbon_added_g_cm2',materialC.manure, ...
        'manure_nitrogen_added_g_cm2',materialN.manure, ...
        'root_residue_carbon_added_g_cm2',materialC.root_residue, ...
        'root_residue_nitrogen_added_g_cm2',materialN.root_residue, ...
        'incorporated_layer_carbon_g_cm2',layerC, ...
        'incorporated_layer_nitrogen_g_cm2',layerN, ...
        'tan_layer_nitrogen_g_cm2',layerTanN, ...
        'direct_pom1_carbon_g_cm2',directC(1), ...
        'direct_pom2_carbon_g_cm2',directC(2), ...
        'direct_doc_carbon_g_cm2',directC(3), ...
        'direct_pon1_nitrogen_g_cm2',directN(1), ...
        'direct_pon2_nitrogen_g_cm2',directN(2), ...
        'direct_don_nitrogen_g_cm2',directN(3), ...
        'manure_fast_doc_C_g_cm3',manureFastC, ...
        'manure_fast_don_N_g_cm3',manureFastN, ...
        'manure_aom_C_added_g_cm3',manureAomC, ...
        'manure_aom_organic_N_added_g_cm3',manureAomN, ...
        'manure_tan_N_added_g_cm3',manureTanN, ...
        'manure_surface_aom_C_added_g_cm2',surfaceManureC, ...
        'manure_surface_aom_N_added_g_cm2',surfaceManureN, ...
        'carbon_residual_g_cm2',carbonResidual, ...
        'nitrogen_residual_g_cm2',nitrogenResidual);
    tolerance=1e-14*max(1,max(expectedC,expectedN));
    assert(abs(diagnostic.carbon_residual_g_cm2)<=tolerance ...
        &&abs(diagnostic.nitrogen_residual_g_cm2)<=tolerance, ...
        'stageOrganicFertilizer:massClosure', ...
        'Organic management event staging failed C/N closure.');
end

function [carbon,nitrogen]=aomTotals(aom,dzBottom)
% Integrate each sole material owner; component totals are caches, not extras.
% Editor: Wenjie Liu; Date: 2026-10-17.
    carbon=0;nitrogen=0;
    for material={'litter','straw','manure'}
        s=aom.(material{1});
        assert(all(isfinite([s.surface_C;s.surface_N;s.element_C(:);s.element_N(:)])) ...
            &&all([s.surface_C;s.surface_N;s.element_C(:);s.element_N(:)]>=0), ...
            'stageOrganicFertilizer:badAom','AOM inventories must be finite and nonnegative.');
        carbon=carbon+s.surface_C+sum(s.element_C(:).*dzBottom);
        nitrogen=nitrogen+s.surface_N+sum(s.element_N(:).*dzBottom);
    end
end

function diagnostic=emptyDiagnostic(n)
% Stable no-event diagnostic contract.
% Editor: Jizixing Zhu; Date: 2026-08-28.
    diagnostic=struct('triggered',false,'event_index',0, ...
        'event_count',0,'event_indices',[],'material_type',{{}}, ...
        'carbon_added_g_cm2',0,'nitrogen_added_g_cm2',0, ...
        'organic_nitrogen_added_g_cm2',0,'tan_nitrogen_added_g_cm2',0, ...
        'straw_carbon_added_g_cm2',0,'straw_nitrogen_added_g_cm2',0, ...
        'manure_carbon_added_g_cm2',0,'manure_nitrogen_added_g_cm2',0, ...
        'root_residue_carbon_added_g_cm2',0,'root_residue_nitrogen_added_g_cm2',0, ...
        'incorporated_layer_carbon_g_cm2',[], ...
        'incorporated_layer_nitrogen_g_cm2',[], ...
        'tan_layer_nitrogen_g_cm2',[], ...
        'direct_pom1_carbon_g_cm2',0,'direct_pom2_carbon_g_cm2',0, ...
        'direct_doc_carbon_g_cm2',0,'direct_pon1_nitrogen_g_cm2',0, ...
        'direct_pon2_nitrogen_g_cm2',0,'direct_don_nitrogen_g_cm2',0, ...
        'manure_fast_doc_C_g_cm3',zeros(n,1), ...
        'manure_fast_don_N_g_cm3',zeros(n,1), ...
        'manure_aom_C_added_g_cm3',zeros(n,1), ...
        'manure_aom_organic_N_added_g_cm3',zeros(n,1), ...
        'manure_tan_N_added_g_cm3',zeros(n,1), ...
        'manure_surface_aom_C_added_g_cm2',0, ...
        'manure_surface_aom_N_added_g_cm2',0, ...
        'carbon_residual_g_cm2',0,'nitrogen_residual_g_cm2',0);
end

function enabled=fastManureEnabled(p)
% Optional immutable physical switch; absence exactly retains legacy routing.
% Editor: Wenjie Liu; Date: 2026-10-30.
    enabled=false;
    if isempty(p)||~isstruct(p)||~isfield(p,'manure_fast_input_enabled'),return;end
    flag=p.manure_fast_input_enabled;
    assert(isscalar(flag)&&(islogical(flag) ...
        ||(isnumeric(flag)&&isreal(flag)&&isfinite(flag)&&ismember(flag,[0 1]))), ...
        'stageOrganicFertilizer:fastManureSwitch', ...
        'manure_fast_input_enabled must be a scalar logical or zero/one.');
    enabled=logical(flag);
end

function fraction=manureFastFraction(p)
% Third configured fraction splits paired event C/organic N, never total N.
% Editor: Wenjie Liu; Date: 2026-10-30.
    assert(isfield(p,'manure_input_fraction') ...
        &&isnumeric(p.manure_input_fraction)&&isreal(p.manure_input_fraction) ...
        &&isvector(p.manure_input_fraction)&&numel(p.manure_input_fraction)==3 ...
        &&all(isfinite(p.manure_input_fraction(:))) ...
        &&all(p.manure_input_fraction(:)>=0&p.manure_input_fraction(:)<=1) ...
        &&abs(sum(p.manure_input_fraction)-1)<1e-12 ...
        &&sum(p.manure_input_fraction(1:2))>0, ...
        'stageOrganicFertilizer:fastManureFraction', ...
        'Fast manure input requires normalized POM1/POM2/DOM fractions and a positive particulate share.');
    fraction=double(p.manure_input_fraction(3));
end

function validateDomCandidate(dom,n)
% Sole DOM owner stores bottom-up total phase mass; sorbed phases stay intact.
% Editor: Wenjie Liu; Date: 2026-10-30.
    assert(isstruct(dom)&&isscalar(dom), ...
        'stageOrganicFertilizer:fastManureDom','Fast manure input requires the live DOM owner.');
    for key={'doc_total_g_cm3','don_total_g_cm3'}
        name=key{1};
        assert(isfield(dom,name)&&isnumeric(dom.(name))&&isreal(dom.(name)) ...
            &&isequal(size(dom.(name)),[n,1])&&all(isfinite(dom.(name))) ...
            &&all(dom.(name)>=0), ...
            'stageOrganicFertilizer:fastManureDom', ...
            'Fast manure input requires nonnegative bottom-up NL-by-1 DOM totals.');
    end
end

function [weights,fractions]=rootProfile(context,n)
% Immutable root allocation, bottom-up layer mass fractions, not rate density.
% Editor: Wenjie Liu; Date: 2026-10-30.
    assert(isstruct(context)&&isscalar(context) ...
        &&isfield(context,'weights_bottom_up')&&isfield(context,'fractions'), ...
        'stageOrganicFertilizer:rootContext','Root events require explicit Jackson weights and component fractions.');
    weights=context.weights_bottom_up;fractions=context.fractions;
    assert(isnumeric(weights)&&isreal(weights)&&isvector(weights)&&numel(weights)==n ...
        &&all(isfinite(weights(:)))&&all(weights(:)>=0)&&abs(sum(weights(:))-1)<=1e-12, ...
        'stageOrganicFertilizer:rootWeights','Root weights must be finite normalized bottom-up layer fractions.');
    assert(isnumeric(fractions)&&isreal(fractions)&&isvector(fractions)&&numel(fractions)==3 ...
        &&all(isfinite(fractions(:)))&&all(fractions(:)>=0)&&abs(sum(fractions(:))-1)<=1e-12, ...
        'stageOrganicFertilizer:rootFractions','Root component fractions must be finite, nonnegative and sum to one.');
    weights=weights(:);fractions=fractions(:);
end

function cumulative=addCumulative(cumulative,name,amount)
% Add a material-specific ledger field while accepting older checkpoints.
% Editor: Codex; Date: 2026-08-28.
    if ~isfield(cumulative,name), cumulative.(name)=0; end
    cumulative.(name)=cumulative.(name)+amount;
end

function material=eventMaterial(schedule,event)
% Resolve the explicit native material identity; default organic events are manure.
% Editor: Jizixing Zhu; Date: 2026-08-28.
    material='manure';
    if isfield(schedule,'material_type')&&numel(schedule.material_type)>=event
        material=lower(strtrim(char(string(schedule.material_type{event}))));
    end
    if strcmp(material,'organic_fertilizer'), material='manure'; end
end

function weights=uniformDepthWeights(dzTop,depthStart,depthEnd)
% Exact layer-overlap fractions for uniform input per unit event depth.
% Editor: Jizixing Zhu; Date: 2026-08-28.
% Change: 2026-10-10 Wenjie Liu -- replace density-weighted mixing by depth overlap.
    assert(isfinite(depthStart)&&isfinite(depthEnd)&&depthStart>=0 ...
        &&depthEnd>depthStart&&depthEnd<=sum(dzTop)+1e-10, ...
        'stageOrganicFertilizer:eventDepth', ...
        'Incorporation depths must lie inside the COMPASS soil column.');
    top=[0;cumsum(dzTop(1:end-1))]; bottom=cumsum(dzTop);
    overlap=max(0,min(bottom,depthEnd)-max(top,depthStart));
    weights=overlap./sum(overlap);
    assert(abs(sum(weights)-1)<=1e-12, ...
        'stageOrganicFertilizer:depthClosure', ...
        'Event-to-layer overlap weights do not sum to one.');
end
