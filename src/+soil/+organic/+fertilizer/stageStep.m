function [candidate, diagnostic, candidateApplied] = stageStep( ...
        accepted, schedule, applied, stepStart, dtSeconds, zplus, bulkDensityTop)
%STAGESTEP Stage organic management events into native material AOM stores.
%
%   Event masses [g C/N cm^-2] retain their material identity. Manure TAN is
%   excluded from AOM and returned as a placement-resolved NH4 input ledger.
%   Surface events
%   remain in the surface AOM scalar; incorporated events are distributed to
%   COMPASS elements by exact depth overlap and stored as [g cm^-3 bulk soil].
%   No event enters POC, MOC, DOC, or DON directly. APPLIED is copied and is
%   committed only after the complete outer transaction accepts.
%
% Editor: Jizixing Zhu; Date: 2026-08-28.
% Change: 2026-08-30 Codex -- stage manure TAN separately from organic AOM N.

    candidate=accepted;
    candidateApplied=logical(applied(:).');
    diagnostic=emptyDiagnostic();
    if schedule.n_events==0, return; end
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
        'Top-down positive bulk density is required for soil-mass mixing.');
    addedC=0; addedN=0; addedOrganicN=0; addedTanN=0;
    materials=cell(1,numel(events));
    materialC=struct('straw',0,'manure',0);
    materialN=struct('straw',0,'manure',0);
    layerC=zeros(numel(dzTop),1); layerN=zeros(numel(dzTop),1);
    layerTanN=zeros(numel(dzTop),1);
    for q=1:numel(events)
        event=events(q);
        material=eventMaterial(schedule,event);
        if ~isfield(candidate.aom,material)
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
        assert(abs(organicNitrogen+tanNitrogen-nitrogen)<= ...
            1e-14*max(1,nitrogen), 'stageOrganicFertilizer:tanClosure', ...
            'Organic N plus TAN must equal event total N.');
        if strcmpi(schedule.placement{event},'surface')
            candidate.aom.(material).surface_C = ...
                candidate.aom.(material).surface_C+carbon;
            candidate.aom.(material).surface_N = ...
                candidate.aom.(material).surface_N+organicNitrogen;
            layerTanN(1)=layerTanN(1)+tanNitrogen;
        else
            eta=1;
            if isfield(schedule,'mixing_efficiency') ...
                    &&numel(schedule.mixing_efficiency)>=event
                eta=schedule.mixing_efficiency(event);
            end
            assert(isfinite(eta)&&eta>=0&&eta<=1, ...
                'stageOrganicFertilizer:mixingEfficiency', ...
                'Mixing efficiency must lie within [0,1].');
            weights=soilMassWeights(dzTop,bulkDensityTop, ...
                schedule.depth_start_cm(event),schedule.depth_end_cm(event));
            incorporatedC=eta.*carbon;
            incorporatedN=eta.*organicNitrogen;
            surfaceC=(1-eta).*carbon;
            surfaceN=(1-eta).*organicNitrogen;
            cTop=incorporatedC.*weights./dzTop;
            nTop=incorporatedN.*weights./dzTop;
            candidate.aom.(material).element_C = ...
                candidate.aom.(material).element_C+flip(cTop);
            candidate.aom.(material).element_N = ...
                candidate.aom.(material).element_N+flip(nTop);
            candidate.aom.(material).surface_C = ...
                candidate.aom.(material).surface_C+surfaceC;
            candidate.aom.(material).surface_N = ...
                candidate.aom.(material).surface_N+surfaceN;
            layerC=layerC+cTop.*dzTop;
            layerN=layerN+nTop.*dzTop;
            layerTanN=layerTanN+eta.*tanNitrogen.*weights;
            layerTanN(1)=layerTanN(1)+(1-eta).*tanNitrogen;
        end
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
        'incorporated_layer_carbon_g_cm2',layerC, ...
        'incorporated_layer_nitrogen_g_cm2',layerN, ...
        'tan_layer_nitrogen_g_cm2',layerTanN, ...
        'carbon_residual_g_cm2',addedC-expectedC, ...
        'nitrogen_residual_g_cm2',addedN-expectedN);
    tolerance=1e-14*max(1,max(expectedC,expectedN));
    assert(abs(diagnostic.carbon_residual_g_cm2)<=tolerance ...
        &&abs(diagnostic.nitrogen_residual_g_cm2)<=tolerance, ...
        'stageOrganicFertilizer:massClosure', ...
        'Organic management event staging failed C/N closure.');
end

function diagnostic=emptyDiagnostic()
% Stable no-event diagnostic contract.
% Editor: Jizixing Zhu; Date: 2026-08-28.
    diagnostic=struct('triggered',false,'event_index',0, ...
        'event_count',0,'event_indices',[],'material_type',{{}}, ...
        'carbon_added_g_cm2',0,'nitrogen_added_g_cm2',0, ...
        'organic_nitrogen_added_g_cm2',0,'tan_nitrogen_added_g_cm2',0, ...
        'straw_carbon_added_g_cm2',0,'straw_nitrogen_added_g_cm2',0, ...
        'manure_carbon_added_g_cm2',0,'manure_nitrogen_added_g_cm2',0, ...
        'incorporated_layer_carbon_g_cm2',[], ...
        'incorporated_layer_nitrogen_g_cm2',[], ...
        'tan_layer_nitrogen_g_cm2',[], ...
        'carbon_residual_g_cm2',0,'nitrogen_residual_g_cm2',0);
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

function weights=soilMassWeights(dzTop,bulkDensityTop,depthStart,depthEnd)
% Conservative SWAT-style tillage mixing weights based on mixed soil mass.
% Editor: Jizixing Zhu; Date: 2026-08-28.
    assert(isfinite(depthStart)&&isfinite(depthEnd)&&depthStart>=0 ...
        &&depthEnd>depthStart&&depthEnd<=sum(dzTop)+1e-10, ...
        'stageOrganicFertilizer:eventDepth', ...
        'Incorporation depths must lie inside the COMPASS soil column.');
    top=[0;cumsum(dzTop(1:end-1))]; bottom=cumsum(dzTop);
    overlap=max(0,min(bottom,depthEnd)-max(top,depthStart));
    mixedSoilMass=overlap.*bulkDensityTop;
    weights=mixedSoilMass./sum(mixedSoilMass);
    assert(abs(sum(weights)-1)<=1e-12, ...
        'stageOrganicFertilizer:depthClosure', ...
        'Event-to-layer overlap weights do not sum to one.');
end
