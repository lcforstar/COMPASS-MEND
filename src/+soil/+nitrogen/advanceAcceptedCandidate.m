function nitrogen=advanceAcceptedCandidate(varargin)
%ADVANCEACCEPTEDCANDIDATE Couple detached MEND C/N through explicit sources.
%   One source MEND candidate supplies NH4/NO3 exchange and respiratory O2
%   demand. Replay the reference oxygen equations on this hour's accepted
%   carrier tiles; infeasible demands reject before the atomic commit.
% Change: 2026-10-23 Wenjie Liu -- restore reference oxygen/nitrification
%   and remove whole-MEND oxygen scaling; retain source/sink and mass ledgers.
% Editor: Jizixing Zhu; Date: 2026-09-14.
% Editor: Wenjie Liu; Date: 2026-10-17.
% Change: 2026-10-17 Wenjie Liu -- pair AOM products with configured POM2
%   C:N; align DOM water drivers with outer DIN while preserving particle tiles.
% Change: 2026-10-30 Wenjie Liu -- opt-in material-specific manure inputs and
%   independent passive organic/TAN labels on the same accepted gross maps.
% Change: 2026-11-01 Wenjie Liu -- optional clay-scaled NH4 accessibility
%   reserves the same nodal mass for source MEND and accepted O2/N replay.
% Change: 2026-11-03 Wenjie Liu -- drive free DOM with chronological liquid
%   fluxes and actual liquid storage endpoints from the accepted water tiles.
    assert(numel(varargin)>=7&&numel(varargin)<=12);
    args=varargin;
    if numel(args)<12,args(end+1:12)=repmat({[]},1,12-numel(args));end
    nitrogen=advanceWithSourceSinks(args{:});
    nitrogen.diagnostics.oxygen_reaction_factor_top_down=ones(args{1}.domain.NL,1);
    nitrogen.diagnostics.oxygen_coupling_iterations=1;
end

function nitrogen = advanceWithSourceSinks( ...
        built, root_request, soil_advance, surface_water_request, fertilizer, ...
        kt, dt_s, surface_water_candidate, organic_accepted, ...
        dom_accepted, microbe_accepted, organic_management)
%ADVANCEACCEPTEDCANDIDATE  Form one detached accepted soil-N endpoint.
%
%   This adapter preserves the existing reaction, freeze-partition, O2, and
%   ADE kernels.  It supplies their legacy arrays as two time levels only:
%   slot 1 is accepted t^n and slot 2 is the candidate t^(n+1).  Therefore N
%   history has explicit soil ownership in Y.soil_nitrogen rather than being
%   hidden in full-run diagnostic vectors.  No field of BUILT is mutated.
%
%   The soil water/heat/gas solve is already accepted when this routine runs.
%   Its accepted-time-mean liquid flux drives N advection; mean absolute
%   liquid flux drives dispersion, preserving the reference averaging order.
%   N is not allowed to rerun or alter the coupled h/U/P_g solution.
%
%   Inputs:
%     built [in] struct Accepted state, immutable parameters, and current forcing.
%     root_request [in] struct NH4/NO3 sink request [ug cm^-3 s^-1], soil-to-root positive.
%     soil_advance [in] struct Accepted water fluxes and final coupled diagnostics.
%     surface_water_request [in] struct Current top water boundary request.
%     fertilizer [in] struct Accept-only fertilizer ledger event.
%     kt [in] [-] Outer forcing index; dt_s [in] [s] accepted outer duration.
%
%   Returns:
%     nitrogen [out] struct Detached t^(n+1) soil-N state, pools, O2, and ledger.
%
%   Discrete closure:
%     Liquid and ice storage at t^n are carried in Y.soil_nitrogen; all
%     reaction, transport, boundary, fertilizer, and root terms advance once
%     over dt_s.  The candidate is committed only with the accepted soil step.
%
%   Reference: retained STEMMUS-N solveNitrogen kernel and its finite-volume
%   NH4/NO3 mass ledger (src/+soil/+nitrogen/solveNitrogen.m).
%
%   Editor: Jizixing Zhu; Date: 2026-07-20.
% Change: 2026-07-20 Jizixing Zhu -- treat disabled root-N requests as strict zero soil sources.
% Change: 2026-07-20 Jizixing Zhu — preserve pre-fertilizer accepted pools
%   and pass management inputs into the N ledger before residual evaluation.
% Change: 2026-07-20 Jizixing Zhu -- use the accepted outer candidate for surface
%   pond storage without prematurely writing surface state.
% Change: 2026-08-29 Codex -- stage native AOM/DOM/MEMC/MEND proposals and
%   pass signed NH4/NO3 exchange into the retained transport solver.
% Change: 2026-08-30 Codex -- inject event TAN through the native NH4 source
%   interface at its conservative placement depth.
% Change: 2026-09-05 Jizixing Zhu -- pass accepted endpoint water and NH4
%   sorption capacity with the matching accepted liquid DIN to the C/N ledger.
% Change: 2026-09-12 Codex -- retain the post-root DIN shadow inventory in
%   the accepted diagnostic result for direct plant--microbe competition audit.
% Change: 2026-09-12 Codex -- retain MEND's post-reaction DOC when applying
%   same-stage organic-N transfers to the sole-owner DOM candidate.
% Change: 2026-09-12 Codex -- drive MEND carbon reactions with the same
%   accepted soil temperature and water potential used by MEND DIN kinetics.
% Change: 2026-09-12 Codex -- advance organic N on the accepted MEND carbon
%   RK4 substep ledger instead of replaying the integrated hourly flux total.
% Change: 2026-09-13 Codex -- apply enabled BIOPARTICLE C/N transport after
%   MEND reaction and N limitation, before the atomic soil-state commit.
% Change: 2026-09-13 Codex -- provide the DOM-owned DON Kd profile to the
%   MEND RK-stage aqueous organic-N transfer.
% Change: 2026-09-13 Codex -- pass source kinetic DIN uptake through the
%   root-first inventory ledger and return post-uptake excess MBA N to NH4.
% Change: 2026-09-14 Codex -- evaluate MEND kinetics on free NH4 after root
%   reservation while retaining total NH4 for the conservative inventory cap.
% Change: 2026-09-14 Codex -- rebuild MEND C/N from accepted owners under
%   oxygen supply; separate surface respiration from soil-pore oxygen demand.
% Change: 2026-09-14 Jizixing Zhu -- select paired source-hour dynamic DOM/C/N
%   before the same root-reserved inorganic solve and oxygen retry gate.
% Change: 2026-09-21 Wenjie Liu -- drive DOM/BIOPARTICLE transport with the
%   chronological accepted soil-tile QL/theta schedule and drive the outer-step
%   inorganic-N solve with accepted-dt-weighted liquid-flux components.
% Change: 2026-09-21 Wenjie Liu -- retain the exact liquid-flux driver passed
%   to the inorganic-N ADE so production regressions can distinguish the
%   accepted-dt mean from the legacy final-tile fallback.
% Change: 2026-09-21 Wenjie Liu -- advance a detached mineral-fertilizer N
%   source label through accepted gross MEND and frozen transport operators.
% Change: 2026-10-12 Wenjie Liu -- bound source-hour MEND allocation by
%   shared-node DIN while preserving physical/kinetic inventories and C/N feedback.
    if nargin < 8 || isempty(surface_water_candidate)
        surface_water_candidate = surface.currentWaterState(built);
    end
    assert(isstruct(root_request) && isfield(root_request, 'sink_nh4') ...
        && isfield(root_request, 'sink_no3') && isfield(root_request, 'enabled'), ...
        'soil:nitrogen:missingRootRequest', ...
        'Accepted soil-N advance requires the typed vegetation sink request.');
    assert(isstruct(soil_advance) && isfield(soil_advance, 'solver_diagnostics'), ...
        'soil:nitrogen:missingSoilDiagnostics', ...
        'Accepted soil-N advance requires final coupled-soil flux diagnostics.');
    assert(isscalar(kt) && isfinite(kt) && kt >= 1 && ...
        isscalar(dt_s) && isfinite(dt_s) && dt_s > 0, ...
        'soil:nitrogen:badStep', 'kt and dt_s must be finite positive scalars.');

    [Y, ~] = soil.nitrogen.materializeAcceptedState( ...
        built.Y, built.workspace.derived, built.domain);
    state = Y.soil_nitrogen;
    state.variables = updateTotalNitrogen( ...
        state.variables, Y.nh4, Y.no3, state.theta_liquid_prev, ...
        built.parameters.soil_n, built.domain);
    cfg = nitrogenCfg(built);
    tracerOn=isfield(built.control.switches,'lFertilizerSourceTracer') ...
        &&built.control.switches.lFertilizerSourceTracer==1;
    manureTracerOn=isfield(built.control.switches,'lManureSourceTracer') ...
        &&built.control.switches.lManureSourceTracer==1;
    anyTracerOn=tracerOn||manureTracerOn;
    soil_variables = soilVariables(built, Y, state);
    nitrogen_variables = nitrogenVariables(built, state, root_request, kt);
    nitrogen_flux = fluxWork(state.flux);
    nitrogen_flux = seedInitialStorage( ...
        nitrogen_flux, built, Y, state, nitrogen_variables);
    nitrogen_boundary = nitrogenBoundary( ...
        built, soil_advance, surface_water_request, ...
        surface_water_candidate, kt);
    [air_variables,liquidTransportDriver] = ...
        airVariables(soil_advance,built.domain.NL);
    gas_variables = struct();
    oxygenScale=ones(built.domain.NL,1); % source MEND map is never oxygen-rescaled

    externalInputs = fertilizer;
    tanRate=zeros(built.domain.NL,1);
    if nargin>=12&&isstruct(organic_management) ...
            &&isfield(organic_management,'tan_layer_nitrogen_g_cm2') ...
            &&~isempty(organic_management.tan_layer_nitrogen_g_cm2)
        tanLayer=organic_management.tan_layer_nitrogen_g_cm2(:);
        assert(numel(tanLayer)==built.domain.NL&&all(isfinite(tanLayer)) ...
            &&all(tanLayer>=0),'advanceAcceptedCandidate:tanLayer', ...
            'Organic-fertilizer TAN layer masses must be finite NL-by-1.');
        tanRate=flip(tanLayer./flip(built.domain.dz(:))).*1e6./dt_s;
    end
    organicCandidate = [];
    organicExchange = [];
    organicDiagnostic = [];
    aomCandidate=[];manureAomBefore=[];
    manurePartition=[];manurePartitionDonor=[];
    domCandidate=[]; microbeCandidate=[];
    tracerTotalBeforeReaction=[];tracerTotalAfterReaction=[];
    domOn=isfield(built.control.switches,'lSoilDOM') ...
        &&built.control.switches.lSoilDOM==1;
    microbeOn=isfield(built.control.switches,'lExplicitSoilMicrobes') ...
        &&built.control.switches.lExplicitSoilMicrobes==1;
    thresholdEnabled=isfield(built.parameters.soil_n,'nh4_threshold_enabled') ...
        &&built.parameters.soil_n.nh4_threshold_enabled;
    protectedProfile=[];
    if thresholdEnabled
        % Nitrification keeps its own threshold and conservative node budget.
        % Earlier positive-threshold sites also protected MEND by default.
        if ~isfield(built.parameters.soil_n,'nh4_protection_scope') ...
                ||strcmp(built.parameters.soil_n.nh4_protection_scope,'nitrification_and_mend')
            protectedProfile=built.parameters.soil_n.S_min_nh4_mg_kg;
        end
        if microbeOn
            assert(strcmp(built.parameters.soil_microbes.reaction_method,'source_hourly'), ...
                'soil:nitrogen:unsupportedProtectionMethod', ...
                'Nonzero NH4 protection supports source_hourly microbes only.');
        end
    end
    microbeTransportOn=isfield(built.control.switches,'lMicrobialTransport') ...
        &&built.control.switches.lMicrobialTransport==1;
    thetaElement=mean(built.workspace.derived.Theta_LL,2);
    particleDrivers=[];
    if microbeTransportOn
        [particleTheta,particleQ,particleDt]=acceptedTransportDrivers( ...
            soil_advance,thetaElement,dt_s,built.domain.NL);
        particleDrivers=struct('theta',particleTheta,'q',particleQ,'dt',particleDt);
    end
    organicOn=isfield(built.control.switches,'lSoilOrganic') ...
        &&built.control.switches.lSoilOrganic==1;
    if organicOn
        if nargin<9||isempty(organic_accepted)
            organic_accepted=built.Y.soil_organic;
        end
        trial=built.workspace.step.trial;
        temperatureElement=0.5.*(trial.T(1:end-1)+trial.T(2:end));
        headElement=0.5.*(trial.h(1:end-1)+trial.h(2:end));
        if manureTracerOn,manureAomBefore=organic_accepted.aom.manure;end
        aomCandidate=soil.aom.advanceAOM(organic_accepted, ...
            temperatureElement,headElement,thetaElement,built.domain.dz(:), ...
            built.parameters.soil_organic.aom,dt_s);
        organicCandidate=aomCandidate.state;
        organicExchange=struct('aom',aomCandidate, ...
            'nh4_ug_cm3_s',aomCandidate.tan_N_rate_ug_cm3_s, ...
            'no3_ug_cm3_s',zeros(built.domain.NL,1), ...
            'tan_nh4_ug_cm3_s',tanRate, ...
            'co2_c_ug_cm3_s',aomCandidate.co2_C_rate_ug_cm3_s);
        organicDiagnostic=struct('aom_carbon_residual_g_cm2', ...
            aomCandidate.carbon_residual_g_cm2, ...
            'aom_nitrogen_residual_g_cm2', ...
            aomCandidate.nitrogen_residual_g_cm2, ...
            'aom_surface_tan_N_g_cm2',sum( ...
                aomCandidate.tan_N_g_cm3.*built.domain.dz(:)));
        if isempty(dom_accepted), dom_accepted=built.Y.soil_dom; end
        dom_accepted.doc_total_g_cm3=dom_accepted.doc_total_g_cm3 ...
            +aomCandidate.doc_C_g_cm3;
        dom_accepted.don_total_g_cm3=dom_accepted.don_total_g_cm3 ...
            +aomCandidate.don_N_g_cm3;
        if isempty(microbe_accepted)
            microbe_accepted=built.Y.soil_microbes;
        end
        bdTop=flip(built.parameters.soil.bulk_density(:));
        p=built.parameters.soil_microbes;
        if isfield(p,'manure_fast_input_enabled')&&p.manure_fast_input_enabled
            % The slow manure branch retains its own particulate allocation.
            % Litter/straw remain on the existing generic mapping. Never
            % apply manure fractions to all coincident AOM materials.
            routes=aomCandidate.material_routes;
            otherC=routes.litter.poc_C+routes.straw.poc_C;
            otherN=routes.litter.pon_N+routes.straw.pon_N;
            other=soil.aom.partitionMendProducts(otherC,otherN,p.input_fraction,p.CN_POM2);
            manurePartition=soil.aom.partitionMendProducts(routes.manure.poc_C, ...
                routes.manure.pon_N,p.manure_input_fraction,p.CN_POM2);
            aomPartition=struct('C_POM1',other.C_POM1+manurePartition.C_POM1, ...
                'C_POM2',other.C_POM2+manurePartition.C_POM2, ...
                'N_PON1',other.N_PON1+manurePartition.N_PON1, ...
                'N_PON2',other.N_PON2+manurePartition.N_PON2);
            manurePartitionDonor=routes.manure.pon_N;
        else
            % Preserve the complete legacy arithmetic with the option off.
            aomPartition=soil.aom.partitionMendProducts(aomCandidate.poc_C_g_cm3, ...
                aomCandidate.pon_N_g_cm3,p.input_fraction,p.CN_POM2);
            manurePartition=aomPartition;
            manurePartitionDonor=aomCandidate.pon_N_g_cm3;
        end
        microbe_accepted.POM1=microbe_accepted.POM1 ...
            +flip(aomPartition.C_POM1)./bdTop.*1000;
        microbe_accepted.POM2=microbe_accepted.POM2 ...
            +flip(aomPartition.C_POM2)./bdTop.*1000;
        microbe_accepted.POC=microbe_accepted.POM1+microbe_accepted.POM2;
        organicCandidate.PON1=organicCandidate.PON1 ...
            +aomPartition.N_PON1;
        organicCandidate.PON2=organicCandidate.PON2 ...
            +aomPartition.N_PON2;
        organicCandidate.PON=organicCandidate.PON1+organicCandidate.PON2;
        organic_accepted=organicCandidate;
    end
    if domOn
        if nargin<10||isempty(dom_accepted), dom_accepted=built.Y.soil_dom; end
        % Change: 2026-09-29 Wenjie Liu -- pass water-owned bottom faces on
        % the same accepted chronology, including every oxygen retry.
        [domWater,domHourly]=soil.dom.waterDriverSchedule(soil_advance,thetaElement, ...
            dt_s,liquidTransportDriver,built.parameters.soil_dom);
        transportTheta=domWater.theta;transportQ=domWater.QL;transportDt=domWater.dt_s;
        domCandidate=soil.dom.advanceAcceptedCandidate(dom_accepted, ...
            transportTheta,built.parameters.soil.bulk_density(:),built.domain.dz(:), ...
            transportQ,built.parameters.soil_dom,transportDt,domWater.bottom_faces,domWater.theta_initial);
    end
    [rootReservation,rootRemaining] = rootDinReservation( ...
        root_request,dt_s,Y.nh4,Y.no3,state.theta_liquid_prev, ...
        built.parameters.soil_n.RHOKG,built.parameters.soil.bulk_density(:));
    if microbeOn
        if nargin<11||isempty(microbe_accepted)
            microbe_accepted=built.Y.soil_microbes;
        end
        if isempty(organic_accepted)
            organic_accepted=built.Y.soil_organic;
        end
        domBeforeMicrobe=domCandidate.state;
        microbeParameters=built.parameters.soil_microbes;
        microbeParameters.reaction_scale_top_down=oxygenScale;
        sourceHourly=isfield(microbeParameters,'reaction_method') ...
            &&strcmp(microbeParameters.reaction_method,'source_hourly');
        if sourceHourly
            [rootRemaining,mendNodalInventory]=soil.nitrogen.limitMendDinInventoryAtNodes( ...
                rootRemaining,root_request,Y.nh4,Y.no3,state.theta_liquid_prev, ...
                built.parameters.soil_n.RHOKG,built.parameters.soil.bulk_density(:), ...
                built.domain.dz(:),dt_s,protectedProfile);
            tracerTotalBeforeReaction=tracerMendBundle( ...
                organic_accepted,domBeforeMicrobe,built.parameters.soil_dom);
            dinAfterRoot=soil.nitrogen.mendKineticDinAfterRoot( ...
                rootRemaining,rootReservation,Y.nh4,Y.no3,state.theta_liquid_prev, ...
                built.parameters.soil_n.RHOKG,built.parameters.soil.bulk_density(:),protectedProfile);
            [microbeCandidate,organicCandidate,domCandidate.state,mendN,nFeedback]= ...
                soil.microbe.advanceSourceCandidate(microbe_accepted,organic_accepted, ...
                domBeforeMicrobe,thetaElement,built.parameters.soil.bulk_density(:), ...
                built.parameters.soil_dom,microbeParameters,dt_s, ...
                temperatureElement,headElement,dinAfterRoot,rootRemaining,oxygenScale);
            tracerTotalAfterReaction=tracerMendBundle( ...
                organicCandidate,domCandidate.state,built.parameters.soil_dom);
            organicExchange.mend=mendN;
            mendAllocation=mendN.din_allocation;mendKinetic=mendN.din_kinetic_request;
            organicDiagnostic.mend_nitrogen_residual_mg_g=mendN.nitrogen_residual_mg_g;
            organicDiagnostic.mend_nitrogen_feedback=nFeedback;
            organicDiagnostic.mend_din_kinetic_request=mendKinetic;
        else
        microbeCandidate=soil.microbe.advanceAcceptedCandidate( ...
            microbe_accepted,domBeforeMicrobe,thetaElement, ...
            built.parameters.soil.bulk_density(:), ...
            built.parameters.soil_dom,microbeParameters,dt_s, ...
            temperatureElement,headElement);
        domCandidate.state.doc_total_g_cm3=microbeCandidate.doc_total_g_cm3;
        [organicCandidate,domNitrogenCandidate,mendN] = ...
            soil.microbe.mendOrganicNitrogenCandidate( ...
            organic_accepted,domBeforeMicrobe,microbe_accepted, ...
            microbeCandidate.state,microbeCandidate.substep_ledger, ...
            thetaElement,built.parameters.soil.bulk_density(:), ...
            built.parameters.soil_microbes,built.parameters.soil_dom);
        % The N ledger reads DOC from each accepted carbon RK stage, but owns
        % only the DON update. Keep MEND's post-reaction DOC as the sole
        % accepted carbon state instead of restoring domBeforeMicrobe.
        assert(isequaln(domNitrogenCandidate.doc_total_g_cm3, ...
            domBeforeMicrobe.doc_total_g_cm3), ...
            'advanceAcceptedCandidate:mendNitrogenChangedDoc', ...
            'MEND organic-N transfer must not modify the DOC owner.');
        domCandidate.state.don_total_g_cm3= ...
            domNitrogenCandidate.don_total_g_cm3;
        organicExchange.mend = mendN;
        dinAfterRoot=soil.nitrogen.mendKineticDinAfterRoot( ...
            rootRemaining,rootReservation,Y.nh4,Y.no3, ...
            state.theta_liquid_prev,built.parameters.soil_n.RHOKG, ...
            built.parameters.soil.bulk_density(:));
        mendEnvironment=soil.microbe.mendEnvironmentalParameters( ...
            built.parameters.soil_microbes,flip(temperatureElement), ...
            flip(headElement).*9.80665e-5, ...
            built.parameters.soil_microbes.pH_top_down);
        activeNitrogen=flip(organicCandidate.MB_N(:)) ...
            ./flip(built.parameters.soil.bulk_density(:)).*1000;
        mendEnvironment.VNup_MB_h=mendEnvironment.VNup_MB_h.*oxygenScale;
        mendKinetic=soil.microbe.mendDinKineticRequest( ...
            dinAfterRoot,microbeCandidate.state.MB,activeNitrogen, ...
            mendN.din_request_mg_g,mendEnvironment,dt_s);
        mendAllocation=soil.nitrogen.allocateMendDinStages( ...
            mendKinetic.requested_mg_g,mendN.din_request_mg_g, ...
            rootRemaining,mendN.din_release_mg_g);
        organicExchange.mend.din_kinetic_request=mendKinetic;
        organicExchange.mend.substep_coupling_policy=...
            'converged_frozen_carbon_trajectory_then_atomic_CN_DIN_ledger';
        organicExchange.mend.dt_s=dt_s;
        mendActual=mendAllocation.nh4_actual_mg_g+...
            mendAllocation.no3_actual_mg_g;
        [microbeCandidate.state,organicCandidate,nFeedback]=...
            soil.microbe.applyNitrogenLimitation(microbeCandidate.state, ...
            organicCandidate,mendActual,built.parameters.soil.bulk_density(:), ...
            built.parameters.soil_microbes,dt_s);
        mendAllocation.nh4_release_mg_g=...
            mendAllocation.nh4_release_mg_g+...
            nFeedback.nh4_post_uptake_release_mg_g;
        mendAllocation.net_nh4_sink_mg_g=...
            mendAllocation.nh4_actual_mg_g-mendAllocation.nh4_release_mg_g;
        organicExchange.mend.din_allocation=mendAllocation;
        microbeCandidate.integrated_flux_mg_g.N_overflow_CO2=...
            nFeedback.overflow_co2_mg_g;
        microbeCandidate.co2_c_rate_ug_cm3_s=...
            microbeCandidate.co2_c_rate_ug_cm3_s+...
            nFeedback.co2_c_rate_ug_cm3_s;
        microbeCandidate.carbon_residual_mg_g=...
            microbeCandidate.carbon_residual_mg_g+...
            nFeedback.carbon_residual_mg_g;
        organicDiagnostic.mend_nitrogen_residual_mg_g = ...
            mendN.nitrogen_residual_mg_g;
        organicDiagnostic.mend_nitrogen_feedback=nFeedback;
        organicDiagnostic.mend_din_kinetic_request=mendKinetic;
        end
        if microbeTransportOn
            % Change: 2026-09-29 Wenjie Liu -- particle policy is independent
            % of DOM; both consume identical water-owned accepted bottom faces.
            particleFaces=soil.dom.acceptedBottomDrivers(soil_advance, ...
                built.parameters.soil_microbes.bioparticle,particleDt);
            particle=soil.microbe.bioparticleTransport( ...
                microbeCandidate.state,organicCandidate,particleTheta, ...
                built.parameters.soil.bulk_density(:),built.domain.dz(:), ...
                particleQ,built.parameters.soil_microbes.bioparticle,particleDt,[],particleFaces);
            microbeCandidate.state=particle.microbes;
            organicCandidate=particle.organic_n;
            microbeCandidate.bioparticle_transport=particle;
            organicDiagnostic.bioparticle_carbon_residual_g_cm2= ...
                particle.carbon_residual_g_cm2;
            organicDiagnostic.bioparticle_nitrogen_residual_g_cm2= ...
                particle.nitrogen_residual_g_cm2;
        end
        % MEND owns the complete organic-carbon/microbial reaction.  Its
        % same-stage donor-N transfers update the paired organic pools and DON;
        % DIN requests remain a separate, root-prioritized solver concern.
        externalInputs.organic_enabled=true;
        externalInputs.organic_nh4_rate=domCandidate.nh4_source_ug_cm3_s ...
            +aomCandidate.tan_N_rate_ug_cm3_s+tanRate ...
            -flip(mendAllocation.net_nh4_sink_mg_g ...
            .*flip(built.parameters.soil.bulk_density(:)).*1000./dt_s);
        externalInputs.organic_no3_rate=...
            -flip(mendAllocation.no3_actual_mg_g ...
            .*flip(built.parameters.soil.bulk_density(:)).*1000./dt_s);
        externalInputs.organic_co2_c_rate= ...
            microbeCandidate.co2_c_rate_ug_cm3_s;
        % AOM's direct respiration is surface-only and uses atmospheric O2.
        % Incorporated AOM has no direct CO2 branch; MEND owns soil respiration.
    elseif organicOn
        error('advanceAcceptedCandidate:soilCarbonWithoutMemc', ...
            'The native soil-carbon path requires DOM and MEND.');
    end
    % O2 is replayed from t^n on the exact accepted water/heat/gas tiles.
    % No previous-hour respiration or separately writable O2 state is used.
    respiration=zeros(built.domain.NL,1);
    if microbeOn,respiration=externalInputs.organic_co2_c_rate;end
    assert(isfield(soil_advance.solver_diagnostics,'oxygen_replay_schedule'), ...
        'soil:nitrogen:missingOxygenReplay', ...
        'Source/sink coupling requires every accepted carrier tile.');
    if thresholdEnabled
        nodeBudget=protectedNitrificationBudget( ...
            built,Y,state,root_request,externalInputs,dt_s);
        oxygenHour=soil.gas.replayOxygenInterval(built, ...
            soil_advance.solver_diagnostics.oxygen_replay_schedule,respiration,[],nodeBudget);
    else
        oxygenHour=soil.gas.replayOxygenInterval(built, ...
            soil_advance.solver_diagnostics.oxygen_replay_schedule,respiration);
    end
    assert(oxygenHour.feasible,'soil:nitrogen:oxygenSourceInfeasible', ...
        'Reference oxygen transport cannot realize this hour''s fixed C/N sources: %s', ...
        oxygenHour.failure_reason);
    externalInputs.oxygen_hour=oxygenHour;
    if anyTracerOn
        nitrogen_variables.fertilizer_tracer_total_gross= ...
            fertilizerTracerTotalGross(externalInputs,mendAllocation, ...
            built.parameters.soil.bulk_density(:),dt_s);
    end

    [nitrogen_variables, soil_variables, nitrogen_flux, diagnostics] = ...
        soil.nitrogen.solveNitrogen(nitrogen_variables, nitrogen_flux, ...
            nitrogen_boundary, soil_variables, air_variables, gas_variables, ...
            cfg, 2, dt_s, externalInputs);
    if microbeOn&&sourceHourly
        diagnostics.mend_nodal_inventory=mendNodalInventory;
    end
    if thresholdEnabled,diagnostics.nh4_protection=oxygenHour.nh4_protection;end
    nitrogen_variables = updateTotalNitrogen( ...
        nitrogen_variables, soil_variables.C_NH4_new, ...
        soil_variables.C_NO3_new, soil_variables.Theta_LL_for_N_prev, ...
        built.parameters.soil_n, built.domain);

    nitrogen = struct();
    nitrogen.state = stateFromWork(state, nitrogen_variables, nitrogen_flux, ...
        soil_variables);
    nitrogen.trial = struct( ...
        'nh4', soil_variables.C_NH4_new, ...
        'no3', soil_variables.C_NO3_new, ...
        'urea', soil_variables.C_Urea_new, ...
        'C_NH4_granule', soil_variables.C_NH4_granule);
    nitrogen.ice = struct('nh4', soil_variables.C_NH4_ice, ...
        'no3', soil_variables.C_NO3_ice);
    nitrogen.oxygen = nitrogen_variables.C_O2; % derived concentration
    nitrogen.oxygen_y = oxygenHour.final_y_O2;
    nitrogen.oxygen_phase = oxygenHour.phase_new;
    nitrogen.diagnostics = diagnostics;
    nitrogen.diagnostics.liquid_transport_driver=liquidTransportDriver;
    if domOn
        % Existing consumers read these four outer-step summaries as vectors
        % and scalars. The actual chronological operator input is in .tiles.
        domHourly.tiles=domWater;
        nitrogen.diagnostics.dom_transport_driver=domHourly;
    end
    nitrogen.boundary = nitrogen_boundary;
    nitrogen.root_request = root_request;
    nitrogen.root_reservation = rootReservation;
    nitrogen.root_remaining = rootRemaining;
    sinkActual = nitrogen_variables.SinkS;
    % SinkS is layer x endpoint x species [ug N cm^-3 s^-1].
    % Integrate each layer after averaging endpoints and adding NH4/NO3;
    % collapsing the layer dimension before dz weighting both overcounts
    % uptake and leaves a two-species vector instead of a scalar.
    % Change: 2026-09-30 Wenjie Liu -- actual plant inflow is one physical
    % ground-area amount, shared by soil and winter-wheat candidate replay.
    actualRootRate = sum(sum(mean(sinkActual,2),3).*built.domain.dz(:))*0.01;
    nitrogen.root_actual = struct('sink_nh4',sinkActual(:,:,1), ...
        'sink_no3',sinkActual(:,:,2), ...
        'U_N_total',actualRootRate);
    nitrogen.organic_exchange = organicExchange;
    nitrogen.organic_diagnostics = organicDiagnostic;
    if microbeOn && isfield(organicDiagnostic,'mend_nitrogen_residual_mg_g')
        actual=mendAllocation.nh4_actual_mg_g+mendAllocation.no3_actual_mg_g;
        nitrogen.organic_diagnostics.mend_din_actual_mg_g=actual;
        nitrogen.organic_diagnostics.mend_din_unmet_mg_g=...
            mendAllocation.unmet_mg_g;
    end
    nitrogen.organic = organicCandidate;
    nitrogen.dom = domCandidate;
    nitrogen.microbes = microbeCandidate;
    if tracerOn&&~manureTracerOn
        nitrogen.fertilizer_tracer=soil.nitrogen.tracer.advanceCandidate( ...
            built,Y,fertilizer,soil_variables,nitrogen_variables,diagnostics, ...
            rootReservation,rootRemaining,tracerTotalBeforeReaction, ...
            tracerTotalAfterReaction,organicCandidate,domCandidate, ...
            microbeCandidate,organicExchange.mend.source_hour_ledger, ...
            transportTheta,transportQ,transportDt,dt_s,particleDrivers);
    end
    if manureTracerOn
        assert(organicOn&&microbeOn&&sourceHourly ...
            &&isstruct(organic_management)&&isstruct(manureAomBefore), ...
            'soil:nitrogen:manureTracerContext', ...
            'Manure tracing requires exact management, AOM and source-hour ledgers.');
        context=struct('source_kind','manure_organic', ...
            'accepted',Y.soil_manure_n_tracer.organic, ...
            'management',organic_management, ...
            'aom_before',manureAomBefore, ...
            'aom_after',aomCandidate.state.aom.manure, ...
            'manure_routes',aomCandidate.material_routes.manure, ...
            'partition_pon_N_g_cm3',manurePartitionDonor, ...
            'aom_partition',manurePartition);
        % Change: 2026-11-02 Wenjie Liu -- stage all mutually exclusive
        % sources before sharing the same frozen mineral-N storage capacity.
        contexts={};
        if tracerOn,contexts={struct('source_kind','mineral')};end
        for cohort={'organic','tan'}
            key=cohort{1};context.source_kind=['manure_' key];
            context.accepted=Y.soil_manure_n_tracer.(key);
            contexts{end+1}=context; %#ok<AGROW>
        end
        result=soil.nitrogen.tracer.advanceCandidate( ...
            built,Y,fertilizer,soil_variables,nitrogen_variables,diagnostics, ...
            rootReservation,rootRemaining,tracerTotalBeforeReaction, ...
            tracerTotalAfterReaction,organicCandidate,domCandidate, ...
            microbeCandidate,organicExchange.mend.source_hour_ledger, ...
            transportTheta,transportQ,transportDt,dt_s,particleDrivers,contexts);
        if tracerOn,nitrogen.fertilizer_tracer=result{1};end
        candidates=struct('organic',result{end-1},'tan',result{end});
        nitrogen.manure_tracer=struct( ...
            'state',struct('schema_version',1,'organic',candidates.organic.state, ...
                'tan',candidates.tan.state), ...
            'diagnostics',struct('organic',candidates.organic.diagnostics, ...
                'tan',candidates.tan.diagnostics));
        tags={candidates.organic.state,candidates.tan.state};
        if tracerOn,tags{end+1}=nitrogen.fertilizer_tracer.state;end
        total=struct('nh4',nitrogen.trial.nh4,'no3',nitrogen.trial.no3, ...
            'urea',nitrogen.trial.urea,'nh4_granule',nitrogen.trial.C_NH4_granule, ...
            'nh4_ice',nitrogen.ice.nh4,'no3_ice',nitrogen.ice.no3, ...
            'organic_n',organicCandidate,'dom',domCandidate.state, ...
            'aom',organicCandidate.aom);
        soil.nitrogen.tracer.assertCombinedBounds(tags,total);
    end
end

function bundle=tracerMendBundle(organic,dom,domParameters)
% Package sole total-N owners for passive accepted-gross replay.
% Change: 2026-10-28 Wenjie Liu -- retain the dynamic DON phase subpool.
    bundle=struct('organic_n',organic, ...
        'don_total_g_cm3',dom.don_total_g_cm3, ...
        'sorption_model',domParameters.sorption_model);
    if strcmp(domParameters.sorption_model,'mend_dynamic')
        bundle.don_sorbed_g_cm3=dom.don_sorbed_g_cm3;
    end
end

function gross=fertilizerTracerTotalGross(external,mend,bdBottom,dt_s)
% Recover simultaneous accepted gross organic release and uptake while the
% physical ADE continues to consume only their signed net exchange.
    bdTop=flip(bdBottom(:));
    uptake4Bottom=flip(mend.nh4_actual_mg_g(:).*bdTop.*1000./dt_s);
    uptake3Bottom=flip(mend.no3_actual_mg_g(:).*bdTop.*1000./dt_s);
    source4Bottom=external.organic_nh4_rate(:)+uptake4Bottom;
    source3Bottom=external.organic_no3_rate(:)+uptake3Bottom;
    scale=max(1,max(abs([source4Bottom;source3Bottom;uptake4Bottom;uptake3Bottom])));
    tolerance=1024*eps(scale)*scale;
    assert(all(source4Bottom>=-tolerance)&&all(source3Bottom>=-tolerance), ...
        'advanceAcceptedCandidate:tracerOrganicGross', ...
        'Accepted organic gross source reconstructed as negative.');
    source4Bottom=max(source4Bottom,0);source3Bottom=max(source3Bottom,0);
    gross=struct( ...
        'nh4',struct( ...
            'sources',struct('organic_release',repmat(source4Bottom,1,2)), ...
            'sinks',struct('organic_immobilization',repmat(uptake4Bottom,1,2))), ...
        'no3',struct( ...
            'sources',struct('organic_release',repmat(source3Bottom,1,2)), ...
            'sinks',struct('organic_immobilization',repmat(uptake3Bottom,1,2))));
end

function [reservation,remaining]=rootDinReservation(request,dt_s,nh4Node,no3Node,theta,rhokg,bd)
%ROOTDINRESERVATION Reserve root NH4/NO3 uptake before organic C/N allocation.
%   Request sinks are [NL x nD] endpoint concentrations [ug cm^-3 s^-1].
%   Outputs are top-down [NL x 1] integrated inventories [mg N g^-1].
%   The reservation is a shadow inventory only; solveNitrogen applies the
%   accepted root sink exactly once later in the transaction.
%   Editor: Jizixing Zhu; Date: 2026-09-12.
    [availableNh4,availableNo3]=soil.organic_n.availableDinMgG( ...
        nh4Node,no3Node,theta,rhokg,bd);
    z=zeros(size(availableNh4));
    if ~request.enabled
        reservation=struct('nh4_mg_g',z,'no3_mg_g',z,'dt_s',dt_s, ...
            'source','disabled');
        remaining=struct('nh4_mg_g',availableNh4,'no3_mg_g',availableNo3);
        return
    end
    assert(isequal(size(request.sink_nh4),[numel(bd),2]) ...
        &&isequal(size(request.sink_no3),[numel(bd),2]), ...
        'soil:nitrogen:rootReservationShape','Root reservation requires NL-by-2 endpoint sinks.');
    rNh4=flip(mean(request.sink_nh4,2).*dt_s./bd)*1e-3;
    rNo3=flip(mean(request.sink_no3,2).*dt_s./bd)*1e-3;
    assert(all(isfinite([rNh4;rNo3]))&&all([rNh4;rNo3]>=0), ...
        'soil:nitrogen:rootReservationValues','Root reservation must be finite and nonnegative.');
    assert(all(rNh4<=availableNh4+1e-12) ...
        &&all(rNo3<=availableNo3+1e-12), ...
        'soil:nitrogen:rootReservationOverdraw','Root uptake exceeds available DIN inventory.');
    reservation=struct('nh4_mg_g',rNh4,'no3_mg_g',rNo3,'dt_s',dt_s, ...
        'source','COMPASS_root_priority_shadow');
    remaining=struct('nh4_mg_g',availableNh4-rNh4,'no3_mg_g',availableNo3-rNo3);
end

function budget=protectedNitrificationBudget(built,Y,state,root,inputs,dt)
% Fixed NN budget [ug N/cm2/s], bottom-up, before final accepted O2 replay.
% Keep the signed stock-minus-protection difference until positive same-hour
% sources have filled any deficit. Empirical mineralization/urea positives
% are conservatively omitted; their new inventory is usable next hour.
% No phase state, water driver or accepted inventory is modified.
% Editor: Wenjie Liu; Date: 2026-11-01.
    p=built.parameters.soil_n;dz=built.domain.dz(:);n=numel(dz);
    lump=@(v)[v(:,1).*dz/2;0]+[0;v(:,2).*dz/2];
    load=@(v)[(2*v(:,1)+v(:,2)).*dz/6;0]+[0;(v(:,1)+2*v(:,2)).*dz/6];
    stock=lump(state.theta_liquid_prev+repmat(p.RHOKG(:),1,2)).*Y.nh4(:);
    protected=lump(repmat(p.M_min_nit_ug_cm3(:),1,2));
    other=zeros(n,2);
    if isfield(inputs,'organic_enabled')&&inputs.organic_enabled
        other=repmat(inputs.organic_nh4_rate(:),1,2);
    end
    if root.enabled,other=other-root.sink_nh4;end
    budget=max((stock-protected)/dt+load(other),0);
    assert(all(isfinite(budget))&&all(budget>=0),'soil:nitrogen:protectionBudget');
end

function cfg = nitrogenCfg(built)
% Rebuild the small legacy view required by retained N kernels.
% Change: 2026-11-03 Wenjie Liu -- route the same effective transport bits
% used by the restart signature to the legacy soil nitrogen ADE kernels.
    cfg = init.makeLegacyRuntimeCfg(built);
    flags=soil.nitrogen.transport.numericalOptions(built.control);
    cfg.lUpW=flags.lUpW;
    cfg.lArtD=flags.lArtD;
    cfg.Grid = struct('NL', built.domain.NL, 'NN', built.domain.NN, ...
        'nD', built.domain.grid.nD, 'nS', built.domain.grid.nS, ...
        'DeltZ', built.domain.dz);
    for name = {'lNtotal','lThmrlefc','lBNF'}
        if ~isfield(cfg, name{1}), cfg.(name{1}) = 0; end
    end
end

function sv = soilVariables(built, Y, state)
% Assemble only the physical fields read by solveNitrogen and its kernels.
    tr = built.workspace.step.trial;
    d = built.workspace.derived;
    sv = struct( ...
        'hh', tr.h(:), 'h', Y.h(:), 'TT', tr.T(:), ...
        'P_g', tr.p_gas(:), 'Theta_LL', d.Theta_LL, ...
        'Theta_L', state.theta_liquid_prev, 'Theta_II', d.Theta_II, ...
        'Theta_g', d.Theta_g, 'POR', built.parameters.soil.porosity(:), ...
        'Theta_s', built.parameters.soil.theta_s(:), ...
        'Srt', built.workspace.flux.Srt, ...
        'C_NH4', Y.nh4(:), 'C_NH4_new', tr.nh4(:), ...
        'C_NO3', Y.no3(:), 'C_NO3_new', tr.no3(:), ...
        'C_Urea', Y.urea(:), 'C_Urea_new', tr.urea(:), ...
        'C_NH4_granule', tr.C_NH4_granule(:), ...
        'C_NH4_granule_old', Y.C_NH4_granule(:), ...
        'C_NH4_ice', Y.C_NH4_ice, 'C_NO3_ice', Y.C_NO3_ice, ...
        'Theta_LL_for_N_prev', state.theta_liquid_prev, ...
        'Theta_II_for_N_prev', state.theta_ice_prev);

end

function nv = nitrogenVariables(built, state, root_request, kt)
% Static parameters plus dynamic accepted N memory and this-step root sink.
    p = built.parameters.soil_n;
    nv = p;
    dyn = state.variables;
    names = fieldnames(dyn);
    for i = 1:numel(names)
        nv.(names{i}) = dyn.(names{i});
    end
    NL = built.domain.NL;
    nD = built.domain.grid.nD;
    nS = built.domain.grid.nS;
    assert(isequal(size(root_request.sink_nh4), [NL, nD]) ...
        && isequal(size(root_request.sink_no3), [NL, nD]), ...
        'soil:nitrogen:rootSinkShape', ...
        'Root N request must be [NL x nD] for both NH4 and NO3.');
    nv.SinkS = zeros(NL, nD, nS);
    if logical(root_request.enabled)
        nv.SinkS(:, :, 1) = root_request.sink_nh4;
        nv.SinkS(:, :, 2) = root_request.sink_no3;
    else
        assert(all(root_request.sink_nh4(:) == 0) ...
            && all(root_request.sink_no3(:) == 0), ...
            'soil:nitrogen:disabledRootSinkNonzero', ...
            'A disabled root-N request must be exactly zero for both species.');
    end
    nv.eT = finiteOr(built.workspace.derived.eT, nv.eT);
    nv.eTheta = finiteOr(built.workspace.derived.eTheta, nv.eTheta);
    nv.ePH = finiteOr(built.workspace.derived.ePH, nv.ePH);
    nv.Tm = twoSlot(forcingAt(built.forcing.hourly.Ta, kt, 0));
    nv.u_2 = twoSlot(forcingAt(built.forcing.hourly.U, kt, 0));
    nv.C_O2 = built.Y.soil_gas.C_O2(:);
end

function flux = fluxWork(previous)
% Retain only t^n and t^(n+1); the legacy kernels use positional indexing.
    flux = previous;
    names = fieldnames(flux);
    for i = 1:numel(names)
        name = names{i};
        value = flux.(name);
        if isnumeric(value) && isvector(value) && numel(value) == 2
            flux.(name) = [value(2); value(2)];
        end
    end
end

function flux = seedInitialStorage(flux, built, Y, state, nv)
% The first candidate has no previous time-series row.  Seed the accepted
% t^n storage explicitly so diagnostic and cumulative balances begin at the
% physical initial inventory rather than an artificial zero vector.
    if flux.cvnh4(1) ~= 0 || flux.cvno3(1) ~= 0
        return;
    end
    granule = Y.C_NH4_granule;
    [nh4, no3] = soil.nitrogen.transport.calculateStorage( ...
        state.theta_liquid_prev, Y.nh4, Y.no3, nv.RHOKG, ...
        built.domain.dz, built.domain.NL, state.theta_ice_prev, ...
        Y.C_NH4_ice, Y.C_NO3_ice, granule);
    flux.cvnh4(:) = nh4;
    flux.cvno3(:) = no3;
    flux.cvnh4OLD(:) = nh4;
    flux.cvno3OLD(:) = no3;
end

function boundary = nitrogenBoundary( ...
        built, advance, surface_request, surface_water, kt)
% Current-step surface and groundwater drivers in the retained N convention.
    h = built.forcing.hourly;
    b = built.boundary;
    water = advance.water_flux_accum;
    topHead = 0;
    nbch = 2;
    if isfield(advance, 'surface_boundary_summary') && ...
            advance.surface_boundary_summary.has_data
        topHead = max(advance.surface_boundary_summary.final_psi_top_cm, 0);
        if strcmp(advance.surface_boundary_summary.final_mode, 'finite_ponding')
            nbch = 1;
        end
    end
    boundary = struct( ...
        'NBCN', b.NBC_NH4_top, 'NBCNB', b.NBC_NH4_bottom, ...
        'cT_dep_flux', twoSlot(forcingAt(h.NH4_dep, kt, b.BC_NH4_top)), ...
        'cT1_dep_flux', twoSlot(forcingAt(h.NO3_dep, kt, b.BC_NO3_top)), ...
        'cB_NH4_conc', twoSlot(forcingAt(h.NH4_bot, kt, b.BC_NH4_bottom)), ...
        'cB_NO3_conc', twoSlot(forcingAt(h.NO3_bot, kt, b.BC_NO3_bottom)), ...
        'Precip_cms', twoSlot(forcingAt(h.Precip, kt, 0)), ...
        'DSTMAX', built.workspace.surface.DSTMAX, ...
        'DSTOR0', surface_water.DSTOR0, 'h_target', topHead, ...
        'NBCh', nbch, 'Infiltration', twoSlot(water.Infiltration), ...
        'EVAP', twoSlot(surface_request.EVAP), 'QMB', twoSlot(water.QMB), ...
        'NoTime', twoSlot(kt), 'BNF', 0);
    if isfield(built.control.switches, 'lBNF') && built.control.switches.lBNF == 1
        boundary.BNF = soil.nitrogen.calculateBiologicalNFixation( ...
            built.domain.dz, built.domain.NL);
    end
end

function [air,driver] = airVariables(soilAdvance, NL)
% Use the accepted-dt-weighted outer liquid-flux decomposition for N ADE.
% A legacy final-tile fallback is retained only for direct adapters/tests that
% predate the typed scheduler result. The mean is an outer-step approximation;
% unlike DOM/BIOPARTICLE, inorganic N is not yet replayed per accepted tile.
    diagnostics=soilAdvance.solver_diagnostics;
    selection='legacy_final_accepted_tile';
    if isfield(soilAdvance,'liquid_flux_outer_mean') ...
            &&isstruct(soilAdvance.liquid_flux_outer_mean) ...
            &&~isempty(fieldnames(soilAdvance.liquid_flux_outer_mean))
        diagnostics=soilAdvance.liquid_flux_outer_mean;
        selection='accepted_tile_dt_weighted_outer_mean';
    end
    if strcmp(selection,'legacy_final_accepted_tile')
        diagnostics.QL_abs=abs(diagnostics.QL);
    end
    required = {'QL','QL_abs','QL_h','QL_T','QL_a'};
    for i = 1:numel(required)
        assert(isfield(diagnostics, required{i}) ...
            && numel(diagnostics.(required{i})) == NL ...
            && all(isfinite(diagnostics.(required{i})(:))), ...
            'soil:nitrogen:missingLiquidFlux', ...
            ['Accepted N transport requires coupled-soil %s diagnostics; ', ...
             'do not substitute zero flux.'], required{i});
    end
    air = struct('QL', diagnostics.QL, 'QL_abs',diagnostics.QL_abs, ...
        'QL_h', diagnostics.QL_h, ...
        'QL_T', diagnostics.QL_T, 'QL_A', diagnostics.QL_a);
    driver=struct('schema_version',1,'selection',selection, ...
        'units','cm s-1','spatial_order','bottom_to_surface_elements', ...
        'QL',air.QL(:),'QL_abs',air.QL_abs(:),'QL_h',air.QL_h(:),'QL_T',air.QL_T(:), ...
        'QL_a',air.QL_A(:));
end

function [theta,q,dt]=acceptedTransportDrivers(soilAdvance,thetaFinal,dtOuter,nl)
% Select chronological QL/theta tiles, with a legacy one-tile fallback.
    diagnostics=soilAdvance.solver_diagnostics;
    if isfield(diagnostics,'accepted_transport_schedule')
        schedule=diagnostics.accepted_transport_schedule;
        required={'schema_version','tile_count','elapsed_s','dt_s', ...
            'theta_liquid_endpoint','QL'};
        for i=1:numel(required)
            assert(isfield(schedule,required{i}), ...
                'soil:nitrogen:transportScheduleSchema', ...
                'Accepted transport schedule lacks %s.',required{i});
        end
        dt=schedule.dt_s(:).';theta=schedule.theta_liquid_endpoint;q=schedule.QL;
        assert(schedule.schema_version==1&&schedule.tile_count==numel(dt) ...
            &&abs(schedule.elapsed_s-sum(dt))<=1e-9*max(1,sum(dt)) ...
            &&~isempty(dt)&&size(theta,1)==nl&&size(q,1)==nl ...
            &&size(theta,2)==numel(dt)&&size(q,2)==numel(dt) ...
            &&all(isfinite(theta(:)))&&all(theta(:)>=0)&&all(isfinite(q(:))) ...
            &&all(isfinite(dt))&&all(dt>0) ...
            &&abs(sum(dt)-dtOuter)<=1e-9*max(1,dtOuter), ...
            'soil:nitrogen:transportScheduleDriver', ...
            ['Accepted transport schedule must be finite, NL-by-tile, and ', ...
             'cover the exact outer-step duration.']);
        return
    end
    assert(isfield(diagnostics,'QL')&&numel(diagnostics.QL)==nl ...
        &&all(isfinite(diagnostics.QL(:))), ...
        'soil:nitrogen:missingLiquidFlux', ...
        'Accepted transport requires coupled-soil QL diagnostics.');
    theta=thetaFinal(:);q=diagnostics.QL(:);dt=dtOuter;
end

function state = stateFromWork(previous, nv, flux, sv)
% Scatter only accepted dynamic memory back to its soil-owned state.
    state = previous;
    dynamic = fieldnames(previous.variables);
    for i = 1:numel(dynamic)
        name = dynamic{i};
        if isfield(nv, name)
            state.variables.(name) = nv.(name);
        end
    end
    state.theta_liquid_prev = sv.Theta_LL_for_N_prev;
    state.theta_ice_prev = sv.Theta_II_for_N_prev;
    % Each retained vector is a two-time-level adapter, not a growing global
    % history.  Its accepted endpoint becomes both slots for the next call.
    flux = normalizeTwoSlotFlux(flux);
    state.flux = flux;
    state.schema_version = 1;
end

function flux = normalizeTwoSlotFlux(flux)
    names = fieldnames(flux);
    for i = 1:numel(names)
        name = names{i};
        value = flux.(name);
        if isnumeric(value) && isvector(value) && numel(value) == 2
            flux.(name) = [value(2); value(2)];
        end
    end
end

function value = forcingAt(series, kt, fallback)
    value = fallback;
    if ~isempty(series)
        value = series(min(kt, numel(series)));
    end
    if ~(isscalar(value) && isfinite(value))
        value = fallback;
    end
end

function value = twoSlot(value)
    value = [value; value];
end

function value = finiteOr(candidate, fallback)
    value = candidate;
    if isempty(value) || any(~isfinite(value(:)))
        value = fallback;
    end
end

function variables = updateTotalNitrogen( ...
        variables, C_NH4, C_NO3, thetaLiquid, parameters, domain)
%UPDATETOTALNITROGEN  Materialize accepted node total-N mass concentrations.
%
%   C_NH4/C_NO3 [NN x 1] are liquid concentrations [ug cm^-3 water].
%   thetaLiquid [NL x 2] is effective liquid content [-]. RHOKG [-] adds the
%   adsorbed NH4 buffer. Stored totals are [ug cm^-3 soil] and use two equal
%   accepted slots; diagnostics divide by dry bulk density to obtain [ug g^-1].
%
%   The endpoint mapping matches the retained bottom-up finite-element state:
%   node 1:NL use each element's first endpoint and node NN uses the last
%   element's second endpoint.
%
%   Reference: legacy accepted-state update in src/COMPASS.m:492-503.
% Editor: Jizixing Zhu; Date: 2026-07-20.
    NL = domain.NL;
    NN = domain.NN;
    assert(isequal(size(thetaLiquid), [NL, 2]) ...
            && numel(C_NH4) == NN && numel(C_NO3) == NN ...
            && numel(parameters.RHOKG) == NL, ...
        'soil:nitrogen:badTotalNitrogenShape', ...
        'Total-N materialization requires NN pools and NL-by-2 liquid water.');
    thetaNode = [thetaLiquid(:, 1); thetaLiquid(end, 2)];
    rhokg = parameters.RHOKG(:);
    rhokgNode = [rhokg; rhokg(end)];
    nh4Total = C_NH4(:) .* (thetaNode + rhokgNode);
    no3Total = C_NO3(:) .* thetaNode;
    assert(all(isfinite(nh4Total)) && all(nh4Total >= 0) ...
            && all(isfinite(no3Total)) && all(no3Total >= 0), ...
        'soil:nitrogen:badTotalNitrogen', ...
        'Accepted total-N node concentrations must be finite and nonnegative.');
    variables.NH4_total = repmat(nh4Total, 1, 2);
    variables.NO3_total = repmat(no3Total, 1, 2);
end
