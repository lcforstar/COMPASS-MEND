function nitrogen=advanceAcceptedCandidate(varargin)
%ADVANCEACCEPTEDCANDIDATE Rebuild detached C/N proposals until oxygen-feasible.
%   Inputs retain the core adapter contract below. Trials always read the
%   identical accepted owners; no transport, fertilizer or ledger is committed
%   twice. A per-layer common multiplier preserves MEND reaction stoichiometry.
% Editor: Jizixing Zhu; Date: 2026-09-14.
    assert(numel(varargin)>=7&&numel(varargin)<=12);
    args=varargin;
    if numel(args)<12,args(end+1:12)=repmat({[]},1,12-numel(args));end
    factor=ones(args{1}.domain.NL,1);
    for iteration=1:24
        nitrogen=advanceAtOxygenScale(args{:},factor);
        if ~isfield(nitrogen,'oxygen_retry_scale')
            nitrogen.diagnostics.oxygen_reaction_factor_top_down=factor;
            nitrogen.diagnostics.oxygen_coupling_iterations=iteration;
            return
        end
        factor=factor.*nitrogen.oxygen_retry_scale;
        factor(factor<1e-12)=0;
    end
    error('soil:nitrogen:oxygenCoupling', ...
        'Could not rebuild an oxygen-feasible C/N candidate before commit.');
end

function nitrogen = advanceAtOxygenScale( ...
        built, root_request, soil_advance, surface_water_request, fertilizer, ...
        kt, dt_s, surface_water_candidate, organic_accepted, ...
        dom_accepted, microbe_accepted, organic_management, oxygenScale)
%ADVANCEACCEPTEDCANDIDATE  Form one detached accepted soil-N endpoint.
%
%   This adapter preserves the existing reaction, freeze-partition, O2, and
%   ADE kernels.  It supplies their legacy arrays as two time levels only:
%   slot 1 is accepted t^n and slot 2 is the candidate t^(n+1).  Therefore N
%   history has explicit soil ownership in Y.soil_nitrogen rather than being
%   hidden in full-run diagnostic vectors.  No field of BUILT is mutated.
%
%   The soil water/heat/gas solve is already accepted when this routine runs.
%   Its final, converged liquid-flux decomposition is the N transport driver;
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
    soil_variables = soilVariables(built, Y, state);
    nitrogen_variables = nitrogenVariables(built, state, root_request, kt);
    nitrogen_flux = fluxWork(state.flux);
    nitrogen_flux = seedInitialStorage( ...
        nitrogen_flux, built, Y, state, nitrogen_variables);
    nitrogen_boundary = nitrogenBoundary( ...
        built, soil_advance, surface_water_request, ...
        surface_water_candidate, kt);
    air_variables = airVariables(soil_advance.solver_diagnostics, built.domain.NL);
    gas_variables = gasVariables(built, kt);

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
    aomCandidate=[];
    domCandidate=[]; microbeCandidate=[];
    domOn=isfield(built.control.switches,'lSoilDOM') ...
        &&built.control.switches.lSoilDOM==1;
    microbeOn=isfield(built.control.switches,'lExplicitSoilMicrobes') ...
        &&built.control.switches.lExplicitSoilMicrobes==1;
    microbeTransportOn=isfield(built.control.switches,'lMicrobialTransport') ...
        &&built.control.switches.lMicrobialTransport==1;
    thetaElement=mean(built.workspace.derived.Theta_LL,2);
    organicOn=isfield(built.control.switches,'lSoilOrganic') ...
        &&built.control.switches.lSoilOrganic==1;
    if organicOn
        if nargin<9||isempty(organic_accepted)
            organic_accepted=built.Y.soil_organic;
        end
        trial=built.workspace.step.trial;
        temperatureElement=0.5.*(trial.T(1:end-1)+trial.T(2:end));
        headElement=0.5.*(trial.h(1:end-1)+trial.h(2:end));
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
        particulateFraction=built.parameters.soil_microbes.input_fraction(1:2);
        particulateFraction=particulateFraction./sum(particulateFraction);
        addedPoc=flip(aomCandidate.poc_C_g_cm3)./bdTop.*1000;
        microbe_accepted.POM1=microbe_accepted.POM1 ...
            +particulateFraction(1).*addedPoc;
        microbe_accepted.POM2=microbe_accepted.POM2 ...
            +particulateFraction(2).*addedPoc;
        microbe_accepted.POC=microbe_accepted.POM1+microbe_accepted.POM2;
        organicCandidate.PON1=organicCandidate.PON1 ...
            +particulateFraction(1).*aomCandidate.pon_N_g_cm3;
        organicCandidate.PON2=organicCandidate.PON2 ...
            +particulateFraction(2).*aomCandidate.pon_N_g_cm3;
        organicCandidate.PON=organicCandidate.PON1+organicCandidate.PON2;
        organic_accepted=organicCandidate;
    end
    if domOn
        if nargin<10||isempty(dom_accepted), dom_accepted=built.Y.soil_dom; end
        q=soil_advance.solver_diagnostics.QL(:);
        domCandidate=soil.dom.advanceAcceptedCandidate(dom_accepted, ...
            thetaElement,built.parameters.soil.bulk_density(:), ...
            built.domain.dz(:),q,built.parameters.soil_dom,dt_s);
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
            dinAfterRoot=soil.nitrogen.mendKineticDinAfterRoot( ...
                rootRemaining,rootReservation,Y.nh4,Y.no3,state.theta_liquid_prev, ...
                built.parameters.soil_n.RHOKG,built.parameters.soil.bulk_density(:));
            [microbeCandidate,organicCandidate,domCandidate.state,mendN,nFeedback]= ...
                soil.microbe.advanceSourceCandidate(microbe_accepted,organic_accepted, ...
                domBeforeMicrobe,thetaElement,built.parameters.soil.bulk_density(:), ...
                built.parameters.soil_dom,microbeParameters,dt_s, ...
                temperatureElement,headElement,dinAfterRoot,rootRemaining,oxygenScale);
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
            particle=soil.microbe.bioparticleTransport( ...
                microbeCandidate.state,organicCandidate,thetaElement, ...
                built.parameters.soil.bulk_density(:),built.domain.dz(:), ...
                q,built.parameters.soil_microbes.bioparticle,dt_s);
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
        demand=soil.gas.calcO2Sink(zeros(built.domain.NL,2), ...
            repmat(externalInputs.organic_co2_c_rate(:),1,2),built.domain.NL);
        [~,~,~,~,oxygenProposal]=soil.gas.solveO2Transport( ...
            nitrogen_variables.C_O2,soil_variables.Theta_g,soil_variables.POR, ...
            soil_variables.TT,gas_variables.V_A,demand,built.domain.dz,dt_s, ...
            built.domain.NL,built.domain.NN,gas_variables.O2a_molcm3(2),[], ...
            soil_variables.oxygen_phase);
        if any(oxygenProposal.additional_scale<1-1e-10)
            retry=ones(built.domain.NL,1);
            limited=oxygenProposal.additional_scale<1-1e-10;
            retry(limited)=0.8*oxygenProposal.additional_scale(limited);
            nitrogen=struct('oxygen_retry_scale', ...
                flip(retry));
            return
        end
    elseif organicOn
        error('advanceAcceptedCandidate:soilCarbonWithoutMemc', ...
            'The native soil-carbon path requires DOM and MEND.');
    end
    externalInputs.root_reservation_mg_g = [rootReservation.nh4_mg_g, ...
        rootReservation.no3_mg_g];
    if microbeOn
        externalInputs.pre_reaction_reservation_mg_g=...
            externalInputs.root_reservation_mg_g+...
            [max(mendAllocation.net_nh4_sink_mg_g,0), ...
             mendAllocation.net_no3_sink_mg_g];
    else
        externalInputs.pre_reaction_reservation_mg_g=...
            externalInputs.root_reservation_mg_g;
    end
    externalInputs.root_remaining_mg_g = [rootRemaining.nh4_mg_g, ...
        rootRemaining.no3_mg_g];

    [nitrogen_variables, soil_variables, nitrogen_flux, diagnostics] = ...
        soil.nitrogen.solveNitrogen(nitrogen_variables, nitrogen_flux, ...
            nitrogen_boundary, soil_variables, air_variables, gas_variables, ...
            cfg, 2, dt_s, externalInputs);
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
    nitrogen.oxygen = nitrogen_variables.C_O2;
    nitrogen.diagnostics = diagnostics;
    nitrogen.boundary = nitrogen_boundary;
    nitrogen.root_request = root_request;
    nitrogen.root_reservation = rootReservation;
    nitrogen.root_remaining = rootRemaining;
    sinkActual = nitrogen_variables.SinkS;
    nitrogen.root_actual = struct('sink_nh4',sinkActual(:,:,1), ...
        'sink_no3',sinkActual(:,:,2), ...
        'U_N_total',sum(sum(mean(sinkActual,2),1).*built.domain.dz(:))*0.01);
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

function cfg = nitrogenCfg(built)
% Rebuild the small legacy view required by retained N kernels.
    cfg = init.makeLegacyRuntimeCfg(built);
    cfg.Grid = struct('NL', built.domain.NL, 'NN', built.domain.NN, ...
        'nD', built.domain.grid.nD, 'nS', built.domain.grid.nS, ...
        'DeltZ', built.domain.dz);
    for name = {'lNtotal','lArtD','lUpW','lThmrlefc','lBNF'}
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
    assert(isfield(Y.soil_gas,'oxygen_phase'), ...
        'soil:gas:missingOxygenHistory','O2 requires phase history; rebuild the old checkpoint.');
    old=Y.soil_gas.oxygen_phase;
    sv.oxygen_phase=struct('theta_liquid_new',d.Theta_LL, ...
        'theta_liquid_old',old.theta_liquid,'theta_gas_old',old.theta_gas, ...
        'temperature_old',old.temperature_C);
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

function air = airVariables(diagnostics, NL)
% The final accepted coupled soil tile owns the transport decomposition.
    required = {'QL','QL_h','QL_T','QL_a'};
    for i = 1:numel(required)
        assert(isfield(diagnostics, required{i}) ...
            && numel(diagnostics.(required{i})) == NL ...
            && all(isfinite(diagnostics.(required{i})(:))), ...
            'soil:nitrogen:missingLiquidFlux', ...
            ['Accepted N transport requires final coupled-soil %s diagnostics; ', ...
             'do not substitute zero flux.'], required{i});
    end
    air = struct('QL', diagnostics.QL, 'QL_h', diagnostics.QL_h, ...
        'QL_T', diagnostics.QL_T, 'QL_A', diagnostics.QL_a);
end

function gas = gasVariables(built, kt)
% O2 boundary is an atmospheric forcing; V_A is a final soil transport field.
    gas = struct('V_A', built.workspace.flux.V_A, ...
        'O2a_molcm3', twoSlot(forcingAt( ...
            built.forcing.gas.O2a_molcm3, kt, NaN)));
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
