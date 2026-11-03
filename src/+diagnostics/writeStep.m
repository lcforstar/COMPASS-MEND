function writeStep(model, kt, report)
%WRITESTEP  Write accepted or failed outer-step diagnostics to the registry.
%
%   Inputs:
%     model  [in] Model shell containing runtime registry and accepted ledgers.
%     kt     [-] One-based outer-step index.
%     report [in] Fixed step report containing status and duration [s].
%
%   This output-only function never updates physical state and physical code
%   never reads values back from the registry.
%   Change: 2026-09-29 Wenjie Liu -- organic input/top-output diagnostics and
%     schema-4 terminal source-label routes; mineral-N diagnostics unchanged.
%   Editor: Jizixing Zhu; Date: 2026-07-16.
%   Change: 2026-09-29 Wenjie Liu -- publish DIN dry-mass observables from
%     accepted half-element capacities, without changing kinetic total-N slots.
%   Change: 2026-07-16 Jizixing Zhu -- migrate the accepted SCOPE, plant,
%     root-uptake, surface-water, and frozen-soil outputs from COMPASS.m.
%   Change: 2026-07-17 Jizixing Zhu -- export accepted finite-ponding
%     complementarity states, fluxes, and worst substep residuals.
%   Change: 2026-07-17 Jizixing Zhu -- publish one canonical rate/amount
%     boundary schema and the combined surface-soil water ledger.
%   Change: 2026-07-17 Jizixing Zhu -- write the accepted ClimaLand snow
%     state, rain/snow partition, transfers, and conservation residuals.
%   Change: 2026-07-17 Jizixing Zhu -- distinguish the SCOPE radiative skin
%     temperature from the soil heat boundary imposed under snow.
%   Change: 2026-07-17 Jizixing Zhu -- publish realized soil evaporation
%     separately from the upstream SCOPE evaporation demand at the dry limit.
%   Change: 2026-07-17 Jizixing Zhu -- publish conservative terminal snow meltout.
%   Change: 2026-07-18 Jizixing Zhu -- publish terminal meltwater enthalpy transfer.
%   Change: 2026-07-18 Jizixing Zhu -- publish the first accepted
%     surface-plus-soil energy ledger.
%   Change: 2026-07-18 Jizixing Zhu -- publish the before/after bottom-temperature
%     initialization mismatch without the discarded A1 correction fields.
%   Change: 2026-07-18 Jizixing Zhu -- publish positive ecological ET,
%     signed net vapor exchange, and separate dew/deposition pathways.
%   Change: 2026-07-19 Jizixing Zhu -- publish the generalized prognostic
%     surface-snow-soil total-energy ledger.
%   Change: 2026-07-19 Jizixing Zhu -- publish outer interface closure audit.
%   Change: 2026-07-20 Jizixing Zhu -- restore accepted soil-N registry writes.
%   Change: 2026-07-20 Jizixing Zhu — write accepted root sinks and convert
%     two-slot total-N state to finite node diagnostics with explicit units.
%   Change: 2026-07-20 Jizixing Zhu — materialize accepted vegetation C/N,
%     LUNA/Cab, turnover, and complete plant C/N closure diagnostics.
%   Change: 2026-08-29 Codex -- write accepted native organic pools,
%     mineral-N/CO2 exchange, conservation residuals, and organic fertilizer.
%   Change: 2026-08-16 Codex -- publish the continuous carbon and nitrogen
%     inputs used by the evidence-based dynamic spin-up.
%   Change: 2026-08-28 Codex -- publish native AOM, MEMC/MEND, DOM, paired-N,
%     and BIOPARTICLE pool diagnostics for layerwise closure inspection.
%   Change: 2026-09-12 Codex -- publish complete accepted MEND C/N pool,
%     internal-flux, root-first DIN competition, and residual audit ledgers.
%   Change: 2026-09-13 Jizixing Zhu -- publish accepted source-MEND direct
%     DOC/DON inputs from PlantCN litter.
%   Change: 2026-10-09 Wenjie Liu -- publish accepted DOC/DON soil-water
%     concentrations in mg/L while retaining all existing inventory outputs.
%   Change: 2026-10-12 Wenjie Liu -- expose MAOC=MOC+QOC and D+Q explicitly
%     as diagnostic aggregates; DOC remains aqueous D and no owner is added.
%   Change: 2026-11-01 Wenjie Liu -- publish NH4 protection reductions and
%     nodal rate budgets from accepted diagnostic caches only when present.
    runtime = model.runtime;
    if ~isfield(runtime, 'diagnostics') ...
            || ~isfield(runtime.diagnostics, 'registry')
        return;
    end
    registry = runtime.diagnostics.registry;
    registry.write('step_failed', double(report.step_failed), kt);
    if report.step_failed
        return;
    end
    diagnostics.writeScopeStep(registry, model.built, kt);
    diagnostics.writePlantStep(registry, model.built, kt, report.dt);
    diagnostics.writePlantNitrogenStep(registry, model.built, kt);
    diagnostics.writeSoilStateStep(registry, model.built, kt);
    writeSoilNitrogen(registry, model.built, kt, report);
    writeSoilOrganic(registry, model.built, kt, report);
    diagnostics.writeSolverInternalStep(registry, model.built, kt);
    writeSnow(registry, model.built, kt, report.dt);
    writeWaterBalance(registry, model.built, kt, report.dt);
    writeSoilEnergyBalance(registry, model.built, kt);
    writeSurfaceSoilEnergyBalance(registry, model.built, kt);
    writeSoilInterfaceClosure(registry, model.built, kt);
    writeDryAirMassBalance(registry, model.built, kt);
    writeSoilBoundary(registry, model.built, kt);
end

function writeSoilOrganic(registry, built, kt, report)
%WRITESOILORGANIC Export accepted native AOM--MEMC--DOM C/N state.
% Editor: Jizixing Zhu; Date: 2026-08-16.
% Change: 2026-08-28 Codex -- replace retired-pool diagnostics with native
%   three-AOM, MEMC carbon, DOM, and paired organic-N ledgers.
% Change: 2026-08-30 Codex -- publish manure TAN and MIMICS-CN closure terms.
% Change: 2026-09-01 Codex -- distinguish surface-AOM mineralized TAN from
%   event-supplied manure TAN and expose material-resolved TAN production.
% Change: 2026-09-05 Jizixing Zhu -- export accepted MEMC mortality routing.
% Change: 2026-09-12 Codex -- export accepted MEND microbial DIN competition
%   and nitrogen-limitation results from the committed transaction ledger.
% Change: 2026-09-13 Jizixing Zhu -- export source-MEND direct PlantCN
%   litter inputs to DOC and DON.
% Change: 2026-09-13 Codex -- export accepted microbial C/N phase pairs.
% Change: 2026-09-14 Jizixing Zhu -- read dynamic DOM phases from their owner
%   and include intracellular enzyme C/N in storage and profile diagnostics.
% Change: 2026-09-21 Wenjie Liu -- include accepted BIOPARTICLE transport
%   C/N residuals from their production candidate field exactly once.
    if ~isfield(built.control.switches,'lSoilOrganic') ...
            || built.control.switches.lSoilOrganic ~= 1 ...
            || ~isfield(built.Y,'soil_organic')
        return;
    end
    state=built.Y.soil_organic; n=built.domain.NL;
    dzBottom=built.domain.dz(:); dzTop=flip(dzBottom);
    bdTop=flip(built.parameters.soil.bulk_density(:));
    microbes=built.Y.soil_microbes;
    convert=bdTop./1000;
    poc=microbes.POC(:).*convert; moc=microbes.MOC(:).*convert;
    activeC=microbes.MB(:).*convert; dormantC=microbes.BD(:).*convert;
    smbC=activeC+dormantC;
    enzyme=(microbes.EP(:)+microbes.EM(:)).*convert;
    intracellularC=zeros(n,1);intracellularN=zeros(n,1);
    if isfield(microbes,'ENZN')
        intracellularC=microbes.ENZN(:).*convert;
        intracellularN=flip(state.ENZN_N(:));
        enzyme=enzyme+intracellularC;
    end
    docTotalBottom=built.Y.soil_dom.doc_total_g_cm3(:);
    donTotalBottom=built.Y.soil_dom.don_total_g_cm3(:);
    thetaBottom=mean(built.workspace.derived.Theta_LL,2);
    docPhase=soil.dom.phaseState(built.Y.soil_dom,'doc',thetaBottom, ...
        built.parameters.soil.bulk_density(:),built.parameters.soil_dom);
    donPhase=soil.dom.phaseState(built.Y.soil_dom,'don',thetaBottom, ...
        built.parameters.soil.bulk_density(:),built.parameters.soil_dom);
    doc=flip(thetaBottom.*docPhase.aqueous_g_cm3_water);
    qoc=flip(built.parameters.soil.bulk_density(:) ...
        .*docPhase.sorbed_g_g_soil);
    don=flip(thetaBottom.*donPhase.aqueous_g_cm3_water);
    qon=flip(built.parameters.soil.bulk_density(:) ...
        .*donPhase.sorbed_g_g_soil);
    docTotal=flip(docTotalBottom);
    aomC=zeros(n,1); aomN=zeros(n,1);
    aomProfiles=struct();
    names=fieldnames(state.aom);
    for j=1:numel(names)
        a=state.aom.(names{j});
        aomC=aomC+flip(a.element_C(:));
        aomN=aomN+flip(a.element_N(:));
        aomC(1)=aomC(1)+a.surface_C/dzTop(1);
        aomN(1)=aomN(1)+a.surface_N/dzTop(1);
        materialC=flip(a.element_C(:)); materialN=flip(a.element_N(:));
        materialC(1)=materialC(1)+a.surface_C/dzTop(1);
        materialN(1)=materialN(1)+a.surface_N/dzTop(1);
        aomProfiles.([names{j} '_C'])=materialC;
        aomProfiles.([names{j} '_N'])=materialN;
    end
    iom=flip(state.IOM_C(:)); nec=flip(state.NEC_C(:));
    assert(max(abs(doc+qoc-docTotal))<=1e-12*max(1,max(docTotal)), ...
        'diagnostics:writeStep:docPartitionClosure', ...
        'Equilibrium DOC/QOC diagnostic partition failed closure.');
    assert(max(abs(don+qon-flip(donTotalBottom))) ...
            <=1e-12*max(1,max(donTotalBottom)), ...
        'diagnostics:writeStep:donPartitionClosure', ...
        'Equilibrium DON/QON diagnostic partition failed closure.');
    somC=poc+moc+docTotal+enzyme+iom+nec;
    % Observation operator for dichromate-heated bulk-soil SOC. Unreacted
    % surface or incorporated AOM is excluded; its carbon remains in totalC
    % and therefore in the complete conservation ledger. DOC/QOC, microbial
    % biomass, necromass, enzymes and IOM are internal soil-carbon owners.
    socObservationC=smbC+somC;
    totalC=socObservationC+aomC;
    totalN=flip(state.PON+state.MON+state.QON+state.MB_N+state.BD_N ...
        +state.EP_N+state.EM_N+state.NEC_N ...
        +built.Y.soil_dom.don_total_g_cm3)+aomN+intracellularN;
    assert(all(isfinite(totalC)) && all(totalC>=0) ...
            && all(isfinite(totalN)) && all(totalN>=0), ...
        'diagnostics:writeStep:badSoilOrganicState', ...
        'Accepted soil-organic C/N profiles must be finite and nonnegative.');
    registry.write('organic_carbon_profile',totalC,kt);
    registry.write('soc_observation_profile',socObservationC,kt);
    registry.write('organic_smb_carbon_profile',smbC,kt);
    registry.write('organic_som_carbon_profile',somC,kt);
    registry.write('organic_aom_carbon_profile',aomC,kt);
    for material={'litter','straw','manure'}
        key=material{1};
        registry.write(['aom_' key '_carbon_profile'],aomProfiles.([key '_C']),kt);
        registry.write(['aom_' key '_nitrogen_profile'],aomProfiles.([key '_N']),kt);
    end
    registry.write('poc_carbon_profile',poc,kt);
    registry.write('pom1_carbon_profile',microbes.POM1(:).*convert,kt);
    registry.write('pom2_carbon_profile',microbes.POM2(:).*convert,kt);
    registry.write('moc_carbon_profile',moc,kt);
    registry.write('maoc_carbon_profile',moc+qoc,kt);
    registry.write('doc_carbon_profile',doc,kt);
    registry.write('doc_plus_qoc_carbon_profile',docTotal,kt);
    registry.write('doc_aqueous_concentration_profile',flip(docPhase.aqueous_mg_L),kt);
    registry.write('qoc_carbon_profile',qoc,kt);
    registry.write('microbe_active_carbon_profile',activeC,kt);
    registry.write('microbe_dormant_carbon_profile',dormantC,kt);
    registry.write('microbe_total_carbon_profile',smbC,kt);
    registry.write('enzyme_pom1_carbon_profile',microbes.ENZP1(:).*convert,kt);
    registry.write('enzyme_pom2_carbon_profile',microbes.ENZP2(:).*convert,kt);
    registry.write('enzyme_mom_carbon_profile',microbes.EM(:).*convert,kt);
    registry.write('microbial_necromass_carbon_profile',nec,kt);
    registry.write('organic_nitrogen_profile',totalN,kt);
    registry.write('pon1_nitrogen_profile',flip(state.PON1(:)),kt);
    registry.write('pon2_nitrogen_profile',flip(state.PON2(:)),kt);
    registry.write('mon_nitrogen_profile',flip(state.MON(:)),kt);
    registry.write('don_nitrogen_profile',don,kt);
    registry.write('don_aqueous_concentration_profile',flip(donPhase.aqueous_mg_L),kt);
    registry.write('qon_nitrogen_profile',qon,kt);
    registry.write('microbe_active_nitrogen_profile',flip(state.MB_N(:)),kt);
    registry.write('microbe_dormant_nitrogen_profile',flip(state.BD_N(:)),kt);
    registry.write('enzyme_pom1_nitrogen_profile',flip(state.ENZP1_N(:)),kt);
    registry.write('enzyme_pom2_nitrogen_profile',flip(state.ENZP2_N(:)),kt);
    registry.write('enzyme_mom_nitrogen_profile',flip(state.EM_N(:)),kt);
    registry.write('enzyme_intracellular_carbon_profile',intracellularC,kt);
    registry.write('enzyme_intracellular_nitrogen_profile',intracellularN,kt);
    registry.write('microbial_necromass_nitrogen_profile', ...
        flip(state.NEC_N(:)),kt);
    for phase={'mobile','attached','retained'}
        key=phase{1};
        registry.write(['microbe_active_' key '_carbon_profile'], ...
            microbes.(['MB_' key])(:).*convert,kt);
        registry.write(['microbe_dormant_' key '_carbon_profile'], ...
            microbes.(['BD_' key])(:).*convert,kt);
        registry.write(['microbe_' key '_carbon_profile'], ...
            (microbes.(['MB_' key])(:)+microbes.(['BD_' key])(:)).*convert,kt);
        if isfield(microbes,['ENZN_' key])
            registry.write(['enzyme_intracellular_' key '_carbon_profile'], ...
                microbes.(['ENZN_' key])(:).*convert,kt);
        else
            registry.write(['enzyme_intracellular_' key '_carbon_profile'],zeros(n,1),kt);
        end
        activePhaseN=state.(['MB_' key '_N']);
        dormantPhaseN=state.(['BD_' key '_N']);
        registry.write(['microbe_active_' key '_nitrogen_profile'], ...
            flip(activePhaseN(:)),kt);
        registry.write(['microbe_dormant_' key '_nitrogen_profile'], ...
            flip(dormantPhaseN(:)),kt);
    end
    registry.write('organic_carbon_storage',sum(totalC.*dzTop),kt);
    registry.write('organic_nitrogen_storage',sum(totalN.*dzTop),kt);
    registry.write('soc_observation_storage',sum(socObservationC.*dzTop),kt);
    registry.write('doc_carbon_storage',sum(doc.*dzTop),kt);
    registry.write('qoc_carbon_storage',sum(qoc.*dzTop),kt);
    registry.write('moc_carbon_storage',sum(moc.*dzTop),kt);
    registry.write('maoc_carbon_storage',sum((moc+qoc).*dzTop),kt);
    registry.write('doc_plus_qoc_carbon_storage',sum(docTotal.*dzTop),kt);
    if isfield(built.workspace,'soil_organic')
        work=built.workspace.soil_organic;
        x=work.exchange; d=work.diagnostics;
        aom=x.aom;
        registry.write('aom_respiration_co2_profile',flip(aom.co2_C_g_cm3(:)),kt);
        if isfield(built.workspace,'soil_dom')&&isfield(built.workspace.soil_dom,'transport')
            transport=built.workspace.soil_dom.transport;
            registry.write('doc_transport_net_change_profile', ...
                flip(transport.doc.net_layer_change_g_cm3(:)),kt);
            registry.write('don_transport_net_change_profile', ...
                flip(transport.don.net_layer_change_g_cm3(:)),kt);
            for key={'doc','don'}
                element=key{1}; ledger=transport.(element);
                registry.write([element '_transport_storage_before'], ...
                    ledger.storage_before_g_cm2,kt);
                registry.write([element '_transport_storage_after'], ...
                    ledger.storage_after_g_cm2,kt);
                registry.write([element '_transport_closure_residual'], ...
                    ledger.closure_residual_g_cm2,kt);
                % Change: 2026-09-29 Wenjie Liu -- accepted bottom exports,
                % distinct from respiration and internal C/N exchange.
                registry.write([element '_bottom_outflow_step'],ledger.bottom_outflow_g_cm2,kt);
                registry.write([element '_bottom_outflow_flux'], ...
                    ledger.bottom_outflow_g_cm2/ledger.driver_elapsed_s,kt);
                cumulative=0;
                if isfield(built.Y.soil_dom,'bottom_export_cumulative')
                    cumulative=built.Y.soil_dom.bottom_export_cumulative.([element '_g_cm2']);
                end
                registry.write([element '_bottom_outflow_cumulative'],cumulative,kt);
            end
        end
        registry.write('microbial_transport_enabled',double( ...
            built.parameters.soil_microbes.transport_enabled),kt);
        diagnostics.writeBioparticleExports(registry,built,kt);
        diagnostics.writeOrganicBoundaryExchange(registry,built,kt);
        registry.write('organic_abiotic_factor', ...
            flip(aom.temperature_scalar.*aom.water_scalar.*aom.depth_scalar),kt);
        aomCo2Rate=sum(aom.co2_C_rate_ug_cm3_s(:).*dzBottom).*1e-6;
        co2Rate=aomCo2Rate;
        registry.write('aom_respiration_co2_flux',aomCo2Rate,kt);
        if isfield(built.workspace,'soil_microbes')
            co2Rate=co2Rate+sum(built.workspace.soil_microbes. ...
                co2_c_rate_ug_cm3_s(:).*dzBottom).*1e-6;
        end
        registry.write('organic_co2_c_flux',co2Rate,kt);
        nh4Profile=zeros(n,1); no3Profile=zeros(n,1);
        tanProfile=zeros(n,1); aomTanProfile=zeros(n,1);
        if isfield(x,'nh4_ug_cm3_s')
            aomTanProfile=x.nh4_ug_cm3_s(:);
            nh4Profile=nh4Profile+aomTanProfile;
        end
        if isfield(x,'tan_nh4_ug_cm3_s')
            tanProfile=x.tan_nh4_ug_cm3_s(:);
            nh4Profile=nh4Profile+tanProfile;
        end
        if isfield(built.workspace,'soil_dom')
            nh4Profile=nh4Profile+built.workspace.soil_dom.nh4_source_ug_cm3_s(:);
        end
        if isfield(x,'mend')&&isfield(x.mend,'din_allocation')
            ma=x.mend.din_allocation;
            stepSeconds=x.mend.dt_s;
            mk=x.mend.din_kinetic_request;
            mendNh4Rate=flip(-ma.net_nh4_sink_mg_g.*bdTop.*1000./stepSeconds);
            mendNo3Rate=flip(-ma.net_no3_sink_mg_g.*bdTop.*1000./stepSeconds);
            nh4Profile=nh4Profile+mendNh4Rate;
            no3Profile=no3Profile+mendNo3Rate;
            nh4Request=ma.first.requested_mg_g(:,1)+ma.second.requested_mg_g(:,1);
            no3Request=ma.first.requested_mg_g(:,2)+ma.second.requested_mg_g(:,2);
            registry.write('mend_nh4_request_profile',nh4Request,kt);
            registry.write('mend_no3_request_profile',no3Request,kt);
            registry.write('mend_nh4_actual_profile',ma.nh4_actual_mg_g,kt);
            registry.write('mend_no3_actual_profile',ma.no3_actual_mg_g,kt);
            registry.write('mend_din_unmet_profile',ma.unmet_mg_g,kt);
            registry.write('mend_n_limitation_factor', ...
                x.mend.din_kinetic_request.phi_n_limitation,kt);
            registry.write('mend_stoichiometric_demand_profile', ...
                mk.stoichiometric_demand_mg_g,kt);
            registry.write('mend_kinetic_nh4_request_profile', ...
                mk.raw_kinetic_request_mg_g(:,1),kt);
            registry.write('mend_kinetic_no3_request_profile', ...
                mk.raw_kinetic_request_mg_g(:,2),kt);
            registry.write('mend_first_nh4_actual_profile', ...
                ma.first.nh4_actual_mg_g,kt);
            registry.write('mend_first_no3_actual_profile', ...
                ma.first.no3_actual_mg_g,kt);
            registry.write('mend_second_nh4_actual_profile', ...
                ma.second.nh4_actual_mg_g,kt);
            registry.write('mend_second_no3_actual_profile', ...
                ma.second.no3_actual_mg_g,kt);
            registry.write('mend_nh4_release_profile',ma.nh4_release_mg_g,kt);
            registry.write('mend_nh4_net_sink_profile',ma.net_nh4_sink_mg_g,kt);
            registry.write('mend_no3_net_sink_profile',ma.net_no3_sink_mg_g,kt);
            assert(isfield(built.workspace,'nitrogen') ...
                &&isfield(built.workspace.nitrogen,'root_reservation') ...
                &&isfield(built.workspace.nitrogen,'root_remaining'), ...
                'diagnostics:writeStep:missingDinCompetitionLedger', ...
                'Accepted MEND output requires root reservation and remaining DIN.');
            rootReserved=built.workspace.nitrogen.root_reservation;
            rootRemaining=built.workspace.nitrogen.root_remaining;
            registry.write('root_reserved_nh4_profile',rootReserved.nh4_mg_g,kt);
            registry.write('root_reserved_no3_profile',rootReserved.no3_mg_g,kt);
            registry.write('post_root_nh4_available_profile', ...
                rootRemaining.nh4_mg_g,kt);
            registry.write('post_root_no3_available_profile', ...
                rootRemaining.no3_mg_g,kt);
        end
        registry.write('organic_nh4_source_profile',flip(nh4Profile),kt);
        registry.write('organic_no3_source_profile',flip(no3Profile),kt);
        registry.write('organic_tan_nh4_source_profile',flip(tanProfile),kt);
        registry.write('aom_surface_tan_nh4_source_profile', ...
            flip(aomTanProfile),kt);
        nh4Rate=sum(nh4Profile.*dzBottom);
        registry.write('organic_nh4_exchange_flux',nh4Rate,kt);
        registry.write('organic_no3_exchange_flux',sum(no3Profile.*dzBottom),kt);
        carbonResidual=d.aom_carbon_residual_g_cm2;
        nitrogenResidual=d.aom_nitrogen_residual_g_cm2;
        mendCarbonResidual=0;
        mendNitrogenResidual=0;
        if isfield(built.workspace,'soil_microbes')
            mw=built.workspace.soil_microbes;
            if isfield(mw,'carbon_residual_mg_g')
                mendCarbonResidual=sum(mw.carbon_residual_mg_g(:).* ...
                    bdTop./1000.*dzTop);
                carbonResidual=carbonResidual+mendCarbonResidual;
                registry.write('mend_carbon_residual_profile', ...
                    mw.carbon_residual_mg_g,kt);
            end
            if isfield(d,'mend_nitrogen_residual_mg_g')
                mendNitrogenResidual=sum(d.mend_nitrogen_residual_mg_g(:) ...
                    .*bdTop./1000.*dzTop);
                nitrogenResidual=nitrogenResidual+mendNitrogenResidual;
                registry.write('mend_nitrogen_residual_profile', ...
                    d.mend_nitrogen_residual_mg_g,kt);
            end
            if isfield(mw,'integrated_flux_mg_g')
                names=mendCarbonFluxNames();
                for i=1:numel(names)
                    name=names{i}; amount=mw.integrated_flux_mg_g.(name);
                    amount=amount(:);
                    key=['mend_' lower(name)];
                    registry.write([key '_profile'],amount,kt);
                    registry.write([key '_flux'], ...
                        sum(amount.*bdTop./1000.*dzTop)./x.mend.dt_s,kt);
                end
                respiration=mw.integrated_flux_mg_g.CO2_growth ...
                    +mw.integrated_flux_mg_g.CO2_maintenance_MBA ...
                    +mw.integrated_flux_mg_g.CO2_maintenance_MBD ...
                    +mw.integrated_flux_mg_g.N_overflow_CO2;
                registry.write('mend_microbe_respiration_co2_profile',respiration,kt);
                registry.write('mend_microbe_respiration_co2_flux', ...
                    sum(respiration(:).*bdTop./1000.*dzTop)./x.mend.dt_s,kt);
            end
            registry.write('mend_reaction_substep_count',mw.substep_count,kt);
            registry.write('mend_rejected_substep_count',mw.rejected_substeps,kt);
            registry.write('mend_min_reaction_substep',mw.min_substep_s,kt);
            if isfield(mw,'organic_n')
                nitrogenResidual=nitrogenResidual ...
                    +mw.organic_n.nitrogen_residual_g_cm2;
            end
            if isfield(mw,'bioparticle_transport')
                carbonResidual=carbonResidual ...
                    +mw.bioparticle_transport.carbon_residual_g_cm2;
                nitrogenResidual=nitrogenResidual ...
                    +mw.bioparticle_transport.nitrogen_residual_g_cm2;
            end
        end
        registry.write('organic_carbon_residual',carbonResidual,kt);
        registry.write('organic_nitrogen_residual',nitrogenResidual,kt);
        registry.write('mend_carbon_residual',mendCarbonResidual,kt);
        registry.write('mend_nitrogen_residual',mendNitrogenResidual,kt);
        registry.write('organic_spinup_carbon_added',0,kt);
        registry.write('organic_spinup_nitrogen_added',0,kt);
    end
    registry.write('organic_fertilizer_carbon_added', ...
        report.organic_fertilizer_carbon_added_g_cm2,kt);
    registry.write('organic_fertilizer_nitrogen_added', ...
        report.organic_fertilizer_nitrogen_added_g_cm2,kt);
    registry.write('litter_carbon_added',report.litter_carbon_added_g_cm2,kt);
    registry.write('litter_nitrogen_added',report.litter_nitrogen_added_g_cm2,kt);
    registry.write('litter_direct_doc_carbon_added', ...
        report.litter_direct_doc_carbon_added_g_cm2,kt);
    registry.write('litter_direct_don_nitrogen_added', ...
        report.litter_direct_don_nitrogen_added_g_cm2,kt);
    registry.write('straw_carbon_added', ...
        report.straw_return_carbon_added_g_cm2,kt);
    registry.write('straw_nitrogen_added', ...
        report.straw_return_nitrogen_added_g_cm2,kt);
    % Change: 2026-10-30 Wenjie Liu -- external root input has disjoint
    % event and cumulative diagnostics, without crediting PlantCN litter.
    registry.write('root_residue_carbon_added',report.root_residue_carbon_added_g_cm2,kt);
    registry.write('root_residue_nitrogen_added',report.root_residue_nitrogen_added_g_cm2,kt);
    registry.write('root_residue_carbon_input_cumulative', ...
        cumulativeValue(state.cumulative,'input_root_residue_C_g_cm2'),kt);
    registry.write('root_residue_nitrogen_input_cumulative', ...
        cumulativeValue(state.cumulative,'input_root_residue_N_g_cm2'),kt);
    registry.write('manure_carbon_added',report.manure_carbon_added_g_cm2,kt);
    registry.write('manure_nitrogen_added',report.manure_nitrogen_added_g_cm2,kt);
    registry.write('manure_tan_nitrogen_added', ...
        report.organic_fertilizer_tan_added_g_cm2,kt);
    for material={'litter','straw','manure'}
        name=material{1};
        registry.write([name '_carbon_input_cumulative'], ...
            cumulativeValue(state.cumulative,['input_' name '_C_g_cm2']),kt);
        registry.write([name '_nitrogen_input_cumulative'], ...
            cumulativeValue(state.cumulative,['input_' name '_N_g_cm2']),kt);
        for product={'poc_c','doc_c','co2_c','pon_n','don_n','tan_n'}
            key=product{1}; stepValue=0;
            if exist('aom','var')&&isfield(aom,'materials') ...
                    &&isfield(aom.materials,name)
                sourceKey=diagnosticProductKey(key);
                stepValue=aom.materials.(name).(sourceKey);
            end
            registry.write([name '_to_' key],stepValue,kt);
            registry.write([name '_to_' key '_cumulative'], ...
                cumulativeValue(state.cumulative, ...
                    [name '_to_' cumulativeProductKey(key)]),kt);
        end
    end
end

function names=mendCarbonFluxNames()
%MENDCARBONFLUXNAMES Accepted integrated MEND carbon-ledger field names.
% Editor: Jizixing Zhu; Date: 2026-09-12.
    names={'POM1_decomp','POM2_decomp','POM1_to_DOM','POM2_to_DOM', ...
        'POM1_to_MOM','POM2_to_MOM','MOM_to_DOM','DOM_to_MBA', ...
        'CO2_growth','CO2_maintenance_MBA','CO2_maintenance_MBD', ...
        'MBA_mortality','MBA_to_DOM','MBA_to_POM1','MBA_to_POM2', ...
        'MBA_to_ENZP1','MBA_to_ENZP2','MBA_to_ENZM', ...
        'ENZP1_to_DOM','ENZP2_to_DOM','ENZM_to_DOM', ...
        'MBA_to_MBD','MBD_to_MBA','DOM_to_QOM','QOM_to_DOM', ...
        'DOM_sorption_net','ENZN_production','ENZN_to_MBA','N_overflow_CO2'};
end

function value=cumulativeValue(cumulative,name)
% Read a new material ledger safely from legacy checkpoints.
% Editor: Codex; Date: 2026-08-28.
    value=0;
    if isfield(cumulative,name), value=cumulative.(name); end
end

function key=diagnosticProductKey(product)
% Map output-safe product names to the AOM diagnostic contract.
% Editor: Codex; Date: 2026-08-28.
    map=struct('poc_c','poc_C_g_cm2','doc_c','doc_C_g_cm2', ...
        'co2_c','co2_C_g_cm2','pon_n','pon_N_g_cm2', ...
        'don_n','don_N_g_cm2','tan_n','tan_N_g_cm2');
    key=map.(product);
end

function key=cumulativeProductKey(product)
% Map output-safe names to cumulative AOM-ledger suffixes.
% Editor: Codex; Date: 2026-08-28.
    map=struct('poc_c','POC_C_g_cm2','doc_c','DOC_C_g_cm2', ...
        'co2_c','CO2_C_g_cm2','pon_n','PON_N_g_cm2', ...
        'don_n','DON_N_g_cm2','tan_n','TAN_N_g_cm2');
    key=map.(product);
end

function writeSoilNitrogen(registry, built, kt, report)
%WRITESOILNITROGEN  Export accepted soil-N state, sources, and ledger terms.
%
%   The N candidate is committed only after the accepted h/U/P_g endpoint.
%   This writer consequently reads the soil-owned accepted state and its
%   accepted solver diagnostic cache; it must not read a staged root request
%   or recompute a transport/reaction term for output.
%
%   Editor: Jizixing Zhu; Date: 2026-07-20.
% Change: 2026-09-14 Codex -- publish accepted O2 mass closure and reaction factor.
% Change: 2026-09-21 Wenjie Liu -- publish raw upward-positive mineral-N
%   boundary fluxes separately from nonnegative physical inputs and losses.
% Change: 2026-09-21 Wenjie Liu -- export the accepted passive
%   mineral-fertilizer source-label state, route ledgers, and closure.
% Change: 2026-10-23 Wenjie Liu -- fill the existing soil_oxygen group from
%   accepted O2 state and tile-mean kinetics; skip it when gas is disabled.
% Change: 2026-11-01 Wenjie Liu -- write accepted protection reductions
%   and node budgets without changing physical state or recomputing reactions.
    if ~isfield(built.control, 'switches') || ...
            ~isfield(built.control.switches, 'lNefc') || ...
            built.control.switches.lNefc ~= 1 || ...
            ~isfield(built.Y, 'soil_nitrogen')
        return;
    end
    state = built.Y.soil_nitrogen;
    flux = state.flux;
    variables = state.variables;
    diagnostic = struct();
    if isfield(built.workspace, 'nitrogen')
        diagnostic = built.workspace.nitrogen;
    end
    if isfield(diagnostic,'nh4_protection')
        protection=diagnostic.nh4_protection;
        reduction=protection.rate_reduction_ug_cm3_s;
        registry.write('nitrification_protection_reduction_profile',reduction,kt);
        registry.write('nitrification_protection_reduction_flux', ...
            sum(mean(reduction,2).*built.domain.dz(:)),kt);
        writeNode(registry,'nh4_protection_budget_profile', ...
            protection.node_budget_ug_cm2_s,kt);
    end
    if isfield(diagnostic,'oxygen')
        writeNode(registry,'oxygen_concentration',built.Y.soil_gas.C_O2,kt);
        for key={'storage_before','storage_after','consumed','top_input','bottom_input','balance_residual'}
            name=key{1};
            registry.write(['oxygen_' name],diagnostic.oxygen.([name '_mol_cm2']),kt);
        end
        registry.write('oxygen_reaction_factor',diagnostic.oxygen_reaction_factor_top_down,kt);
        if built.control.switches.lSoilairefc==1
            % The registered soil_oxygen group uses SI concentration units;
            % accepted state/ledger use mol per cm3. Registry maps layer
            % endpoints to top-down nodes for fractions and sinks.
            writeNode(registry,'o2_concentration',built.Y.soil_gas.C_O2*1e6,kt);
            writeLayerEndpoints(registry,'aerobic_fraction',variables,'f_aerobic',kt);
            writeLayerEndpoints(registry,'anaerobic_fraction',variables,'f_anaerobic',kt);
            registry.write('soil_gas_phase_fraction',built.Y.soil_gas.oxygen_phase.theta_gas,kt);
            appliedSink=diagnostic.oxygen.actual_sink_mol_cm3_s(:)*1e6;
            registry.write('o2_sink_total',repmat(appliedSink,1,2),kt);
        end
    end

    writeNode(registry, 'nh4_liquid', built.Y.nh4, kt);
    writeNode(registry, 'no3_liquid', built.Y.no3, kt);
    writeNode(registry, 'urea_liquid', built.Y.urea, kt);
    if isfield(built.control.switches, 'lNtotal') && ...
            built.control.switches.lNtotal == 1
        [nh4Observed,no3Observed]=soil.nitrogen.observationProfile( ...
            built.Y,built.parameters.soil_n,built.domain.dz);
        writeNode(registry,'nh4_total',nh4Observed,kt); % [ug g^-1]
        writeNode(registry,'no3_total',no3Observed,kt); % [ug g^-1]
    end
    writeLayerEndpoints(registry, 'nitrification_profile', variables, 'Snit', kt);
    writeLayerEndpoints(registry, 'denitrification_profile', variables, 'Sden', kt);
    writeLayerEndpoints(registry, 'mineralization_profile', variables, 'Smin', kt);
    if isfield(variables, 'SinkS') && ndims(variables.SinkS) == 3
        writeLayerEndpoints(registry, 'root_sink_profile', ...
            struct('value', variables.SinkS(:, :, 1) + variables.SinkS(:, :, 2)), ...
            'value', kt);
    end
    writeLayerEndpoints(registry, 'nit_mm_factor', variables, 'nit_mm_factor', kt);
    % Change: 2026-10-28 Wenjie Liu -- accepted-time dry-response diagnostic.
    writeLayerEndpoints(registry, 'nit_moisture_factor', variables, 'nit_moisture_factor', kt);
    writeLayerEndpoints(registry, 'nh4_substrate_total', variables, 'nit_M_total', kt);

    fluxMap = { ...
        'nitrification_flux', 'cvSnit'; 'denitrification_flux', 'cvSden'; ...
        'mineralization_flux', 'cvSmin'; 'ammonia_volatilization', 'cvSvol'; ...
        'root_sink_nh4', 'cvSinkSnh4'; 'root_sink_no3', 'cvSinkSno3'; ...
        'bot_flux_nh4', 'cvBot'; 'bot_flux_no3', 'cvBot1'; ...
        'top_flux_nh4', 'cvTop'; 'top_flux_no3', 'cvTop1'; ...
        'n2o_from_nitrification', 'Cvnit_N2o'; ...
        'n2o_from_denitrification', 'Cvdenit_N2o'; ...
        'n2_from_denitrification', 'Cvdenit_N2'; ...
        'nh4_storage', 'cvnh4'; 'no3_storage', 'cvno3'};
    for i = 1:size(fluxMap, 1)
        writeFlux(registry, fluxMap{i,1}, flux, fluxMap{i,2}, kt);
    end
    boundarySpecies={ ...
        'nh4','cvTop','cvBot'; ...
        'no3','cvTop1','cvBot1'};
    directionalRateMap={ ...
        'top_input','top_%s_input_flux'; ...
        'top_loss','top_%s_loss_flux'; ...
        'bottom_input','bottom_%s_input_flux'; ...
        'bottom_loss','bottom_%s_loss_flux'};
    for speciesIndex=1:size(boundarySpecies,1)
        species=boundarySpecies{speciesIndex,1};
        diagnosticField=['boundary_' species];
        if isfield(diagnostic,diagnosticField)
            boundaryLedger=diagnostic.(diagnosticField);
        else
            boundaryLedger=soil.nitrogen.boundaryFluxLedger( ...
                fluxValue(flux,boundarySpecies{speciesIndex,2}), ...
                fluxValue(flux,boundarySpecies{speciesIndex,3}));
        end
        for directionIndex=1:size(directionalRateMap,1)
            registry.write(sprintf(directionalRateMap{directionIndex,2},species), ...
                boundaryLedger.(directionalRateMap{directionIndex,1}),kt);
        end
    end
    if isfield(variables, 'Shys')
        shys = variables.Shys;
        if isequal(size(shys), [built.domain.NL, 2])
            registry.write('urea_hydrolysis_flux', ...
                sum(built.domain.dz(:) .* mean(shys, 2)), kt);
        end
    end
    writeScalar(registry, 'nh4_cum_residual', diagnostic, 'cum_residual_nh4', kt);
    writeScalar(registry, 'no3_cum_residual', diagnostic, 'cum_residual_no3', kt);
    registry.write('fert_nh4_added', report.fert_nh4_added, kt);
    registry.write('fert_no3_added', report.fert_no3_added, kt);
    registry.write('fert_urea_added', report.fert_urea_added, kt);
    cumulativeMap = { ...
        'fert_nh4_added_cumulative','cum_Fert_NH4'; ...
        'fert_no3_added_cumulative','cum_Fert_NO3'; ...
        'fert_urea_added_cumulative','cum_Fert_Urea'; ...
        'total_n_urea_hydrolysis_cumulative','cum_cvShys'; ...
        'total_n_mineralization_cumulative','cum_cvSmin'; ...
        'total_n_nitrification_cumulative','cum_cvSnit'; ...
        'total_n_nitrification_n2o_cumulative','cum_Cvnit_N2o'; ...
        'total_n_denitrification_cumulative','cum_cvSden'; ...
        'total_n_ammonia_volatilization_cumulative','cum_cvSvol'; ...
        'total_n_root_nh4_uptake_cumulative','cum_cvSinkSnh4'; ...
        'total_n_root_no3_uptake_cumulative','cum_cvSinkSno3'; ...
        'total_n_top_nh4_transfer_cumulative','cum_cvTop_NH4'; ...
        'total_n_bottom_nh4_transfer_cumulative','cum_cvBot_NH4'; ...
        'total_n_top_no3_transfer_cumulative','cum_cvTop_NO3'; ...
        'total_n_bottom_no3_transfer_cumulative','cum_cvBot_NO3'; ...
        'total_n_top_nh4_input_cumulative','cum_top_input_NH4'; ...
        'total_n_top_nh4_loss_cumulative','cum_top_loss_NH4'; ...
        'total_n_bottom_nh4_input_cumulative','cum_bottom_input_NH4'; ...
        'total_n_bottom_nh4_loss_cumulative','cum_bottom_loss_NH4'; ...
        'total_n_top_no3_input_cumulative','cum_top_input_NO3'; ...
        'total_n_top_no3_loss_cumulative','cum_top_loss_NO3'; ...
        'total_n_bottom_no3_input_cumulative','cum_bottom_input_NO3'; ...
        'total_n_bottom_no3_loss_cumulative','cum_bottom_loss_NO3'};
    for i=1:size(cumulativeMap,1)
        writeScalar(registry,cumulativeMap{i,1},flux,cumulativeMap{i,2},kt);
    end
    fertilizerTotal=fluxValue(flux,'cum_Fert_NH4') ...
        +fluxValue(flux,'cum_Fert_NO3')+fluxValue(flux,'cum_Fert_Urea');
    registry.write('fert_total_n_added_cumulative',fertilizerTotal,kt);
    if isfield(built.control.switches,'lFertilizerSourceTracer') ...
            &&built.control.switches.lFertilizerSourceTracer==1
        writeFertilizerSourceTracer(registry,built,kt);
    end
    % Change: 2026-10-30 Wenjie Liu -- write only accepted manure cohorts.
    if isfield(built.control.switches,'lManureSourceTracer') ...
            &&built.control.switches.lManureSourceTracer==1
        diagnostics.writeManureSourceTracer(registry,built,kt);
    end
    if isfield(built.diagnostics, 'root_demand')
        demand = built.diagnostics.root_demand;
        if numel(demand.nh4) >= kt, registry.write('root_demand_nh4', demand.nh4(kt), kt); end
        if numel(demand.no3) >= kt, registry.write('root_demand_no3', demand.no3(kt), kt); end
        if numel(demand.nh4) >= kt && numel(demand.no3) >= kt
            demandTotal = demand.nh4(kt) + demand.no3(kt);
            if demandTotal > 0
                uptake = fluxValue(flux, 'cvSinkSnh4') ...
                    + fluxValue(flux, 'cvSinkSno3');
                registry.write('nitrogen_supply_limit_factor', ...
                    max(0, min(1, uptake * 1e-2 / demandTotal)), kt);
            end
        end
    end
end

function writeFertilizerSourceTracer(registry,built,kt)
% Write only the committed fertilizer-source state and accepted ledgers.
% Editor: Wenjie Liu; Date: 2026-09-21.
% Change: 2026-09-22 Wenjie Liu -- retain detailed source routes and phase pools.
% Change: 2026-10-28 Wenjie Liu -- also accept dynamic DON tracer schema 5.
    assert(isfield(built.Y,'soil_fertilizer_n_tracer') ...
        &&isfield(built.workspace,'fertilizer_n_tracer'), ...
        'diagnostics:writeStep:missingFertilizerTracer', ...
        ['Enabled fertilizer source tracer requires committed state and ', ...
        'accepted diagnostics.']);
    state=built.Y.soil_fertilizer_n_tracer;
    diagnostic=built.workspace.fertilizer_n_tracer;
    % Change: 2026-09-29 Wenjie Liu -- schema 3 includes particle terminal loss.
    assert(isfield(state,'schema_version')&&ismember(state.schema_version,[4 5]) ...
        &&isfield(diagnostic,'schema_version')&&diagnostic.schema_version==1, ...
        'diagnostics:writeStep:fertilizerTracerSchema', ...
        'Unsupported fertilizer source-tracer diagnostic schema.');

    nodeMap={ ...
        'fertilizer_tracer_nh4_liquid','nh4'; ...
        'fertilizer_tracer_no3_liquid','no3'; ...
        'fertilizer_tracer_urea_liquid','urea'; ...
        'fertilizer_tracer_nh4_granule','nh4_granule'};
    for i=1:size(nodeMap,1)
        assert(isfield(state,nodeMap{i,2}), ...
            'diagnostics:writeStep:fertilizerTracerNode', ...
            'Fertilizer source tracer is missing node field %s.',nodeMap{i,2});
        writeNode(registry,nodeMap{i,1},state.(nodeMap{i,2}),kt);
    end

    assert(isfield(state,'dom')&&isfield(state.dom,'don_total_g_cm3') ...
        &&isfield(state,'organic_n')&&isstruct(state.organic_n), ...
        'diagnostics:writeStep:fertilizerTracerOrganic', ...
        'Fertilizer source tracer lacks DON or organic-N owners.');
    writeTracerLayer(registry,'fertilizer_tracer_don_profile', ...
        state.dom.don_total_g_cm3,built.domain.NL,kt);
    layerMap={ ...
        'fertilizer_tracer_pon1_profile','PON1'; ...
        'fertilizer_tracer_pon2_profile','PON2'; ...
        'fertilizer_tracer_mon_profile','MON'; ...
        'fertilizer_tracer_active_microbe_n_profile','MB_N'; ...
        'fertilizer_tracer_dormant_microbe_n_profile','BD_N'; ...
        'fertilizer_tracer_pom1_enzyme_n_profile','ENZP1_N'; ...
        'fertilizer_tracer_pom2_enzyme_n_profile','ENZP2_N'; ...
        'fertilizer_tracer_mom_enzyme_n_profile','EM_N'; ...
        'fertilizer_tracer_intracellular_enzyme_n_profile','ENZN_N'};
    for i=1:size(layerMap,1)
        field=layerMap{i,2};
        assert(isfield(state.organic_n,field), ...
            'diagnostics:writeStep:fertilizerTracerOrganic', ...
            'Fertilizer source tracer is missing organic-N field %s.',field);
        writeTracerLayer(registry,layerMap{i,1},state.organic_n.(field), ...
            built.domain.NL,kt);
    end

    assert(isfield(diagnostic,'storage')&&isstruct(diagnostic.storage) ...
        &&isfield(diagnostic,'step')&&isstruct(diagnostic.step) ...
        &&isfield(state,'cumulative')&&isstruct(state.cumulative), ...
        'diagnostics:writeStep:fertilizerTracerLedger', ...
        'Fertilizer source tracer lacks accepted storage or route ledgers.');
    for key={'nh4','no3','urea','organic','total'}
        name=key{1};field=[name '_g_cm2'];
        registry.write(['fertilizer_tracer_' name '_storage'], ...
            tracerScalar(diagnostic.storage,field),kt);
    end
    routeFields=fertilizerTracerRouteFields();
    for i=1:numel(routeFields)
        name=routeFields{i};field=[name '_g_cm2'];
        registry.write(['fertilizer_tracer_' name '_step'], ...
            tracerScalar(diagnostic.step,field),kt);
        registry.write(['fertilizer_tracer_' name '_cumulative'], ...
            tracerScalar(state.cumulative,field),kt);
    end
    registry.write('fertilizer_tracer_mend_closure_residual_step', ...
        tracerScalar(diagnostic.step,'mend_closure_residual_g_cm2'),kt);
    registry.write('fertilizer_tracer_din_closure_residual_step', ...
        tracerScalar(diagnostic.step,'din_closure_residual_g_cm2'),kt);
    registry.write('fertilizer_tracer_cumulative_closure_residual', ...
        tracerScalar(diagnostic,'cumulative_closure_residual_g_cm2'),kt);
    diagnostics.writeFertilizerTracerDetails(registry,built,kt);
end

function writeTracerLayer(registry,name,value,nl,kt)
% Reverse one bottom-up accepted tracer layer profile exactly once.
% Editor: Wenjie Liu; Date: 2026-09-21.
    assert(isnumeric(value)&&isreal(value)&&isvector(value) ...
        &&numel(value)==nl&&all(isfinite(value(:))), ...
        'diagnostics:writeStep:fertilizerTracerLayer', ...
        'Accepted fertilizer tracer %s must be one finite NL-layer vector.',name);
    registry.write(name,flip(value(:)),kt);
end

function value=tracerScalar(source,field)
% Read one required finite scalar from an accepted tracer ledger.
% Editor: Wenjie Liu; Date: 2026-09-21.
    assert(isstruct(source)&&isfield(source,field) ...
        &&isnumeric(source.(field))&&isreal(source.(field)) ...
        &&isscalar(source.(field))&&isfinite(source.(field)), ...
        'diagnostics:writeStep:fertilizerTracerScalar', ...
        'Accepted fertilizer tracer ledger lacks finite scalar %s.',field);
    value=source.(field);
end

function fields=fertilizerTracerRouteFields()
% Route field names shared by accepted step and cumulative tracer ledgers.
% Editor: Wenjie Liu; Date: 2026-09-21.
    fields={ ...
        'fertilizer_input','root_nh4','root_no3','nh3_volatilization', ...
        'nitrification_n2o','denitrification_n2o','denitrification_n2', ...
        'nh4_top_input','nh4_top_loss','nh4_bottom_input','nh4_bottom_loss', ...
        'no3_top_input','no3_top_loss','no3_bottom_input','no3_bottom_loss', ...
        'don_bottom_loss','bioparticle_mb_bottom_loss','bioparticle_bd_bottom_loss', ...
        'bioparticle_enzn_bottom_loss','don_top_loss','bioparticle_mb_top_loss', ...
        'bioparticle_bd_top_loss','bioparticle_enzn_top_loss', ...
        'organic_entry_gross','organic_release_gross','numerical_correction'};
end

function writeNode(registry, name, value, kt)
    assert(isnumeric(value) && isvector(value) ...
            && numel(value) == registry.NN && all(isfinite(value(:))), ...
        'diagnostics:writeStep:badNitrogenNode', ...
        'Accepted diagnostic %s must be one finite NN-node vector.', name);
    registry.write(name, flip(value(:)), kt);
end

function writeLayerEndpoints(registry, name, source, field, kt)
    if isfield(source, field) && ~isempty(source.(field))
        registry.write(name, source.(field), kt);
    end
end

function writeFlux(registry, name, source, field, kt)
    if isfield(source, field)
        registry.write(name, fluxValue(source, field), kt);
    end
end

function writeScalar(registry, name, source, field, kt)
    if isfield(source, field) && isscalar(source.(field)) ...
            && isfinite(source.(field))
        registry.write(name, source.(field), kt);
    end
end

function value = fluxValue(flux, field)
    value = flux.(field);
    if ~isscalar(value), value = value(end); end
end

function writeSoilInterfaceClosure(registry, built, kt)
% Export one accepted outer strict-Dirichlet temperature closure.
% Editor: Jizixing Zhu; Date: 2026-07-19.
% Change: 2026-07-20 Jizixing Zhu -- export atmospheric closure fallback mode,
%         bounded flag, applied threshold, and refined time-tile cap.
    if ~isfield(built.workspace, 'surface') ...
            || ~isfield(built.workspace.surface, ...
                'soil_interface_temperature_closure')
        return;
    end
    c = built.workspace.surface.soil_interface_temperature_closure;
    mappings = { ...
        'soil_interface_closure_initial_temperature', ...
            'initial_temperature_C'; ...
        'soil_interface_closure_final_temperature', 'final_temperature_C'; ...
        'soil_interface_closure_initial_residual', 'initial_residual_W_m2'; ...
        'soil_interface_closure_final_residual', 'final_residual_W_m2'; ...
        'soil_interface_closure_soil_flux', 'soil_flux_W_m2'; ...
        'soil_interface_closure_surface_flux', 'surface_flux_W_m2'; ...
        'soil_interface_closure_iterations', 'iterations'; ...
        'soil_interface_closure_evaluations', 'evaluations'; ...
        'soil_interface_closure_converged', 'converged'; ...
        'soil_interface_closure_acceptance_mode', 'acceptance_mode_code'; ...
        'soil_interface_closure_accepted_bounded_residual', ...
            'accepted_bounded_residual'; ...
        'soil_interface_closure_bounded_tolerance', ...
            'bounded_tolerance_W_m2'; ...
        'soil_interface_closure_substep_retry_used', ...
            'substep_retry_used'; ...
        'soil_interface_closure_chosen_substep_cap', ...
            'chosen_substep_cap_s'};
    for i = 1:size(mappings, 1)
        registry.write(mappings{i,1}, double(c.(mappings{i,2})), kt);
    end
end

function writeSurfaceSoilEnergyBalance(registry, built, kt)
% Export the accepted surface-snow-soil total-energy ledger.
% Editor: Jizixing Zhu; Date: 2026-07-18.
% Change: 2026-07-19 Jizixing Zhu -- include prognostic snow storage.
    if ~isfield(built.workspace, 'flux') ...
            || ~isfield(built.workspace.flux, ...
                'surface_soil_energy_balance')
        return;
    end
    b = built.workspace.flux.surface_soil_energy_balance;
    if ~isfield(b, 'has_data') || ~b.has_data
        return;
    end
    mappings = { ...
        'ground_net_radiation', 'ground_net_radiation_W_m2'; ...
        'ground_sensible_heat', 'ground_sensible_heat_W_m2'; ...
        'ground_latent_heat_demand', 'ground_latent_heat_demand_W_m2'; ...
        'ground_latent_heat_actual', 'ground_latent_heat_actual_W_m2'; ...
        'diagnosed_bare_liquid_infiltration_sensible_energy_flux', ...
            'diagnosed_bare_liquid_infiltration_sensible_energy_flux_W_m2'; ...
        'diagnosed_dry_air_sensible_energy_transport', ...
            'diagnosed_dry_air_sensible_energy_transport_W_m2'; ...
        'bare_ground_external_energy_flux_actual', ...
            'bare_ground_external_energy_flux_actual_W_m2'; ...
        'bare_ground_external_energy_flux_reconstructed', ...
            'bare_ground_external_energy_flux_reconstructed_W_m2'; ...
        'bare_ground_external_energy_flux_reconstruction_error', ...
            'bare_ground_external_energy_flux_reconstruction_error_W_m2'; ...
        'snow_external_energy_flux', 'snow_external_energy_flux_W_m2'; ...
        'snow_to_soil_internal_energy_flux', ...
            'snow_to_soil_internal_energy_flux_W_m2'; ...
        'surface_external_energy_flux', 'surface_external_energy_flux_W_m2'; ...
        'scope_ground_heat_flux_empirical', ...
            'scope_ground_heat_flux_empirical_W_m2'; ...
        'scope_ground_ebal_residual', 'scope_ground_ebal_residual_W_m2'; ...
        'soil_top_energy_flux_total_actual', ...
            'soil_top_energy_flux_total_actual_W_m2'; ...
        'surface_to_soil_interface_power_residual', ...
            'surface_to_soil_interface_power_residual_W_m2'; ...
        'snow_internal_energy_change', 'snow_internal_energy_change_J_m2'; ...
        'surface_soil_total_energy_change', ...
            'surface_soil_total_energy_change_J_m2'; ...
        'surface_external_energy', 'surface_external_energy_J_m2'; ...
        'surface_to_soil_interface_energy_residual', ...
            'surface_to_soil_interface_energy_residual_J_m2'; ...
        'surface_soil_total_energy_residual', ...
            'surface_soil_total_energy_residual_J_m2'; ...
        'surface_soil_interface_cancellation_error', ...
            'surface_soil_interface_cancellation_error_J_m2'; ...
        'surface_energy_precipitation_temperature', ...
            'precipitation_temperature_C'};
    for i = 1:size(mappings, 1)
        registry.write(mappings{i,1}, b.(mappings{i,2}), kt);
    end
    if kt == 1 && isfield(built, 'p') ...
            && isfield(built.p, 'initialization') ...
            && isfield(built.p.initialization, 'thermal_bottom_alignment')
        initReport = built.p.initialization.thermal_bottom_alignment;
        registry.write('initial_bottom_temperature_mismatch_before_alignment', ...
            readFinite(initReport, 'initial_mismatch_K'), kt);
        registry.write('initial_bottom_temperature_mismatch_after_alignment', ...
            readFinite(initReport, 'final_mismatch_K'), kt);
    end
end

function value = readFinite(container, name)
% Read one finite scalar or NaN for an inapplicable diagnostic.
% Editor: Jizixing Zhu; Date: 2026-07-18.
    value = NaN;
    if isstruct(container) && isfield(container, name) ...
            && isscalar(container.(name)) && isfinite(container.(name))
        value = double(container.(name));
    end
end

function writeSnow(registry, built, kt, dt_s)
% Export the accepted ClimaLand-compatible snow state and flux transaction.
% Editor: Jizixing Zhu; Date: 2026-07-17.
% Change: 2026-07-17 Jizixing Zhu -- write terminal ice-free meltwater drainage.
% Change: 2026-07-18 Jizixing Zhu -- write its conservative enthalpy transfer.
    if ~isfield(built, 'diagnostics') || ~isfield(built.diagnostics, 'snow')
        return;
    end
    f = built.diagnostics.snow;
    if isempty(f) || ~isfield(f, 'state'), return; end
    d = f.state;
    registry.write('snow_water_equivalent', d.SWE_cm, kt);
    registry.write('snow_liquid_water', d.liquid_water_cm, kt);
    registry.write('snow_ice_water_equivalent', d.ice_water_equivalent_cm, kt);
    registry.write('snow_internal_energy', d.energy_J_cm2 * 1e4, kt);
    registry.write('snow_temperature', d.temperature_C, kt);
    registry.write('snow_depth', d.depth_m, kt);
    registry.write('snow_density', d.density_kg_m3, kt);
    registry.write('snow_cover_fraction', d.cover_fraction, kt);
    writeRateAmount(registry, 'snowfall', f.snowfall_cm_s, dt_s, kt);
    writeRateAmount(registry, 'rainfall', f.rainfall_cm_s, dt_s, kt);
    writeRateAmount(registry, 'snow_runoff', f.snow_runoff_cm_s, dt_s, kt);
    writeRateAmount(registry, 'snow_complete_meltout', ...
        readSnowFlux(f, 'complete_meltout_cm_s', 0), dt_s, kt);
    meltoutEnergyFlux = readSnowFlux( ...
        f, 'complete_meltout_energy_W_m2', 0);
    registry.write('snow_complete_meltout_energy_flux', ...
        meltoutEnergyFlux, kt);
    registry.write('snow_complete_meltout_energy', ...
        meltoutEnergyFlux * dt_s, kt);
    writeRateAmount(registry, 'snow_vapor_loss', f.snow_vapor_flux_cm_s, dt_s, kt);
    registry.write('frozen_surface_water_to_snow', ...
        f.frozen_pond_transfer_cm * 10, kt);
    registry.write('snow_soil_heat_flux', f.snow_soil_heat_flux_W_m2, kt);
    registry.write('snow_surface_energy_flux', f.total_energy_flux_W_m2, kt);
    registry.write('snow_applied_energy_flux', f.applied_energy_flux_W_m2, kt);
    registry.write('snow_phase_change_flux', f.phase_change_flux_m_s * 100, kt);
    registry.write('snow_water_balance_residual', ...
        f.water_balance_residual_cm * 10, kt);
    registry.write('snow_energy_balance_residual', ...
        f.energy_balance_residual_J_cm2 * 1e4, kt);
end

function value = readSnowFlux(flux, name, fallback)
% Read an optional snow diagnostic while preserving older checkpoints.
% Editor: Jizixing Zhu; Date: 2026-07-17.
    value = fallback;
    if isfield(flux, name) && ~isempty(flux.(name))
        value = flux.(name);
    end
end

function writeRateAmount(registry, prefix, rate, dt_s, kt)
% Write one snow rate and its accepted-step integral.
% Editor: Jizixing Zhu; Date: 2026-07-17.
    registry.write([prefix '_flux'], rate, kt);
    registry.write([prefix '_amount'], rate * dt_s * 10, kt);
end

function writeDryAirMassBalance(registry, built, kt)
% Export the accepted dry-air ledger in [g m^-2].
    if ~isfield(built.workspace, 'flux') ...
            || ~isfield(built.workspace.flux, ...
                'soil_dry_air_mass_balance')
        return;
    end
    balance = built.workspace.flux.soil_dry_air_mass_balance;
    if ~isfield(balance, 'has_data') || ~balance.has_data
        return;
    end
    to_g_m2 = 1e4;
    registry.write('soil_dry_air_mass', ...
        balance.storage_after_g_cm2 * to_g_m2, kt);
    registry.write('soil_dry_air_mass_change', ...
        balance.storage_change_g_cm2 * to_g_m2, kt);
    registry.write('soil_dry_air_top_input', ...
        balance.top_input_g_cm2 * to_g_m2, kt);
    registry.write('soil_dry_air_bottom_input', ...
        balance.bottom_input_g_cm2 * to_g_m2, kt);
    registry.write('soil_dry_air_boundary_net_input', ...
        balance.boundary_net_input_g_cm2 * to_g_m2, kt);
    registry.write('soil_dry_air_mass_balance_residual', ...
        balance.residual_g_cm2 * to_g_m2, kt);
end

function writeSoilBoundary(registry, built, kt)
% Export accepted heat and dry-air boundary results.
% Editor: Jizixing Zhu; Date: 2026-07-16.
% Change: 2026-07-17 Jizixing Zhu -- write snow-aware heat mode and request metadata.
% Change: 2026-07-17 Jizixing Zhu -- write bottom heat mode and requested-value metadata.
% Change: 2026-07-19 Jizixing Zhu -- write affine top heat-response metadata.
    if ~isfield(built.workspace.flux, 'soil_boundary_hourly')
        return;
    end
    boundary = built.workspace.flux.soil_boundary_hourly;
    if ~boundary.has_data
        return;
    end
    if boundary.heat.enabled
        registry.write('soil_temperature', flip(built.Y.T(:)), kt);
        if isfield(built.boundary, 'heat') ...
                && strcmp(built.boundary.heat.bottom.kind, ...
                    'seasonal_bottom_temperature') ...
                && isfield(built.workspace, 'boundary') ...
                && isfield(built.workspace.boundary, 'thermal') ...
                && numel(built.workspace.boundary.thermal.Tbh) >= kt
            registry.write('bottom_temperature', ...
                built.workspace.boundary.thermal.Tbh(kt), kt);
        end
        registry.write('surface_heat_flux', ...
            boundary.heat.top_flux_into_domain * 1e4, kt);
        registry.write('bottom_heat_flux', ...
            boundary.heat.bottom_flux_into_domain * 1e4, kt);
        registry.write('soil_bottom_heat_boundary_mode', ...
            boundary.heat.bottom_boundary_mode_code, kt);
        registry.write('soil_bottom_boundary_temperature', ...
            boundary.heat.bottom_requested_temperature_C, kt);
        registry.write('soil_bottom_requested_heat_flux', ...
            boundary.heat.bottom_requested_flux_into_domain * 1e4, kt);
        registry.write('soil_top_heat_boundary_mode', ...
            boundary.heat.top_boundary_mode_code, kt);
        registry.write('soil_top_boundary_temperature', ...
            boundary.heat.top_requested_temperature_C, kt);
        registry.write('soil_top_requested_heat_flux', ...
            boundary.heat.top_requested_flux_into_domain * 1e4, kt);
        registry.write('soil_top_energy_response_reference_temperature', ...
            boundary.heat.top_reference_temperature_C, kt);
        registry.write('soil_top_energy_flux_temperature_sensitivity', ...
            boundary.heat.top_flux_temperature_derivative * 1e4, kt);
        registry.write( ...
            'soil_top_energy_response_evaluation_temperature', ...
            boundary.heat.top_response_evaluation_temperature_C, kt);
        registry.write('soil_top_energy_response_evaluated_flux', ...
            boundary.heat.top_response_evaluated_flux_into_domain * 1e4, kt);
        registry.write('surface_skin_temperature', ...
            boundary.heat.surface_skin_temperature_C, kt);
    end
    if boundary.dry_air.enabled
        registry.write('air_pressure', flip(built.Y.p_gas(:)) / 10, kt);
        registry.write('top_dry_air_mass_flux', ...
            boundary.dry_air.top_flux_into_domain, kt);
        registry.write('bottom_dry_air_mass_flux', ...
            boundary.dry_air.bottom_flux_into_domain, kt);
    end
end

function writeSoilEnergyBalance(registry, built, kt)
% Export the accepted soil-energy ledger in [J m^-2].
    if ~isfield(built.workspace, 'flux') ...
            || ~isfield(built.workspace.flux, 'soil_energy_balance')
        return;
    end
    balance = built.workspace.flux.soil_energy_balance;
    if ~isfield(balance, 'has_data') || ~balance.has_data
        return;
    end
    to_j_m2 = 1e4;
    registry.write('soil_internal_energy', ...
        balance.storage_after_J_cm2 * to_j_m2, kt);
    registry.write('soil_internal_energy_change', ...
        balance.storage_change_J_cm2 * to_j_m2, kt);
    registry.write('soil_energy_top_input', ...
        balance.top_input_J_cm2 * to_j_m2, kt);
    registry.write('soil_energy_bottom_input', ...
        balance.bottom_input_J_cm2 * to_j_m2, kt);
    registry.write('soil_energy_boundary_net_input', ...
        balance.boundary_net_input_J_cm2 * to_j_m2, kt);
    registry.write('soil_energy_root_sink', ...
        balance.root_sink_J_cm2 * to_j_m2, kt);
    registry.write('soil_energy_balance_residual', ...
        balance.residual_J_cm2 * to_j_m2, kt);
end

function writeWaterBalance(registry, built, kt, dt_s)
% Export the accepted water ledger and hydrologic boundary rates.
% Editor: Jizixing Zhu; Date: 2026-07-18.
    if ~isfield(built.workspace, 'flux') ...
            || ~isfield(built.workspace.flux, 'water_balance')
        return;
    end
    balance = built.workspace.flux.water_balance;
    if ~isfinite(balance.residual_mm)
        return;
    end
    hourly = built.workspace.flux.water_hourly;
    registry.write('water_balance_residual', balance.residual_mm, kt);
    registry.write('water_storage_column', ...
        balance.soil_storage_after_mm, kt);
    registry.write('water_storage_change', ...
        balance.soil_storage_change_mm, kt);

    evaporationDemand_cm_s = built.workspace.flux.EVAP(1);
    evaporation_cm_s = evaporationDemand_cm_s;
    if isfield(built.workspace.flux, 'surface_boundary_hourly') ...
            && built.workspace.flux.surface_boundary_hourly.has_data
        surfaceSummary = built.workspace.flux.surface_boundary_hourly;
        evaporation_cm_s = surfaceSummary.mean_evap_actual_cm_s;
        evaporationDemand_cm_s = surfaceSummary.mean_evap_demand_cm_s;
    end
    transpiration_cm_s = waterOnlyTranspiration(built);
    evaporationDeficit_cm_s = max(evaporationDemand_cm_s - evaporation_cm_s, 0);
    canopyNet_cm_s = 0;
    if isfield(built.workspace.flux, 'CanopyInterception')
        canopyNet_cm_s = built.workspace.flux.CanopyInterception(1);
    end
    wetCanopyEvaporation_cm_s = max(canopyNet_cm_s, 0);
    canopyDew_cm_s = max(-canopyNet_cm_s, 0);
    bareVapor_cm_s = scopeSignedMassRateCmS( ...
        built, 'bare_ground_upward_water_vapor_flux_kg_m2_s');
    if isfinite(bareVapor_cm_s)
        bareDew_cm_s = max(-bareVapor_cm_s, 0);
    else
        bareDew_cm_s = scopeMassRateCmS(built, 'bare_dew_kg_m2_s');
    end
    snowVapor_cm_s = 0;
    if isfield(built, 'diagnostics') && isfield(built.diagnostics, 'snow')
        snowVapor_cm_s = readSnowFlux( ...
            built.diagnostics.snow, 'snow_vapor_flux_cm_s', 0);
    end
    snowSublimation_cm_s = max(snowVapor_cm_s, 0);
    snowDeposition_cm_s = max(-snowVapor_cm_s, 0);
    evapotranspirationTotal_cm_s = evaporation_cm_s ...
        + wetCanopyEvaporation_cm_s + max(transpiration_cm_s, 0) ...
        + snowSublimation_cm_s;
    netSurfaceVaporExchange_cm_s = evapotranspirationTotal_cm_s ...
        - bareDew_cm_s - canopyDew_cm_s - snowDeposition_cm_s;
    registry.write('soil_evaporation_flux', evaporation_cm_s, kt);
    registry.write('soil_evaporation_amount', ...
        evaporation_cm_s * dt_s * 10, kt);
    registry.write('soil_evaporation_demand_flux', ...
        evaporationDemand_cm_s, kt);
    registry.write('soil_evaporation_demand_amount', ...
        evaporationDemand_cm_s * dt_s * 10, kt);
    registry.write('soil_evaporation_deficit_flux', ...
        evaporationDeficit_cm_s, kt);
    registry.write('soil_evaporation_deficit_amount', ...
        evaporationDeficit_cm_s * dt_s * 10, kt);
    registry.write('transpiration_flux', transpiration_cm_s, kt);
    registry.write('transpiration_amount', ...
        transpiration_cm_s * dt_s * 10, kt);
    writeRateAmount(registry, 'bare_surface_dew', bareDew_cm_s, dt_s, kt);
    writeRateAmount(registry, 'canopy_dew', canopyDew_cm_s, dt_s, kt);
    writeRateAmount(registry, 'snow_sublimation', ...
        snowSublimation_cm_s, dt_s, kt);
    writeRateAmount(registry, 'snow_deposition', ...
        snowDeposition_cm_s, dt_s, kt);
    writeRateAmount(registry, 'evapotranspiration_total', ...
        evapotranspirationTotal_cm_s, dt_s, kt);
    writeRateAmount(registry, 'net_surface_vapor_exchange', ...
        netSurfaceVaporExchange_cm_s, dt_s, kt);
    registry.write('soil_top_water_flux', -hourly.QMT, kt);
    registry.write('soil_top_water_amount', -hourly.QMT * dt_s * 10, kt);
    registry.write('soil_bottom_water_flux', hourly.QMB, kt);
    registry.write('soil_bottom_water_amount', hourly.QMB * dt_s * 10, kt);
    writeSurfaceWater(registry, built, hourly, kt, dt_s);
    writeSurfaceSoilWaterBalance(registry, built, kt);

    if size(built.Y.h, 2) == 1
        zero_head = soil.water.computeZeroHeadDiagnostics( ...
            built.Y.h, built.domain.z_node);
        registry.write('soil_matric_head', flip(built.Y.h(:)), kt);
        registry.write('soil_saturated_node_mask', ...
            double(zero_head.node_saturated), kt);
        registry.write('zero_head_crossing_depth', ...
            zero_head.crossing_depth_cm, kt);
        registry.write('zero_head_crossing_count', ...
            zero_head.crossing_count, kt);
        registry.write('shallowest_zero_head_depth', ...
            zero_head.shallowest_crossing_depth_cm, kt);
        registry.write('saturated_profile_fraction', ...
            zero_head.saturated_fraction, kt);
    end
end

function writeSurfaceWater(registry, built, hourly, kt, dt_s)
% Export accepted precipitation, interception, and runoff rates and amounts.
    precipitation = 0;
    if isfield(built.forcing, 'hourly') ...
            && isfield(built.forcing.hourly, 'Precip') ...
            && numel(built.forcing.hourly.Precip) >= kt
        precipitation = built.forcing.hourly.Precip(kt);
    end
    interception = 0;
    if isfield(built.workspace.flux, 'CanopyInterception')
        interception = max(built.workspace.flux.CanopyInterception(1), 0);
    end
    supply = hourly.Throughfall;
    runoff = hourly.runoffHort;
    storageChangeMm = NaN;
    if isfield(built.workspace.flux, 'surface_boundary_hourly') ...
            && built.workspace.flux.surface_boundary_hourly.has_data
        summary = built.workspace.flux.surface_boundary_hourly;
        supply = summary.mean_liquid_supply_cm_s;
        runoff = summary.mean_runoff_cm_s;
        storageChangeMm = summary.storage_tendency_integral_cm * 10;
    end
    registry.write('precipitation_flux', precipitation, kt);
    registry.write('precipitation_amount', precipitation * dt_s * 10, kt);
    registry.write('surface_liquid_water_supply_flux', supply, kt);
    registry.write('surface_liquid_water_supply_amount', ...
        supply * dt_s * 10, kt);
    registry.write('wet_canopy_evaporation_flux', interception, kt);
    registry.write('wet_canopy_evaporation_amount', ...
        interception * dt_s * 10, kt);
    registry.write('surface_runoff_flux', runoff, kt);
    registry.write('surface_runoff_amount', runoff * dt_s * 10, kt);
    if isfield(built.Y, 'surface_water') ...
            && isfield(built.Y.surface_water, 'DSTOR0')
        registry.write('surface_water_storage', ...
            built.Y.surface_water.DSTOR0, kt);
    end
    if isfinite(storageChangeMm)
        registry.write('surface_water_storage_change', storageChangeMm, kt);
    end
    writeSurfaceComplementarity(registry, built, kt);
end

function writeSurfaceComplementarity(registry, built, kt)
% Export the accepted outer-step summary of the two surface KKT pairs.
% Editor: Jizixing Zhu; Date: 2026-07-17.
    if ~isfield(built.workspace.flux, 'surface_boundary_hourly')
        return;
    end
    summary = built.workspace.flux.surface_boundary_hourly;
    if ~isfield(summary, 'has_data') || ~summary.has_data
        return;
    end
    registry.write('surface_pressure_head', summary.mean_psi_top_cm, kt);
    registry.write('surface_water_available_flux', ...
        summary.mean_q_available_cm_s, kt);
    registry.write('surface_boundary_mode', summary.final_mode_code, kt);
    registry.write('surface_mass_residual_max', ...
        summary.max_abs_mass_residual_cm, kt);
    registry.write('surface_lower_complementarity_residual_max', ...
        summary.max_lower_complementarity_cm2, kt);
    registry.write('surface_upper_complementarity_residual_max', ...
        summary.max_upper_complementarity_cm2_s, kt);
end

function writeSurfaceSoilWaterBalance(registry, built, kt)
% Export the post-canopy/post-snow combined surface-soil water ledger.
% Editor: Jizixing Zhu; Date: 2026-07-17.
    if ~isfield(built.workspace.flux, 'surface_soil_water_balance')
        return;
    end
    balance = built.workspace.flux.surface_soil_water_balance;
    if ~balance.has_data
        return;
    end
    registry.write('surface_soil_liquid_input', ...
        balance.liquid_supply_mm, kt);
    registry.write('surface_soil_evaporation_loss', ...
        balance.surface_evaporation_mm, kt);
    registry.write('surface_soil_runoff_loss', balance.runoff_mm, kt);
    registry.write('surface_soil_bottom_loss', balance.bottom_loss_mm, kt);
    registry.write('surface_soil_root_sink', balance.root_sink_mm, kt);
    registry.write('surface_soil_storage_change', ...
        balance.total_storage_change_mm, kt);
    registry.write('surface_soil_water_balance_residual', ...
        balance.residual_mm, kt);
end

function transpiration_cm_s = waterOnlyTranspiration(built)
% Read current ebal transpiration [cm/s], or zero when not staged.
% Editor: Jizixing Zhu; Date: 2026-07-16.
    transpiration_cm_s = 0;
    if isfield(built.workspace, 'scope_step') ...
            && isfield(built.workspace.scope_step, 'hydraulic') ...
            && isfield(built.workspace.scope_step.hydraulic, 'Tp_t') ...
            && ~isempty(built.workspace.scope_step.hydraulic.Tp_t)
        transpiration_cm_s = built.workspace.scope_step.hydraulic.Tp_t;
    end
end

function rate_cm_s = scopeMassRateCmS(built, fieldName)
% Convert one accepted SCOPE grid vapor mass contribution to cm s-1.
% Editor: Jizixing Zhu; Date: 2026-07-18.
    rate_cm_s = 0;
    if ~isfield(built.workspace, 'scope_step') ...
            || ~isfield(built.workspace.scope_step, 'fluxes')
        return;
    end
    fluxes = built.workspace.scope_step.fluxes;
    if isfield(fluxes, fieldName) && isscalar(fluxes.(fieldName)) ...
            && isfinite(fluxes.(fieldName))
        rate_cm_s = max(double(fluxes.(fieldName)), 0) * 0.1;
    end
end

function rate_cm_s = scopeSignedMassRateCmS(built, fieldName)
% Convert one signed SCOPE grid vapor mass flux to cm s-1 without clipping.
% Positive is upward loss; negative is condensation toward the surface.
% Editor: Jizixing Zhu; Date: 2026-07-18.
    rate_cm_s = NaN;
    if ~isfield(built.workspace, 'scope_step') ...
            || ~isfield(built.workspace.scope_step, 'fluxes')
        return;
    end
    fluxes = built.workspace.scope_step.fluxes;
    if isfield(fluxes, fieldName) && isscalar(fluxes.(fieldName)) ...
            && isfinite(fluxes.(fieldName))
        rate_cm_s = double(fluxes.(fieldName)) * 0.1;
    end
end
