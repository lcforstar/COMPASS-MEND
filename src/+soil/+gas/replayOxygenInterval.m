function oxygen = replayOxygenInterval(built, schedule, respiratoryC, nitrogenCap, nodeBudget)
%REPLAYOXYGENINTERVAL Replay original O2 equations on accepted carrier tiles.
% Inputs: built accepted y_O2 and static soil-N parameters; a value-only
% accepted tile schedule; same-hour frozen respiratoryC [ug C cm^-3 s^-1,
% NL-by-1 bottom-up]. Optional nitrogenCap [ug N cm^-3 s^-1, NL-by-2] is an
% explicit upper bound, absent in the production original-chemistry path.
% Returns a detached original hour accumulator, y_O2 and derived C_O2. No
% state is written. Infeasible respiratory demand or failed original Picard
% convergence returns feasible=false; no reaction is silently scaled.
% Source: COMPASS-COMPASS_MATLAB advanceOxygenSubstep (30 iterations, 1e-6)
% and unchanged solveDryMoleFractionTransport storage/transport/boundaries.
% Editor: Wenjie Liu; Date: 2026-10-23.
% Change: 2026-10-23 Wenjie Liu -- pass the optional biological pH-response
% switch in gas-on and gas-off kinetics; actual soil pH remains unchanged.
% Change: 2026-10-28 Wenjie Liu -- use each frozen tile's matric head for
% dry-side nitrification and publish its accepted-time-mean response.
% Change: 2026-11-01 Wenjie Liu -- opt-in NH4 protection uses one fixed
% NN-by-1 nodal budget [ug N/cm2/s] before the shared N/O2 reaction is applied.
% Change: 2026-11-03 Wenjie Liu -- retain the linearized nitrification sink
% separately when its magnitude is below the precision of total oxygen demand.
    n=built.domain.NL; nn=built.domain.NN; dz=built.domain.dz(:);
    assert(isstruct(schedule)&&schedule.schema_version==1 ...
        &&schedule.tile_count==numel(schedule.tiles)&&schedule.tile_count>0, ...
        'replayOxygenInterval:missingSchedule','Accepted O2 carrier schedule is required.');
    assert(iscolumn(respiratoryC)&&numel(respiratoryC)==n ...
        &&all(isfinite(respiratoryC))&&all(respiratoryC>=0), ...
        'replayOxygenInterval:respiration','Same-hour respiratory C must be a nonnegative NL column.');
    if nargin<4||isempty(nitrogenCap), nitrogenCap=inf(n,2); end
    assert(isequal(size(nitrogenCap),[n,2])&&all(nitrogenCap(:)>=0), ...
        'replayOxygenInterval:nitrogenCap','Nitrogen cap must be nonnegative NL-by-2.');
    if nargin<5,nodeBudget=[];end
    protection=[];
    if isfield(built.parameters.soil_n,'nh4_threshold_enabled') ...
            &&built.parameters.soil_n.nh4_threshold_enabled
        protection=built.parameters.soil_n.M_min_nit_ug_cm3;
        assert(iscolumn(nodeBudget)&&numel(nodeBudget)==nn ...
            &&all(isfinite(nodeBudget))&&all(nodeBudget>=0), ...
            'replayOxygenInterval:protectionBudget', ...
            'Enabled NH4 protection requires the fixed hourly nodal budget.');
    else
        assert(isempty(nodeBudget),'replayOxygenInterval:unexpectedProtectionBudget');
    end
    if isfield(schedule,'oxygen_coupling_enabled')&&~schedule.oxygen_coupling_enabled
        oxygen=neutralHour(built,schedule,nitrogenCap,protection,nodeBudget);
        return;
    end
    y=built.Y.soil_gas.y_O2(:);
    assert(numel(y)==nn&&all(isfinite(y))&&all(y>=0&y<=1), ...
        'replayOxygenInterval:oldState','Accepted dry-mole-fraction state is required.');
    accumulator=soil.gas.initOxygenAccumulator(n,nn);
    sinkIntegral=zeros(n,1); requestedIntegral=zeros(n,1);
    respiratoryActualIntegral=zeros(n,1); nitrogenActualIntegral=zeros(n,1);
    nitMMIntegral=zeros(n,2); nitMIntegral=zeros(n,2); nitTemperatureIntegral=zeros(n,2);
    nitMoistureIntegral=zeros(n,2);
    protectedReduction=zeros(n,2);
    for i=1:schedule.tile_count
        driver=schedule.tiles{i};
        if i>1
            previous=schedule.tiles{i-1};
            assert(abs(driver.t_start_s-previous.t_start_s-previous.dt_s)<=1e-8 ...
                &&samePhase(driver.phase_old,previous.phase_new), ...
                'replayOxygenInterval:discontinuousCarrier','Accepted carrier tiles must be contiguous.');
        end
        [tile,failure]=solveTile(built,driver,y,respiratoryC,nitrogenCap,protection,nodeBudget);
        if ~isempty(failure)
            oxygen=struct('has_data',false,'feasible',false,'failure_reason',failure, ...
                'failed_tile',i,'elapsed_s',accumulator.elapsed_s,'failed_tile_diagnostic',tile);
            return;
        end
        accumulator=soil.gas.accumulateOxygenSubstep(accumulator,tile,driver.dt_s);
        sinkIntegral=sinkIntegral+tile.actual_O2_sink*driver.dt_s;
        requestedIntegral=requestedIntegral+tile.R_O2_sink*driver.dt_s;
        respiratoryActualIntegral=respiratoryActualIntegral+tile.respiratory_actual_sink*driver.dt_s;
        nitrogenActualIntegral=nitrogenActualIntegral+tile.nitrification_actual_sink*driver.dt_s;
        nitMMIntegral=nitMMIntegral+tile.nit_mm_factor*driver.dt_s;
        nitMIntegral=nitMIntegral+tile.nit_M_total*driver.dt_s;
        nitTemperatureIntegral=nitTemperatureIntegral+tile.eT_nit*driver.dt_s;
        nitMoistureIntegral=nitMoistureIntegral+tile.nit_moisture_factor*driver.dt_s;
        if ~isempty(protection)
            protectedReduction=protectedReduction+tile.nh4_protection.rate_reduction_ug_cm3_s*driver.dt_s;
        end
        y=tile.y_O2;
    end
    assert(abs(accumulator.elapsed_s-schedule.elapsed_s)<=1e-9*max(1,schedule.elapsed_s), ...
        'replayOxygenInterval:duration','O2 replay must cover the exact accepted interval.');
    oxygen=soil.gas.finalizeOxygenAccumulator(accumulator);
    oxygen.feasible=true; oxygen.failure_reason=''; oxygen.oxygen_coupling_enabled=true;
    oxygen.oxygen_respiration_coupled=true;
    oxygen.nit_mm_factor=nitMMIntegral/accumulator.elapsed_s;
    oxygen.nit_M_total=nitMIntegral/accumulator.elapsed_s;
    oxygen.eT_nit=nitTemperatureIntegral/accumulator.elapsed_s;
    oxygen.nit_moisture_factor=nitMoistureIntegral/accumulator.elapsed_s;
    oxygen.nitrification_diagnostic_aggregation='accepted_tile_dt_mean';
    oxygen.y_O2=y; oxygen.final_y_O2=y;
    phase=driver.phase_new;
    oxygen.C_O2=soil.gas.species.dryMoleFractionConcentration( ...
        phase.gas_pressure_Pa,phase.temperature_C,y,phase.vapor_pressure_Pa);
    oxygen.R_O2_sink=sinkIntegral/accumulator.elapsed_s;
    oxygen.ledger.consumed_mol_cm2=oxygen.ledger.reaction_loss_mol_cm2;
    oxygen.ledger.actual_sink_mol_cm3_s=oxygen.R_O2_sink;
    oxygen.ledger.requested_sink_mol_cm3_s=requestedIntegral/accumulator.elapsed_s;
    oxygen.ledger.phase_new=struct('theta_liquid',phase.theta_l, ...
        'theta_gas',phase.theta_g,'temperature_C',phase.temperature_C);
    oxygen.phase_new=oxygen.ledger.phase_new;
    oxygen.ledger.original_phase_new=phase;
    oxygen.ledger.respiratory_demand_mol_cm2=sum(respiratoryC.*dz)*1e-6/12*accumulator.elapsed_s;
    oxygen.ledger.respiratory_actual_by_layer_mol_cm3=respiratoryActualIntegral;
    oxygen.ledger.nitrification_actual_by_layer_mol_cm3=nitrogenActualIntegral;
    oxygen.ledger.respiratory_actual_mol_cm2=sum(dz.*respiratoryActualIntegral);
    oxygen.ledger.nitrification_actual_mol_cm2=sum(dz.*nitrogenActualIntegral);
    oxygen.ledger.nitrification_demand_mol_cm2=sum(dz.*mean(oxygen.Snit,2))*1e-6/14*2*accumulator.elapsed_s;
    oxygen.ledger.demand_residual_mol_cm2=oxygen.ledger.reaction_loss_mol_cm2 ...
        -oxygen.ledger.respiratory_demand_mol_cm2-oxygen.ledger.nitrification_demand_mol_cm2;
    if ~isempty(protection)
        oxygen.nh4_protection=struct('node_budget_ug_cm2_s',nodeBudget, ...
            'rate_reduction_ug_cm3_s',protectedReduction/accumulator.elapsed_s);
    end
end

function yes=samePhase(a,b)
% Rebuilt gas fraction can differ by floating-point subtraction roundoff.
    yes=true;
    for key=fieldnames(a)'
        x=a.(key{1}); z=b.(key{1});
        yes=yes&&isequal(size(x),size(z)) ...
            &&all(abs(x(:)-z(:))<=1e-12*max(1,max(abs([x(:);z(:)]))));
    end
end

function [tile,failure]=solveTile(built,driver,yOld,respiratoryC,cap,protection,nodeBudget)
% Original fixed-old-state Picard, with a separately named respiratory sink.
    n=built.domain.NL; dz=built.domain.dz(:); p=built.parameters.soil_n;
    phResponseEnabled=~isfield(p,'ph_response_enabled')||p.ph_response_enabled;
    phys=equations.physicalConstants();
    snit=zeros(n,2);
    if isfield(built.Y,'soil_nitrogen')&&isfield(built.Y.soil_nitrogen,'variables') ...
            &&isfield(built.Y.soil_nitrogen.variables,'Snit') ...
            &&~isempty(built.Y.soil_nitrogen.variables.Snit)
        snit=min(built.Y.soil_nitrogen.variables.Snit,cap);
    end
    if ~isempty(protection)
        snit=soil.nitrogen.limitNitrificationByNodeBudget(snit,nodeBudget,dz);
    end
    reference=yOld;
    [sinkReference,nitrogenReference]=soil.gas.calcO2Sink(snit,respiratoryC,n);
    converged=false; failure=''; tile=struct();
    for iteration=1:30
        linearSink=sinkReference;
        oldMean=mean([reference(1:n),reference(2:n+1)],2);
        if any(oldMean<=0 & sinkReference>0)
            failure='positive_respiration_without_oxygen'; return;
        end
        [candidate,dissolved,aerobic,anaerobic,ledger]=soil.gas.solveO2Transport( ...
            yOld,reference,driver.phase_old,driver.phase_new,driver.transport, ...
            sinkReference,dz,driver.dt_s,driver.boundary_applied);
        [snit,~,~,nitMM,nitM,nitTemperature,~,nitMoisture]=soil.nitrogen.chemistry.nitrification( ...
            driver.phase_new.temperature_C,driver.nh4,n,driver.phase_new.theta_l, ...
            p.Theta_f,built.parameters.soil.theta_s(:),p.RHOKG,p.Vmax_nit,p.Km_nit, ...
            zeros(n,2),phys.R,phys.Tr,ones(n,2),0,p.pH,zeros(n,2),1,dz,aerobic,phResponseEnabled,driver.matric_head_cm,protection);
        snit=min(snit,cap);
        if ~isempty(protection)
            [snit,protectionAudit]=soil.nitrogen.limitNitrificationByNodeBudget(snit,nodeBudget,dz);
        end
        snitN2O=0.02*snit;
        [sinkCandidate,nitrogenCandidate]=soil.gas.calcO2Sink(snit,respiratoryC,n);
        yScale=max([1e-12;abs(candidate);abs(reference)]);
        sinkScale=max([1e-30;abs(sinkCandidate);abs(sinkReference)]);
        yRelative=max(abs(candidate-reference))/yScale;
        sinkRelative=max(abs(sinkCandidate-sinkReference))/sinkScale;
        newMean=mean([candidate(1:n),candidate(2:n+1)],2);
        actual=zeros(n,1); active=oldMean>0;
        actual(active)=linearSink(active).*newMean(active)./oldMean(active);
        % Layer-local agreement prevents positive/negative errors cancelling
        % in the column ledger. Both reactions share this original implicit
        % linearization, so report their separate actually applied amounts.
        ratio=zeros(n,1);ratio(active)=newMean(active)./oldMean(active);
        respiratoryActual=respiratoryC*1e-6/12.*ratio;
        nitrogenActual=nitrogenReference.*ratio;
        demandRelative=max(abs(actual-sinkCandidate)./max(sinkCandidate,1e-24));
        respiratoryRelative=max(abs(respiratoryActual-respiratoryC*1e-6/12) ...
            ./max(respiratoryC*1e-6/12,1e-24));
        nitrogenRequested=mean(snit,2)*1e-6/14*2;
        nitrogenRelative=max(abs(nitrogenActual-nitrogenRequested)./max(nitrogenRequested,1e-24));
        reference=candidate; sinkReference=sinkCandidate; nitrogenReference=nitrogenCandidate;
        if yRelative<=1e-6&&sinkRelative<=1e-6&&demandRelative<=1e-6 ...
                &&respiratoryRelative<=1e-6&&nitrogenRelative<=1e-6
            converged=true; break;
        end
    end
    consistency=carrierConsistency(driver.phase_new,candidate, ...
        driver.dry_air_storage_after,ledger.reaction_loss_mol_cm2,dz);
    admissible=consistency.minimum_background_storage_g_cm3>=-1e-15 ...
        &&abs(ledger.balance_residual_mol_cm2)<=1e-12+1e-8*max(ledger.storage_after_mol_cm2,1e-20);
    tile=struct('converged',converged,'admissible',admissible,'iterations',iteration, ...
        'y_relative_change',yRelative,'sink_relative_change',sinkRelative, ...
        'respiration_relative_mismatch',respiratoryRelative, ...
        'nitrification_relative_mismatch',nitrogenRelative, ...
        'demand_relative_mismatch',demandRelative,'y_O2',candidate,'C_O2_liq',dissolved, ...
        'f_aerobic',aerobic,'f_anaerobic',anaerobic,'Snit',snit,'Snit_N2O',snitN2O, ...
        'nit_mm_factor',nitMM,'nit_M_total',nitM,'eT_nit',nitTemperature, ...
        'nit_moisture_factor',nitMoisture, ...
        'R_O2_sink',sinkCandidate,'actual_O2_sink',actual, ...
        'respiratory_actual_sink',respiratoryActual,'nitrification_actual_sink',nitrogenActual, ...
        'QL',driver.QL,'QL_h',driver.QL_h,'QL_T',driver.QL_T,'QL_a',driver.QL_a, ...
        'ledger',ledger,'boundary_request',driver.boundary_request, ...
        'boundary_applied',driver.boundary_applied,'carrier_consistency',consistency);
    if ~converged
        failure='same_hour_oxygen_demand_not_converged';
    elseif ~admissible
        failure='oxygen_carrier_or_mass_balance_inadmissible';
    end
    if ~isempty(protection),tile.nh4_protection=protectionAudit;end
end

function oxygen=neutralHour(built,schedule,cap,protection,nodeBudget)
% Original accumulateWithoutAir policy: explicit neutral kinetics, no O2 ledger.
% Source: COMPASS-COMPASS_MATLAB soil.nitrogen.accumulateWithoutAir.
    n=built.domain.NL;p=built.parameters.soil_n;phys=equations.physicalConstants();
    phResponseEnabled=~isfield(p,'ph_response_enabled')||p.ph_response_enabled;
    neutral=ones(n,2); elapsed=0;
    protectedReduction=zeros(n,2);
    oxygen=struct('has_data',true,'feasible',true,'failure_reason','', ...
        'oxygen_coupling_enabled',false,'oxygen_respiration_coupled',false, ...
        'f_aerobic',neutral,'f_anaerobic',neutral,'Snit',zeros(n,2), ...
        'nit_mm_factor',zeros(n,2),'nit_M_total',zeros(n,2),'eT_nit',zeros(n,2), ...
        'nit_moisture_factor',zeros(n,2), ...
        'nitrification_diagnostic_aggregation','accepted_tile_dt_mean', ...
        'Snit_N2O',zeros(n,2),'QL',zeros(n,1),'QL_abs',zeros(n,1), ...
        'QL_h',zeros(n,1),'QL_T',zeros(n,1),'QL_a',zeros(n,1));
    for i=1:schedule.tile_count
        d=schedule.tiles{i};
        [snit,~,~,nitMM,nitM,nitTemperature,~,nitMoisture]=soil.nitrogen.chemistry.nitrification( ...
            d.temperature_C,d.nh4,n,d.theta_liquid,p.Theta_f,built.parameters.soil.theta_s(:), ...
            p.RHOKG,p.Vmax_nit,p.Km_nit,zeros(n,2),phys.R,phys.Tr,neutral,0,p.pH,zeros(n,2), ...
            1,built.domain.dz,neutral,phResponseEnabled,d.matric_head_cm,protection);
        snit=min(snit,cap);
        if ~isempty(protection)
            [snit,protectionAudit]=soil.nitrogen.limitNitrificationByNodeBudget(snit,nodeBudget,built.domain.dz);
            protectedReduction=protectedReduction+protectionAudit.rate_reduction_ug_cm3_s*d.dt_s;
        end
        oxygen.Snit=oxygen.Snit+snit*d.dt_s;
        oxygen.Snit_N2O=oxygen.Snit_N2O+0.02*snit*d.dt_s;
        oxygen.nit_mm_factor=oxygen.nit_mm_factor+nitMM*d.dt_s;
        oxygen.nit_M_total=oxygen.nit_M_total+nitM*d.dt_s;
        oxygen.eT_nit=oxygen.eT_nit+nitTemperature*d.dt_s;
        oxygen.nit_moisture_factor=oxygen.nit_moisture_factor+nitMoisture*d.dt_s;
        for key={'QL','QL_h','QL_T','QL_a'}
            k=key{1}; oxygen.(k)=oxygen.(k)+d.(k)*d.dt_s;
        end
        oxygen.QL_abs=oxygen.QL_abs+abs(d.QL)*d.dt_s;
        elapsed=elapsed+d.dt_s;
    end
    assert(abs(elapsed-schedule.elapsed_s)<=1e-9*max(1,elapsed));
    for key={'Snit','Snit_N2O','QL','QL_abs','QL_h','QL_T','QL_a','nit_mm_factor','nit_M_total','eT_nit','nit_moisture_factor'}
        k=key{1};oxygen.(k)=oxygen.(k)/elapsed;
    end
    oxygen.elapsed_s=elapsed;oxygen.accepted_tiles=schedule.tile_count;
    if ~isempty(protection)
        oxygen.nh4_protection=struct('node_budget_ug_cm2_s',nodeBudget, ...
            'rate_reduction_ug_cm3_s',protectedReduction/elapsed);
    end
    oxygen.final_y_O2=built.Y.soil_gas.y_O2;oxygen.y_O2=oxygen.final_y_O2;
    phase=soil.gas.captureOxygenPhase(built,true);
    oxygen.phase_new=struct('theta_liquid',phase.theta_l, ...
        'theta_gas',phase.theta_g,'temperature_C',phase.temperature_C);
    oxygen.C_O2=soil.gas.species.dryMoleFractionConcentration( ...
        phase.gas_pressure_Pa,phase.temperature_C,oxygen.y_O2,phase.vapor_pressure_Pa);
end

function consistency=carrierConsistency(phase,y,dryAirStorage,reactionMol,dz)
% Original dilute-species carrier check, not feedback into h/U/P_g.
    p=soil.gas.species.oxygenProperties(); n=numel(dz);
    T=[phase.temperature_C(1:n),phase.temperature_C(2:n+1)]+273.15;
    Pg=[phase.gas_pressure_Pa(1:n),phase.gas_pressure_Pa(2:n+1)];
    Pv=[phase.vapor_pressure_Pa(1:n),phase.vapor_pressure_Pa(2:n+1)];
    cd=(Pg-Pv)./(8.31446261815324.*T)*1e-6;
    KH=p.henry_reference_mol_L_atm.*exp(p.henry_enthalpy_over_R_K.*(1./T-1/p.henry_reference_temperature_K));
    beta=KH.*p.gas_constant_cm3_atm_mol_K.*T/1000;
    storage=(phase.theta_g.*cd+phase.theta_l.*beta.*cd).*[y(1:n),y(2:n+1)]*p.species_molar_mass_g_mol;
    background=dryAirStorage-storage;
    column=sum(dz.*mean(dryAirStorage,2));
    consistency=struct('minimum_background_storage_g_cm3',min(background(:)), ...
        'omitted_reaction_feedback_g_cm2',reactionMol*p.species_molar_mass_g_mol, ...
        'omitted_reaction_feedback_ratio',reactionMol*p.species_molar_mass_g_mol/max(column,realmin));
end
