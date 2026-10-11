function [C_O2,C_O2_liq,f_aerobic,f_anaerobic,ledger]=solveO2Transport( ...
        C_old,gasNew,porosity,Tnew,velocity,requested,dz,dt,NL,NN,top,bottom,phase)
%SOLVEO2TRANSPORT Conservative positive gas/dissolved oxygen interval.
%   Bottom-up node C [mol/cm3 gas], endpoint liquid/gas fractions [NL,2],
%   requested layer sink [mol/cm3 bulk/s], velocity [cm/s] upward positive.
%   PHASE supplies theta_liquid_new, theta_liquid_old, theta_gas_old, and
%   temperature_old [NN,1]. Optional reserved_sink_mol_cm3_s gets priority.
%   Without PHASE only fixed, unfrozen phase tests are supported. Production
%   always supplies actual phase history. Top is Dirichlet; bottom defaults
%   to closed. Lumped Mnew*Cnew-Mold*Cold, existing molecular diffusion and
%   donor-upwind gas advection form a positive operator. An active set limits
%   demand and returns actual removal, unmet demand and boundary mass.
%   An infeasible reserved demand fails before commit: caller rebuilds C/N.
%   Source: retained COMPASS Henry/diffusion constants; repair contract at
%   docs/soil/biogeochemistry/OXYGEN_CONSERVATIVE_REPAIR_20260914.md.
% Editor: Jizixing Zhu; Date: 2026-09-14.
% Change: 2026-09-14 Codex -- conserve phase-resolved oxygen storage and
%   return applied consumption instead of clipping material negative states.
% Change: 2026-10-05 Wenjie Liu -- retain successful active-set solutions
%   exactly; if node-to-layer demand mapping produces a material negative,
%   preserve feasible reserved respiration and conservatively limit only
%   the remaining layer demand along its common nonnegative feasible ray.
    if nargin<12, bottom=[]; end
    C_old=C_old(:); Tnew=Tnew(:); dz=dz(:); requested=requested(:);
    porosity=porosity(:); velocity=velocity(:);
    if nargin<13||isempty(phase)
        liquid=max(porosity-gasNew,0);
        phase=struct('theta_liquid_new',liquid,'theta_liquid_old',liquid, ...
            'theta_gas_old',gasNew,'temperature_old',Tnew);
    end
    oldGas=phase.theta_gas_old; oldLiquid=phase.theta_liquid_old;
    liquid=phase.theta_liquid_new; Told=phase.temperature_old(:);
    assert(NN==NL+1&&isequal(size(gasNew),[NL,2]) ...
        &&isequal(size(oldGas),[NL,2])&&isequal(size(liquid),[NL,2]) ...
        &&isequal(size(oldLiquid),[NL,2])&&numel(C_old)==NN ...
        &&numel(Tnew)==NN&&numel(Told)==NN&&numel(dz)==NL ...
        &&numel(requested)==NL&&numel(velocity)==NL&&numel(porosity)==NL, ...
        'soil:gas:oxygenShape','Oxygen phase and transport arrays must align.');
    assert(all(isfinite([C_old;gasNew(:);oldGas(:);liquid(:);oldLiquid(:); ...
        Tnew;Told;dz;requested;velocity;porosity])) ...
        &&all([C_old;gasNew(:);oldGas(:);liquid(:);oldLiquid(:);requested]>=0) ...
        &&all(dz>0)&&all(porosity>0)&&all([Tnew;Told]>-273.15) ...
        &&isscalar(dt)&&isfinite(dt)&&dt>0&&isscalar(top)&&isfinite(top)&&top>=0, ...
        'soil:gas:oxygenValues','Oxygen inputs must be finite and physical.');
    boundary=NN; boundaryValue=top;
    if ~isempty(bottom)&&~isnan(bottom)
        assert(isscalar(bottom)&&isfinite(bottom)&&bottom>=0);
        boundary=[1;NN]; boundaryValue=[bottom;top];
    end
    newCapacity=gasNew+henryRatio(Tnew).*liquid;
    oldCapacity=oldGas+henryRatio(Told).*oldLiquid;
    Mnew=accumarray([(1:NL)';(2:NN)'],[newCapacity(:,1);newCapacity(:,2)].*[dz;dz]/2,[NN,1]);
    Mold=accumarray([(1:NL)';(2:NN)'],[oldCapacity(:,1);oldCapacity(:,2)].*[dz;dz]/2,[NN,1]);
    assert(all(Mnew>0),'soil:gas:oxygenCapacity','Oxygen capacity must be positive.');
    tk=[Tnew(1:NL),Tnew(2:NN)]+273.15;
    diffusivity=0.205*(tk/298.15).^1.75.*gasNew.^(10/3)./porosity.^2;
    conductance=mean(diffusivity,2)./dz;
    up=max(velocity,0); down=min(velocity,0);
    lo=conductance+up; hi=conductance-down;
    K=sparse([(1:NL)';(1:NL)';(2:NN)';(2:NN)'], ...
        [(1:NL)';(2:NN)';(1:NL)';(2:NN)'],[lo;-hi;-lo;hi],NN,NN);
    A=spdiags(Mnew/dt,0,NN,NN)+K;
    B=sparse([(1:NL)';(2:NN)'],[(1:NL)';(1:NL)'],[dz;dz]/2,NN,NL);
    b=Mold.*C_old/dt;
    fixedA=A; fixedA(boundary,:)=0;
    fixedA(sub2ind([NN,NN],boundary,boundary))=1;
    fixedB=B; fixedB(boundary,:)=0;
    fixedb=b; fixedb(boundary)=boundaryValue;
    reserved=zeros(NL,1);
    if isfield(phase,'reserved_sink_mol_cm3_s')
        reserved=phase.reserved_sink_mol_cm3_s(:);
        assert(numel(reserved)==NL&&all(isfinite(reserved)) ...
            &&all(reserved>=0)&&all(reserved<=requested+1e-25), ...
            'soil:gas:oxygenReservation','Reserved oxygen must be within demand.');
    end
    base=fixedb-fixedB*reserved;
    passive=fixedA\base;
    tolerance=1e-12*max(1e-12,max(abs(passive)));
    assert(all(passive>=-tolerance),'soil:gas:oxygenReservationInfeasible', ...
        'Organic respiration exceeds oxygen supply; rebuild its C/N candidate.');
    extra=max(requested-reserved,0);
    rhs=base-fixedB*extra;
    c=fixedA\rhs; active=false(NN,1);
    for iteration=1:2*NN+4
        residual=fixedA*c-rhs;
        next=(c< -tolerance)|(active&residual>0);
        next(boundary)=false;
        if isequal(next,active), break; end
        active=next; matrix=fixedA; target=rhs;
        rows=find(active); matrix(rows,:)=0;
        matrix(sub2ind([NN,NN],rows,rows))=1; target(rows)=0;
        c=matrix\target;
    end
    assert(iteration<2*NN+4&&all(c>=-tolerance), ...
        'soil:gas:oxygenActiveSet','Oxygen demand active set did not converge.');
    nodeDemand=fixedB*extra;
    deficit=max(fixedA*c-rhs,0);
    deficit(~active)=0;
    nodeScale=ones(NN,1); has=nodeDemand>0;
    nodeScale(has)=max(0,min(1,1-deficit(has)./nodeDemand(has)));
    nodeScale(boundary)=1;
    scale=min(nodeScale(1:NL),nodeScale(2:NN));
    actual=reserved+scale.*extra;
    C_O2=fixedA\(fixedb-fixedB*actual);
    feasibilityFallback=[];
    if all(isfinite(C_O2))&&any(C_O2< -tolerance)
        [C_O2,actual,scale,feasibilityFallback]=limitAdditionalDemand( ...
            fixedA,fixedB,fixedb,reserved,extra,scale,passive,tolerance,C_O2);
    end
    assert(all(isfinite(C_O2))&&all(C_O2>=-tolerance), ...
        'soil:gas:oxygenNegative','Applied oxygen demand produced a negative state.');
    % Roundoff-only correction remains visible in the final mass residual.
    C_O2(C_O2<0)=0;
    reaction=B*actual;
    rowResidual=A*C_O2-b+reaction;
    boundaryMass=zeros(NN,1); boundaryMass(boundary)=dt*rowResidual(boundary);
    storageBefore=sum(Mold.*C_old); storageAfter=sum(Mnew.*C_O2);
    consumed=dt*sum(actual.*dz);
    balance=storageAfter-storageBefore+consumed-sum(boundaryMass);
    assert(abs(balance)<=1e-10*max(1e-12,storageBefore+storageAfter+consumed ...
        +sum(abs(boundaryMass))), 'soil:gas:oxygenClosure','Oxygen mass balance failed.');
    avgTk=mean(tk,2); KH=1.3e-3*exp(1700*(1./avgTk-1/298.15));
    C_O2_liq=KH.*mean([C_O2(1:NL),C_O2(2:NN)],2).*82.06.*avgTk;
    aerobic=C_O2_liq./(C_O2_liq+20e-6);
    f_aerobic=[aerobic,aerobic]; f_anaerobic=1-f_aerobic;
    ledger=struct('requested_sink_mol_cm3_s',requested, ...
        'actual_sink_mol_cm3_s',actual,'unmet_sink_mol_cm3_s',requested-actual, ...
        'additional_scale',scale,'storage_before_mol_cm2',storageBefore, ...
        'storage_after_mol_cm2',storageAfter,'consumed_mol_cm2',consumed, ...
        'top_input_mol_cm2',boundaryMass(end),'bottom_input_mol_cm2',boundaryMass(1), ...
        'balance_residual_mol_cm2',balance,'active_set_iterations',iteration, ...
        'phase_new',struct('theta_liquid',liquid,'theta_gas',gasNew,'temperature_C',Tnew));
    if ~isempty(feasibilityFallback)
        ledger.feasibility_fallback=feasibilityFallback;
    end
end

function [concentration,actual,scale,audit]=limitAdditionalDemand( ...
        operator,layerToNode,boundaryRhs,reserved,extra,originalScale,passive,tolerance,failed)
%LIMITADDITIONALDEMAND Feasible reserved stock plus a bounded extra-demand ray.
% All concentrations are node mol/cm3 gas; layer sinks are mol/cm3 bulk/s.
% The unchanged positive transport operator gives C(lambda)=passive-lambda*H,
% H=A\(B*(originalScale.*extra)). Restricting lambda to [0,1] can only reduce
% the original additional removal; reserved respiration remains exact.
% The maximal feasible scalar along this already-limited demand direction is
% min(passive_i/H_i) over H_i>0. This is not a global allocation optimum.
% Re-solve with the actual layer removal so concentration and the mass ledger
% use identical sinks. A bounded bisection is used only if roundoff at the
% analytical boundary fails the original nonnegativity tolerance. The known
% feasible zero-extra endpoint remains available; no concentration is clipped
% here, and the original material-negative and mass-balance guards still run.
% Editor: Wenjie Liu; Date: 2026-10-05.
    response=operator\(layerToNode*(originalScale.*extra));
    assert(all(isfinite(response))&&all(response>=-tolerance), ...
        'soil:gas:oxygenFeasibilityOperator','Extra-demand response must retain the positive transport operator.');
    positive=response>0;
    assert(any(positive),'soil:gas:oxygenFeasibilityDirection', ...
        'A negative mapped solution requires a finite additional-demand direction.');
    analyticCap=max(0,min(1,min(max(passive(positive),0)./response(positive))));
    lambda=analyticCap*(1-32*eps(1));
    scale=originalScale*lambda;actual=reserved+scale.*extra;
    concentration=operator\(boundaryRhs-layerToNode*actual);
    solves=1;usedBisection=false;
    if ~all(isfinite(concentration))||any(concentration< -tolerance)
        usedBisection=true;lower=0;upper=lambda;
        concentration=passive;actual=reserved;scale=zeros(size(originalScale));
        for iteration=1:40
            trialLambda=(lower+upper)/2;
            trialScale=originalScale*trialLambda;
            trialActual=reserved+trialScale.*extra;
            trialConcentration=operator\(boundaryRhs-layerToNode*trialActual);
            solves=solves+1;
            if all(isfinite(trialConcentration))&&all(trialConcentration>=-tolerance)
                lower=trialLambda;concentration=trialConcentration;
                actual=trialActual;scale=trialScale;
            else
                upper=trialLambda;
            end
        end
        lambda=lower;
    end
    assert(all(actual>=reserved)&&all(actual<=reserved+originalScale.*extra), ...
        'soil:gas:oxygenFeasibilitySink','Fallback may reduce only nonreserved additional demand.');
    audit=struct('method','common_extra_demand_feasible_ray', ...
        'trigger_min_concentration_mol_cm3_gas',min(failed), ...
        'reserved_preserved',true,'analytic_lambda_cap',analyticCap,'applied_lambda',lambda, ...
        'linear_solves',solves,'roundoff_bisection_used',usedBisection, ...
        'original_additional_scale',originalScale,'original_nonnegative_tolerance',tolerance);
end

function ratio=henryRatio(temperature)
% Dissolved/gas concentration ratio at each endpoint; ice is not a solvent.
% Editor: Jizixing Zhu; Date: 2026-09-14.
    tk=temperature(:)+273.15;
    alpha=1.3e-3*exp(1700*(1./tk-1/298.15)).*82.06.*tk/1000;
    ratio=[alpha(1:end-1),alpha(2:end)];
end
