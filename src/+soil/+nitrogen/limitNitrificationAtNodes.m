function [limited,info]=limitNitrificationAtNodes(rate,otherRate,cOld,thetaOld,rhoKd,dz,dt,protectedNodeMass_ug_cm2)
%LIMITNITRIFICATIONATNODES Bound NH4 oxidation by the assembled nodal inventory.
% RATE and OTHER_RATE are bottom-up [NL x 2] endpoint rates
% [ug N cm-3 soil s-1]. OTHER_RATE contains all unchanged NH4 exchanges:
% Shys+Smin+Sdis-Svol-root_NH4. COLD is the post-event, step-start nodal
% dissolved concentration [NN x 1, ug N cm-3 water]; THETAOLD [NL x 2]
% and RHOKD [NL x 1] define the same liquid-plus-adsorbed lumped storage
% used by transport.assembleMatrixCoefficients. DZ is [NL x 1, cm]; DT [s].
% Optional PROTECTEDNODEMASS_UG_CM2 is a finite nonnegative NN-vector of
% bottom-up protected NH4 area masses [ug N/cm2]. Empty or scalar/vector zero
% executes the original availability arithmetic. A nonzero threshold is
% subtracted BEFORE adding OTHER_RATE, so positive sources first fill any
% initial protection deficit. Independent root/transport/volatilization
% exchanges retain their physical equations and can cross the threshold.
%
% Layer inventory rationing and the oxygen cap must run first. Consistent
% linear FE source quadrature can otherwise withdraw more N from one node
% than that node holds, even when its incident layers have sufficient N.
% Each element receives the minimum allowance at its two nodes, scaling
% both endpoints together. No accepted concentration, external exchange,
% boundary or transport coefficient is changed. Unused nitrification is
% returned in INFO and must propagate to NO3, N2O, O2 and allocation ledgers.
% If OTHER_RATE alone overdraws a node, its negative allowance is retained;
% this helper cannot repair that independent inconsistency. The unchanged
% transport positivity/conservation gates still decide acceptance.
%
% Outputs: LIMITED [NL x 2] nonnegative rates no greater than RATE; INFO
% contains nodal available/requested/applied loads [ug N cm-2 s-1] and
% element_scale [NL x 1]. This is an availability cap, not a general ADE
% positivity proof. A feasible input remains bitwise unchanged.
%
% Editor: Wenjie Liu; Date: 2026-10-10.
% Change: 2026-11-01 Wenjie Liu -- optionally reserve protected node masses;
% keep the signed protection deficit and the original zero-threshold path.
    dz=dz(:);cOld=cOld(:);rhoKd=rhoKd(:);nl=numel(dz);
    assert(isequal(size(rate),[nl,2],size(otherRate),size(thetaOld)) ...
        &&numel(cOld)==nl+1&&numel(rhoKd)==nl ...
        &&isscalar(dt)&&isfinite(dt)&&dt>0 ...
        &&all(isfinite([rate(:);otherRate(:);thetaOld(:);rhoKd;cOld;dz])) ...
        &&all(rate>=0,'all')&&all(thetaOld>=0,'all') ...
        &&all(rhoKd>=0)&&all(cOld>=0)&&all(dz>0), ...
        'soil:nitrogen:nodalNitrificationInputs', ...
        'Invalid NH4 nodal reaction allocation inputs.');
    thresholdEnabled=false;
    if nargin>=8&&~isempty(protectedNodeMass_ug_cm2)
        assert(isnumeric(protectedNodeMass_ug_cm2)&&isreal(protectedNodeMass_ug_cm2) ...
            &&isvector(protectedNodeMass_ug_cm2) ...
            &&(numel(protectedNodeMass_ug_cm2)==nl+1 ...
                ||(isscalar(protectedNodeMass_ug_cm2)&&protectedNodeMass_ug_cm2==0)) ...
            &&all(isfinite(protectedNodeMass_ug_cm2(:))) ...
            &&all(protectedNodeMass_ug_cm2(:)>=0), ...
            'soil:nitrogen:nodalProtectionInputs', ...
            'Protected NH4 mass must be a finite nonnegative NN-vector or zero scalar.');
        thresholdEnabled=any(protectedNodeMass_ug_cm2(:)>0);
    end
    load=@(v)[(2*v(:,1)+v(:,2)).*dz/6;0] ...
        +[0;(v(:,1)+2*v(:,2)).*dz/6];
    stock=[(thetaOld(:,1)+rhoKd).*dz/2;0] ...
        +[0;(thetaOld(:,2)+rhoKd).*dz/2];
    if thresholdEnabled
        protectedNodeMass_ug_cm2=protectedNodeMass_ug_cm2(:);
        available=(stock.*cOld-protectedNodeMass_ug_cm2)/dt+load(otherRate);
    else
        available=stock.*cOld/dt+load(otherRate);
    end
    demand=load(rate);
    factor=ones(nl+1,1);
    active=demand>max(available,0);
    factor(active)=max(available(active),0)./demand(active)*(1-8*eps);
    element=min(factor(1:end-1),factor(2:end));
    limited=rate.*element;
    info=struct('node_available_ug_cm2_s',available, ...
        'node_requested_ug_cm2_s',demand, ...
        'node_applied_ug_cm2_s',load(limited),'element_scale',element, ...
        'storage_basis','step-start liquid plus equilibrium adsorbed NH4', ...
        'source_quadrature','consistent linear finite element');
    if thresholdEnabled
        info.protected_node_mass_ug_cm2=protectedNodeMass_ug_cm2;
    end
end
