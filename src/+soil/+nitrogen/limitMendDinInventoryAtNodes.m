function [remaining,audit]=limitMendDinInventoryAtNodes( ...
        remaining,root,nh4,no3,theta,rhoKd,bd,dz,dt,S_min_nh4_mg_kg)
%LIMITMENDDININVENTORYATNODES Bound uniform layer uptake by shared-node DIN.
% Physical post-root inventories REMAINING remain top-down mg N/g dry soil.
% New *_accessible_mg_g fields are allocation ceilings, not changed states
% or kinetic concentrations. The uniform layer sink uses the ADE's consistent
% linear FE source quadrature; node stocks use lumped liquid-plus-sorbed NH4
% or liquid NO3. ROOT supplies bottom-up NL-by-2 ug N/cm3 soil/s; NH4/NO3 are
% bottom-up NN concentrations [ug N/cm3 water]; THETA is NL-by-2; BD g/cm3,
% DZ cm and DT seconds. Each layer receives the minimum of its two incident
% node allowance factors. The combined assembled demand then cannot exceed
% nonnegative post-root nodal storage. Independent root overdraw stays visible
% in the audit and remains subject to the unchanged transport gates.
% MEND must use ceilings in both uptake stages and rebuild C/N feedback;
% this function never truncates an already committed organic sink.
% Optional S_MIN_NH4_MG_KG is a finite nonnegative bottom-up NL-vector of
% protected NH4 [mg N/kg dry soil]. Empty or scalar/vector zero retains the
% original arithmetic. Nonzero protection first reduces the post-root layer
% budget by 1e-3*flip(S_MIN_NH4_MG_KG) [mg/g], then reserves the matching
% half-element projected node masses in the shared-node FE limiter. Physical
% inventories and NO3 remain unchanged; same-hour internal release may be
% recaptured, so the ceiling bounds old-inventory contribution, not gross
% immobilization including that release.
% Source: limitNitrificationAtNodes uses the identical ADE quadrature.
% Editor: Wenjie Liu; Date: 2026-10-12.
% Change: 2026-11-01 Wenjie Liu -- combine optional layer protection with
% projected shared-node NH4 reservations without changing physical stocks.
    bd=bd(:);nl=numel(bd);
    assert(all(isfinite(bd))&&all(bd>0),'soil:nitrogen:mendNodeDensity');
    thresholdEnabled=false;
    if nargin>=10&&~isempty(S_min_nh4_mg_kg)
        assert(isnumeric(S_min_nh4_mg_kg)&&isreal(S_min_nh4_mg_kg) ...
            &&isvector(S_min_nh4_mg_kg) ...
            &&(numel(S_min_nh4_mg_kg)==nl ...
                ||(isscalar(S_min_nh4_mg_kg)&&S_min_nh4_mg_kg==0)) ...
            &&all(isfinite(S_min_nh4_mg_kg(:)))&&all(S_min_nh4_mg_kg(:)>=0), ...
            'soil:nitrogen:mendNodeProtection', ...
            'S_min_nh4_mg_kg must be a finite nonnegative NL-vector or zero scalar.');
        thresholdEnabled=any(S_min_nh4_mg_kg(:)>0);
    end
    if thresholdEnabled
        S_min_nh4_mg_kg=S_min_nh4_mg_kg(:);
        protectedLayerMgG=1e-3*flip(S_min_nh4_mg_kg);
        protectedHalfMass=bd.*S_min_nh4_mg_kg.*dz(:)/2;
        protectedNodeMass_ug_cm2=[protectedHalfMass;0]+[0;protectedHalfMass];
    end
    species={'nh4','no3'};concentrations={nh4,no3};storage={rhoKd,zeros(nl,1)};
    audit=struct('policy','shared_node_FE_accessible_inventory', ...
        'physical_inventory_unchanged',true,'kinetic_concentration_unchanged',true);
    for j=1:2
        key=species{j};budget=remaining.([key '_mg_g'])(:);
        assert(numel(budget)==nl&&all(isfinite(budget))&&all(budget>=0), ...
            'soil:nitrogen:mendNodeInventory','Invalid physical post-root DIN.');
        if j==1&&thresholdEnabled
            budget=max(budget-protectedLayerMgG,0);
        end
        rate=repmat(flip(budget).*bd.*1000./dt,1,2);
        rootRate=zeros(nl,2);
        if root.enabled,rootRate=root.(['sink_' key]);end
        if j==1&&thresholdEnabled
            [~,cap]=soil.nitrogen.limitNitrificationAtNodes( ...
                rate,-rootRate,concentrations{j},theta,storage{j},dz,dt, ...
                protectedNodeMass_ug_cm2);
        else
            [~,cap]=soil.nitrogen.limitNitrificationAtNodes( ...
                rate,-rootRate,concentrations{j},theta,storage{j},dz,dt);
        end
        remaining.([key '_accessible_mg_g'])=budget.*flip(cap.element_scale);
        audit.(key)=cap;
    end
end
