function result=mendDinKineticRequest( ...
        dinRemaining,activeCarbon,activeNitrogen,stoichiometricDemand, ...
        parameters,dt_s)
%MENDDINKINETICREQUEST Evaluate source MEND microbial NH4/NO3 uptake demand.
%   DINREMAINING is top-down NL-by-2 *dissolved* [mg N g-1 dry soil] after
%   COMPASS root reservation. The separate allocator caps uptake by total
%   post-root NH4 (including equilibrium sorption) and liquid NO3 inventory.
%   ACTIVE pools are top-down [mg element g-1 dry soil]. The
%   MEND MM request uses VNup_MB, KsNH4_MB and KsNO3_MB on an hourly clock.
%   MEND vegetation terms are deliberately excluded: COMPASS plant radial
%   transport owns plant uptake and has already reserved its applied demand.
%   The first-stage source kinetic request is independent of the later
%   microbial C:N deficit; post-root DIN inventory limits its application.
%   Reference: MOD_MEND.F90 commit 92323c7, microbial N immobilization.
%   Editor: Jizixing Zhu; Date: 2026-09-12.
% Change: 2026-09-13 Codex -- share the source NUE closure with MEND
%   organic-N mineralization, preventing divergent C:N response shapes.
% Change: 2026-09-13 Codex -- retain the source first-stage kinetic request
%   even when the separate post-reaction C:N deficit is zero.
% Change: 2026-09-14 Codex -- document the free-NH4 source substrate,
%   distinct from the total-NH4 allocation stock.
    dinRemaining=dinRemaining(:,:); activeCarbon=activeCarbon(:);
    activeNitrogen=activeNitrogen(:); stoichiometricDemand=stoichiometricDemand(:);
    n=numel(activeCarbon);
    assert(isequal(size(dinRemaining),[n,2]) ...
        &&isequal(size(activeNitrogen),[n,1]) ...
        &&isequal(size(stoichiometricDemand),[n,1]) ...
        &&all(isfinite([dinRemaining(:);activeCarbon;activeNitrogen;stoichiometricDemand])) ...
        &&all([dinRemaining(:);activeCarbon;activeNitrogen;stoichiometricDemand]>=0) ...
        &&isscalar(dt_s)&&isfinite(dt_s)&&dt_s>0, ...
        'soil:microbe:mendDinKinetics','Invalid MEND DIN kinetic inputs.');
    v=expand(parameters.VNup_MB_h,n);
    k=[expand(parameters.KsNH4_MB_mg_g,n), ...
        expand(parameters.KsNO3_MB_mg_g,n)];
    cnMin=expand(parameters.CN_MB_min,n);
    cnMax=expand(parameters.CN_MB_max,n);
    w=parameters.wdorm;
    assert(all(v>=0)&&all(k(:)>0)&&all(cnMax>cnMin)&&all(cnMin>0) ...
        &&isscalar(w)&&isfinite(w)&&w>=1, ...
        'soil:microbe:mendDinKinetics','Invalid MEND DIN kinetic parameters.');
    phi=soil.microbe.mendNitrogenUseEfficiency( ...
        activeCarbon,activeNitrogen,cnMin,cnMax,w);
    totalDin=sum(dinRemaining,2); share=zeros(n,2); hasDin=totalDin>0;
    share(hasDin,:)=dinRemaining(hasDin,:)./totalDin(hasDin);
    competition=1+sum(dinRemaining./k,2);
    raw=zeros(n,2); dt_h=dt_s/3600;
    for species=1:2
        raw(:,species)=dt_h.*v.*share(:,species).*phi ...
            .*dinRemaining(:,species).*activeNitrogen ...
            ./(k(:,species).*(competition+activeNitrogen./k(:,species)));
    end
    raw(~isfinite(raw))=0;
    result=struct('requested_mg_g',raw, ...
        'raw_kinetic_request_mg_g',raw,'phi_n_limitation',phi, ...
        'stoichiometric_demand_mg_g',stoichiometricDemand, ...
        'plant_terms_in_denominator',false, ...
        'competition_owner','COMPASS_root_reservation_then_MEND_microbes');
end

function value=expand(value,n)
    value=value(:); if isscalar(value), value=repmat(value,n,1); end
    assert(numel(value)==n&&all(isfinite(value)), ...
        'soil:microbe:mendDinKinetics','MEND kinetic profile must align with NL.');
end
