function [row,cumulative]=acceptedOrganicAudit(built,previousY,initialY,previousDIN,kt,startFresh,firstStep)
%ACCEPTEDORGANICAUDIT Read-only interval/cumulative C/N, DIN and O2 gates.
% Inputs: accepted Y endpoints, previous DIN residual [ug/cm2], KT and fresh
% flags. Row: KT,C,N[g/cm2],DIN interval/cumulative[ug/cm2],O2[mol/cm2],
% SOC/free DOC [mg/kg,0--20 cm]. Gates: 1e-10 C,1e-12 N,1e-6 DIN,1e-14 O2.
% No state/physical ledger mutation. SOC uses historical independent solid
% pools plus DOM/IOM/NEC over dry mass; QOC is excluded from that observable.
% Editor: Wenjie Liu; Date: 2026-10-17.
% Change: 2026-10-23 Wenjie Liu -- report O2 as not applicable (NaN) when
%   the reference air-off mode is selected; retain all enabled-process gates.
    options=struct('carbon_tolerance_g_cm2',1e-10,'nitrogen_tolerance_g_cm2',1e-12, ...
        'start_is_fresh',startFresh&&firstStep);
    interval=diagnostics.biogeochemistryBudget(built,previousY,built.Y,options);
    options.start_is_fresh=startFresh;
    cumulative=diagnostics.biogeochemistryBudget(built,initialY,built.Y,options);
    assert(interval.passed&&cumulative.passed,'acceptedAudit:CN','C/N conservation failed at KT=%d.',kt);
    cn=[interval.carbon_physical_residual_g_cm2,interval.nitrogen_residual_g_cm2];
    n=built.workspace.nitrogen;din=[n.cum_residual_nh4,n.cum_residual_no3];change=din-previousDIN;
    oxygenEnabled=~isfield(n,'oxygen_coupling_enabled')||n.oxygen_coupling_enabled;
    o2=NaN;oxygenPassed=true;
    if oxygenEnabled
        o2=n.oxygen.balance_residual_mol_cm2;
        oxygenPassed=isfinite(o2)&&abs(o2)<1e-14;
    end
    assert(all(isfinite([cn din change]))&&all(abs(din)<1e-6) ...
        &&all(abs(change)<1e-6)&&oxygenPassed,'acceptedAudit:DINO2','DIN/O2 conservation failed at KT=%d.',kt);
    c=zeros(built.domain.NL,1);
    for key={'POM1','POM2','MOC','MB','BD','ENZP1','ENZP2','EM','ENZN'}
        c=c+built.Y.soil_microbes.(key{1});
    end
    rho=built.parameters.soil.bulk_density(:);z=built.domain.z_node(:);
    mass=flip(c).*rho/1000+built.Y.soil_dom.doc_total_g_cm3 ...
        +built.Y.soil_organic.IOM_C+built.Y.soil_organic.NEC_C;
    width=max(0,min(z(1:end-1),20)-max(z(2:end),0));
    soc=1e6*sum(mass.*width)/sum(rho.*width);
    row=[kt cn change din o2 soc diagnostics.docFreeLayerMean(built,0,20)];
end
