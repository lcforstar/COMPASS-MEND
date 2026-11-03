function storage=columnStorage(state,thetaLiquid,thetaIce,rhokg,dz)
%COLUMNSTORAGE Integrate the non-overlapping fertilizer-label owners.
%
% Mineral node pools are converted from ug N cm-2 to g N cm-2. Organic and
% DON element pools already use g N cm-3 bulk. Derived aggregates and
% microbial phase subpools are intentionally excluded from the sum.
% Editor: Wenjie Liu; Date: 2026-09-21.
% Change: 2026-10-30 Wenjie Liu -- retain manure-source surface and element
% AOM as separate owners without adding caches twice.
    dz=dz(:);nl=numel(dz);
    [din4,din3]=soil.nitrogen.transport.calculateStorage( ...
        thetaLiquid,state.nh4,state.no3,rhokg,dz,nl,thetaIce, ...
        state.nh4_ice,state.no3_ice,state.nh4_granule);
    capacity=[dz;0]./2+[0;dz]./2;
    urea=capacity'*state.urea(:);
    fields={'PON1','PON2','MON','MB_N','BD_N', ...
        'ENZP1_N','ENZP2_N','EM_N','ENZN_N'};
    profile=state.dom.don_total_g_cm3(:);
    for i=1:numel(fields)
        profile=profile+state.organic_n.(fields{i})(:);
    end
    organic=sum(profile.*dz);
    storage=struct('nh4_g_cm2',din4.*1e-6, ...
        'no3_g_cm2',din3.*1e-6,'urea_g_cm2',urea.*1e-6, ...
        'organic_g_cm2',organic, ...
        'total_g_cm2',(din4+din3+urea).*1e-6+organic);
    if isfield(state,'aom')
        assert(isfield(state.aom,'manure'), ...
            'soil:nitrogen:tracer:aomStorage','Manure AOM label owner is missing.');
        aom=state.aom.manure;
        storage.aom_g_cm2=aom.surface_N+sum(aom.element_N(:).*dz);
        storage.total_g_cm2=storage.total_g_cm2+storage.aom_g_cm2;
    end
end
