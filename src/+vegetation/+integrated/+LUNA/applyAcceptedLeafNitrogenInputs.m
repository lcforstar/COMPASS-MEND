function model=applyAcceptedLeafNitrogenInputs(model)
%APPLYACCEPTEDLEAFNITROGENINPUTS Bound staged green-leaf optics and enzymes.
% Call after the final physical green LAI is staged, before optical RT.
% Accepted C.leaf/N.leaf [g element/m2 ground] are green leaves only. Their
% area density is formed with the actual exchange LAI; standing dead foliage
% and seed pools cannot contribute. This auxiliary projection changes no
% accepted state, pool, ledger or daily optimization memory.
% Source: LUNA N-to-capacity constants and nLeafToCab conversion; initial
% Cab is the explicit site optical prior bounded by available functional N.
% Editor: Wenjie Liu; Date: 2026-09-30.
    built=model.built;pa=built.p.plant_acclimation;
    if ~pa.luna_enabled||~strcmp(pa.luna_leaf_n_source,'plant_pool'),return;end
    pc=built.Y.plant_cn;
    area=built.workspace.scope_step.canopy.LAI;
    params=built.parameters.plant_cn.resolved;
    p=vegetation.integrated.LUNA.lunaParams(params.Jmaxb1);
    lnc=0;slatop=params.SLA0;
    if area>0&&pc.C.leaf>0&&pc.N.leaf>0
        lnc=pc.N.leaf/area;slatop=area/pc.C.leaf;
    end
    state=built.Y.plant_acclimation.luna_state;
    if state.first_day||isempty(state.Vcmax25_z)
        state.Vcmax25_z=p.Vcmax25_fallback;
        state.Jmax25_z=p.Jmax25_fallback;
        state.PNlc_z=0;state.PNresp_z=0;
        state=vegetation.integrated.LUNA.boundLeafNitrogen(state,lnc,slatop,p);
        functional=state.FNCa_z(1);
        if functional>0
            cab_prior=built.workspace.scope_step.leafbio.Cab;
            cab_for_all_n=vegetation.integrated.LUNA.nLeafToCab(functional);
            state.PNlc_z=min(1,max(0,cab_prior)/cab_for_all_n);
        end
    end
    state=vegetation.integrated.LUNA.boundLeafNitrogen(state,lnc,slatop,p);
    bio=built.workspace.scope_step.leafbio;
    bio.use_luna=1;
    bio.Vcmax25_LUNA=state.Vcmax25_z(1);
    bio.Jmax25_LUNA=state.Jmax25_z(1);
    bio.Vcmo=state.Vcmax25_z(1);
    built.workspace.scope_step.leafbio=bio;
    built.workspace.scope_step.accepted_leaf_n_cab= ...
        vegetation.integrated.LUNA.nLeafToCab(state.PNlc_z(1)*state.FNCa_z(1));
    built.workspace.scope_step.leaf_n_budget=struct( ...
        'lnc_g_m2_green_leaf',lnc,'green_lai',area, ...
        'available_functional_N_g_m2_green_leaf',state.FNCa_z(1), ...
        'used_functional_N_g_m2_green_leaf',state.N_functional_used_z(1), ...
        'capacity_scale',state.capacity_n_scale_z(1));
    model.built=built;
end
