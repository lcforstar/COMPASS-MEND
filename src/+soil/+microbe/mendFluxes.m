function f = mendFluxes(state,p)
%MENDFLUXES Evaluate selected MEND carbon flux rates [mg C g-1 h-1].
%
% This is the carbon-only flux topology of MOD_MEND.F90 commit 92323c7,
% iKinetics=0/1/2. Dynamic QOM sorption is deliberately absent because COMPASS
% total DOM is the sole dissolved/sorbed owner.
%
% Editor: Jizixing Zhu; Date: 2026-09-12.
% Change: 2026-09-12 Codex -- expose source MM, first-order, and
%   second-order POM/MOM kinetic branches; DOM uptake remains source MM.
    required={'POM1','POM2','MOC','MB','BD','DOC','ENZP1','ENZP2','EM'};
    for i=1:numel(required)
        v=state.(required{i});
        assert(iscolumn(v)&&all(isfinite(v))&&all(v>=-1e-12), ...
            'soil:microbe:mendState','Invalid MEND state %s.',required{i});
    end
    f.POM1_decomp=kinetic(p.kinetics,p.VdPOM_h(:,1),state.POM1,state.ENZP1, ...
        p.KsPOM_mg_g(:,1));
    f.POM2_decomp=kinetic(p.kinetics,p.VdPOM_h(:,2),state.POM2,state.ENZP2, ...
        p.KsPOM_mg_g(:,2));
    f.POM1_to_DOM=p.frPOM2DOM.*f.POM1_decomp;
    f.POM2_to_DOM=p.frPOM2DOM.*f.POM2_decomp;
    f.POM1_to_MOM=(1-p.frPOM2DOM).*f.POM1_decomp;
    f.POM2_to_MOM=(1-p.frPOM2DOM).*f.POM2_decomp;
    f.MOM_to_DOM=kinetic(p.kinetics,p.VdMOM_h,state.MOC,state.EM,p.KsMOM_mg_g);
    f.DOM_to_MBA=mm((p.Vg_h+p.Vm_h)./p.Yg,state.DOC,state.MB, ...
        p.KsDOM_mg_g);
    f.CO2_growth=mm(p.Vg_h.*(1./p.Yg-1),state.DOC,state.MB, ...
        p.KsDOM_mg_g);
    f.CO2_maintenance_MBA=mm(p.Vm_h.*(1./p.Yg-1),state.DOC, ...
        state.MB,p.KsDOM_mg_g);
    f.CO2_maintenance_MBD=p.VmD_h.*state.BD;
    f.MBA_mortality=p.rMORT_h.*state.MB;
    f.MBA_to_DOM=p.frMB2DOM.*f.MBA_mortality;
    particulate=(1-p.frMB2DOM).*f.MBA_mortality;
    f.MBA_to_POM1=0.1.*particulate;
    f.MBA_to_POM2=0.9.*particulate;
    totalPOM=state.POM1+state.POM2;
    share1=zeros(size(totalPOM)); active=totalPOM>0;
    share1(active)=state.POM1(active)./totalPOM(active);
    % Source divides by sumPOM. The zero-POM extension sends production to
    % hydrolytic enzyme, avoiding source NaN while preserving total production.
    f.MBA_to_ENZP1=share1.*p.pENZP.*p.Vm_h.*state.MB;
    f.MBA_to_ENZP2=(1-share1).*p.pENZP.*p.Vm_h.*state.MB;
    f.MBA_to_ENZM=p.pENZM.*p.Vm_h.*state.MB;
    f.ENZP1_to_DOM=p.rENZP_h(1).*state.ENZP1;
    f.ENZP2_to_DOM=p.rENZP_h(2).*state.ENZP2;
    f.ENZM_to_DOM=p.rENZM_h.*state.EM;
    phi=zeros(size(state.DOC)); denom=state.DOC+p.KsDOM_mg_g;
    positive=denom>0; phi(positive)=state.DOC(positive)./denom(positive);
    f.substrate_saturation=phi;
    f.MBA_to_MBD=(1-phi).*p.VmA2D_h.*state.MB;
    f.MBD_to_MBA=phi.*p.VmD2A_h.*state.BD;
end

function value=kinetic(kind,vm,substrate,enzyme,km)
    switch string(kind)
        case "MM"
            value=mm(vm,substrate,enzyme,km);
        case "first_order"
            value=max(vm.*substrate,0);
        case "second_order"
            value=max(vm.*enzyme.*substrate,0);
        otherwise
            error('soil:microbe:mendKinetics', ...
                'Unsupported MEND kinetics branch: %s.',string(kind));
    end
end

function value=mm(vm,substrate,enzyme,km)
% Source fV_MM, interpreted as an hourly derivative without its dt=1 cap.
    denom=km+substrate;
    value=zeros(size(substrate)); active=denom>0;
    value(active)=vm(active).*substrate(active).*enzyme(active)./denom(active);
    value=max(value,0);
end
