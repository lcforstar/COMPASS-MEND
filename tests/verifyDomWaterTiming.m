function receipt=verifyDomWaterTiming()
% Verify uniform concentration under locally conservative liquid reallocation.
% All columns and fluxes follow the native bottom-to-surface convention.
% Editor: Wenjie Liu; Date: 2026-11-03.
    dz=5*ones(6,1);rho=1.3*ones(6,1);start=[.20;.23;.24;.27;.29;.28];
    waterFlux=-[0;.10;.30;.55;.10;0]/3600;duration=[1800 1800];
    q=[waterFlux,-0.6*waterFlux];
    liquidFaces=[zeros(1,2);0.5*(q(1:end-1,:)+q(2:end,:));zeros(1,2)];
    delta=-diff(liquidFaces,1,1).*duration./dz;
    ends=start+cumsum(delta,2); c0=100e-6; n0=10e-6;
    p=struct('sorption_model','mend_dynamic','diffusion_cm2_s',0, ...
        'dispersion_cm',0,'transport_substep_s',300);
    sorbedC=(1:6)'*1e-5;sorbedN=sorbedC/15;
    initial=struct('schema_version',2,'sorption_model',p.sorption_model, ...
        'doc_total_g_cm3',start*c0+sorbedC,'don_total_g_cm3',start*n0+sorbedN, ...
        'doc_sorbed_g_cm3',sorbedC,'don_sorbed_g_cm3',sorbedN);
    fixed=soil.dom.advanceAcceptedCandidate(initial,ends,rho,dz,q,p,duration,{});
    sol=soil.dom.advanceAcceptedCandidate(initial,ends,rho,dz,q,p,duration,{},start);
    C=(sol.state.doc_total_g_cm3-sorbedC)./ends(:,end);
    N=(sol.state.don_total_g_cm3-sorbedN)./ends(:,end);
    fixedC=(fixed.state.doc_total_g_cm3-sorbedC)./ends(:,end);
    rel=max(abs([C/c0-1;N/n0-1]));
    assert(rel<2e-13 && max(abs(fixedC/c0-1))>1e-3);
    assert(isequal(initial.doc_sorbed_g_cm3,sol.state.doc_sorbed_g_cm3) ...
        &&isequal(initial.don_sorbed_g_cm3,sol.state.don_sorbed_g_cm3));
    assert(sol.transport.doc.driver_segments==2 && sol.transport.don.accepted_substeps==12);
    [tag,~]=soil.nitrogen.tracer.replayTransportLedger(.37*start*n0,sol.transport.don,dz);
    assert(max(abs(tag-.37*ends(:,end)*n0))<1e-18);
    % Under instantaneous linear equilibrium, uniform aqueous C also
    % remains exact with constant (depth-varying) sorbed inventories.
    ep=p;ep.sorption_model='linear_equilibrium';ep.kd_doc_cm3_g=(1:6)'/10;ep.kd_don_cm3_g=ep.kd_doc_cm3_g;
    eq=struct('doc_total_g_cm3',(start+rho.*ep.kd_doc_cm3_g)*c0, ...
        'don_total_g_cm3',(start+rho.*ep.kd_don_cm3_g)*n0);
    e=soil.dom.advanceAcceptedCandidate(eq,ends,rho,dz,q,ep,duration,{},start);
    eqError=max(abs(e.state.doc_total_g_cm3./(ends(:,end)+rho.*ep.kd_doc_cm3_g)/c0-1));
    assert(eqError<2e-13);
    % Boundary concentration forces a diffusive exchange even at zero q.
    % Two explicit steps independently evaluate that conductance at t=0 and half-time.
    b=struct('kind','configured','q_liquid_cm_s',0,'q_top_cm_s',0, ...
        'bottom',struct('type',1,'concentration_g_cm3_water',2*c0), ...
        'top',struct('type',1,'concentration_g_cm3_water',0));
    old=start*c0;last=ends(:,1);mid=.5*(start+last);q0=zeros(6,1);
    [a,ledger]=soil.dom.transportConservative(old,last,rho,0,dz,q0,5e-5,5,600,300,[],b,start);
    [half,h]=soil.dom.transportConservative(old,start,rho,0,dz,q0,5e-5,5,300,300,[],b);
    [ref,k]=soil.dom.transportConservative(half,mid,rho,0,dz,q0,5e-5,5,300,300,[],b);
    assert(max(abs(a-ref))<1e-17);
    assert(abs(ledger.bottom_inflow_g_cm2-h.bottom_inflow_g_cm2-k.bottom_inflow_g_cm2)<1e-17);
    assert(abs(ledger.top_outflow_g_cm2-h.top_outflow_g_cm2-k.top_outflow_g_cm2)<1e-17);
    % Type-2 top evaporates water while retaining dissolved carbon. The
    % porewater concentration must increase as water volume falls.
    b.top.type=2;b.bottom.type=3;b.bottom.concentration_g_cm3_water=0;
    b.q_top_cm_s=.05/3600;thetaOld=.30;thetaDry=thetaOld-.05/5;
    [remaining,ev]=soil.dom.transportConservative(thetaOld*c0,thetaDry, ...
        1.3,0,5,0,0,0,3600,300,[],b,thetaOld);
    assert(remaining==thetaOld*c0&&remaining/thetaDry>c0 ...
        &&ev.top_outflow_g_cm2==0&&ev.closure_residual_g_cm2==0);
    % At the 5:1 interface the face lies five-sixths of the distance from
    % the lower to the upper centre. The simple mean would reverse flow.
    sizes=[5;1]; centres=[-.8;.4]/3600; exact=.2/3600;
    unequalError=0;
    for direction=[-1 1]
        s=[.3;.25]; qFace=direction*exact;
        finish=s+3600*qFace*[-1/sizes(1);1/sizes(2)];
        [mass,led]=soil.dom.transportConservative(s*c0,finish,[1;1],0,sizes, ...
            direction*centres,0,0,3600,300,[],struct('kind','closed'),s);
        unequalError=max(unequalError,max(abs(mass./finish/c0-1)));
        assert(abs(led.net_layer_change_g_cm3(2)*sizes(2)-3600*qFace*c0)<1e-15);
    end
    s=[0;.3];finish=[.04;.1];
    mass=soil.dom.transportConservative(s*c0,finish,[1;1],0,sizes,-centres, ...
        0,0,3600,300,[],struct('kind','closed'),s);
    assert(unequalError<2e-13&&max(abs(mass./finish/c0-1))<2e-13);
    receipt=struct('status','PASS','dynamic_uniform_relative_error',rel, ...
        'unequal_grid_uniform_relative_error',unequalError, ...
        'linear_uniform_relative_error',eqError,'control_fixed_theta_error',max(abs(fixedC/c0-1)), ...
        'boundary_reference_max_difference_g_cm3',max(abs(a-ref)), ...
        'closed_carbon_residual_g_cm2',sol.transport.doc.closure_residual_g_cm2, ...
        'evaporation_retained_C_g_cm3_bulk',remaining);
end
