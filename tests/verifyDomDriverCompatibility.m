function v=verifyDomDriverCompatibility()
% Verify legacy/marked adapter input, strict missing-data failures, rewetting
% with empty free pools, and rejection before an undefined dry donor is used.
    active=struct('sorption_model','mend_dynamic','diffusion_cm2_s',0,'dispersion_cm',0,'transport_substep_s',300);
    legacy=struct('solver_diagnostics',struct('QL',[0;0]),'water_flux_accum',struct('QMB',0));
    th=[.2;.3];lq=struct('QL',[0;0]);
    [a,b]=soil.dom.waterDriverSchedule(legacy,th,3600,lq,active);
    assert(strcmp(a.policy,'din_outer_mean_endpoint_v1')&&isempty(a.theta_initial)&&isequal(b.QL,[0;0])&&b.dt_s==3600);
    spec=struct('schema_version',1,'bottom_kind','advective_outflow');
    open=active; open.boundary=spec; legacy.water_flux_accum.QMB=-1e-6;
    legacy.solver_diagnostics.accepted_transport_schedule=struct('schema_version',1);
    [a,b]=soil.dom.waterDriverSchedule(legacy,th,3600,lq,open);
    assert(a.bottom_faces{1}.q_liquid_cm_s==-1e-6&&isequal(a.bottom_faces,b.bottom_faces));
    water=struct('schema_version',1,'tile_count',2,'elapsed_s',3600,'dt_s',[900,2700], ...
        'theta_liquid_initial',[.25;.25],'theta_liquid_endpoint',[.24 .2;.26 .3],'QL',[4e-6 0;0 4e-6]);
    lq.QL=water.QL*water.dt_s.'/3600;legacy.solver_diagnostics.dom_water_driver_version=2;
    legacy.solver_diagnostics.accepted_transport_schedule=water;
    [a,b]=soil.dom.waterDriverSchedule(legacy,th,3600,lq,active);
    assert(isequal(a.QL,water.QL)&&isequal(a.theta_initial,water.theta_liquid_initial) ...
        &&isequal(b.QL,lq.QL)&&isequal(b.theta,th)&&isscalar(b.dt_s)&&isequal(b.bottom_faces,{}));
    wrong=legacy;
    wrong.solver_diagnostics=rmfield(wrong.solver_diagnostics,'accepted_transport_schedule');
    fails(@() soil.dom.waterDriverSchedule(wrong,th,3600,lq,active),'soil:dom:waterDriverVersion');
    wrong=legacy;wrong.solver_diagnostics.accepted_transport_schedule=rmfield(water,'theta_liquid_initial');
    fails(@() soil.dom.waterDriverSchedule(wrong,th,3600,lq,active),'soil:dom:waterDriverSchema');
    wrong=legacy;wrong.solver_diagnostics=rmfield(wrong.solver_diagnostics,'dom_water_driver_version');
    fails(@() soil.dom.waterDriverSchedule(wrong,th,3600,lq,active),'soil:dom:waterDriverVersion');
    fails(@() soil.dom.waterDriverSchedule(legacy,th+1e-4,3600,lq,active),'soil:dom:waterDriverCoverage');
    fails(@() soil.dom.waterDriverSchedule(legacy,th,3600,struct('QL',[0;0]),active),'soil:dom:waterDriverMean');
    name={'doc','don'};closed=struct('kind','closed');
    % Liquid transfers 0.05 cm from one 1-cm layer to a previously dry
    % layer. A uniform incoming concentration must produce the same final C.
    old=[0;.2];last=[.05;.15];q=-.05/3600*ones(2,1);rho=1.3*ones(2,1);dz=ones(2,1);
    sorb=[.001;.002]; y=struct('schema_version',2,'sorption_model','mend_dynamic', ...
        'doc_sorbed_g_cm3',sorb,'don_sorbed_g_cm3',sorb/10);
    for i=1:2
        key=name{i};y.([key '_total_g_cm3'])=y.([key '_sorbed_g_cm3'])+old*1e-4/(10^(i-1));
    end
    wet=soil.dom.advanceAcceptedCandidate(y,last,rho,dz,q,active,3600,{},old);
    assert(max(abs(wet.state.doc_total_g_cm3-sorb-last*1e-4))<2e-18);
    assert(isequal(wet.state.doc_sorbed_g_cm3,sorb)&&abs(wet.transport.doc.closure_residual_g_cm2)<1e-17);
    [wetLabel,~]=soil.nitrogen.tracer.replayTransportLedger(.37*(y.don_total_g_cm3-y.don_sorbed_g_cm3),wet.transport.don,dz);
    assert(max(abs(wetLabel-.37*last*1e-5))<1e-18);
    le=active;le.sorption_model='linear_equilibrium';le.kd_doc_cm3_g=0;le.kd_don_cm3_g=0;
    ly=struct('doc_total_g_cm3',old*1e-4,'don_total_g_cm3',old*1e-5);
    eq=soil.dom.advanceAcceptedCandidate(ly,last,rho,dz,q,le,3600,{},old);
    assert(max(abs(eq.state.doc_total_g_cm3-last*1e-4))<1e-17);
    fails(@() soil.dom.advanceAcceptedCandidate(y,last,rho,dz,-q,active,3600,{},old),'soil:dom:dryDonor');
    yy=y;yy.doc_total_g_cm3(1)=sorb(1)+1e-6;
    fails(@() soil.dom.advanceAcceptedCandidate(yy,last,rho,dz,q,active,3600,{},old),'soil:dom:dynamicDryFreeMass');
    fails(@() soil.dom.advanceAcceptedCandidate(y,[0;.15],rho,dz,q,active,3600,{},old),'soil:dom:dynamicDryTransport');
    % Boundary case. Type 2 infiltration into an empty 1-cm dry cell carries
    % the specified external DOC concentration with exactly its liquid mass.
    bc=struct('kind','configured','q_liquid_cm_s',0,'q_top_cm_s',-.05/3600, ...
        'top',struct('type',2,'concentration_g_cm3_water',1e-4), ...
        'bottom',struct('type',3,'concentration_g_cm3_water',0));
    [s,ld]=soil.dom.transportConservative(0,.05,1.3,0,1,0,0,0,3600,300,[],bc,0);
    assert(abs(s-5e-6)<1e-18&&abs(ld.top_inflow_g_cm2-s)<1e-18);
    fails(@() soil.dom.transportConservative(0,0,1.3,0,1,0,0,0,3600,300,[],closed),'soil:dom:zeroCapacity');
    v=struct('status','PASS','old_schedule_and_missing_schedule_supported',true, ...
        'malformed_modern_rejected',true,'partial_wet_dynamic_and_linear_verified',true, ...
        'rewet_boundary_input_g_cm2',ld.top_inflow_g_cm2,'rewet_label_replay',true,'hours',0);
end

function fails(fn,identifier)
    try
        fn();
    catch err
        assert(strcmp(err.identifier,identifier),'domCheck:wrongError','Expected %s, got %s.',identifier,err.identifier);
        return
    end
    error('domCheck:unexpectedSuccess','Expected %s.',identifier);
end
