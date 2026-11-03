function writeFertilizerTracerDetails(reg,built,kt,prefix,state,diagnostic)
%WRITEFERTILIZERTRACERDETAILS Export existing accepted source/phase ledgers.
% BUILT supplies accepted Y, workspace, NL-layer dz [cm] and bulk density
% [g dry soil cm-3 bulk]. Source replay routes are top-down [mg N g-1]
% amounts for this accepted step. State owners are bottom-up [g N cm-3].
% This function only writes diagnostics; physical state is never changed.
% Editor: Wenjie Liu; Date: 2026-09-22.
% Change: 2026-09-29 Wenjie Liu -- particle bottom faces match terminal exports.
% Change: 2026-09-29 Wenjie Liu -- also verify organic top exterior faces
% against accepted terminal losses; external organic inputs are unlabelled.
% Change: 2026-10-30 Wenjie Liu -- allow explicit passive-source cohorts;
% the original mineral-fertilizer caller and all its names remain unchanged.
    if nargin<4,prefix='fertilizer_tracer_';end
    if nargin<5,state=built.Y.soil_fertilizer_n_tracer;end
    if nargin<6,diagnostic=built.workspace.fertilizer_n_tracer;end
    schema=diagnostics.fertilizerTracerDetailSchema();
    replay=diagnostic.mend_routes;
    assert(strcmp(replay.grid_order,'top_down') ...
        &&strcmp(replay.units,'mg_N_g_dry_soil'), ...
        'diagnostics:tracer:routeContract','Unsupported accepted source-route units/order.');
    n=built.domain.NL;dz=flip(built.domain.dz(:));
    conversion=flip(built.parameters.soil.bulk_density(:))./1000;
    for i=1:size(schema.groups,1)
        for j=1:numel(schema.phases)
            phase=schema.phases{j};field=[schema.groups{i,2} '_' phase '_N'];
            value=flip(layer(state.organic_n.(field),n));
            write([schema.groups{i,1} '_' phase '_n_profile'],value);
        end
    end
    for i=1:size(schema.routes,1)
        value=layer(replay.nitrogen_flux.(schema.routes{i,2}),n).*conversion;
        key=['route_' schema.routes{i,1}];
        write([key '_profile_step'],value);
        write([key '_step'],sum(value.*dz));
    end
    for i=1:size(schema.phase_routes,1)
        matrix=replay.(schema.phase_routes{i,2}).(schema.phase_routes{i,3});
        assert(isequal(size(matrix),[n,3])&&all(isfinite(matrix),'all'), ...
            'diagnostics:tracer:phaseRoute','Accepted phase route must be finite NL-by-3.');
        for j=1:3
            write(['route_' schema.phase_routes{i,1} '_' schema.phases{j} ...
                '_profile_step'],matrix(:,j).*conversion);
        end
    end
    theta=mean(built.workspace.derived.Theta_LL,2);
    phase=soil.dom.phaseState(state.dom,'don',theta, ...
        built.parameters.soil.bulk_density(:),built.parameters.soil_dom);
    write('don_free_profile',flip(layer(phase.free_g_cm3_bulk,n)));
    write('don_sorbed_profile',flip(layer(phase.sorbed_g_cm3_bulk,n)));
    write('don_aqueous_profile',flip(layer(phase.aqueous_g_cm3_water,n)));
    pools={'pon1','PON1';'pon2','PON2';'mon','MON'; ...
        'active_microbe_n','MB_N';'dormant_microbe_n','BD_N'; ...
        'pom1_enzyme_n','ENZP1_N';'pom2_enzyme_n','ENZP2_N'; ...
        'mom_enzyme_n','EM_N';'intracellular_enzyme_n','ENZN_N'};
    write('don_storage',sum(flip(layer(state.dom.don_total_g_cm3,n)).*dz));
    for i=1:size(pools,1)
        write([pools{i,1} '_storage'], ...
            sum(flip(layer(state.organic_n.(pools{i,2}),n)).*dz));
    end
    writeFaces('don',diagnostic.don_transport);
    transportEnabled=isfield(built.control.switches,'lMicrobialTransport') ...
        &&built.control.switches.lMicrobialTransport==1;
    for i=1:size(schema.groups,1)
        key=schema.groups{i,1};group=schema.groups{i,2};
        if transportEnabled
            transport=diagnostic.bioparticle_transport.(group);
        else
            % A disabled process has a physical zero transfer, not missing data.
            transport=struct('upward_face_amount_g_cm2',zeros(n+1,1), ...
                'downward_face_amount_g_cm2',zeros(n+1,1), ...
                'net_upward_face_amount_g_cm2',zeros(n+1,1), ...
                'attachment_g_cm3',zeros(n,1),'detachment_g_cm3',zeros(n,1), ...
                'straining_g_cm3',zeros(n,1),'bottom_outflow_g_cm2',0,'top_outflow_g_cm2',0);
        end
        writeFaces(key,transport);
        for process={'attachment','detachment','straining'}
            value=flip(layer(transport.([process{1} '_g_cm3']),n));
            name=['transport_' key '_' process{1}];
            write([name '_profile_step'],value);
            write([name '_step'],sum(value.*dz));
        end
    end
    function write(key,value)
        reg.write([prefix key],value,kt);
    end
    function writeFaces(key,transport)
        for direction={'upward','downward','net_upward'}
            value=flip(layer(transport.([direction{1} '_face_amount_g_cm2']),n+1));
            % External organic input is unlabelled; labelled export is terminal.
            expected=0;expectedTop=0;
            if ~strcmp(direction{1},'downward'),expectedTop=transport.top_outflow_g_cm2;end
            if strcmp(direction{1},'downward'),expected=transport.bottom_outflow_g_cm2;end
            if strcmp(direction{1},'net_upward'),expected=-transport.bottom_outflow_g_cm2;end
            assert(value(1)==expectedTop&&value(end)==expected, ...
                'diagnostics:tracer:boundaryFace','Source-label exterior face disagrees with export ledger.');
            write(['transport_' key '_' direction{1} '_face_step'],value);
        end
    end
end

function value=layer(value,n)
    assert(isnumeric(value)&&isreal(value)&&isvector(value) ...
        &&numel(value)==n&&all(isfinite(value)), ...
        'diagnostics:tracer:layer','Accepted tracer route/pool must be finite NL-vector.');
    value=value(:);
end
