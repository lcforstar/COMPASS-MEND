function report=verifyExclusiveDinTracer(fixtureDir)
%VERIFYEXCLUSIVEDINTRACER Check exclusive nodal bounds, origin mass and locality.
% Optional fixture: three producer-captured original frozen ADE requests from
% Changwu cand0030 KT140575. No physical model step is executed here.
% Editor: Wenjie Liu; Date: 2026-11-02.
    total=[2;1;1;2];w=[2;1;3;1];dz=[1;1;3];
    input=[0.1,0,0;0.8,0.4,0.3;0.1,0,0;0.5,0.3,0];
    [out,a]=soil.nitrogen.transport.coProjectExclusiveTracers(input,total,w,dz);
    assert(a.active&&all(sum(out,2)<=total+4*eps(total))&&all(out>=0,'all'));
    assert(norm(w.'*(out-input),inf)<1e-14&&all(a.distance_cm==1));
    assert(isequal(out(4,:),input(4,:))&&isequal(sort(a.to_node),[1;3]));
    [permuted,~]=soil.nitrogen.transport.coProjectExclusiveTracers(input(:,[3,1,2]),total,w,dz);
    assert(max(abs(permuted-out(:,[3,1,2])),[],'all')<=16*eps);
    [already,b]=soil.nitrogen.transport.coProjectExclusiveTracers(input/2,total,w,dz);
    assert(~b.active&&isequal(already,input/2));
    zeroW=[2;0;3;1];zero=soil.nitrogen.transport.coProjectExclusiveTracers(input,total,zeroW,dz);
    assert(isequal(zero([1,3,4],:),input([1,3,4],:))&&isequal(zeroW.'*zero,zeroW.'*input));
    tiny=[0,0;8e-20,8e-20;0,0];
    small=soil.nitrogen.transport.coProjectExclusiveTracers(tiny,[1e-19;1e-19;1e-19],ones(3,1),ones(2,1));
    assert(all(sum(small)>0)&&max(abs(sum(small)-sum(tiny)))<1e-34);
    rejected=false;
    try
        soil.nitrogen.transport.coProjectExclusiveTracers(.7*ones(4,2),ones(4,1),w,dz);
    catch err
        rejected=strcmp(err.identifier,'soil:nitrogen:exclusiveTracerInfeasible');
    end
    assert(rejected,'A materially infeasible total labelled inventory must reject.');
    report=struct('status','PASS','manufactured',true,'captured_failure_replayed',false, ...
        'coupled_denitrification_iterations',requireDenitrificationFeedback());
    if nargin==0||isempty(fixtureDir),return;end
    names={'mineral','manure_organic','manure_tan'};requests=cell(1,3);old=cell(1,3);
    for i=1:3
        record=load(fullfile(fixtureDir,[names{i} '.mat']),'frozen','request','din','adeLedger');
        if i==1,frozen=record.frozen;else,assert(isequaln(frozen,record.frozen));end
        requests{i}=record.request;old{i}=record.din;
        [single,~]=soil.nitrogen.tracer.solveExclusiveDIN(frozen,{record.request});
        used=record.adeLedger.used_request;
        original=soil.nitrogen.transport.solveFrozenFertilizerTracer(frozen,used);
        fromUsed=soil.nitrogen.tracer.solveExclusiveDIN(frozen,{used});
        for key={'nh4_ug_cm3','no3_ug_cm3'}
            field=key{1};remote=record.din.(field);
            assert(isequal(fromUsed{1}.(field),original.(field)), ...
                'A single origin must exactly retain the existing native ADE result for the identical used request.');
            assert(all(abs(single{1}.(field)-remote)<=64*eps(max(1,max(abs(remote))))), ...
                'Independent replay of the server-captured source differs beyond floating-point roundoff.');
        end
    end
    getNo3=@(d)d.no3_ug_cm3;
    before=cell2mat(cellfun(getNo3,old,'UniformOutput',false));
    total=frozen.no3.projected_total_ug_cm3;
    assert(sum(before(9,:))-total(9)>4e-7,'Captured case must reproduce the actual old failure.');
    [dins,ledgers]=soil.nitrogen.tracer.solveExclusiveDIN(frozen,requests);
    after=cell2mat(cellfun(getNo3,dins,'UniformOutput',false));
    assert(all(sum(after,2)<=total+16*eps(total))&&all(after>=0,'all'));
    weights=frozen.no3.projection_map.storage_weights_cm;
    for i=1:3
        a=ledgers{i}.no3.exclusive_projection;
        assert(all(a.from_node==9)&&all(ismember(a.to_node,[8,10]))&&all(a.distance_cm==25));
        assert(ledgers{i}.denitrification_label_residual<=1e-11*max(1,max(frozen.gross_rates.no3.sinks.denitrification,[],'all')));
        assert(abs(weights.'*(after(:,i)-dins{i}.raw_no3_ug_cm3))<1e-10);
        assert(abs(ledgers{i}.no3.residual_ug_cm2)<1e-8&&abs(ledgers{i}.combined_din.residual_ug_cm2)<1e-8);
    end
    report.captured_failure_replayed=true;
    report.source_order=names;
    report.initial_excess_ug_cm3=sum(before(9,:))-total(9);
    report.final_excess_ug_cm3=sum(after(9,:))-total(9);
    report.corrected_source_concentrations_ug_cm3=after(9,:);
    report.shifted_source_mass_ug_cm2=cellfun(@(l)l.no3.exclusive_projection.redistributed_ug_cm2,ledgers);
    report.mass_residual_ug_cm2=cellfun(@(l)l.no3.exclusive_projection.inventory_residual_ug_cm2,ledgers);
    report.denit_rate_residual=cellfun(@(l)l.denitrification_label_residual,ledgers);
    report.joint_iterations=cellfun(@(l)l.exclusive_source_iterations,ledgers);
end

function count=requireDenitrificationFeedback()
% Three closed nodes, lumped end storage, no transport. The physical DIN
% endpoint includes a known FE gross immobilization load and denitrification.
% A high-labelled middle node overfills after subtracting the unlabeled
% immobilization. Its coupled endpoint fraction must alter the actual gas
% attribution. An endpoint-only repair without feedback fails this test.
    z=zeros(2,2);dz=[1;1];dt=10;w=[.5;1;.5];old=ones(3,1);
    immob=[0,.023;.023,0];denit=1e-6*ones(2,2);sink=immob+denit;
    feLoad=[(2*sink(1,1)+sink(1,2))/6; ...
        (sink(1,1)+2*sink(1,2)+2*sink(2,1)+sink(2,2))/6; ...
        (sink(2,1)+2*sink(2,2))/6];
    final=old-dt*feLoad./w;
    frozen=struct('schema_version',1,'dt_s',dt,'NL',2,'NN',3,'dz_cm',dz, ...
        'ql_nodes_cm_s',zeros(3,1),'boundary',struct('NBCN',2,'NBCNB',2,'lThmrlefc',1,'QMB_now',0));
    gross=struct('nh4',struct('sources',struct('organic_release',z),'sinks',struct('nitrification',z)), ...
        'no3',struct('sources',struct('nitrification',z),'sinks',struct('organic_immobilization',immob,'denitrification',denit)));
    frozen.gross_rates=gross;
    matrix=[w./dt,zeros(3,1)];
    for species={'nh4','no3'}
        s=species{1};if strcmp(s,'nh4'),total=old;else,total=final;end
        map=struct('raw_total_ug_cm3',total,'projected_total_ug_cm3',total, ...
            'storage_weights_cm',w,'projected_total_mass_ug_cm2',w.'*total, ...
            'mass_tolerance_ug_cm2',256*eps(w.'*total));
        frozen.(s)=struct('old_total_ug_cm3',old,'projected_total_ug_cm3',total, ...
            'implicit_reaction_coefficient_s_1',z,'old_storage_lumped',w,'matrix_lower',zeros(2,1), ...
            'matrix_diag_upper',matrix,'pre_boundary_matrix_diag_upper',matrix,'projection_map',map);
    end
    requests=cell(1,2);concentration=[.06,.14;.50,.45;.12,.08];
    for i=1:2
        f=concentration(:,i);rates=gross;
        rates.no3.sinks.organic_immobilization=z;
        rates.no3.sinks.denitrification=denit.*[f(1:2),f(2:3)];
        requests{i}=struct('old',struct('nh4_ug_cm3',zeros(3,1),'no3_ug_cm3',f),'rates',rates);
    end
    [dins,ledgers]=soil.nitrogen.tracer.solveExclusiveDIN(frozen,requests);
    summed=dins{1}.no3_ug_cm3+dins{2}.no3_ug_cm3;
    assert(all(summed<=final+16*eps(final)));
    for i=1:2
        fractions=dins{i}.no3_ug_cm3./final;
        sinkAmount=dt*sum(mean(denit.*[fractions(1:2),fractions(2:3)],2));
        assert(abs(sinkAmount-ledgers{i}.no3.sinks_ug_cm2.denitrification)<=sum(dz)*dt*1e-11);
        assert(abs(ledgers{i}.no3.residual_ug_cm2)<1e-14);
        assert(ledgers{i}.denitrification_label_residual<=1e-11);
    end
    count=ledgers{1}.exclusive_source_iterations;
    assert(count>1,'This fixture must require correction of source-specific denitrification rates.');
end
