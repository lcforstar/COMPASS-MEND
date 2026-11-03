function [label,audit]=coProjectExclusiveTracers(label,total,weights,dz)
%COPROJECTEXCLUSIVETRACERS Share one frozen total-solute capacity across origins.
% label: NN x nSource individually bounded concentrations, total: NN x 1,
% both [ug N cm^-3 water], bottom-up. weights: exact frozen ADE end-storage
% weights [cm]; dz: bottom-up element widths [cm]. No total state is changed.
% At an overfull node move all origins in their local mixing fractions. Fill
% the closest nodes with unlabelled capacity first. Moving a source's mass
% (not its concentration) preserves its column inventory and external routes.
% Zero-storage nodes can be bounded but cannot store redistributed mass.
% Materially infeasible profiles reject the candidate. Feasible inputs are
% returned unchanged. The displacement is numerical, not a physical N flux.
% Editor: Wenjie Liu; Date: 2026-11-02.
    total=total(:);weights=weights(:);dz=dz(:);
    nn=numel(total);
    assert(size(label,1)==nn&&size(label,2)>=2&&numel(weights)==nn ...
        &&numel(dz)==nn-1&&all(isfinite([label(:);total;weights;dz])) ...
        &&all(weights>=0)&&all(dz>0)&&all(total>=0) ...
        &&all(label>=0,'all')&&all(label<=total,'all'), ...
        'soil:nitrogen:exclusiveTracerDomain','Require individually bounded finite source labels on the frozen ADE grid.');
    originalMass=weights.'*label;
    sumLabel=sum(label,2);
    donors=find(sumLabel>total);
    moved=zeros(size(originalMass));
    audit=struct('active',~isempty(donors),'initial_max_excess_ug_cm3',max(sumLabel-total), ...
        'redistributed_ug_cm2',moved,'inventory_residual_ug_cm2',moved, ...
        'from_node',zeros(0,1),'to_node',zeros(0,1),'distance_cm',zeros(0,1), ...
        'transfer_ug_cm2',zeros(0,size(label,2)));
    if isempty(donors),return;end

    capacityMass=weights.'*total;
    roundoff=32*nn*eps(max(capacityMass,sum(originalMass)));
    assert(sum(originalMass)<=capacityMass+roundoff, ...
        'soil:nitrogen:exclusiveTracerInfeasible','Combined source inventory exceeds the physical column capacity.');
    z=[0;cumsum(dz)];
    % Determine all initially admissible recipients before visiting a donor.
    room=weights.*max(total-sumLabel,0);
    for j=1:numel(donors)
        d=donors(j);
        excess=sumLabel(d)-total(d);
        fraction=label(d,:)./sumLabel(d);
        delta=excess.*fraction;
        label(d,:)=max(0,label(d,:)-delta);
        pending=weights(d)*excess;
        if pending==0,continue;end
        recipients=find(room>0);
        distances=abs(z(recipients)-z(d));
        for distance=reshape(unique(distances),1,[])
            % At equal distances share by the available unlabelled mass;
            % a shallower/deeper array ordering must not break the tie.
            nodes=recipients(distances==distance);
            available=sum(room(nodes));
            take=min(available,pending);
            amounts=take.*(room(nodes)./available);
            for k=1:numel(nodes)
                node=nodes(k);transfer=amounts(k).*fraction;
                label(node,:)=label(node,:)+transfer./weights(node);
                moved=moved+transfer;
                edge=numel(audit.from_node)+1;
                audit.from_node(edge,1)=d;
                audit.to_node(edge,1)=node;
                audit.distance_cm(edge,1)=distance;
                audit.transfer_ug_cm2(edge,:)=transfer;
            end
            room(nodes)=max(0,room(nodes)-amounts);
            pending=pending-take;
            if pending==0,break;end
        end
        assert(pending<=roundoff,'soil:nitrogen:exclusiveTracerInfeasible', ...
            'No positive-storage node has room for the remaining labelled inventory.');
    end
    residual=weights.'*label-originalMass;
    sourceRoundoff=32*nn*eps(originalMass);
    localRoundoff=16*eps(total);
    assert(all(label>=0,'all')&&all(sum(label,2)<=total+localRoundoff) ...
        &&all(abs(residual)<=sourceRoundoff), ...
        'soil:nitrogen:exclusiveTracerClosure', ...
        'Joint projection failed a nodal capacity or an individual source inventory.');
    audit.redistributed_ug_cm2=moved;
    audit.inventory_residual_ug_cm2=residual;
end
