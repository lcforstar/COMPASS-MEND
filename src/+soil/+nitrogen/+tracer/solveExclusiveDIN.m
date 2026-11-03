function [dins,ledgers]=solveExclusiveDIN(frozen,requests)
%SOLVEEXCLUSIVEDIN Replay one or more exclusive N origins on a frozen ADE.
% REQUESTS is an ordered cell of old nodal concentrations [ug N cm^-3 water]
% and explicit gross rates [ug N cm^-3 soil s^-1] on the identical operator.
% Output cells retain that source order. The raw ADE owns boundary exchange;
% final DIN storage and denitrification fractions use the admissible endpoint.
% Original single-source fixed points are kept when jointly admissible.
% Editor: Wenjie Liu; Date: 2026-11-02.
% Preserve the independent solver and its original stopping policy whenever
% its endpoints already satisfy the shared admissibility contract.
    dins=cell(size(requests));ledgers=cell(size(requests));
    for i=1:numel(requests)
        [dins{i},ledgers{i},requests{i}]=iterateDenitrification(frozen,requests{i});
    end
    if numel(requests)==1,return;end
    for iteration=1:20
        if iteration>1
            for i=1:numel(requests)
                [dins{i},ledgers{i}]=soil.nitrogen.transport.solveFrozenFertilizerTracer(frozen,requests{i});
            end
        end
        anyProjection=false;
        for species={'nh4','no3'}
            name=species{1};field=[name '_ug_cm3'];
            values=cellfun(@(c)c.(field),dins,'UniformOutput',false);
            values=cat(2,values{:});
            map=frozen.(name).projection_map;
            [joint,audit]=soil.nitrogen.transport.coProjectExclusiveTracers( ...
                values,map.projected_total_ug_cm3,map.storage_weights_cm,frozen.dz_cm);
            if ~audit.active,continue;end
            anyProjection=true;
            for i=1:numel(requests)
                dins{i}.(field)=joint(:,i);
                % Source/sink integrals and boundary amounts are properties
                % of the raw frozen ADE solve. Only its final storage uses
                % the jointly admissible endpoint. No new flux is booked.
                item=ledgers{i}.(name);
                item.final_storage_ug_cm2=map.storage_weights_cm(:).'*joint(:,i);
                item.storage_change_ug_cm2=item.final_storage_ug_cm2-item.initial_storage_ug_cm2;
                item.residual_ug_cm2=item.storage_change_ug_cm2-item.expected_change_ug_cm2;
                item.exclusive_projection=struct('initial_max_excess_ug_cm3',audit.initial_max_excess_ug_cm3, ...
                    'redistributed_ug_cm2',audit.redistributed_ug_cm2(i), ...
                    'inventory_residual_ug_cm2',audit.inventory_residual_ug_cm2(i), ...
                    'from_node',audit.from_node,'to_node',audit.to_node, ...
                    'distance_cm',audit.distance_cm,'transfer_ug_cm2',audit.transfer_ug_cm2(:,i));
                ledgers{i}.(name)=item;
            end
        end
        if iteration==1&&~anyProjection,return;end
        converged=true;
        for i=1:numel(requests)
            fraction=soil.nitrogen.tracer.nodeFraction(dins{i}.no3_ug_cm3,frozen.no3.projected_total_ug_cm3,'joint NO3');
            target=frozen.gross_rates.no3.sinks.denitrification.* ...
                [fraction(1:frozen.NL),fraction(2:frozen.NN)];
            used=requests{i}.rates.no3.sinks.denitrification;
            change=max(abs(target-used),[],'all');
            converged=converged&&(change<=1e-11*max(1,max(abs(target),[],'all')));
            requests{i}.rates.no3.sinks.denitrification=0.5*used+0.5*target;
            ledger=ledgers{i};
            ledger.combined_din.storage_change_ug_cm2=ledger.nh4.storage_change_ug_cm2+ledger.no3.storage_change_ug_cm2;
            ledger.combined_din.residual_ug_cm2=ledger.nh4.residual_ug_cm2+ledger.no3.residual_ug_cm2;
            ledger.denitrification_label_iterations=iteration;
            ledger.denitrification_label_relaxation=0.5;
            ledger.denitrification_label_residual=change;
            ledger.exclusive_source_iterations=iteration;
            ledgers{i}=ledger;
        end
        if converged,return;end
    end
    error('soil:nitrogen:tracer:denitrificationIteration', ...
        'Joint passive-source denitrification attribution did not converge.');
end

function [candidate,ledger,request]=iterateDenitrification(frozen,request)
    for iteration=1:20
        [next,nextLedger]=soil.nitrogen.transport.solveFrozenFertilizerTracer( ...
            frozen,request);
        fraction=soil.nitrogen.tracer.nodeFraction(next.no3_ug_cm3, ...
            frozen.no3.projected_total_ug_cm3,'accepted NO3');
        updated=frozen.gross_rates.no3.sinks.denitrification.* ...
            [fraction(1:frozen.NL),fraction(2:frozen.NN)];
        change=max(abs(updated-request.rates.no3.sinks.denitrification),[],'all');
        scale=max(1,max(abs(updated),[],'all'));
        candidate=next;ledger=nextLedger;
        if change<=1e-11*scale
            ledger.denitrification_label_iterations=iteration;
            ledger.denitrification_label_relaxation=0.5;
            ledger.denitrification_label_residual=change;
            return
        end
        % Endpoint attribution is a negative-feedback fixed point. A convex
        % half update damps alternating iterates without changing the fixed
        % point or the admissible [0,total gross rate] interval. Convergence
        % above uses the full, undamped residual, never the smaller update.
        request.rates.no3.sinks.denitrification= ...
            0.5*request.rates.no3.sinks.denitrification+0.5*updated;
    end
    error('soil:nitrogen:tracer:denitrificationIteration', ...
        'Passive fertilizer-N denitrification attribution did not converge.');
end
