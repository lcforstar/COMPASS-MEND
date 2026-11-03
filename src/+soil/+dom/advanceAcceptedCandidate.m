function candidate = advanceAcceptedCandidate( ...
        accepted, theta, bulkDensity, dz_cm, q_cm_s, parameters, dt_s, bottomFaces,thetaInitial)
%ADVANCEACCEPTEDCANDIDATE Conservatively transport DOC and DON.
%
%   This DOM owner transports bottom-up element total DOC/DON [g/cm3 bulk
%   soil]. Linear-equilibrium states repartition by Kd; explicitly tagged
%   dynamic states retain their sorbed inventory and transport only free mass.
%   DON reactions belong to the paired microbial C/N ledger called from
%   soil.nitrogen.advanceAcceptedCandidate. This routine has no independent
%   DON-to-NH4 or DON-to-NO3 reaction and never mutates the accepted input.
%   theta, q_cm_s and dt_s may also carry a chronological accepted-tile
%   schedule (NL-by-K, NL-by-K and 1-by-K). Each QL column is then replayed
%   with the liquid-water endpoint from the same accepted soil tile.
%   With thetaInitial [NL x 1, liquid cm3/cm3] the explicit solute steps use
%   interpolated liquid storage from the start to the end of each water tile.
%   Omitting it retains the prescribed piecewise-constant theta interface.
%
%   Editor: Jizixing Zhu; Date: 2026-08-21.
% Change: 2026-09-01 Codex -- remove the standalone DON ammonification
%   bypass so all soil DON processing is owned by the explicit microbial C-N
%   ledger; MEMC/MEND source equations are unchanged.
% Change: 2026-09-14 Jizixing Zhu -- support explicitly supplied dynamic
%   phase inventories; transport free DOM only and retain sorbed C/N exactly.
% Change: 2026-09-21 Wenjie Liu -- replay every accepted water tile with its
%   matching endpoint theta instead of applying the final tile over the outer
%   duration; aggregate transport diagnostics over the complete outer step.
% Change: 2026-09-29 Wenjie Liu -- accepted-tile bottom aqueous outflow;
% cumulative exports are detached candidate memory, installed only by commit.
% Change: 2026-09-29 Wenjie Liu -- separate DOC/DON numbered face policies;
% explicit external concentrations and all directed boundary cumulative routes.
% Editor: Wenjie Liu; Date: 2026-10-17.
% Change: 2026-10-17 Wenjie Liu -- production can use one DIN-aligned outer
%   driver segment; reject every negative DON inventory without a final clamp.
% Change: 2026-11-03 Wenjie Liu -- pass both liquid storage endpoints to
%   the free DOC/DON operator; dynamic sorbed matter remains immobile.
    if nargin<8,bottomFaces={};end
    if nargin<9,thetaInitial=[];end
    bulkDensity=bulkDensity(:);dz_cm=dz_cm(:);
    [thetaSchedule,qSchedule,dtSchedule]=normalizeSchedule( ...
        theta,q_cm_s,dt_s,numel(bulkDensity));
    openBottom=isfield(parameters,'boundary') ...
        &&~strcmp(parameters.boundary.bottom_kind,'closed');
    if openBottom
        assert(iscell(bottomFaces)&&numel(bottomFaces)==numel(dtSchedule), ...
            'soil:dom:missingBottomSchedule','Open DOM requires every accepted bottom-face tile.');
    else
        assert(isempty(bottomFaces),'soil:dom:closedBottomDriver', ...
            'Do not pass open bottom drivers to a closed DOM boundary.');
        bottomFaces=repmat({struct('kind','closed')},size(dtSchedule));
    end
    next=accepted;docTransport=[];donTransport=[];
    for segment=1:numel(dtSchedule)
        thetaEnd=thetaSchedule(:,segment);
        if isempty(thetaInitial)
            thetaStart=thetaEnd;
        elseif segment==1
            thetaStart=thetaInitial(:);
        else
            thetaStart=thetaSchedule(:,segment-1);
        end
        [next,docPart,donPart]=advanceSegment(next,thetaStart,thetaEnd, ...
            bulkDensity,dz_cm,qSchedule(:,segment),parameters,dtSchedule(segment),bottomFaces{segment});
        docTransport=accumulateTransport(docTransport,docPart);
        donTransport=accumulateTransport(donTransport,donPart);
    end
    docTransport=finalizeTransport(docTransport, ...
        next.doc_total_g_cm3-accepted.doc_total_g_cm3,dtSchedule);
    donTransport=finalizeTransport(donTransport, ...
        next.don_total_g_cm3-accepted.don_total_g_cm3,dtSchedule);
    zeroRate=zeros(size(next.don_total_g_cm3));
    if openBottom
        assert(isfield(accepted,'bottom_export_cumulative'), ...
            'soil:dom:missingBoundaryMemory','Open DOM state lacks cumulative export memory.');
        next.bottom_export_cumulative.doc_g_cm2=accepted.bottom_export_cumulative.doc_g_cm2 ...
            +docTransport.bottom_outflow_g_cm2;
        next.bottom_export_cumulative.don_g_cm2=accepted.bottom_export_cumulative.don_g_cm2 ...
            +donTransport.bottom_outflow_g_cm2;
        if parameters.boundary.schema_version==2
            for carrier={'doc','don'}
                key=carrier{1};ledger=docTransport;if strcmp(key,'don'),ledger=donTransport;end
                for route={'top_inflow_g_cm2','top_outflow_g_cm2','bottom_inflow_g_cm2'}
                    field=route{1};next.boundary_exchange_cumulative.(key).(field)= ...
                        accepted.boundary_exchange_cumulative.(key).(field)+ledger.(field);
                end
            end
        end
    end
    candidate=struct('state',next, ...
        'nh4_source_ug_cm3_s',zeroRate, ...
        'no3_source_ug_cm3_s',zeroRate, ...
        'mineralized_don_g_cm3',zeroRate, ...
        'transport',struct('doc',docTransport,'don',donTransport), ...
        'direct_don_to_no3',false, ...
        'reaction_ownership','explicit_microbial_cn_ledger');
end

function [next,docTransport,donTransport]=advanceSegment( ...
        accepted,thetaStart,theta,bulkDensity,dz_cm,q_cm_s,parameters,dt_s,bottom)
% Advance one liquid flow interval from thetaStart to theta (both physical liquid).
    dynamic=strcmp(parameters.sorption_model,'mend_dynamic');
    rewetting=thetaStart==0 & theta>0;
    docPhase=soil.dom.phaseState(accepted,'doc',thetaStart,bulkDensity,parameters,rewetting);
    donPhase=soil.dom.phaseState(accepted,'don',thetaStart,bulkDensity,parameters,rewetting);
    if dynamic
        % The existing free-solute ADE needs positive liquid capacity.
        % Do not change dry/frozen transport physics through a water floor.
        assert(all(theta(:)>0),'soil:dom:dynamicDryTransport', ...
            'Dynamic DOM transport requires positive liquid water at the end of the water tile.');
        docInput=docPhase.free_g_cm3_bulk;donInput=donPhase.free_g_cm3_bulk;
        kdDoc=0;kdDon=0;
    else
        docInput=accepted.doc_total_g_cm3;donInput=accepted.don_total_g_cm3;
        kdDoc=parameters.kd_doc_cm3_g;kdDon=parameters.kd_don_cm3_g;
    end
    docBoundary=bottom;donBoundary=bottom;
    if isfield(parameters,'boundary')
        docBoundary=soil.dom.materializeTransportBoundary(parameters.boundary,bottom,'doc');
        donBoundary=soil.dom.materializeTransportBoundary(parameters.boundary,bottom,'don');
    end
    [doc,docTransport]=soil.dom.transportConservative( ...
        docInput,theta,bulkDensity, ...
        kdDoc,dz_cm,q_cm_s, ...
        parameters.diffusion_cm2_s,parameters.dispersion_cm,dt_s, ...
        parameters.transport_substep_s,[],docBoundary,thetaStart);
    [don,donTransport]=soil.dom.transportConservative( ...
        donInput,theta,bulkDensity, ...
        kdDon,dz_cm,q_cm_s, ...
        parameters.diffusion_cm2_s,parameters.dispersion_cm,dt_s, ...
        parameters.transport_substep_s,[],donBoundary,thetaStart);
    if dynamic
        doc=doc+accepted.doc_sorbed_g_cm3;
        don=don+accepted.don_sorbed_g_cm3;
        % Public storage diagnostics retain their total-inventory meaning.
        sorbedC=sum(accepted.doc_sorbed_g_cm3.*dz_cm(:));
        sorbedN=sum(accepted.don_sorbed_g_cm3.*dz_cm(:));
        docTransport.storage_before_g_cm2=docTransport.storage_before_g_cm2+sorbedC;
        docTransport.storage_after_g_cm2=docTransport.storage_after_g_cm2+sorbedC;
        donTransport.storage_before_g_cm2=donTransport.storage_before_g_cm2+sorbedN;
        donTransport.storage_after_g_cm2=donTransport.storage_after_g_cm2+sorbedN;
    end
    assert(all(don>=0)&&all(isfinite(don)), ...
        'soil:dom:badDonCandidate','DON candidate is invalid.');
    next=accepted;
    next.doc_total_g_cm3=doc;
    next.don_total_g_cm3=don;
    next.direct_don_to_no3=false;
end

function [theta,q,dt]=normalizeSchedule(theta,q,dt,nl)
% Normalize legacy one-tile inputs or an explicit accepted-tile schedule.
    assert(isnumeric(dt)&&isvector(dt)&&~isempty(dt) ...
        &&all(isfinite(dt(:)))&&all(dt(:)>0), ...
        'soil:dom:transportScheduleTime', ...
        'DOM transport segment durations must be finite and positive.');
    dt=double(dt(:).');segments=numel(dt);
    if segments==1
        theta=theta(:);q=q(:);
    else
        assert(ismatrix(theta)&&ismatrix(q)&&size(theta,1)==nl ...
            &&size(q,1)==nl&&size(theta,2)==segments&&size(q,2)==segments, ...
            'soil:dom:transportScheduleShape', ...
            'Scheduled DOM theta and QL must be NL-by-number-of-segments.');
    end
    assert(size(theta,1)==nl&&size(q,1)==nl ...
        &&size(theta,2)==segments&&size(q,2)==segments ...
        &&all(isfinite(theta(:)))&&all(theta(:)>=0)&&all(isfinite(q(:))), ...
        'soil:dom:transportScheduleDriver', ...
        'Scheduled DOM theta and QL drivers are invalid.');
end

function total=accumulateTransport(total,part)
% Combine closed-transport diagnostics without discarding earlier segments.
    if isempty(total)
        total=part;
        return
    end
    total.net_layer_change_g_cm3=total.net_layer_change_g_cm3 ...
        +part.net_layer_change_g_cm3;
    total.storage_after_g_cm2=part.storage_after_g_cm2;
    total.bottom_outflow_g_cm2=total.bottom_outflow_g_cm2+part.bottom_outflow_g_cm2;
    for route={'top_inflow_g_cm2','top_outflow_g_cm2','bottom_inflow_g_cm2'}
        f=route{1};total.(f)=total.(f)+part.(f);
    end
    total.closure_residual_g_cm2=total.closure_residual_g_cm2 ...
        +part.closure_residual_g_cm2;
    total.accepted_substeps=total.accepted_substeps+part.accepted_substeps;
    total.min_substep_s=min(total.min_substep_s,part.min_substep_s);
    total.accepted_dt_s=[total.accepted_dt_s,part.accepted_dt_s];
    % Change: 2026-09-22 Wenjie Liu -- retain accepted source-label maps.
    total.source_label_tiles=[total.source_label_tiles,part.source_label_tiles];
    total.frozen_schedule_replay=total.frozen_schedule_replay ...
        &&part.frozen_schedule_replay;
end

function total=finalizeTransport(total,endpointChange,dt)
% Make the outer-step layer change authoritative and audit segment summation.
    summed=total.net_layer_change_g_cm3;
    scale=max(1,max(abs(endpointChange(:))));
    assert(max(abs(summed(:)-endpointChange(:)))<=1e-12*scale, ...
        'soil:dom:transportScheduleClosure', ...
        'Segment DOM changes do not match the accepted outer-step endpoints.');
    total.net_layer_change_g_cm3=endpointChange;
    total.driver_segments=numel(dt);
    total.driver_elapsed_s=sum(dt);
    total.segment_dt_s=dt;
end
