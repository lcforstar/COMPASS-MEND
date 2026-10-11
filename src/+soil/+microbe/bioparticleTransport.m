function candidate=bioparticleTransport(microbes,organicN,thetaBottom, ...
        bulkDensityBottom,dzBottom,qBottom,parameters,dt_s,frozenTransport,bottomFaces)
%BIOPARTICLETRANSPORT Transport active/dormant microbial C/N between phases.
%
%   Mobile biomass follows an advection-dispersion operator (closed by default);
%   attached/mobile exchange and irreversible straining follow the PFLOTRAN
%   BIOPARTICLE transport-only source structure. Growth, respiration, death,
%   and dormancy are excluded. C and N use the identical linear operator and
%   are independently conserved. State units remain MEMC [mg g^-1] for C and
%   COMPASS [g cm^-3] for N.
%   The irreversible retained store follows the first-order straining term of
%   Gargiulo et al. (2007), J. Contam. Hydrol. 92, 255--273,
%   https://doi.org/10.1016/j.jconhyd.2007.01.009. PFLOTRAN BIOPARTICLE is
%   retained as the attachment/detachment transport interface reference.
%   thetaBottom, qBottom and dt_s may describe a chronological accepted-tile
%   schedule (NL-by-K, NL-by-K and 1-by-K); each flux is paired with theta at
%   the same converged soil-tile endpoint.
%   Aggregated carbon_transport/nitrogen_transport storage fields cover the
%   complete mobile+attached+retained inventory. closure_residual_g_cm2 is the
%   total-phase endpoint change plus bottom export; the summed segment-only ADE
%   residual is retained separately as segment_ade_closure_residual_g_cm2.
%
% Editor: Jizixing Zhu; Date: 2026-08-28.
% Change: 2026-08-30 Codex -- document the separate literature ownership of
%   irreversible straining and enable its nonzero site-calibration path.
% Change: 2026-09-14 Jizixing Zhu -- carry explicitly stored intracellular
%   enzyme phases with the same cell-particle C/N transport/exchange operator.
% Change: 2026-09-21 Wenjie Liu -- replay all accepted water tiles in order
%   and aggregate ADE/exchange diagnostics instead of extending the final QL
%   tile over the full outer duration.
% Change: 2026-09-21 Wenjie Liu -- define aggregate transport storage on the
%   complete mobile+attached+retained phase inventory and refresh the derived
%   Tot cache for both current and legacy (no ENZN) microbial states.
% Change: 2026-09-21 Wenjie Liu -- permit passive source labels to replay
%   the accepted total-N ADE tiles without selecting a second time grid.
% Change: 2026-09-22 Wenjie Liu -- expose chronological total-N transport
%   and phase-exchange ledgers; the physical C/N calculation is unchanged.

% Change: 2026-09-29 Wenjie Liu -- optional accepted bottom liquid faces for
% mobile C/N and intracellular cargo, with paired monotone clock and exports.
% Attached/retained phases have no direct external flux. C is top-down mg/g,
% N bottom-up g/cm3; all export ledgers are positive outward g C/N per cm2.
% Change: 2026-09-29 Wenjie Liu -- numbered organic boundaries with explicit
% paired external C/N; all input/output routes are ledgered without DIN edits.
    if nargin<9,frozenTransport=[];end
    if nargin<10,bottomFaces={};end
    n=numel(dzBottom); bd=bulkDensityBottom(:); dz=dzBottom(:);
    [theta,q,dt]=normalizeSchedule(thetaBottom,qBottom,dt_s,n);
    openBottom=isfield(parameters,'boundary') ...
        &&~strcmp(parameters.boundary.bottom_kind,'closed');
    if openBottom
        assert(iscell(bottomFaces)&&numel(bottomFaces)==numel(dt), ...
            'bioparticle:missingBottomSchedule','Every accepted particle tile requires its bottom liquid face.');
        assert(isfield(microbes,'bottom_export_cumulative') ...
            &&microbes.bottom_export_cumulative.schema_version==1, ...
            'bioparticle:missingExportMemory','Open particle state lacks cumulative export memory.');
    else
        assert(isempty(bottomFaces),'bioparticle:closedBottomDriver', ...
            'Closed particle boundary must not consume open drivers.');
        bottomFaces=repmat({struct('kind','closed')},size(dt));
    end
    assert(numel(bd)==n, ...
        'bioparticleTransport:shape','BIOPARTICLE driver shapes must equal NL.');
    % Change: 2026-09-11 Jizixing Zhu -- reaction owner must already close
    % phase sums. Transport cannot repair them by inventing phase transfers.
    groups={'MB','BD'};
    if isfield(microbes,'ENZN')
        assert(isfield(organicN,'ENZN_N'),'bioparticleTransport:enzymeOwner', ...
            'Intracellular enzyme C requires its paired N owner.');
        groups{end+1}='ENZN';
    end
    for key=groups
        k=key{1};
        csum=microbes.([k '_mobile'])+microbes.([k '_attached'])+microbes.([k '_retained']);
        nsum=organicN.([k '_mobile_N'])+organicN.([k '_attached_N'])+organicN.([k '_retained_N']);
        assert(all(abs(csum-microbes.(k))<=1e-11.*max(1,abs(microbes.(k)))) ...
            &&all(abs(nsum-organicN.([k '_N']))<=1e-13.*max(1,abs(organicN.([k '_N'])))), ...
            'bioparticleTransport:phaseClosure','Reaction phase sums must match C/N aggregates.');
    end
    beforeC=columnCarbon(microbes,bd,dz);
    beforeN=columnNitrogen(organicN,dz);
    diagnostics=struct();exportC=0;exportN=0;commonClock=[];otherNet=zeros(1,2);
    for prefix=groups
        name=prefix{1};
        c=struct('mobile',flip(microbes.([name '_mobile'])).*bd./1000, ...
            'attached',flip(microbes.([name '_attached'])).*bd./1000, ...
            'retained',flip(microbes.([name '_retained'])).*bd./1000);
        nn=struct('mobile',organicN.([name '_mobile_N']), ...
            'attached',organicN.([name '_attached_N']), ...
            'retained',organicN.([name '_retained_N']));
        frozen=[];
        if ~isempty(frozenTransport)
            assert(isstruct(frozenTransport)&&isfield(frozenTransport,name) ...
                &&isfield(frozenTransport.(name),'nitrogen_transport') ...
                &&isfield(frozenTransport.(name).nitrogen_transport,'accepted_dt_s'), ...
                'bioparticleTransport:frozenSchedule', ...
                'Frozen BIOPARTICLE schedule lacks %s nitrogen tiles.',name);
            frozen=frozenTransport.(name).nitrogen_transport.accepted_dt_s;
        end
        groupParameters=parameters;groupParameters.boundary_carrier=name;
        [c,nn,d]=advancePairSchedule(c,nn,theta,bd,dz,q,groupParameters,dt,frozen,bottomFaces);
        if openBottom
            clock=d.carbon_transport.accepted_dt_s;
            if isempty(commonClock),commonClock=clock;end
            assert(isequal(clock,commonClock),'bioparticle:commonClock', ...
                'All mobile cell and enzyme cargo must share an accepted clock.');
            for element={'carbon','nitrogen'}
                e=element{1};field=[e '_g_cm2'];
                microbes.bottom_export_cumulative.(name).(field)= ...
                    microbes.bottom_export_cumulative.(name).(field) ...
                    +d.([e '_transport']).bottom_outflow_g_cm2;
                if parameters.boundary.schema_version==2
                    for route={'top_inflow_g_cm2','top_outflow_g_cm2','bottom_inflow_g_cm2'}
                        r=route{1};key=[name '_' e];
                        microbes.boundary_exchange_cumulative.(key).(r)= ...
                            microbes.boundary_exchange_cumulative.(key).(r)+d.([e '_transport']).(r);
                    end
                end
            end
        end
        elements={'carbon','nitrogen'};
        for e=1:2
            ledger=d.([elements{e} '_transport']);
            otherNet(e)=otherNet(e)+ledger.top_outflow_g_cm2 ...
                -ledger.top_inflow_g_cm2-ledger.bottom_inflow_g_cm2;
        end
        exportC=exportC+d.carbon_transport.bottom_outflow_g_cm2;
        exportN=exportN+d.nitrogen_transport.bottom_outflow_g_cm2;
        microbes.([name '_mobile'])=flip(c.mobile./bd.*1000);
        microbes.([name '_attached'])=flip(c.attached./bd.*1000);
        microbes.([name '_retained'])=flip(c.retained./bd.*1000);
        organicN.([name '_mobile_N'])=nn.mobile;
        organicN.([name '_attached_N'])=nn.attached;
        organicN.([name '_retained_N'])=nn.retained;
        microbes.(name)=microbes.([name '_mobile']) ...
            +microbes.([name '_attached'])+microbes.([name '_retained']);
        organicN.([name '_N'])=organicN.([name '_mobile_N']) ...
            +organicN.([name '_attached_N'])+organicN.([name '_retained_N']);
        diagnostics.(name)=d;
    end
    microbes=refreshTotalCache(microbes);
    afterC=columnCarbon(microbes,bd,dz); afterN=columnNitrogen(organicN,dz);
    residualC=afterC-beforeC+exportC+otherNet(1);
    residualN=afterN-beforeN+exportN+otherNet(2);
    assert(abs(residualC)<=1e-12*max(1,beforeC) ...
        &&abs(residualN)<=1e-12*max(1,beforeN), ...
        'bioparticleTransport:closure', ...
        'BIOPARTICLE C/N storage plus boundary outputs minus inputs does not close.');
    candidate=struct('microbes',microbes,'organic_n',organicN, ...
        'diagnostics',diagnostics,'carbon_residual_g_cm2',residualC, ...
        'nitrogen_residual_g_cm2',residualN, ...
        'bottom_carbon_outflow_g_cm2',exportC,'bottom_nitrogen_outflow_g_cm2',exportN, ...
        'reaction_terms_disabled',true,'elapsed_s',sum(dt));
end

function [c,n,d]=advancePairSchedule(c,n,theta,bd,dz,q,p,dt,frozen,bottomFaces)
% Replay one phase pair through chronological accepted soil-water tiles.
    if nargin<9,frozen=[];end
    cBefore=c.mobile+c.attached+c.retained;
    nBefore=n.mobile+n.attached+n.retained;
    dc=[];dn=[];exchangeSubsteps=0;frozenOffset=0;labelSegments=cell(1,0);
    for segment=1:numel(dt)
        frozenPart=[];
        if ~isempty(frozen)
            [frozenPart,frozenOffset]=takeFrozenSegment( ...
                frozen,frozenOffset,dt(segment));
        end
        [c,n,part]=advancePair(c,n,theta(:,segment),bd,dz, ...
            q(:,segment),p,dt(segment),frozenPart,bottomFaces{segment});
        dc=accumulateTransport(dc,part.carbon_transport);
        dn=accumulateTransport(dn,part.nitrogen_transport);
        exchangeSubsteps=exchangeSubsteps+part.exchange_substeps;
        labelSegments{end+1}=struct('transport',part.nitrogen_transport, ...
            'exchange',{part.source_label_exchange}); %#ok<AGROW>
    end
    assert(isempty(frozen)||frozenOffset==numel(frozen), ...
        'bioparticleTransport:frozenSchedule', ...
        'Frozen BIOPARTICLE schedule contains unused transport tiles.');
    dc=finalizePhaseTransport(dc,cBefore, ...
        c.mobile+c.attached+c.retained,dz,dt,'carbon');
    dn=finalizePhaseTransport(dn,nBefore, ...
        n.mobile+n.attached+n.retained,dz,dt,'nitrogen');
    d=struct('carbon_transport',dc,'nitrogen_transport',dn, ...
        'exchange_substeps',exchangeSubsteps,'active_motility',false, ...
        'driver_segments',numel(dt),'segment_dt_s',dt);
    d.source_label_segments=labelSegments;
end

function [c,n,d]=advancePair(c,n,theta,bd,dz,q,p,dt,frozen,bottom)
% Advance one active/dormant phase pair with separate C and N ledgers.
% Editor: Jizixing Zhu; Date: 2026-08-28.
    if nargin<9,frozen=[];end
    maxStep=p.transport_substep_s;
    cBoundary=bottom;nBoundary=bottom;
    if isfield(p,'boundary')&&p.boundary.schema_version==2
        cBoundary=soil.dom.materializeTransportBoundary(p.boundary,bottom,[p.boundary_carrier '_carbon']);
        nBoundary=soil.dom.materializeTransportBoundary(p.boundary,bottom,[p.boundary_carrier '_nitrogen']);
    end
    if ~strcmp(bottom.kind,'closed')
        maxStep=soil.microbe.mobileTransportStep(theta,dz,q,p,cBoundary);
    end
    [c.mobile,dc]=soil.dom.transportConservative(c.mobile,theta,bd,0,dz,q, ...
        p.diffusion_cm2_s,p.dispersion_cm,dt,maxStep,frozen,cBoundary);
    [n.mobile,dn]=soil.dom.transportConservative(n.mobile,theta,bd,0,dz,q, ...
        p.diffusion_cm2_s,p.dispersion_cm,dt,maxStep,frozen,nBoundary);
    if ~strcmp(bottom.kind,'closed')
        assert(isequal(dc.accepted_dt_s,dn.accepted_dt_s),'bioparticle:commonClock', ...
            'Paired mobile C/N must use the same accepted transport operator.');
    end
    elapsed=0; count=0;labelExchange=cell(1,0);
    while elapsed<dt-1e-12
        step=min(p.transport_substep_s,dt-elapsed);
        attach=p.attachment_rate_s.*c.mobile.*step;
        detach=p.detachment_rate_s.*c.attached.*step;
        strain=p.straining_rate_s.*c.mobile.*step;
        mobileOut=attach+strain;
        scale=min(1,c.mobile./max(mobileOut,realmin));
        attach=attach.*scale; strain=strain.*scale;
        detach=min(detach,c.attached);
        ratioMobile=safeRatio(n.mobile,c.mobile);
        ratioAttached=safeRatio(n.attached,c.attached);
        nAttach=min(n.mobile,attach.*ratioMobile);
        nStrain=min(n.mobile-nAttach,strain.*ratioMobile);
        nDetach=min(n.attached,detach.*ratioAttached);
        exchangeBefore=n;
        c.mobile=c.mobile-attach-strain+detach;
        c.attached=c.attached+attach-detach; c.retained=c.retained+strain;
        n.mobile=n.mobile-nAttach-nStrain+nDetach;
        n.attached=n.attached+nAttach-nDetach; n.retained=n.retained+nStrain;
        labelExchange{end+1}=struct('before',exchangeBefore,'after',n, ...
            'attach',nAttach,'detach',nDetach,'strain',nStrain); %#ok<AGROW>
        elapsed=elapsed+step; count=count+1;
    end
    d=struct('carbon_transport',dc,'nitrogen_transport',dn, ...
        'exchange_substeps',count,'active_motility',false);
    d.source_label_exchange=labelExchange;
end

function [part,nextOffset]=takeFrozenSegment(schedule,offset,target)
% Extract the accepted inner tiles belonging to one outer water tile.
    schedule=double(schedule(:).');part=zeros(1,0);elapsed=0;
    while elapsed<target-1e-12
        index=offset+numel(part)+1;
        assert(index<=numel(schedule), ...
            'bioparticleTransport:frozenSchedule', ...
            'Frozen BIOPARTICLE schedule ended before an outer tile.');
        part(end+1)=schedule(index); %#ok<AGROW>
        elapsed=elapsed+part(end);
        assert(elapsed<=target+1e-10*max(1,target), ...
            'bioparticleTransport:frozenSchedule', ...
            'Frozen BIOPARTICLE schedule crosses an outer-tile boundary.');
    end
    assert(abs(elapsed-target)<=1e-10*max(1,target), ...
        'bioparticleTransport:frozenSchedule', ...
        'Frozen BIOPARTICLE schedule does not cover the outer tile exactly.');
    nextOffset=offset+numel(part);
end

function [theta,q,dt]=normalizeSchedule(theta,q,dt,nl)
% Normalize legacy one-tile inputs or an explicit accepted-tile schedule.
    assert(isnumeric(dt)&&isvector(dt)&&~isempty(dt) ...
        &&all(isfinite(dt(:)))&&all(dt(:)>0), ...
        'bioparticleTransport:scheduleTime', ...
        'BIOPARTICLE transport segment durations must be finite and positive.');
    dt=double(dt(:).');segments=numel(dt);
    if segments==1
        theta=theta(:);q=q(:);
    else
        assert(ismatrix(theta)&&ismatrix(q)&&size(theta,1)==nl ...
            &&size(q,1)==nl&&size(theta,2)==segments&&size(q,2)==segments, ...
            'bioparticleTransport:scheduleShape', ...
            'Scheduled BIOPARTICLE theta and QL must be NL-by-segments.');
    end
    assert(size(theta,1)==nl&&size(q,1)==nl ...
        &&size(theta,2)==segments&&size(q,2)==segments ...
        &&all(isfinite(theta(:)))&&all(theta(:)>=0)&&all(isfinite(q(:))), ...
        'bioparticleTransport:scheduleDriver', ...
        'Scheduled BIOPARTICLE theta and QL drivers are invalid.');
end

function total=accumulateTransport(total,part)
% Preserve all ADE segment changes and endpoint/storage audit fields.
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
    total.source_label_tiles=[total.source_label_tiles,part.source_label_tiles];
    total.frozen_schedule_replay=total.frozen_schedule_replay ...
        &&part.frozen_schedule_replay;
end

function total=finalizePhaseTransport(total,beforeProfile,afterProfile,dz,dt,kind)
% Define aggregate storage on all phases; exchange is internal to each layer.
    summedChange=total.net_layer_change_g_cm3;
    endpointChange=afterProfile-beforeProfile;
    scale=max(1,max(abs([summedChange(:);endpointChange(:)])));
    assert(max(abs(summedChange-endpointChange))<=1e-12*scale, ...
        'bioparticleTransport:aggregateLayerBalance', ...
        'Aggregated %s ADE change does not match the total-phase endpoint.',kind);
    segmentClosure=total.closure_residual_g_cm2;
    total.net_layer_change_g_cm3=endpointChange;
    total.storage_before_g_cm2=sum(beforeProfile.*dz);
    total.storage_after_g_cm2=sum(afterProfile.*dz);
    total.closure_residual_g_cm2=total.storage_after_g_cm2 ...
        -total.storage_before_g_cm2+total.bottom_outflow_g_cm2 ...
        +total.top_outflow_g_cm2-total.top_inflow_g_cm2-total.bottom_inflow_g_cm2;
    total.segment_ade_closure_residual_g_cm2=segmentClosure;
    total.layer_balance_residual_g_cm2=sum((summedChange-endpointChange).*dz);
    total.storage_basis='mobile_attached_retained_total';
    total.phase_exchange_is_internal=true;
    total.driver_segments=numel(dt);
    total.driver_elapsed_s=sum(dt);
    total.segment_dt_s=dt;
end

function microbes=refreshTotalCache(microbes)
% Refresh derived total C without requiring ENZN in legacy state schemas.
    if ~isfield(microbes,'Tot')&&~isfield(microbes,'ENZN'),return;end
    required={'POC','MOC','MB','BD','EP','EM'};
    for i=1:numel(required)
        assert(isfield(microbes,required{i}), ...
            'bioparticleTransport:totalCacheSchema', ...
            'Microbial Tot cache requires pool %s.',required{i});
    end
    total=microbes.POC+microbes.MOC+microbes.MB+microbes.BD ...
        +microbes.EP+microbes.EM;
    if isfield(microbes,'ENZN'),total=total+microbes.ENZN;end
    microbes.Tot=total;
end


function value=columnCarbon(s,bdBottom,dz)
% Active plus dormant microbial C column [g C cm^-2].
% Editor: Jizixing Zhu; Date: 2026-08-28.
% Change: 2026-09-14 Jizixing Zhu -- include intracellular enzyme cargo.
    value=sum(flip(s.MB+s.BD).*bdBottom(:)./1000.*dz(:));
    if isfield(s,'ENZN'),value=value+sum(flip(s.ENZN).*bdBottom(:)./1000.*dz(:));end
end

function value=columnNitrogen(s,dz)
% Active plus dormant microbial N column [g N cm^-2].
% Editor: Jizixing Zhu; Date: 2026-08-28.
% Change: 2026-09-14 Jizixing Zhu -- include intracellular enzyme cargo N.
    value=sum((s.MB_N+s.BD_N).*dz(:));
    if isfield(s,'ENZN_N'),value=value+sum(s.ENZN_N.*dz(:));end
end

function ratio=safeRatio(n,c)
% Phase N:C with exact zero-carbon behavior.
% Editor: Jizixing Zhu; Date: 2026-08-28.
    ratio=zeros(size(c)); active=c>0; ratio(active)=n(active)./c(active);
end
