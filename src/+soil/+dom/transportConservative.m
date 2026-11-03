function [totalNext, diagnostic] = transportConservative( ...
        total, theta, bulkDensity, kd, dz_cm, q_cm_s, ...
        diffusion_cm2_s, dispersion_cm, dt_s, maxSubstep_s, frozenSubsteps_s, bottom,thetaStart)
%TRANSPORTCONSERVATIVE Finite-volume DOM ADE with explicit boundary ledger.
%
%   Element order is COMPASS bottom-up. q_cm_s is the accepted liquid flux at
%   each element centre, positive deeper-to-shallower. Internal-face advection
%   is upwind; dispersion/diffusion is centred. Only internal faces exchange
%   mass by default (closed); optional bottom outflow is specified below.
%   q_cm_s and theta are the flux and final liquid content of one accepted
%   water tile. Optional thetaStart [NL x 1, cm3 liquid/cm3 bulk] makes
%   liquid storage evolve linearly within the tile. Explicit donor and
%   diffusive concentrations use the beginning of each accepted substep.
%   A prescribed constant theta is supported by omitting thetaStart.
%
%   Editor: Jizixing Zhu; Date: 2026-08-21.
%   Change: 2026-09-02 Codex -- accept depth-resolved Kd and dispersivity;
%   internal-face dispersivity is the arithmetic mean of adjacent elements.
%   Change: 2026-09-08 Jizixing Zhu -- expose accepted per-layer DOM transport
%     change [g/cm3 bulk soil] for independent reaction/transport budget audits.
%   Change: 2026-09-14 Jizixing Zhu -- reject negative trial inventories
%     instead of clipping them; absolute clipping breaks paired C/N transport.
%   Change: 2026-09-21 Wenjie Liu -- document the matched accepted-endpoint
%     time-level contract used by chronological outer-step replay.
%   Change: 2026-09-21 Wenjie Liu -- expose every accepted internal time
%     tile and permit a passive source tracer to replay that frozen schedule.
%     The replay can neither choose a different tile nor feed back to the
%     total-DOM transport solution.
%   Change: 2026-09-22 Wenjie Liu -- export accepted total endpoints and
%     directed exchange masses for bounded one-way label attribution.
%   Change: 2026-09-29 Wenjie Liu -- optional aqueous bottom advection with
%     zero external diffusion; loss participates in positivity and the ledger.
%     bottom.kind=closed (default) or advective_outflow with accepted
%     q_liquid_cm_s<=0 (upward/into-column positive). Top stays closed.
%   Change: 2026-09-29 Wenjie Liu -- configured top/bottom concentration or
%     inflow/outflow boundaries; all four directed mass routes enter closure.
%   Change: 2026-11-03 Wenjie Liu -- synchronize free concentration and
%     diffusive conductance with liquid capacity on the explicit flux step.
%   Change: 2026-11-03 Codex -- interpolate centre water flux to the physical
%     internal face on unequal cells, including the rewetting donor check.
    if nargin<12||isempty(bottom),bottom=struct('kind','closed');end
    assert(isstruct(bottom)&&isscalar(bottom)&&isfield(bottom,'kind') ...
        &&ismember(bottom.kind,{'closed','advective_outflow','configured'}), ...
        'soil:dom:transportBoundary','Invalid DOM transport boundary.');
    openBottom=strcmp(bottom.kind,'advective_outflow');qBottom=0;
    if openBottom
        assert(isfield(bottom,'q_liquid_cm_s')&&isscalar(bottom.q_liquid_cm_s) ...
            &&isfinite(bottom.q_liquid_cm_s),'soil:dom:bottomWaterFlux', ...
            'Open DOM requires an explicit accepted liquid bottom-face flux.');
        qBottom=bottom.q_liquid_cm_s;
        assert(qBottom<=0,'soil:dom:bottomInflowUnsupported', ...
            'External DOC/DON concentration is required for bottom inflow.');
    end
    bottomOut=0;bottomIn=0;topOut=0;topIn=0;
    if nargin<11
        frozenSubsteps_s=[];
    end
    total=total(:); theta=theta(:); bulkDensity=bulkDensity(:);
    if nargin<13,thetaStart=theta;else,thetaStart=thetaStart(:);end
    dz_cm=dz_cm(:); q_cm_s=q_cm_s(:);
    n=numel(total);
    kd=expandProfile(kd,n,'kd');
    dispersion_cm=expandProfile(dispersion_cm,n,'dispersion_cm');
    assert(n>=1 && numel(theta)==n && numel(thetaStart)==n && numel(bulkDensity)==n ...
        && numel(dz_cm)==n && numel(q_cm_s)==n, ...
        'soil:dom:transportSize','DOM transport arrays must have NL values.');
    assert(all(isfinite(total))&&all(total>=0)&&all(isfinite(theta))&&all(theta>=0) ...
        &&all(isfinite(thetaStart))&&all(thetaStart>=0) ...
        && all(bulkDensity>0)&&all(isfinite(dz_cm)&dz_cm>0)&&all(isfinite(q_cm_s)) ...
        && diffusion_cm2_s>=0&&all(dispersion_cm>=0) ...
        &&dt_s>0&&maxSubstep_s>0, ...
        'soil:dom:transportValues','DOM transport inputs are invalid.');
    beforeProfile=total;
    waterChange=theta-thetaStart;
    rewetting=thetaStart==0 & theta>0;
    % The face is half a cell width away from either centre. Preserve the
    % original arithmetic mean exactly wherever neighbouring widths agree.
    qLower=q_cm_s(1:end-1);qUpper=q_cm_s(2:end);
    dzLower=dz_cm(1:end-1);dzUpper=dz_cm(2:end);
    qInterface=(qLower+qUpper)/2;
    unequal=dzLower~=dzUpper;
    fraction=dzLower(unequal)./(dzLower(unequal)+dzUpper(unequal));
    qInterface(unequal)=qLower(unequal)+fraction.*(qUpper(unequal)-qLower(unequal));
    if any(rewetting)
        outgoing=[max(qInterface,0);0]+[0;max(-qInterface,0)];
        if ~strcmp(bottom.kind,'closed')
            outgoing(1)=outgoing(1)+max(-bottom.q_liquid_cm_s,0);
            if isfield(bottom,'q_top_cm_s')
                outgoing(end)=outgoing(end)+max(bottom.q_top_cm_s,0);
            end
        end
        assert(all(outgoing(rewetting)==0),'soil:dom:dryDonor', ...
            'An initially dry cell cannot supply outgoing liquid on the same fixed-flux tile.');
    end
    assert(~openBottom||qBottom==0||(thetaStart(1)>0&&theta(1)>0),'soil:dom:dryBottomOutflow', ...
        'Nonzero liquid drainage cannot advect solute from a dry bottom cell.');
    before=sum(total.*dz_cm);
    if isempty(frozenSubsteps_s)
        replay=false;
        requestedSchedule=[];
    else
        replay=true;
        requestedSchedule=double(frozenSubsteps_s(:).');
        assert(~isempty(requestedSchedule)&&all(isfinite(requestedSchedule)) ...
            &&all(requestedSchedule>0) ...
            &&abs(sum(requestedSchedule)-dt_s)<=1e-10*max(1,dt_s), ...
            'soil:dom:transportFrozenSchedule', ...
            'Frozen DOM transport substeps must be positive and sum to dt_s.');
    end
    elapsed=0; acceptedSubsteps=0; minSubstep=Inf; acceptedDt=zeros(1,0);
    labelTiles=cell(1,0);
    while elapsed < dt_s-1e-12
        theta=thetaStart+(elapsed/dt_s)*waterChange;
        [boundaryInput,boundaryLoss]=soil.dom.boundaryCoefficients( ...
            bottom,theta,dz_cm,diffusion_cm2_s,dispersion_cm,rewetting);
        if replay
            assert(acceptedSubsteps<numel(requestedSchedule), ...
                'soil:dom:transportFrozenSchedule', ...
                'Frozen DOM transport schedule ended before dt_s.');
            trialDt=requestedSchedule(acceptedSubsteps+1);
        else
            trialDt=min(maxSubstep_s,dt_s-elapsed);
        end
        while true
            phase=soil.dom.equilibrate(total,theta,bulkDensity,kd,rewetting);
            concentration=phase.aqueous_g_cm3_water;
            flux=zeros(max(n-1,0),1);
            mixing=zeros(max(n-1,0),1);
            for face=1:n-1
                qFace=qInterface(face);
                if qFace>=0, upwind=concentration(face); ...
                else, upwind=concentration(face+1); end
                thetaFace=max(0.5*(theta(face)+theta(face+1)),0);
                distance=0.5*(dz_cm(face)+dz_cm(face+1));
                dispersionFace=0.5*(dispersion_cm(face)+dispersion_cm(face+1));
                dispersive=thetaFace*diffusion_cm2_s ...
                    +dispersionFace*abs(qFace);
                flux(face)=qFace*upwind ...
                    +dispersive*(concentration(face)-concentration(face+1))/distance;
                grossForward=(max(qFace,0)+dispersive/distance)*concentration(face);
                grossBackward=(max(-qFace,0)+dispersive/distance)*concentration(face+1);
                mixing(face)=min(grossForward,grossBackward);
            end
            candidate=total;
            if n>1
                transfer=flux*trialDt;
                candidate(1:end-1)=candidate(1:end-1)-transfer./dz_cm(1:end-1);
                candidate(2:end)=candidate(2:end)+transfer./dz_cm(2:end);
            end
            loss=boundaryLoss(1)*concentration(1)*trialDt;
            gain=boundaryInput(1)*trialDt;
            topLoss=boundaryLoss(2)*concentration(end)*trialDt;
            topGain=boundaryInput(2)*trialDt;
            if ~strcmp(bottom.kind,'closed')
                candidate(1)=candidate(1)+(gain-loss)/dz_cm(1);
                candidate(end)=candidate(end)+(topGain-topLoss)/dz_cm(end);
            end
            if all(candidate>=0) && all(isfinite(candidate))
                labelTiles{end+1}=struct('before',total,'after',candidate, ...
                    'forward',(max(flux,0)+mixing)*trialDt, ...
                    'backward',(max(-flux,0)+mixing)*trialDt, ...
                    'bottom_outflow_g_cm2',loss,'bottom_inflow_g_cm2',gain, ...
                    'top_outflow_g_cm2',topLoss,'top_inflow_g_cm2',topGain); %#ok<AGROW>
                bottomOut=bottomOut+loss;
                bottomIn=bottomIn+gain;topOut=topOut+topLoss;topIn=topIn+topGain;
                total=candidate;
                break;
            end
            if replay
                error('soil:dom:transportFrozenReplayNegative', ...
                    ['Passive-tracer replay became negative on an accepted ', ...
                     'total-DOM transport tile.']);
            end
            trialDt=trialDt/2;
            if trialDt<1e-6
                error('soil:dom:transportNonnegative', ...
                    'DOM transport could not find a nonnegative substep.');
            end
        end
        elapsed=elapsed+trialDt;
        acceptedSubsteps=acceptedSubsteps+1;
        minSubstep=min(minSubstep,trialDt);
        acceptedDt(end+1)=trialDt; %#ok<AGROW>
    end
    if replay
        assert(acceptedSubsteps==numel(requestedSchedule), ...
            'soil:dom:transportFrozenSchedule', ...
            'Frozen DOM transport schedule contains unused substeps.');
    end
    after=sum(total.*dz_cm);
    closure=after-before+bottomOut+topOut-bottomIn-topIn;
    closureTolerance=512*eps(max([before,after,bottomOut,topOut,bottomIn,topIn,realmin]))*max(1,acceptedSubsteps);
    assert(abs(closure)<=closureTolerance,'soil:dom:transportMassClosure', ...
        'DOM storage change plus outputs minus inputs does not close.');
    totalNext=total;
    diagnostic=struct('net_layer_change_g_cm3',total-beforeProfile, ...
        'storage_before_g_cm2',before, ...
        'storage_after_g_cm2',after, ...
        'closure_residual_g_cm2',closure, ...
        'bottom_outflow_g_cm2',bottomOut, ...
        'bottom_inflow_g_cm2',bottomIn,'top_outflow_g_cm2',topOut,'top_inflow_g_cm2',topIn, ...
        'accepted_substeps',acceptedSubsteps, ...
        'min_substep_s',minSubstep, ...
        'accepted_dt_s',acceptedDt, ...
        'frozen_schedule_replay',replay, ...
        'boundary_condition','closed_zero_mass_flux');
    diagnostic.source_label_tiles=labelTiles;
    if openBottom,diagnostic.boundary_condition='aqueous_advective_outflow_zero_diffusion';end
    if strcmp(bottom.kind,'configured'),diagnostic.boundary_condition='configured_organic_FV_faces_v2';end
end

function value=expandProfile(value,n,name)
% Normalize a scalar or element profile without changing its physical value.
    assert(isnumeric(value)&&all(isfinite(value(:)))&&all(value(:)>=0) ...
        &&(isscalar(value)||numel(value)==n), ...
        'soil:dom:transportProfile', ...
        '%s must be a finite nonnegative scalar or NL-element profile.',name);
    if isscalar(value), value=repmat(double(value),n,1); ...
    else, value=double(value(:)); end
end
