function [tag,diagnostic]=replayTransportLedger(tag,ledger,dz)
%REPLAYTRANSPORTLEDGER Bounded implicit mixing on immutable total exchanges.
% tag and total endpoints: bottom-up g N/cm3 bulk; dz: cm. Directed
% face exchanges: g N/cm2 per accepted tile (including symmetric dispersion).
% Solve (diag(m_new+outgoing)-incoming)*f_new = m_tag_old, then
% m_tag_new=m_new*f_new. With no external input, row sums equal m_old.
% Unlabelled external input adds row dominance but contributes no labelled RHS.
% Thus 0<=f<=1 and label mass plus external export is conserved, even when an explicit
% total operator is not monotone. This is an attribution discretization,
% not an identical explicit-operator replay. It never changes total states.
% No material clipping, global redistribution, or tracer-selected tiles.
% Editor: Wenjie Liu; Date: 2026-09-22.
% Change: 2026-09-22 Wenjie Liu -- retain already-solved implicit-fraction
% directed label exchanges on all layer faces, without changing label state.
    tag=tag(:);dz=dz(:);initial=sum(tag.*dz);
    assert(~isempty(tag)&&numel(tag)==numel(dz)&&all(isfinite(dz))&&all(dz>0), ...
        'tracer:transportShape','Positive layer thicknesses must align with labels.');
    assert(isfield(ledger,'source_label_tiles'), ...
        'tracer:missingTransportLedger','Accepted total transport ledger required.');
% Change: 2026-09-29 Wenjie Liu -- allow an immutable total bottom export;
% its implicit donor fraction defines the labelled export and terminal ledger.
% Change: 2026-09-29 Wenjie Liu -- organic external input is UNLABELLED;
% it dilutes the implicit fractions, never creates mineral-fertilizer labels.
    maxResidual=0;roundoffMass=0;bottomTagOut=0;topTagOut=0;
    upwardFace=zeros(numel(dz)+1,1);downwardFace=upwardFace;
    for k=1:numel(ledger.source_label_tiles)
        tile=ledger.source_label_tiles{k};
        old=tile.before(:).*dz;next=tile.after(:).*dz;
        mass=tag.*dz;n=numel(mass);fwd=tile.forward(:);back=tile.backward(:);
        assert(all(isfinite([old;next;mass;fwd;back]))&&all([old;next;fwd;back]>=0), ...
            'tracer:transportLedgerDomain','Invalid accepted mass/exchange.');
        tol=512*eps(max([old;next;realmin]));
        assert(all(mass>=-tol&mass<=old+tol), ...
            'tracer:transportDonorBound','Input label exceeds total donor.');
        incoming=[0;fwd]+[back;0];outgoing=[fwd;0]+[0;back];
        bottomOut=0;
        if isfield(tile,'bottom_outflow_g_cm2'),bottomOut=tile.bottom_outflow_g_cm2;end
        assert(isscalar(bottomOut)&&isfinite(bottomOut)&&bottomOut>=0, ...
            'tracer:transportBoundary','Invalid total bottom export.');
        outgoing(1)=outgoing(1)+bottomOut;
        boundaryIn=zeros(n,1);topOut=0;
        if isfield(tile,'top_outflow_g_cm2')
            topOut=tile.top_outflow_g_cm2;
            values=[topOut,tile.bottom_inflow_g_cm2,tile.top_inflow_g_cm2];
            assert(all(isfinite(values))&&all(values>=0),'tracer:transportBoundary', ...
                'Invalid organic external input/output ledger.');
            boundaryIn(1)=boundaryIn(1)+tile.bottom_inflow_g_cm2;
            boundaryIn(end)=boundaryIn(end)+tile.top_inflow_g_cm2;
            outgoing(end)=outgoing(end)+topOut;
        end
        incoming=incoming+boundaryIn;
        residual=next-old-incoming+outgoing;
        assert(all(abs(residual)<=tol), ...
            'tracer:transportLedgerClosure','Accepted total face exchanges do not close.');
        % old+incoming includes unlabelled input and gives row dominance.
        % Its difference from
        % next+outgoing is only the independently checked total roundoff.
        diagonal=old+incoming;
        empty=diagonal==0;diagonal(empty)=1;
        A=sparse(1:n,1:n,diagonal,n,n) ...
            +sparse(2:n,1:n-1,-fwd,n,n)+sparse(1:n-1,2:n,-back,n,n);
        fraction=A\mass;
        assert(all(isfinite(fraction))&&all(fraction>=-512*eps&fraction<=1+512*eps), ...
            'tracer:transportFraction','Implicit label fraction outside [0,1].');
        bounded=min(max(fraction,0),1);
        roundoffMass=roundoffMass+sum(next.*(bounded-fraction));
        result=next.*bounded;
        tagOut=bottomOut*bounded(1);
        tagTopOut=topOut*bounded(end);topTagOut=topTagOut+tagTopOut;
        bottomTagOut=bottomTagOut+tagOut;
        errorMass=sum(result)-sum(mass)+tagOut+tagTopOut;
        massTol=2048*eps(max([sum(abs(mass));sum(next);realmin]));
        assert(abs(errorMass)<=massTol, ...
            'tracer:transportMass','Implicit label transport does not conserve mass.');
        maxResidual=max(maxResidual,abs(errorMass));
        % Internal memory is bottom-up: forward i->i+1 is upward. These are
        % effective directed exchanges of the accepted implicit attribution,
        % not new explicit fluxes. The total ledger owns external export.
        upwardFace(2:n)=upwardFace(2:n)+fwd.*bounded(1:n-1);
        downwardFace(2:n)=downwardFace(2:n)+back.*bounded(2:n);
        downwardFace(1)=downwardFace(1)+tagOut;
        upwardFace(end)=upwardFace(end)+tagTopOut;
        tag=result./dz;
    end
    diagnostic=struct('method','implicit fractions on frozen directed total exchanges', ...
        'accepted_substeps',numel(ledger.source_label_tiles), ...
        'storage_before_g_cm2',initial,'storage_after_g_cm2',sum(tag.*dz), ...
        'closure_residual_g_cm2',sum(tag.*dz)-initial+bottomTagOut+topTagOut, ...
        'bottom_outflow_g_cm2',bottomTagOut, ...
        'top_outflow_g_cm2',topTagOut, ...
        'max_tile_mass_residual_g_cm2',maxResidual, ...
        'fraction_roundoff_mass_g_cm2',roundoffMass, ...
        'upward_face_amount_g_cm2',upwardFace, ...
        'downward_face_amount_g_cm2',downwardFace, ...
        'net_upward_face_amount_g_cm2',upwardFace-downwardFace, ...
        'face_grid_order','bottom_up_including_external_faces', ...
        'face_exchange_definition','effective_implicit_fraction_directed_exchange');
end
