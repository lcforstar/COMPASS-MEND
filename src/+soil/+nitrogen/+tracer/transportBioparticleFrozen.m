function [tagNext,diagnostic]=transportBioparticleFrozen(tag,~,~,~,~,dz,~,totalTransport,totalAfter)
%TRANSPORTBIOPARTICLEFROZEN Attribute accepted microbial N transport/exchange.
% Label/total: bottom-up g N/cm3 bulk; dz: cm. Exact total-N phase transfers
% use donor fractions, interleaved with bounded implicit mobile transport.
% Editor: Wenjie Liu; Date: 2026-09-21.
% Change: 2026-09-22 Wenjie Liu -- consume immutable accepted total ledgers;
% do not apply a nonpositive explicit operator to a synthetic carbon carrier.
% Change: 2026-09-22 Wenjie Liu -- retain existing labelled attachment,
% detachment, straining and effective transport face amounts for diagnostics.
% Change: 2026-09-29 Wenjie Liu -- mobile bottom export is a terminal label loss.
% Change: 2026-09-29 Wenjie Liu -- organic top exports are terminal too;
% external unlabelled input dilution is owned by the frozen transport replay.
    groups={'MB','BD','ENZN'};phases={'mobile','attached','retained'};
    tagNext=tag;diagnostic=struct();dz=dz(:);
    for i=1:numel(groups)
        group=groups{i};state=struct();
        for j=1:3
            state.(phases{j})=tag.([group '_' phases{j} '_N']);
        end
        initial=sum((state.mobile+state.attached+state.retained).*dz);
        assert(isfield(totalTransport.(group),'source_label_segments'), ...
            'tracer:missingParticleLedger','Accepted particle ledger required.');
        segments=totalTransport.(group).source_label_segments;parts=cell(size(segments));
        attachLabel=zeros(size(dz));detachLabel=attachLabel;strainLabel=attachLabel;
        upwardFace=zeros(numel(dz)+1,1);downwardFace=upwardFace;bottomOut=0;topOut=0;
        for k=1:numel(segments)
            segment=segments{k};
            [state.mobile,parts{k}]=soil.nitrogen.tracer.replayTransportLedger( ...
                state.mobile,segment.transport,dz);
            upwardFace=upwardFace+parts{k}.upward_face_amount_g_cm2;
            downwardFace=downwardFace+parts{k}.downward_face_amount_g_cm2;
            bottomOut=bottomOut+parts{k}.bottom_outflow_g_cm2;
            topOut=topOut+parts{k}.top_outflow_g_cm2;
            for step=1:numel(segment.exchange)
                e=segment.exchange{step};
                fm=fraction(state.mobile,e.before.mobile);
                fa=fraction(state.attached,e.before.attached);
                a=e.attach.*fm;s=e.strain.*fm;d=e.detach.*fa;
                attachLabel=attachLabel+a;detachLabel=detachLabel+d;
                strainLabel=strainLabel+s;
                % Multiply the accepted remaining total, rather than
                % subtracting almost equal labelled amounts at depletion.
                state.mobile=(e.before.mobile-e.attach-e.strain).*fm+d;
                state.attached=(e.before.attached-e.detach).*fa+a;
                state.retained=state.retained+s;
                for j=1:3
                    p=phases{j};upper=e.after.(p);tol=512*eps(max([upper(:);realmin]));
                    assert(all(state.(p)>=-tol&state.(p)<=upper+tol), ...
                        'tracer:particleExchangeBound','Label exchange violates bounds.');
                end
            end
        end
        aggregate=zeros(size(dz));
        for j=1:3
            field=[group '_' phases{j} '_N'];value=state.(phases{j});
            upper=totalAfter.(field);tol=512*eps(max([upper(:);realmin]));
            assert(all(value>=-tol&value<=upper+tol), ...
                'tracer:particleEndpointBound','Particle label exceeds final total.');
            tagNext.(field)=value;aggregate=aggregate+value;
        end
        tagNext.([group '_N'])=aggregate;
        final=sum(aggregate.*dz);
        assert(abs(final-initial+bottomOut+topOut)<=4096*eps(max([initial;final;realmin])), ...
            'tracer:particleClosure','Particle transport/exchange does not close.');
        diagnostic.(group)=struct('transport',{parts}, ...
            'storage_before_g_cm2',initial,'storage_after_g_cm2',final, ...
            'closure_residual_g_cm2',final-initial+bottomOut+topOut, ...
            'bottom_outflow_g_cm2',bottomOut, ...
            'top_outflow_g_cm2',topOut, ...
            'attachment_g_cm3',attachLabel,'detachment_g_cm3',detachLabel, ...
            'straining_g_cm3',strainLabel, ...
            'upward_face_amount_g_cm2',upwardFace, ...
            'downward_face_amount_g_cm2',downwardFace, ...
            'net_upward_face_amount_g_cm2',upwardFace-downwardFace, ...
            'grid_order','bottom_up');
    end
end

function f=fraction(tag,total)
    tol=512*eps(max([total(:);realmin]));
    assert(all(tag>=-tol&tag<=total+tol),'tracer:particleDonorBound', ...
        'Particle label outside total donor.');
    f=zeros(size(total));active=total>0;f(active)=tag(active)./total(active);
    assert(all(f>=-512*eps&f<=1+512*eps),'tracer:particleFraction', ...
        'Invalid particle donor fraction.');
    f=min(max(f,0),1);
end
