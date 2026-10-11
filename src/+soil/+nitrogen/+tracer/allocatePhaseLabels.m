function [tagOut,diagnostic]=allocatePhaseLabels(totalIn,tagIn,totalOut,tagTarget)
%ALLOCATEPHASELABELS Bounded phase attribution of an aggregate reaction label.
% Rows are soil layers, columns mobile/attached/retained; all inputs use the
% same N concentration unit. Total owner scales phase inventories while
% gross reaction accounting supplies the aggregate labelled target.
% First carry phase fractions with the total resize. A decrease of labelled
% aggregate then removes labels proportionally; an increase occupies each
% phase's unlabelled capacity proportionally. Complement labels receive the
% complementary allocation. This is an attribution remap, not physical
% interphase mass transfer. It preserves the prescribed labelled aggregate
% without allowing any phase to exceed its immutable total endpoint.
% Editor: Wenjie Liu; Date: 2026-09-22.
    current=sum(totalIn,2);target=sum(totalOut,2);tagTarget=tagTarget(:);
    assert(isequal(size(totalIn),size(tagIn),size(totalOut)) ...
        &&size(totalIn,1)==numel(tagTarget), ...
        'tracer:phaseAllocationShape','Phase allocation shapes differ.');
    assert(all(isfinite([totalIn(:);tagIn(:);totalOut(:);tagTarget])) ...
        &&all([totalIn(:);tagIn(:);totalOut(:);tagTarget]>=0), ...
        'tracer:phaseAllocationDomain','Invalid phase allocation input.');
    tol=512*eps(max([totalIn(:);totalOut(:);realmin]));
    assert(all(tagIn<=totalIn+tol,'all')&&all(tagTarget<=target+3*tol), ...
        'tracer:phaseAllocationBound','Aggregate or initial phase label exceeds total.');
    carried=zeros(size(tagIn));active=current>0;
    carried(active,:)=tagIn(active,:).*(target(active)./current(active));
    tagOut=zeros(size(tagIn));
    for i=1:numel(target)
        base=carried(i,:);baseMass=sum(base);wanted=tagTarget(i);
        if wanted==0,continue;end
        % Summing resized phases can land one ULP below the aggregate
        % target. The bound check above permits only roundoff here, and
        % the residual below records/audits that endpoint discrepancy.
        if wanted>=target(i)
            tagOut(i,:)=totalOut(i,:);
        elseif wanted<=baseMass
            tagOut(i,:)=base*(wanted/baseMass);
        else
            room=totalOut(i,:)-base;
            assert(all(room>=-tol),'tracer:phaseAllocationCapacity','Negative phase capacity.');
            % Only endpoint arithmetic roundoff can make a capacity negative.
            room=max(room,0);roomMass=sum(room);
            assert(roomMass>0&&wanted-baseMass<=roomMass+3*tol, ...
                'tracer:phaseAllocationCapacity','Insufficient unlabelled phase capacity.');
            tagOut(i,:)=base+room*((wanted-baseMass)/roomMass);
        end
    end
    residual=sum(tagOut,2)-tagTarget;
    assert(all(tagOut>=0&tagOut<=totalOut+tol,'all') ...
        &&all(abs(residual)<=2048*eps(max([tagTarget;realmin]))), ...
        'tracer:phaseAllocationClosure','Bounded phase attribution did not close.');
    diagnostic=struct('method','total-resize carry then complementary capacity allocation', ...
        'aggregate_residual',residual,'carried_label',carried,'label_change',tagOut-carried);
end
