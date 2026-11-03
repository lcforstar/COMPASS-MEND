function [tagNext,diagnostic]=transportDonFrozen(tag,~,~,~,~,dz,parameters,totalDiagnostic,totalAfter,tagSorbed,totalSorbed)
%TRANSPORTDONFROZEN Attribute accepted total-DON exchanges without feedback.
% Inputs/outputs: bottom-up g N/cm3 bulk; dz: cm. Accepted transport supplies
% chronological total-DON (equilibrium) or free-DON (dynamic) exchanges.
% Editor: Wenjie Liu; Date: 2026-09-21.
% Change: 2026-09-22 Wenjie Liu -- bounded implicit fractions replace unsafe
% explicit label replay. The accepted total calculation is unchanged.
% Change: 2026-10-28 Wenjie Liu -- dynamic transport replays only the free
% donor ledger; the immobile sorbed label remains within total DON.
    dynamic=strcmp(parameters.sorption_model,'mend_dynamic');
    assert(dynamic||strcmp(parameters.sorption_model,'linear_equilibrium'), ...
        'soil:nitrogen:tracer:donTopology','Unsupported DON topology.');
    if dynamic
        assert(nargin==11&&isequal(size(tag),size(tagSorbed),size(totalSorbed)) ...
            &&all(isfinite([tagSorbed(:);totalSorbed(:)])) ...
            &&all(tagSorbed>=0&tagSorbed<=tag) ...
            &&all(tagSorbed<=totalSorbed)&&all(totalSorbed<=totalAfter), ...
            'soil:nitrogen:tracer:dynamicSchema', ...
            'Dynamic DON requires explicit bounded total and sorbed labels.');
        tag=tag-tagSorbed;
    end
    [tagNext,diagnostic]=soil.nitrogen.tracer.replayTransportLedger(tag,totalDiagnostic,dz);
    if dynamic
        tol=512*eps(max([totalAfter(:);realmin]));
        assert(all(tagNext<=totalAfter-totalSorbed+tol), ...
            'soil:nitrogen:tracer:donFreeBound','Free DON label exceeds free total.');
        tagNext=tagNext+tagSorbed;
        fixed=sum(tagSorbed(:).*dz(:));
        diagnostic.storage_before_g_cm2=diagnostic.storage_before_g_cm2+fixed;
        diagnostic.storage_after_g_cm2=diagnostic.storage_after_g_cm2+fixed;
    end
    tol=512*eps(max([totalAfter(:);realmin]));
    assert(all(tagNext>=0&tagNext<=totalAfter(:)+tol), ...
        'soil:nitrogen:tracer:donBound','DON label outside total endpoint.');
end
