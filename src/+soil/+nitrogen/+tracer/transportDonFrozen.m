function [tagNext,diagnostic]=transportDonFrozen(tag,~,~,~,~,dz,parameters,totalDiagnostic,totalAfter)
%TRANSPORTDONFROZEN Attribute accepted total-DON exchanges without feedback.
% Inputs/outputs: bottom-up g N/cm3 bulk; dz: cm. Linear-equilibrium total
% transport supplies chronological accepted endpoints and directed face masses.
% Editor: Wenjie Liu; Date: 2026-09-21.
% Change: 2026-09-22 Wenjie Liu -- bounded implicit fractions replace unsafe
% explicit label replay. The accepted total calculation is unchanged.
    assert(strcmp(parameters.sorption_model,'linear_equilibrium'), ...
        'soil:nitrogen:tracer:donTopology','Linear equilibrium is required.');
    [tagNext,diagnostic]=soil.nitrogen.tracer.replayTransportLedger(tag,totalDiagnostic,dz);
    tol=512*eps(max([totalAfter(:);realmin]));
    assert(all(tagNext>=0&tagNext<=totalAfter(:)+tol), ...
        'soil:nitrogen:tracer:donBound','DON label outside total endpoint.');
end
