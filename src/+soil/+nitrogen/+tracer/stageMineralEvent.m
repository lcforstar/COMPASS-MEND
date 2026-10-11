function [state,input_g_cm2]=stageMineralEvent(accepted,fertilizer,totalTrial)
%STAGEMINERALEVENT Add the exact accepted mineral-fertilizer nodal label.
%
% This function is detached: it changes only its returned tracer candidate.
% Event flags and the candidate are committed together after outer acceptance.
% Editor: Wenjie Liu; Date: 2026-09-21.
    state=accepted;input_g_cm2=0;
    if ~fertilizer.applied
        return
    end
    assert(isfield(fertilizer,'tracer_increment') ...
        &&isstruct(fertilizer.tracer_increment), ...
        'soil:nitrogen:tracer:missingEventIncrement', ...
        'Applied fertilizer lacks exact nodal source-label increments.');
    delta=fertilizer.tracer_increment;
    mapping={'nh4','nh4_liquid';'no3','no3_liquid'; ...
        'urea','urea';'nh4_granule','nh4_granule'};
    for i=1:size(mapping,1)
        owner=mapping{i,1};increment=mapping{i,2};
        value=delta.(increment)(:);
        assert(isequal(size(value),size(state.(owner))) ...
            &&all(isfinite(value))&&all(value>=0), ...
            'soil:nitrogen:tracer:eventIncrement', ...
            'Invalid fertilizer source-label increment for %s.',owner);
        state.(owner)=state.(owner)+value;
        total=totalTrial.(owner)(:);
        tolerance=2e-12.*max(1,abs(total));
        assert(all(state.(owner)<=total+tolerance), ...
            'soil:nitrogen:tracer:eventBound', ...
            'Staged fertilizer label exceeds total %s.',owner);
    end
    input_g_cm2=1e-6.*(fertilizer.nh4_added+fertilizer.no3_added ...
        +fertilizer.urea_added);
end
