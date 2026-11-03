function memory=emptyBoundaryExchange(carriers)
%EMPTYBOUNDARYEXCHANGE Additional organic cumulative routes [g element/cm2].
% Bottom outflow remains owned by existing bottom_export_cumulative; no copy.
% This memory is candidate-owned until the enclosing KT commits.
% Editor: Wenjie Liu; Date: 2026-09-29.
    memory=struct('schema_version',1);
    for k=1:numel(carriers)
        memory.(carriers{k})=struct('top_inflow_g_cm2',0,'top_outflow_g_cm2',0, ...
            'bottom_inflow_g_cm2',0);
    end
end
