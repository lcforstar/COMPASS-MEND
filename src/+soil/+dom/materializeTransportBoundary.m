function boundary=materializeTransportBoundary(spec,face,carrier)
%MATERIALIZETRANSPORTBOUNDARY Freeze one organic carrier's accepted tile.
% face QMB and top QL are cm/s upward; external values already g/cm3 water.
% Legacy bottom-only faces are passed through unchanged. No state mutation.
% Editor: Wenjie Liu; Date: 2026-09-29.
    boundary=face;
    if ~isfield(spec,'schema_version')||spec.schema_version~=2,return;end
    assert(strcmp(spec.bottom_kind,'configured')&&isfield(spec.carriers,carrier) ...
        &&isfield(face,'q_top_cm_s'),'soil:dom:missingBoundaryDriver', ...
        'Configured organic boundary requires matching accepted top and bottom water drivers.');
    boundary.kind='configured';
    boundary.top=spec.carriers.(carrier).top;
    boundary.bottom=spec.carriers.(carrier).bottom;
end
