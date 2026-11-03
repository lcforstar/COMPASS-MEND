function spec=normalizeBoundary(boundary)
%NORMALIZEBOUNDARY Immutable DOC/DON column boundary policy from site.toml.
% dom_bottom_kind: closed (legacy default), or advective_outflow.
% The latter explicitly requires dom_bottom_water_phase=liquid_only and
% water_bottom_kind=free_drainage. This is a zero-vapor bottom closure, not
% a general decomposition of total water flux. Top external DOM flux and
% external dispersive flux are zero. No external DOM concentration is assumed.
% Editor: Wenjie Liu; Date: 2026-09-29.
    assert(isstruct(boundary)&&isscalar(boundary),'soil:dom:boundaryConfig', ...
        'DOM boundary input must be a scalar struct.');
    % Change: 2026-09-29 Wenjie Liu -- numbered DOC/DON boundaries, organic only.
    keys=fieldnames(boundary);
    if any(startsWith(keys,'doc_')|startsWith(keys,'don_'))
        spec=soil.dom.normalizeColumnBoundary(boundary,{'doc','don'},false);return
    end
    names=fieldnames(boundary);selected=names(startsWith(names,'dom_'));
    assert(isempty(setdiff(selected,{'dom_bottom_kind','dom_bottom_water_phase'})), ...
        'soil:dom:boundaryConfig','Unknown [boundary] dom_* field.');
    kind='closed';phase='none';
    if isfield(boundary,'dom_bottom_kind'),kind=char(string(boundary.dom_bottom_kind));end
    if isfield(boundary,'dom_bottom_water_phase'),phase=char(string(boundary.dom_bottom_water_phase));end
    assert(ismember(kind,{'closed','advective_outflow'}), ...
        'soil:dom:boundaryConfig','dom_bottom_kind must be closed or advective_outflow.');
    if strcmp(kind,'advective_outflow')
        assert(strcmp(phase,'liquid_only')&&isfield(boundary,'water_bottom_kind') ...
            &&strcmp(boundary.water_bottom_kind,'free_drainage'), ...
            'soil:dom:boundaryWaterPolicy', ...
            ['Open DOM requires free_drainage and explicit ', ...
             'dom_bottom_water_phase=liquid_only (zero bottom vapor exchange).']);
    else
        assert(ismember(phase,{'none','liquid_only'}),'soil:dom:boundaryWaterPolicy', ...
            'Unsupported DOM bottom water phase.');
        phase='none'; % inactive alternatives do not alter the closed signature
    end
    spec=struct('schema_version',1,'bottom_kind',kind,'water_phase',phase, ...
        'top_kind','closed','external_diffusive_flux','zero', ...
        'reverse_flow_policy','error_requires_external_concentration', ...
        'advected_concentration','aqueous_g_cm3_water', ...
        'ledger_version',1);
end
