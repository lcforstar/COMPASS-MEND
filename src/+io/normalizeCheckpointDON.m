function [Y,audit]=normalizeCheckpointDON(Y)
%NORMALIZECHECKPOINTDON Enforce one conserved DON store at restart.
% Input/output Y: detached checkpoint value, never overwrites the MAT file.
% soil_dom.don_total_g_cm3 and soil_organic.QON are NLx1 bottom-up [g N/cm3].
% Legacy unmarked zero QON is safe. Nonzero QON requires an explicit audited
% soil_organic.qon_ownership: 'independent_legacy' or 'included_in_dom_total'.
% No inference from numerical equality or schema number is permitted.
% Independent legacy mass is added once; an included diagnostic is not added.
% Compatibility QON is then zero; equilibrium diagnostics are rebuilt by DOM.
% Ambiguous, negative, nonfinite or incompatible arrays reject before resume.
% Policy: user-approved single DON-total owner; not a MEND rate equation.
% Editor: Jizixing Zhu; Date: 2026-09-12.
    audit=struct('applied',false,'reason','no_paired_dom_state');
    if ~isfield(Y,'soil_dom')||~isfield(Y,'soil_organic'),return;end
    d=Y.soil_dom;o=Y.soil_organic;
    if isempty(d)&&isempty(o),return;end
    assert(isstruct(d)&&isstruct(o)&&isfield(d,'don_total_g_cm3')&& ...
        isfield(o,'QON'),'normalizeCheckpointDON:missingState','Missing paired DON/QON state.');
    total=d.don_total_g_cm3;q=o.QON;
    assert(iscolumn(total)&&~isempty(total)&&isequal(size(total),size(q))&& ...
        all(isfinite([total;q]))&&all([total;q]>=0), ...
        'normalizeCheckpointDON:invalidState','DON/QON must be matching nonnegative columns.');
    ownership='unmarked';
    if isfield(o,'qon_ownership'),ownership=o.qon_ownership;end
    assert(ischar(ownership)&&ismember(ownership, ...
        {'unmarked','independent_legacy','included_in_dom_total','dom_total_sole_owner'}), ...
        'normalizeCheckpointDON:unknownOwnership','Unknown QON provenance.');
    added=zeros(size(q));
    switch ownership
        case 'independent_legacy'
            added=q;
        case 'included_in_dom_total'
            assert(all(q<=total),'normalizeCheckpointDON:phaseExceedsTotal', ...
                'An included QON phase cannot exceed total DON.');
        otherwise
            assert(all(q==0),'normalizeCheckpointDON:ambiguousQON', ...
                'Nonzero QON has no verified owner. Audit provenance before migration.');
    end
    after=total+added;
    assert(all(isfinite(after)),'normalizeCheckpointDON:invalidState','DON mass overflow.');
    Y.soil_dom.don_total_g_cm3=after;
    Y.soil_organic.QON=zeros(size(q));
    Y.soil_organic.qon_ownership='dom_total_sole_owner';
    audit=struct('applied',true,'reason',ownership,'legacy_QON_g_cm3',q, ...
        'added_to_DOM_g_cm3',added,'DON_before_g_cm3',total, ...
        'DON_after_g_cm3',after,'residual_g_cm3',after-total-added);
end
