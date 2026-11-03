function out=toCompass(mend,bulkDensityBottomUp)
%TOCOMPASS Convert accepted native MEND pools to COMPASS diagnostic views.
%   MEND authoritative pools remain in native [mg g-1] top-down units. The
%   returned aggregates are reconstructable diagnostics, never extra owners.
%   Editor: Jizixing Zhu; Date: 2026-09-12.
    assert(isstruct(mend)&&isfield(mend,'POM1')&&isfield(mend,'POM2') ...
        &&isfield(mend,'MOC')&&isfield(mend,'MB')&&isfield(mend,'BD'), ...
        'soil:mend:adapterState','Incomplete MEND state.');
    bd=bulkDensityBottomUp(:); %#ok<NASGU>
    out=struct();
    out.POC=mend.POM1+mend.POM2;
    if isfield(mend,'QOM'), out.MAOC=mend.MOC+mend.QOM; else, out.MAOC=mend.MOC; end
    out.MBC=mend.MB+mend.BD;
    if isfield(mend,'DOMT'), out.DOC=mend.DOMT; elseif isfield(mend,'DOC'), out.DOC=mend.DOC; end
    out.authoritative_fields={'POM1','POM2','MOC','MB','BD','ENZP1','ENZP2','EM'};
end
