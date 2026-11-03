function [organic,microbes,dom,diagnostic]=stageExternalInput( ...
        organic,microbes,dom,carbonBottom,nitrogenBottom,fractions, ...
        parameters,bulkDensityBottom)
%STAGEEXTERNALINPUT Add one external C/N transfer to detached MEND owners.
%   CARBONBOTTOM/NITROGENBOTTOM are bottom-up layer columns [g C/N cm-3
%   bulk soil], BULKDENSITYBOTTOM is bottom-up [g dry soil cm-3]. FRACTIONS
%   gives [POM1 POM2 DOC] carbon fractions. PARAMETERS supplies CN_POM2 and
%   input_cn_lignin_to_labile for the conservative external N partition.
%   POM owners are top-down [mg C g-1 dry soil]; PON and total DOC/DON owners
%   are bottom-up [g element cm-3 bulk soil]. This boundary performs the only
%   unit/grid conversion. DIAGNOSTIC contains the six partition increments
%   in bottom-up [g element cm-3 bulk soil], plus the supplied totals.
%   Existing AOM and other physical owners are unchanged. POC, PON and Tot
%   are refreshed as diagnostic sums; Tot excludes the separately owned DOM.
%   No transaction is committed, no reaction is run, and no DIN is borrowed.
%   Editor: Wenjie Liu; Date: 2026-10-10.
    assert(isnumeric(carbonBottom)&&iscolumn(carbonBottom) ...
        &&isnumeric(nitrogenBottom)&&iscolumn(nitrogenBottom) ...
        &&isequal(size(carbonBottom),size(nitrogenBottom)) ...
        &&isnumeric(bulkDensityBottom)&&isreal(bulkDensityBottom) ...
        &&iscolumn(bulkDensityBottom) ...
        &&isequal(size(carbonBottom),size(bulkDensityBottom)) ...
        &&all(isfinite(bulkDensityBottom))&&all(bulkDensityBottom>0), ...
        'soil:mend:externalInputGrid', ...
        'External C/N and bulk density must be aligned bottom-up columns.');
    assert(isstruct(parameters)&&isfield(parameters,'CN_POM2') ...
        &&isfield(parameters,'input_cn_lignin_to_labile'), ...
        'soil:mend:externalInputParameters', ...
        'External MEND C:N parameters are required.');
    diagnostic=soil.mend.partitionExternalInput(carbonBottom,nitrogenBottom, ...
        fractions,parameters.CN_POM2,parameters.input_cn_lignin_to_labile);
    count=numel(carbonBottom);
    validateOwner(organic,{'PON1','PON2'},count);
    validateOwner(microbes, ...
        {'POM1','POM2','MOC','MB','BD','ENZP1','ENZP2','EM'},count);
    validateOwner(dom,{'doc_total_g_cm3','don_total_g_cm3'},count);
    scale=1000./flip(bulkDensityBottom);
    microbes.POM1=microbes.POM1+flip(diagnostic.C_POM1).*scale;
    microbes.POM2=microbes.POM2+flip(diagnostic.C_POM2).*scale;
    organic.PON1=organic.PON1+diagnostic.N_PON1;
    organic.PON2=organic.PON2+diagnostic.N_PON2;
    dom.doc_total_g_cm3=dom.doc_total_g_cm3+diagnostic.C_DOC;
    dom.don_total_g_cm3=dom.don_total_g_cm3+diagnostic.N_DON;
    microbes.POC=microbes.POM1+microbes.POM2;
    microbes.Tot=microbes.POC+microbes.MOC+microbes.MB+microbes.BD ...
        +microbes.ENZP1+microbes.ENZP2+microbes.EM;
    if isfield(microbes,'ENZN')
        validateOwner(microbes,{'ENZN'},count);
        microbes.Tot=microbes.Tot+microbes.ENZN;
    end
    organic.PON=organic.PON1+organic.PON2;
    diagnostic.total_C_g_cm3=carbonBottom;
    diagnostic.total_N_g_cm3=nitrogenBottom;
    diagnostic.grid_order='bottom_up';
    diagnostic.units='g_element_cm3_bulk_soil';
end

function validateOwner(owner,names,count)
% Require the active source-MEND column owners without changing their data.
% Editor: Wenjie Liu; Date: 2026-10-10.
    for j=1:numel(names)
        name=names{j};
        assert(isstruct(owner)&&isfield(owner,name) ...
            &&isnumeric(owner.(name))&&isreal(owner.(name)) ...
            &&iscolumn(owner.(name))&&numel(owner.(name))==count ...
            &&all(isfinite(owner.(name))), ...
            'soil:mend:externalInputOwner', ...
            'External MEND input requires aligned finite owner %s.',name);
    end
end
