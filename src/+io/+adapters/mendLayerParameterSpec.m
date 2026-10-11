function spec=mendLayerParameterSpec()
%MENDLAYERPARAMETERSPEC Canonical layer-resolved MEND carbon kinetic inputs.
% Public soil.csv values retain MEND's hourly clock and dry-soil mass basis.
% The four-column group is atomic: soil owns all four, or legacy [microbe]
% scalars own all four. POM2/MOM sharing and fixed physiology are unchanged.
% Editor: Wenjie Liu; Date: 2026-09-29.
    spec=struct('names',{{'mend_vd_pom_h','mend_ks_pom1_mg_g', ...
        'mend_vg_h','mend_ks_dom_mg_g'}}, ...
        'units',{{'h^-1','mg C/g dry soil','h^-1','mg C/g dry soil'}});
end
