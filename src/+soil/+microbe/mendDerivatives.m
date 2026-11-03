function [d,f] = mendDerivatives(state,p)
%MENDDERIVATIVES Assemble conservative MEND carbon derivatives per hour.
%
% Editor: Jizixing Zhu; Date: 2026-09-12.
% Change: 2026-09-14 Codex -- apply one oxygen factor to the full reaction
%   ledger, preserving stoichiometry while rebuilding a supply-limited candidate.
    f=soil.microbe.mendFluxes(state,p);
    if isfield(p,'reaction_scale_top_down')
        scale=p.reaction_scale_top_down(:);
        assert(isequal(size(scale),size(state.MB))&&all(isfinite(scale)) ...
            &&all(scale>=0&scale<=1),'soil:microbe:oxygenScale','Invalid oxygen reaction factor.');
        names=fieldnames(f);
        for j=1:numel(names)
            key=names{j};
            if strcmp(key,'substrate_saturation'),continue;end
            f.(key)=f.(key).*scale;
        end
    end
    d.POM1=f.MBA_to_POM1-f.POM1_decomp;
    d.POM2=f.MBA_to_POM2-f.POM2_decomp;
    d.MOC=f.POM1_to_MOM+f.POM2_to_MOM-f.MOM_to_DOM;
    d.DOC=f.POM1_to_DOM+f.POM2_to_DOM+f.MOM_to_DOM+f.MBA_to_DOM ...
        +f.ENZP1_to_DOM+f.ENZP2_to_DOM+f.ENZM_to_DOM-f.DOM_to_MBA;
    production=f.MBA_to_ENZP1+f.MBA_to_ENZP2+f.MBA_to_ENZM;
    d.MB=f.DOM_to_MBA-f.CO2_growth-f.CO2_maintenance_MBA ...
        -f.MBA_mortality-production-f.MBA_to_MBD+f.MBD_to_MBA;
    d.BD=-f.CO2_maintenance_MBD+f.MBA_to_MBD-f.MBD_to_MBA;
    d.ENZP1=f.MBA_to_ENZP1-f.ENZP1_to_DOM;
    d.ENZP2=f.MBA_to_ENZP2-f.ENZP2_to_DOM;
    d.EM=f.MBA_to_ENZM-f.ENZM_to_DOM;
    d.IC=f.CO2_growth+f.CO2_maintenance_MBA+f.CO2_maintenance_MBD;
    residual=d.POM1+d.POM2+d.MOC+d.DOC+d.MB+d.BD ...
        +d.ENZP1+d.ENZP2+d.EM+d.IC;
    assert(max(abs(residual))<=5e-13.*max(1,max(abs(d.IC))), ...
        'soil:microbe:mendDerivativeClosure', ...
        'MEND carbon derivative does not close.');
end
