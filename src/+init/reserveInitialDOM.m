function built=reserveInitialDOM(built,reserveDon)
% Editor: Wenjie Liu; Date: 2026-10-17.
%RESERVEINITIALDOM Reserve new total-mg/kg DOM within native C/N inventory.
% Editor: Jizixing Zhu; Date: 2026-09-09.
% MOC C [g/cm3, top-down] and MON N [g/cm3, bottom-up] are donors.
% This changes initial partition only, not source reaction equations.
% Change: 2026-10-12 Wenjie Liu -- the historical *_total_mgkg_basis flag
% also covers explicit free-D input: reserve its reconstructed D+Q, not D.
% Change: 2026-10-12 Wenjie Liu -- optional reserveDon=false funds deferred
% free DOC after final water closure without reserving the existing DON twice.
    % Change: 2026-10-17 Wenjie Liu -- detailed CSV DOM is already partitioned.
    if nargin<2,reserveDon=true;end
    if ~isfield(built.Y,'soil_organic'),return;end
    if isfield(built.parameters.soil_organic,'initial_dom'),return;end
    d=built.Y.soil_dom;
    if d.initial_doc_total_mgkg_basis
        c=flip(d.doc_total_g_cm3);
        old=built.parameters.soil_organic.initial_carbon_top_down.MOC;
        assert(all(c<=old),'reserveInitialDOM:carbonBudget', ...
            'Initial DOC exceeds available MOC; revise initial partition.');
        built.parameters.soil_organic.initial_carbon_top_down.MOC=old-c;
        built.initialization.soil_organic.dom_C_reserved_g_cm2= ...
            sum(d.doc_total_g_cm3.*built.domain.dz(:));
    end
    if reserveDon&&d.initial_don_total_mgkg_basis
        old=built.Y.soil_organic.MON;n=d.don_total_g_cm3;
        assert(all(n<=old),'reserveInitialDOM:nitrogenBudget', ...
            'Initial DON exceeds available MON; revise initial partition.');
        built.Y.soil_organic.MON=old-n;
    end
end
