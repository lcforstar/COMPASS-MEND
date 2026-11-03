function report=exportOrganicInitialPoolsCsv(sourceStateOrMat,domain,bulkDensity,targetCsvFile)
%EXPORTORGANICINITIALPOOLSCSV Export a state as resolved SOC-relative CSV priors.
% Source is accepted Y, built/model, or a MAT containing initialY (offline
% conversion only). domain supplies bottom-up z_node/dz, bulkDensity is
% bottom-up NLx1 [g/cm3] or struct.layer. Empty domain/rho use MAT metadata.
% C mg/g top-down -> g/cm3 bottom-up. SOC excludes AOM, includes Tot, total
% DOM, IOM. Dynamic DOC/QOC and DON/QON are exported as independent phases;
% explicit legacy equilibrium exports its historical total-DOM prior.
% Independent C pools sum to one; N/SOC and parent-relative C/N phases
% preserve stoichiometry. AOM totals/components exported separately in their
% native density units; cumulative fluxes are never exported as initial mass.
% Output CSV uses %.17g, ascending positive centers. Throws on invalid mass,
% aliases, phase closure, unsupported NEC, or missing spatial information.
% Returns provenance and closure diagnostics. No model state is modified.
% Editor: Wenjie Liu; Date: 2026-10-17.
% Change: 2026-10-23 Wenjie Liu -- preserve dynamic D/Q and DON/QON in CSV;
%   total DOM must not be mislabeled as D or counted in both DOC and QOC.
    source=sourceStateOrMat;sourceFile='provided accepted state';
    if ischar(source)||isstring(source)
        sourceFile=char(source);source=load(sourceFile);
    end
    if isfield(source,'built'),source=source.built;end
    if isempty(domain)&&isfield(source,'domain'),domain=source.domain;end
    if isempty(bulkDensity)&&isfield(source,'parameters'),bulkDensity=source.parameters.soil.bulk_density;end
    if isstruct(bulkDensity),bulkDensity=bulkDensity.layer;end
    if isfield(source,'initialY'),y=source.initialY;
    elseif isfield(source,'Y'),y=source.Y;else,y=source;end
    rho=bulkDensity(:);n=numel(rho);spec=io.adapters.organicPoolSpec();
    assert(n==domain.NL&&all(isfinite(rho))&&all(rho>0),'organicExport:domain','Valid matching domain and bulk density required.');
    m=y.soil_microbes;o=y.soil_organic;d=y.soil_dom;
    assert(all(o.NEC_C==0)&&all(o.NEC_N==0),'organicExport:nec', ...
        'Nonzero NEC is not represented by this twelve-pool format.');
    C=struct();
    for key=spec.carbon(1:10),k=key{1};C.(k)=flip(m.(k)(:)).*rho/1000;end
    C.IOM=o.IOM_C(:);C.DOC=d.doc_total_g_cm3(:);
    soc=flip(m.Tot(:)).*rho/1000+C.DOC+C.IOM;
    dynamic=isfield(d,'sorption_model')&&strcmp(d.sorption_model,'mend_dynamic');
    if dynamic
        assert(d.schema_version>=2&&all(C.QOC==0)&&all(o.QON==0), ...
            'organicExport:duplicateQ','Dynamic QOC/QON must be owned only by soil_dom.');
        for kind={'doc','don'}
            k=kind{1};qkey=[k '_sorbed_g_cm3'];totalDOM=d.([k '_total_g_cm3']);
            assert(isfield(d,qkey)&&isequal(size(d.(qkey)),size(totalDOM)) ...
                &&all(isfinite(d.(qkey)))&&all(d.(qkey)>=0)&&all(d.(qkey)<=totalDOM), ...
                'organicExport:dynamicInventory','Invalid accepted dynamic %s Q.',k);
        end
        C.QOC=d.doc_sorbed_g_cm3(:);C.DOC=C.DOC-C.QOC;
        domDefinition='DOC=D; QOC=Q; DON=DN; QON=QN (independent dynamic inventories)';
    else
        assert(~isfield(d,'doc_sorbed_g_cm3')&&~isfield(d,'don_sorbed_g_cm3'), ...
            'organicExport:phaseMode','Sorbed owners require an explicit mend_dynamic state tag.');
        domDefinition='legacy equilibrium: DOC=D+Q; DON=DN+QN (total priors)';
    end
    total=zeros(n,1);for key=spec.carbon,total=total+C.(key{1});end
    assert(all(soc>0)&&all(abs(total-soc)<=1e-13.*soc),'organicExport:closure','SOC alias must match independent pools.');
    z=.5*(domain.z_node(1:end-1)+domain.z_node(2:end));
    names={'depth_cm'};values=flip(z(:));
    for key=spec.carbon,append(['frac_' key{1}],flip(C.(key{1})./soc));end
    % Make the largest positive pool absorb only the floating-point summation
    % remainder, never normalize physical pools or silently repair bad states.
    carbonColumns=find(startsWith(names,'frac_'));
    for row=1:n
        [~,largest]=max(values(row,carbonColumns));idx=carbonColumns(largest);
        other=carbonColumns(carbonColumns~=idx);values(row,idx)=1-sum(values(row,other),2);
    end
    for key=spec.nitrogen
        k=key{1};if strcmp(k,'DON'),v=d.don_total_g_cm3(:);else,v=o.(k)(:);end
        if dynamic
            if strcmp(k,'DON'),v=v-d.don_sorbed_g_cm3(:);
            elseif strcmp(k,'QON'),v=d.don_sorbed_g_cm3(:);end
        end
        append(['nfrac_' k],flip(v./soc));
    end
    for group={'MB','BD','ENZN'}
        g=group{1};parentC=m.(g)(:);parentN=o.([g '_N'])(:);
        for phase=spec.phases
            k=[g '_' phase{1}];
            append(['frac_' k],safeRatio(m.(k)(:),parentC));
            append(['frac_' k '_N'],flip(safeRatio(o.([k '_N'])(:),parentN)));
        end
    end
    for material={'litter','straw','manure'}
        g=material{1};a=o.aom.(g);addAom(a,['aom_' g]);
        if isfield(a,'components')
            for c=fieldnames(a.components)',addAom(a.components.(c{1}),['aom_' g '_' c{1}]);end
        end
    end
    assert(all(isfinite(values),'all')&&all(values>=0,'all'),'organicExport:mass','Export values must be finite and nonnegative.');
    fid=fopen(targetCsvFile,'w','n','UTF-8');assert(fid>=0,'organicExport:file','Cannot create target CSV.');
    cleanup=onCleanup(@()fclose(fid));
    fprintf(fid,'# Resolved C/SOC and N/SOC priors; phase fractions are relative to parent pools.\n');
    fprintf(fid,'# AOM is outside SOC; surface g/cm2 is present in the first row only.\n');
    fprintf(fid,'# Editor: Wenjie Liu; Change: 2026-10-17 export complete accepted organic t0 as CSV.\n');
    fprintf(fid,'# Change: 2026-10-23 Wenjie Liu -- %s.\n',domDefinition);
    fprintf(fid,'%s\n',strjoin(names,','));format=[repmat('%.17g,',1,numel(names)-1) '%.17g\n'];
    fprintf(fid,format,values.');clear cleanup
    report=struct('source',sourceFile,'target',char(targetCsvFile),'layers',n, ...
        'columns',numel(names),'SOC_definition','Tot_mg_g*rho/1000 + total_DOM_g_cm3 + IOM_g_cm3', ...
        'DOM_export_definition',domDefinition, ...
        'max_export_fraction_residual',max(abs(sum(values(:,2:13),2)-1)), ...
        'runtime_mat_dependency',false);
    function append(name,v),names{end+1}=name;values(:,end+1)=v(:);end
    function r=safeRatio(v,parent)
        assert(all(isfinite(v))&&all(v>=0)&&all(v(parent==0)==0),'organicExport:phase','Invalid phase inventory.');
        r=zeros(n,1);nonzero=parent>0;r(nonzero)=v(nonzero)./parent(nonzero);
        r(~nonzero)=1; % zero-parent closure is set below by the mobile phase
        if ~strcmp(phase{1},'mobile'),r(~nonzero)=0;end
    end
    function addAom(a,prefix)
        for element={'C','N'}
            e=element{1};append([prefix '_' e '_g_cm3'],flip(a.(['element_' e])(:)));
            v=zeros(n,1);v(1)=a.(['surface_' e]);append([prefix '_surface_' e '_g_cm2'],v);
        end
    end
end
