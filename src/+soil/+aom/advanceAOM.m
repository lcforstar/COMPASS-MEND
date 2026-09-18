function candidate=advanceAOM(accepted,temperatureC,headCm,theta, ...
        dzCm,parameters,dt_s)
%ADVANCEAOM Advance litter, straw, and manure AOM with the CoLM rate kernel.
%
%   The temperature, matric-potential, and depth scalars are an algebraically
%   equivalent MATLAB transcription of CoLM
%   MOD_BGC_Soil_BiogeochemDecompCascadeBGC.F90. Surface AOM retains the
%   implicit-decomposer respiratory split: decomposed C becomes CO2-C plus
%   POC-C, excess decomposed organic N becomes TAN, and water-driven leaching
%   enters DOC/DON. Incorporated AOM uses
%   the same frozen rate/environment kernel but has no direct respiratory
%   branch: decomposed material enters POC/PON and leached material enters
%   DOC/DON. Its later CO2-C is owned solely by explicit MEMC/MEND microbes.
%   Organic N follows donor C:N and never accompanies CO2-C.
%
%   Element stores and outputs are bottom-up [g cm^-3 bulk soil]; surface
%   stores are [g cm^-2]. Flux diagnostics are step-integrated masses.
%
% Editor: Jizixing Zhu; Date: 2026-08-28.
% Change: 2026-08-30 Codex -- separate surface implicit respiration from the
%   incorporated AOM precursor path and prohibit direct soil-AOM CO2.
% Change: 2026-09-01 Codex -- close surface implicit decomposition against
%   the POC product C:N; excess decomposed organic N becomes TAN for the
%   COMPASS NH4/NH3 owner, never direct NH3.
% Change: 2026-09-08 Jizixing Zhu -- persist CoLM litter chemical components.

    validateInputs(accepted,temperatureC,headCm,theta,dzCm,dt_s);
    n=numel(dzCm); dz=dzCm(:);
    t=temperatureC(:); h=headCm(:); water=theta(:); %#ok<NASGU>
    [fT,fW]=colmScalars(t,h,parameters);
    dzTop=flip(dz);
    centers=flip(cumsum(dzTop)-0.5.*dzTop);
    fDepth=exp(-(centers./100)./parameters.depth_efolding_m);
    materials={'litter','straw','manure'};
    k=[parameters.k_litter_ref_s,parameters.k_straw_ref_s, ...
        parameters.k_manure_ref_s];
    state=accepted; poc=zeros(n,1); pon=zeros(n,1);
    doc=zeros(n,1); don=zeros(n,1); tan=zeros(n,1); co2=zeros(n,1);
    materialDiag=struct();
    for i=1:numel(materials)
        name=materials{i}; a=state.aom.(name);
        if strcmp(name,'litter')&&isfield(parameters,'litter_component_rates_s')
            assert(isfield(a,'components'),'advanceAOM:missingComponents', ...
                'Component litter requires a new initialization or compatible checkpoint.');
            [a,flux]=advanceComponents(a,fT,fW,fDepth,dz,parameters,dt_s);
        else
            [a,flux]=advanceOne(a,k(i),fT,fW,fDepth,dz,parameters,dt_s);
        end
        state.aom.(name)=a;
        poc=poc+flux.poc_C; pon=pon+flux.pon_N;
        doc=doc+flux.doc_C; don=don+flux.don_N;
        tan=tan+flux.tan_N; co2=co2+flux.co2_C;
        materialDiag.(name)=flux.diagnostic;
        for product={'poc_C','doc_C','co2_C','pon_N','don_N','tan_N'}
            key=product{1};
            state.cumulative=addCumulative(state.cumulative, ...
                [name '_to_' cumulativeSuffix(key)], ...
                flux.diagnostic.([key '_g_cm2']));
        end
    end
    state.PON=state.PON+pon;
    state.cumulative.co2_C_g_cm2=state.cumulative.co2_C_g_cm2+sum(co2.*dz);
    candidate=struct('state',state,'poc_C_g_cm3',poc, ...
        'pon_N_g_cm3',pon,'doc_C_g_cm3',doc,'don_N_g_cm3',don, ...
        'tan_N_g_cm3',tan,'tan_N_rate_ug_cm3_s',tan.*1e6./dt_s, ...
        'co2_C_g_cm3',co2,'co2_C_rate_ug_cm3_s',co2.*1e6./dt_s, ...
        'temperature_scalar',fT,'water_scalar',fW, ...
        'depth_scalar',fDepth,'materials',materialDiag, ...
        'surface_rate_equations_unchanged',true, ...
        'soil_rate_equations_unchanged',true, ...
        'coupling_partition_adapted',true);
    candidate.carbon_residual_g_cm2=sumResidual(materialDiag,'carbon');
    candidate.nitrogen_residual_g_cm2=sumResidual(materialDiag,'nitrogen');
    tolerance=1e-12*max(1,totalAomCarbon(accepted,dz));
    assert(abs(candidate.carbon_residual_g_cm2)<=tolerance ...
        &&abs(candidate.nitrogen_residual_g_cm2)<=tolerance, ...
        'advanceAOM:massClosure','AOM C/N step failed mass closure.');
end

function [a,flux]=advanceComponents(a,fT,fW,fDepth,dz,p,dt)
% CoLM per-component rates; existing COMPASS POC/DOM coupling unchanged.
% CWD fragments internally, without respiration or direct soil-pool input.
% Editor: Jizixing Zhu; Date: 2026-09-08.
    a=soil.litter.componentTotals(a,true); before=a;
    names={'metabolic','cellulose','lignin'};
    for j=1:3
        pj=p; pj.respiration_fraction=p.litter_component_respiration(j);
        [a.components.(names{j}),fj]=advanceOne(a.components.(names{j}), ...
            p.litter_component_rates_s(j),fT,fW,fDepth,dz,pj,dt);
        if j==1
            flux=fj;
        else
            for key={'poc_C','pon_N','doc_C','don_N','tan_N','co2_C'}
                flux.(key{1})=flux.(key{1})+fj.(key{1});
            end
            keys=fieldnames(flux.diagnostic);
            for q=1:numel(keys)
                flux.diagnostic.(keys{q})=flux.diagnostic.(keys{q})+fj.diagnostic.(keys{q});
            end
        end
    end
    fracElement=-expm1(-p.litter_component_rates_s(4).*fT.*fW.*fDepth.*dt);
    fracSurface=-expm1(-p.litter_component_rates_s(4).*fT(end).*fW(end).*dt);
    for key={'surface_C','surface_N','element_C','element_N'}
        field=key{1}; fraction=fracElement;
        if startsWith(field,'surface'), fraction=fracSurface; end
        moved=a.components.cwd.(field).*fraction;
        a.components.cwd.(field)=a.components.cwd.(field)-moved;
        a.components.cellulose.(field)=a.components.cellulose.(field)+p.litter_cwd_fractions(1).*moved;
        a.components.lignin.(field)=a.components.lignin.(field)+p.litter_cwd_fractions(2).*moved;
    end
    a=soil.litter.componentTotals(a,false);
    flux.diagnostic.carbon_residual_g_cm2=a.surface_C-before.surface_C ...
        +sum((a.element_C-before.element_C+flux.poc_C+flux.doc_C+flux.co2_C).*dz);
    flux.diagnostic.nitrogen_residual_g_cm2=a.surface_N-before.surface_N ...
        +sum((a.element_N-before.element_N+flux.pon_N+flux.don_N+flux.tan_N).*dz);
end

function [a,flux]=advanceOne(a,k,fT,fW,fDepth,dz,p,dt)
% Advance one material while retaining surface versus incorporated identity.
% Editor: Jizixing Zhu; Date: 2026-08-28.
    c=a.element_C(:); n=a.element_N(:);
    rate=k.*fT.*fW.*fDepth;
    decompC=c.*(1-exp(-rate.*dt));
    leachC=(c-decompC).*(1-exp(-p.doc_leach_rate_s.*fW.*dt));
    ratio=safeRatio(n,c);
    decompN=min(n,decompC.*ratio);
    remainingN=n-decompN;
    leachN=min(remainingN,leachC.*ratio);
    cNew=c-decompC-leachC; nNew=n-decompN-leachN;

    surfaceRate=k.*fT(end).*fW(end);
    surfaceDecompC=a.surface_C.*(1-exp(-surfaceRate.*dt));
    surfaceLeachC=(a.surface_C-surfaceDecompC) ...
        .*(1-exp(-p.doc_leach_rate_s.*fW(end).*dt));
    surfaceRatio=safeRatio(a.surface_N,a.surface_C);
    surfaceDecompN=min(a.surface_N,surfaceDecompC.*surfaceRatio);
    surfaceRemainingN=a.surface_N-surfaceDecompN;
    surfaceLeachN=min(surfaceRemainingN,surfaceLeachC.*surfaceRatio);
    a.surface_C=max(a.surface_C-surfaceDecompC-surfaceLeachC,0);
    a.surface_N=max(a.surface_N-surfaceDecompN-surfaceLeachN,0);
    a.element_C=max(cNew,0); a.element_N=max(nNew,0);

    % Incorporated AOM is a soil precursor: no direct CO2 branch. Surface AOM
    % uses the CoLM implicit-decomposer respiratory fraction and deposits its
    % non-respired particulate product plus soluble leachate in the top cell.
    poc=decompC;
    pon=decompN;
    doc=leachC;
    don=leachN;
    co2=zeros(size(decompC));
    co2(end)=p.respiration_fraction.*surfaceDecompC./dz(end);
    poc(end)=poc(end)+(1-p.respiration_fraction).*surfaceDecompC./dz(end);
    surfacePocC=(1-p.respiration_fraction).*surfaceDecompC;
    surfacePonN=min(surfaceDecompN,surfacePocC./p.surface_product_cn);
    surfaceTanN=surfaceDecompN-surfacePonN;
    pon(end)=pon(end)+surfacePonN./dz(end);
    tan=zeros(size(decompN));
    tan(end)=surfaceTanN./dz(end);
    doc(end)=doc(end)+surfaceLeachC./dz(end);
    don(end)=don(end)+surfaceLeachN./dz(end);
    beforeC=sum(c.*dz)+surfaceDecompC+surfaceLeachC+a.surface_C;
    afterC=sum(a.element_C.*dz)+a.surface_C;
    beforeN=sum(n.*dz)+surfaceDecompN+surfaceLeachN+a.surface_N;
    afterN=sum(a.element_N.*dz)+a.surface_N;
    carbonResidual=afterC+sum((co2+poc+doc).*dz)-beforeC;
    nitrogenResidual=afterN+sum((pon+don+tan).*dz)-beforeN;
    flux=struct('poc_C',poc,'pon_N',pon, ...
        'doc_C',doc,'don_N',don,'tan_N',tan, ...
        'co2_C',co2,'diagnostic',struct( ...
            'poc_C_g_cm2',sum(poc.*dz), ...
            'doc_C_g_cm2',sum(doc.*dz), ...
            'co2_C_g_cm2',sum(co2.*dz), ...
            'pon_N_g_cm2',sum(pon.*dz), ...
            'don_N_g_cm2',sum(don.*dz), ...
            'tan_N_g_cm2',sum(tan.*dz), ...
            'surface_tan_N_g_cm2',surfaceTanN, ...
            'surface_co2_C_g_cm2',p.respiration_fraction.*surfaceDecompC, ...
            'soil_co2_C_g_cm2',0, ...
            'carbon_residual_g_cm2',carbonResidual, ...
            'nitrogen_residual_g_cm2',nitrogenResidual));
end

function suffix=cumulativeSuffix(product)
% Stable cumulative-ledger suffix for one AOM product.
% Editor: Codex; Date: 2026-08-28.
    map=struct('poc_C','POC_C_g_cm2','doc_C','DOC_C_g_cm2', ...
        'co2_C','CO2_C_g_cm2','pon_N','PON_N_g_cm2', ...
        'don_N','DON_N_g_cm2','tan_N','TAN_N_g_cm2');
    suffix=map.(product);
end

function cumulative=addCumulative(cumulative,name,amount)
% Add a material-resolved output while accepting older checkpoints.
% Editor: Codex; Date: 2026-08-28.
    if amount==0&&~isfield(cumulative,name), return; end
    if ~isfield(cumulative,name), cumulative.(name)=0; end
    cumulative.(name)=cumulative.(name)+amount;
end

function [fT,fW]=colmScalars(tC,hCm,p)
% Exact CoLM Q10 and log-matric-potential scalar transcription.
% Editor: Jizixing Zhu; Date: 2026-08-28.
    fT=zeros(size(tC)); warm=tC>=0;
    fT(warm)=p.q10.^((tC(warm)-25)./10);
    fT(~warm)=p.q10.^(-2.5).*p.frozen_q10.^(tC(~warm)./10);
    catan=@(x)11.75+(29.7/pi).*atan(pi.*0.031.*(x-15.4));
    normalization=(catan(15)./catan(30))./(p.q10.^(-1));
    fT=fT.*normalization;
    psi=min(hCm(:).*10,p.smpmax_mm);
    fW=repmat(0.001,size(psi)); active=psi>p.smpmin_mm;
    fW(active)=log(p.smpmin_mm./psi(active)) ...
        ./log(p.smpmin_mm./p.smpmax_mm);
    fW=min(max(fW,0.001),1);
end

function ratio=safeRatio(n,c)
% Donor N:C with an exact zero-carbon limit.
% Editor: Jizixing Zhu; Date: 2026-08-28.
    ratio=zeros(size(c)); active=c>0; ratio(active)=n(active)./c(active);
end

function value=sumResidual(d,name)
% Sum material ledger residuals.
% Editor: Jizixing Zhu; Date: 2026-08-28.
    fields=fieldnames(d); value=0;
    key=[name '_residual_g_cm2'];
    for i=1:numel(fields), value=value+d.(fields{i}).(key); end
end

function value=totalAomCarbon(state,dzBottom)
% Total AOM carbon column for tolerance scaling.
% Editor: Jizixing Zhu; Date: 2026-08-28.
    names=fieldnames(state.aom); value=0;
    for i=1:numel(names)
        a=state.aom.(names{i});
        value=value+a.surface_C+sum(a.element_C(:).*dzBottom(:));
    end
end

function validateInputs(state,t,h,theta,dz,dt)
% Validate shapes and finite environmental drivers.
% Editor: Jizixing Zhu; Date: 2026-08-28.
    n=numel(dz);
    assert(isfield(state,'aom')&&state.schema_version>=2, ...
        'advanceAOM:state','Native AOM state is required.');
    assert(numel(t)==n&&numel(h)==n&&numel(theta)==n ...
        &&all(isfinite([t(:);h(:);theta(:);dz(:)])) ...
        &&all(theta(:)>=0)&&all(dz(:)>0)&&isscalar(dt)&&isfinite(dt)&&dt>0, ...
        'advanceAOM:drivers','AOM environmental drivers are invalid.');
end
