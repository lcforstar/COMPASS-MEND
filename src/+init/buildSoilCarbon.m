function [state, parameters, report] = buildSoilCarbon(raw, built, bd)
% Editor: Wenjie Liu; Date: 2026-10-17.
%BUILDSOILCARBON Initialize COMPASS organic N and MEND carbon priors.
%
%   Total SOC is constrained by soil_initial_profile.csv. The existing
%   depth-resolved native fractions initialize POC, MOC, MB and IOM directly.
%   Carbon reaction pools are handed to MEND by
%   PARAMETERS.initial_carbon_*; this accepted state owns AOM and all paired
%   organic-N ledgers without duplicating MEMC carbon storage.
%
%   Pool concentrations use [g C or N cm^-3 bulk soil], bottom-up unless a
%   field explicitly ends in _top_down. Surface AOM stores use [g cm^-2].
%
% Editor: Jizixing Zhu; Date: 2026-08-28.
% Change: 2026-09-08 Jizixing Zhu -- optional CoLM component litter state.
% Change: 2026-09-12 Jizixing Zhu -- make MEND POM/enzyme N subpools
%   authoritative while retaining checked compatibility aggregates.
% Change: 2026-10-17 Wenjie Liu -- resolve cold or detailed CSV priors once;
%   retain full C/N phases and AOM without an external state import.
% Change: 2026-10-17 Wenjie Liu -- source-hour D/Q priors fund one DOM owner;
%   preserve layer C/N and let configured sorption determine phase partition.
% Change: 2026-10-23 Wenjie Liu -- restore independent source-MEND D/Q and
%   DON/QON priors for mend_dynamic; transfer Q ownership without equilibrating
%   or debiting MOC/MON. Explicit linear_equilibrium retains its total prior.

    assert(istable(raw.organic_initial) && ~isempty(raw.organic_initial), ...
        'buildSoilCarbon:missingFractions', ...
        'organic_initial_pools.csv is required for first native initialization.');
    assert(isfield(raw.config,'organic_n') && isstruct(raw.config.organic_n), ...
        'buildSoilCarbon:missingOrganicN', ...
        '[organic_n] initialization and target C:N parameters are required.');

    n = built.domain.NL;
    dzTop = flip(built.domain.dz(:));
    zplus = -cumsum(dzTop);
    centersTop = -0.5.*([0;zplus(1:end-1)]+zplus);
    bdTop = flip(bd.layer(:));
    soc = initialSoc(raw.initial,centersTop,bdTop);
    prior=init.mapOrganicInitialPools(raw.organic_initial,centersTop,soc, ...
        raw.config.bioparticle.mobile_fraction_initial);
    if prior.detailed && strcmp(string(raw.config.microbe.reaction_method),"source_hourly")
        mode="mend_dynamic";
        if isfield(raw.config,'dom')&&isfield(raw.config.dom,'sorption_model')
            mode=string(raw.config.dom.sorption_model);
        end
        assert(isscalar(mode)&&ismember(mode,["mend_dynamic","linear_equilibrium"]), ...
            'buildSoilCarbon:sorption','Unknown dom.sorption_model.');
        if mode=="mend_dynamic"
            % Source initialization assigns D and Q independently (MEND
            % MOD_MEND.F90:143-151). Store Q before transferring ownership;
            % total-Q then recovers the supplied D without theta or Kd.
            prior.dom.doc_sorbed_g_cm3=flip(prior.carbon.QOC);
            prior.dom.don_sorbed_g_cm3=prior.nitrogen.QON;
            prior.dom.sorption_model='mend_dynamic';
        end
        % Both modes keep one DOM inventory. Dynamic also owns independent
        % Q; equilibrium keeps only total. Compatibility Q fields stay zero.
        prior.carbon.DOC=prior.carbon.DOC+prior.carbon.QOC;
        prior.nitrogen.DON=prior.nitrogen.DON+prior.nitrogen.QON;
        prior.carbon.QOC=zeros(n,1);prior.nitrogen.QON=zeros(n,1);
        prior.dom.doc_total_g_cm3=flip(prior.carbon.DOC);
        prior.dom.don_total_g_cm3=prior.nitrogen.DON;
    end
    if prior.detailed,fractions=[];
    else,fractions=initializationFractions(raw.organic_initial,centersTop,dzTop);end

    if ~prior.detailed
    pocC=soc.*fractions(:,1);
    mocC=soc.*fractions(:,2);
    mbC=soc.*fractions(:,3);
    iomC=soc.*fractions(:,4);
    qocC=zeros(n,1); epC=zeros(n,1); emC=zeros(n,1);

    cn = buildCnParameters(raw.config.organic_n,raw.config.microbe);
    zero=zeros(n,1);
    activeFraction=raw.config.microbe.active_fraction_initial;
    assert(isnumeric(activeFraction)&&isscalar(activeFraction) ...
        &&isfinite(activeFraction)&&activeFraction>0&&activeFraction<=1, ...
        'buildSoilCarbon:activeFraction', ...
        'microbe.active_fraction_initial must be in (0,1].');
    mobileFraction=raw.config.bioparticle.mobile_fraction_initial;
    activeN=flip(activeFraction.*mbC./cn.MB);
    dormantN=flip((1-activeFraction).*mbC./cn.MB);
    lci=raw.config.microbe.mend_lci_initial;
    pom1C=lci.*pocC; pom2C=(1-lci).*pocC;
    pon1=flip(pom1C./cn.POM1); pon2=flip(pom2C./cn.POM2);
    enzp1N=flip(lci.*epC./cn.ENZP1);
    enzp2N=flip((1-lci).*epC./cn.ENZP2);
    nitrogen=struct( ...
        'PON1',pon1,'PON2',pon2,'PON',pon1+pon2,'MON',flip(mocC./cn.MOC), ...
        'QON',flip(qocC./cn.QOC),'MB_N',activeN,'BD_N',dormantN, ...
        'MB_mobile_N',mobileFraction.*activeN,'MB_attached_N',(1-mobileFraction).*activeN, ...
        'MB_retained_N',zero,'BD_mobile_N',mobileFraction.*dormantN, ...
        'BD_attached_N',(1-mobileFraction).*dormantN,'BD_retained_N',zero, ...
        'ENZP1_N',enzp1N,'ENZP2_N',enzp2N,'EP_N',enzp1N+enzp2N, ...
        'EM_N',flip(emC./cn.ENZM),'NEC_C',zero,'NEC_N',zero,'IOM_C',flip(iomC));
    initialCarbon=struct('POM1',pom1C,'POM2',pom2C,'POC',pocC,'MOC',mocC, ...
        'QOC',qocC,'MB',mbC,'EP',epC,'EM',emC,'IOM',iomC);
    else
        cn=buildCnParameters(raw.config.organic_n,raw.config.microbe);
        zero=zeros(n,1);nitrogen=rmfield(prior.nitrogen,'DON');initialCarbon=prior.carbon;
    end
    state=struct( ...
        'schema_version',4,'model','COMPASS_MEND_CN', ...
        'grid_order','bottom_up','volume_units','g_per_cm3_bulk_soil', ...
        'surface_units','g_per_cm2_ground', ...
        'aom',struct( ...
            'litter',emptyAom(zero), ...
            'straw',emptyAom(zero), ...
            'manure',emptyAom(zero)), ...
        'cumulative',struct('input_C_g_cm2',0,'input_N_g_cm2',0, ...
            'input_litter_C_g_cm2',0,'input_litter_N_g_cm2',0, ...
            'input_straw_C_g_cm2',0,'input_straw_N_g_cm2',0, ...
            'input_manure_C_g_cm2',0,'input_manure_N_g_cm2',0, ...
            'co2_C_g_cm2',0,'unmet_N_g_cm2',0));

    % External root input has its own ledger, not the PlantCN litter ledger.
    % Change: 2026-10-30 Wenjie Liu -- initialize both root ledgers before
    %   the first event; configurations without root events retain old state.
    if isfield(built,'sources')&&isfield(built.sources,'organic_fertilization_schedule') ...
            &&any(strcmp(built.sources.organic_fertilization_schedule.material_type,'root_residue'))
        state.cumulative.input_root_residue_C_g_cm2=0;
        state.cumulative.input_root_residue_N_g_cm2=0;
    end

    parameters=struct( ...
        'grid_zplus_cm',zplus,'grid_order','top_down', ...
        'state_units','g_per_cm3_bulk_soil', ...
        'initial_carbon_top_down',initialCarbon, ...
        'target_cn',cn,'aom',buildAomParameters( ...
            raw.config.organic,cn.POC), ...
        'native_compass_state',true);
    for key=fieldnames(nitrogen)',state.(key{1})=nitrogen.(key{1});end
    if prior.detailed
        parameters.initialization_mode='detailed_soc_fractions';parameters.initial_dom=prior.dom;
        for material={'litter','straw','manure'},state.aom.(material{1})=prior.aom.(material{1});end
    end

    if isfield(parameters.aom,'litter_component_rates_s')
        state.schema_version=max(state.schema_version,4);
        for name={'metabolic','cellulose','lignin','cwd'}
            if ~isfield(state.aom.litter,'components')||~isfield(state.aom.litter.components,name{1})
                assert(~prior.detailed||all(state.aom.litter.element_C==0) ...
                    &&state.aom.litter.surface_C==0,'organicInitial:litterComponents', ...
                    'Nonzero component-mode litter requires explicit component priors.');
                state.aom.litter.components.(name{1})=emptyAom(zero);
            end
        end
    end
    initialC=sum(soc.*dzTop);
    if prior.detailed
        spec=io.adapters.organicPoolSpec();total=zeros(n,1);
        for key=spec.carbon,total=total+prior.carbon.(key{1});end
        partitionC=sum(total.*dzTop);
    else,partitionC=sum((pocC+mocC+qocC+mbC+epC+emC+iomC).*dzTop);end
    initialN=columnNitrogen(state,built.domain.dz(:));
    report=struct('initial_C_g_cm2',initialC, ...
        'state_C_g_cm2',partitionC, ...
        'carbon_residual_g_cm2',partitionC-initialC, ...
        'state_N_g_cm2',initialN, ...
        'initialization_mode','native_soc_priors', ...
        'total_soc_source','soil_initial_profile.csv:soc_total_mgkg', ...
        'fraction_source','organic_initial_pools.csv:native_pool_priors', ...
        'n_cells',n,'uses_compass_grid',true);
    if prior.detailed,report.initialization_mode='detailed_soc_fractions';end
    tolerance=1e-12*max(1,initialC);
    assert(abs(report.carbon_residual_g_cm2)<=tolerance, ...
        'buildSoilCarbon:carbonClosure', ...
        'Native SOC initialization residual exceeds tolerance.');
    validateState(state,n);
end

function p=buildAomParameters(cfg,surfaceProductCn)
% Build the frozen CoLM/CLM decomposition-scalar and AOM partition contract.
% Editor: Jizixing Zhu; Date: 2026-08-28.
% Change: 2026-09-02 Codex -- read site-global AOM kinetics from [organic]
%   and convert public d^-1 rates once at the configuration boundary.
    rateKeys={'k_litter_1d','k_straw_1d','k_manure_1d','doc_leach_1d'};
    rateFields={'k_litter_ref_s','k_straw_ref_s','k_manure_ref_s', ...
        'doc_leach_rate_s'};
    for i=1:numel(rateKeys)
        assert(isfield(cfg,rateKeys{i})&&isnumeric(cfg.(rateKeys{i})) ...
            &&isscalar(cfg.(rateKeys{i}))&&isfinite(cfg.(rateKeys{i})) ...
            &&cfg.(rateKeys{i})>=0,'buildSoilCarbon:badAomParameter', ...
            'organic.%s must be a finite nonnegative scalar.',rateKeys{i});
        p.(rateFields{i})=double(cfg.(rateKeys{i}))/86400;
    end
    names={'respiration_fraction','q10','frozen_q10', ...
        'smpmax_mm','smpmin_mm','depth_efolding_m'};
    for i=1:numel(names)
        assert(isfield(cfg,names{i})&&isnumeric(cfg.(names{i})) ...
            &&isscalar(cfg.(names{i}))&&isfinite(cfg.(names{i})), ...
            'buildSoilCarbon:badAomParameter', ...
            'organic.%s must be a finite scalar.',names{i});
        p.(names{i})=double(cfg.(names{i}));
    end
    assert(all([p.k_litter_ref_s,p.k_straw_ref_s,p.k_manure_ref_s, ...
        p.doc_leach_rate_s]>=0)&&p.respiration_fraction>=0 ...
        &&p.respiration_fraction<=1&&p.q10>0&&p.frozen_q10>0 ...
        &&p.smpmin_mm<p.smpmax_mm&&p.smpmax_mm<0 ...
        &&p.depth_efolding_m>0, ...
        'buildSoilCarbon:invalidAomParameters', ...
        'The CoLM AOM parameter contract is physically invalid.');
    assert(isscalar(surfaceProductCn)&&isfinite(surfaceProductCn) ...
        &&surfaceProductCn>0,'buildSoilCarbon:surfaceProductCn', ...
        'The surface AOM particulate product C:N must be positive.');
    p.surface_product_cn=double(surfaceProductCn);
    % Change: 2026-09-08 Jizixing Zhu -- validate source component priors.
    if isfield(cfg,'litter_component_rates_1d')
        rates=cfg.litter_component_rates_1d(:);
        rf=cfg.litter_component_respiration(:);
        leaf=cfg.litter_leaf_fractions(:); root=cfg.litter_root_fractions(:);
        cwd=cfg.litter_cwd_fractions(:);
        assert(numel(rates)==4&&numel(rf)==3&&numel(leaf)==3 ...
            &&numel(root)==3&&numel(cwd)==2 ...
            &&all(isfinite([rates;rf;leaf;root;cwd])) ...
            &&all(rates>=0)&&all(rf>=0&rf<=1) ...
            &&all([leaf;root;cwd]>=0) ...
            &&all(abs([sum(leaf),sum(root),sum(cwd)]-1)<1e-12), ...
            'buildSoilCarbon:litterComponents','Invalid litter component parameters.');
        p.litter_component_rates_s=rates/86400;
        p.litter_component_respiration=rf;
        p.litter_leaf_fractions=leaf;
        p.litter_root_fractions=root;
        p.litter_cwd_fractions=cwd;
    end
    p.source='CoLM_MOD_BGC_Soil_BiogeochemDecompCascadeBGC';
    p.nitrogen_partition='COMPASS_stoichiometric_interface_adaptation';
end

function aom=emptyAom(zero)
% Allocate one material-specific surface plus incorporated AOM C/N store.
% Editor: Jizixing Zhu; Date: 2026-08-28.
    aom=struct('surface_C',0,'surface_N',0, ...
        'element_C',zero,'element_N',zero);
end

function soc=initialSoc(tbl,centersTop,bdTop)
% Convert measured SOC [mg C kg^-1] to [g C cm^-3 bulk soil].
% Editor: Jizixing Zhu; Date: 2026-08-28.
    required={'depth_cm','soc_total_mgkg'};
    assert(istable(tbl)&&all(ismember(required,tbl.Properties.VariableNames)), ...
        'buildSoilCarbon:missingInitialSoc', ...
        'soil_initial_profile.csv requires depth_cm and soc_total_mgkg.');
    fraction=init.interpLinearWithExtrapolation(tbl.depth_cm(:), ...
        tbl.soc_total_mgkg(:)./1e6,centersTop(:));
    assert(all(isfinite(fraction))&&all(fraction>0), ...
        'buildSoilCarbon:badInitialSoc', ...
        'Interpolated SOC mass fraction must be finite and positive.');
    soc=fraction.*bdTop;
end

function f=initializationFractions(tbl,centersTop,dzTop)
% Interpolate native POC/MOC/MB/IOM initialization fractions.
% Editor: Jizixing Zhu; Date: 2026-08-28.
% Change: 2026-08-30 Codex -- hold endpoint fractions outside sampled
%   depths; linear extrapolation created an artificial surface POC excess.
    names={'frac_POC','frac_MOC','frac_MB','frac_IOM'};
    assert(all(ismember([{'depth_cm'},names],tbl.Properties.VariableNames)), ...
        'buildSoilCarbon:badFractionColumns', ...
        'The native migration prior table is missing required columns.');
    f=zeros(numel(centersTop),numel(names));
    % Change: 2026-09-08 Jizixing Zhu -- explicit horizon priors are
    % thickness weighted, while legacy point profiles retain interpolation.
    if ismember('depth_top_cm',tbl.Properties.VariableNames)
        a=tbl.depth_top_cm(:); b=tbl.depth_cm(:); v=tbl{:,names};
        assert(~isempty(a)&&all(isfinite([a;b;v(:)])) ...
            &&all(a>=0)&&all(b>a)&&a(1)==0 ...
            &&all(abs(a(2:end)-b(1:end-1))<1e-10) ...
            &&all(v(:)>=0&v(:)<=1)&&all(abs(sum(v,2)-1)<1e-10), ...
            'buildSoilCarbon:badHorizons','Invalid contiguous pool-fraction horizons.');
        edges=[0;cumsum(dzTop(:))];
        assert(b(end)>=edges(end)-1e-10,'buildSoilCarbon:horizonCoverage', ...
            'Pool-fraction horizons must cover the model column.');
        for j=1:numel(centersTop)
            overlap=max(0,min(b,edges(j+1))-max(a,edges(j)));
            assert(abs(sum(overlap)-dzTop(j))<1e-9, ...
                'buildSoilCarbon:horizonCoverage','Incomplete layer coverage.');
            f(j,:)=overlap'*v/dzTop(j);
        end
        return;
    end
    for i=1:numel(names)
        f(:,i)=interp1(tbl.depth_cm(:),tbl.(names{i})(:), ...
            centersTop(:),'linear','extrap');
        f(centersTop(:)<tbl.depth_cm(1),i)=tbl.(names{i})(1);
        f(centersTop(:)>tbl.depth_cm(end),i)=tbl.(names{i})(end);
    end
    f=min(max(f,0),1);
    assert(all(abs(sum(f,2)-1)<=1e-10), ...
        'buildSoilCarbon:fractionClosure', ...
        'Native initialization fractions must sum to measured SOC.');
end

function cn=buildCnParameters(cfg,microbe)
% Validate pool-specific initial/target C:N ratios [g C g^-1 N].
% Editor: Jizixing Zhu; Date: 2026-08-28.
% Change: 2026-09-12 Jizixing Zhu -- add source-MEND POM/enzyme and
%   microbial C:N bounds.
    names={'POC','MOC','QOC','MB','EP','EM','NEC'};
    cn=struct();
    for i=1:numel(names)
        key=['CN_' names{i} '0'];
        assert(isfield(cfg,key)&&isnumeric(cfg.(key))&&isscalar(cfg.(key)) ...
            &&isfinite(cfg.(key))&&cfg.(key)>0, ...
            'buildSoilCarbon:badCnParameter', ...
            'organic_n.%s must be a finite positive scalar.',key);
        cn.(names{i})=double(cfg.(key));
    end
    required={'mend_lci_initial','mend_cn_pom2','mend_cn_mb_min', ...
        'mend_cn_mb_max','mend_cn_enzp','mend_cn_enzm'};
    for i=1:numel(required)
        key=required{i};
        assert(isfield(microbe,key)&&isnumeric(microbe.(key)) ...
            &&isscalar(microbe.(key))&&isfinite(microbe.(key)) ...
            &&microbe.(key)>0,'buildSoilCarbon:badMendCn', ...
            'microbe.%s must be a finite positive scalar.',key);
    end
    lci=double(microbe.mend_lci_initial);
    cn.POM2=double(microbe.mend_cn_pom2);
    denominator=1/cn.POC-(1-lci)/cn.POM2;
    assert(lci>0&&lci<1&&denominator>0, ...
        'buildSoilCarbon:incompatiblePocCn', ...
        'MEND POM split and CN_POC0 imply a nonpositive POM1 C:N.');
    cn.POM1=lci/denominator;
    cn.MBA_min=double(microbe.mend_cn_mb_min);
    cn.MBA_max=double(microbe.mend_cn_mb_max);
    assert(cn.MBA_max>cn.MBA_min,'buildSoilCarbon:badMendCn', ...
        'MEND microbial maximum C:N must exceed minimum C:N.');
    cn.ENZP1=double(microbe.mend_cn_enzp);
    cn.ENZP2=double(microbe.mend_cn_enzp);
    cn.ENZM=double(microbe.mend_cn_enzm);
end

function value=columnNitrogen(state,dzBottom)
% Total initialized organic N column [g N cm^-2].
% Editor: Jizixing Zhu; Date: 2026-08-28.
    profile=state.PON1+state.PON2+state.MON+state.QON ...
        +state.MB_N+state.BD_N+state.ENZP1_N+state.ENZP2_N ...
        +state.EM_N+state.NEC_N;
    value=sum(profile(:).*dzBottom(:));
end

function validateState(state,n)
% Enforce finite nonnegative paired state shapes before checkpoint ownership.
% Editor: Jizixing Zhu; Date: 2026-08-28.
    fields={'PON1','PON2','PON','MON','QON','MB_N','BD_N','MB_mobile_N', ...
        'MB_attached_N','MB_retained_N','BD_mobile_N','BD_attached_N', ...
        'BD_retained_N','ENZP1_N','ENZP2_N','EP_N','EM_N', ...
        'NEC_C','NEC_N','IOM_C'};
    values=[];
    for i=1:numel(fields)
        value=state.(fields{i});
        assert(isequal(size(value),[n,1]), ...
            'buildSoilCarbon:stateShape','%s must be NL-by-1.',fields{i});
        values=[values;value(:)]; %#ok<AGROW>
    end
    materials={'litter','straw','manure'};
    for i=1:numel(materials)
        a=state.aom.(materials{i});
        values=[values;a.surface_C;a.surface_N;a.element_C(:);a.element_N(:)]; %#ok<AGROW>
    end
    assert(all(isfinite(values))&&all(values>=0), ...
        'buildSoilCarbon:invalidState', ...
        'Native organic C/N state must be finite and nonnegative.');
    assert(all(abs(state.PON-state.PON1-state.PON2)<=1e-14) ...
        &&all(abs(state.EP_N-state.ENZP1_N-state.ENZP2_N)<=1e-14), ...
        'buildSoilCarbon:mendAggregate', ...
        'MEND organic-N compatibility aggregates are inconsistent.');
end
