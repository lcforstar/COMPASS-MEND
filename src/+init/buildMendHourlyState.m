function built=buildMendHourlyState(built,config)
%BUILDMENDHOURLYSTATE Initialize explicit source-hour owners without new mass.
% Called after native organic, DOM and microbial builders. Dynamic phase
% stores are bottom-up g/cm3; ENZN C is top-down mg/g, paired N bottom-up
% g/cm3. Partition initial total DOM once; zero ENZN adds no C/N. Schemas
% differ so old equilibrium checkpoints cannot silently enter this mode.
% Source base sorption parameters are explicit site inputs, not fitted here.
% Editor: Jizixing Zhu; Date: 2026-09-14.
% Change: 2026-09-16 Jizixing Zhu -- explicit DOM sorption selection lets
% source-hour reactions retain the COMPASS linear-equilibrium total owner.
    if ~isfield(config,'microbe')||~isfield(built.Y,'soil_microbes'),return;end
    cfg=config.microbe;method='rk4_equilibrium';
    if isfield(cfg,'reaction_method'),method=char(string(cfg.reaction_method));end
    assert(ismember(string(method),["rk4_equilibrium","source_hourly"]), ...
        'buildMendHourlyState:method','Unknown microbe.reaction_method.');
    built.parameters.soil_microbes.reaction_method=method;
    mode='mend_dynamic';
    if strcmp(method,'rk4_equilibrium'),mode='linear_equilibrium';end
    if isfield(config,'dom')&&isfield(config.dom,'sorption_model')
        mode=char(string(config.dom.sorption_model));
    end
    assert(ismember(string(mode),["linear_equilibrium","mend_dynamic"]), ...
        'buildMendHourlyState:sorption','Unknown dom.sorption_model.');
    if strcmp(method,'rk4_equilibrium')
        assert(strcmp(mode,'linear_equilibrium'),'buildMendHourlyState:sorption', ...
            'rk4_equilibrium requires linear_equilibrium DOM.');
        return;
    end
    p=built.parameters.soil_microbes;
    assert(strcmp(p.kinetics,'MM')&&p.CN_MB_min==2&&p.CN_MB_max==14 ...
        &&all(p.CN_ENZP==3)&&p.CN_ENZM==3,'buildMendHourlyState:subset', ...
        'Source-hour mode requires MM, microbial CN 2--14 and enzyme CN 3.');
    names={'mend_qmax_mg_g','mend_kba_g_mg','mend_kdes_mg_g_h','mend_fp_enzn','mend_ihr'};
    for j=1:numel(names)
        key=names{j};assert(isfield(cfg,key)&&isnumeric(cfg.(key)) ...
            &&isscalar(cfg.(key))&&isfinite(cfg.(key))&&cfg.(key)>=0, ...
            'buildMendHourlyState:parameter','microbe.%s must be supplied and nonnegative.',key);
    end
    assert(cfg.mend_qmax_mg_g>0&&ismember(cfg.mend_ihr,[0 1]), ...
        'buildMendHourlyState:parameter','Positive Qmax and iHR=0/1 required.');
    p.hourly=struct('Qmax',cfg.mend_qmax_mg_g,'Kba',cfg.mend_kba_g_mg, ...
        'Kdes_ref',cfg.mend_kdes_mg_g_h,'fpENZN',cfg.mend_fp_enzn,'iHR',cfg.mend_ihr);
    p.hourly.no2_policy='unresolved_species_excluded_from_enzyme_allocation';
    built.parameters.soil_microbes=p;
    d=built.Y.soil_dom;dp=built.parameters.soil_dom;
    theta=.5*(built.Y.theta(1:end-1)+built.Y.theta(2:end));bd=built.parameters.soil.bulk_density(:);
    if strcmp(mode,'mend_dynamic')
        for kind={'doc','don'}
            k=kind{1};phase=soil.dom.phaseState(d,k,theta,bd,dp);
            built.Y.soil_dom.([k '_sorbed_g_cm3'])=phase.sorbed_g_cm3_bulk;
        end
        built.Y.soil_dom.schema_version=2;
    else
        % Only total DOC/DON are prognostic; phases are algebraic at the
        % current theta and independently configured DOC/DON Kd profiles.
        built.parameters.soil_microbes.hourly.sorption_policy='linear_equilibrium_stage_split';
    end
    built.Y.soil_dom.sorption_model=mode;
    built.parameters.soil_dom.sorption_model=mode;
    z=zeros(built.domain.NL,1);
    built.Y.soil_microbes.ENZN=z;built.Y.soil_organic.ENZN_N=z;
    for phase={'mobile','attached','retained'}
        built.Y.soil_microbes.(['ENZN_' phase{1}])=z;
        built.Y.soil_organic.(['ENZN_' phase{1} '_N'])=z;
    end
    built.Y.soil_microbes.schema_version=4;built.Y.soil_organic.schema_version=5;
end
