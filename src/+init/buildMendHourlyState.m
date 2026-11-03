function built=buildMendHourlyState(built,config)
% Editor: Wenjie Liu; Date: 2026-10-17.
%BUILDMENDHOURLYSTATE Initialize explicit source-hour owners without new mass.
% Called after native organic, DOM and microbial builders. Dynamic phase
% stores are bottom-up g/cm3; ENZN C is top-down mg/g, paired N bottom-up
% g/cm3. Preserve independent initial D/Q; zero ENZN adds no C/N. Schemas
% differ so old equilibrium checkpoints cannot silently enter this mode.
% Source base sorption parameters are explicit site inputs, not fitted here.
% Editor: Jizixing Zhu; Date: 2026-09-14.
% Change: 2026-09-16 Jizixing Zhu -- explicit DOM sorption selection lets
% source-hour reactions retain the COMPASS linear-equilibrium total owner.
% Change: 2026-10-17 Wenjie Liu -- preserve resolved ENZN C/N and phases;
%   allocate zero only for missing cold-start fields.
% Change: 2026-10-17 Wenjie Liu -- support configured microbial C:N bounds
%   and a common positive enzyme C:N, retaining source defaults.
% Change: 2026-10-23 Wenjie Liu -- require prebuilt dynamic Q owners from
%   independent source priors; remove the former theta/Kd initial repartition.
% Change: 2026-10-31 Wenjie Liu -- use resolved soil-or-site sorption references,
%   preserving scalars and top-down NL-by-1 profiles without duplicate defaults.
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
    assert(strcmp(p.kinetics,'MM')&&p.CN_MB_min>0&&p.CN_MB_max>p.CN_MB_min ...
        &&all(p.CN_ENZP==p.CN_ENZM)&&p.CN_ENZM>0,'buildMendHourlyState:subset', ...
        'Source-hour mode requires MM, ordered microbial C:N bounds and common enzyme C:N.');
    names={'mend_fp_enzn','mend_ihr'};
    for j=1:numel(names)
        key=names{j};assert(isfield(cfg,key)&&isnumeric(cfg.(key)) ...
            &&isscalar(cfg.(key))&&isfinite(cfg.(key))&&cfg.(key)>=0, ...
            'buildMendHourlyState:parameter','microbe.%s must be supplied and nonnegative.',key);
    end
    assert(ismember(cfg.mend_ihr,[0 1]), ...
        'buildMendHourlyState:parameter','iHR=0/1 required.');
    refs={'Qmax_ref_mg_g','Kba_ref_g_mg','Kdes_ref_mg_g_h'};
    for j=1:numel(refs)
        key=refs{j};
        assert(isfield(p,key)&&(isscalar(p.(key)) ...
            ||isequal(size(p.(key)),[built.domain.NL,1])), ...
            'buildMendHourlyState:parameter','Resolved %s must be scalar or top-down NL-by-1.',key);
    end
    p.hourly=struct('Qmax',p.Qmax_ref_mg_g,'Kba',p.Kba_ref_g_mg, ...
        'Kdes_ref',p.Kdes_ref_mg_g_h,'fpENZN',cfg.mend_fp_enzn,'iHR',cfg.mend_ihr);
    p.hourly.no2_policy='unresolved_species_excluded_from_enzyme_allocation';
    built.parameters.soil_microbes=p;
    d=built.Y.soil_dom;dp=built.parameters.soil_dom;
    if strcmp(mode,'mend_dynamic')
        assert(strcmp(dp.sorption_model,'mend_dynamic') ...
            &&isfield(d,'sorption_model')&&strcmp(d.sorption_model,'mend_dynamic') ...
            &&d.schema_version>=2,'buildMendHourlyState:dynamicPriors', ...
            'Dynamic DOM must already own independent initial D/Q and DON/QON.');
        for kind={'doc','don'}
            k=kind{1};qkey=[k '_sorbed_g_cm3'];total=d.([k '_total_g_cm3']);
            assert(isfield(d,qkey)&&isequal(size(d.(qkey)),size(total)) ...
                &&all(isfinite(d.(qkey)))&&all(d.(qkey)>=0)&&all(d.(qkey)<=total), ...
                'buildMendHourlyState:dynamicPriors','Invalid independent initial %s Q.',k);
        end
        assert(all(built.Y.soil_microbes.QOC==0)&&all(built.Y.soil_organic.QON==0), ...
            'buildMendHourlyState:duplicateQ','Only soil_dom may own source-hour sorbed C/N.');
        built.parameters.soil_microbes.hourly.sorption_policy='source_dynamic_independent_phases';
    else
        % Only total DOC/DON are prognostic; phases are algebraic at the
        % current theta and independently configured DOC/DON Kd profiles.
        built.parameters.soil_microbes.hourly.sorption_policy='linear_equilibrium_stage_split';
    end
    built.Y.soil_dom.sorption_model=mode;
    built.parameters.soil_dom.sorption_model=mode;
    z=zeros(built.domain.NL,1);
    if ~isfield(built.Y.soil_microbes,'ENZN'),built.Y.soil_microbes.ENZN=z;end
    if ~isfield(built.Y.soil_organic,'ENZN_N'),built.Y.soil_organic.ENZN_N=z;end
    for phase={'mobile','attached','retained'}
        ck=['ENZN_' phase{1}];nk=[ck '_N'];
        if ~isfield(built.Y.soil_microbes,ck),built.Y.soil_microbes.(ck)=z;end
        if ~isfield(built.Y.soil_organic,nk),built.Y.soil_organic.(nk)=z;end
    end
    built.Y.soil_microbes.schema_version=4;built.Y.soil_organic.schema_version=5;
end
