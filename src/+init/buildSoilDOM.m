function [state, parameters, freeTarget] = buildSoilDOM(raw, built, ~, deferFree)
% Editor: Wenjie Liu; Date: 2026-10-17.
%BUILDSOILDOM Build accepted DOC/DON total-storage state on COMPASS nodes.
%
%   DOC and DON accepted states are bottom-up element total mass
%   concentrations per bulk-soil volume [g cm^-3 bulk soil]. Aqueous
%   concentrations are diagnostics. mend_dynamic preserves independently
%   supplied sorbed Q and recovers free D=total-Q. Explicit linear_equilibrium
%   instead uses total=(theta_l+rho_b*Kd)*C_aq. Missing legacy doc_mgL/don_mgL
%   columns initialize zero only in the equilibrium cold-start path.
%
%   Editor: Jizixing Zhu; Date: 2026-08-21.
%   Change: 2026-09-01 Codex -- remove the retired independent DON Michaelis-
%   Menten sink; DON consumption is the C:N companion of MEMC DOC uptake.
%   Change: 2026-09-02 Codex -- read depth-resolved sorption and dispersion
%   from soil.csv while retaining numerical transport controls in [dom].
%   Change: 2026-09-29 Wenjie Liu -- initialize numbered organic boundary
%     policy and disjoint input/top-export cumulative memory for restart.
%   Change: 2026-10-12 Wenjie Liu -- explicit doc_free_mgkg means D only;
%     preserve historical doc_mgkg=D+Q and reserve reconstructed total C.
%   Change: 2026-10-12 Wenjie Liu -- the facade defers explicit free-D
%     reconstruction until its final liquid-water state. freeTarget is a
%     bottom-up mg C/kg dry-soil initialization request, not accepted storage.
%     Direct three-argument calls retain the supplied-theta algebra contract.
% Change: 2026-10-17 Wenjie Liu -- resolved CSV DOM is already SOC-funded;
%   build once without a free-D request or another MOC/MON debit.
% Change: 2026-10-23 Wenjie Liu -- initialize dynamic total and independent
%   sorbed inventories directly from detailed priors. Reject missing Q priors
%   before the legacy Kd/free-D reconstruction; no inferred equilibrium Q.
    if nargin<4,deferFree=false;end
    cfg = requiredSection(raw.config, 'dom');
    mode='linear_equilibrium';
    if isfield(raw.config,'microbe')&&isfield(raw.config.microbe,'reaction_method') ...
            &&strcmp(string(raw.config.microbe.reaction_method),"source_hourly")
        mode='mend_dynamic';
    end
    if isfield(cfg,'sorption_model'),mode=char(string(cfg.sorption_model));end
    assert(ismember(string(mode),["linear_equilibrium","mend_dynamic"]), ...
        'buildSoilDOM:sorption','Unknown dom.sorption_model.');
    names = {'diffusion_cm2_s','transport_substep_s'};
    for i = 1:numel(names)
        assert(isfield(cfg, names{i}), 'buildSoilDOM:missingParameter', ...
            'dom.%s is required.', names{i});
    end
    soilParameters=built.parameters.soil;
    parameters = struct( ...
        'kd_doc_cm3_g', soilProfile(soilParameters,'dom_kd_doc_cm3_g',true), ...
        'kd_don_cm3_g', soilProfile(soilParameters,'dom_kd_don_cm3_g',true), ...
        'diffusion_cm2_s', positiveOrZero(cfg.diffusion_cm2_s, 'diffusion_cm2_s'), ...
        'dispersion_cm', soilProfile(soilParameters,'dom_dispersion_cm',false), ...
        'transport_substep_s', positiveScalar(cfg.transport_substep_s, 'transport_substep_s'), ...
        'sorption_model', mode, ...
        'direct_don_to_no3', false);

    % Change: 2026-09-29 Wenjie Liu -- boundary policy is site-wide; only
    % enabled open boundaries own checkpointed cumulative export memory.
    parameters.boundary=soil.dom.normalizeBoundary(raw.config.boundary);

    theta = 0.5 * (built.Y.theta(1:end-1) + built.Y.theta(2:end));
    bulkDensity = built.parameters.soil.bulk_density(:);
    centerDepthBottomUp = 0.5 * (built.domain.z_node(1:end-1) ...
        + built.domain.z_node(2:end));
    if ismember('doc_free_mgkg',raw.initial.Properties.VariableNames) ...
            &&~all(ismissing(raw.initial.doc_free_mgkg)) ...
            &&isfield(raw.config,'switches')&&isfield(raw.config.switches,'lFreezing') ...
            &&raw.config.switches.lFreezing==1
        assert(isfield(built.Y,'T')&&all(isfinite(built.Y.T))&&all(built.Y.T>0), ...
            'buildSoilDOM:frozenFreeInitial', ...
            ['doc_free_mgkg currently requires an unfrozen (>0 C) cold start. ' ...
             'Frozen free-D targets require explicit post-freezing initialization.']);
    end
    % Change: 2026-09-09 Jizixing Zhu -- mg/kg dry-soil TOTAL DOM input.
    % Legacy mg/L remains aqueous; never interpret the same number as both.
    detailed=isfield(built.parameters,'soil_organic') ...
        &&isfield(built.parameters.soil_organic,'initial_dom');
    if detailed
        d=built.parameters.soil_organic.initial_dom;
        docTotal=d.doc_total_g_cm3;donTotal=d.don_total_g_cm3;
        docBlank=false;donBlank=false;docMassBasis=true;donMassBasis=true;freeTarget=[];
    else
        assert(~strcmp(mode,'mend_dynamic'),'buildSoilDOM:missingDynamicPriors', ...
            ['mend_dynamic requires detailed independent DOC/QOC and DON/QON priors; ', ...
             'total or free DOC alone cannot determine source Q without an additional assumption.']);
        [docTotal,docBlank,docMassBasis,freeTarget]=initialMass(raw.initial,centerDepthBottomUp, ...
            'doc',theta,bulkDensity,parameters.kd_doc_cm3_g,deferFree);
        [donTotal,donBlank,donMassBasis]=initialMass(raw.initial,centerDepthBottomUp, ...
            'don',theta,bulkDensity,parameters.kd_don_cm3_g,false);
    end
    state = struct( ...
        'schema_version', 1, ...
        'doc_total_g_cm3', ...
            docTotal, ...
        'don_total_g_cm3', ...
            donTotal, ...
        'initial_doc_total_mgkg_basis',docMassBasis, ...
        'initial_don_total_mgkg_basis',donMassBasis, ...
        'initialized_from_blank', docBlank || donBlank, ...
        'initial_doc_blank', docBlank, ...
        'initial_don_blank', donBlank, ...
        'direct_don_to_no3', false);
    if strcmp(mode,'mend_dynamic')
        assert(isfield(d,'sorption_model')&&strcmp(d.sorption_model,'mend_dynamic') ...
            &&isfield(d,'doc_sorbed_g_cm3')&&isfield(d,'don_sorbed_g_cm3'), ...
            'buildSoilDOM:missingDynamicPriors', ...
            'Dynamic initial DOC/QOC and DON/QON must be explicitly supplied.');
        state.sorption_model=mode;state.schema_version=2;
        for kind={'doc','don'}
            k=kind{1};q=d.([k '_sorbed_g_cm3']);total=state.([k '_total_g_cm3']);
            assert(isequal(size(q),size(total))&&all(isfinite(q)) ...
                &&all(q>=0)&&all(q<=total),'buildSoilDOM:dynamicPrior', ...
                'Initial %s sorbed inventory must lie between zero and total.',k);
            state.([k '_sorbed_g_cm3'])=q;
        end
    end
    if ~strcmp(parameters.boundary.bottom_kind,'closed')
        state.bottom_export_cumulative=struct('schema_version',1,'doc_g_cm2',0,'don_g_cm2',0);
    end
    if parameters.boundary.schema_version==2
        state.boundary_exchange_cumulative=soil.dom.emptyBoundaryExchange({'doc','don'});
    end
end

function [mass,blank,massBasis,freeTarget]=initialMass(tbl,z,prefix,theta,bd,kd,deferFree)
% Total mg/kg * rho[g/cm3] * 1e-6 -> g/cm3 bulk soil.
    freeTarget=[];
    modern=[prefix '_mgkg'];legacy=[prefix '_mgL'];
    [v,blank]=initialConcentration(tbl,z,modern);
    [old,oldBlank]=initialConcentration(tbl,z,legacy);
    if strcmp(prefix,'doc')
        [free,freeBlank]=initialConcentration(tbl,z,'doc_free_mgkg');
        assert(freeBlank||(blank&&oldBlank),'buildSoilDOM:ambiguousUnits', ...
            'Use only one of doc_free_mgkg (D), doc_mgkg (D+Q), or doc_mgL (water).');
        if ~freeBlank
            freeTarget=free;
            if deferFree
                % Unfunded scaffold only: buildState records the pending
                % request; initialize must finalize it before returning.
                mass=zeros(size(free));
            else
                mass=soil.dom.totalFromFree(free.*bd*1e-6,theta,bd,kd);
            end
            blank=false;massBasis=true; % Reserve reconstructed D+Q, not D alone.
            return;
        end
    end
    assert(blank||oldBlank,'buildSoilDOM:ambiguousUnits', ...
        'Do not populate both %s and %s.',modern,legacy);
    massBasis=~blank;
    if massBasis
        mass=v.*bd*1e-6;
    else
        blank=oldBlank;
        mass=(theta+bd.*kd).*old*1e-6;
    end
end

function value=soilProfile(soil,name,mustBePositive)
% Read one depth-resolved bottom-up parameter profile from soil.csv.
    assert(isfield(soil,name),'buildSoilDOM:missingSoilProfile', ...
        'soil.csv parameter %s is required.',name);
    value=soil.(name)(:);
    assert(all(isfinite(value))&&all(value>=0), ...
        'buildSoilDOM:badSoilProfile', ...
        'soil.csv:%s must be finite and nonnegative.',name);
    if mustBePositive
        assert(all(value>0),'buildSoilDOM:badSoilProfile', ...
            'soil.csv:%s must be positive in every layer.',name);
    end
end

function section = requiredSection(config, name)
    if ~isfield(config, name) || ~isstruct(config.(name))
        error('buildSoilDOM:missingSection', '[%s] is required.', name);
    end
    section = config.(name);
end

function [value, blank] = initialConcentration(tbl, zBottomUp, name)
    blank = ~ismember(name, tbl.Properties.VariableNames) ...
        || all(ismissing(tbl.(name)));
    if blank
        value = zeros(size(zBottomUp));
        return;
    end
    samples = tbl.(name)(:);
    if any(~isfinite(samples)) || any(samples < 0)
        error('buildSoilDOM:badInitialProfile', ...
            '%s must be finite and nonnegative when supplied.', name);
    end
    valueTopDown = init.interpLinearWithExtrapolation( ...
        tbl.depth_cm(:), samples, flip(zBottomUp));
    value = flip(valueTopDown(:));
    assert(all(isfinite(value))&&all(value>=0),'buildSoilDOM:badInitialProfile', ...
        'Interpolated %s must remain finite and nonnegative.',name);
end

function value = positiveOrZero(value, name)
    assert(isnumeric(value) && isscalar(value) && isfinite(value) && value >= 0, ...
        'buildSoilDOM:badParameter', 'dom.%s must be finite and >= 0.', name);
    value = double(value);
end

function value = positiveScalar(value, name)
    assert(isnumeric(value) && isscalar(value) && isfinite(value) && value > 0, ...
        'buildSoilDOM:badParameter', 'dom.%s must be finite and > 0.', name);
    value = double(value);
end
