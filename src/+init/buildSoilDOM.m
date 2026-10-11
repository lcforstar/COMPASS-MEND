function [state, parameters] = buildSoilDOM(raw, built, ~)
%BUILDSOILDOM Build accepted DOC/DON total-storage state on COMPASS nodes.
%
%   DOC and DON accepted states are bottom-up element total mass
%   concentrations per bulk-soil volume [g cm^-3 bulk soil]. Aqueous
%   concentrations are diagnostic closure
%   variables recovered by instantaneous linear equilibrium sorption:
%       S_X = (theta_l + rho_b Kd_X) C_X.
%   Thus M_X/S_X is the prognostic store and C_X is never a second accepted
%   state. Missing doc_mgL/don_mgL columns initialize exact zero and are marked
%   for the spin-up guard in beginRun.
%
%   Editor: Jizixing Zhu; Date: 2026-08-21.
%   Change: 2026-09-01 Codex -- remove the retired independent DON Michaelis-
%   Menten sink; DON consumption is the C:N companion of MEMC DOC uptake.
%   Change: 2026-09-02 Codex -- read depth-resolved sorption and dispersion
%   from soil.csv while retaining numerical transport controls in [dom].
%   Change: 2026-09-29 Wenjie Liu -- initialize numbered organic boundary
%     policy and disjoint input/top-export cumulative memory for restart.
    cfg = requiredSection(raw.config, 'dom');
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
        'sorption_model', 'linear_equilibrium', ...
        'direct_don_to_no3', false);

    % Change: 2026-09-29 Wenjie Liu -- boundary policy is site-wide; only
    % enabled open boundaries own checkpointed cumulative export memory.
    parameters.boundary=soil.dom.normalizeBoundary(raw.config.boundary);

    theta = 0.5 * (built.Y.theta(1:end-1) + built.Y.theta(2:end));
    bulkDensity = built.parameters.soil.bulk_density(:);
    centerDepthBottomUp = 0.5 * (built.domain.z_node(1:end-1) ...
        + built.domain.z_node(2:end));
    % Change: 2026-09-09 Jizixing Zhu -- mg/kg dry-soil TOTAL DOM input.
    % Legacy mg/L remains aqueous; never interpret the same number as both.
    [docTotal,docBlank,docMassBasis]=initialMass(raw.initial,centerDepthBottomUp, ...
        'doc',theta,bulkDensity,parameters.kd_doc_cm3_g);
    [donTotal,donBlank,donMassBasis]=initialMass(raw.initial,centerDepthBottomUp, ...
        'don',theta,bulkDensity,parameters.kd_don_cm3_g);
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
    if ~strcmp(parameters.boundary.bottom_kind,'closed')
        state.bottom_export_cumulative=struct('schema_version',1,'doc_g_cm2',0,'don_g_cm2',0);
    end
    if parameters.boundary.schema_version==2
        state.boundary_exchange_cumulative=soil.dom.emptyBoundaryExchange({'doc','don'});
    end
end

function [mass,blank,massBasis]=initialMass(tbl,z,prefix,theta,bd,kd)
% Total mg/kg * rho[g/cm3] * 1e-6 -> g/cm3 bulk soil.
    modern=[prefix '_mgkg'];legacy=[prefix '_mgL'];
    [v,blank]=initialConcentration(tbl,z,modern);
    [old,oldBlank]=initialConcentration(tbl,z,legacy);
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
