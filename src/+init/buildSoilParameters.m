function soil = buildSoilParameters(rawSoil, domain, rawConfig)
%BUILDSOILPARAMETERS  Map soil.csv rows onto the model layers (stage B).
%
%   soil = init.buildSoilParameters(rawSoil, domain, rawConfig) maps soil.csv
%   values onto model-layer centers according to
%   rawConfig.soil_input.vertical_parameter_mode:
%     horizon_bounds  depth_cm is the lower boundary of a piecewise-constant
%                     horizon (legacy COMPASS behavior).
%     sample_linear   depth_cm is a measurement coordinate; primitive values
%                     are linearly interpolated to layer centers and endpoint
%                     values are held outside the sampled range (default).
%
%   Inputs:
%     rawSoil  table   raw.soil from io.loadSite (one row per horizon boundary
%                      or sampled depth; canonical column names from
%                      io.adapters.tableSpec('soil'))
%     domain   struct  built.domain from init.buildState (uses domain.sumdz, the
%                      top-down cumulative layer depth [cm], == legacy SUMDELTZ)
%     rawConfig struct  raw.config; optional [soil_input] mode selector
%
%   Returns soil: struct of [NL x 1] per-layer arrays, BOTTOM-UP (layer 1 =
%   deepest). Units (match legacy SoilProperties):
%     theta_s,theta_r,msoc_volume [m3/m3]; alpha [1/cm]; n,m [-];  (porosity is
%       NOT set here -- it follows the bulk-density decision in buildBulkDensity)
%     ks [cm/s]; clay,sand,silt [-]; msoc [kg/kg]; bulk_density [g/cm3]; ph [-];
%     dif_nh4,dif_no3 [cm2/s]; disp [cm]; kg [cm3/g];
%     k_mineralization,kurea,vmax_nit,vmax_denit,k_nh4_dissolve [1/s];
%     km_nit,km_denit [ugN/cm3].
%     fmax [-] (TOPMODEL max saturated fraction, 0=disabled); ks0 [cm/s] top-layer Ks.
%
%   References:
%     legacy SoilProperties.load_from_soilproperty_xlsx and distribute_to_nodes
%       (unit factors and horizon_bounds mapping).
%     ClimaLand docs/src/tutorials/standalone/Soil/layered_soil.jl:46-71
%       (sampled hydraulic parameters mapped with SpaceVaryingInput).
% Editor: Jizixing Zhu; Date: 2026-06-27.
% Change: 2026-07-11 Jizixing Zhu -- keep distributed layer parameters column-shaped for a single homogeneous horizon.
% Change: 2026-06-28 Jizixing Zhu — add surface-runoff params soil.fmax/soil.ks0 (from setBoundaryCondition).
% Change: 2026-07-12 Jizixing Zhu -- add canonical ClimaLand-style sample-linear soil.csv mapping.
% Change: 2026-07-16 Jizixing Zhu -- document the descriptive soil initialization filename.
% Change: 2026-09-02 Codex -- keep soil.csv depth-resolved: site-global
%   AOM/DOM/MEMC/MEND/BIOPARTICLE parameters are owned by site.toml.
% Change: 2026-09-29 Wenjie Liu -- map four MEND kinetic columns by the same
%   public depth contract; retain h^-1 and mg C/g (no day/second conversion).
% Change: 2026-10-12 Wenjie Liu -- add an independent MEND group-constant
%   mapping so hydraulic and inorganic-N profiles retain their configured map.

    if nargin < 3 || isempty(rawConfig)
        rawConfig = struct();
    end
    verticalMode = readVerticalParameterMode(rawConfig);
    mendMode = readMendVerticalParameterMode(rawConfig,verticalMode);
    cols = rawSoil.Properties.VariableNames;
    has  = @(c) ismember(c, cols);

    % --- per-row primitive values with unit conversion (1 x NS samples) ---
    T = struct();
    % hydraulic (van Genuchten)
    T.theta_s = rawSoil.theta_s(:)';            % [m3/m3]
    T.theta_r = rawSoil.theta_r(:)';            % [m3/m3]
    T.alpha   = rawSoil.vg_alpha_1cm(:)';        % [1/cm]
    T.n       = rawSoil.vg_n(:)';                % [-]
    T.ks      = rawSoil.ks_cmd(:)' / 86400;      % [cm/s]   <- cm/d
    % texture / structure
    T.clay         = rawSoil.clay_pct(:)' / 100;          % [-] fraction <- %
    T.sand         = rawSoil.sand_pct(:)' / 100;          % [-]          <- %
    T.bulk_density = rawSoil.bulk_density_gcm3(:)';        % [g/cm3]
    if has('msoc_mgkg'); T.msoc = rawSoil.msoc_mgkg(:)' / 1e6; end   % [kg/kg] <- mg/kg
    if has('ph');        T.ph   = rawSoil.ph(:)'; end                % [-]
    % N-cycle transport + rates (cm2/d -> cm2/s, 1/d -> 1/s)
    if has('difnh4_cm2d');         T.dif_nh4          = rawSoil.difnh4_cm2d(:)'         / 86400; end
    if has('difno3_cm2d');         T.dif_no3          = rawSoil.difno3_cm2d(:)'         / 86400; end
    if has('disp_cm');             T.disp             = rawSoil.disp_cm(:)';                     end
    if has('kg_cm3g');             T.kg               = rawSoil.kg_cm3g(:)';                     end
    if has('kmineralization_1d');  T.k_mineralization = rawSoil.kmineralization_1d(:)'  / 86400; end
    if has('vmax_nit_ugn_cm3d');   T.vmax_nit         = rawSoil.vmax_nit_ugn_cm3d(:)'   / 86400; end
    if has('km_nit_ugn_cm3');      T.km_nit           = rawSoil.km_nit_ugn_cm3(:)';              end
    if has('kurea_1d');            T.kurea            = rawSoil.kurea_1d(:)'            / 86400; end
    if has('vmax_denit_ugn_cm3d'); T.vmax_denit       = rawSoil.vmax_denit_ugn_cm3d(:)' / 86400; end
    if has('km_denit_ugn_cm3');    T.km_denit         = rawSoil.km_denit_ugn_cm3(:)';            end
    if has('k_nh4_dissolve_1d');   T.k_nh4_dissolve   = rawSoil.k_nh4_dissolve_1d(:)'   / 86400; end

    % Depth-resolved DOM retention/dispersion and porous-medium microbial
    % transport parameters. Public first-order rates use d^-1.
    if has('kd_doc_cm3_g');        T.dom_kd_doc_cm3_g = rawSoil.kd_doc_cm3_g(:)'; end
    if has('kd_don_cm3_g');        T.dom_kd_don_cm3_g = rawSoil.kd_don_cm3_g(:)'; end
    if has('dom_disp_cm');          T.dom_dispersion_cm = rawSoil.dom_disp_cm(:)'; end
    if has('attachment_1d');       T.bioparticle_attachment_rate_s = rawSoil.attachment_1d(:)'/86400; end
    if has('detachment_1d');       T.bioparticle_detachment_rate_s = rawSoil.detachment_1d(:)'/86400; end
    if has('straining_1d');        T.bioparticle_straining_rate_s = rawSoil.straining_1d(:)'/86400; end
    if has('microbe_disp_cm');     T.bioparticle_dispersion_cm = rawSoil.microbe_disp_cm(:)'; end
    mend=io.adapters.mendLayerParameterSpec();
    for j=1:numel(mend.names)
        key=mend.names{j};
        if has(key),T.(key)=rawSoil.(key)(:)';end
    end

    % NOTE: porosity is NOT set here. It depends on which bulk density wins
    % (soil.csv vs soil_initial_profile.csv), decided once in init.buildBulkDensity, after
    % which buildState sets soil.bulk_density / soil.porosity. bulk_density below
    % is the soil.csv value, used as the 'soil' source by buildBulkDensity.

    sampleDepth = rawSoil.depth_cm(:);            % [cm] top-down coordinates
    validateSampleDepth(sampleDepth);

    % --- map primitive inputs to layers, then flip to bottom-up ---
    soil = struct();
    fns = fieldnames(T);
    jmap = mapLayerToHorizon(domain.sumdz, sampleDepth'); % legacy lower-face map
    layerThicknessTopDown = flip(domain.dz(:));
    centerDepth = domain.sumdz(:) - layerThicknessTopDown / 2; % [cm]
    for i = 1:numel(fns)
        name = fns{i};
        fieldMode = verticalMode;
        samples = T.(name);
        if ismember(name,mend.names),fieldMode=mendMode;end
        switch fieldMode
            case 'horizon_bounds'
                mapped = samples(jmap);
            case 'sample_linear'
                % ClimaLand layered_soil.jl:46-71 uses sampled coordinates in
                % SpaceVaryingInput; endpoint holding prevents unsupported
                % extrapolation beyond the first and last observations.
                mapped = interpolateHeld(sampleDepth, samples(:), centerDepth);
            case 'group_constant'
                mapped = mapMendGroups(rawSoil,name,samples(:),centerDepth);
        end
        soil.(name) = flip(mapped(:));  % [NL x 1] bottom-up
    end
    soil.mend_layer_mapping=buildMendMappingAudit(rawSoil,mendMode,centerDepth);

    % --- derive nonlinear/dependent quantities after primitive interpolation ---
    soil.m    = 1 - 1 ./ soil.n;                 % [-] van Genuchten m
    soil.silt = 1 - soil.clay - soil.sand;       % [-]
    if isfield(soil, 'msoc')
        soil.msoc_volume = (soil.msoc * 2650) ./ ...
            ((soil.msoc * 2650) + (1 - soil.msoc) * 1300); % [m3/m3]
    end

    % surface-runoff parameters (legacy setBoundaryCondition.m:118-119 from
    % SoilProperties): fmax disabled (run_stemmus_initialize.m:234 sets 0); ks0 is
    % the top-layer saturated conductivity (legacy SoilProperties.m:790 = Ks(end)).
    soil.fmax = 0;                 % [-]    TOPMODEL max saturated fraction (disabled)
    soil.ks0  = soil.ks(end);      % [cm/s] top-layer (surface) saturated conductivity
end

% ---------------------------------------------------------------------------
function mode=readMendVerticalParameterMode(rawConfig,defaultMode)
% MEND-only mapping; omission exactly inherits the established soil mapping.
    mode=defaultMode;
    if ~isfield(rawConfig,'soil_input') ...
            ||~isfield(rawConfig.soil_input,'mend_vertical_parameter_mode')
        return;
    end
    rawMode=rawConfig.soil_input.mend_vertical_parameter_mode;
    if ~(ischar(rawMode)||(isstring(rawMode)&&isscalar(rawMode)))
        error('buildSoilParameters:badMendVerticalMode', ...
            'soil_input.mend_vertical_parameter_mode must be scalar text.');
    end
    mode=lower(strtrim(char(rawMode)));
    if strcmp(mode,'inherit'),mode=defaultMode;end
    if ~ismember(mode,{'horizon_bounds','sample_linear','group_constant'})
        error('buildSoilParameters:badMendVerticalMode', ...
            ['unknown MEND vertical mode "%s" (expected inherit|', ...
             'horizon_bounds|sample_linear|group_constant)'],mode);
    end
end

% ---------------------------------------------------------------------------
function mapped=mapMendGroups(rawSoil,name,samples,queryDepth)
% Each contiguous group uses one value and its deepest row as lower bound.
    [ids,groups,bounds]=validateMendGroups(rawSoil);
    values=zeros(numel(groups),1);
    for j=1:numel(groups)
        selected=ids==groups(j);
        assert(all(samples(selected)==samples(find(selected,1))), ...
            'buildSoilParameters:nonconstantMendGroup', ...
            '%s must be exactly constant within each mend_group_id.',name);
        bounds(j)=max(rawSoil.depth_cm(selected));values(j)=samples(find(selected,1));
    end
    mapped=zeros(size(queryDepth));
    for j=1:numel(queryDepth)
        index=find(queryDepth(j)<=bounds,1);
        if isempty(index),index=numel(bounds);end
        mapped(j)=values(index);
    end
end

% ---------------------------------------------------------------------------
function audit=buildMendMappingAudit(rawSoil,mode,centerDepth)
    audit=struct('schema_version',1,'mode',mode, ...
        'independent_of_general_soil_mapping',true);
    if strcmp(mode,'group_constant')
        [~,groups,bounds]=validateMendGroups(rawSoil);
        mapped=zeros(size(centerDepth));
        for j=1:numel(centerDepth)
            index=find(centerDepth(j)<=bounds,1);
            if isempty(index),index=numel(bounds);end
            mapped(j)=groups(index);
        end
        audit.group_id=groups(:);audit.group_bottom_cm=bounds(:);
        audit.layer_group_id_bottom_up=flip(mapped(:));
    end
end

% ---------------------------------------------------------------------------
function [ids,groups,bounds]=validateMendGroups(rawSoil)
    assert(ismember('mend_group_id',rawSoil.Properties.VariableNames), ...
        'buildSoilParameters:missingMendGroup', ...
        'group_constant requires soil.csv mend_group_id.');
    ids=double(rawSoil.mend_group_id(:));
    assert(all(isfinite(ids)&ids>=1&ids==fix(ids))&&all(diff(ids)>=0) ...
        &&isequal(unique(ids,'stable'),(1:max(ids))'), ...
        'buildSoilParameters:badMendGroup', ...
        'mend_group_id must be nondecreasing contiguous integers starting at one.');
    groups=unique(ids,'stable');
    bounds=arrayfun(@(g)max(rawSoil.depth_cm(ids==g)),groups(:));
end

% ---------------------------------------------------------------------------
function mode = readVerticalParameterMode(rawConfig)
% Read the public soil.csv vertical-coordinate contract.
% Editor: Jizixing Zhu; Date: 2026-07-12.
    mode = 'sample_linear';
    if ~isfield(rawConfig, 'soil_input') ...
            || ~isfield(rawConfig.soil_input, 'vertical_parameter_mode')
        return;
    end
    rawMode = rawConfig.soil_input.vertical_parameter_mode;
    if ~(ischar(rawMode) || (isstring(rawMode) && isscalar(rawMode)))
        error('buildSoilParameters:badVerticalMode', ...
            'soil_input.vertical_parameter_mode must be scalar text.');
    end
    mode = lower(strtrim(char(rawMode)));
    if ~ismember(mode, {'horizon_bounds', 'sample_linear'})
        error('buildSoilParameters:badVerticalMode', ...
            ['unknown soil_input.vertical_parameter_mode "%s" ' ...
             '(expected horizon_bounds|sample_linear)'], mode);
    end
end

% ---------------------------------------------------------------------------
function validateSampleDepth(depth)
% Require one unambiguous increasing coordinate for both mapping modes.
% Editor: Jizixing Zhu; Date: 2026-07-12.
    if isempty(depth) || any(~isfinite(depth)) || any(depth < 0)
        error('buildSoilParameters:badDepth', ...
            'soil.csv depth_cm must be finite, non-negative, and non-empty.');
    end
    if any(diff(depth) <= 0)
        error('buildSoilParameters:badDepthOrder', ...
            'soil.csv depth_cm must be strictly increasing with no duplicates.');
    end
end

% ---------------------------------------------------------------------------
function values = interpolateHeld(depth, samples, query)
% Linear sampled-depth mapping with flat endpoint extension.
% Editor: Jizixing Zhu; Date: 2026-07-12.
    if isscalar(depth)
        values = repmat(samples(1), size(query));
        return;
    end
    queryHeld = min(max(query(:), depth(1)), depth(end));
    values = interp1(depth, samples, queryHeld, 'linear');
end

% ---------------------------------------------------------------------------
function jmap = mapLayerToHorizon(sumdz, boundaries)
%LAYERHORIZONMAP  Step-function depth->horizon index from lower boundaries.
%   Row k is selected when depth <= boundaries(k). Depths deeper than the last
%   boundary use the last row, so a soil table with NH rows can use all NH rows.
% Editor: Jizixing Zhu; Date: 2026-06-27.
    NL = numel(sumdz);
    jmap = ones(NL, 1);
    nb = numel(boundaries);
    if nb == 0
        error('buildSoilParameters:noBoundaries', 'rawSoil.depth_cm is empty.');
    end
    for i = 1:NL
        d = sumdz(i);
        idx = nb;
        for k = 1:nb
            if d <= boundaries(k)
                idx = k;
                break;
            end
        end
        jmap(i) = idx;
    end
end
