function bd = buildBulkDensity(rawInit, soil, domain, cfg)
%BUILDBULKDENSITY  Decide the single bulk-density source up front (before C).
%
%   bd = init.buildBulkDensity(rawInit, soil, domain, cfg) picks ONE bulk density
%   and derives every representation the rest of the build needs, so soil
%   parameters, porosity, thermal params, Y.bd and the nitrogen conversion all
%   consume the same numbers (no post-hoc reconciliation).
%
%   switches.bulk_density_source:
%     'initial' (default): soil_initial_profile.csv per-depth bulk density wins. Node values
%                are interpolated onto the grid; layer values are the node mean,
%                clamped to [0.3, 2.0] g/cm3 (legacy run_stemmus [7b]).
%     'soil'   : the soil.csv per-horizon bulk density wins (soil.bulk_density);
%                node values are the per-layer value mapped to nodes.
%     Falls back to 'soil' when soil_initial_profile.csv has no usable bulk_density column.
%
%   Inputs:
%     rawInit  table   raw.initial (uses depth_cm, bulk_density_gcm3 if present)
%     soil     struct  built.parameters.soil (uses soil.bulk_density and soil.theta_s)
%     domain   struct  built.domain (z_node, NL, NN)
%     cfg      struct  raw.config (switches.bulk_density_source)
%
%   Returns bd (all bottom-up):
%     bd.node     [g/cm3, NN x 1]  node bulk density (drives Y.bd + N conversion)
%     bd.layer    [g/cm3, NL x 1]  layer bulk density (drives soil bd + thermal)
%     bd.porosity [m3/m3, NL x 1]  1 - bd.layer/2.65, clamped to [0.20, 0.85] and
%                                  raised to theta_s + 0.03 where bulk density would
%                                  otherwise leave too little dry-air headroom
%     bd.source   [char]           'initial' or 'soil'
%
%   Reference: run_stemmus_initialize.m [7b] (node->layer mean, [0.3,2.0]/
%   [0.20,0.85] clamps, porosity = 1 - bd/2.65) and InitialProfile.m bd fallback.
%   The POR >= theta_s + 0.03 self-consistency follows the legacy
%   init.setVanGenuchtenParameters (POR vs VG.Theta_s strong constraint).
% Editor: Jizixing Zhu; Date: 2026-06-27.
% Change: 2026-06-28 Jizixing Zhu — restore POR >= theta_s + 0.03 raise + init:PorosityVGConflict warning.
% Change: 2026-07-16 Jizixing Zhu -- document the descriptive soil initialization filename.

    NL = domain.NL;
    NN = domain.NN;

    source = lower(strtrim(getSwitchString(cfg, 'bulk_density_source', 'initial')));
    if ~any(strcmp(source, {'initial', 'soil'}))
        error('buildBulkDensity:badSource', ...
            'unknown switches.bulk_density_source "%s" (expected initial|soil)', source);
    end
    hasIC = istable(rawInit) ...
        && ismember('bulk_density_gcm3', rawInit.Properties.VariableNames) ...
        && any(rawInit.bulk_density_gcm3 ~= 0) ...
        && ~all(isnan(rawInit.bulk_density_gcm3));

    if strcmp(source, 'initial') && hasIC
        bd.source = 'initial';
        % node bulk density: interpolate IC profile onto node depths (top-down),
        % then flip to bottom-up.
        depths_target = flip(domain.z_node(:));                 % top-down (surface=0)
        node_td = init.interpLinearWithExtrapolation(rawInit.depth_cm(:), ...
            rawInit.bulk_density_gcm3(:), depths_target);
        bd.node = flip(node_td);                                % [NN x 1] bottom-up
        % layer bulk density = mean of bounding nodes, clamped (legacy [7b])
        layer = (bd.node(1:NL) + bd.node(2:NL+1)) / 2;
        bd.layer = clampWithWarn(layer, 0.3, 2.0, 'bulk density');
    else
        bd.source = 'soil';
        bd.layer = soil.bulk_density(:);                        % [NL x 1] soil.csv
        bd.node  = mapLayerToNodeMean(bd.layer, NN, NL);       % [NN x 1] bottom-up
    end

    % porosity from the chosen layer bulk density, clamped (legacy [7b])
    bd.porosity = clampWithWarn(1 - bd.layer / 2.65, 0.20, 0.85, 'porosity');

    % POR vs theta_s self-consistency (legacy init.setVanGenuchtenParameters):
    % geometric porosity must sit >= trapped_air_margin above the VG saturated
    % water content so Theta_V = POR - Theta_L - Theta_II stays positive. Raise POR
    % where the bulk-density value violates this and warn the user to check inputs.
    if isfield(soil, 'theta_s') && ~isempty(soil.theta_s)
        trapped_air_margin = 0.03;                          % [m3/m3] dry-air headroom
        por_required = soil.theta_s(:) + trapped_air_margin;
        n_bad = sum(bd.porosity(:) < por_required);
        if n_bad > 0
            max_gap = max(por_required - bd.porosity(:));
            warning('init:PorosityVGConflict', ...
                ['porosity (from bulk density) < theta_s + %.2f in %d/%d layer(s), ', ...
                 'max gap %.4f m3/m3. Raising POR = max(POR, theta_s + %.2f); ', ...
                 'check soil.csv bulk_density vs theta_s.'], ...
                trapped_air_margin, n_bad, numel(bd.porosity), max_gap, trapped_air_margin);
            bd.porosity = max(bd.porosity(:), por_required);
        end
    end
end
