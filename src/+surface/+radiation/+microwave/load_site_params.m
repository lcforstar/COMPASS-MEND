function site = load_site_params(site_dir, depth_node_cm)
%LOAD_SITE_PARAMS  Load canonical V2 site inputs for the microwave operator.
%
%   site = surface.radiation.microwave.load_site_params(site_dir, depth_node_cm)
%   reads one self-contained V2 site package. It maps the primitive texture
%   samples in soil.csv to the accepted-output node coordinates using the
%   site-selected sample-linear contract with held endpoints.
%
%   Inputs:
%     site_dir       [char|string]  Site-package directory containing site.toml
%                                    and soil.csv (for example data/sites/changwu)
%     depth_node_cm  [cm, NN x 1]   Output-node depth below surface, top-down
%
%   Returns:
%     site.sand_prof [--, NN x 1]      Sand fraction at output nodes
%     site.clay_prof [--, NN x 1]      Clay fraction at output nodes
%     site.bd_prof   [g cm-3, NN x 1]  Active bulk density at output nodes
%     site.msoc_prof [--, NN x 1]      Organic-matter mass fraction at output nodes
%     site.depth_cm  [cm, NN x 1]      Returned top-down node coordinates
%
%   The active bulk-density source follows switches.bulk_density_source. The
%   default initial-profile path uses init.interpLinearWithExtrapolation, the
%   same node mapping used by init.buildBulkDensity. Microwave diagnostics are
%   post-processing only: this loader never modifies model state or parameters.
%
% Editor: Jizixing Zhu; Date: 2026-07-21.
% Change: 2026-07-21 Jizixing Zhu — replace V1 config/XLSX reads with canonical V2 site inputs.

    site_dir = char(site_dir);
    if ~isfolder(site_dir)
        error('load_site_params:SiteDirectory', ...
            'Missing V2 site directory: %s', site_dir);
    end
    depth_node_cm = depth_node_cm(:);
    if isempty(depth_node_cm) || any(~isfinite(depth_node_cm)) || ...
            any(depth_node_cm < 0) || any(diff(depth_node_cm) < 0)
        error('load_site_params:DepthNodes', ...
            'depth_node_cm must be finite, non-negative, and top-down ordered.');
    end

    config = io.loadSiteConfig(site_dir);
    soil_path = fullfile(site_dir, 'soil.csv');
    if ~isfile(soil_path)
        error('load_site_params:SoilInput', 'Missing canonical soil input: %s', soil_path);
    end
    soil = io.adapters.readCsvTable(soil_path);
    requireColumns(soil, {'depth_cm', 'sand_pct', 'clay_pct', ...
        'bulk_density_gcm3', 'msoc_mgkg'}, soil_path);
    validateSampleDepth(soil.depth_cm);

    site.sand_prof = interpolateHeld(soil.depth_cm, soil.sand_pct / 100, depth_node_cm);
    site.clay_prof = interpolateHeld(soil.depth_cm, soil.clay_pct / 100, depth_node_cm);
    site.msoc_prof = interpolateHeld(soil.depth_cm, soil.msoc_mgkg / 1e6, depth_node_cm);
    site.bd_prof = readActiveBulkDensity(site_dir, config, soil, depth_node_cm);
    site.depth_cm = depth_node_cm;
    site.sample_depth_cm = soil.depth_cm(:);

    site.Z_km = requiredScalar(config.site, 'elevation', 'site.elevation') / 1000;
    site.Latitude = requiredScalar(config.site, 'latitude', 'site.latitude');
    site.Longitude = requiredScalar(config.site, 'longitude', 'site.longitude');
    site.hc = requiredScalar(config.vegetation, 'canopy_height', ...
        'vegetation.canopy_height');
    site.Cw = requiredScalar(config.leaf_optical, 'Cw', 'leaf_optical.Cw');
    site.Cdm = requiredScalar(config.leaf_optical, 'Cdm', 'leaf_optical.Cdm');
    site.N_leaf = requiredScalar(config.leaf_optical, 'N_mesophyll', ...
        'leaf_optical.N_mesophyll');
    site.LIDFa = requiredScalar(config.vegetation, 'LIDFa', 'vegetation.LIDFa');
    site.LIDFb = requiredScalar(config.vegetation, 'LIDFb', 'vegetation.LIDFb');
    site.m_d = site.Cdm / max(site.Cdm + site.Cw, 1e-6);  % [--]

    if site.LIDFa <= -0.5 && site.LIDFb <= -0.5
        site.tt_H_prior = 1.0;
        site.tt_V_prior = 0.0;
        site.lidf_label = 'planophile';
    elseif site.LIDFa >= 0.5
        site.tt_H_prior = 0.0;
        site.tt_V_prior = 2.0;
        site.lidf_label = 'erectophile';
    else
        site.tt_H_prior = 1.0;
        site.tt_V_prior = 1.0;
        site.lidf_label = 'spherical';
    end
end

% -------------------------------------------------------------------------
function bd = readActiveBulkDensity(site_dir, config, soil, depth_node_cm)
%READACTIVEBULKDENSITY  Reproduce the V2 node-level bulk-density source.
% Editor: Jizixing Zhu; Date: 2026-07-21.
    source = 'initial';
    if isfield(config, 'switches') && isfield(config.switches, 'bulk_density_source')
        source = lower(strtrim(char(config.switches.bulk_density_source)));
    end
    if ~ismember(source, {'initial', 'soil'})
        error('load_site_params:BulkDensitySource', ...
            'Unsupported switches.bulk_density_source: %s', source);
    end

    initial_path = fullfile(site_dir, 'soil_initial_profile.csv');
    if strcmp(source, 'initial') && isfile(initial_path)
        initial = io.adapters.readCsvTable(initial_path);
        has_bd = ismember('bulk_density_gcm3', initial.Properties.VariableNames) ...
            && any(isfinite(initial.bulk_density_gcm3)) ...
            && any(initial.bulk_density_gcm3 ~= 0);
        if has_bd
            bd = init.interpLinearWithExtrapolation(initial.depth_cm, ...
                initial.bulk_density_gcm3, depth_node_cm);
            return;
        end
    end

    bd = interpolateHeld(soil.depth_cm, soil.bulk_density_gcm3, depth_node_cm);
end

% -------------------------------------------------------------------------
function values = interpolateHeld(sample_depth_cm, samples, query_depth_cm)
%INTERPOLATEHELD  Linear depth interpolation with constant endpoint extension.
% Editor: Jizixing Zhu; Date: 2026-07-21.
    sample_depth_cm = sample_depth_cm(:);
    samples = samples(:);
    if numel(sample_depth_cm) == 1
        values = repmat(samples(1), size(query_depth_cm));
        return;
    end
    query_held = min(max(query_depth_cm(:), sample_depth_cm(1)), sample_depth_cm(end));
    values = interp1(sample_depth_cm, samples, query_held, 'linear');
end

% -------------------------------------------------------------------------
function validateSampleDepth(depth_cm)
%VALIDATESAMPLEDEPTH  Require one increasing canonical soil-coordinate axis.
% Editor: Jizixing Zhu; Date: 2026-07-21.
    depth_cm = depth_cm(:);
    if isempty(depth_cm) || any(~isfinite(depth_cm)) || any(depth_cm < 0) ...
            || any(diff(depth_cm) <= 0)
        error('load_site_params:SoilDepth', ...
            'soil.csv depth_cm must be finite, non-negative, and strictly increasing.');
    end
end

% -------------------------------------------------------------------------
function requireColumns(tbl, names, path)
%REQUIRECOLUMNS  Fail fast if a V2 microwave input lacks a required column.
% Editor: Jizixing Zhu; Date: 2026-07-21.
    missing = names(~ismember(names, tbl.Properties.VariableNames));
    if ~isempty(missing)
        error('load_site_params:SoilColumns', ...
            'Missing columns in %s: %s', path, strjoin(missing, ', '));
    end
end

% -------------------------------------------------------------------------
function value = requiredScalar(group, field, label)
%REQUIREDSCALAR  Read one finite scalar from the canonical site configuration.
% Editor: Jizixing Zhu; Date: 2026-07-21.
    if ~isstruct(group) || ~isfield(group, field)
        error('load_site_params:SiteConfig', 'Missing required %s.', label);
    end
    value = group.(field);
    if ~(isnumeric(value) && isscalar(value) && isfinite(value))
        error('load_site_params:SiteConfig', '%s must be a finite scalar.', label);
    end
end
