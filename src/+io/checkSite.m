function report = checkSite(raw, plan)
%CHECKSITE  Stage-1.5 validator: check raw containers against the column catalogue.
%
%   report = io.checkSite(raw) verifies the legacy full raw site. The optional
%   report = io.checkSite(raw, plan) form makes only plan-enabled inputs
%   mandatory and skips validators for inactive components.
%   BEFORE any grid construction or interpolation (io.buildState). It is the
%   "check" step between assembly and build. Validation is BY NAME against
%   io.adapters.tableSpec, so reordered or newly added columns never break the
%   read: a required column that is absent is a hard error; a column not in the
%   catalogue is reported as "extra" and left untouched (e.g. a user-added
%   'soil_type' column flows through to buildState unchanged).
%
%   Inputs:
%     raw   struct  output of io.loadSite (fields: config, soil, initial,
%                   plant, root, forcing, water_boundary, ...). Optional tables
%                   may be [] if the file was absent.
%
%   Returns:
%     report  struct  per-table outcome with fields:
%                       report.<kind>.present   [logical]  table was read
%                       report.<kind>.extra     {cellstr}  columns/keys not in catalogue
%                       report.<kind>.missing   {cellstr}  required but absent (-> also errors)
%                       report.ok               [logical]  true if no hard errors
%
%   Errors (hard): a required column/pool is missing, or a mandatory file
%   (soil, initial) was not read at all. Nitrogen columns are required only when
%   raw.config.switches.solve_nitrogen is true.
%
%   No physics, no unit conversion, no interpolation -- those live in buildState.
%
%   Editor: Jizixing Zhu; Date: 2026-06-26.
%   Change: 2026-07-11 Jizixing Zhu -- validate optional water_boundary.csv against typed series sources.
%   Change: 2026-07-11 Jizixing Zhu -- require an active prescribed series to cover the full configured simulation interval.
%   Change: 2026-07-11 Jizixing Zhu -- report the resolved configured water-boundary path in validation errors.
%   Change: 2026-07-11 Jizixing Zhu -- derive fixed top_value/bottom_value CSV requirements from variable_head/variable_flux modes.
%   Change: 2026-07-11 Jizixing Zhu -- require CSV for all prescribed modes and enforce constant columns for constant_head/constant_flux.
%   Change: 2026-07-11 Jizixing Zhu -- validate only the inputs and components selected by the initialization plan.
%   Change: 2026-07-12 Jizixing Zhu -- validate soil.csv vertical mapping mode and sample coordinates.
%   Change: 2026-07-16 Jizixing Zhu -- report the descriptive soil and plant initialization filenames.
%   Change: 2026-07-22 Jizixing Zhu -- reject the retired [advanced] and [debug] sections instead of silently ignoring them.
%   Change: 2026-08-29 Codex -- validate measured SOC, native pool fractions,
%     and ordered soil.csv calibration priors when soil carbon is active.

    if nargin < 2 || isempty(plan)
        plan = fullValidationPlan();
    end
    report = struct('ok', true);
    problems = {};

    % --- which optional groups become mandatory, from config switches ---
    nitrogenOn = getSwitch(raw, 'solve_nitrogen', false);

    % --- wide tables: soil (mandatory), initial (mandatory), root (optional file) ---
    [report.soil,    problems] = checkWide(raw, 'soil',    'soil',    inputEnabled(plan, 'soil'), problems);
    [report.initial, problems] = checkWide(raw, 'initial', 'initial', inputEnabled(plan, 'initial'), problems);
    [report.root,    problems] = checkWide(raw, 'root',    'root',    inputEnabled(plan, 'root'), problems);
    [report.forcing, problems] = checkWide(raw, 'forcing', 'forcing', inputEnabled(plan, 'forcing'), problems);
    [report.water_boundary, problems] = checkWide(raw, 'water_boundary', ...
        'water_boundary', false, problems);
    organicOn = getSoilCarbonMaster(raw);
    if inputEnabled(plan, 'organic_initial')
        [report.organic_initial, problems] = checkOrganicInitial(raw, problems);
    else
        report.organic_initial = inactiveReport();
    end
    if organicOn
        [report.organic_config, problems] = checkOrganicConfiguration( ...
            raw, problems);
    else
        report.organic_config = inactiveReport();
    end

    if report.soil.present && isfield(raw, 'config') && isfield(raw.config, 'grid') ...
            && isfield(raw.config.grid, 'NS')
        expectedSoilRows = raw.config.grid.NS;
        if height(raw.soil) ~= expectedSoilRows
            report.soil.missing = [report.soil.missing, {'row_count_matches_grid_NS'}];
            report.ok = false;
            problems{end+1} = sprintf('soil.csv: row count (%d) must match grid.NS (%d)', ...
                height(raw.soil), expectedSoilRows);
        end
    end
    if report.soil.present && ismember('depth_cm', raw.soil.Properties.VariableNames)
        soilDepth = raw.soil.depth_cm;
        if ~isnumeric(soilDepth) || any(~isfinite(soilDepth)) ...
                || any(soilDepth < 0) || any(diff(soilDepth) <= 0)
            report.soil.missing = [report.soil.missing, {'strictly_increasing_depth_cm'}];
            report.ok = false;
            problems{end+1} = ...
                'soil.csv: depth_cm must be finite, non-negative, and strictly increasing';
        end
    end

    % soil and initial nitrogen columns are required only when the N cycle is on
    if nitrogenOn && report.soil.present
        nCols = {'difnh4_cm2d', 'difno3_cm2d', 'kg_cm3g', 'vmax_nit_ugn_cm3d', ...
                 'km_nit_ugn_cm3', 'kmineralization_1d'};
        miss = nCols(~ismember(nCols, raw.soil.Properties.VariableNames));
        if ~isempty(miss)
            report.soil.missing = [report.soil.missing, miss];
            report.ok = false;
            problems{end+1} = sprintf('soil.csv: nitrogen is ON but missing N columns: %s', strjoin(miss, ', '));
        end
    end
    if nitrogenOn && report.initial.present
        nInitCols = {'nh4', 'no3'};
        miss = nInitCols(~ismember(nInitCols, raw.initial.Properties.VariableNames));
        if ~isempty(miss)
            report.initial.missing = [report.initial.missing, miss];
            report.ok = false;
            problems{end+1} = sprintf('soil_initial_profile.csv: nitrogen is ON but missing initial N column(s): %s', ...
                strjoin(miss, ', '));
        end
    end

    % --- long table: plant pools ---
    if inputEnabled(plan, 'plant')
        [report.plant, problems] = checkPlant(raw, problems);
    else
        report.plant = inactiveReport();
    end

    % --- required config presence (minimal smoke check, not full schema) ---
    problems = checkConfigKeys(raw, problems);
    problems = checkWaterBoundarySeries(raw, report.water_boundary, problems);

    % --- SCOPE/static-library inputs (paths, NetCDF variables, light ranges) ---
    if inputEnabled(plan, 'scope')
        [report.scope, problems] = checkScopeInputs(raw, problems);
    else
        report.scope = inactiveReport();
    end

    % --- raise a single, aggregated error if anything is wrong ---
    if ~isempty(problems)
        error('checkSite:invalid', 'Site check failed:\n  - %s', strjoin(problems, '\n  - '));
    end
end

% ---------------------------------------------------------------------------
function value = getSoilCarbonMaster(raw)
% Read the sole supported soil-carbon master.
% Editor: Codex; Date: 2026-08-29.
    value = getSwitch(raw, 'soil_carbon_model', false);
end

% ---------------------------------------------------------------------------
function plan = fullValidationPlan()
% Backward-compatible validation policy for direct checkSite callers.
% Editor: Jizixing Zhu; Date: 2026-07-11.
    plan.inputs = struct('soil', true, 'initial', true, 'forcing', true, ...
        'plant', true, 'root', false, 'water_boundary', true, ...
        'vegetation', true, 'scope', true, 'organic_initial', false);
end

% ---------------------------------------------------------------------------
function [r, problems] = checkOrganicInitial(raw, problems)
% Validate native POC/MOC/MB/IOM initialization fractions.
% Editor: Codex; Date: 2026-08-29.
    r = struct('present', false, 'extra', {{}}, 'missing', {{}});
    tbl = getfieldOr(raw, 'organic_initial', []);
    required = {'depth_cm','frac_POC','frac_MOC','frac_MB','frac_IOM'};
    if ~istable(tbl) || isempty(tbl)
        r.missing = required;
        problems{end+1} = ...
            'organic_initial_pools.csv: required table missing or empty';
        return;
    end
    r.present = true;
    columns = tbl.Properties.VariableNames;
    r.missing = required(~ismember(required, columns));
    r.extra = columns(~ismember(columns, required));
    if ~isempty(r.missing)
        problems{end+1} = sprintf( ...
            'organic_initial_pools.csv: missing required column(s): %s', ...
            strjoin(r.missing, ', '));
        return;
    end
    values = tbl{:, required};
    fractions = values(:, 2:end);
    if any(~isfinite(values), 'all') || any(values(:,1) < 0) ...
            || any(diff(values(:,1)) <= 0)
        problems{end+1} = ['organic_initial_pools.csv: depth_cm must be ', ...
            'finite, nonnegative, and strictly increasing'];
    end
    if any(fractions < 0, 'all') || any(fractions > 1, 'all') ...
            || any(abs(sum(fractions, 2)-1) > 1e-10)
        problems{end+1} = ['organic_initial_pools.csv: fractions must be ', ...
            'within [0,1] and each row must sum to one'];
    end
end

% ---------------------------------------------------------------------------
function [r, problems] = checkOrganicConfiguration(raw, problems)
% Validate native SOC ownership and site-global organic parameter sections.
% Editor: Codex; Date: 2026-08-29.
% Change: 2026-09-02 Codex -- validate site-global parameters in site.toml;
%   soil.csv is reserved for quantities that may vary with depth.
    r = struct('present', true, 'extra', {{}}, 'missing', {{}});
    if ~isfield(raw, 'initial') || ~istable(raw.initial) ...
            || ~ismember('soc_total_mgkg', raw.initial.Properties.VariableNames)
        r.missing{end+1} = 'soil_initial_profile.csv:soc_total_mgkg';
        problems{end+1} = ['soil carbon model requires measured initial ', ...
            'soc_total_mgkg in soil_initial_profile.csv'];
    else
        soc = raw.initial.soc_total_mgkg;
        if ~isnumeric(soc) || any(~isfinite(soc)) || any(soc <= 0)
            problems{end+1} = ['soil_initial_profile.csv:soc_total_mgkg ', ...
                'must contain finite positive values'];
        end
    end
    sections={ ...
        'organic',{'k_litter_1d','k_straw_1d','k_manure_1d','doc_leach_1d'}; ...
        'dom',{'diffusion_cm2_s','transport_substep_s'}; ...
        'microbe',{'reaction_substep_s','active_fraction_initial', ...
            'mend_lci_initial','mend_vd_pom_h','mend_ks_pom1_mg_g', ...
            'mend_ks_pom_ratio','mend_enzyme_turnover_h','mend_penzp', ...
            'mend_penzm_ratio','mend_frac_pom_to_dom', ...
            'mend_frac_mortality_to_dom','mend_vg_h','mend_alpha', ...
            'mend_ks_dom_mg_g','mend_cue_20','mend_cue_slope_per_c', ...
            'mend_mortality_multiplier','mend_dormant_maintenance_ratio', ...
            'mend_swp_a2d_mpa','mend_tau','mend_wdorm','mend_vnup_mb_h', ...
            'mend_ks_nh4_mb_mg_g','mend_ks_no3_mb_mg_g', ...
            'mend_cn_mb_min','mend_cn_mb_max','mend_cn_pom2', ...
            'mend_cn_enzp','mend_cn_enzm','epsilon_n_mg_g'}; ...
        'bioparticle',{'mobile_fraction_initial','diffusion_cm2_s', ...
            'transport_substep_s'}};
    for i=1:size(sections,1)
        sectionName=sections{i,1}; keys=sections{i,2};
        if ~isfield(raw.config,sectionName)||~isstruct(raw.config.(sectionName))
            r.missing{end+1}=['[' sectionName ']'];
            problems{end+1}=sprintf('site.toml: [%s] is required',sectionName); %#ok<AGROW>
            continue;
        end
        section=raw.config.(sectionName);
        for j=1:numel(keys)
            key=keys{j}; label=[sectionName '.' key];
            if ~isfield(section,key)
                r.missing{end+1}=label;
                problems{end+1}=sprintf('site.toml: %s is required',label); %#ok<AGROW>
            elseif ~isnumeric(section.(key))||~isscalar(section.(key)) ...
                    ||~isfinite(section.(key))||section.(key)<0
                problems{end+1}=sprintf( ...
                    'site.toml: %s must be a finite nonnegative scalar',label); %#ok<AGROW>
            end
        end
    end
    positive={ ...
        'organic','k_litter_1d';'organic','k_straw_1d'; ...
        'organic','k_manure_1d';'microbe','reaction_substep_s'; ...
        'microbe','mend_vd_pom_h';'microbe','mend_ks_pom1_mg_g'; ...
        'microbe','mend_ks_pom_ratio';'microbe','mend_vg_h'; ...
        'microbe','mend_ks_dom_mg_g';'microbe','mend_cue_20'; ...
        'microbe','mend_cn_mb_min';'microbe','mend_cn_mb_max'; ...
        'microbe','mend_cn_pom2';'microbe','mend_cn_enzp'; ...
        'microbe','mend_cn_enzm';'microbe','epsilon_n_mg_g'};
    for i=1:size(positive,1)
        section=positive{i,1}; key=positive{i,2};
        if isfield(raw.config,section)&&isstruct(raw.config.(section)) ...
                &&isfield(raw.config.(section),key) ...
                &&isnumeric(raw.config.(section).(key)) ...
                &&isscalar(raw.config.(section).(key)) ...
                &&isfinite(raw.config.(section).(key)) ...
                &&raw.config.(section).(key)<=0
            problems{end+1}=sprintf('site.toml: %s.%s must be positive', ...
                section,key); %#ok<AGROW>
        end
    end
    soilProfiles={'kd_doc_cm3_g','kd_don_cm3_g','dom_disp_cm', ...
        'attachment_1d','detachment_1d','straining_1d','microbe_disp_cm'};
    if ~isfield(raw,'soil')||~istable(raw.soil)
        missing=soilProfiles;
    else
        missing=soilProfiles(~ismember(soilProfiles, ...
            raw.soil.Properties.VariableNames));
    end
    if ~isempty(missing)
        r.missing=[r.missing,missing];
        problems{end+1}=sprintf( ...
            'soil.csv: missing depth-resolved organic parameter columns: %s', ...
            strjoin(missing,', '));
    else
        for i=1:numel(soilProfiles)
            value=raw.soil.(soilProfiles{i});
            if any(~isfinite(value))||any(value<0)
                problems{end+1}=sprintf( ...
                    'soil.csv:%s must be finite and nonnegative', ...
                    soilProfiles{i}); %#ok<AGROW>
            end
        end
        if any(raw.soil.kd_doc_cm3_g<=0)||any(raw.soil.kd_don_cm3_g<=0)
            problems{end+1}=['soil.csv: kd_doc_cm3_g and kd_don_cm3_g ', ...
                'must be positive at every sampled depth'];
        end
    end
end

% ---------------------------------------------------------------------------
function tf = inputEnabled(plan, name)
% Return whether one input belongs to the active component plan.
% Editor: Jizixing Zhu; Date: 2026-07-11.
    tf = true;
    if isfield(plan, 'inputs') && isfield(plan.inputs, name)
        tf = logical(plan.inputs.(name));
    end
end

% ---------------------------------------------------------------------------
function r = inactiveReport()
% Stable report shape for a component intentionally absent by configuration.
% Editor: Jizixing Zhu; Date: 2026-07-11.
    r = struct('present', false, 'extra', {{}}, 'missing', {{}}, ...
        'inactive_by_configuration', true);
end

% ---------------------------------------------------------------------------
function [r, problems] = checkWide(raw, field, kind, mandatory, problems)
% Check one wide raw table against its tableSpec column catalogue.
% Editor: Jizixing Zhu; Date: 2026-06-26.
    spec = io.adapters.tableSpec(kind);
    r = struct('present', false, 'extra', {{}}, 'missing', {{}});
    tbl = getfieldOr(raw, field, []);
    if ~istable(tbl) || isempty(tbl)
        if mandatory
            problems{end+1} = sprintf('%s.csv: required table missing or empty', kind);
        end
        return;
    end
    r.present = true;
    cols = tbl.Properties.VariableNames;
    r.missing = spec.required(~ismember(spec.required, cols));
    known     = [spec.required, spec.optional];
    r.extra   = cols(~ismember(cols, known));   % user-added columns: allowed, carried through
    if ~isempty(r.missing)
        problems{end+1} = sprintf('%s.csv: missing required column(s): %s', ...
            kind, strjoin(r.missing, ', '));
    end
end

% ---------------------------------------------------------------------------
function [r, problems] = checkPlant(raw, problems)
% Check the long-form initial plant pool table.
% Editor: Jizixing Zhu; Date: 2026-06-26.
    spec = io.adapters.tableSpec('plant');
    r = struct('present', false, 'extra', {{}}, 'missing', {{}});
    tbl = getfieldOr(raw, 'plant', []);
    if ~istable(tbl) || isempty(tbl)
        problems{end+1} = 'plant_initial_pools.csv: required table missing or empty';
        return;
    end
    if ~ismember(spec.keyCol, tbl.Properties.VariableNames)
        problems{end+1} = sprintf('plant_initial_pools.csv: missing key column "%s"', spec.keyCol);
        return;
    end
    if ~ismember(spec.valCol, tbl.Properties.VariableNames)
        problems{end+1} = sprintf('plant_initial_pools.csv: missing value column "%s"', spec.valCol);
        return;
    end
    r.present = true;
    pools = cellstr(string(tbl.(spec.keyCol)));
    r.missing = spec.required(~ismember(spec.required, pools));
    r.extra   = pools(~ismember(pools, spec.required))';
    if ~isempty(r.missing)
        problems{end+1} = sprintf('plant_initial_pools.csv: missing required pool(s): %s', ...
            strjoin(r.missing, ', '));
    end
    for i = 1:numel(spec.required)
        pool = spec.required{i};
        row = strcmp(pools, pool);
        val = tbl.(spec.valCol)(row);
        if isempty(val) || ~isnumeric(val) || any(~isfinite(val))
            problems{end+1} = sprintf('plant_initial_pools.csv: pool "%s" has missing or non-finite value', pool); %#ok<AGROW>
        end
    end
end

% ---------------------------------------------------------------------------
function problems = checkConfigKeys(raw, problems)
% Check minimal required site.toml sections and enum strings.
% Editor: Jizixing Zhu; Date: 2026-06-26.
    if ~isfield(raw, 'config') || ~isstruct(raw.config)
        problems{end+1} = 'site.toml: config missing';
        return;
    end
    needSections = {'site', 'time', 'grid', 'switches'};
    miss = needSections(~isfield(raw.config, needSections));
    if ~isempty(miss)
        problems{end+1} = sprintf('site.toml: missing section(s): [%s]', strjoin(miss, '], ['));
    end

    problems = checkUnsupportedConfigSection(raw.config, 'advanced', ...
        'use the typed [switches] process controls instead', problems);
    problems = checkUnsupportedConfigSection(raw.config, 'debug', ...
        'use explicit runner diagnostics rather than scientific site configuration', problems);

    % enum switches: catch typos that would otherwise silently fall back to a default
    if isfield(raw.config, 'switches')
        problems = checkEnum(raw.config.switches, 'nitrogen_input_phase', ...
            {'total', 'dissolved'}, problems);
        problems = checkEnum(raw.config.switches, 'bulk_density_source', ...
            {'initial', 'soil'}, problems);
        problems = checkEnum(raw.config.switches, 'diurnal_temp_method', ...
            {'legacy_sine', 'sunrise_anchored'}, problems);
    end
    if isfield(raw.config, 'time')
        problems = checkEnum(raw.config.time, 'clock_basis', {'local_standard', 'utc'}, problems);
    end
    if isfield(raw.config, 'soil_input')
        problems = checkEnum(raw.config.soil_input, 'vertical_parameter_mode', ...
            {'sample_linear', 'horizon_bounds'}, problems);
    end
end

% ---------------------------------------------------------------------------
function problems = checkUnsupportedConfigSection(config, sectionName, guidance, problems)
% Reject a retired section before any of its keys can be mistaken for controls.
% Editor: Jizixing Zhu; Date: 2026-07-22.
    if ~isfield(config, sectionName)
        return;
    end
    problems{end+1} = sprintf( ...
        'site.toml: unsupported configuration section [%s]; %s', ...
        sectionName, guidance);
end

% ---------------------------------------------------------------------------
function problems = checkWaterBoundarySeries(raw, report, problems)
% Require and validate the two fixed external value columns selected by modes.
% Editor: Jizixing Zhu; Date: 2026-07-11.
% Change: 2026-07-11 Jizixing Zhu -- read the single flat water_top_*/water_bottom_* site schema.
    if ~isfield(raw, 'config') || ~isstruct(raw.config) ...
            || ~isfield(raw.config, 'boundary') ...
            || ~isstruct(raw.config.boundary)
        return;
    end

    boundary = raw.config.boundary;
    fileLabel = waterBoundaryFileLabel(raw);
    positions = {'top','bottom'};
    requiredColumns = {};
    constantColumns = {};
    for i = 1:numel(positions)
        position = positions{i};
        kindField = ['water_' position '_kind'];
        if isfield(boundary, kindField)
            mode = lower(strtrim(char(string(boundary.(kindField)))));
            if ismember(mode, {'constant_head','variable_head', ...
                    'constant_flux','variable_flux'})
                requiredColumns{end+1} = [position '_value']; %#ok<AGROW>
            end
            if ismember(mode, {'constant_head','constant_flux'})
                constantColumns{end+1} = [position '_value']; %#ok<AGROW>
            end
        end
    end
    requiredColumns = unique(requiredColumns, 'stable');
    if isempty(requiredColumns)
        return;
    end

    if ~report.present || ~isfield(raw, 'water_boundary') ...
            || ~istable(raw.water_boundary) || isempty(raw.water_boundary)
        problems{end+1} = sprintf( ...
            '%s: required by a prescribed head or flux water boundary', ...
            fileLabel);
        return;
    end

    tableData = raw.water_boundary;
    names = tableData.Properties.VariableNames;
    missing = requiredColumns(~ismember(requiredColumns, names));
    if ~isempty(missing)
        problems{end+1} = sprintf( ...
            '%s: missing configured series column(s): %s', ...
            fileLabel, strjoin(missing, ', '));
    end

    if ismember('time_s', names)
        time_s = tableData.time_s;
        timeValid = isnumeric(time_s) && isvector(time_s) && numel(time_s) >= 2 ...
            && all(isfinite(time_s)) && all(diff(double(time_s(:))) > 0);
        if ~timeValid
            problems{end+1} = sprintf( ...
                '%s: time_s must be finite and strictly increasing with at least two rows', ...
                fileLabel);
        else
            duration_s = simulationDurationSeconds(raw.config);
            time_s = double(time_s(:));
            if isfinite(duration_s) ...
                    && (time_s(1) > 0 || time_s(end) < duration_s)
                problems{end+1} = sprintf( ...
                    ['%s: time_s coverage [%.15g, %.15g] s ' ...
                     'must include the simulation interval [0, %.15g] s'], ...
                    fileLabel, time_s(1), time_s(end), duration_s);
            end
        end
    end

    presentColumns = requiredColumns(ismember(requiredColumns, names));
    for i = 1:numel(presentColumns)
        name = presentColumns{i};
        values = tableData.(name);
        if ~isnumeric(values) || ~isvector(values) ...
                || numel(values) ~= height(tableData) || any(~isfinite(values))
            problems{end+1} = sprintf( ...
                '%s: column %s must be a finite numeric vector', ...
                fileLabel, name); %#ok<AGROW>
        elseif ismember(name, constantColumns) && ~isConstantColumn(values)
            problems{end+1} = sprintf( ...
                '%s: constant boundary mode requires column %s to remain unchanged', ...
                fileLabel, name); %#ok<AGROW>
        end
    end
end

% ---------------------------------------------------------------------------
function tf = isConstantColumn(values)
% Compare one numeric CSV column against its first value at roundoff tolerance.
% Editor: Jizixing Zhu; Date: 2026-07-11.
    values = double(values(:));
    reference = values(1);
    tolerance = 8 * eps(max(1, abs(reference)));
    tf = all(abs(values - reference) <= tolerance);
end

% ---------------------------------------------------------------------------
function label = waterBoundaryFileLabel(raw)
% Prefer the resolved configured path while remaining compatible with raw fixtures.
% Editor: Jizixing Zhu; Date: 2026-07-11.
    label = 'water_boundary.csv';
    if isfield(raw, 'water_boundary_path') ...
            && (ischar(raw.water_boundary_path) ...
            || (isstring(raw.water_boundary_path) && isscalar(raw.water_boundary_path)))
        candidate = strtrim(char(raw.water_boundary_path));
        if ~isempty(candidate)
            label = candidate;
        end
    end
end

% ---------------------------------------------------------------------------
function duration_s = simulationDurationSeconds(config)
% Return the configured simulation span in seconds, or NaN if unavailable.
% Editor: Jizixing Zhu; Date: 2026-07-11.
    duration_s = NaN;
    if ~isfield(config, 'time') || ~isstruct(config.time) ...
            || ~all(isfield(config.time, {'start','end'}))
        return;
    end
    try
        format = 'yyyy-MM-dd''T''HH:mm:ss';
        startTime = datetime(char(string(config.time.start)), ...
            'InputFormat', format);
        endTime = datetime(char(string(config.time.end)), ...
            'InputFormat', format);
        value = seconds(endTime - startTime);
        if isfinite(value) && value >= 0
            duration_s = double(value);
        end
    catch
        % The general time-schema validator owns malformed clock strings.
    end
end

% ---------------------------------------------------------------------------
function [r, problems] = checkScopeInputs(raw, problems)
% Validate SCOPE-related config and NetCDF static-library contracts.
% Editor: Jizixing Zhu; Date: 2026-06-26.
    r = struct('present', false, 'missing', {{}});
    if ~isfield(raw, 'config') || ~isstruct(raw.config)
        r.missing = {'config'};
        return;
    end

    cfg = raw.config;
    need = {'scope', 'leaf_optical', 'soil_optical', 'geometry', 'model_library'};
    miss = need(~isfield(cfg, need));
    if ~isempty(miss)
        r.missing = miss;
        problems{end+1} = sprintf('site.toml: missing SCOPE section(s): [%s]', strjoin(miss, '], ['));
        return;
    end
    r.present = true;

    problems = checkScopeEnums(cfg.scope, problems);
    soilspectrum = getNumericOrDefault(cfg.scope, 'soilspectrum', 0);
    problems = checkLeafOpticalRanges(cfg.leaf_optical, problems);
    problems = checkSoilOpticalRanges(soilspectrum, cfg.soil_optical, problems);
    problems = checkGeometryRanges(cfg.geometry, problems);
    problems = checkModelLibrary(raw, cfg, soilspectrum, problems);
end

% ---------------------------------------------------------------------------
function problems = checkScopeEnums(scope, problems)
% Check SCOPE enum-like options before builders consume them.
% Editor: Jizixing Zhu; Date: 2026-06-26.
    problems = checkIntegerEnum(scope, 'simulation', [0 1 2], 'scope', problems);
    problems = checkIntegerEnum(scope, 'soilspectrum', [0 1], 'scope', problems);
    problems = checkIntegerEnum(scope, 'rt_thermal', [0 1], 'scope', problems);
    problems = checkIntegerEnum(scope, 'fluorescence_model', [0 1], 'scope', problems);
    problems = checkIntegerEnum(scope, 'soil_heat_method', [0 1 2], 'scope', problems);
end

% ---------------------------------------------------------------------------
function problems = checkLeafOpticalRanges(leaf, problems)
% Light numeric range checks for leaf optical baselines.
% Editor: Jizixing Zhu; Date: 2026-06-26.
    nonneg = {'Cab','Cca','Cw','Cdm','Cant','fluor_kNPQs'};
    for i = 1:numel(nonneg)
        problems = checkFiniteLower(leaf, nonneg{i}, 0, 'leaf_optical', problems);
    end
    fractions = {'Cs','fqe','V2Z','fluor_beta','fluor_qLs','fluor_stressfactor'};
    for i = 1:numel(fractions)
        problems = checkFiniteRange(leaf, fractions{i}, 0, 1, 'leaf_optical', problems);
    end
    problems = checkFiniteLower(leaf, 'N_mesophyll', eps, 'leaf_optical', problems);
    problems = checkFiniteRange(leaf, 'leaf_rho_thermal', 0, 1, 'leaf_optical', problems);
    problems = checkFiniteRange(leaf, 'leaf_tau_thermal', 0, 1, 'leaf_optical', problems);
    if isfield(leaf, 'leaf_rho_thermal') && isfield(leaf, 'leaf_tau_thermal') ...
            && isfinite(leaf.leaf_rho_thermal) && isfinite(leaf.leaf_tau_thermal) ...
            && leaf.leaf_rho_thermal + leaf.leaf_tau_thermal >= 1
        problems{end+1} = 'leaf_optical: leaf_rho_thermal + leaf_tau_thermal must be < 1';
    end
end

% ---------------------------------------------------------------------------
function problems = checkSoilOpticalRanges(soilspectrum, soil, problems)
% Check soil optical baselines and BSM/file-mode requirements.
% Editor: Jizixing Zhu; Date: 2026-06-26.
    for nm = {'rss','rbs','rs_thermal'}
        problems = checkFiniteLower(soil, nm{1}, 0, 'soil_optical', problems);
    end
    if soilspectrum == 1
        problems = checkFiniteRange(soil, 'BSM_brightness', 0, 0.9, 'soil_optical', problems);
        problems = checkFiniteRange(soil, 'BSM_lat', 20, 40, 'soil_optical', problems);
        problems = checkFiniteRange(soil, 'BSM_lon', 40, 60, 'soil_optical', problems);
    else
        problems = checkFiniteLower(soil, 'soil_spectrum_index', 1, 'soil_optical', problems);
    end
end

% ---------------------------------------------------------------------------
function problems = checkGeometryRanges(geometry, problems)
% Check viewing geometry ranges.
% Editor: Jizixing Zhu; Date: 2026-06-26.
    problems = checkFiniteRange(geometry, 'observer_zenith', 0, 90, 'geometry', problems);
    problems = checkFiniteRange(geometry, 'relative_azimuth', 0, 360, 'geometry', problems);
end

% ---------------------------------------------------------------------------
function problems = checkModelLibrary(raw, cfg, soilspectrum, problems)
% Check [model_library] paths and required NetCDF variables/shapes.
% Editor: Jizixing Zhu; Date: 2026-06-26.
    ml = cfg.model_library;
    format = 'netcdf';
    if isfield(ml, 'format')
        format = lower(strtrim(char(ml.format)));
    end
    if ~strcmp(format, 'netcdf')
        problems{end+1} = sprintf('model_library.format must be "netcdf", got "%s"', format);
        return;
    end

    rootDir = fullfile(fileparts(fileparts(raw.siteDir)), 'model_library_nc');
    [leafFile, problems] = checkLibraryFile(rootDir, ml, 'leaf_optics_file', problems);
    [atmFile,  problems] = checkLibraryFile(rootDir, ml, 'atmosphere_file', problems);
    [soilFile, problems] = checkLibraryFile(rootDir, ml, 'soil_spectrum_file', problems);

    if ~isempty(leafFile)
        requiredLeaf = {'wl','nr','Kab','Kca','Ks','Kw','Kdm','Kant','nw','phi','phiI','phiII','KcaV','KcaZ'};
        problems = checkNcVariables(leafFile, requiredLeaf, 'model_library.leaf_optics_file', problems);
    end
    if ~isempty(atmFile)
        problems = checkNcMatrixWidth(atmFile, 20, 'model_library.atmosphere_file', problems);
    end
    if soilspectrum == 0 && ~isempty(soilFile)
        [problems, ncol] = checkNcMatrixWidth(soilFile, 2, 'model_library.soil_spectrum_file', problems);
        if isfinite(ncol)
            idx = cfg.soil_optical.soil_spectrum_index + 1;
            if idx < 1 || idx > ncol
                problems{end+1} = sprintf(['soil_optical.soil_spectrum_index selects column %d, ' ...
                    'but soil_spectrum_file has %d columns'], idx, ncol);
            end
        end
    end
end

% ---------------------------------------------------------------------------
function problems = checkIntegerEnum(s, name, allowed, group, problems)
% Check an integer enum field.
% Editor: Jizixing Zhu; Date: 2026-06-26.
    if ~isfield(s, name)
        problems{end+1} = sprintf('%s.%s is required', group, name);
        return;
    end
    val = s.(name);
    if ~isnumeric(val) || ~isscalar(val) || ~ismember(val, allowed)
        problems{end+1} = sprintf('%s.%s must be one of {%s}', ...
            group, name, sprintf('%g ', allowed));
    end
end

% ---------------------------------------------------------------------------
function problems = checkFiniteField(s, name, group, problems)
% Check that a numeric field exists and is finite.
% Editor: Jizixing Zhu; Date: 2026-06-26.
    if ~isfield(s, name) || ~isnumeric(s.(name)) || ~isscalar(s.(name)) || ~isfinite(s.(name))
        problems{end+1} = sprintf('%s.%s must be a finite numeric scalar', group, name);
    end
end

% ---------------------------------------------------------------------------
function problems = checkFiniteLower(s, name, lower, group, problems)
% Check a finite scalar lower bound.
% Editor: Jizixing Zhu; Date: 2026-06-26.
    oldCount = numel(problems);
    problems = checkFiniteField(s, name, group, problems);
    if numel(problems) == oldCount && s.(name) < lower
        problems{end+1} = sprintf('%s.%s must be >= %.15g', group, name, lower);
    end
end

% ---------------------------------------------------------------------------
function problems = checkFiniteRange(s, name, lo, hi, group, problems)
% Check a finite scalar closed interval.
% Editor: Jizixing Zhu; Date: 2026-06-26.
    oldCount = numel(problems);
    problems = checkFiniteField(s, name, group, problems);
    if numel(problems) == oldCount && (s.(name) < lo || s.(name) > hi)
        problems{end+1} = sprintf('%s.%s must be in [%.15g, %.15g]', group, name, lo, hi);
    end
end

% ---------------------------------------------------------------------------
function val = getNumericOrDefault(s, name, dflt)
% Read a numeric scalar field with a default fallback.
% Editor: Jizixing Zhu; Date: 2026-06-26.
    val = dflt;
    if isfield(s, name) && isnumeric(s.(name)) && isscalar(s.(name)) && isfinite(s.(name))
        val = s.(name);
    end
end

% ---------------------------------------------------------------------------
function [filePath, problems] = checkLibraryFile(rootDir, ml, key, problems)
% Resolve and validate a relative NetCDF library file path.
% Editor: Jizixing Zhu; Date: 2026-06-26.
    filePath = '';
    if ~isfield(ml, key) || isempty(ml.(key))
        problems{end+1} = sprintf('model_library.%s is required', key);
        return;
    end
    relPath = char(ml.(key));
    if startsWith(relPath, filesep) || ~isempty(regexp(relPath, '^[A-Za-z]:', 'once'))
        problems{end+1} = sprintf('model_library.%s must be relative to data/model_library_nc', key);
        return;
    end
    parts = regexp(relPath, '[\\/]', 'split');
    filePath = fullfile(rootDir, parts{:});
    if ~endsWith(lower(filePath), '.nc')
        problems{end+1} = sprintf('model_library.%s must point to a .nc file', key);
        filePath = '';
        return;
    end
    if ~isfile(filePath)
        problems{end+1} = sprintf('model_library.%s does not exist: %s', key, filePath);
        filePath = '';
    end
end

% ---------------------------------------------------------------------------
function problems = checkNcVariables(ncPath, required, label, problems)
% Check required NetCDF variable names without loading full arrays.
% Editor: Jizixing Zhu; Date: 2026-06-26.
    try
        info = ncinfo(ncPath);
    catch ME
        problems{end+1} = sprintf('%s cannot be inspected as NetCDF: %s', label, ME.message);
        return;
    end
    names = {info.Variables.Name};
    miss = required(~ismember(required, names));
    if ~isempty(miss)
        problems{end+1} = sprintf('%s missing NetCDF variable(s): %s', label, strjoin(miss, ', '));
    end
end

% ---------------------------------------------------------------------------
function [problems, ncol] = checkNcMatrixWidth(ncPath, minWidth, label, problems)
% Check NetCDF variable "data" exists and has enough logical columns.
% Editor: Jizixing Zhu; Date: 2026-06-26.
    ncol = NaN;
    try
        info = ncinfo(ncPath);
    catch ME
        problems{end+1} = sprintf('%s cannot be inspected as NetCDF: %s', label, ME.message);
        return;
    end
    names = {info.Variables.Name};
    id = find(strcmp(names, 'data'), 1);
    if isempty(id)
        problems{end+1} = sprintf('%s missing NetCDF variable "data"', label);
        return;
    end
    sz = info.Variables(id).Size;
    if numel(sz) < 2
        problems{end+1} = sprintf('%s data variable must be a 2-D matrix', label);
        return;
    end
    ncol = min(sz(1), sz(2));  % converted libraries may be stored transposed
    if ncol < minWidth
        problems{end+1} = sprintf('%s data variable has %d columns; need at least %d', ...
            label, ncol, minWidth);
    end
end

% ---------------------------------------------------------------------------
function problems = checkEnum(s, name, allowed, problems)
% Check a string enum field.
% Editor: Jizixing Zhu; Date: 2026-06-26.
    if isfield(s, name)
        val = char(s.(name));
        if ~ismember(val, allowed)
            problems{end+1} = sprintf('switches.%s must be one of {%s}, got "%s"', ...
                name, strjoin(allowed, ', '), val);
        end
    end
end

% ---------------------------------------------------------------------------
function tf = getSwitch(raw, name, dflt)
% Read a logical switch with a default fallback.
% Editor: Jizixing Zhu; Date: 2026-06-26.
    tf = dflt;
    if isfield(raw, 'config') && isfield(raw.config, 'switches') ...
            && isfield(raw.config.switches, name)
        tf = logical(raw.config.switches.(name));
    end
end

% ---------------------------------------------------------------------------
function v = getfieldOr(s, name, dflt)
% Read a struct field with a default fallback.
% Editor: Jizixing Zhu; Date: 2026-06-26.
    if isfield(s, name); v = s.(name); else; v = dflt; end
end
