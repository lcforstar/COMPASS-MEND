function plan = resolveSimulationPlan(config)
%RESOLVESIMULATIONPLAN  Convert site context and switches into an input plan.
%
%   plan = init.resolveSimulationPlan(config) is the single initialization
%   decision point for point-field, point-laboratory, and future regional
%   tasks. A false component switch allows that component's input file or TOML
%   section to be absent; a true component switch never follows a missing-input
%   fallback. The plan controls only loading, validation, and assembly. It does
%   not alter solver equations or boundary values.
%
%   Current W1 support:
%     point + field       full coupled Changwu-compatible initialization
%     point + laboratory scalar standalone Richards initialization
%     regional            full initialization through initializeRegional
%
%   Editor: Jizixing Zhu; Date: 2026-07-11.
%   Change: 2026-08-16 Codex -- load soil-organic initialization only when
%     enabled and require the mineral-N solver for atomic C/N closure.
%   Change: 2026-08-21 Jizixing Zhu -- resolve the DOM--MEMC--litter ownership
%     contract at initialization and expose it in the component plan.
%   Change: 2026-08-29 Codex -- use soil_carbon_model as the sole master and
%     require native POC/MOC/MB/IOM initialization priors.

    assert(isstruct(config) && isscalar(config), ...
        'resolveSimulationPlan:badConfig', 'config must be a scalar struct.');

    spatialMode = 'point';
    pointContext = 'field';
    if isfield(config, 'simulation')
        simulation = config.simulation;
        spatialMode = readEnum(simulation, 'spatial_mode', spatialMode, ...
            {'point','regional'}, 'simulation.spatial_mode');
        pointContext = readEnum(simulation, 'point_context', pointContext, ...
            {'field','laboratory'}, 'simulation.point_context');
    end
    if strcmp(spatialMode, 'regional') && ~strcmp(pointContext, 'field')
        error('resolveSimulationPlan:regionalLaboratory', ...
            'simulation.point_context must be "field" when spatial_mode is "regional".');
    end

    plan = struct();
    plan.spatial_mode = spatialMode;
    plan.point_context = pointContext;
    plan.is_laboratory = strcmp(spatialMode, 'point') && strcmp(pointContext, 'laboratory');
    plan.driver = 'runOneStep';

    if plan.is_laboratory
        validateLaboratorySwitches(config);
        plan.driver = 'solver.runWaterOnlyExperimentStep';
        % Field canopy and plant hydraulics still consume the configured root
        % profile when lRwuef=0; only the laboratory context disables that
        % whole component and may omit root.csv.
        plan.inputs = struct( ...
            'soil', true, 'initial', true, 'water_boundary', true, ...
            'forcing', false, 'plant', false, 'root', false, ...
            'vegetation', false, 'scope', false, ...
            'organic_initial', false);
        plan.components = struct( ...
            'soil_water', true, 'soil_heat', false, 'soil_air', false, ...
            'soil_nitrogen', false, 'plant_cn', false, ...
            'root_water_uptake', false, 'root_nitrogen_uptake', false, ...
            'soil_organic', false, ...
            'soil_dom', false, 'explicit_soil_microbes', false, ...
            'litter_input', false, ...
            'canopy', false, 'scope', false);
    else
        plantCnOn = readSwitch(config, 'plant_cn_cycle');
        rootUptakeOn = readSwitch(config, 'root_water_uptake');
        soilNitrogenOn = readSwitch(config, 'solve_nitrogen');
        soilOrganicOn = readSoilCarbonMaster(config);
        biogeochemistry = init.resolveBiogeochemistryMode(config);
        rootNitrogenOn = readSwitch(config, 'root_nitrogen_uptake');
        if rootNitrogenOn && ~soilNitrogenOn
            error('resolveSimulationPlan:rootNitrogenNeedsSoilNitrogen', ...
                ['switches.root_nitrogen_uptake=true requires ' ...
                 'switches.solve_nitrogen=true.']);
        end
        if soilOrganicOn && ~soilNitrogenOn
            error('resolveSimulationPlan:soilOrganicNeedsSoilNitrogen', ...
                ['switches.soil_carbon_model=true requires ' ...
                 'switches.solve_nitrogen=true for atomic C/N closure.']);
        end
        plan.inputs = struct( ...
            'soil', true, 'initial', true, 'water_boundary', true, ...
            'forcing', true, 'plant', plantCnOn, 'root', true, ...
            'vegetation', true, 'scope', true, ...
            'organic_initial', soilOrganicOn);
        plan.components = struct( ...
            'soil_water', true, ...
            'soil_heat', readSwitch(config, 'solve_soil_heat'), ...
            'soil_air', readSwitch(config, 'solve_soil_air'), ...
            'soil_nitrogen', soilNitrogenOn, ...
            'soil_organic', soilOrganicOn, ...
            'soil_dom', biogeochemistry.dom_enabled, ...
            'explicit_soil_microbes', ...
                biogeochemistry.explicit_microbes_enabled, ...
            'litter_input', biogeochemistry.litter_enabled, ...
            'plant_cn', plantCnOn, ...
            'root_water_uptake', rootUptakeOn, ...
            'root_nitrogen_uptake', rootNitrogenOn, ...
            'canopy', true, 'scope', true);
        plan.biogeochemistry = biogeochemistry;
    end
end

% -------------------------------------------------------------------------
function value = readSoilCarbonMaster(config)
% Read the sole supported soil-carbon switch.
% Editor: Codex; Date: 2026-08-29.
    value = readOptionalSwitch(config, 'soil_carbon_model', false);
end

% -------------------------------------------------------------------------
function value = readOptionalSwitch(config, name, defaultValue)
% Read one optional boolean module switch without changing legacy sites.
% Editor: Jizixing Zhu; Date: 2026-08-16.
    value = defaultValue;
    if isfield(config, 'switches') && isfield(config.switches, name)
        raw = config.switches.(name);
        if ~(islogical(raw) || (isnumeric(raw) && isscalar(raw) && isfinite(raw)))
            error('resolveSimulationPlan:badSwitch', ...
                'switches.%s must be boolean.', name);
        end
        value = logical(raw);
    end
end

% -------------------------------------------------------------------------
function validateLaboratorySwitches(config)
% Enforce the component set supported by the standalone water driver.
% Editor: Jizixing Zhu; Date: 2026-07-11.
    disabled = {'solve_soil_heat','solve_soil_air','solve_nitrogen', ...
        'plant_cn_cycle','soil_freezing','root_water_uptake', ...
        'root_nitrogen_uptake'};
    active = {};
    for i = 1:numel(disabled)
        if readSwitch(config, disabled{i})
            active{end+1} = disabled{i}; %#ok<AGROW>
        end
    end
    if ~isempty(active)
        error('resolveSimulationPlan:laboratoryComponentOn', ...
            ['The W1 laboratory context is standalone soil water only. ' ...
             'Set these switches to false or use point_context="field": %s'], ...
            strjoin(active, ', '));
    end
    if ~isfield(config, 'numerics')
        error('resolveSimulationPlan:missingNumerics', ...
            'site.toml requires [numerics] for a laboratory water calculation.');
    end
    for name = {'hydraulic_kernel','matric_kernel'}
        key = name{1};
        if ~isfield(config.numerics, key) ...
                || ~strcmpi(strtrim(char(string(config.numerics.(key)))), 'scalar')
            error('resolveSimulationPlan:laboratoryKernel', ...
                'numerics.%s must be "scalar" for the W1 laboratory context.', key);
        end
    end
end

% -------------------------------------------------------------------------
function value = readSwitch(config, name)
% Read one required module switch without inferring state from missing input.
% Editor: Jizixing Zhu; Date: 2026-07-11.
    if ~isfield(config, 'switches') || ~isfield(config.switches, name)
        error('resolveSimulationPlan:missingSwitch', ...
            'site.toml requires switches.%s; missing input never disables a module.', name);
    end
    raw = config.switches.(name);
    if ~(islogical(raw) || (isnumeric(raw) && isscalar(raw) && isfinite(raw)))
        error('resolveSimulationPlan:badSwitch', 'switches.%s must be boolean.', name);
    end
    value = logical(raw);
end

% -------------------------------------------------------------------------
function value = readEnum(s, name, fallback, allowed, label)
% Read one optional lower-case enum with a validated default.
% Editor: Jizixing Zhu; Date: 2026-07-11.
    value = fallback;
    if isfield(s, name) && ~isempty(s.(name))
        value = lower(strtrim(char(string(s.(name)))));
    end
    if ~ismember(value, allowed)
        error('resolveSimulationPlan:badEnum', ...
            '%s must be one of {%s}, got "%s".', label, strjoin(allowed, ', '), value);
    end
end
