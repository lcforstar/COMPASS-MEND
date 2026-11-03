function sw = buildSwitches(rawConfig)
%BUILDSWITCHES  Normalize [switches] into legacy flags + conditional cascade.
%
%   sw = init.buildSwitches(rawConfig) maps the friendly site.toml switch names
%   onto the legacy code-name flags the main loop reads (lNefc, lThmrlefc, ...),
%   converts enum strings to their legacy integers, and applies the
%   conditional-switch cascade. Pure config translation -- no physics.
%
%   Inputs:
%     rawConfig  struct  raw.config (uses .switches)
%
%   Returns sw  struct  legacy switch flags (0/1 or int) with cascade applied.
%
%   Enum integers and the cascade are verified against the legacy code
%   (setInitialThermal.m, legacy ModelSettings.setup_conditional_switches).
% Editor: Jizixing Zhu; Date: 2026-06-27.
% Change: 2026-07-13 Jizixing Zhu -- default finite-rate phase relaxation off so strict SFCC remains the sole freeze-thaw closure.
% Change: 2026-07-11 Jizixing Zhu -- restore lRsuef as the independent master switch for all root nitrogen uptake.
% Change: 2026-07-03 Jizixing Zhu — add sw.lEbalAnderson (ebal_anderson_latent, default true = legacy ebal_clm45 implicit default); makes the ebal latent-heat closure config-driven instead of relying on the omitted-field fallback in runStepEbal.
% Change: 2026-07-18 Jizixing Zhu -- remove the retired Fourier-G ebal feedback switch.
% Change: 2026-07-20 Jizixing Zhu -- document independent soil-N, PlantCN, and passive root-N switches.
% Change: 2026-08-16 Codex -- add the optional soil organic-cycle master switch.
% Change: 2026-08-21 Jizixing Zhu -- add independently gated DOM, explicit
%   microbial, and live litter-coupling modules; all default off.
% Change: 2026-08-28 Codex -- make soil_carbon_model the master switch;
%   organic N follows it and all dependent soil-C modules are forced off when
%   the master is false. Only soil_carbon_model is accepted.
% Change: 2026-09-21 Wenjie Liu -- add the one-way mineral-fertilizer N
%   source tracer without coupling it into any total-N process switch.
% Change: 2026-10-23 Wenjie Liu -- expose the soil-N kinetic pH response;
%   absent configuration preserves the reference response, false means fPH=1.

    s = rawConfig.switches;
    sw = struct();

    boolMap = { ...
        'solve_soil_heat',               'lThmrlefc'; ...
        'solve_soil_air',                'lSoilairefc'; ...
        'solve_nitrogen',                'lNefc'; ...
        'plant_cn_cycle',                'lPlantCN'; ...          % Plant C/N inventory; not the root-N master switch
        'soil_freezing',                 'lFreezing'; ...
        'root_water_uptake',             'lRwuef'; ...
        'root_nitrogen_uptake',          'lRsuef'; ...
        'soil_freezing_curve',           'lSFCC'; ...
        'temperature_affects_retention', 'lhThmrl'; ...
        'effective_thermal_properties',  'lThmrlCondCap'; ...
        'temperature_driven_liquid_flow','lKLT_Switch'; ...
        'pressure_driven_vapor_diffusion','lDVa_Switch'; ...
        'pressure_driven_liquid_flow',   'lKLa_Switch'; ...
        'dual_porosity',                 'lDualPor'; ...
        'track_total_nitrogen',          'lNtotal'; ...
        'biological_n_fixation',         'lBNF'};
    for i = 1:size(boolMap, 1)
        sw.(boolMap{i, 2}) = double(logical(s.(boolMap{i, 1})));
    end

    % enum switches: friendly string -> legacy integer (verified in code)
    sw.lSWCC      = enumToInt(s.water_retention_model, ...
        {'van_genuchten', 'clapp_hornberger'}, [1 0], 'water_retention_model');
    sw.lW_Chg     = enumToInt(s.wetting_heat_method, ...
        {'prunty', 'milly'}, [1 0], 'wetting_heat_method');
    sw.lThermCond = enumToInt(s.thermal_conductivity_model, ...
        {'de_vries', 'tian', 'peters_lidard', 'johansen', 'farouki'}, [1 3 2 2 4], 'thermal_conductivity_model');

    % ebal latent-heat closure: Anderson formulation. Legacy ebal_clm45 defaults
    % true when cfg lacks the field (no config path existed); made explicit here.
    sw.lEbalAnderson = readOptionalBool(s, 'ebal_anderson_latent', true);
    nitrogenPH=true;
    if isfield(s,'nitrogen_ph_response')
        nitrogenPH=s.nitrogen_ph_response;
        assert(isscalar(nitrogenPH)&&(islogical(nitrogenPH) ...
            ||(isnumeric(nitrogenPH)&&isreal(nitrogenPH)&&isfinite(nitrogenPH) ...
            &&ismember(nitrogenPH,[0 1]))), ...
            'buildSwitches:invalidNitrogenPHResponse', ...
            'nitrogen_ph_response must be a scalar logical or zero/one.');
    end
    sw.lNitrogenPH=double(logical(nitrogenPH));
    sw.lSoilOrganic = readSoilCarbonMaster(s);
    sw.lOrganicNitrogen = sw.lSoilOrganic;
    sw.lSoilDOM = readOptionalBool(s, 'dissolved_organic_matter', false);
    sw.lExplicitSoilMicrobes = readOptionalBool( ...
        s, 'explicit_soil_microbes', false);
    sw.lLitterInput = readOptionalBool(s, 'litter_input', false);
    sw.lMicrobialDormancy = readOptionalBool(s, 'microbial_dormancy', false);
    sw.lMicrobialTransport = readOptionalBool(s, 'microbial_transport', false);
    sw.lFertilizerSourceTracer = readOptionalBool( ...
        s, 'fertilizer_source_tracer', false);
    % Change: 2026-10-30 Wenjie Liu -- separate manure-organic and TAN labels
    % from the existing mineral-fertilizer cohort; optional and off by default.
    if isfield(s,'manure_source_tracer')
        value=s.manure_source_tracer;
        assert(isscalar(value)&&(islogical(value)||(isnumeric(value) ...
            &&isreal(value)&&isfinite(value)&&ismember(value,[0 1]))), ...
            'buildSwitches:manureSourceTracer','manure_source_tracer must be boolean.');
    end
    sw.lManureSourceTracer = readOptionalBool(s,'manure_source_tracer',false);
    if sw.lSoilOrganic == 0
        sw.lSoilDOM = 0;
        sw.lExplicitSoilMicrobes = 0;
        sw.lLitterInput = 0;
        sw.lMicrobialDormancy = 0;
        sw.lMicrobialTransport = 0;
    end

    % cascade-only flags: legacy defaults (1); no direct TOML switch
    sw.lHystrs = 1;       % HYDRUS/Lenhard wetting-history state, legacy default
    sw.lDallAmicoPhaseChange = readOptionalBool(s, 'dallamico_phase_change', false);
    sw.lDVT_Switch = 1;   % temperature-driven vapor diffusion
    sw.lKaT_Switch = 1;   % temperature-driven dry-air diffusion
    sw.lKaa_Switch = 1;   % pressure-driven dry-air diffusion

    % conditional cascade (legacy ModelSettings.setup_conditional_switches)
    if sw.lSoilairefc == 0
        sw.lKaT_Switch = 0; sw.lKaa_Switch = 0;
        sw.lDVa_Switch = 0; sw.lKLa_Switch = 0;
    end
    if sw.lThmrlefc == 0
        sw.lKLT_Switch = 0; sw.lDVT_Switch = 0;
    end

    % forcing diurnal-temperature disaggregation method (new path; no legacy code).
    % 'legacy_sine' = fixed cos peak 15:00 (MeteoData); 'sunrise_anchored' = sunrise/
    % sunset-anchored sine daytime + linear night cooling (init.buildForcing).
    if isfield(s, 'diurnal_temp_method')
        sw.diurnal_temp_method = lower(strtrim(char(s.diurnal_temp_method)));
    else
        sw.diurnal_temp_method = 'legacy_sine';
    end
end

% ---------------------------------------------------------------------------
function value = readSoilCarbonMaster(s)
% Read the sole supported soil-carbon switch.
% Editor: Codex; Date: 2026-08-29.
    value = readOptionalBool(s, 'soil_carbon_model', false);
end

% ---------------------------------------------------------------------------
function v = enumToInt(s, keys, vals, name)
% Map a friendly enum string to its legacy integer; error on an unknown value.
% Editor: Jizixing Zhu; Date: 2026-06-27.
    s = lower(strtrim(char(s)));
    idx = find(strcmp(s, keys), 1);
    if isempty(idx)
        error('buildSwitches:badEnum', '%s = "%s" not in {%s}', name, s, strjoin(keys, ', '));
    end
    v = vals(idx);
end

% ---------------------------------------------------------------------------
function v = readOptionalBool(s, name, defaultValue)
% Read an optional friendly boolean switch.
% Editor: Jizixing Zhu; Date: 2026-06-28.
    if isfield(s, name)
        v = double(logical(s.(name)));
    else
        v = double(logical(defaultValue));
    end
end
