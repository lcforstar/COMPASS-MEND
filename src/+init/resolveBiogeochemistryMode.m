function mode = resolveBiogeochemistryMode(config)
%RESOLVEBIOGEOCHEMISTRYMODE Validate revised soil C/N module ownership.
%
%   MODE = init.resolveBiogeochemistryMode(CONFIG) is the single fail-fast
%   decision point for the optional soil-biogeochemistry extensions.  With
%   When the master is false, all dependent modules are inactive and the
%   existing COMPASS empirical mineralization is the sole fallback.
%
%   MEMC requires DOM because its DOC state is owned by Y.soil_dom, not copied
%   into a second accepted store. Live litter input requires PlantCN.
%
%   Input/return units: switches are dimensionless booleans; MODE contains
%   only ownership flags and the configured litter source string.
%
%   Editor: Jizixing Zhu; Date: 2026-08-21.
%   Change: 2026-08-29 Codex -- cascade all organic submodules from the sole
%   soil_carbon_model master.
    assert(isstruct(config) && isscalar(config), ...
        'resolveBiogeochemistryMode:badConfig', ...
        'config must be a scalar struct.');
    s = config.switches;
    organicOn = readSoilCarbonMaster(s);
    nitrogenOn = readBool(s, 'solve_nitrogen', false);
    plantCnOn = readBool(s, 'plant_cn_cycle', false);
    domOn = organicOn && readBool(s, 'dissolved_organic_matter', false);
    microbeOn = organicOn && readBool(s, 'explicit_soil_microbes', false);
    litterOn = organicOn && readBool(s, 'litter_input', false);
    dormancyOn = organicOn && readBool(s, 'microbial_dormancy', false);
    transportOn = organicOn && readBool(s, 'microbial_transport', false);

    if domOn && ~nitrogenOn
        error('resolveBiogeochemistryMode:domDependencies', ...
            'dissolved_organic_matter=true requires solve_nitrogen=true.');
    end
    if microbeOn && ~domOn
        error('resolveBiogeochemistryMode:microbeNeedsDom', ...
            'explicit_soil_microbes=true requires dissolved_organic_matter=true.');
    end
    if organicOn && ~microbeOn
        error('resolveBiogeochemistryMode:soilCarbonNeedsMemc', ...
            ['soil_carbon_model=true requires explicit_soil_microbes=true ', ...
             'in the native implementation.']);
    end
    if dormancyOn && ~microbeOn
        error('resolveBiogeochemistryMode:dormancyNeedsMicrobes', ...
            'microbial_dormancy=true requires explicit_soil_microbes=true.');
    end
    if transportOn && ~microbeOn
        error('resolveBiogeochemistryMode:transportNeedsMicrobes', ...
            'microbial_transport=true requires explicit_soil_microbes=true.');
    end
    litterSource = 'plant_cn';
    if isfield(config, 'litter') && isfield(config.litter, 'source')
        litterSource = lower(strtrim(char(string(config.litter.source))));
    end
    if litterOn
        if ~organicOn
            error('resolveBiogeochemistryMode:litterNeedsOrganic', ...
                'litter_input=true requires soil_carbon_model=true.');
        end
        if ~strcmp(litterSource, 'plant_cn')
            error('resolveBiogeochemistryMode:badLitterSource', ...
                'litter.source must be "plant_cn".');
        end
        if ~plantCnOn
            error('resolveBiogeochemistryMode:litterNeedsPlantCN', ...
                'litter_input=true with source=plant_cn requires plant_cn_cycle=true.');
        end
    end
    strawOwnership='none';
    if isfield(config,'straw_return')&&hasStrawEvents(config.straw_return)
        if ~isfield(config.straw_return,'ownership')
            error('resolveBiogeochemistryMode:missingStrawOwnership', ...
                ['straw_return.ownership must explicitly be ', ...
                 '"external_addition" or "crop_residue".']);
        end
        strawOwnership=lower(strtrim(char(string( ...
            config.straw_return.ownership))));
        if ~ismember(strawOwnership,{'external_addition','crop_residue'})
            error('resolveBiogeochemistryMode:badStrawOwnership', ...
                ['straw_return.ownership must be "external_addition" ', ...
                 'or "crop_residue".']);
        end
        if litterOn&&strcmp(strawOwnership,'crop_residue')
            error('resolveBiogeochemistryMode:duplicateStrawResidue', ...
                ['PlantCN litter_input already owns crop-residue litter; ', ...
                 'a management straw event marked crop_residue would ', ...
                 'double count the same biomass.']);
        end
    end

    mode = struct( ...
        'organic_enabled', organicOn, ...
        'dom_enabled', domOn, ...
        'explicit_microbes_enabled', microbeOn, ...
        'litter_enabled', litterOn, ...
        'organic_nitrogen_enabled', organicOn, ...
        'microbial_dormancy_enabled', dormancyOn, ...
        'microbial_transport_enabled', transportOn, ...
        'litter_source', litterSource, ...
        'straw_return_ownership',strawOwnership, ...
        'memc_reaction_owner', microbeOn, ...
        'legacy_empirical_mineralization', ~organicOn, ...
        'strict_fallback', ~organicOn);
end

function value=hasStrawEvents(section)
% True only when the management section contains at least one event.
% Editor: Codex; Date: 2026-08-28.
    value=isstruct(section)&&isfield(section,'events_time') ...
        &&~isempty(section.events_time);
end

function value = readSoilCarbonMaster(s)
% Read the sole supported soil-carbon switch.
% Editor: Codex; Date: 2026-08-29.
    value = readBool(s, 'soil_carbon_model', false);
end

function value = readBool(container, name, fallback)
% Read one optional scalar boolean without inferring an active process.
    value = fallback;
    if ~isfield(container, name), return; end
    raw = container.(name);
    if ~(islogical(raw) || (isnumeric(raw) && isscalar(raw) && isfinite(raw)))
        error('resolveBiogeochemistryMode:badSwitch', ...
            'switches.%s must be boolean.', name);
    end
    value = logical(raw);
end
