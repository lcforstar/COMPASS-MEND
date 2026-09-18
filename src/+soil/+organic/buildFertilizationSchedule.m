function schedule = buildFertilizationSchedule(config, enabled)
%BUILDFERTILIZATIONSCHEDULE Normalize native organic-fertilizer AOM events.
%
%   Public carbon [kg C ha^-1] is converted to [g C cm^-2] with
%   1 kg ha^-1 = 1e-5 g cm^-2. Event C:N initializes an independent nitrogen
%   mass. Manure total N is separated into organic N and TAN; only organic N
%   enters AOM, while TAN is staged to the existing COMPASS NH4 pool at the
%   same placement depth. Material identity selects the native AOM store; no
%   retired slow/fast/direct pool partition is accepted.
%   Editor: Codex; Date: 2026-08-29.
%   Change: 2026-08-30 Codex -- add conservative manure organic-N/TAN split.
    schedule = emptySchedule();
    if ~enabled || ~isfield(config, 'fertilization'), return; end
    f = config.fertilization;
    if ~isfield(f, 'organic_events_time') || isempty(f.organic_events_time), return; end
    fields = {'organic_carbon_kg_ha','organic_cn_ratio','organic_placement', ...
        'organic_depth_start_cm','organic_depth_end_cm'};
    missing = fields(~isfield(f, fields));
    if ~isempty(missing)
        error('buildFertilizationSchedule:missingField', ...
            'Missing organic fertilizer field(s): %s.', strjoin(missing, ', '));
    end
    times = datetime(string(f.organic_events_time), ...
        'InputFormat','yyyy-MM-dd''T''HH:mm:ss');
    n = numel(times);
    carbon = eventRow(f.organic_carbon_kg_ha,n,'organic_carbon_kg_ha');
    totalCn = eventRow(f.organic_cn_ratio,n,'organic_cn_ratio');
    depthStart = sharedDepth(f.organic_depth_start_cm,n,'organic_depth_start_cm');
    depthEnd = sharedDepth(f.organic_depth_end_cm,n,'organic_depth_end_cm');
    placement = lower(string(f.organic_placement));
    if numel(placement) ~= n || any(~ismember(placement, ...
            ["surface","incorporated"]))
        error('buildFertilizationSchedule:placement', ...
            'organic_placement must contain one surface/incorporated value per event.');
    end
    incorporated=placement=="incorporated";
    validDepth=(incorporated&depthStart>=0&depthEnd>depthStart) ...
        |(~incorporated&depthStart==0&depthEnd==0);
    if any(carbon<=0)||any(totalCn<=0)||any(~validDepth)
        error('buildFertilizationSchedule:values', ...
            'Organic fertilizer masses, C:N, placement, or depths are invalid.');
    end
    materialType=repmat("manure",1,n);
    if isfield(f,'organic_material_type')
        materialType=lower(string(f.organic_material_type));
        if numel(materialType)~=n
            error('buildFertilizationSchedule:materialLength', ...
                'organic_material_type must contain one value per event.');
        end
    end
    if any(~ismember(materialType,["straw","manure"]))
        error('buildFertilizationSchedule:materialValue', ...
            'organic_material_type must contain only straw or manure.');
    end
    carbonGcm2=carbon.*1e-5;
    tanFraction=zeros(1,n);
    if isfield(f,'organic_tan_fraction')
        tanFraction=eventRow(f.organic_tan_fraction,n, ...
            'organic_tan_fraction');
    end
    if any(tanFraction<0|tanFraction>1)
        error('buildFertilizationSchedule:tanFraction', ...
            'organic_tan_fraction must be within [0,1] for every event.');
    end
    totalN=carbonGcm2./totalCn;
    tanN=totalN.*tanFraction;
    mixingEfficiency=ones(1,n);
    if isfield(f,'organic_mixing_efficiency')
        mixingEfficiency=eventRow(f.organic_mixing_efficiency,n, ...
            'organic_mixing_efficiency');
    end
    if any(mixingEfficiency<0|mixingEfficiency>1)
        error('buildFertilizationSchedule:mixingEfficiency', ...
            'organic_mixing_efficiency must be within [0,1].');
    end
    schedule=struct('n_events',n,'time',times(:).', ...
        'material_type',{cellstr(materialType(:)).'}, ...
        'placement',{cellstr(placement(:)).'}, ...
        'depth_start_cm',depthStart,'depth_end_cm',depthEnd, ...
        'carbon_g_cm2',carbonGcm2,'nitrogen_g_cm2',totalN, ...
        'organic_nitrogen_g_cm2',totalN-tanN, ...
        'tan_nitrogen_g_cm2',tanN,'tan_fraction',tanFraction, ...
        'mixing_efficiency',mixingEfficiency, ...
        'mass_units','g_per_cm2');
end

function value=sharedDepth(value,n,name)
% Expand one management depth scalar to all events; an event array is retained
% only for backward-compatible reads during this configuration migration.
    value=value(:).';
    if isscalar(value), value=repmat(value,1,n); end
    if numel(value)~=n||any(~isfinite(value))
        error('buildFertilizationSchedule:depth', ...
            '%s must be one finite scalar (or one value per event).',name);
    end
end

function value = eventRow(value,n,name)
% Validate one numeric event array and return a row.
    value = value(:).';
    if numel(value) ~= n || any(~isfinite(value))
        error('buildFertilizationSchedule:length', ...
            '%s must contain %d finite values.',name,n);
    end
end

function schedule = emptySchedule()
% Stable disabled/no-event schedule.
    schedule = struct('n_events',0,'time',datetime.empty, ...
        'material_type',{{}},'placement',{{}}, ...
        'depth_start_cm',[],'depth_end_cm',[], ...
        'carbon_g_cm2',[],'nitrogen_g_cm2',[], ...
        'organic_nitrogen_g_cm2',[],'tan_nitrogen_g_cm2',[], ...
        'tan_fraction',[],'mixing_efficiency',[], ...
        'mass_units','g_per_cm2');
end
