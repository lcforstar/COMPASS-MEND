function schedule=appendRootResidueSchedule(schedule,config,enabled)
%APPENDROOTRESIDUESCHEDULE Append prescribed paired root-residue C/N events.
% The optional [root_residue] has aligned events_time, events_carbon_kg_ha
% and events_nitrogen_kg_ha arrays. C>0 and N>=0 are independently prescribed;
% 1 kg/ha = 1e-5 g/cm2. Empty arrays leave the original schedule unchanged.
% Events use material_type=root_residue and placement=root_profile. Their
% spatial allocation is supplied separately from the unchanged Jackson RAI;
% depth_start/end are inactive placeholders, never a uniform mixing band.
% This external crop-residue route requires soil organic C/N on, PlantCN and
% live litter off. No plant stock, TAN, or prognostic root pool is created.
% Invalid array lengths, dates, masses or ownership fail before appending.
% Editor: Wenjie Liu; Date: 2026-10-30.
    if ~isfield(config,'root_residue'), return; end
    s=config.root_residue;
    required={'events_time','events_carbon_kg_ha','events_nitrogen_kg_ha'};
    assert(isstruct(s)&&isscalar(s)&&all(isfield(s,required)), ...
        'appendRootResidueSchedule:missingField', ...
        '[root_residue] requires events_time, events_carbon_kg_ha and events_nitrogen_kg_ha.');
    n=numel(s.events_time);
    carbon=eventMass(s.events_carbon_kg_ha,n,'events_carbon_kg_ha');
    nitrogen=eventMass(s.events_nitrogen_kg_ha,n,'events_nitrogen_kg_ha');
    if n==0, return; end
    assert(enabled&&isfield(config,'switches') ...
        &&readFlag(config.switches,'soil_carbon_model') ...
        &&~readFlag(config.switches,'plant_cn_cycle') ...
        &&~readFlag(config.switches,'litter_input'), ...
        'appendRootResidueSchedule:ownership', ...
        'Root residue requires soil_carbon_model=true, plant_cn_cycle=false and litter_input=false.');
    assert(isfield(config,'vegetation') ...
        &&isfield(config.vegetation,'root_dist_mode') ...
        &&config.vegetation.root_dist_mode==0, ...
        'appendRootResidueSchedule:profile', ...
        'Root residue currently requires the unchanged configured Jackson root profile.');
    assert(all(carbon>0)&&all(nitrogen>=0), ...
        'appendRootResidueSchedule:mass', ...
        'Root carbon must be positive and root nitrogen nonnegative.');
    times=datetime(string(s.events_time),'InputFormat','yyyy-MM-dd''T''HH:mm:ss');
    assert(all(~isnat(times(:))),'appendRootResidueSchedule:time', ...
        'Root residue event times must be finite datetimes.');
    before=numel(schedule.time);
    if ~isfield(schedule,'mixing_efficiency'),schedule.mixing_efficiency=ones(1,before);end
    if ~isfield(schedule,'organic_nitrogen_g_cm2')
        schedule.organic_nitrogen_g_cm2=schedule.nitrogen_g_cm2;
    end
    if ~isfield(schedule,'tan_nitrogen_g_cm2'),schedule.tan_nitrogen_g_cm2=zeros(1,before);end
    if ~isfield(schedule,'tan_fraction'),schedule.tan_fraction=zeros(1,before);end
    schedule.time=[schedule.time,times(:).'];
    schedule.material_type=[schedule.material_type,repmat({'root_residue'},1,n)];
    schedule.placement=[schedule.placement,repmat({'root_profile'},1,n)];
    schedule.depth_start_cm=[schedule.depth_start_cm,zeros(1,n)];
    schedule.depth_end_cm=[schedule.depth_end_cm,zeros(1,n)];
    schedule.carbon_g_cm2=[schedule.carbon_g_cm2,carbon*1e-5];
    schedule.nitrogen_g_cm2=[schedule.nitrogen_g_cm2,nitrogen*1e-5];
    schedule.organic_nitrogen_g_cm2=[schedule.organic_nitrogen_g_cm2,nitrogen*1e-5];
    schedule.tan_nitrogen_g_cm2=[schedule.tan_nitrogen_g_cm2,zeros(1,n)];
    schedule.tan_fraction=[schedule.tan_fraction,zeros(1,n)];
    schedule.mixing_efficiency=[schedule.mixing_efficiency,ones(1,n)];
    schedule.n_events=numel(schedule.time);
end

function value=eventMass(value,n,name)
% Validate public per-event masses [kg element/ha], preserving zero N.
% Editor: Wenjie Liu; Date: 2026-10-30.
    assert(isnumeric(value)&&isreal(value)&&numel(value)==n ...
        &&all(isfinite(value(:))),'appendRootResidueSchedule:length', ...
        '[root_residue].%s must contain one finite value per event.',name);
    value=double(value(:).');
end

function on=readFlag(s,key)
% Absent optional PlantCN/live-litter flags are inactive.
% Editor: Wenjie Liu; Date: 2026-10-30.
    on=false;
    if isfield(s,key)
        assert(isscalar(s.(key))&&ismember(s.(key),[0 1]), ...
            'appendRootResidueSchedule:switch','%s must be boolean.',key);
        on=logical(s.(key));
    end
end
