function schedule=appendStrawReturnSchedule(schedule,config,enabled)
%APPENDSTRAWRETURNSCHEDULE Append site.toml straw-return management events.
%
%   [straw_return] follows the indexed-array convention used by
%   [fertilization]. Event carbon is public [kg C ha^-1] and converted with
%   1 kg ha^-1 = 1e-5 g cm^-2. Nitrogen is independently initialized from the
%   event C:N. Surface events require zero depths; incorporated events require
%   a positive interval. All events retain material_type='straw'.
%
% Editor: Jizixing Zhu; Date: 2026-08-28.
% Change: 2026-08-30 Codex -- append explicit zero-TAN fields for straw.
% Change: 2026-10-10 Wenjie Liu -- append unit mixing entries under the
%   uniform-depth placement contract shared with direct manure events.

    if ~enabled||~isfield(config,'straw_return'), return; end
    s=config.straw_return;
    if ~isfield(s,'events_time')||isempty(s.events_time), return; end
    required={'events_carbon_kg_ha','events_cn_ratio','placement', ...
        'depth_start_cm','depth_end_cm'};
    missing=required(~isfield(s,required));
    if ~isempty(missing)
        error('appendStrawReturnSchedule:missingField', ...
            'Missing [straw_return] field(s): %s.',strjoin(missing,', '));
    end
    times=datetime(string(s.events_time), ...
        'InputFormat','yyyy-MM-dd''T''HH:mm:ss');
    n=numel(times);
    carbonKg=eventRow(s.events_carbon_kg_ha,n,'events_carbon_kg_ha');
    cn=eventRow(s.events_cn_ratio,n,'events_cn_ratio');
    depthStart=sharedDepth(s.depth_start_cm,n,'depth_start_cm');
    depthEnd=sharedDepth(s.depth_end_cm,n,'depth_end_cm');
    placement=lower(string(s.placement));
    assert(numel(placement)==n&&all(ismember(placement, ...
        ["surface","incorporated"])), ...
        'appendStrawReturnSchedule:placement', ...
        'placement must contain one surface/incorporated value per event.');
    incorporated=placement=="incorporated";
    validDepth=(incorporated&depthStart>=0&depthEnd>depthStart) ...
        |(~incorporated&depthStart==0&depthEnd==0);
    assert(all(carbonKg>0)&&all(cn>0)&&all(validDepth), ...
        'appendStrawReturnSchedule:values', ...
        'Straw carbon, C:N, placement, or depths are invalid.');

    carbon=carbonKg.*1e-5; nitrogen=carbon./cn;
    if ~isfield(schedule,'mixing_efficiency')
        schedule.mixing_efficiency=ones(1,numel(schedule.time));
    end
    schedule.mixing_efficiency=[schedule.mixing_efficiency,ones(1,n)];
    % The staging adapter ignores legacy three-pool splits: material identity
    % selects the native straw AOM rate kernel after the event is accepted.
    schedule.time=[schedule.time,times(:).'];
    schedule.material_type=[schedule.material_type,repmat({'straw'},1,n)];
    schedule.placement=[schedule.placement,cellstr(placement(:)).'];
    schedule.depth_start_cm=[schedule.depth_start_cm,depthStart];
    schedule.depth_end_cm=[schedule.depth_end_cm,depthEnd];
    schedule.carbon_g_cm2=[schedule.carbon_g_cm2,carbon];
    schedule.nitrogen_g_cm2=[schedule.nitrogen_g_cm2,nitrogen];
    if ~isfield(schedule,'organic_nitrogen_g_cm2')
        schedule.organic_nitrogen_g_cm2=schedule.nitrogen_g_cm2(1:end-n);
    end
    if ~isfield(schedule,'tan_nitrogen_g_cm2')
        schedule.tan_nitrogen_g_cm2=zeros(1,numel(schedule.time)-n);
    end
    if ~isfield(schedule,'tan_fraction')
        schedule.tan_fraction=zeros(1,numel(schedule.time)-n);
    end
    schedule.organic_nitrogen_g_cm2=[schedule.organic_nitrogen_g_cm2,nitrogen];
    schedule.tan_nitrogen_g_cm2=[schedule.tan_nitrogen_g_cm2,zeros(1,n)];
    schedule.tan_fraction=[schedule.tan_fraction,zeros(1,n)];
    schedule.n_events=numel(schedule.time);
end

function value=sharedDepth(value,n,name)
% Expand one site-level incorporation depth to every straw event.
    value=value(:).';
    if isscalar(value), value=repmat(value,1,n); end
    assert(numel(value)==n&&all(isfinite(value)), ...
        'appendStrawReturnSchedule:depth', ...
        '[straw_return].%s must be one finite scalar (or one value per event).',name);
end

function value=eventRow(value,n,name)
% Validate one finite site event array and return row orientation.
% Editor: Jizixing Zhu; Date: 2026-08-28.
    value=value(:).';
    assert(numel(value)==n&&all(isfinite(value)), ...
        'appendStrawReturnSchedule:length', ...
        '[straw_return].%s must contain %d finite values.',name,n);
end
