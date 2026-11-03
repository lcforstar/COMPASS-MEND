function phen = stageCandidate(model, kt, dt_s)
%STAGECANDIDATE Freeze current-step crop events before canopy exchange.
% Inputs are an accepted model shell, outer forcing index KT [-], and dt_s
% [s]. Return is an uncommitted phenology candidate; no accepted Y is written.
% Thermal phenology advances once per forcing-clock day. Event pulses are
% cleared on every outer step, including an hourly retry from accepted state.
% Observed_lai_zero uses only the current precompiled forcing edge: a strictly
% positive observed green LAI followed by exact zero is the user-selected
% harvest proxy. It is not a physiological-maturity or field-date observation.
% Sources: unchanged CoLM daily kernel (winterWheatPhenology); COMPASS observed
% management adaptation authorized 2026-09-30. No future runtime lookup.
% Editor: Wenjie Liu; Date: 2026-09-30.
    built = model.built;
    phen = built.Y.vegetation.phenology;
    pulseNames = {'crop_planting_event','crop_emergence_event', ...
        'crop_harvest_event','crop_grainfill_event','crop_maturity_event'};
    for k = 1:numel(pulseNames), phen.state.(pulseNames{k}) = false; end
    time = built.control.time;
    elapsed = step.readAcceptedTime(model.runtime);
    timestamp = time.start_time + seconds(elapsed);
    pd = built.vegetation_state.phen_driver;
    forcingStart = time.start_time;
    if isfield(pd,'forcing_start_time'), forcingStart = pd.forcing_start_time; end
    forcingDay = floor(days(dateshift(timestamp,'start','day') - ...
        dateshift(forcingStart,'start','day'))) + 1;
    lastDay = 0;
    if isfield(phen,'last_advance_day_index'), lastDay = phen.last_advance_day_index; end
    if forcingDay ~= lastDay
        phen = vegetation.integrated.phenology.advanceDay( ...
            phen,pd,built.parameters.veg.structure.PFT_type,timestamp, ...
            time.start_time,built.domain.site.latitude,built.parameters.veg.phenology);
        phen.last_advance_day_index = forcingDay;
    end
    if isfield(built.vegetation_state,'crop_lifecycle')
        lifecycle = built.vegetation_state.crop_lifecycle;
        if strcmp(lifecycle.harvest_rule,'observed_lai_zero')
            assert(kt >= 1 && kt <= numel(lifecycle.harvest_event), ...
                'phenology:stageCandidate:forcingIndex', ...
                'Observed crop lifecycle forcing does not cover this KT.');
            if lifecycle.harvest_event(kt) && phen.state.crop_live
                phen.state.crop_live = false;
                phen.state.crop_phase = 4;
                phen.state.crop_harvest_event = true;
                phen.state.offset_flag = true;
                phen.state.harvested_flag = true;
                phen.state.harvest_doy = day(timestamp,'dayofyear');
                phen.state.crop_a_leaf = 0; phen.state.crop_a_root = 0;
                phen.state.crop_a_stem = 0; phen.state.crop_a_repr = 0;
                phen.state.dormant_flag = true;
                phen.state.phen = 0;
                phen.now = 0;
            end
        end
    end
    assert(isscalar(dt_s) && isfinite(dt_s) && dt_s > 0, ...
        'phenology:stageCandidate:duration','Outer duration must be positive seconds.');
end
