function candidates = stageDailyVegetation(model, kt, dt_outer_s)
%STAGEDAILYVEGETATION  Stage phenology, root-N uptake, and plant C/N state.
%
%   Inputs:
%     model      [in] Model shell after post-ebal water and hydraulic staging.
%     kt         [-] One-based outer-step index.
%     dt_outer_s [s] Actual outer-step duration.
%
%   Return:
%     candidates [struct] Pending phenology, plant C/N, root-N sink,
%                diagnostics, and optional NEE/Resp values.
%
%   Phenology advances only on the first outer step of a day. Root nitrogen
%   uptake remains independent from PlantCN and is controlled by its own
%   switch. No candidate is committed by this function.
%   Editor: Jizixing Zhu; Date: 2026-07-16.
% Change: 2026-07-20 Jizixing Zhu — expose the typed vegetation-to-soil nitrogen request source.
    phenology_next = advanceDailyPhenology(model, kt, dt_outer_s);
    nitrogen_demand = vegetation.integrated.preparePlantNitrogenDemand( ...
        model, kt, dt_outer_s);
    root_nitrogen = vegetation.root.stepNitrogenUptake(model, nitrogen_demand);

    plant_cn_next = [];
    plant_nee = [];
    if readField(model.built.control.switches, 'lPlantCN', 0) == 1
        [plant_cn_next, plant_nee] = vegetation.integrated.stepPlantCN( ...
            model, kt, dt_outer_s, phenology_next, ...
            nitrogen_demand, root_nitrogen);
    end
    candidates = struct( ...
        'phenology', phenology_next, ...
        'plant_cn', plant_cn_next, ...
        'plant_nee', plant_nee, ...
        'SinkS', root_nitrogen.SinkS, ...
        'root_nitrogen', root_nitrogen, ...
        'root_demand', struct( ...
            'nh4', root_nitrogen.demand_nh4, ...
            'no3', root_nitrogen.demand_no3), ...
        'root_uptake_diagnostics', root_nitrogen.uptake_diag);
end

function phenology_next = advanceDailyPhenology(model, kt, dt_outer_s)
% Advance the pending phenology state at the first step of each model day.
    phenology_next = model.built.Y.vegetation.phenology;
    time_control = model.built.control.time;
    steps_per_day = max(1, round(86400 / time_control.Delt_t));
    if mod(kt - 1, steps_per_day) ~= 0
        return;
    end
    timestamp = time_control.start_time + seconds( ...
        step.readAcceptedTime(model.runtime) + dt_outer_s);
    phenology_next = vegetation.integrated.phenology.advanceDay( ...
        phenology_next, model.built.vegetation_state.phen_driver, ...
        model.built.parameters.veg.structure.PFT_type, timestamp, ...
        time_control.start_time, model.built.domain.site.latitude, ...
        model.built.parameters.veg.phenology);
end

function value = readField(container, name, fallback)
% Read a nonempty struct field or return a local fallback.
    value = fallback;
    if isfield(container, name) && ~isempty(container.(name))
        value = container.(name);
    end
end
