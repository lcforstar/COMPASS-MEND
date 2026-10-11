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
% Change: 2026-09-30 Wenjie Liu -- consume the pre-energy phenology candidate;
% daily thermal progression and one-shot management events are staged once.
% Change: 2026-09-30 Wenjie Liu -- carry the ephemeral accepted-input
% record needed to reconcile annual crop growth with actual soil N uptake.
    if isfield(model.built.workspace.scope_step,'phenology_candidate')
        phenology_next = model.built.workspace.scope_step.phenology_candidate;
    else
        phenology_next = vegetation.integrated.phenology.stageCandidate(model,kt,dt_outer_s);
    end
    nitrogen_demand = vegetation.integrated.preparePlantNitrogenDemand( ...
        model, kt, dt_outer_s);
    root_nitrogen = vegetation.root.stepNitrogenUptake(model, nitrogen_demand);

    plant_cn_next = [];
    plant_nee = [];
    plant_cn_replay = [];
    if readField(model.built.control.switches, 'lPlantCN', 0) == 1
        [plant_cn_next, plant_nee, plant_cn_debug] = vegetation.integrated.stepPlantCN( ...
            model, kt, dt_outer_s, phenology_next, ...
            nitrogen_demand, root_nitrogen);
        if isfield(plant_cn_debug,'plant_cn_replay')
            plant_cn_replay = plant_cn_debug.plant_cn_replay;
        end
    end
    candidates = struct( ...
        'phenology', phenology_next, ...
        'plant_cn', plant_cn_next, ...
        'plant_cn_replay', plant_cn_replay, ...
        'plant_nee', plant_nee, ...
        'SinkS', root_nitrogen.SinkS, ...
        'root_nitrogen', root_nitrogen, ...
        'root_demand', struct( ...
            'nh4', root_nitrogen.demand_nh4, ...
            'no3', root_nitrogen.demand_no3), ...
        'root_uptake_diagnostics', root_nitrogen.uptake_diag);
end

function value = readField(container, name, fallback)
% Read a nonempty struct field or return a local fallback.
    value = fallback;
    if isfield(container, name) && ~isempty(container.(name))
        value = container.(name);
    end
end
