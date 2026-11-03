function pending = buildPending(input_candidates, post_ebal_candidates, daily_candidates, fertilizer_applied)
%BUILDPENDING  Assemble producer candidates into the fixed step transaction.
%
%   Inputs:
%     input_candidates      [in] LAI-input and leaf-temperature candidates.
%     post_ebal_candidates  [in] Hydraulic, surface, water, snow, and LUNA
%                           candidates produced after the finalized ordinary
%                           or single-correction SCOPE/energy pass.
%     daily_candidates      [in] Phenology, root-N, and plant C/N candidates.
%     fertilizer_applied    [in] Accept-only fertilizer event mask/state.
%
%   Return:
%     pending [struct] Destination-keyed transaction from step.emptyPending.
%
%   This is the single producer-to-destination mapping. It prevents the
%   runOneStep orchestrator from growing a scattered sequence of field writes.
%   Required-field completeness is checked by step.validatePending at commit.
%   Editor: Jizixing Zhu; Date: 2026-07-16.
%   Change: 2026-07-17 Jizixing Zhu -- map staged ClimaLand snow state and
%     diagnostics into the accepted-step transaction.
%   Change: 2026-07-17 Jizixing Zhu -- document that post-ebal candidates
%     belong to the finalized ordinary or bounded-correction surface pass.
%   Change: 2026-07-18 Jizixing Zhu -- stage finalized force-restore history as an
%     accept-only prognostic candidate.
% Change: 2026-07-20 Jizixing Zhu — materialize root nitrogen as a typed soil request.
% Change: 2026-07-20 Jizixing Zhu -- eliminate duplicated workspace flux fields;
%   legacy output values are projected from typed requests after acceptance.
    pending = step.emptyPending();
    pending.candidate.vegetation = vegetation.materializeCandidate( ...
        input_candidates, post_ebal_candidates, daily_candidates);
    pending.candidate.surface = ...
        surface.materializeCandidate(post_ebal_candidates);

    pending.request.vegetation_to_soil.water = ...
        post_ebal_candidates.root_soil_exchange;
    pending.request.vegetation_to_soil.nitrogen = ...
        coupling.solveRootSoilNitrogenExchange( ...
        daily_candidates.root_nitrogen, size(daily_candidates.SinkS, 1));
    pending.request.surface_to_soil.water = ...
        post_ebal_candidates.surface_soil_water;

    pending.diagnostics.canopy_water = ...
        post_ebal_candidates.canopy_water_diagnostics;
    pending.diagnostics.snow = post_ebal_candidates.snow_diagnostics;
    pending.diagnostics.root_demand = daily_candidates.root_demand;
    pending.diagnostics.root_uptake_M1_diag = ...
        daily_candidates.root_uptake_diagnostics;
    if ~isempty(post_ebal_candidates.plant_hydraulic_W)
        pending.diagnostics.v2_hydraulic = ...
            post_ebal_candidates.hydraulic_diagnostics;
    end
    if ~isempty(daily_candidates.plant_nee)
        if isfield(daily_candidates.plant_nee, 'NEE')
            pending.ledger.vegetation_to_atmosphere. ...
                carbon_result.NEE = daily_candidates.plant_nee.NEE;
        end
        if isfield(daily_candidates.plant_nee, 'Resp')
            pending.ledger.vegetation_to_atmosphere. ...
                carbon_result.Resp = daily_candidates.plant_nee.Resp;
        end
    end
    pending.management.fertilizer_applied = fertilizer_applied;
end
