function [model, fertilizer, applied_next] = stageStep(model, ~, dt_outer_s)
%STAGESTEP  Stage one scheduled fertilizer event without consuming it.
%
%   Inputs:
%     model      [inout] Model shell with warm-started nitrogen trial pools.
%     kt         [-] One-based outer-step index.
%     dt_outer_s [s] Actual outer-step duration.
%
%   Returns:
%     model        [out] Model with matched fertilizer injected into trial pools.
%     fertilizer   [struct] Per-pool applied-mass diagnostics.
%     applied_next [logical array] Candidate once-only event-consumption flags.
%
%   The event flags remain pending until the outer step is accepted. A failed
%   step therefore neither consumes the event nor advances nitrogen pools.
%   Editor: Jizixing Zhu; Date: 2026-07-16.
% Change: 2026-09-16 Jizixing Zhu -- match [step start, step end) like organic inputs.
    fertilizer = struct( ...
        'applied', false, 'nh4_added', 0, ...
        'no3_added', 0, 'urea_added', 0, ...
        'tracer_increment',[]);
    built = model.built;
    applied_next = model.runtime.management.fertilizer.applied;
    if ~isfield(built.control.switches, 'lNefc') ...
            || built.control.switches.lNefc ~= 1
        return;
    end
    schedule = built.sources.fertilization_schedule;
    if ~isfield(schedule, 'n_events') || schedule.n_events == 0
        return;
    end
    timestamp = built.control.time.start_time + seconds( ...
        step.readAcceptedTime(model.runtime));
    [should_apply, ~, ~, applied_next, pools] = ...
        soil.nitrogen.fertilizer.matchFertilizerEvent( ...
        schedule, timestamp, applied_next, dt_outer_s);
    if ~should_apply
        return;
    end

    trial = built.workspace.step.trial;
    soil_variables = struct( ...
        'C_NH4_new', trial.nh4, ...
        'C_NO3_new', trial.no3, ...
        'C_Urea_new', trial.urea, ...
        'C_NH4_granule', trial.C_NH4_granule, ...
        'Theta_LL', built.workspace.derived.Theta_LL);
    nitrogen_variables = struct('RHOKG', built.parameters.soil_n.RHOKG);
    grid = built.domain;
    [soil_variables, fertilizer] = ...
        soil.nitrogen.fertilizer.injectFertilizerPools( ...
        soil_variables, pools, schedule.DEPTHSTART, schedule.DEPTHEND, ...
        grid.sumdz, grid.NL, grid.NN, grid.dz, nitrogen_variables);
    model.built.workspace.step.trial.nh4 = soil_variables.C_NH4_new;
    model.built.workspace.step.trial.no3 = soil_variables.C_NO3_new;
    model.built.workspace.step.trial.urea = soil_variables.C_Urea_new;
    model.built.workspace.step.trial.C_NH4_granule = ...
        soil_variables.C_NH4_granule;
end
