function step_workspace = warmStartTrial(step_workspace, accepted_state)
%WARMSTARTTRIAL  Copy accepted prognostic fields into the soil trial state.
%
%   Inputs:
%     step_workspace [inout] Soil step workspace containing trial/previous.
%     accepted_state [in] Accepted model prognostic state Y.
%
%   Return:
%     step_workspace [out] Workspace with trial initialized from Y.
%
%   The previous layer is intentionally untouched because it remains the
%   t_(n-1) extrapolation and rollback anchor until acceptance.
%   Editor: Jizixing Zhu; Date: 2026-07-16.
%   Change: 2026-07-20 Jizixing Zhu -- warm-start soil-mode candidates and clear the
%     trial-only dual-porosity exchange rate.
    names = { ...
        'h','T','p_gas','U','nh4','no3','urea','C_NH4_granule', ...
        'canopy_air','scope_soil','soil_gas','soil_modes','surface_water'};
    for i = 1:numel(names)
        name = names{i};
        if isfield(accepted_state, name)
            step_workspace.trial.(name) = accepted_state.(name);
        end
    end
    if isfield(step_workspace, 'trial_aux') ...
            && isfield(step_workspace.trial_aux, 'dual_porosity')
        step_workspace.trial_aux.dual_porosity.SinkIm_current(:) = 0;
    end
    if isfield(step_workspace, 'trial_aux') ...
            && isfield(step_workspace.trial_aux, 'surface_water')
        step_workspace.trial_aux.surface_water.active = false;
    end
end
