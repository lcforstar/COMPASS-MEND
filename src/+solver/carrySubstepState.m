function built = carrySubstepState(built)
%CARRYSUBSTEPSTATE  Commit only the accepted soil/solve state layers.
%
%   This is the shared pending-free carry ritual used by the natural outer
%   step and the standalone water experiment driver:
%     previous <- accepted Y; Y solve fields <- trial; water-content time
%     layers <- refreshed accepted endpoint.
%   It does not advance the runtime clock or commit canopy, plant, management,
%   checkpoint, or diagnostic side effects.
%
%   Editor: Jizixing Zhu; Date: 2026-07-11.
%   Change: 2026-07-14 Jizixing Zhu -- carry prognostic U through previous/Y/trial.
%   Change: 2026-07-20 Jizixing Zhu -- carry soil-mode candidates only when a
%     substep/outer transaction accepts; trial exchange rates remain auxiliary.

    built.workspace.step = commitPreviousFromY(built.workspace.step, built.Y);
    built.Y = acceptTrialState(built.Y, built.workspace.step.trial);
    built.Y = carryThetaLayer(built.Y, built.workspace.derived);
end

% -------------------------------------------------------------------------
function step = commitPreviousFromY(step, Y)
% Preserve the accepted pre-trial state as the extrapolation/rollback anchor.
    names = {'h','T','p_gas','U','nh4','no3','urea', ...
        'C_NH4_granule','soil_gas','soil_modes'};
    for i = 1:numel(names)
        name = names{i};
        if isfield(Y, name)
            step.previous.(name) = Y.(name);
        end
    end
end

% -------------------------------------------------------------------------
function Y = acceptTrialState(Y, trial)
% Copy only prognostic fields owned by the coupled solve trial.
    names = {'h','T','p_gas','U','nh4','no3','urea', ...
        'C_NH4_granule','soil_gas','soil_modes'};
    for i = 1:numel(names)
        name = names{i};
        if isfield(trial, name)
            Y.(name) = trial.(name);
        end
    end
end

% -------------------------------------------------------------------------
function Y = carryThetaLayer(Y, derived)
% Advance the accepted water-content layer to the refreshed trial endpoint.
    if isfield(Y, 'Theta_L') && isfield(derived, 'Theta_LL')
        Y.Theta_L = derived.Theta_LL;
    end
    if isfield(Y, 'Theta_U') && isfield(derived, 'Theta_UU')
        Y.Theta_U = derived.Theta_UU;
    end
    if isfield(Y, 'Theta_I') && isfield(derived, 'Theta_II')
        Y.Theta_I = derived.Theta_II;
    end
end
