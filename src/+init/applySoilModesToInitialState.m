function [Y, state] = applySoilModesToInitialState(Y, state, built, rawConfig)
%APPLYSOILMODESTOINITIALSTATE  Apply switch-dependent soil-mode IC transforms.
%
%   [Y, state] = init.applySoilModesToInitialState(Y, state, built, rawConfig)
%   applies initialization transforms that need both the accepted initial soil
%   state and the soil-mode state container. At present this covers HYDRUS-style
%   dual porosity. The wetting-history arrays are initialized in
%   init.buildSoilModes because they only need soil parameters.
%
%   Inputs:
%     Y          struct  accepted initial soil state from init.buildInitialState
%     state      struct  Y.soil_modes candidate from init.buildSoilModes
%     built      struct  buildState container with domain, switches and soil params
%     rawConfig  struct  raw.config; uses [dualporosity] theta_s_immobile and
%                        theta_r_immobile when switches.dual_porosity is true.
%
%   Returns:
%     Y      struct  with accepted Theta_L, Theta_U, and Theta_I phase state;
%                    when lDualPor=1, Theta_L and theta are mobile-domain water.
%     state  struct  dual_porosity fields filled with immobile-domain IC values.
%
%   Closure:
%     HYDRUS InitW split is conservative at t=0:
%       theta_mobile + theta_immobile = theta_observed.
%     Van Genuchten theta_s is not changed; the user-provided theta_s remains
%     the mobile-domain saturation, matching the legacy dualPorIni contract.
%
%   Editor: Jizixing Zhu; Date: 2026-06-28.
%   Change: 2026-07-12 Jizixing Zhu -- initialize accepted total-water and ice phase
%           state after all soil-mode transforms.
%   Change: 2026-07-20 Jizixing Zhu -- read immobile-domain material constants from
%     immutable parameters; accepted state now contains prognostic memory only.

    sw = built.control.switches;
    if ~getFlag(sw, 'lDualPor', 0)
        Y = initializeAcceptedPhaseState(Y);
        return;
    end

    NL = built.domain.NL;
    nD = getGridDim(built.domain.grid, 'nD', 2);
    theta_s_mobile = built.parameters.soil.theta_s(:);
    if any(theta_s_mobile <= 0)
        error('applySoilModesToInitialState:badThetaS', ...
            'soil.theta_s must be positive for dual-porosity initialization.');
    end

    dpParameters = built.parameters.soil.modes.dual_porosity;
    theta_s_im = dpParameters.theta_s_im(:);
    theta_r_im = dpParameters.theta_r_im(:);

    theta_obs = Y.Theta_L;
    theta_m  = zeros(NL, nD);
    theta_im = zeros(NL, nD);
    for i = 1:NL
        denom = max(theta_s_mobile(i) + theta_s_im(i), 1e-12);
        for j = 1:nD
            tobs = max(0, min(theta_s_mobile(i) + theta_s_im(i), theta_obs(i, j)));
            theta_m(i, j)  = tobs * theta_s_mobile(i) / denom;
            theta_im(i, j) = tobs - theta_m(i, j);
        end
    end

    Y.Theta_L = theta_m;
    Y.theta = [theta_m(:, 1); theta_m(end, nD)];

    state.dual_porosity.Theta_im       = theta_im;
    state.dual_porosity.Theta_im_old   = theta_im;
    state.dual_porosity.Theta_m_old    = theta_m;
    Y = initializeAcceptedPhaseState(Y);
end

% ---------------------------------------------------------------------------
function Y = initializeAcceptedPhaseState(Y)
% Initialize phase state after the final initial liquid-water transform.
    Y.Theta_U = Y.Theta_L;
    Y.Theta_I = zeros(size(Y.Theta_L));
end

% ---------------------------------------------------------------------------
function value = getFlag(sw, name, fallback)
% Read one switch flag with a scalar fallback.
% Editor: Jizixing Zhu; Date: 2026-06-28.
    value = fallback;
    if isfield(sw, name) && ~isempty(sw.(name))
        value = sw.(name);
    end
end

% ---------------------------------------------------------------------------
function value = getGridDim(gridCfg, name, fallback)
% Read one grid dimension with a scalar fallback.
% Editor: Jizixing Zhu; Date: 2026-06-28.
    value = fallback;
    if isfield(gridCfg, name) && ~isempty(gridCfg.(name))
        value = gridCfg.(name);
    end
end
