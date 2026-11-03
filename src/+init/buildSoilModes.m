function [soil, state] = buildSoilModes(built, rawConfig)
%BUILDSOILMODES  Split legacy SoilVariables support fields by ownership.
%
%   [soil, state] = init.buildSoilModes(built) recreates the non-primary
%   fields that init.defineSoilVariables used to add to the old monolithic
%   SoilVariables struct. It does not duplicate the accepted soil state
%   (h/T/theta/N/gas), which is already owned by built.Y, and it does not
%   allocate per-step derived scratch, which is owned by init.allocateWorkspace.
%
%   Inputs:
%     built  struct  buildState container with built.domain,
%                    built.control.switches and built.parameters.soil
%
%   Returns:
%     soil   struct  built.parameters.soil with .modes added:
%                    .modes.freeze.imped [NN x 1],
%                    .modes.retention.{Ss,h_trans}. Material values stay in
%                    their canonical soil fields (xk, wilting_point,
%                    field_capacity, porosity, sand/silt/clay, ks) and are not
%                    duplicated here.
%     state   struct  cross-step soil-mode memory:
%                    .hysteresis.* and .dual_porosity.* initialised like legacy
%                    defineSoilVariables.
%
%   Legacy provenance:
%     init.defineSoilVariables.m after the state-field copy block. The legacy
%     IH/XWRE wetting-history branch is currently disabled in COMPASS.m, but the
%     state placeholders are still required for mode-compatible initialization
%     and future runOneStep adapters.
%
%   Editor: Jizixing Zhu; Date: 2026-06-28.
%   Change: 2026-07-20 Jizixing Zhu -- move dual-porosity material constants out of
%     prognostic state and into immutable soil-mode parameters.

    if nargin < 2 || isempty(rawConfig)
        rawConfig = struct();
    end

    soil = built.parameters.soil;
    NL = built.domain.NL;
    NN = built.domain.NN;
    nD = getGridDim(built.domain.grid, 'nD', 2);
    sw = built.control.switches;

    soil.modes = struct();
    soil.modes.freeze = struct();
    soil.modes.freeze.imped = 7 * ones(NN, 1);       % [-] legacy ice impedance factor
    soil.modes.retention = struct();
    soil.modes.retention.Ss = 1.0e-5;                % [1/cm] specific storage
    soil.modes.retention.h_trans = 0.1;              % [cm] saturated-transition half width
    soil.modes.dual_porosity = struct( ...
        'theta_s_im', readDualPorosityScalar( ...
            rawConfig, 'theta_s_immobile', 0.05) * ones(NL, 1), ...
        'theta_r_im', readDualPorosityScalar( ...
            rawConfig, 'theta_r_immobile', 0.0) * ones(NL, 1), ...
        'OmegaExch', readDualPorosityScalar( ...
            rawConfig, 'omega_exch', 0.01));

    state = struct();
    state.freeze = struct();
    state.freeze.ISFT = 0;                       % [-] legacy initial soil-freezing type flag
    state.hysteresis = buildHysteresisState(NL, nD, soil, getFlag(sw, 'lHystrs', 1));
    state.dual_porosity = buildDualPorosityState(NL, nD);
end

% ---------------------------------------------------------------------------
function h = buildHysteresisState(NL, nD, soil, enabled)
% Build wetting-history state used by the HYDRUS/Lenhard branch.
% Editor: Jizixing Zhu; Date: 2026-06-28.
    IPATH = 7;
    h = struct();
    h.IH     = ones(NL, 1);
    h.XWRE   = zeros(NL, nD);
    h.RHSW   = zeros(NL, IPATH);
    h.RASW   = zeros(NL, IPATH);
    h.IPSW   = ones(NL, 1);
    h.JJH    = zeros(NL, 1);
    h.MPSW   = 3 * ones(NL, 1);
    h.PHSW   = zeros(NL, 1);
    h.SARW   = zeros(NL, 1);
    h.RSWAW  = ones(NL, 1);
    h.AlphaW = zeros(NL, 1);
    h.SARWI  = zeros(NL, 1);
    if ~logical(enabled)
        return;
    end

    h.RASW(:, 1) = 1.0;
    h.RHSW(:, 1) = 0.0;
    for ij = 4:2:IPATH
        h.RASW(:, ij) = 1.0;
    end
    for ij = 3:2:(IPATH - 1)
        h.RASW(:, ij) = 0.0;
    end
    if mod(IPATH, 2) == 1
        h.RASW(:, IPATH) = 0.0;
    end
    h.RHSW(:, 2) = 0.0;
    h.RASW(:, 2) = 1.0;
    h.AlphaW = 2.0 * soil.alpha(:);
    h.SARWI  = 0.05 * ones(NL, 1);
    h.SARW   = h.SARWI;
end

% ---------------------------------------------------------------------------
function d = buildDualPorosityState(NL, nD)
% Build immobile-domain water-content memory for the dual-porosity branch.
% Editor: Jizixing Zhu; Date: 2026-06-28.
    d = struct();
    d.Theta_im     = zeros(NL, nD);
    d.Theta_im_old = zeros(NL, nD);
    d.Theta_m_old  = zeros(NL, nD);
end

% ---------------------------------------------------------------------------
function value = readDualPorosityScalar(rawConfig, name, fallback)
% Read one immutable dual-porosity material parameter.
% Editor: Jizixing Zhu; Date: 2026-07-20.
    value = fallback;
    if isfield(rawConfig, 'dualporosity') ...
            && isfield(rawConfig.dualporosity, name) ...
            && ~isempty(rawConfig.dualporosity.(name))
        value = rawConfig.dualporosity.(name);
    end
    assert(isnumeric(value) && isscalar(value) && isfinite(value), ...
        'buildSoilModes:badDualPorosityParameter', ...
        'dualporosity.%s must be a finite scalar.', name);
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

% ---------------------------------------------------------------------------
function value = getFlag(sw, name, fallback)
% Read one switch flag with a scalar fallback.
% Editor: Jizixing Zhu; Date: 2026-06-28.
    value = fallback;
    if isfield(sw, name) && ~isempty(sw.(name))
        value = sw.(name);
    end
end
