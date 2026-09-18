function [Y, wasMaterialized] = materializeAcceptedState(Y, derived, domain)
%MATERIALIZEACCEPTEDSTATE  Complete accepted N memory without changing pools.
%
%   [Y, wasMaterialized] = soil.nitrogen.materializeAcceptedState(Y, derived,
%   domain) fills only missing nitrogen solver memory from an already accepted
%   soil endpoint.  It never changes dissolved/ice N pools or O2.  Fresh runs
%   obtain their exact t=0 effective-water history here; old checkpoints that
%   predate Y.soil_nitrogen receive a clearly marked compatibility seed.
%
%   This is a state-schema operation, not an N reaction or transport step.
%
%   Inputs:
%     Y [inout] struct Accepted prognostic state; dissolved/ice N pools are unchanged.
%     derived [in] struct Accepted soil liquid/ice water layers [m3 m^-3].
%     domain [in] struct Soil grid dimensions and layer thickness [cm].
%
%   Returns:
%     Y [out] struct Same physical pools plus complete Y.soil_nitrogen memory.
%     wasMaterialized [out] logical True only when compatibility/history fields were seeded.
%
%   Editor: Jizixing Zhu; Date: 2026-07-20.
%   Change: 2026-07-20 Jizixing Zhu — retain the accepted root-N sink field
%           required by diagnostics and restart-equivalent solver memory.
    wasMaterialized = false;
    NL = domain.NL;
    NN = domain.NN;
    nD = domain.grid.nD;
    nS = domain.grid.nS;
    if ~isfield(Y, 'soil_nitrogen') || ~isstruct(Y.soil_nitrogen)
        Y.soil_nitrogen = soil.nitrogen.emptyAcceptedState();
        Y.soil_nitrogen.initialized_from_legacy_checkpoint = true;
        wasMaterialized = true;
    end
    state = Y.soil_nitrogen;
    if ~isfield(state, 'schema_version')
        state.schema_version = 1;
        state.initialized_from_legacy_checkpoint = true;
        wasMaterialized = true;
    end
    if ~isfield(state, 'theta_liquid_prev') || isempty(state.theta_liquid_prev)
        state.theta_liquid_prev = derived.Theta_LL;
        wasMaterialized = true;
    end
    if ~isfield(state, 'theta_ice_prev') || isempty(state.theta_ice_prev)
        state.theta_ice_prev = derived.Theta_II;
        wasMaterialized = true;
    end
    if ~isfield(state, 'variables') || isempty(fieldnames(state.variables))
        state.variables = defaultVariables(NL, NN, nD, nS);
        wasMaterialized = true;
    end
    if ~isfield(state, 'flux') || isempty(fieldnames(state.flux))
        state.flux = defaultFlux(NL);
        wasMaterialized = true;
    end
    if ~isfield(state, 'initialized_from_legacy_checkpoint')
        state.initialized_from_legacy_checkpoint = false;
    end
    Y.soil_nitrogen = state;
end

function variables = defaultVariables(NL, NN, nD, nS)
% Dynamic, accepted N-solver memory only; static parameters are not copied.
    z2 = zeros(NL, nD);
    z3 = zeros(NL, nD, nS);
    variables = struct( ...
        'Retard', ones(NL, nD), 'Dsh', z3, 'Dif', z3, ...
        'Dsh_h', z3, 'Dsh_T', z3, 'Dsh_A', z3, ...
        'Smin', z2, 'Sden', z2, 'SVden', z2, 'Snit', z2, ...
        'SVnit', z2, 'Svol', z2, 'SVvol', z2, 'Shys', z2, ...
        'Vpore', z2, 'WFPS', z2, 'eT', zeros(NL, 1), ...
        'eTheta', zeros(NL, 1), 'ePH', zeros(NL, 1), ...
        'Sdis', z2, 'f_aerobic', ones(NL, nD), ...
        'f_anaerobic', zeros(NL, nD), 'N_SAVE', zeros(3, 3, 5), ...
        'SinkS', z3, ...
        'NH4_total', zeros(NN, 2), 'NO3_total', zeros(NN, 2), ...
        'CHK2', inf(NN, 1), 'CHK3', inf(NN, 1), ...
        'R_O2_sink', zeros(NL, 1), 'nit_mm_factor', zeros(NL, nD), ...
        'nit_M_total', zeros(NL, nD), 'eT_nit', zeros(NL, nD), ...
        'Snit_N2o', zeros(NL, nD));
end

function flux = defaultFlux(NL)
% Two slots: slot 1 is accepted n, slot 2 is the candidate n+1.
    flux = struct( ...
        'RetardBAR', zeros(NL, 1), ...
        'cTop_NH4', zeros(2, 1), 'cBot_NH4', zeros(2, 1), ...
        'cTop_NO3', zeros(2, 1), 'cBot_NO3', zeros(2, 1), ...
        'cvTop', zeros(2, 1), 'cvBot', zeros(2, 1), ...
        'cvTop1', zeros(2, 1), 'cvBot1', zeros(2, 1), ...
        'cvSnit', zeros(2, 1), 'cvSmin', zeros(2, 1), ...
        'cvSden', zeros(2, 1), 'cvSvol', zeros(2, 1), ...
        'cvSinkSnh4', zeros(2, 1), 'cvSinkSno3', zeros(2, 1), ...
        'Cvnit_N2o', zeros(2, 1), 'Cvdenit_N2o', zeros(2, 1), ...
        'Cvdenit_N2', zeros(2, 1), 'cvnh4OLD', zeros(2, 1), ...
        'cvno3OLD', zeros(2, 1), 'cvnh4', zeros(2, 1), ...
        'cvno3', zeros(2, 1));
end
