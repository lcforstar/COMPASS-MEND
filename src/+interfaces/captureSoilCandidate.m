function candidate = captureSoilCandidate(built)
%CAPTURESOILCANDIDATE Capture the uncommitted soil endpoint of one outer step.
%
%   candidate = interfaces.captureSoilCandidate(built) copies the accepted
%   trial endpoint required by the coupled h/U/P_g discretization.  The soil
%   candidate is not yet an accepted state: step.commitTransaction is the
%   only public installer after all cross-domain candidates validate.
%
%   Inputs:
%     built     [struct] Legacy model shell containing solver trial/derived data.
%
%   Returns:
%     candidate [struct] Soil-owned trial state and derived endpoint needed to
%               carry h, U, p_gas, solutes, and phase water layers.
%
%   Physical contract:
%     The coupled soil solve has already formed a discrete candidate at t^(n+1).
%     This function neither solves nor changes it; it makes that candidate
%     explicit so the transaction can install it atomically with surface and
%     vegetation candidates.  `derived` is included because theta layers are
%     diagnosed from the accepted nonlinear endpoint during the carry step.
%
%   Editor: Jizixing Zhu; Date: 2026-07-20.
% Change: 2026-07-20 Jizixing Zhu -- reserve the soil-owned accepted N candidate slot.
% Change: 2026-07-20 Jizixing Zhu -- project an explicit soil-owned field set instead
%   of copying surface and vegetation trial payloads into the soil candidate.
    assert(isstruct(built) && isfield(built, 'workspace') ...
        && isfield(built.workspace, 'step') ...
        && isfield(built.workspace.step, 'trial'), ...
        'interfaces:captureSoilCandidate:missingTrial', ...
        'An accepted outer soil candidate requires built.workspace.step.trial.');
    assert(isfield(built.workspace, 'derived'), ...
        'interfaces:captureSoilCandidate:missingDerived', ...
        'An accepted outer soil candidate requires built.workspace.derived.');

    trial = built.workspace.step.trial;
    soilNames = {'h','T','p_gas','U','nh4','no3','urea', ...
        'C_NH4_granule','soil_gas','soil_modes'};
    soilTrial = struct();
    for i = 1:numel(soilNames)
        name = soilNames{i};
        if isfield(trial, name)
            soilTrial.(name) = trial.(name);
        end
    end
    candidate = struct( ...
        'trial', soilTrial, ...
        'derived', built.workspace.derived, ...
        'nitrogen', []);
end
