function built = initializeFreezingHead(built)
%INITIALIZEFREEZINGHEAD  Apply one-time initial freezing head correction.
%
%   built = soil.water.initializeFreezingHead(built) applies the initial
%   freezing-pressure head correction after workspace allocation and before the
%   first time step. It updates accepted state built.Y.h once and synchronizes
%   workspace.step.{trial,previous}.h and .h_frez.
%
%   Inputs:
%     built  [inout] struct canonical model container with built.Y,
%                     built.control.switches, and built.workspace allocated
%
%   Output:
%     built  [out] struct with initialized h/h_frez and Gamma_hh cache
%
%   Role:
%     This is the canonical replacement for the initialization part of
%     init.applySoilHeteroWithInitialFreezing. It is not a file reader and it
%     does not build legacy SoilVariables/cfg containers.
%
%   Editor: Jizixing Zhu; Date: 2026-06-28.
%   Change: Added canonical soil-water initialization entry for initial
%   freezing-pressure correction.

    requireWorkspace(built);
    phys = equations.physicalConstants();
    L_f = fusionLatentHeat(built.control.switches, phys);
    nCell = readDomainField(built.domain, 'nCell', 1);

    for c = 1:nCell
        h = soil.water.getNodeColumn(built.Y.h, c);
        T = soil.water.getNodeColumn(built.Y.T, c);
        h_frez = zeros(size(h));

        freezingNode = T <= 0;
        h_frez(freezingNode) = L_f * 1e4 .* T(freezingNode) ./ (phys.g * phys.T_freeze);
        tooLow = h_frez <= h + 1e-6;
        h_frez(tooLow) = h(tooLow) + 1e-6;
        h = h - h_frez;

        built.Y.h = soil.water.setNodeColumn(built.Y.h, c, h);
        built.workspace.step.trial.h = soil.water.setNodeColumn(built.workspace.step.trial.h, c, h);
        built.workspace.step.previous.h = soil.water.setNodeColumn(built.workspace.step.previous.h, c, h);
        built.workspace.step.trial.h_frez = soil.water.setNodeColumn( ...
            built.workspace.step.trial.h_frez, c, h_frez);
        built.workspace.step.previous.h_frez = soil.water.setNodeColumn( ...
            built.workspace.step.previous.h_frez, c, h_frez);
        built.workspace.derived.Gamma_hh = soil.water.setNodeColumn( ...
            built.workspace.derived.Gamma_hh, c, ...
            soil.water.computeGamma(h, built.parameters.soil.constants));
    end
end

function requireWorkspace(built)
% Guard that the caller allocated workspace before this state synchronization.
    if ~isfield(built, 'workspace') || ~isfield(built.workspace, 'derived') ...
            || ~isfield(built.workspace, 'step')
        error('initializeFreezingHead:missingWorkspace', ...
            'Call init.allocateWorkspace before soil.water.initializeFreezingHead.');
    end
end

function L_f = fusionLatentHeat(switches, phys)
% Return fusion latent heat only when freezing physics is enabled.
    if isfield(switches, 'lFreezing') && switches.lFreezing == 1
        L_f = phys.lambda_fus;
    else
        L_f = 0;
    end
end

function value = readDomainField(domain, name, defaultValue)
% Read a domain field with fallback.
    if isfield(domain, name) && ~isempty(domain.(name))
        value = domain.(name);
    else
        value = defaultValue;
    end
end
