function model = stageWaterUptake(model)
%STAGEWATERUPTAKE  Stage the current root-water sink under lRwuef.
%
%   Input:
%     model [inout] Model shell after the canopy energy balance.
%
%   Return:
%     model [out] Model with current nodal RWU and layer-endpoint Srt staged.
%
%   Srt [1/s] may be negative under hydraulic redistribution and is never
%   clipped. With lRwuef disabled, both current RWU and Srt are exact zero.
%   The staged `root_soil_exchange` is the typed vegetation-to-soil request;
%   it is not an independently owned soil-water state.
%   Editor: Jizixing Zhu; Date: 2026-07-16.
    scope_step = model.built.workspace.scope_step;
    n_layers = model.built.domain.NL;
    n_nodes = model.built.domain.NN;
    switches = model.built.control.switches;
    if readField(switches, 'lRwuef', 0) ~= 1
        scope_step.RWU = zeros(n_nodes, 1);
        exchange = coupling.solveRootSoilExchange( ...
            scope_step.RWU, n_layers);
        scope_step.Srt = exchange.element_endpoint_rate_1_s;
        scope_step.root_soil_exchange = exchange;
        model.built.workspace.scope_step = scope_step;
        return;
    end
    if ~isfield(scope_step, 'RWU') || isempty(scope_step.RWU)
        scope_step.RWU = zeros(n_nodes, 1);
        exchange = coupling.solveRootSoilExchange( ...
            scope_step.RWU, n_layers);
        scope_step.Srt = exchange.element_endpoint_rate_1_s;
        scope_step.root_soil_exchange = exchange;
        model.built.workspace.scope_step = scope_step;
        return;
    end
    root_uptake = scope_step.RWU(:);
    assert(numel(root_uptake) >= n_nodes, ...
        'runOneStep:badRwuShape', ...
        'Current RWU must contain at least NN=%d nodal values.', n_nodes);
    exchange = coupling.solveRootSoilExchange(root_uptake, n_layers);
    scope_step.Srt = exchange.element_endpoint_rate_1_s;
    scope_step.root_soil_exchange = exchange;
    model.built.workspace.scope_step = scope_step;
end

function value = readField(container, name, fallback)
% Read a nonempty struct field or return a local fallback.
    value = fallback;
    if isfield(container, name) && ~isempty(container.(name))
        value = container.(name);
    end
end
