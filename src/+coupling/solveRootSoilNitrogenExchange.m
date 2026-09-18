function exchange = solveRootSoilNitrogenExchange(root_nitrogen, n_layers)
%SOLVEROOTSOILNITROGENEXCHANGE Materialize the vegetation-to-soil N-sink request.
%
%   exchange = coupling.solveRootSoilNitrogenExchange(root_nitrogen,n_layers)
%   preserves the independently staged NH4 and NO3 root sinks at the two
%   endpoints of each soil element.  It does not transport N or modify soil
%   concentrations; soil nitrogen owns the realized removal when its accepted
%   transaction bridge is enabled.
%
%   Inputs:
%     root_nitrogen [struct] Output of vegetation.root.stepNitrogenUptake.
%     n_layers      [-] Number of finite elements, NL.
%
%   Returns:
%     exchange [struct] Immutable vegetation-to-soil request. Positive values
%              remove mineral nitrogen from soil [ug cm^-3 s^-1].
%
%   Physical contract:
%     The plant module owns uptake demand and the root-scale kinetic estimate.
%     The soil N transport equation owns the actual concentration change and
%     conservation ledger. Keeping this request typed prevents the current
%     `SinkS` diagnostic field from becoming a hidden cross-domain boundary.
%
%   Editor: Jizixing Zhu; Date: 2026-07-20.
% Change: 2026-07-20 Jizixing Zhu -- retain PlantCN-inventory versus passive-sink accounting explicitly.
    assert(isscalar(n_layers) && isfinite(n_layers) ...
            && n_layers >= 1 && n_layers == floor(n_layers), ...
        'coupling:solveRootSoilNitrogenExchange:badLayerCount', ...
        'n_layers must be a positive integer.');
    assert(isstruct(root_nitrogen) && isfield(root_nitrogen, 'SinkS') ...
            && isnumeric(root_nitrogen.SinkS), ...
        'coupling:solveRootSoilNitrogenExchange:missingSink', ...
        'root_nitrogen.SinkS is required.');
    sink = root_nitrogen.SinkS;
    assert(isequal(size(sink), [n_layers, 2, 2]) ...
            && all(isfinite(sink(:))), ...
        'coupling:solveRootSoilNitrogenExchange:badSinkShape', ...
        'SinkS must be finite [NL x 2 x 2] for NH4 and NO3.');

    enabled = readLogical(root_nitrogen, 'enabled', false);
    if ~enabled
        % A disabled request is a strict zero boundary.  This prevents stale
        % workspace SinkS values from removing soil N after the switch closes.
        sink(:) = 0;
        depthIntegrated = zeros(1, 2);
        totalUptake = 0;
        scheme = 'off';
        accountingMode = 'disabled';
    else
        depthIntegrated = readNumeric(root_nitrogen, 'SPUptake', zeros(1, 2));
        totalUptake = readNumeric(root_nitrogen, 'U_N_total', 0);
        scheme = readText(root_nitrogen, 'scheme', 'passive_mass_flow');
        accountingMode = readText(root_nitrogen, 'accounting_mode', ...
            'passive_sink_untracked');
    end
    exchange = struct( ...
        'schema_version', 1, ...
        'producer', 'vegetation', ...
        'receiver', 'soil', ...
        'kind', 'root_mineral_nitrogen_sink_request', ...
        'unit', 'ug cm^-3 s^-1', ...
        'positive_direction', 'soil_to_vegetation', ...
        'enabled', enabled, ...
        'sink_nh4', sink(:, :, 1), ...
        'sink_no3', sink(:, :, 2), ...
        'depth_integrated_uptake_ug_cm2_s', depthIntegrated, ...
        'uptake_total_g_m2_s', totalUptake, ...
        'scheme', scheme, ...
        'accounting_mode', accountingMode);
end

function value = readLogical(container, name, fallback)
%READLOGICAL Read one optional scalar switch without accepting non-scalars.
    value = fallback;
    if isfield(container, name) && isscalar(container.(name))
        value = logical(container.(name));
    end
end

function value = readNumeric(container, name, fallback)
%READNUMERIC Read one optional finite numeric payload.
    value = fallback;
    if isfield(container, name) && isnumeric(container.(name)) ...
            && all(isfinite(container.(name)), 'all')
        value = container.(name);
    end
end

function value = readText(container, name, fallback)
%READTEXT Read one optional request metadata label.
    value = fallback;
    if isfield(container, name) && (ischar(container.(name)) ...
            || isstring(container.(name)))
        value = char(string(container.(name)));
    end
end
