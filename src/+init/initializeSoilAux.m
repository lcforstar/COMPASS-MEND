function built = initializeSoilAux(built)
%INITIALIZESOILAUX  Pre-step soil auxiliary initialization.
%
%   built = init.initializeSoilAux(built) performs the one-time soil-water
%   auxiliary refresh that must happen after workspace allocation and before the
%   first time step. It also evaluates the initial internal-energy mirror from
%   (h,T,P_g) through the SAME unified kernel used by the solver. When heat is
%   disabled, this is a one-time finite diagnostic state; it is not advanced.
%
%   Editor: Jizixing Zhu; Date: 2026-06-28.
%   Change: Added a separate pre-step soil auxiliary stage between
%   allocateWorkspace and the main time loop.
%   Change: 2026-07-14 Jizixing Zhu -- run full B1+B2+equations refresh and
%           compute Y.U from the unified internal-energy kernel.
%   Change: 2026-07-14 Jizixing Zhu -- persist the returned U state, initialize
%           regional U column by column, and use the configured outer time step.
%   Change: 2026-07-14 Jizixing Zhu -- initialize a finite U mirror even when
%           the heat equation is disabled, without adding runtime heat solves.
%   Change: 2026-07-20 Jizixing Zhu -- seed accepted N solver history from the initialized soil endpoint.
% Change: 2026-09-14 Codex -- retain actual O2 liquid/gas capacity at t=0.

    built = soil.water.initializeFreezingHead(built);
    built = soil.water.updateWaterContentAuxiliary(built);

    sw = built.control.switches;
    nCell = size(built.Y.T, 2);
    if nCell > 1
        built = initializeRegionalInternalEnergy(built, nCell);
        return;
    end

    % Use the same switch-controlled latent heat and configured outer time step
    % that seed the runtime soil solve. Tmin=0 follows refreshSoilCoefficients'
    % runtime default; it is not the [-60,90] temperature-inversion bracket.
    L_f = resolveInitialLatentHeat(sw);
    built.workspace.step.L_f = L_f;
    opts = struct('dt_s', built.control.time.Delt_t, 'L_f', L_f, 'Tmin', 0, ...
        'DSTMAX', 0, 'KT', 1, 'NoTime_KT', 1, ...
        'Srt', built.workspace.flux.Srt, 'airHeatImpedance', 0);

    % Run the full B1+B2+equations refresh the solver uses.
    [built, equations] = soil.refreshSoilCoefficients(built, opts);

    % N uses its own effective-liquid/ice old-step snapshots.  Seed them only
    % after the authoritative initial soil endpoint has been refreshed.
    [built.Y, ~] = soil.nitrogen.materializeAcceptedState( ...
        built.Y, built.workspace.derived, built.domain);
    built.Y.soil_gas.oxygen_phase=struct( ...
        'theta_liquid',built.workspace.derived.Theta_LL, ...
        'theta_gas',built.workspace.derived.Theta_g, ...
        'temperature_C',built.Y.T(:));

    % A heat-disabled model still owns a finite initial U mirror. Build only
    % the coefficient-level heat state on a value copy; do not enable the heat
    % equation, refresh a heat boundary, or add a runtime heat solve.
    if ~isfield(equations, 'heat') || ~isfield(equations.heat, 'energyStorage') ...
            || isempty(equations.heat.energyStorage)
        energyBuilt = built;
        energyBuilt.control.switches.lThmrlefc = 1;
        energyEquations = soil.buildEquationCoefficients(energyBuilt, opts);
        equations.heat = energyEquations.heat;
    end
    assert(isfield(equations.heat, 'energyStorage') ...
        && ~isempty(equations.heat.energyStorage), ...
        'initializeSoilAux:missingInitialEnergy', ...
        'The unified internal-energy kernel did not return initial storage.');
    built = init.setPrognosticInternalEnergy( ...
        built, equations.heat.energyStorage);
end

function built = initializeRegionalInternalEnergy(built, nCell)
% Compute regional initial U by reusing the validated single-column pathway.
% The retained B2 compatibility adapter is column based, so each regional cell
% is selected without flattening, refreshed independently, and gathered only
% after its endpoint energy has been assembled to nodes.
% Editor: Jizixing Zhu; Date: 2026-07-14.
    U = zeros(built.domain.NN, nCell);
    for cellIndex = 1:nCell
        column = selectRegionalSoilColumn(built, cellIndex, nCell);
        column = init.allocateWorkspace(column);
        column = init.initializeSoilAux(column);
        U(:, cellIndex) = column.Y.U;
    end

    built.Y.U = U;
    built.workspace.step.trial.U = U;
    built.workspace.step.previous.U = U;
end

function column = selectRegionalSoilColumn(built, cellIndex, nCell)
% Select the state, soil parameters, and forcing for one regional soil column.
% Editor: Jizixing Zhu; Date: 2026-07-14.
    column = built;
    column.domain.kind = 'point';
    column.domain.nCell = 1;
    column.Y = selectRegionalLeaves(built.Y, cellIndex, nCell);
    column.parameters.soil = selectRegionalLeaves( ...
        built.parameters.soil, cellIndex, nCell);
    column.forcing = selectRegionalLeaves(built.forcing, cellIndex, nCell);

    if isfield(column, 'workspace')
        column = rmfield(column, 'workspace');
    end
    if isfield(column, 'diagnostics')
        column = rmfield(column, 'diagnostics');
    end
end

function value = selectRegionalLeaves(value, cellIndex, nCell)
% Recursively remove only the trailing regional-cell dimension from a tree.
% Element-endpoint arrays are [NL x nD x nCell], while node/layer/time fields
% are [N x nCell]; both become their normal single-column shapes here.
% Editor: Jizixing Zhu; Date: 2026-07-14.
    if isstruct(value)
        names = fieldnames(value);
        for k = 1:numel(names)
            name = names{k};
            value.(name) = selectRegionalLeaves( ...
                value.(name), cellIndex, nCell);
        end
        return;
    end

    if ~(isnumeric(value) || islogical(value)) || isempty(value)
        return;
    end

    if ~ismatrix(value) && size(value, ndims(value)) == nCell
        indices = repmat({':'}, 1, ndims(value));
        indices{end} = cellIndex;
        value = value(indices{:});
    elseif ismatrix(value) && size(value, 2) == nCell
        value = value(:, cellIndex);
    end
end

function L_f = resolveInitialLatentHeat(switches)
% Resolve initialization latent heat exactly like runOneStep's staged state.
% Editor: Jizixing Zhu; Date: 2026-07-14.
    if isfield(switches, 'lFreezing') && switches.lFreezing == 1
        L_f = equations.getPhysicalConstant('lambda_fus');
    else
        L_f = 0;
    end
    override = str2double(io.getEnv('LFUSION_OVERRIDE'));
    if isfinite(override)
        L_f = override;
    end
end
