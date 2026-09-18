function built = alignInitialThermalBoundary(built)
%ALIGNINITIALTHERMALBOUNDARY  Make the initial bottom node BC-compatible.
%
%   A Dirichlet temperature acts at the bottom computational node.  Carrying a
%   different value at that same node creates an artificial first-step storage
%   impulse.  After the initial thermal auxiliaries have been built, this
%   routine aligns only that boundary node, then rebuilds h/T/P_g-dependent
%   phase, transport, and internal-energy auxiliaries.  Interior temperatures
%   remain the user-provided profile.  Flux and free-drainage energy boundaries
%   are not modified because they prescribe no bottom temperature.
%
%   The seasonal boundary depends weakly on the initialized thermal properties.
%   A bounded fixed-point loop therefore refreshes the auxiliary state until the
%   node and the recomputed boundary agree.  The audit record is stored under
%   built.p.initialization.thermal_bottom_alignment.
%
%   Editor: Jizixing Zhu; Date: 2026-07-18.
%   Change: 2026-07-18 Jizixing Zhu -- update the accepted, trial, and previous
%           temperature mirrors atomically before rebuilding internal energy.

    report = struct( ...
        'applicable', false, ...
        'applied', false, ...
        'boundary_kind', '', ...
        'initial_temperature_C', NaN, ...
        'initial_boundary_temperature_C', NaN, ...
        'initial_mismatch_K', NaN, ...
        'final_temperature_C', NaN, ...
        'final_boundary_temperature_C', NaN, ...
        'final_mismatch_K', NaN, ...
        'refresh_count', 0);

    if ~heatEnabled(built) || ~isfield(built.boundary, 'heat')
        built = storeReport(built, report);
        return;
    end
    kind = built.boundary.heat.bottom.kind;
    report.boundary_kind = kind;
    if ~any(strcmp(kind, {'seasonal_bottom_temperature','fixed_temperature'}))
        built = storeReport(built, report);
        return;
    end

    report.applicable = true;
    report.initial_temperature_C = built.Y.T(1);
    target = readBoundaryTemperature(built, kind);
    report.initial_boundary_temperature_C = target;
    report.initial_mismatch_K = built.Y.T(1) - target;

    toleranceK = 1e-8;
    maxRefresh = 3;
    for iteration = 1:maxRefresh
        target = readBoundaryTemperature(built, kind);
        mismatch = built.Y.T(1) - target;
        if abs(mismatch) <= toleranceK
            break;
        end
        built = setBottomTemperatureMirrors(built, target);
        built = init.initializeSoilAux(built);
        report.applied = true;
        report.refresh_count = iteration;
    end

    target = readBoundaryTemperature(built, kind);
    finalMismatch = built.Y.T(1) - target;
    assert(abs(finalMismatch) <= toleranceK, ...
        'alignInitialThermalBoundary:notConverged', ...
        ['Initial bottom temperature did not converge to the %s boundary: ', ...
         'mismatch %.6g K after %d refreshes.'], ...
        kind, finalMismatch, maxRefresh);
    report.final_temperature_C = built.Y.T(1);
    report.final_boundary_temperature_C = target;
    report.final_mismatch_K = finalMismatch;
    built = storeReport(built, report);
end

function built = setBottomTemperatureMirrors(built, target)
% Keep the t0 temperature time layers consistent with the rebuilt U mirrors.
% Editor: Jizixing Zhu; Date: 2026-07-18.
    built.Y.T(1) = target;
    if isfield(built, 'workspace') && isfield(built.workspace, 'step')
        if isfield(built.workspace.step, 'trial') ...
                && isfield(built.workspace.step.trial, 'T')
            built.workspace.step.trial.T(1) = target;
        end
        if isfield(built.workspace.step, 'previous') ...
                && isfield(built.workspace.step.previous, 'T')
            built.workspace.step.previous.T(1) = target;
        end
    end
end

function target = readBoundaryTemperature(built, kind)
% Read the already materializable t0 bottom temperature.
% Editor: Jizixing Zhu; Date: 2026-07-18.
    switch kind
        case 'seasonal_bottom_temperature'
            thermal = built.workspace.boundary.thermal;
            assert(isfield(thermal, 'Tbh') && ~isempty(thermal.Tbh) ...
                    && isfinite(thermal.Tbh(1)), ...
                'alignInitialThermalBoundary:missingSeasonalBoundary', ...
                'The t0 seasonal bottom temperature must be initialized first.');
            target = thermal.Tbh(1);
        case 'fixed_temperature'
            target = built.boundary.heat.bottom.configured_value;
        otherwise
            error('alignInitialThermalBoundary:unsupportedKind', ...
                'Boundary kind %s is not a temperature Dirichlet condition.', kind);
    end
end

function tf = heatEnabled(built)
% Read the normalized soil-heat switch.
% Editor: Jizixing Zhu; Date: 2026-07-18.
    tf = isfield(built.control, 'switches') ...
        && isfield(built.control.switches, 'lThmrlefc') ...
        && built.control.switches.lThmrlefc == 1;
end

function built = storeReport(built, report)
% Preserve initialization provenance outside the prognostic state.
% Editor: Jizixing Zhu; Date: 2026-07-18.
    if ~isfield(built, 'p') || ~isstruct(built.p)
        built.p = struct();
    end
    if ~isfield(built.p, 'initialization') || ~isstruct(built.p.initialization)
        built.p.initialization = struct();
    end
    built.p.initialization.thermal_bottom_alignment = report;
end
