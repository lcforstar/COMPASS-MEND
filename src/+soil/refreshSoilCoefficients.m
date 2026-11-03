function [built, equations] = refreshSoilCoefficients(built, opts)
%REFRESHSOILCOEFFICIENTS  Refresh soil auxiliary values and build equation coefficients.
%
%   [built, equations] = soil.refreshSoilCoefficients(built, opts) refreshes
%   water/phase and thermal-gas auxiliary fields, then calls
%   soil.buildEquationCoefficients once. The returned equations struct contains
%   water, dry-air, and heat coefficient banks only. This function performs no
%   matrix assembly, RHS construction, boundary application, or solve and does
%   not write accepted state in built.Y.
%
%   Inputs:
%     built [inout] model container; workspace.derived/flux are refreshed
%     opts  [in] coefficient options, including dt_s [s], L_f [J/kg],
%                Tmin [degC], DSTMAX [cm], Srt [1/s], airHeatImpedance [-],
%                KT/NoTime_KT [-] for the runtime bottom-temperature boundary
%
%   Returns:
%     built     [out] model container with refreshed auxiliary fields only
%     equations [out] struct with .water, .dry_air, and .heat coefficients
%
%   Reference: docs/superpowers/specs/2026-07-07-stage4-refresh-soil-coefficients-orchestrator-design.md.
%   Oracle: tests/test_refreshSoilCoefficients_orchestrator.m.
%   Editor: Jizixing Zhu; Date: 2026-07-07.
%   Change: 2026-07-12 Jizixing Zhu -- make the conservative water coefficient
%           bank the sole water path while retaining the dry-air and thermal
%           refresh implementations.
%   Change: 2026-07-12 Jizixing Zhu -- replace B3a/B3b adapters with the
%           coefficient-only soil.buildEquationCoefficients interface.
%   Change: 2026-07-12 Jizixing Zhu -- refresh Tbh only when the soil-heat
%           equation is enabled; equation banks obey their master switches.
%   Change: 2026-07-15 Jizixing Zhu -- allow the transactional solver to skip
%           Tbh refresh after a heat boundary request has been materialized;
%           direct callers retain the previous default behavior.
    if nargin < 2 || isempty(opts), opts = struct(); end
    opts = normalizeRefreshOptions(built, opts);

    built = soil.refreshWaterAux(built, opts.dt_s, opts.L_f, opts.Tmin);            % B1
    built = soil.refreshThermalGasTransport(built, opts.dt_s);                      % B2
    if switchEnabled(built.control.switches, 'lThmrlefc') ...
            && opts.refreshThermalBoundary
        built = soil.refreshThermalBoundaryAuxiliary(built, opts.KT, opts.NoTime_KT);
    end
    equations = soil.buildEquationCoefficients(built, opts);
end

function enabled = switchEnabled(switches, name)
% Read one master equation switch without inventing an enabled default.
    enabled = isfield(switches, name) && switches.(name) == 1;
end

function opts = normalizeRefreshOptions(built, opts)
% Fill defaults only; physical coefficients are built after B1/B2 refresh.
    d = struct('dt_s', 0, 'L_f', built.workspace.step.L_f, 'Tmin', 0, 'DSTMAX', 0, ...
        'Srt', built.workspace.flux.Srt, 'airHeatImpedance', 0, ...
        'KT', 1, 'NoTime_KT', 1, 'refreshThermalBoundary', true);
    f = fieldnames(d);
    for i = 1:numel(f)
        if ~isfield(opts, f{i}) || isempty(opts.(f{i}))
            opts.(f{i}) = d.(f{i});
        end
    end
end
