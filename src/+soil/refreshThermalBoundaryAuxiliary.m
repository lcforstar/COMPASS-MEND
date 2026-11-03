function built = refreshThermalBoundaryAuxiliary(built, KT, NoTime_KT)
%REFRESHTHERMALBOUNDARYAUXILIARY  Refresh the runtime bottom-temperature boundary.
%
%   built = soil.refreshThermalBoundaryAuxiliary(built, KT, NoTime_KT)
%   evaluates the retained deep-soil harmonic boundary model after thermal
%   transport properties have been refreshed. It writes only the recomputable
%   runtime boundary cache and performs no equation-coefficient or matrix
%   assembly.
%
%   Inputs:
%     built     [inout] model container with static thermal-bottom parameters
%                       and refreshed Lambda_eff/c_unsat auxiliary fields
%     KT        [-] current model-step index
%     NoTime_KT [-] elapsed hourly-step index used by the daily harmonic model
%
%   Returns:
%     built.workspace.boundary.thermal.Tbh [degC, tS x nCell]
%         hourly bottom-temperature boundary values
%     built.workspace.boundary.thermal.Tb [degC, nDay x nCell]
%         daily deep-soil temperature values
%
%   Method:
%     Static Tav/At remain in built.boundary.thermal_bottom. Runtime Tbh/Tb
%     belong to workspace.boundary because they are recomputable boundary
%     auxiliaries, not equation coefficients or accepted prognostic state.
%
%   Reference:
%     soil.thermal.calculateBottomTempBoundary
%   Editor: Jizixing Zhu; Date: 2026-07-12.

    if nargin < 2 || isempty(KT), KT = 1; end
    if nargin < 3 || isempty(NoTime_KT), NoTime_KT = KT; end

    nCell = readDomainDim(built.domain, 'nCell', 1);
    thermal = built.workspace.boundary.thermal;
    nDay = numel(equations.generateJDSeries( ...
        built.control.time.start_time, built.control.time.end_time));
    thermal = ensureCacheShape(thermal, built.control.time.tS, nDay, nCell);
    staticBoundary = readStaticBoundary(built);

    for cellIndex = 1:nCell
        SV = struct('Tb', thermal.Tb(:, cellIndex));
        if ~isempty(thermal.Tb_jd_series)
            SV.Tb_jd_series = thermal.Tb_jd_series;
            SV.Tb_JD0 = thermal.Tb_JD0;
        end
        TC = struct( ...
            'Lambda_eff', endpointColumn(built.workspace.derived.Lambda_eff, cellIndex), ...
            'c_unsat', endpointColumn(built.workspace.derived.c_unsat, cellIndex));
        cfg = struct();
        cfg.Grid = struct('Tot_Depth', built.domain.grid.Tot_Depth);
        cfg.Env = struct( ...
            'Tav', scalarColumn(staticBoundary.Tav, cellIndex), ...
            'At', scalarColumn(staticBoundary.At, cellIndex));
        cfg.Time = struct( ...
            'start_time', built.control.time.start_time, ...
            'end_time', built.control.time.end_time, ...
            'Delt_t', built.control.time.Delt_t);

        SV = soil.thermal.calculateBottomTempBoundary(SV, TC, cfg, KT, NoTime_KT);
        thermal.Tbh(KT, cellIndex) = SV.Tbh(KT);
        thermal.Tb(:, cellIndex) = SV.Tb(:);
        thermal.Tb_jd_series = SV.Tb_jd_series;
        thermal.Tb_JD0 = SV.Tb_JD0;
    end
    built.workspace.boundary.thermal = thermal;
end

function boundary = readStaticBoundary(built)
% Read static harmonic parameters, with an inert fallback when heat is disabled.
    if isfield(built.boundary, 'thermal_bottom')
        boundary = built.boundary.thermal_bottom;
        return;
    end
    heatEnabled = isfield(built.control.switches, 'lThmrlefc') ...
        && built.control.switches.lThmrlefc == 1;
    assert(~heatEnabled, ...
        'refreshThermalBoundaryAuxiliary:missingThermalBottom', ...
        'built.boundary.thermal_bottom is required when soil heat is enabled.');
    boundary = struct('Tav', mean(built.Y.T, 1), 'At', zeros(1, size(built.Y.T, 2)));
end

function thermal = ensureCacheShape(thermal, tS, nDay, nCell)
% Allocate fixed-shape runtime boundary arrays when initialization left them empty.
    if ~isequal(size(thermal.Tbh), [tS, nCell])
        thermal.Tbh = zeros(tS, nCell);
    end
    if ~isequal(size(thermal.Tb), [nDay, nCell])
        thermal.Tb = zeros(nDay, nCell);
    end
end

function A = endpointColumn(A, cellIndex)
% Select one regional endpoint column; point mode is already two-dimensional.
    if ndims(A) >= 3
        A = A(:, :, cellIndex);
    end
end

function value = scalarColumn(value, cellIndex)
% Select one regional scalar parameter while preserving point-site scalars.
    if ~isscalar(value)
        value = value(cellIndex);
    end
end

function value = readDomainDim(domain, name, defaultValue)
% Read one optional domain dimension.
    value = defaultValue;
    if isfield(domain, name) && ~isempty(domain.(name))
        value = domain.(name);
    end
end
