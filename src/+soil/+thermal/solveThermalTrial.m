function [built, CHKHEAT, heatFlux, linearSystem, boundaryApplied] = solveThermalTrial( ...
        built, coefficients, matrix, energyAccepted, dt_s, kt, ...
        waterFlux, waterBoundaryFlux, assembleOnly, boundaryRequest)
%SOLVETHERMALTRIAL  Assemble and solve one heat equation on trial state.
%
%   [built, CHKHEAT, heatFlux] = soil.thermal.solveThermalTrial(built,
%       coefficients, matrix, energyAccepted, dt_s, kt, waterFlux,
%       waterBoundaryFlux)
%   consumes the already refreshed heat coefficient and FEM matrix banks,
%   assembles the heat RHS, applies heat-specific boundaries, solves with the
%   common Thomas kernel, and updates only trial temperature and SAVE.
%
%   Inputs:
%     built              [inout] model container; accepted Y remains unchanged
%     coefficients       [heat coefficient, struct] equations.heat
%     matrix             [assembled coefficient, struct] matrices.heat
%     energyAccepted     [J cm^-3, NL x 2] accepted endpoint internal energy
%     dt_s               [s] trial-substep duration
%     kt                 [-] outer forcing index
%     waterFlux          [struct] current water solution, including QMB/Throughfall [cm/s]
%     waterBoundaryFlux  [struct] current forcing, including EVAP [cm/s]
%
%   Returns:
%     built     [out] trial temperature and workspace.solver.SAVE updated
%     CHKHEAT   [degC, NN x 1] temperature iterate increment
%     heatFlux  [struct] QET/QEB [J cm-2 s-1] top/bottom heat flux
%     boundaryApplied [struct] mathematical BC imposed on each endpoint:
%                         fixed_temperature or fixed_total_energy_flux
%
%   Method:
%     Coefficients and FEM matrices are read-only inputs. Heat RHS, heat
%     boundaries, temperature safeguards, and flux recovery remain local to
%     the heat module; only A*x=b is shared with other scalar equations.
%
%   Reference:
%     soil.thermal.assembleCoefficientMatrices
%     soil.thermal.calculateBoundaryConditions
%     soil.thermal.calculateBoundaryFlux
%   Editor: Jizixing Zhu; Date: 2026-07-13.
%   Change: 2026-07-13 Jizixing Zhu -- use the actual internal-energy storage
%   difference and CTT only as the temperature-correction tangent.
%   Change: 2026-07-14 Jizixing Zhu -- expose an assembly-only mode and the
%           post-boundary temperature-coordinate linear system.
%   Change: 2026-07-15 Jizixing Zhu -- consume the immutable trial heat
%           boundary request; retain the local materializer only for direct
%           compatibility callers.
%   Change: 2026-07-15 Jizixing Zhu -- apply materialized surface-energy-
%           balance Gtot as a frozen Neumann flux; the production transaction
%           no longer reconstructs the surface balance during Picard.
%   Change: 2026-07-16 Jizixing Zhu -- return the applied Dirichlet/Neumann
%           heat condition separately from the requested physical policy.
%   Change: 2026-07-17 Jizixing Zhu -- expose free-drainage liquid internal-energy
%           export as a named coupled bottom Robin condition.
%   Change: 2026-07-18 Jizixing Zhu -- apply the invariant bare/snow tile-
%           aggregate top energy flux through the existing Neumann row.
%   Change: 2026-07-19 Jizixing Zhu -- apply the CLM-style affine surface
%           energy response as a stabilizing Robin row in the U-primary map.
%   Change: 2026-07-20 Jizixing Zhu -- rename the natural fixed-flux request
%           atmospheric while retaining the same Neumann matrix row.
%   Change: 2026-07-20 Jizixing Zhu -- use the current surface-domain water
%     candidate for pond heat terms without writing accepted Y.

    if nargin < 9 || isempty(assembleOnly); assembleOnly = false; end
    if nargin < 10; boundaryRequest = []; end

    assert(~isempty(fieldnames(coefficients)) && ~isempty(fieldnames(matrix)), ...
        'solveThermalTrial:missingCoefficients', ...
        'Heat coefficients and matrices are required when soil heat is enabled.');
    tr = built.workspace.step.trial;
    Y = built.Y;
    oldTemperature = tr.T(:);
    cfg = heatConfig(built);
    [boundary, boundaryApplied] = heatBoundary(built, kt, boundaryRequest);
    SV = struct('hh', tr.h(:), 'P_gg', tr.p_gas(:), ...
        'T', Y.T(:), 'Tbh', boundary.runtimeTbh);
    SAVE = built.workspace.solver.SAVE;
    conservativeState = struct('h', tr.h(:), 'p_gas', tr.p_gas(:));
    [RHS, matrix, SAVE] = soil.thermal.assembleConservativeEnergyEquation( ...
        matrix, conservativeState, cfg, dt_s, oldTemperature, ...
        coefficients.energyStorage, energyAccepted, SAVE);

    meteo = struct('Ta', built.forcing.hourly.Ta(kt));
    hBoundaryFlux = struct('QMB', waterFlux.QMB);
    localKt = 1;
    [RHS, matrix] = soil.thermal.calculateBoundaryConditions( ...
        boundary, matrix, hBoundaryFlux, meteo, SV, ...
        waterFlux.Throughfall, waterBoundaryFlux.EVAP, dt_s, RHS, ...
        built.workspace.derived.L, localKt);
    linearSystem = struct('matrix', matrix, 'rhs', RHS(:), 'SAVE', SAVE, ...
        'boundary', boundary);

    if assembleOnly
        solved = oldTemperature;
    else
        solved = solver.solveTridiagonal( ...
            matrix.C5_a, matrix.C5(:, 1), matrix.C5(:, 2), ...
            RHS, [1, built.domain.NN]);
    end

    if any(~isfinite(solved)) || ~isreal(solved) || any(abs(solved) >= 100)
        solved = Y.T(:);
    end
    solved(solved <= -272) = -272;
    CHKHEAT = abs(solved - oldTemperature);
    [QET, QEB] = soil.thermal.calculateBoundaryFlux(SAVE, solved, cfg);
    heatFlux = struct('QET', QET, 'QEB', QEB);
    if ~assembleOnly
        built.workspace.step.trial.T = reshape( ...
            solved, size(built.workspace.step.trial.T));
    end
    built.workspace.solver.SAVE = SAVE;
end

function cfg = heatConfig(built)
% Build the minimal immutable config read by heat RHS and boundary assembly.
% Editor: Jizixing Zhu; Date: 2026-07-14.
    cfg = struct();
    cfg.Grid = struct('NL', built.domain.NL, 'NN', built.domain.NN, ...
        'nD', built.domain.grid.nD, 'DeltZ', built.domain.dz);
    cfg.lThmrlefc = 1;
    cfg.lSoilairefc = readSwitch(built.control.switches, 'lSoilairefc');
end

function [boundary, applied] = heatBoundary(built, kt, request)
% Translate the typed request to the retained heat-kernel boundary frame.
% Editor: Jizixing Zhu; Date: 2026-07-14.
% Change: 2026-07-15 Jizixing Zhu -- use scheduler-materialized values when supplied.
    surfaceWater = surface.currentWaterState(built);
% Change: 2026-07-16 Jizixing Zhu -- expose the mathematical condition
%         actually imposed by the retained matrix adapter.
% Change: 2026-07-17 Jizixing Zhu -- accept the dynamic snow-insulated
%         Neumann request produced from a configured surface-temperature policy.
    if ~isempty(request)
        [boundary, applied] = requestToLegacyBoundary(request, built);
        return;
    end
    b = built.boundary;
    runtime = built.workspace.boundary.thermal;
    assert(size(runtime.Tbh, 1) >= kt && isfinite(runtime.Tbh(kt)), ...
        'solveThermalTrial:missingBottomTemperature', ...
        'Runtime bottom temperature Tbh must be refreshed before the heat solve.');
    topTemperature = b.T_top;
    if isfield(built.workspace, 'scope_step') ...
            && isfield(built.workspace.scope_step, 'thermal') ...
            && isfield(built.workspace.scope_step.thermal, 'Tss') ...
            && isfinite(built.workspace.scope_step.thermal.Tss)
        topTemperature = built.workspace.scope_step.thermal.Tss;
    end
    boundary = struct('NBCT', b.NBCT_top, 'NBCTB', b.NBCT_bottom, ...
        'BCT', topTemperature, 'BCTB', b.T_bottom, ...
        'BCT_ref_temperature', NaN, 'BCT_dflux_dtemperature', 0, ...
        'DSTOR0', surfaceWater.DSTOR0, ...
        'runtimeTbh', runtime.Tbh(kt));
    applied = inferLegacyApplied(boundary);
end

function [boundary, applied] = requestToLegacyBoundary(request, built)
% Convert public inward-positive values to legacy heat-kernel signs.
% Editor: Jizixing Zhu; Date: 2026-07-15.
% Change: 2026-07-16 Jizixing Zhu -- record the applied mathematical BC.
    assert(isstruct(request) && request.enabled, ...
        'solveThermalTrial:badBoundaryRequest', ...
        'Enabled heat solve requires an enabled typed boundary request.');
    switch request.top.kind
        case 'surface_temperature'
            topCode = 1; topValue = request.top.temperature_C;
        case {'specified_heat_flux','snow_insulated_heat_flux'}
            topCode = 2;
            topValue = -request.top.requested_flux_into_domain_J_cm2_s;
            topReferenceTemperature = NaN;
            topFluxDerivative = 0;
        case 'atmospheric'
            topCode = 2;
            topValue = -request.top.requested_flux_into_domain_J_cm2_s;
            topReferenceTemperature = NaN;
            topFluxDerivative = 0;
        case 'surface_energy_response'
            topCode = 3;
            topValue = request.top.requested_flux_into_domain_J_cm2_s;
            topReferenceTemperature = request.top.reference_temperature_C;
            topFluxDerivative = ...
                request.top.flux_temperature_derivative_J_cm2_s_K;
        otherwise
            error('solveThermalTrial:badTopBoundaryKind', ...
                'Unsupported heat top kind %s.', request.top.kind);
    end
    switch request.bottom.kind
        case {'seasonal_bottom_temperature','fixed_temperature'}
            bottomCode = 2; bottomValue = request.bottom.temperature_C;
        case 'specified_heat_flux'
            bottomCode = 3;
            bottomValue = request.bottom.requested_flux_into_domain_J_cm2_s;
        case 'free_drainage_internal_energy'
            bottomCode = 4;
            bottomValue = NaN;
        otherwise
            error('solveThermalTrial:badBottomBoundaryKind', ...
                'Unsupported heat bottom kind %s.', request.bottom.kind);
    end
    if topCode == 1
        topReferenceTemperature = NaN;
        topFluxDerivative = 0;
    end
    surfaceWater = surface.currentWaterState(built);
    boundary = struct('NBCT', topCode, 'NBCTB', bottomCode, ...
        'BCT', topValue, 'BCTB', bottomValue, ...
        'BCT_ref_temperature', topReferenceTemperature, ...
        'BCT_dflux_dtemperature', topFluxDerivative, ...
        'DSTOR0', surfaceWater.DSTOR0, ...
        'runtimeTbh', request.bottom.temperature_C);
    applied = struct( ...
        'top_kind', appliedKind(topCode, 'top'), ...
        'bottom_kind', appliedKind(bottomCode, 'bottom'));
end

function applied = inferLegacyApplied(boundary)
% Describe legacy heat codes using the public mathematical vocabulary.
% Editor: Jizixing Zhu; Date: 2026-07-16.
    applied = struct( ...
        'top_kind', appliedKind(boundary.NBCT, 'top'), ...
        'bottom_kind', appliedKind(boundary.NBCTB, 'bottom'));
end

function kind = appliedKind(code, side)
% Map one retained heat code to the equation imposed on the boundary row.
% Editor: Jizixing Zhu; Date: 2026-07-16.
    if strcmp(side, 'top')
        if code == 1
            kind = 'fixed_temperature';
        elseif code == 3
            kind = 'linearized_surface_energy_flux';
        else
            kind = 'fixed_total_energy_flux';
        end
    else
        if ismember(code, [1,2])
            kind = 'fixed_temperature';
        elseif code == 4
            kind = 'free_drainage_internal_energy';
        else
            kind = 'fixed_total_energy_flux';
        end
    end
end

function value = readSwitch(switches, name)
% Read a master switch without inventing an enabled default.
% Editor: Jizixing Zhu; Date: 2026-07-14.
    value = isfield(switches, name) && switches.(name) == 1;
end
