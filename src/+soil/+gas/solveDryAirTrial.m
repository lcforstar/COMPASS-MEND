function [built, CHKAIR, pressureCorrection, linearSystem, boundaryApplied] = solveDryAirTrial( ...
        built, coefficients, matrix, storageAccepted, dt_s, kt, assembleOnly, boundaryRequest)
%SOLVEDRYAIRTRIAL  Assemble and solve one dry-air equation on trial state.
%
%   [built, CHKAIR, pressureCorrection] = soil.gas.solveDryAirTrial(built, coefficients,
%       matrix, storageAccepted, dt_s, kt) consumes the already refreshed dry-air coefficient
%   and FEM matrix banks. It assembles the dry-air RHS, applies the dry-air
%   boundary policy, calls the common scalar Thomas kernel, and writes only
%   workspace.step.trial.p_gas plus the shared boundary-flux snapshot.
%
%   Inputs:
%     built        [inout] model container; Y is the accepted time layer and
%                  workspace.step.trial is the current Picard iterate
%     coefficients [dry-air coefficient, struct] equations.dry_air
%     matrix       [assembled coefficient, struct] matrices.dry_air
%     storageAccepted [g cm^-3, NL x 2] accepted dry-air storage
%     dt_s         [s] trial-substep duration
%     kt           [-] outer forcing index
%
%   Returns:
%     built   [out] trial gas pressure and workspace.solver.SAVE updated; Y unchanged
%     CHKAIR  [internal pressure unit, NN x 1] pressure iterate increment
%     pressureCorrection [internal pressure unit, NN x 1] signed raw Thomas
%                  fixed-point correction. The coupled Picard driver may
%                  relax this correction without changing the equation solve.
%     boundaryApplied [struct] mathematical BC imposed on each endpoint:
%                  fixed_pressure or fixed_mass_flux
%
%   Method:
%     Physical coefficients are not recomputed here. RHS and boundary physics
%     remain dry-air-specific; only A*x=b uses solver.solveTridiagonal.
%
%   Reference:
%     soil.gas.assembleCoefficientMatrices
%     soil.gas.calculateBoundaryConditions
%     soil.gas.solveDryAirEquations
%   Editor: Jizixing Zhu; Date: 2026-07-13.
%   Change: 2026-07-13 Jizixing Zhu -- use the actual dry-air storage change
%   and Caa only as the pressure-correction tangent.
%   Change: 2026-07-14 Jizixing Zhu -- expose the signed raw pressure
%   correction so the distributed Picard driver can stabilize the gas update
%   while retaining the unrelaxed fixed-point defect as its convergence test.
%   Change: 2026-07-14 Jizixing Zhu -- expose an assembly-only mode and the
%           post-boundary pressure-coordinate linear system.
%   Change: 2026-07-15 Jizixing Zhu -- consume the immutable typed dry-air
%           boundary request and confine Pa-to-internal conversion here.
%   Change: 2026-07-16 Jizixing Zhu -- return the applied pressure/mass-flux
%           condition separately from the requested physical policy.

    if nargin < 7 || isempty(assembleOnly); assembleOnly = false; end
    if nargin < 8; boundaryRequest = []; end

    assert(~isempty(fieldnames(coefficients)) && ~isempty(fieldnames(matrix)), ...
        'solveDryAirTrial:missingCoefficients', ...
        'Dry-air coefficients and matrices are required when soil air is enabled.');
    tr = built.workspace.step.trial;
    oldPressure = tr.p_gas(:);

    cfg = gasConfig(built);
    SAVE = built.workspace.solver.SAVE;
    conservativeState = struct('h', tr.h(:), 'T', tr.T(:));
    [RHS, matrix, SAVE] = soil.gas.assembleConservativeDryAirEquation( ...
        matrix, conservativeState, cfg, dt_s, oldPressure, ...
        coefficients.dryAirStorage, storageAccepted, SAVE);
    [boundary, boundaryApplied] = gasBoundary(built, kt, boundaryRequest);
    localKt = 1;
    [RHS, matrix] = soil.gas.calculateBoundaryConditions( ...
        boundary, matrix, RHS, localKt);
    linearSystem = struct('matrix', matrix, 'rhs', RHS(:), 'SAVE', SAVE, ...
        'boundary', boundary);
    if assembleOnly
        solved = oldPressure;
    else
        solved = solver.solveTridiagonal( ...
            matrix.C6_a(1:end-1), matrix.C6(:, 1), matrix.C6(:, 2), ...
            RHS, [1, built.domain.NN]);
    end

    invalid = ~isfinite(solved) | ~isreal(solved);
    solved(invalid) = oldPressure(invalid);
    pressureCorrection = solved - oldPressure;
    CHKAIR = abs(pressureCorrection);
    if ~assembleOnly
        built.workspace.step.trial.p_gas = reshape( ...
            solved, size(built.workspace.step.trial.p_gas));
    end
    built.workspace.solver.SAVE = SAVE;
end

function cfg = gasConfig(built)
% Build the minimal immutable config read by dry-air RHS assembly.
% Editor: Jizixing Zhu; Date: 2026-07-14.
    cfg = struct();
    cfg.Grid = struct('NL', built.domain.NL, 'NN', built.domain.NN, ...
        'nD', built.domain.grid.nD, 'DeltZ', built.domain.dz);
    cfg.lThmrlefc = readSwitch(built.control.switches, 'lThmrlefc');
end

function [boundary, applied] = gasBoundary(built, kt, request)
% Translate a typed trial request to the retained dry-air-kernel frame.
% Editor: Jizixing Zhu; Date: 2026-07-14.
% Change: 2026-07-15 Jizixing Zhu -- use scheduler-materialized values when supplied.
% Change: 2026-07-16 Jizixing Zhu -- expose the mathematical condition
%         actually imposed by the retained matrix adapter.
    if ~isempty(request)
        [boundary, applied] = requestToLegacyBoundary(request, built);
        return;
    end
    b = built.boundary;
    assert(isfield(built.forcing, 'hourly') && isfield(built.forcing.hourly, 'Pg') ...
        && numel(built.forcing.hourly.Pg) >= kt, ...
        'solveDryAirTrial:missingPressureForcing', ...
        'Hourly atmospheric pressure is required by the dry-air boundary.');
    topPg = built.forcing.hourly.Pg(kt) * 10;             % [g cm-1 s-2] Pa -> internal pressure
    boundary = struct('NBCP', b.NBCP_top, 'NBCPB', b.NBCP_bottom, ...
        'BCP', b.P_top, 'BCPB', b.P_bottom, 'TopPg', topPg, ...
        'BtmPg', b.BtmPg, 'UseHydrostaticAirBottom', b.UseHydrostaticAirBottom, ...
        'Tot_Depth', built.domain.grid.Tot_Depth, 'GroundwaterCoupling', 0);
    applied = struct( ...
        'top_kind', appliedKind(boundary.NBCP), ...
        'bottom_kind', appliedKind(boundary.NBCPB));
end

function [boundary, applied] = requestToLegacyBoundary(request, built)
% Convert Pa and inward-positive mass flux to legacy pressure/flux signs.
% Editor: Jizixing Zhu; Date: 2026-07-15.
% Change: 2026-07-16 Jizixing Zhu -- record the applied mathematical BC.
    assert(isstruct(request) && request.enabled, ...
        'solveDryAirTrial:badBoundaryRequest', ...
        'Enabled dry-air solve requires an enabled typed boundary request.');
    [topCode, topValue] = endpointToLegacy(request.top, 'top');
    [bottomCode, bottomValue] = endpointToLegacy(request.bottom, 'bottom');
    boundary = struct('NBCP', topCode, 'NBCPB', bottomCode, ...
        'BCP', topValue, 'BCPB', bottomValue, ...
        'TopPg', topValue, 'BtmPg', bottomValue, ...
        'UseHydrostaticAirBottom', false, ...
        'Tot_Depth', built.domain.grid.Tot_Depth, 'GroundwaterCoupling', 0);
    applied = struct( ...
        'top_kind', appliedKind(topCode), ...
        'bottom_kind', appliedKind(bottomCode));
end

function [code, value] = endpointToLegacy(endpoint, side)
% Map one typed endpoint to the legacy gas boundary convention.
    switch endpoint.kind
        case {'fixed_pressure','atmospheric_pressure','hydrostatic_pressure'}
            code = 1;
            value = endpoint.pressure_Pa * 10;  % Pa -> g cm^-1 s^-2
        case 'specified_mass_flux'
            code = 2;
            value = endpoint.requested_flux_into_domain_g_cm2_s;
            if strcmp(side, 'top'), value = -value; end
        case 'zero_mass_flux'
            code = 2; value = 0;
        otherwise
            error('solveDryAirTrial:badBoundaryKind', ...
                'Unsupported dry-air %s kind %s.', side, endpoint.kind);
    end
end

function kind = appliedKind(code)
% Map one retained gas code to the equation imposed on the boundary row.
% Editor: Jizixing Zhu; Date: 2026-07-16.
    if code == 1
        kind = 'fixed_pressure';
    else
        kind = 'fixed_mass_flux';
    end
end

function value = readSwitch(switches, name)
% Read a master switch without inventing an enabled default.
% Editor: Jizixing Zhu; Date: 2026-07-14.
    value = isfield(switches, name) && switches.(name) == 1;
end
