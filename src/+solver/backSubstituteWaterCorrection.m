function [deltaH, coupledRhs, linearResidual] = ...
        backSubstituteWaterCorrection(waterRow, deltaThermal, deltaP, activeRange)
%BACKSUBSTITUTEWATERCORRECTION  Re-solve the frozen water row after cross updates.
%
%   [deltaH, coupledRhs, linearResidual] = ...
%       solver.backSubstituteWaterCorrection(waterRow, deltaThermal, ...
%       deltaP, activeRange)
%
%   solves the total water-head correction from the frozen correction row
%
%       J_wh*deltaH = rhs_w - J_wq*deltaThermal - J_wc*deltaCross,
%
%   after the sequential dry-air and thermal coordinates have been updated.
%   The thermal coordinate is deltaU when waterRow contains an U bank and
%   deltaT otherwise. This function performs algebra only: it does not apply
%   head safeguards, boundary materialization, state mutation, or commit.
%
%   Inputs:
%     waterRow     correction row with h, P, rhs, and U or T banks
%     deltaThermal [U or degC, NN x 1] total thermal-coordinate correction
%     deltaP       [P_g or dry-air-mass coordinate, NN x 1] total cross
%                  correction selected by the row's P or M bank
%     activeRange  [-, 1 x 2] inclusive Richards solve range
%
%   Returns:
%     deltaH        [cm, NN x 1] total back-substituted head correction;
%                   zero outside activeRange
%     coupledRhs    [water row unit, NN x 1] right-hand side after cross terms
%     linearResidual[water row unit, NN x 1] coupledRhs-J_wh*deltaH
%
%   MATLAB V2.0 source:
%     src/+solver/attemptSolveSoilEquations.m,
%     solveInternalEnergyPrimaryMap correction rows.
%   Mathematical contract:
%     docs/soil/compass_soil_water_heat_physics.md, section 5.9.
%   Editor: Jizixing Zhu; Date: 2026-08-12.
%   Change: 2026-08-12 Jizixing Zhu -- accept the conservative dry-air mass
%           coordinate in addition to the legacy pressure coordinate.

    assert(isstruct(waterRow) && isfield(waterRow, 'h') ...
            && isfield(waterRow, 'rhs') ...
            && (isfield(waterRow, 'P') || isfield(waterRow, 'M')), ...
        'backSubstituteWaterCorrection:badRow', ...
        'waterRow must contain h, rhs, and either a P or M bank.');
    n = numel(waterRow.rhs);
    deltaThermal = deltaThermal(:);
    deltaP = deltaP(:);
    assert(numel(deltaThermal) == n && numel(deltaP) == n ...
            && all(isfinite(deltaThermal)) && all(isfinite(deltaP)), ...
        'backSubstituteWaterCorrection:badCorrection', ...
        'Cross-coordinate corrections must be finite NN-vectors.');
    if nargin < 4 || isempty(activeRange)
        activeRange = [1, n];
    end
    assert(numel(activeRange) == 2 ...
            && activeRange(1) >= 1 && activeRange(2) <= n ...
            && activeRange(1) <= activeRange(2) ...
            && all(activeRange == floor(activeRange)), ...
        'backSubstituteWaterCorrection:badActiveRange', ...
        'activeRange must be an inclusive integer range within 1:NN.');

    if isfield(waterRow, 'U')
        thermalBank = waterRow.U;
    elseif isfield(waterRow, 'T')
        thermalBank = waterRow.T;
    else
        error('backSubstituteWaterCorrection:missingThermalBank', ...
            'waterRow must contain either an U or T correction bank.');
    end

    if isfield(waterRow, 'M')
        crossBank = waterRow.M;
    else
        crossBank = waterRow.P;
    end
    coupledRhs = waterRow.rhs(:) ...
        - solver.applyTridiagonalBank(thermalBank, deltaThermal) ...
        - solver.applyTridiagonalBank(crossBank, deltaP);
    solved = solver.solveTridiagonal( ...
        waterRow.h.lower, waterRow.h.matrix(:, 1), ...
        waterRow.h.matrix(:, 2), coupledRhs, activeRange);
    deltaH = zeros(n, 1);
    first = activeRange(1);
    last = activeRange(2);
    deltaH(first:last) = solved(first:last);
    linearResidual = coupledRhs ...
        - solver.applyTridiagonalBank(waterRow.h, deltaH);
end
