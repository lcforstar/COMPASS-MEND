function [RHS, AirMatrices, SAVE] = assembleConservativeDryAirEquation( ...
        AirMatrices, state, cfg, dt, pressureIteration, ...
        storageIteration, storageAccepted, SAVE)
%ASSEMBLECONSERVATIVEDRYAIREQUATION  Build the mixed dry-air mass RHS.
%
%   This routine evaluates the time-level residual with the actual dry-air
%   storage B_a=rho_a*(theta_g+Hc*theta_l). Caa is used only as the tangent
%   for P_g^(m+1)-P_g^m. Changes caused by h and T are already contained in
%   B_a^m-B_a^n, so Cah and CaT are not applied a second time in the time term.
%   Spatial h/T/P_g operators remain those assembled from the documented
%   frozen-coefficient Galerkin form.
%
%   Inputs use the STEMMUS cm-g-s system. storageIteration/storageAccepted
%   are dry-air mass densities [g cm^-3] at element endpoints; the returned
%   RHS is dry-air mass per area per time [g cm^-2 s^-1].
%
%   Reference: Zeng (2013), Coupled Dynamics in Soil, Eqs. 4.11-4.13 and 4.56.
%   Editor: Jizixing Zhu; Date: 2026-07-13.
%   Change: 2026-07-15 Jizixing Zhu -- snapshot the top-row strict-lower
%           coefficient; the preceding superdiagonal differs when Vaa is nonzero.
%   Change: 2026-07-16 Jizixing Zhu -- remove the artificial pressure
%           diagonal because it acted as an unphysical dry-air mass sink.

    NN = cfg.Grid.NN;
    dz = cfg.Grid.DeltZ(:);
    pressureIteration = pressureIteration(:);
    h = state.h(:);
    temperature = state.T(:);
    assert(isequal(size(storageIteration), size(storageAccepted)) ...
        && isequal(size(storageIteration), [NN - 1, 2]) ...
        && dt > 0 && all(isfinite(storageIteration(:))) ...
        && all(isfinite(storageAccepted(:))), ...
        'assembleConservativeDryAirEquation:badStorage', ...
        'Dry-air endpoint storage must be finite NL-by-2 arrays and dt positive.');
    assert(all(AirMatrices.C3(:,2) == 0), ...
        'assembleConservativeDryAirEquation:notLumped', ...
        'The conservative dry-air correction requires lumped C3 storage.');

    storageResidual = zeros(NN, 1);
    deltaStorage = storageIteration - storageAccepted;
    for e = 1:NN-1
        storageResidual(e) = storageResidual(e) ...
            + deltaStorage(e,1) * dz(e) / 2;
        storageResidual(e+1) = storageResidual(e+1) ...
            + deltaStorage(e,2) * dz(e) / 2;
    end

    spatialCross = applyTridiagonalBank(AirMatrices.C4, AirMatrices.C4_a, h);
    if cfg.lThmrlefc == 1
        spatialCross = spatialCross + applyTridiagonalBank( ...
            AirMatrices.C5, AirMatrices.C5_a, temperature);
    end
    RHS = -AirMatrices.C7 - spatialCross ...
        + AirMatrices.C3(:,1) .* pressureIteration / dt ...
        - storageResidual / dt;
    AirMatrices.C6 = AirMatrices.C6 + AirMatrices.C3 / dt;

    SAVE(1,1,3) = RHS(1);
    SAVE(1,2,3) = AirMatrices.C6(1,1);
    SAVE(1,3,3) = AirMatrices.C6(1,2);
    SAVE(2,1,3) = RHS(NN);
    SAVE(2,2,3) = AirMatrices.C6_a(NN-1);
    SAVE(2,3,3) = AirMatrices.C6(NN,1);
end

function y = applyTridiagonalBank(matrix, lower, x)
% Apply the repository's main/super plus strict-lower matrix storage.
    n = numel(x);
    y = matrix(:,1) .* x;
    y(1:n-1) = y(1:n-1) + matrix(1:n-1,2) .* x(2:n);
    y(2:n) = y(2:n) + lower(1:n-1) .* x(1:n-1);
end
