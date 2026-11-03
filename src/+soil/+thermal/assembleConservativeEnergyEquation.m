function [RHS, EnergyMatrices, SAVE] = assembleConservativeEnergyEquation( ...
        EnergyMatrices, state, cfg, dt, temperatureIteration, ...
        energyIteration, energyAccepted, SAVE)
%ASSEMBLECONSERVATIVEENERGYEQUATION  Build the mixed internal-energy RHS.
%
%   [RHS, EnergyMatrices, SAVE] =
%   soil.thermal.assembleConservativeEnergyEquation(...) forms one Picard
%   correction for the heat equation while retaining temperature as the local
%   scalar Thomas unknown. The nonlinear storage residual is evaluated from
%   the actual endpoint internal energies:
%
%     R_U^m = integral(N_i*(U^m-U^n) dz).
%
%   CTT is used only as the local tangent multiplying T^(m+1)-T^m. CTh and
%   CTa are therefore not used to approximate the time-level energy change;
%   their effects are already present in U^m-U^n. This is the energy analogue
%   of a modified-Picard storage residual and prevents phase energy from being
%   counted once through U and again through h/T cross-capacities.
%
%   Inputs:
%     EnergyMatrices [struct] assembled C1-C9 heat matrix bank
%     state.h [cm] and state.p_gas [internal pressure], NN-by-1 trial state
%     cfg.Grid.DeltZ [cm], cfg.lSoilairefc [-]
%     dt [s]
%     temperatureIteration [degC, NN-by-1] T^m
%     energyIteration/energyAccepted [J cm^-3, NL-by-2] U^m/U^n
%     SAVE [legacy boundary snapshot]
%
%   Returns RHS in J cm^-2 s^-1 and adds C2/dt to the T-column matrix C5.
%
%   Reference: docs/compass_soil_water_heat_physics.md, Eqs. 4.3.13 and 5.18.
%   Editor: Jizixing Zhu; Date: 2026-07-13.
%   Change: 2026-07-15 Jizixing Zhu -- snapshot the top-row strict-lower
%           coefficient; the neighboring superdiagonal is unequal when VTT is nonzero.

    NN = cfg.Grid.NN;
    dz = cfg.Grid.DeltZ(:);
    temperatureIteration = temperatureIteration(:);
    h = state.h(:);
    pGas = state.p_gas(:);
    assert(isequal(size(energyIteration), size(energyAccepted)) ...
        && isequal(size(energyIteration), [NN - 1, 2]) ...
        && dt > 0 && all(isfinite(energyIteration(:))) ...
        && all(isfinite(energyAccepted(:))), ...
        'assembleConservativeEnergyEquation:badStorage', ...
        'Energy endpoint arrays must be finite NL-by-2 values and dt must be positive.');
    assert(all(EnergyMatrices.C2(:,2) == 0), ...
        'assembleConservativeEnergyEquation:notLumped', ...
        'The conservative energy correction requires lumped C2 storage.');

    energyResidual = zeros(NN, 1);
    deltaEnergy = energyIteration - energyAccepted;
    for e = 1:NN-1
        energyResidual(e) = energyResidual(e) ...
            + deltaEnergy(e,1) * dz(e) / 2;
        energyResidual(e+1) = energyResidual(e+1) ...
            + deltaEnergy(e,2) * dz(e) / 2;
    end

    spatialCross = applyTridiagonalBank( ...
        EnergyMatrices.C4, EnergyMatrices.C4_a, h);
    if cfg.lSoilairefc == 1
        spatialCross = spatialCross + applyTridiagonalBank( ...
            EnergyMatrices.C6, EnergyMatrices.C6_a, pGas);
    end

    RHS = -EnergyMatrices.C7 - EnergyMatrices.C9 - spatialCross ...
        + EnergyMatrices.C2(:,1) .* temperatureIteration / dt ...
        - energyResidual / dt;
    EnergyMatrices.C5 = EnergyMatrices.C5 + EnergyMatrices.C2 / dt;

    SAVE(1,1,2) = RHS(1);
    SAVE(1,2,2) = EnergyMatrices.C5(1,1);
    SAVE(1,3,2) = EnergyMatrices.C5(1,2);
    SAVE(2,1,2) = RHS(NN);
    SAVE(2,2,2) = EnergyMatrices.C5_a(NN-1);
    SAVE(2,3,2) = EnergyMatrices.C5(NN,1);
end

function y = applyTridiagonalBank(matrix, lower, x)
% Apply the repository's main/super plus strict-lower matrix storage.
    n = numel(x);
    y = matrix(:,1) .* x;
    y(1:n-1) = y(1:n-1) + matrix(1:n-1,2) .* x(2:n);
    y(2:n) = y(2:n) + lower(1:n-1) .* x(1:n-1);
end
