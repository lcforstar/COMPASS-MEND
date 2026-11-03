function [WaterMatrix, RHS, SAVE] = calculateTimeDerivatives( ...
        WaterMatrix, SoilVariables, P_gg, Delt_t, cfg, SAVE, thetaEqN)
%CALCULATETIMEDERIVATIVES  Assemble the mixed-form water-equation RHS.
%
%   [WaterMatrix, RHS, SAVE] = soil.water.calculateTimeDerivatives(
%   WaterMatrix, SoilVariables, P_gg, Delt_t, cfg, SAVE, thetaEqN)
%   assembles the common water equation first, then adds thermal and gas
%   coupling terms under their respective switches.
%
%   Inputs:
%     WaterMatrix   [inout] assembled finite-element coefficients C1-C9
%     SoilVariables [in] current Picard and accepted water/thermal states
%     P_gg          [internal pressure unit, NN x 1] current gas pressure
%     Delt_t        [s] positive implicit time-step length
%     cfg           [in] grid definition and heat/gas equation switches
%     SAVE          [inout] boundary equations used for flux recovery
%     thetaEqN      [-, NL x 2] accepted total water-equivalent content at
%                   both endpoints of every element
%
%   Returns:
%     WaterMatrix   [out] C4 with C1/Delt_t added to the time diagonal
%     RHS           [out, NN x 1] assembled known side of the water equation
%     SAVE          [out] updated bottom/top equation snapshots
%
%   Coefficient roles:
%     C1  water-head storage Jacobian
%     C2  temperature storage Jacobian
%     C4  water-head spatial operator
%     C5  temperature-driven water-flux operator
%     C6  gas-pressure-driven water-flux operator
%     C7  gravity term
%     C9  root/source-sink term
%
%   The common mixed time contribution is
%
%     [C1*h^m - Rtheta^m] / Delt_t,
%
%   where Rtheta^m is assembled from the exact endpoint change between the
%   accepted state n and current Picard state m. When heat coupling is active,
%   the original C2*(T^n-T^m)/Delt_t storage term and C5*T^m spatial term are
%   retained. Gas coupling retains the original C6*Pg^m spatial term.
%
%   References:
%     Zeng 2013, Coupled Dynamics in Soil, Eqs. 4.47-4.60.
%     Celia et al. 1990 Water Resources Research 26:1483-1496,
%     DOI 10.1029/WR026i007p01483 (modified Picard formulation).
%   Editor: Jizixing Zhu; Date: 2026-07-10.
%   Change: 2026-07-12 Jizixing Zhu -- make mixed storage the sole water-head
%           time formulation while retaining all original coupled C2/C5/C6 terms.
%   Change: 2026-07-12 Jizixing Zhu -- move thetaEqM/Rtheta construction into
%           a local function.
%   Change: 2026-07-12 Jizixing Zhu -- organize assembly as common terms,
%           optional thermal terms, and optional gas terms; each block is
%           ordered surface, interior, then bottom.
%   Change: 2026-07-13 Jizixing Zhu -- prevent thermal storage from being
%           counted twice in coupled mixed form: retain the explicit C2 block
%           but remove its identical linearized increment from Rtheta first.
%   Change: 2026-08-11 Jizixing Zhu -- save the actual top-row strict-lower
%           water coefficient C4_a(NL), not the preceding row's superdiagonal,
%           so nonsymmetric vapor transport yields the correct reaction flux.

    assert(isscalar(Delt_t) && isfinite(Delt_t) && Delt_t > 0, ...
        'calculateTimeDerivatives:badTimeStep', ...
        'Delt_t must be a finite positive scalar [s].');
    assert(nargin >= 7 && ~isempty(thetaEqN), ...
        'calculateTimeDerivatives:missingThetaEqN', ...
        'The mixed water equation requires accepted-state thetaEqN.');

    NL = cfg.Grid.NL;
    NN = cfg.Grid.NN;
    indxBotm = cfg.Grid.indxBotmLayer;

    C1 = WaterMatrix.C1;
    C2 = WaterMatrix.C2;
    C4 = WaterMatrix.C4;
    C4_a = WaterMatrix.C4_a;
    C5 = WaterMatrix.C5;
    C6 = WaterMatrix.C6;
    C7 = WaterMatrix.C7;
    C5_a = WaterMatrix.C5_a;
    C9 = WaterMatrix.C9;

    hh = SoilVariables.hh(:);  % [cm] current Picard water head h^m
    TT = SoilVariables.TT(:);  % [K] current Picard temperature T^m
    T = SoilVariables.T(:);    % [K] accepted temperature T^n
    P_gg = P_gg(:);            % [internal pressure unit] current Pg^m

    assert(all(C1(:, 2) == 0), ...
        'calculateTimeDerivatives:notLumped', ...
        'The mixed form requires lumped C1 with C1(:,2)==0.');
    Rtheta = calculateMixedStorageResidual(SoilVariables, thetaEqN, cfg);
    if cfg.lThmrlefc
        assert(all(C2(:, 2) == 0), ...
            'calculateTimeDerivatives:notLumpedThermalStorage', ...
            'The mixed thermal correction requires lumped C2 with C2(:,2)==0.');
        % Rtheta contains the exact total storage change, including the
        % temperature contribution. The explicit C2 block below contributes
        % -C2*(T^m-T^n)/dt. Subtract that same linearized increment from
        % Rtheta here so their sum is exactly -(thetaEq^m-thetaEq^n)/dt,
        % rather than counting temperature storage twice.
        thermalStorageIncrement = C2(:, 1) .* (TT - T);
        Rtheta = Rtheta - thermalStorageIncrement;
    end

    RHS = zeros(NN, 1);

    % ---------------------------------------------------------------------
    % 1. Common water equation: mixed storage + gravity + source/sink
    % Order inside every block: surface, interior, bottom.
    % ---------------------------------------------------------------------

    % Surface node (NN).
    RHS(NN) = -C9(NN) - C7(NN) ...
        + (C1(NN - 1, 2) * hh(NN - 1) ...
        + C1(NN, 1) * hh(NN) - Rtheta(NN)) / Delt_t;

    % Interior nodes, traversed from shallow to deep.
    for ML = NL:-1:(indxBotm + 1)
        RHS(ML) = -C9(ML) - C7(ML) ...
            + (C1(ML - 1, 2) * hh(ML - 1) ...
            + C1(ML, 1) * hh(ML) ...
            + C1(ML, 2) * hh(ML + 1) ...
            - Rtheta(ML)) / Delt_t;
    end

    % Bottom node (indxBotm).
    RHS(indxBotm) = -C9(indxBotm) - C7(indxBotm) ...
        + (C1(indxBotm, 1) * hh(indxBotm) ...
        + C1(indxBotm, 2) * hh(indxBotm + 1) ...
        - Rtheta(indxBotm)) / Delt_t;

    % ---------------------------------------------------------------------
    % 2. Thermal coupling: add C2 temperature storage and C5 thermal flux
    % ---------------------------------------------------------------------
    if cfg.lThmrlefc
        % Surface node (NN).
        RHS(NN) = RHS(NN) ...
            - (C2(NN - 1, 2) / Delt_t + C5_a(NN - 1)) * TT(NN - 1) ...
            - (C2(NN, 1) / Delt_t + C5(NN, 1)) * TT(NN) ...
            + C2(NN - 1, 2) * T(NN - 1) / Delt_t ...
            + C2(NN, 1) * T(NN) / Delt_t;

        % Interior nodes, traversed from shallow to deep.
        for ML = NL:-1:(indxBotm + 1)
            C2_lower = C2(ML - 1, 2) / Delt_t;
            C2_diag = C2(ML, 1) / Delt_t;
            C2_upper = C2(ML, 2) / Delt_t;
            RHS(ML) = RHS(ML) ...
                - (C2_lower + C5_a(ML - 1)) * TT(ML - 1) ...
                - (C2_diag + C5(ML, 1)) * TT(ML) ...
                - (C2_upper + C5(ML, 2)) * TT(ML + 1) ...
                + C2_lower * T(ML - 1) ...
                + C2_diag * T(ML) ...
                + C2_upper * T(ML + 1);
        end

        % Bottom node (indxBotm).
        RHS(indxBotm) = RHS(indxBotm) ...
            - (C2(indxBotm, 1) / Delt_t + C5(indxBotm, 1)) * TT(indxBotm) ...
            - (C2(indxBotm, 2) / Delt_t + C5(indxBotm, 2)) * TT(indxBotm + 1) ...
            + C2(indxBotm, 1) * T(indxBotm) / Delt_t ...
            + C2(indxBotm, 2) * T(indxBotm + 1) / Delt_t;
    end

    % ---------------------------------------------------------------------
    % 3. Gas coupling: add C6 gas-pressure-driven water flux
    % ---------------------------------------------------------------------
    if cfg.lSoilairefc
        % Surface node (NN).
        RHS(NN) = RHS(NN) ...
            - C6(NN - 1, 2) * P_gg(NN - 1) ...
            - C6(NN, 1) * P_gg(NN);

        % Interior nodes, traversed from shallow to deep.
        for ML = NL:-1:(indxBotm + 1)
            RHS(ML) = RHS(ML) ...
                - C6(ML - 1, 2) * P_gg(ML - 1) ...
                - C6(ML, 1) * P_gg(ML) ...
                - C6(ML, 2) * P_gg(ML + 1);
        end

        % Bottom node (indxBotm).
        RHS(indxBotm) = RHS(indxBotm) ...
            - C6(indxBotm, 1) * P_gg(indxBotm) ...
            - C6(indxBotm, 2) * P_gg(indxBotm + 1);
    end

    % Add the fully implicit water-head time term to the left-hand matrix.
    for MN = indxBotm:NN
        C4(MN, 1) = C4(MN, 1) + C1(MN, 1) / Delt_t;
        C4(MN, 2) = C4(MN, 2) + C1(MN, 2) / Delt_t;
    end

    % Save the two boundary rows for water-flux reconstruction.
    SAVE(1, 1, 1) = RHS(indxBotm);
    SAVE(1, 2, 1) = C4(indxBotm, 1);
    SAVE(1, 3, 1) = C4(indxBotm, 2);
    SAVE(2, 1, 1) = RHS(NN);
    SAVE(2, 2, 1) = C4_a(NN - 1);
    SAVE(2, 3, 1) = C4(NN, 1);

    WaterMatrix.C4 = C4;
end

function Rtheta = calculateMixedStorageResidual(SoilVariables, thetaEqN, cfg)
%CALCULATEMIXEDSTORAGERESIDUAL  Assemble accepted-to-current water storage.
%
%   This local function contains both endpoint conversion and lumped residual
%   assembly. Keeping them here keeps the main equation focused on common,
%   thermal, and gas terms without adding small public files to +soil/+water.
%   Editor: Jizixing Zhu; Date: 2026-07-12.

    RHOL = equations.getPhysicalConstant('RHOL');
    Theta_UU = SoilVariables.Theta_UU;
    Theta_V = SoilVariables.Theta_V;
    RHOV = SoilVariables.RHOV(:);
    DeltZ = cfg.Grid.DeltZ(:);
    NL = cfg.Grid.NL;
    NN = cfg.Grid.NN;

    assert(isequal(size(Theta_UU), [NL, 2]) ...
        && isequal(size(Theta_V), [NL, 2]) ...
        && isequal(size(thetaEqN), [NL, 2]), ...
        'calculateTimeDerivatives:badStorageShape', ...
        'Theta_UU, Theta_V, and thetaEqN must be NL-by-2 endpoint arrays.');
    assert(numel(RHOV) == NN && numel(DeltZ) == NL && NN == NL + 1, ...
        'calculateTimeDerivatives:badStorageGrid', ...
        'Mixed storage requires NN=NL+1, NN vapor nodes, and NL layer widths.');
    assert(isfinite(RHOL) && RHOL > 0 ...
        && all(isfinite(Theta_UU(:))) && all(isfinite(Theta_V(:))) ...
        && all(isfinite(thetaEqN(:))) && all(isfinite(RHOV)) ...
        && all(isfinite(DeltZ)), ...
        'calculateTimeDerivatives:badStorageValue', ...
        'Mixed-storage inputs must be finite and RHOL must be positive.');

    % Zeng 2013, Coupled Dynamics in Soil, Eqs. 4.57-4.60.
    % Total water equivalent at both endpoints of each element:
    % thetaEq = condensed water + pore vapor / liquid-water density.
    RHOV_endpoints = [RHOV(1:NL), RHOV(2:NN)];
    thetaEqM = Theta_UU + Theta_V .* RHOV_endpoints / RHOL;

    % Celia et al. 1990 WRR 26:1483-1496, modified Picard formulation.
    % Apply the same lumped dz/2 endpoint weights used to assemble C1.
    endpointChange = thetaEqM - thetaEqN;
    Rtheta = zeros(NN, 1);
    for ML = 1:NL
        halfLayer = DeltZ(ML) / 2;
        Rtheta(ML) = Rtheta(ML) ...
            + halfLayer * endpointChange(ML, 1);
        Rtheta(ML + 1) = Rtheta(ML + 1) ...
            + halfLayer * endpointChange(ML, 2);
    end
end
