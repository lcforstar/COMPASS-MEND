function [RHS, EnergyMatrices, SAVE] = assembleCoefficientMatrices( ...
    EnergyMatrices, SoilVariables, cfg, Delt_t, T, h, P_g, SAVE)
% assembleCoefficientMatrices - 能量方程右端项组装，更新 C5 主对角并保存快照
%                               Assemble RHS for energy PDE, update C5, save SAVE snapshot
%
% 函数定义 / Function Definition:
%   [RHS, EnergyMatrices, SAVE] = assembleCoefficientMatrices( ...
%       EnergyMatrices, SoilVariables, cfg, Delt_t, T, h, P_g, SAVE)
%
% 输入 / Inputs:
%   EnergyMatrices - 矩阵系数结构体 / Matrix coefficient struct (from calculateMatricCoefficients)
%     .C1/.C2/.C3  - 储存矩阵 (NN×2)
%     .C4/.C5/.C6  - 导度+对流矩阵 (NN×2)
%     .C4_a/.C5_a/.C6_a - 下对角 (NL×1)
%     .C7          - 重力源向量 (NN×1)
%   SoilVariables  - 土壤状态结构体 / Soil state struct
%     .hh          - 当前 Newton 步水头 (NN×1) [cm]
%     .P_gg        - 当前 Newton 步气压 (NN×1) [Pa]
%   cfg            - 模型设置 / Model settings (lSoilairefc, lThmrlefc, Grid.NL/NN)
%   Delt_t         - 时步长 [s] / Time step size
%   T              - 上一 Newton 步温度 (NN×1) [K] / Previous Newton-step temperature
%   h              - 上一 Newton 步水头 (NN×1) [cm] / Previous Newton-step head
%   P_g            - 上一时步气压 (NN×1) [Pa] / Previous time-step pressure
%   SAVE           - 系数快照 (3×3×5) / Coefficient snapshot
%
% 输出 / Outputs:
%   RHS            - 右端向量 (NN×1) / Right-hand side vector
%   EnergyMatrices - 更新后结构体（C5 已加 C2/Delt_t）/ Updated struct (C5 modified)
%   SAVE           - 更新后快照 / Updated snapshot
%
% 引用 / References:
%   原 Enrgy_EQ.m；STEMMUS Technical Notes Section 5
%
% 修改者 / Modified by:   朱寄子星 (Zhu Jizixing)
% 修改日期 / Modified date: 2026-04-25
% Editor: Jizixing Zhu; Date: 2026-04-25.
% Change: 2026-07-13 Jizixing Zhu -- subtract the consistently assembled
% root-energy sink from the heat-equation right-hand side.
% Change: 2026-07-15 Jizixing Zhu -- snapshot C5_a(NL) for the top heat row;
% the preceding superdiagonal is not interchangeable when VTT is nonzero.

    NL    = cfg.Grid.NL;
    NN    = NL + 1;

    hh  = SoilVariables.hh(:);
    P_gg = SoilVariables.P_gg(:);
    T    = T(:);
    h    = h(:);
    P_g  = P_g(:);

    C1   = EnergyMatrices.C1;
    C2   = EnergyMatrices.C2;
    C3   = EnergyMatrices.C3;
    C4   = EnergyMatrices.C4;
    C5   = EnergyMatrices.C5;
    C6   = EnergyMatrices.C6;
    C4_a = EnergyMatrices.C4_a;
    C6_a = EnergyMatrices.C6_a;
    C7   = EnergyMatrices.C7;
    if isfield(EnergyMatrices, 'C9')
        C9 = EnergyMatrices.C9;
    else
        % Archived fixtures may predate the explicit root-energy sink.
        C9 = zeros(NN, 1);
    end

    RHS = zeros(NN, 1);

    if cfg.lSoilairefc && cfg.lThmrlefc
        % --- 底部节点 / Bottom node ---
        RHS(1) = -C7(1) - C9(1) + (C2(1,1)*T(1) + C2(1,2)*T(2)) / Delt_t ...
            - (C1(1,1)/Delt_t + C4(1,1))*hh(1) - (C1(1,2)/Delt_t + C4(1,2))*hh(2) ...
            - (C3(1,1)/Delt_t + C6(1,1))*P_gg(1) - (C3(1,2)/Delt_t + C6(1,2))*P_gg(2) ...
            + (C3(1,1)/Delt_t)*P_g(1) + (C3(1,2)/Delt_t)*P_g(2) ...
            + (C1(1,1)/Delt_t)*h(1)   + (C1(1,2)/Delt_t)*h(2);

        % --- 内部节点 / Interior nodes ---
        for ML = 2:NL
            A1 = C3(ML-1,2)/Delt_t;   A2 = C3(ML,1)/Delt_t;   A3 = C3(ML,2)/Delt_t;
            A4 = C1(ML-1,2)/Delt_t;   A5 = C1(ML,1)/Delt_t;   A6 = C1(ML,2)/Delt_t;

            RHS(ML) = -C7(ML) - C9(ML) + (C2(ML-1,2)*T(ML-1) + C2(ML,1)*T(ML) + C2(ML,2)*T(ML+1)) / Delt_t ...
                - (A1 + C6_a(ML-1))*P_gg(ML-1) - (A2 + C6(ML,1))*P_gg(ML) - (A3 + C6(ML,2))*P_gg(ML+1) ...
                - (A4 + C4_a(ML-1))*hh(ML-1)   - (A5 + C4(ML,1))*hh(ML)   - (A6 + C4(ML,2))*hh(ML+1) ...
                + A1*P_g(ML-1) + A2*P_g(ML) + A3*P_g(ML+1) ...
                + A4*h(ML-1)   + A5*h(ML)   + A6*h(ML+1);
        end

        % --- 顶部节点 / Top node ---
        RHS(NN) = -C7(NN) - C9(NN) + (C2(NN-1,2)*T(NN-1) + C2(NN,1)*T(NN)) / Delt_t ...
            - (C3(NN-1,2)/Delt_t + C6_a(NN-1))*P_gg(NN-1) - (C3(NN,1)/Delt_t + C6(NN,1))*P_gg(NN) ...
            - (C1(NN-1,2)/Delt_t + C4_a(NN-1))*hh(NN-1)   - (C1(NN,1)/Delt_t + C4(NN,1))*hh(NN) ...
            + (C3(NN-1,2)/Delt_t)*P_g(NN-1) + (C3(NN,1)/Delt_t)*P_g(NN) ...
            + (C1(NN-1,2)/Delt_t)*h(NN-1)   + (C1(NN,1)/Delt_t)*h(NN);

    elseif ~cfg.lSoilairefc && cfg.lThmrlefc
        % --- 底部节点 / Bottom node ---
        RHS(1) = -C7(1) - C9(1) + (C2(1,1)*T(1) + C2(1,2)*T(2)) / Delt_t ...
            - (C1(1,1)/Delt_t + C4(1,1))*hh(1) - (C1(1,2)/Delt_t + C4(1,2))*hh(2) ...
            + (C1(1,1)/Delt_t)*h(1) + (C1(1,2)/Delt_t)*h(2);

        % --- 内部节点 / Interior nodes ---
        for ML = 2:NL
            A4 = C1(ML-1,2)/Delt_t;   A5 = C1(ML,1)/Delt_t;   A6 = C1(ML,2)/Delt_t;

            RHS(ML) = -C7(ML) - C9(ML) + (C2(ML-1,2)*T(ML-1) + C2(ML,1)*T(ML) + C2(ML,2)*T(ML+1)) / Delt_t ...
                - (A4 + C4_a(ML-1))*hh(ML-1) - (A5 + C4(ML,1))*hh(ML) - (A6 + C4(ML,2))*hh(ML+1) ...
                + A4*h(ML-1) + A5*h(ML) + A6*h(ML+1);
        end

        % --- 顶部节点 / Top node ---
        RHS(NN) = -C7(NN) - C9(NN) + (C2(NN-1,2)*T(NN-1) + C2(NN,1)*T(NN)) / Delt_t ...
            - (C1(NN-1,2)/Delt_t + C4_a(NN-1))*hh(NN-1) - (C1(NN,1)/Delt_t + C4(NN,1))*hh(NN) ...
            + (C1(NN-1,2)/Delt_t)*h(NN-1) + (C1(NN,1)/Delt_t)*h(NN);

    else
        % --- 等温模式（仅温度存储项）/ Isothermal mode (temperature storage only) ---
        RHS(1) = -C7(1) - C9(1) + (C2(1,1)*T(1) + C2(1,2)*T(2)) / Delt_t;
        for ML = 2:NL
            RHS(ML) = -C7(ML) - C9(ML) + (C2(ML-1,2)*T(ML-1) + C2(ML,1)*T(ML) + C2(ML,2)*T(ML+1)) / Delt_t;
        end
        RHS(NN) = -C7(NN) - C9(NN) + (C2(NN-1,2)*T(NN-1) + C2(NN,1)*T(NN)) / Delt_t;
    end

    % --- C5 主对角更新：加入时间存储项 / Update C5 main diagonal: add temporal storage ---
    C5 = C5 + C2 / Delt_t;

    % --- 快照保存（slot 2 = 能量方程）/ Save snapshot (slot 2 = energy equation) ---
    SAVE(1,1,2) = RHS(1);
    SAVE(1,2,2) = C5(1,1);
    SAVE(1,3,2) = C5(1,2);
    SAVE(2,1,2) = RHS(NN);
    SAVE(2,2,2) = C5_a(NN-1);
    SAVE(2,3,2) = C5(NN,1);

    EnergyMatrices.C5 = C5;
end
