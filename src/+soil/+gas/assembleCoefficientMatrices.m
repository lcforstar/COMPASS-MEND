function [RHS, AirMatrices, SAVE] = assembleCoefficientMatrices(AirMatrices, SoilVariables, cfg, Delt_t, P_g, T, h, SAVE)
% assembleCoefficientMatrices - 干空气方程右端项组装及 C6 对角更新
%                               Assemble RHS and update diagonal C6 for dry-air PDE
%                               对应 STEMMUS Technical Notes Eq. 4.32（p.44）
%
% 函数定义 / Function Definition:
%   [RHS, AirMatrices, SAVE] = assembleCoefficientMatrices(AirMatrices, SoilVariables, cfg, Delt_t, P_g, T, h, SAVE)
%
% 输入 / Inputs:
%   AirMatrices   - 矩阵系数结构体（来自 calculateMatricCoefficients）/ Matrix coeff struct
%     .C1-.C7, .C4_a, .C5_a, .C6_a, .C6  - 各分量矩阵 / Component matrices
%   SoilVariables - 土壤状态结构体 / Soil state struct
%     .TT         - 当前 Newton 步温度 (NN×1) [K]
%     .hh         - 当前 Newton 步水头 (NN×1) [cm]
%   cfg           - 模型设置 / Model settings (cfg.Grid.NL/NN/nD, cfg.lThmrlefc)
%   Delt_t        - 时步长 [s] / Time step size
%   P_g           - 上一时步气压 (NN×1) [Pa] / Previous time-step pressure
%   T             - 上一 Newton 步温度 (NN×1) [K] / Previous Newton-step temperature
%   h             - 上一 Newton 步水头 (NN×1) [cm] / Previous Newton-step head
%   SAVE          - 矩阵系数快照（in/out）/ Coefficient snapshot (3×3×5)
%
% 输出 / Outputs:
%   RHS           - 右端向量 (NN×1) / Right-hand side vector
%   AirMatrices   - 更新后结构体（C6 已加 C3/Delt_t）/ Updated struct (C6 += C3/Delt_t)
%   SAVE          - 更新后快照（第3层写入气压方程端点系数）/ Updated snapshot
%
% 引用 / References:
%   STEMMUS Technical Notes, Section 4, Eq. 4.32；原 Air_EQ.m
%
% 注：非冻融分支内部节点用 C4_a（下对角）而非 C4(:,2)（上对角），
%     源代码此处有误，保持与 Air_EQ.m 一致。
%
% 修改者 / Modified by:   朱寄子星 (Zhu Jizixing)
% 修改日期 / Modified date: 2026-04-25
% Change: 2026-07-15 Jizixing Zhu -- snapshot the top-row strict-lower
%         coefficient; the preceding superdiagonal differs when Vaa is nonzero.

    NL       = cfg.Grid.NL;
    NN       = NL + 1;
    nD       = cfg.Grid.nD;

    % 解包矩阵系数 / Unpack matrix coefficients
    C1   = AirMatrices.C1;    C2   = AirMatrices.C2;    C3   = AirMatrices.C3;
    C4   = AirMatrices.C4;    C5   = AirMatrices.C5;    C6   = AirMatrices.C6;
    C4_a = AirMatrices.C4_a;  C5_a = AirMatrices.C5_a;
    C7   = AirMatrices.C7;

    TT = SoilVariables.TT(:);
    hh = SoilVariables.hh(:);
    T  = T(:);
    h  = h(:);

    RHS = zeros(NN, 1);

    % --- 右端项组装 / RHS assembly ---
    if cfg.lThmrlefc
        % 含热方程耦合 / With thermal coupling
        RHS(1) = -C7(1) + (C3(1,1)*P_g(1) + C3(1,2)*P_g(2)) / Delt_t ...
            - (C2(1,1)/Delt_t + C5(1,1))*TT(1) - (C2(1,2)/Delt_t + C5(1,2))*TT(2) ...
            - (C1(1,1)/Delt_t + C4(1,1))*hh(1) - (C1(1,2)/Delt_t + C4(1,2))*hh(2) ...
            + (C2(1,1)/Delt_t)*T(1) + (C2(1,2)/Delt_t)*T(2) ...
            + (C1(1,1)/Delt_t)*h(1) + (C1(1,2)/Delt_t)*h(2);

        for ML = 2:NL
            ARG1 = C2(ML-1,2)/Delt_t;  ARG2 = C2(ML,1)/Delt_t;  ARG3 = C2(ML,2)/Delt_t;
            ARG4 = C1(ML-1,2)/Delt_t;  ARG5 = C1(ML,1)/Delt_t;  ARG6 = C1(ML,2)/Delt_t;
            RHS(ML) = -C7(ML) + (C3(ML-1,2)*P_g(ML-1) + C3(ML,1)*P_g(ML) + C3(ML,2)*P_g(ML+1)) / Delt_t ...
                - (ARG1 + C5_a(ML-1))*TT(ML-1) - (ARG2 + C5(ML,1))*TT(ML) - (ARG3 + C5(ML,2))*TT(ML+1) ...
                - (ARG4 + C4_a(ML-1))*hh(ML-1) - (ARG5 + C4(ML,1))*hh(ML)  - (ARG6 + C4(ML,2))*hh(ML+1) ...
                + ARG1*T(ML-1) + ARG2*T(ML) + ARG3*T(ML+1) ...
                + ARG4*h(ML-1) + ARG5*h(ML) + ARG6*h(ML+1);
        end

        RHS(NN) = -C7(NN) + (C3(NN-1,2)*P_g(NN-1) + C3(NN,1)*P_g(NN)) / Delt_t ...
            - (C2(NN-1,2)/Delt_t + C5_a(NN-1))*TT(NN-1) - (C2(NN,1)/Delt_t + C5(NN,1))*TT(NN) ...
            - (C1(NN-1,2)/Delt_t + C4_a(NN-1))*hh(NN-1) - (C1(NN,1)/Delt_t + C4(NN,1))*hh(NN) ...
            + (C2(NN-1,2)/Delt_t)*T(NN-1) + (C2(NN,1)/Delt_t)*T(NN) ...
            + (C1(NN-1,2)/Delt_t)*h(NN-1) + (C1(NN,1)/Delt_t)*h(NN);
    else
        % 无热耦合 / Without thermal coupling
        RHS(1) = -C7(1) + (C3(1,1)*P_g(1) + C3(1,2)*P_g(2)) / Delt_t ...
            - (C1(1,1)/Delt_t + C4(1,1))*hh(1) - (C1(1,2)/Delt_t + C4(1,2))*hh(2) ...
            + (C1(1,1)/Delt_t)*h(1) + (C1(1,2)/Delt_t)*h(2);

        for ML = 2:NL
            ARG4 = C1(ML-1,2)/Delt_t;  ARG5 = C1(ML,1)/Delt_t;  ARG6 = C1(ML,2)/Delt_t;
            RHS(ML) = -C7(ML) + (C3(ML-1,2)*P_g(ML-1) + C3(ML,1)*P_g(ML) + C3(ML,2)*P_g(ML+1)) / Delt_t ...
                - (ARG4 + C4_a(ML-1))*hh(ML-1) - (ARG5 + C4(ML,1))*hh(ML) - (ARG6 + C4(ML,2))*hh(ML+1) ...
                + ARG4*h(ML-1) + ARG5*h(ML) + ARG6*h(ML+1);
        end

        RHS(NN) = -C7(NN) + (C3(NN-1,2)*P_g(NN-1) + C3(NN,1)*P_g(NN)) / Delt_t ...
            - (C1(NN-1,2)/Delt_t + C4_a(NN-1))*hh(NN-1) - (C1(NN,1)/Delt_t + C4(NN,1))*hh(NN) ...
            + (C1(NN-1,2)/Delt_t)*h(NN-1) + (C1(NN,1)/Delt_t)*h(NN);
    end

    % --- C6 对角更新：加时间导数项及数值稳定项 / Update C6 diagonal: add C3/Delt_t + stabilizer ---
    for MN = 1:NN
        for ND = 1:nD
            C6(MN, ND) = C3(MN, ND) / Delt_t + C6(MN, ND) + 1e-15;
        end
    end
    AirMatrices.C6 = C6;

    % --- SAVE 快照：记录端点 RHS 和 C6 系数 / Snapshot boundary RHS and C6 for flux calc ---
    SAVE(1,1,3) = RHS(1);     SAVE(1,2,3) = C6(1,1);     SAVE(1,3,3) = C6(1,2);
    SAVE(2,1,3) = RHS(NN);    SAVE(2,2,3) = C6_a(NN-1);  SAVE(2,3,3) = C6(NN,1);
end
