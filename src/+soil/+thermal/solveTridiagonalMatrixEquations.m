function [SoilVariables, CHK, RHS, EnergyMatrices] = solveTridiagonalMatrixEquations( ...
    EnergyMatrices, SoilVariables, RHS, cfg)
% solveTridiagonalMatrixEquations - Thomas 算法求解能量三对角方程组
%                                   Solve energy tridiagonal system via Thomas algorithm
%                                   对应 STEMMUS Technical Notes Eq. 4.25（p.41）
%
% 函数定义 / Function Definition:
%   [SoilVariables, CHK, RHS, EnergyMatrices] = solveTridiagonalMatrixEquations( ...
%       EnergyMatrices, SoilVariables, RHS, cfg)
%
% 输入 / Inputs:
%   EnergyMatrices - 矩阵系数结构体 / Matrix coefficient struct
%     .C5          - 主/上对角 (NN×2)
%     .C5_a        - 下对角 (NL×1)
%   SoilVariables  - 土壤状态结构体 / Soil state struct
%     .TT          - 当前 Newton 步温度 (NN×1) / Current Newton-step temperature
%   RHS            - 右端向量 (NN×1，施加 BC 后）/ RHS after BC application
%   cfg            - 模型设置 / Model settings (cfg.Grid.NL, cfg.Grid.NN)
%
% 输出 / Outputs:
%   SoilVariables  - 更新后结构体（TT 已更新）/ Updated struct (TT updated)
%   CHK            - Newton 收敛残差 (NN×1) / Newton convergence residuals
%   RHS            - 消元后向量（即解向量）/ Eliminated vector (= solution)
%   EnergyMatrices - 更新后结构体（C5 已被前向消元修改）/ Updated struct (C5 modified)
%
% 引用 / References:
%   STEMMUS Technical Notes, Section 5；原 Enrgy_Solve.m
%
% 修改者 / Modified by:   朱寄子星 (Zhu Jizixing)
% 修改日期 / Modified date: 2026-04-25

    NL = cfg.Grid.NL;
    NN = NL + 1;

    % --- 前向消元 / Forward elimination ---
    RHS(1) = RHS(1) / EnergyMatrices.C5(1, 1);

    for i = 2:NN
        EnergyMatrices.C5(i, 1) = EnergyMatrices.C5(i, 1) ...
            - EnergyMatrices.C5_a(i-1) * EnergyMatrices.C5(i-1, 2) / EnergyMatrices.C5(i-1, 1);
        RHS(i) = (RHS(i) - EnergyMatrices.C5_a(i-1) * RHS(i-1)) / EnergyMatrices.C5(i, 1);
    end

    % --- 后向回代 / Back substitution ---
    for i = NL:-1:1
        RHS(i) = RHS(i) - EnergyMatrices.C5(i, 2) * RHS(i+1) / EnergyMatrices.C5(i, 1);
    end

    % --- 收敛残差与温度更新 / Convergence residuals and temperature update ---
    CHK = abs(RHS(:) - SoilVariables.TT(:));
    SoilVariables.TT = RHS(:);
end
