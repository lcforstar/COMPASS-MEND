function [AirVariables, P_gg, SAVE] = solveDryAirEquations(SoilVariables, TransportCoefficient, cfg, BoundaryCondition, P_g, KT, Delt_t, SAVE, T, h)
% solveDryAirEquations - 求解干空气方程，更新土壤气压 P_gg
%                        Solve dry-air PDE (Thomas algorithm) to update P_gg
%
% 函数定义 / Function Definition:
%   [AirVariables, P_gg, SAVE] = solveDryAirEquations(SoilVariables, TransportCoefficient, ...
%       cfg, BoundaryCondition, P_g, KT, Delt_t, SAVE, T, h)
%
% 输入 / Inputs:
%   SoilVariables       - 土壤状态结构体 / Soil state struct (hh, TT, KfL_h, ...)
%   TransportCoefficient - 传输系数结构体 / Transport coefficient struct (D_V, D_Vg, ...)
%   cfg                 - 模型设置 / Model settings (Grid, lThmrlefc, ...)
%   BoundaryCondition   - 边界条件结构体 / Boundary condition struct
%   P_g                 - 上一时步气压 (NN×1) [Pa] / Previous time-step pressure
%   KT                  - 外层时步索引 / Outer time step index
%   Delt_t              - 时步长 [s] / Time step size
%   SAVE                - 矩阵系数快照 (3×3×5) / Coefficient snapshot for flux calc
%   T                   - 上一 Newton 步温度 (NN×1) [K] / Previous Newton-step temperature
%   h                   - 上一 Newton 步水头 (NN×1) [cm] / Previous Newton-step head
%
% 输出 / Outputs:
%   AirVariables        - 干空气系数结构体 / Dry-air coefficient struct (from calculateDryAirParameters)
%   P_gg                - 求解后气压 (NN×1) [Pa] / Updated air pressure
%   SAVE                - 更新后快照 / Updated coefficient snapshot
%
% 引用 / References:
%   STEMMUS Technical Notes, Section 4, Eq. 4.32
%
% 修改者 / Modified by:   朱寄子星 (Zhu Jizixing)
% 修改日期 / Modified date: 2026-04-25

    DSTMAX = BoundaryCondition.DSTMAX;

    % --- [1] 干空气方程系数 / Dry-air equation coefficients ---
    AirVariables = soil.gas.calculateDryAirParameters(SoilVariables, TransportCoefficient, cfg, DSTMAX);

    % --- [2] 组装三对角矩阵 / Assemble tridiagonal matrix ---
    AirMatrices = soil.gas.assembleDryAirEquationMatrices(AirVariables, cfg);

    % --- [3] 组装右端项 / Assemble RHS ---
    [RHS, AirMatrices, SAVE] = soil.gas.assembleCoefficientMatrices( ...
        AirMatrices, SoilVariables, cfg, Delt_t, P_g, T, h, SAVE);

    % --- [4] 施加边界条件 / Apply boundary conditions ---
    [RHS, AirMatrices] = soil.gas.calculateBoundaryConditions(BoundaryCondition, AirMatrices, RHS, KT);

    % --- [5] Thomas 算法求解 / Solve with Thomas algorithm ---
    [~, P_candidate, ~] = soil.gas.solveTridiagonalMatrixEquations(RHS, AirMatrices, cfg);

    P_prev = SoilVariables.P_gg(:);
    bad_pressure = isnan(P_candidate) | ~isreal(P_candidate);
    P_candidate(bad_pressure) = P_prev(bad_pressure);
    P_gg = P_candidate(:);
end
