function [QET, QEB] = calculateBoundaryFlux(SAVE, TT, cfg)
% calculateBoundaryFlux - 能量方程顶部/底部边界通量计算（返回当前时步标量）
%                         Calculate top and bottom boundary heat fluxes (scalar for current step)
%
% 函数定义 / Function Definition:
%   [QET, QEB] = calculateBoundaryFlux(SAVE, TT, cfg)
%
% 输入 / Inputs:
%   SAVE  - 系数快照 (3×3×5) / Coefficient snapshot from assembleCoefficientMatrices
%           slot 2 (energy): SAVE(:,:,2)
%             SAVE(1,:,2) = [RHS(1),   C5(1,1),    C5(1,2)   ]  底部 / Bottom
%             SAVE(2,:,2) = [RHS(NN),  C5_a(NN-1),  C5(NN,1)  ]  顶部 / Top
%   TT    - 当前时步温度 (NN×1) / Current temperature vector
%   cfg   - 模型设置 / Model settings (cfg.Grid.NL)
%
% 输出 / Outputs:
%   QET   - 顶部热通量标量（当前时步）/ Top boundary heat flux scalar (current step)
%   QEB   - 底部热通量标量（当前时步）/ Bottom boundary heat flux scalar (current step)
%
% 注：调用方负责写入时序数组：QET_series(KT) = QET
%
% 引用 / References:
%   原 Enrgy_Bndry_Flux.m；STEMMUS Technical Notes Section 5
%
% 修改者 / Modified by:   朱寄子星 (Zhu Jizixing)
% 修改日期 / Modified date: 2026-04-25
% Change: 2026-07-15 Jizixing Zhu -- document the top row's strict-lower
%         coefficient; the neighboring superdiagonal differs when VTT is nonzero.
% Change: 2026-07-15 Jizixing Zhu -- mark this T-only reaction as a legacy
%         compatibility path; production uses calculateCoupledBoundaryFlux.

    NN = cfg.Grid.NL + 1;

    % 顶部通量 / Top boundary flux
    QET =  SAVE(2,1,2) - SAVE(2,2,2)*TT(NN-1) - SAVE(2,3,2)*TT(NN);

    % 底部通量 / Bottom boundary flux
    QEB = -SAVE(1,1,2) + SAVE(1,2,2)*TT(1)    + SAVE(1,3,2)*TT(2);
end
