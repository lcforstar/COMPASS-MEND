function [RAI] = calculateRootAreaIndex1(root_type, Lmax, RL, NN, DeltZ)
%CALCULATEROOTAREAINDEX1  根面积指数（RAI）剖面计算
%
%   与 init.initialRootBiomass 同出一源（Jackson 1996），但计算的物理量不同：
%
%     initialRootBiomass   → 根长密度 Rl [cm/cm³]，从总生物量(RTB)出发，
%                            使用地表覆盖字符串查 β，DeltZ 方向：top-down（层1=地表）
%     calculateRootAreaIndex1 → 根面积指数 RAI [m²/m²]，直接用 biome 最大 RAI 表，
%                            使用 biome 编号(root_type 1-11)查 β，DeltZ 方向：bottom-up（层1=最深）
%
%   两个函数分别供 Plantrootconductivity（土-根导水，用RAI）和
%   ebal（根系吸水，用Rl）调用，不可互换。
%
%   Jackson 1996 β 累积分布：
%     Y(d) = 1 - β^d  ←  d 为距地表深度（cm），Y(d) 为地表到 d 之间的根系累积比例
%     层 [d_top, d_bot] 内的 RAI = RAImax × (β^d_top - β^d_bot)
%
%   参考文献：
%     Jackson et al. (1996) A global analysis of root distributions
%       for terrestrial biomes. Oecologia 108: 389-411.
%
%   输入：
%     root_type : 整数 1-12（Jackson 1996 biome 类型 + 长武苹果园实测）
%                 1=寒温带林  2=荒漠  3=硬叶灌丛  4=温带针叶  5=温带落叶
%                 6=温带草地  7=热带落叶  8=热带常绿  9=热带草原  10=冻原
%                 11=长武苹果园<10a（Li et al.2023, β=0.995, D50=80cm）
%                 12=长武苹果园>10a（Li et al.2023, β=0.998, D50=512cm）← 推荐用于11龄果树
%     Lmax      : 根系最大深度 [cm]
%     RL        : 土壤剖面总深度 [cm]
%     NN        : 节点数（层数 NL = NN-1）
%     DeltZ     : 各层厚度 [cm]，底部向上排列（DeltZ(1)=最深层，DeltZ(NL)=地表层）
%
%   输出：
%     RAI       : (NL×1)，各层根面积指数 [m²/m²]，底部向上排列（与 DeltZ 相同）
%
%   Author: Zhu Jizixing
%   Date:   2025-09-09  (bugfix 2026-04-16, root_type 11-12 added 2026-05-11)

    % biome 参数表（type 1-10: Jackson 1996；type 11-12: 长武苹果园实测 Li et al.2023）
    % RAImax: type11/12 沿用温带落叶林 11.0，β 来自实测
    RAIlist  = [4.6, 10, 5.5, 11.6, 11.0, 9.8, 79.1, 6.3, 7.4, 42.5, 11.0, 11.0];   % 最大 RAI [m²/m²]
    betalist = [0.943, 0.961, 0.970, 0.950, 0.980, 0.967, 0.943, 0.982, 0.972, 0.972, 0.995, 0.998];

    NL      = NN - 1;
    RAI     = zeros(NL, 1);
    RAImax  = RAIlist(root_type);
    beta    = betalist(root_type);

    % 计算各层顶底距地表深度（DeltZ 为 bottom-up，累积方向=从底部向地表）
    % Elmn_Lnth(ML) = 从底部累积到第 ML 层顶端的总厚度
    %   → 第 ML 层顶端距地表 = RL - Elmn_Lnth(ML)
    %   → 第 ML 层底端距地表 = RL - Elmn_Lnth(ML) + DeltZ(ML)  （= RL - Elmn_Lnth(ML-1)）
    Elmn_Lnth = 0;
    for ML = 1:NL
        Elmn_Lnth = Elmn_Lnth + DeltZ(ML);

        d_top = RL - Elmn_Lnth;                   % 本层顶端距地表 [cm]
        d_bot = d_top + DeltZ(ML);                % 本层底端距地表 [cm]

        % 仅统计根区（0 ~ Lmax）内的部分
        d1 = max(0,    d_top);
        d2 = min(Lmax, d_bot);

        if d2 <= d1
            % 本层完全在根区之外（深于 Lmax）
            RAI(ML) = 0;
        else
            % RAI = RAImax × [Y(d2) - Y(d1)] = RAImax × (β^d1 - β^d2)
            RAI(ML) = RAImax * (beta^d1 - beta^d2);
        end
    end

    % 归一化：补偿 Lmax 截断导致的 β^Lmax 损失，使 sum(RAI) = RAImax
    RAI = RAI / (1 - beta^Lmax);
end
