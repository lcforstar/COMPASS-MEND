function nodeVals = interpolateElementToNode(elementVals, DeltZ)
% interpolateElementToNode - 单元心量加权插值到节点
%                            Weighted interpolation from element centres to shared nodes
%
% 函数定义 / Function Definition:
%   nodeVals = interpolateElementToNode(elementVals, DeltZ)
%
% 输入 / Inputs:
%   elementVals - 单元心量 / Element-centred values (NL×1)
%   DeltZ       - 各层厚度 / Layer thicknesses (NL×1) [cm]
%
% 输出 / Outputs:
%   nodeVals    - 节点量 / Node values (NN×1)
%               边界节点直接继承相邻单元值；
%               内部节点：反距离加权（权重 = 对侧单元厚度）
%               Boundary nodes: copy adjacent element; internal nodes: inverse-distance weighted average.
%
% 修改者 / Modified by:   朱寄子星 (Zhu Jizixing)
% 修改日期 / Modified date: 2026-04-25

    elementVals = elementVals(:);
    DeltZ       = DeltZ(:);
    NL          = length(elementVals);
    NN          = NL + 1;

    nodeVals      = zeros(NN, 1);
    nodeVals(1)   = elementVals(1);
    nodeVals(end) = elementVals(end);

    for i = 2:NL
        nodeVals(i) = (elementVals(i-1)*DeltZ(i) + elementVals(i)*DeltZ(i-1)) / ...
                      (DeltZ(i-1) + DeltZ(i));
    end
end
