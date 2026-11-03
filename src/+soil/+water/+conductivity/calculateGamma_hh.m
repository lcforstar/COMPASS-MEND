function gamma_hh = calculateGamma_hh(hh)
% calculateGamma_hh - 计算液态水有效性系数 Gamma_hh
%                     Compute liquid water availability factor Gamma_hh
%
% 函数定义 / Function Definition:
%   gamma_hh = calculateGamma_hh(hh)
%
% 输入 / Inputs:
%   hh       - 节点基质水头 (cm, scalar) / Nodal matric head
%
% 输出 / Outputs:
%   gamma_hh - 液态水有效性系数 ∈ [0, 1] / Liquid water availability factor
%              用于计算水汽扩散增强因子 / Used for vapor diffusion enhancement
%              0 = 极干（液态水对水汽传输无贡献）
%              1 = 足够湿润（液态水充足，无需折减）
%
% 说明 / Notes:
%   三段逻辑（利用 io.getSoilConstants 中的 hd, hm）：
%     |hh| >= |hd|：绝干区，gamma_hh = 0
%     |hm| <= |hh| < |hd|：对数过渡区
%       gamma_hh = log(|hd|/|hh|) / log(|hd|/|hm|)
%     |hh| < |hm|：湿润区，gamma_hh = 1
%   其中 hd = -1e7 cm（绝干参考水头），hm = -9899 cm（≈-10 bar，湿态上界）
%
% 修改者 / Modified by:   朱寄子星 (Zhu Jizixing)
% 修改日期 / Modified date: 2026-04-11

    SoilConstants = io.getSoilConstants();

    if abs(hh) >= abs(SoilConstants.hd)
        % 绝干区：液态水无贡献 / Oven-dry: no liquid water contribution
        gamma_hh = 0;
    elseif abs(hh) >= abs(SoilConstants.hm)
        % 过渡区：对数折减 / Transition zone: log interpolation
        gamma_hh = log(abs(SoilConstants.hd) / abs(hh)) / ...
                   log(abs(SoilConstants.hd) / abs(SoilConstants.hm));
    else
        % 湿润区：液态水充足 / Moist zone: liquid water sufficient
        gamma_hh = 1;
    end
end
