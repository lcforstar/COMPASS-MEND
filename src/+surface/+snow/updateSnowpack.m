function [Throughfall_liquid, SWE_out] = updateSnowpack(Throughfall_in, Tss, SWE_in, Delt_t)
% updateSnowpack  简单温度指数积雪模块
%   Simple temperature-index snowpack: accumulate when frozen, release when thawed.
%
%   输入 / Inputs:
%     Throughfall_in  - 冠层穿透雨 [cm/s]（Gash 截留后）
%     Tss             - 地表温度 [°C]
%     SWE_in          - 当前积雪水当量 [cm]
%     Delt_t          - 时步长度 [s]
%
%   输出 / Outputs:
%     Throughfall_liquid - 实际到达土面的液态水通量 [cm/s]
%     SWE_out            - 更新后的积雪水当量 [cm]

if Tss <= 0
    % 冻结：降水以雪形式存储，液态通量为零
    SWE_out            = SWE_in + Throughfall_in * Delt_t;
    Throughfall_liquid = 0;
else
    % 融化：积雪全部释放（瞬时融化），叠加当步穿透雨
    Throughfall_liquid = Throughfall_in + SWE_in / Delt_t;
    SWE_out            = 0;
end
end
