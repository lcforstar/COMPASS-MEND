function Theta_II = calculateTheta_II(TT, XCAP, hh, Theta_II_prev, cfg)
% calculateTheta_II - 计算冻结含水量（正弦过渡曲线）
%                     Compute ice content via sinusoidal freezing transition curve
%
% 函数定义 / Function Definition:
%   Theta_II = calculateTheta_II(TT, XCAP, hh, Theta_II_prev, cfg)
%
% 输入 / Inputs:
%   TT           - 节点温度 (°C) / Node temperature (°C)
%   XCAP         - 田间持水量（冰含量上限）/ Field capacity (ice content upper limit)
%   hh           - 截断后节点水头 (cm) / Clamped nodal matric head (cm)
%   Theta_II_prev - 上一时步冰含量 / Ice content from previous time step
%   cfg          - 模型配置对象 / Model settings (legacy ModelSettings)
%
% 输出 / Outputs:
%   Theta_II - 冻结含水量（体积比）/ Ice content (volume fraction)
%
% 说明 / Notes:
%   与原版逻辑完全一致：
%     1. 默认保留上步冰含量
%     2. hh <= -1e7 截断为 0（Clapp 路径有效，VG 路径被温度块覆盖）
%     3. lSWCC=1 && lSFCC=0 时走正弦三段简化曲线
%     4. lSFCC=1 时跳过此块，由外部严格 SFCC 计算（返回上步值）
%
% 参考 / References:
%   Lundin (1990) Cold Reg. Sci. Technol., 18, 47-64
%
% 修改者 / Modified by:   朱寄子星 (Zhu Jizixing)
% 修改日期 / Modified date: 2026-04-27

    % 无冻融模式：冰含量恒为零 / No-freezing mode: ice content is always zero
    if cfg.lFreezing ~= 1
        Theta_II = 0;
        return;
    end

    % 默认保留上步冰含量（与原版一致）/ Default: retain previous ice content
    Theta_II = Theta_II_prev;

    % hh 极干截断 / Extremely dry: zero ice (effective for Clapp path)
    if hh <= -1e7
        Theta_II = 0;
    end

    % VG + 非严格SFCC：正弦三段简化曲线 / VG + simplified sinusoidal curve
    if cfg.lSWCC == 1 && ~cfg.lSFCC
        Tf1 = 273.15 + 1;   % 上限，高于此完全液态 / Fully liquid above
        Tf2 = 273.15 - 3;   % 下限，低于此完全冻结 / Fully frozen below
        TK  = TT + 273.15;

        if TK > Tf1
            Theta_II = 0;
        elseif TK >= Tf2
            frac     = 0.5 * (1 - sin(pi() * (TK - 0.5*Tf1 - 0.5*Tf2) / (Tf1 - Tf2)));
            Theta_II = max(0, min(1, frac)) * XCAP;
        else
            Theta_II = XCAP;
        end
    end
end
