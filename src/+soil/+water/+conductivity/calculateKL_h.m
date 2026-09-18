function KL_h = calculateKL_h(MU_W, MU_WN, Se, Ks, m_vg)
% calculateKL_h - 计算液态导水率（Mualem 模型 + 粘滞温度修正）
%                 Compute liquid hydraulic conductivity (Mualem model + viscosity correction)
%
% 函数定义 / Function Definition:
%   KL_h = calculateKL_h(MU_W, MU_WN, Se, Ks, m_vg)
%
% 输入 / Inputs:
%   MU_W  - 当前温度下液水动力粘度 g/(cm·s) / Dynamic viscosity at current T
%   MU_WN - 参考温度（20°C）液水动力粘度 / Reference viscosity at 20°C
%   Se    - 有效饱和度 ∈ [0, 1] / Effective saturation
%   Ks    - 饱和导水率 (cm/s) / Saturated hydraulic conductivity
%   m_vg  - van Genuchten m 参数 / VG m parameter
%
% 输出 / Outputs:
%   KL_h  - 液态导水率 (cm/s) / Liquid hydraulic conductivity
%
% 参考 / References:
%   Mualem (1976) Water Resour. Res., 12, 513-522
%   van Genuchten (1980) Soil Sci. Soc. Am. J., 44, 892-898
%
% 修改者 / Modified by:   朱寄子星 (Zhu Jizixing)
% 修改日期 / Modified date: 2026-04-10

    % 粘滞温度修正系数 / Viscosity temperature correction factor
    CKT = MU_WN / MU_W;

    if Se == 0
        KL_h = 0;
        return;
    end

    % Mualem-van Genuchten 导水率公式 / Mualem-van Genuchten conductivity formula
    KL_h = CKT * Ks * (Se^0.5) * (1 - (1 - Se^(1 / m_vg))^m_vg)^2;

    % 数值下限 / Numerical lower bound
    if KL_h <= 1e-20
        KL_h = 1e-20;
    end

    % 异常值处理 / Handle NaN or complex
    if isnan(KL_h)
        KL_h = 0;
        warning('calculateKL_h: isnan(KL_h), set to 0');
    end
    if ~isreal(KL_h)
        warning('calculateKL_h: ~isreal(KL_h)');
    end
end
