function Theta_UU = calculateTheta_UU(Theta_m, Gamma_hh, SV, VanGenuchten, i, cfg)
% calculateTheta_UU - 计算总含水量（第 i 层，VG 模型）
%                     Compute total volumetric water content for layer i (VG model)
%
% 函数定义 / Function Definition:
%   Theta_UU = calculateTheta_UU(Theta_m, Gamma_hh, SV, VanGenuchten, i, cfg)
%
% 输入 / Inputs:
%   Theta_m      - 毛管含水量上限 / Capillary water content upper bound
%   Gamma_hh     - 冻结分区因子 / Freezing partition factor
%   SV           - 单节点状态结构体，含 .hh / Per-node state struct with .hh
%   VanGenuchten - VG 参数结构体（各字段 NL×1）/ VG parameter struct (fields NL×1)
%   i            - 土层索引 / Layer index
%   cfg          - 模型配置对象 / Model settings (legacy ModelSettings)
%
% 输出 / Outputs:
%   Theta_UU     - 总含水量（标量）/ Total water content (scalar)
%
% 修改者 / Modified by:   朱寄子星 (Zhu Jizixing)
% 修改日期 / Modified date: 2026-04-10

    hh      = SV.hh;
    theta_s = VanGenuchten.Theta_s(i);
    theta_r = VanGenuchten.Theta_r(i);
    alpha   = VanGenuchten.Alpha(i);
    n       = VanGenuchten.n(i);
    m       = VanGenuchten.m(i);

    % 占位（Clapp-Hornberger 路径预留，当前不启用）
    lamda = 0;
    phi_s = 0;

    if cfg.lSWCC == 1
        if cfg.lFreezing == 1 && cfg.lSFCC == 1
            % lSFCC=1 路径：hh >= -1 时饱和，否则解析计算
            SoilConstants = io.getSoilConstants();
            if hh >= -1
                Theta_UU = theta_s;
            elseif hh <= SoilConstants.hd
                % VC dry 端钉死（hd=-3e4）
                Theta_UU = theta_r;
            elseif cfg.lThmrlefc
                if Gamma_hh == 0
                    Theta_UU = 0;
                else
                    Theta_UU = soil.water.conductivity.calculateTheta(0, Theta_m, hh, Gamma_hh, theta_s, theta_r, lamda, phi_s, alpha, n, m);
                end
            else
                % lThmrlefc=0：无热计算，Gamma_hh 无物理意义，退化为标准 VG
                Theta_UU = soil.water.conductivity.calculateTheta(2, Theta_m, hh, 0, theta_s, theta_r, lamda, phi_s, alpha, n, m);
            end
        else
            % lSFCC=0 路径：原版三段（Theta_UU 含冰，不做弹性储水修正）
            % lSFCC=0 path: original 3-case VG (no Ss term — Theta_UU includes ice)
            if hh >= -1e-6
                Theta_UU = theta_s;
            elseif hh <= -1e7
                Theta_UU = theta_r;
            else
                Theta_UU = soil.water.conductivity.calculateTheta(2, Theta_m, hh, Gamma_hh, theta_s, theta_r, lamda, phi_s, alpha, n, m);
            end
        end
    else
        % Clapp-Hornberger 路径（预留）
        Theta_UU = soil.water.conductivity.calculateTheta(1, Theta_m, hh, Gamma_hh, theta_s, theta_r, lamda, phi_s, alpha, n, m);
    end
end
