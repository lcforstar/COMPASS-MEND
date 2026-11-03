function [Theta_LL, Se_eff] = calculateTheta_LL(Theta_UU, Theta_II, Theta_m, Gamma_hh, SV, VanGenuchten, i, cfg)
% calculateTheta_LL - 计算液态含水量
%                     Compute liquid water content
%
% 函数定义 / Function Definition:
%   [Theta_LL, Se_eff] = calculateTheta_LL(Theta_UU, Theta_II, Theta_m, Gamma_hh, SV, VanGenuchten, i, cfg)
%
% 输入 / Inputs:
%   Theta_UU, Theta_II - 总含水量/冻结含水量 / Total / ice content
%   Theta_m, Gamma_hh  - 毛管上限/冻结分区因子 / VG upper bound / partition factor
%   SV                 - 单节点状态结构体 / Per-node state struct
%                        .hh, .TT, .hh_frez, .Ss, .h_trans
%   VanGenuchten       - VG 参数结构体（各字段 NL×1）/ VG parameter struct (NL×1)
%   i                  - 土层索引 / Layer index
%   cfg                - 模型配置对象 / Model settings (legacy ModelSettings)
%
% 输出 / Outputs:
%   Theta_LL - 液态含水量 / Liquid water content
%   Se_eff   - 内部有效饱和度（用于 Mualem 导水率公式）/ Internal Se for Mualem K
%
% 路径说明 / Path logic:
%   lSFCC=1                   — 严格 SFCC：hh+hh_frez 代入 VG
%   lSFCC=0 + lFreezing=1     — 简化冻融：温度分区，原版 SFCC==0 逻辑
%   lSFCC=0 + lFreezing=0     — 变饱和流 4 区模型（ClimaLand 增广）
%     Zone 1: hh >= h_trans   — 增广饱和区（积水/承压）
%     Zone 2: hh <= -1e7      — 极干区
%     Zone 3: |hh| < h_trans  — 过渡混合区
%     Zone 4: hh < -h_trans   — 标准 VG 非饱和区
%
% 参考 / References:
%   van Genuchten (1980) Soil Sci. Soc. Am. J., 44, 892-898
%   Harlan (1973) Water Resour. Res., 9, 1314-1323
%   ClimaLand augmented liquid fraction mechanism
%
% 修改者 / Modified by:   朱寄子星 (Zhu Jizixing)
% 修改日期 / Modified date: 2026-04-10

    phys    = equations.physicalConstants();
    hh      = SV.hh;
    TT      = SV.TT;
    theta_s = VanGenuchten.Theta_s(i);
    theta_r = VanGenuchten.Theta_r(i);
    Alpha   = VanGenuchten.Alpha(i);
    n       = VanGenuchten.n(i);
    m       = VanGenuchten.m(i);
    Ss      = SV.Ss;
    h_trans = SV.h_trans;

    % 占位（Clapp-Hornberger 路径预留）
    lamda = 0;
    phi_s = 0;

    % 默认输出
    Theta_LL = Theta_UU;
    Se_eff   = 1;

    % [PATCH 2026-05-08 v12 架构调整]
    % hysteresis 只作用于 KL_h（在 calculateHydraulicConductivity 内部），
    % 不修改 θ_LL — 水方程状态变量保持标准 VG 路径，永远收敛。

    if cfg.lSWCC == 1
        % lSFCC is only active when freeze-thaw physics is enabled.
        if cfg.lFreezing == 1 && cfg.lSFCC == 1
            % ------------------------------------------------------------------
            % lSFCC=1 严格路径：hh+hh_frez 代入 VG（原版 SFCC==1 逻辑）
            % ------------------------------------------------------------------
            SoilConstants = io.getSoilConstants();
            if hh >= -1 && SV.hh_frez >= 0
                Theta_LL = theta_s;
            elseif (hh + SV.hh_frez) <= SoilConstants.hd
                % VC dry 端钉死（hd=-3e4）：避开 VG dry tail 雅可比病态
                Theta_LL = theta_r;
            elseif cfg.lThmrlefc
                heatTerm = hh + SV.hh_frez;
                Theta_LL = soil.water.conductivity.calculateTheta(0, Theta_m, heatTerm, Gamma_hh, ...
                    theta_s, theta_r, lamda, phi_s, Alpha, n, m);
            else
                % lThmrlefc=0：无热方程，hh_frez 无意义，无冻结驱动力，全部视为液态
                Theta_LL = Theta_UU;
            end

        elseif cfg.lFreezing == 1
            % ------------------------------------------------------------------
            % lSFCC=0 + lFreezing=1 简化冻融路径（原版 SFCC==0 逻辑）
            % ------------------------------------------------------------------
            Tf1      = 273.15 + 1;
            heatTerm = hh + SV.hh_frez;
            if hh >= -1e-6 && TT + 273.15 > Tf1
                Theta_LL = theta_s;
            elseif hh <= -1e7
                Theta_LL = theta_r;
            elseif TT + 273.15 > Tf1
                Theta_LL = soil.water.conductivity.calculateTheta(2, Theta_m, heatTerm, Gamma_hh, ...
                    theta_s, theta_r, lamda, phi_s, Alpha, n, m);
            else
                Theta_LL = Theta_UU - Theta_II * phys.RHOI / phys.RHOL;
            end
            if Theta_LL <= 0.06
                Theta_LL = 0.06;
            end

        else
            % ------------------------------------------------------------------
            % lSFCC=0 + lFreezing=0 变饱和流 4 区路径
            % ------------------------------------------------------------------

            % Zone 1: 增广饱和区 (hh >= h_trans)
            if hh >= h_trans
                Theta_LL = theta_s + Ss * hh;
                Se_eff   = 1;

            % Zone 2: 极干区 (hh <= -1e7)
            elseif hh <= -1e7
                Theta_LL = theta_r;
                Se_eff   = 0;

            % Zone 3: 过渡混合区 (|hh| < h_trans)
            elseif abs(hh) < h_trans
                w           = (hh + h_trans) / (2 * h_trans);
                term_u      = 1 + abs(Alpha * hh)^n;
                Se_unsat    = term_u^(-m);
                Theta_unsat = theta_r + (theta_s - theta_r) * Se_unsat;
                Theta_sat   = theta_s + Ss * hh;
                Theta_LL    = (1 - w) * Theta_unsat + w * Theta_sat;
                Se_eff      = (1 - w) * Se_unsat + w * 1;

            % Zone 4: 标准 VG 非饱和区 (hh < -h_trans)
            else
                Theta_LL = soil.water.conductivity.calculateTheta(2, Theta_m, hh, Gamma_hh, theta_s, theta_r, lamda, phi_s, Alpha, n, m);
                Se_eff   = (Theta_LL - theta_r) / (theta_s - theta_r);
            end
        end
    else
        % ------------------------------------------------------------------
        % Clapp-Hornberger 路径（预留）/ Clapp-Hornberger path (reserved)
        % ------------------------------------------------------------------
        if hh >= phi_s
            Theta_LL = soil.water.conductivity.calculateTheta(1, Theta_m, hh, Gamma_hh, ...
                theta_s, theta_r, lamda, phi_s, Alpha, n, m);
        elseif hh <= -1e7
            Theta_LL = theta_r;
        end
    end
end
