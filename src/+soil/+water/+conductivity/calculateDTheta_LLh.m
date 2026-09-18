function DTheta_LLh = calculateDTheta_LLh(DTheta_UUh, Theta_m, Theta_UU, Theta_LL, Gamma_hh, SV, VanGenuchten, i, cfg)
% calculateDTheta_LLh - 计算液态含水量对水头的导数 dTheta_LL/dh（4 区模型）
%                       Compute dTheta_LL/dh using 4-zone variably saturated model
%
% 函数定义 / Function Definition:
%   DTheta_LLh = calculateDTheta_LLh(~, Theta_m, Theta_UU, Theta_LL, Gamma_hh, SV, VanGenuchten, i, cfg)
%
% 输入 / Inputs:
%   ~            - 未使用（保留接口兼容性）/ Unused (kept for interface compatibility)
%   Theta_m      - 毛管含水量上限（lFreezing=1 路径）/ Capillary upper bound
%   Theta_UU     - 当前步总含水量 / Current total water content
%   Theta_LL     - 当前步液态含水量（由 calculateTheta_LL 计算）/ Current liquid content
%   Gamma_hh     - 液态水有效性系数 / Liquid water availability factor
%   SV           - 单节点状态结构体 / Per-node state struct
%                  .hh, .hh_prev, .Theta_L, .hh_frez, .Ss, .h_trans
%   VanGenuchten - VG 参数结构体（各字段 NL×1）/ VG parameter struct (NL×1)
%                  .Theta_s, .Theta_r, .Alpha, .n, .m
%   i            - 土层索引 / Layer index
%   cfg          - 模型配置对象 / Model settings (legacy ModelSettings)
%
% 输出 / Outputs:
%   DTheta_LLh   - dTheta_LL/dh（标量）/ Derivative (scalar)
%
% 修改者 / Modified by:   朱寄子星 (Zhu Jizixing)
% 修改日期 / Modified date: 2026-04-10

    theta_s  = VanGenuchten.Theta_s(i);
    theta_r  = VanGenuchten.Theta_r(i);
    alpha    = VanGenuchten.Alpha(i);
    n        = VanGenuchten.n(i);
    m        = VanGenuchten.m(i);

    theta_u  = SV.Theta_U;
    theta_l  = SV.Theta_L;
    hh       = SV.hh;
    h        = SV.hh_prev;
    hh_frez  = SV.hh_frez;
    h_frez   = SV.h_frez;
    phi_s    = 0;
    lamda    = 0;

    heatTerm = hh + hh_frez;

    % [PATCH 2026-05-08 v12] hysteresis 只改 K_h，不改 θ
    % 此处保持标准 VG 路径，dθ/dh 走原雅可比公式

    if cfg.lSWCC == 1
        if cfg.lFreezing == 1 && cfg.lSFCC == 1
            % ------------------------------------------------------------------
            % lSFCC=1 严格路径（原版 SFCC==1 逻辑）
            % ------------------------------------------------------------------
            SoilConstants = io.getSoilConstants();
            if hh >= -1
                % [PATCH 2026-05-08] 改回公开版 STEMMUS / 马铃薯版的 = 0
                % Bug B 当时改成 SV.Ss 是为暴雨期 27mm 残差救场，但水量仍不闭合
                % 公开版 + 马铃薯都是 = 0，配合 Theta_LL 钉死 θ_s 物理自洽
                DTheta_LLh = 0;
            elseif (hh + hh_frez) <= SoilConstants.hd
                % VC dry 端钉死（hd=-3e4）：dθ/dh = 0
                DTheta_LLh = 0;
            elseif cfg.lThmrlefc
                if Gamma_hh == 0
                    DTheta_LLh = 0;
                else
                    DTheta_LLh = soil.water.conductivity.calculateDTheta(0, hh, theta_s, theta_r, Theta_m, Gamma_hh, Theta_LL, theta_l, Theta_UU, theta_u, hh, h, hh_frez, h_frez, phi_s, lamda, alpha, n, m);
                end
            elseif abs(hh - h) < 1e-3
                DTheta_LLh = soil.water.conductivity.calculateDTheta(1, hh, theta_s, theta_r, Theta_m, Gamma_hh, Theta_LL, theta_l, Theta_UU, theta_u, hh, h, hh_frez, h_frez, phi_s, lamda, alpha, n, m);
            else
                DTheta_LLh = soil.water.conductivity.calculateDTheta(3, hh, theta_s, theta_r, Theta_m, Gamma_hh, Theta_LL, theta_l, Theta_UU, theta_u, hh, h, hh_frez, h_frez, phi_s, lamda, alpha, n, m);
            end

        elseif cfg.lFreezing == 1
            % ------------------------------------------------------------------
            % lSFCC=0 + lFreezing=1 简化冻融路径（原版 SFCC==0 逻辑）
            % ------------------------------------------------------------------
            if hh <= -1e7 || Theta_LL <= 0.06
                DTheta_LLh = 0;
            else
                DTheta_LLh = DTheta_UUh;
            end

        else
            % lFreezing=0：无冰，Theta_LL ≡ Theta_UU，导数直接复用
            DTheta_LLh = DTheta_UUh;
        end
    else
        % ------------------------------------------------------------------
        % Clapp-Hornberger 路径（照搬原版）/ Clapp-Hornberger path (original)
        % ------------------------------------------------------------------
        phi_s = SV.Phi_s;
        lamda = SV.Lamda;
        if hh <= -1e7 || hh + hh_frez <= -1e7 || hh + hh_frez >= phi_s
            DTheta_LLh = 0;
        elseif cfg.lThmrlefc
            DTheta_LLh = soil.water.conductivity.calculateDTheta(4, heatTerm, theta_s, theta_r, Theta_m, Gamma_hh, Theta_LL, theta_l, Theta_UU, theta_u, hh, h, hh_frez, h_frez, phi_s, lamda, alpha, n, m);
        elseif abs(hh - h) < 1e-3
            DTheta_LLh = soil.water.conductivity.calculateDTheta(5, heatTerm, theta_s, theta_r, Theta_m, Gamma_hh, Theta_LL, theta_l, Theta_UU, theta_u, hh, h, hh_frez, h_frez, phi_s, lamda, alpha, n, m);
        else
            DTheta_LLh = soil.water.conductivity.calculateDTheta(2, heatTerm, theta_s, theta_r, Theta_m, Gamma_hh, Theta_LL, theta_l, Theta_UU, theta_u, hh, h, hh_frez, h_frez, phi_s, lamda, alpha, n, m);
        end
    end
end
