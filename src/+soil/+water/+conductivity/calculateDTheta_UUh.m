function DTheta_UUh = calculateDTheta_UUh(Theta_UU, Theta_m, Theta_LL, Gamma_hh, SV, VanGenuchten, i, cfg)
% calculateDTheta_UUh - 计算总含水量对水头的导数 dTheta_UU/dh
%                       Compute derivative of total water content w.r.t. matric head
%
% 函数定义 / Function Definition:
%   DTheta_UUh = calculateDTheta_UUh(Theta_UU, Theta_m, Theta_LL, Gamma_hh, SV, VanGenuchten, i, cfg)
%
% 输入 / Inputs:
%   Theta_UU     - 当前步总含水量 / Current step total water content
%   Theta_m      - 毛管含水量上限 / Capillary water content upper bound
%   Theta_LL     - 当前步液态含水量 / Current step liquid water content
%   Gamma_hh     - 冻结分区因子 / Freezing partition factor
%   SV           - 单节点状态结构体（含 .hh, .hh_prev, .Theta_U）/ Per-node state struct
%   VanGenuchten - VG 参数结构体（各字段 NL×1）/ VG parameter struct (fields NL×1)
%   i            - 土层索引 / Layer index
%   cfg          - 模型配置对象 / Model settings (legacy ModelSettings)
%
% 输出 / Outputs:
%   DTheta_UUh   - dTheta_UU/dh（标量）/ Derivative (scalar)
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
    theta_ll = Theta_LL;
    theta_uu = Theta_UU;
    gamma_hh = Gamma_hh;
    hh       = SV.hh;
    h        = SV.hh_prev;
    h_frez   = SV.h_frez;
    hh_frez  = SV.hh_frez;

    heatTerm = hh + hh_frez;
    phi_s    = 0;   % 占位，Clapp 路径内覆盖 / placeholder, overridden in Clapp branch
    lamda    = 0;   % 占位，Clapp 路径内覆盖 / placeholder, overridden in Clapp branch

    % [PATCH 2026-05-08 v12] hysteresis 只改 K_h，不改 θ

    if cfg.lSWCC == 1
        if cfg.lFreezing == 1 && cfg.lSFCC == 1
            % lSFCC=1 严格路径（原版 SFCC==1 逻辑）
            SoilConstants = io.getSoilConstants();
            if hh >= -1
                if hh_frez >= 0
                    % [PATCH 2026-05-08] 改回公开版 = 0（同 calculateDTheta_LLh）
                    DTheta_UUh = 0;
                elseif cfg.lThmrlefc
                    DTheta_UUh = soil.water.conductivity.calculateDTheta(0, heatTerm, theta_s, theta_r, Theta_m, gamma_hh, theta_ll, theta_l, theta_uu, theta_u, hh, h, hh_frez, h_frez, phi_s, lamda, alpha, n, m);
                elseif abs(hh - h) < 1e-3
                    DTheta_UUh = soil.water.conductivity.calculateDTheta(1, heatTerm, theta_s, theta_r, Theta_m, gamma_hh, theta_ll, theta_l, theta_uu, theta_u, hh, h, hh_frez, h_frez, phi_s, lamda, alpha, n, m);
                else
                    DTheta_UUh = soil.water.conductivity.calculateDTheta(2, heatTerm, theta_s, theta_r, Theta_m, gamma_hh, theta_ll, theta_l, theta_uu, theta_u, hh, h, hh_frez, h_frez, phi_s, lamda, alpha, n, m);
                end
            elseif heatTerm <= SoilConstants.hd
                % VC dry 端钉死（hd=-3e4）：dθ/dh = 0
                DTheta_UUh = 0;
            elseif cfg.lThmrlefc
                DTheta_UUh = soil.water.conductivity.calculateDTheta(0, heatTerm, theta_s, theta_r, Theta_m, gamma_hh, theta_ll, theta_l, theta_uu, theta_u, hh, h, hh_frez, h_frez, phi_s, lamda, alpha, n, m);
            else
                % lThmrlefc=0：无热计算，退化为标准 VG 导数
                if abs(hh - h) < 1e-3
                    DTheta_UUh = soil.water.conductivity.calculateDTheta(1, heatTerm, theta_s, theta_r, Theta_m, gamma_hh, theta_ll, theta_l, theta_uu, theta_u, hh, h, hh_frez, h_frez, phi_s, lamda, alpha, n, m);
                else
                    DTheta_UUh = soil.water.conductivity.calculateDTheta(3, heatTerm, theta_s, theta_r, Theta_m, gamma_hh, theta_ll, theta_l, theta_uu, theta_u, hh, h, hh_frez, h_frez, phi_s, lamda, alpha, n, m);
                end
            end

        elseif cfg.lFreezing == 1
            % lSFCC=0 + lFreezing=1 简化冻融路径（原版 SFCC==0 逻辑）
            if hh >= -1e-6 || hh <= -1e7
                DTheta_UUh = 0;
            elseif cfg.lThmrlefc || abs(hh - h) < 1e-3
                DTheta_UUh = soil.water.conductivity.calculateDTheta(1, heatTerm, theta_s, theta_r, Theta_m, gamma_hh, theta_ll, theta_l, theta_uu, theta_u, hh, h, hh_frez, h_frez, phi_s, lamda, alpha, n, m);
            else
                DTheta_UUh = soil.water.conductivity.calculateDTheta(3, heatTerm, theta_s, theta_r, Theta_m, gamma_hh, theta_ll, theta_l, theta_uu, theta_u, hh, h, hh_frez, h_frez, phi_s, lamda, alpha, n, m);
            end

        elseif hh >= SV.h_trans
            % Zone 1: 增广饱和区，dTheta_UU/dh = Ss
            DTheta_UUh = SV.Ss;
        elseif hh <= -1e7
            % Zone 2: 极干区
            DTheta_UUh = 0;
        elseif abs(hh) < SV.h_trans
            % Zone 3: 过渡混合区，对混合公式解析求导
            % d/dh[(1-w)*Theta_unsat + w*Theta_sat]
            % = (1-w)*dTheta_unsat/dh + w*Ss + dw/dh*(Theta_sat - Theta_unsat)
            w           = (hh + SV.h_trans) / (2 * SV.h_trans);
            dw_dh       = 1 / (2 * SV.h_trans);
            term_u      = 1 + abs(alpha * hh)^n;
            Se_u        = term_u^(-m);
            Theta_unsat = theta_r + (theta_s - theta_r) * Se_u;
            Theta_sat   = theta_s + SV.Ss * hh;
            dTheta_unsat = soil.water.conductivity.calculateDTheta(5, hh, theta_s, theta_r, Theta_m, gamma_hh, theta_ll, theta_l, theta_uu, theta_u, hh, h, hh_frez, h_frez, phi_s, lamda, alpha, n, m);
            DTheta_UUh  = (1 - w) * dTheta_unsat + w * SV.Ss + dw_dh * (Theta_sat - Theta_unsat);
        elseif cfg.lThmrlefc || abs(hh - h) < 1e-3
            % Zone 4: 标准VG，解析导数
            DTheta_UUh = soil.water.conductivity.calculateDTheta(1, heatTerm, theta_s, theta_r, Theta_m, gamma_hh, theta_ll, theta_l, theta_uu, theta_u, hh, h, hh_frez, h_frez, phi_s, lamda, alpha, n, m);
        else
            % Zone 4: 标准VG，数值差分
            DTheta_UUh = soil.water.conductivity.calculateDTheta(3, heatTerm, theta_s, theta_r, Theta_m, gamma_hh, theta_ll, theta_l, theta_uu, theta_u, hh, h, hh_frez, h_frez, phi_s, lamda, alpha, n, m);
        end
    else
        % Clapp-Hornberger 路径（预留）/ Clapp-Hornberger path (reserved)
        phi_s = SV.Phi_s;
        lamda = SV.Lamda;
        if hh >= phi_s || hh <= -1e7
            DTheta_UUh = 0;
        elseif cfg.lThmrlefc
            DTheta_UUh = soil.water.conductivity.calculateDTheta(4, heatTerm, theta_s, theta_r, Theta_m, gamma_hh, theta_ll, theta_l, theta_uu, theta_u, hh, h, hh_frez, h_frez, phi_s, lamda, alpha, n, m);
        elseif abs(hh - h) < 1e-3
            DTheta_UUh = soil.water.conductivity.calculateDTheta(5, heatTerm, theta_s, theta_r, Theta_m, gamma_hh, theta_ll, theta_l, theta_uu, theta_u, hh, h, hh_frez, h_frez, phi_s, lamda, alpha, n, m);
        else
            DTheta_UUh = soil.water.conductivity.calculateDTheta(3, heatTerm, theta_s, theta_r, Theta_m, gamma_hh, theta_ll, theta_l, theta_uu, theta_u, hh, h, hh_frez, h_frez, phi_s, lamda, alpha, n, m);
        end
    end
end
