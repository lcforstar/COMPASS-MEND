function Se_eff = calculateSe(Theta_LL, POR, SV, cfg)
% calculateSe - 计算有效饱和度（截断至 [0,1]，用于 Mualem 导水率公式）
%               Compute effective saturation clamped to [0,1] for Mualem K formula
%
% 函数定义 / Function Definition:
%   Se_eff = calculateSe(Theta_LL, POR, SV, cfg)
%
% 输入 / Inputs:
%   Theta_LL - 液态含水量 / Liquid water content
%   POR      - 孔隙度 / Porosity
%   SV       - 单节点状态结构体 / Per-node state struct
%              .hh; lSFCC=1 另需 .hh_frez
%   cfg      - 模型配置对象 / Model settings (legacy ModelSettings)
%
% 输出 / Outputs:
%   Se_eff - 有效饱和度 ∈ [0,1] / Effective saturation in [0,1]
%
% 路径说明 / Path logic:
%   lSFCC=1                   — 严格 SFCC：hh+hh_frez 判断（原版 SFCC==1）
%   lSFCC=0 + lFreezing=1     — 简化冻融：原版 SFCC==0 逻辑
%   lSFCC=0 + lFreezing=0     — 无冻融变饱和流
%
% 修改者 / Modified by:   朱寄子星 (Zhu Jizixing)
% 修改日期 / Modified date: 2026-04-27

    hh           = SV.hh;
    seCalculated = Theta_LL / POR;
    hd           = io.getSoilConstants().hd;

    if cfg.lSWCC == 1
        if cfg.lFreezing == 1 && cfg.lSFCC == 1
            % ------------------------------------------------------------------
            % lSFCC=1 严格路径（原版 SFCC==1 逻辑）
            % ------------------------------------------------------------------
            hh_frez = SV.hh_frez;
            if hh >= -1
                if hh_frez >= 0
                    Se_eff = 1;
                elseif cfg.lThmrlefc
                    if (hh + hh_frez) <= hd
                        Se_eff = 0;
                    else
                        Se_eff = seCalculated;
                    end
                else
                    % lThmrlefc=0：无热方程，退化为液态比例
                    Se_eff = seCalculated;
                end
            elseif cfg.lThmrlefc
                if (hh + hh_frez) <= hd
                    Se_eff = 0;
                else
                    Se_eff = seCalculated;
                end
            else
                % lThmrlefc=0：无热方程，退化为液态比例
                Se_eff = seCalculated;
            end

        elseif cfg.lFreezing == 1
            % ------------------------------------------------------------------
            % lSFCC=0 + lFreezing=1 简化冻融路径（原版 SFCC==0 逻辑）
            % ------------------------------------------------------------------
            if Theta_LL <= 0.06 || hh <= -1e7
                Se_eff = 0;
            else
                Se_eff = seCalculated;
            end

        else
            % ------------------------------------------------------------------
            % lSFCC=0 + lFreezing=0 无冻融变饱和流路径
            % ------------------------------------------------------------------
            if Theta_LL <= 0.06 || hh <= -1e7
                Se_eff = 0;
            else
                Se_eff = seCalculated;
            end
        end

    else
        % ------------------------------------------------------------------
        % Clapp-Hornberger 路径（预留）/ Clapp-Hornberger path (reserved)
        % ------------------------------------------------------------------
        hh_frez = SV.hh_frez;
        if hh + hh_frez <= -1e7 || hh <= -1e7
            Se_eff = 0;
        elseif hh + hh_frez >= SV.Phi_s
            Se_eff = 1;
        else
            Se_eff = seCalculated;
        end
    end

    % 截断至有效范围 / Clamp to [0,1]
    if Se_eff > 1;  Se_eff = 1;  end
    if Se_eff < 0;  Se_eff = 0;  end
end
