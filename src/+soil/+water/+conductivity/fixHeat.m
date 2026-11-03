function hh = fixHeat(hh, cfg)
% fixHeat - 将节点水头截断至数值安全范围
%           Clamp nodal matric head to numerically safe range
%
% 函数定义 / Function Definition:
%   hh = fixHeat(hh, cfg)
%
% 输入 / Inputs:
%   hh  - COR 温度修正后的节点水头 (cm, scalar)
%         COR-corrected nodal matric head
%   cfg - 模型配置对象 / Model settings (legacy ModelSettings)
%         相关开关: lSWCC, lSFCC
%
% 输出 / Outputs:
%   hh  - 截断后水头 (cm) / Clamped head
%
% 说明 / Notes:
%   lSWCC==1 && lSFCC==0（VG+简化路径）：下限截断 -1e7
%     原版上限为 -1e-6，因 4 区模型支持正 hh（Zone 1 积水区），此处去除上限截断
%   lSWCC==1 && lSFCC==1（VG+严格SFCC）：不截断，由 hh_frez+VG 曲线自行处理
%   lSWCC==0（Clapp，预留）：截断至 [−1e7, Phi_s]
%
% 修改者 / Modified by:   朱寄子星 (Zhu Jizixing)
% 修改日期 / Modified date: 2026-04-27

    if cfg.lSWCC == 1
        if ~(cfg.lFreezing == 1 && cfg.lSFCC == 1)
            % 简化路径：仅下限截断，上限开放以支持 4 区模型正水头
            if hh < -1e7;  hh = -1e7;  end
        end
        % lSFCC==1：不截断
    end

    % Clapp-Hornberger 路径预留（当前不启用）
    % if hh >= Phi_s,   hh = Phi_s;   end
    % if hh <= -1e7,    hh = -1e7;    end
end
