function h_frez = updateHfreez(i, L_f, SoilVariables)
% updateHfreez - 计算节点 i 的冻结水头 h_frez
%   替代原 io.getModelSettings() / io.define_constants() 依赖
%
% 输入:
%   i                    - 节点索引
%   L_f                  - 冻结潜热 (J/kg)
%   SoilVariables        - 土壤变量结构体 (含 T, h, Phi_s)
%   cfg                  - legacy ModelSettings 对象

    phys = equations.physicalConstants();

    T      = SoilVariables.T;
    h_frez = SoilVariables.h_frez;
    h      = SoilVariables.h;

    if T(i) <= 0
        h_frez(i) = L_f * 1e4 * T(i) / phys.g / phys.T_freeze;
    else
        h_frez(i) = 0;
    end

    % if cfg.lSWCC == 1
        % Van Genuchten 模式
        if h_frez(i) <= h(i) + 1e-6
            h_frez(i) = h(i) + 1e-6;
        end
    % else
    %     % Clapp-Hornberger 模式，使用节点对应的气泡压 Phi_s
    %     Phi_s_i = SoilVariables.Phi_s(min(i, length(SoilVariables.Phi_s)));
    %     if h_frez(i) <= h(i) - Phi_s_i
    %         h_frez(i) = h(i) - Phi_s_i;
    %     end
    % end

end
