function SoilVariables = updateHyst(SoilVariables, VanGenuchten, cfg, dTheta_II)
% updateHyst - 收敛后更新 hysteresis 反向点 (HYSTER.f UpdateHyst 简化)
%
% 检测每个节点 hh 的方向变化:
%   PHSW (上一步 |hh|) vs 当前 |hh|
%   若方向反转，把当前 (|hh|, ASW) 推入下一个反向点
%
% 冻融钩子: 若该节点 |dTheta_II| > eps，跳过更新（相变扰动非真实方向反转）
%
% 输入:
%   dTheta_II - NL×nD 该步 Theta_II 变化（用于冻融钩子）

    if cfg.lHystrs ~= 1
        return;
    end

    NL = cfg.Grid.NL;
    IPATH = size(SoilVariables.RHSW, 2);

    % 当前 hh (取层 1 列；UpdateSoilWaterContent 已写入 hh)
    hh_now = SoilVariables.hh(1:NL);
    Theta_LL_now = SoilVariables.Theta_LL(:, 1);

    AlphaD = VanGenuchten.Alpha;
    Theta_s = VanGenuchten.Theta_s;
    Theta_r = VanGenuchten.Theta_r;
    n_vg = VanGenuchten.n;
    m_vg = VanGenuchten.m;

    SKIP_DT_II = 1e-4;  % 冻融跳过阈值
    DEBUG_LAYER = 116;  % MATLAB i ~= N10 反算: NL-N+1 → 取 i=10 (假设逐层从顶)

    for i = 1:NL
        % 冻融跳过
        if nargin >= 4 && ~isempty(dTheta_II) && abs(dTheta_II(i, 1)) > SKIP_DT_II
            continue;
        end

        h_now = abs(min(hh_now(i), 0));   % HAW = |hh| 当 hh<0 (饱和时 HAW=0)
        h_prev = SoilVariables.PHSW(i);

        if h_now >= 1e10 || h_now < 0
            SoilVariables.PHSW(i) = h_now;
            continue;
        end

        % 当前 ASW (从 Theta_LL 反算)
        if (Theta_s(i) - Theta_r(i)) > 1e-12
            ASW_now = (Theta_LL_now(i) - Theta_r(i)) / (Theta_s(i) - Theta_r(i));
        else
            ASW_now = 1.0;
        end
        ASW_now = max(0, min(1, ASW_now));

        % 方向判断
        JJJH = SoilVariables.JJH(i);
        IPSW1 = SoilVariables.IPSW(i);

        % 当前方向: HAW 增大 (变干, dh>0) 或减小 (变湿, dh<0)
        dh = h_now - h_prev;

        % 反转检测
        is_reversal = false;
        if JJJH == 0 && dh < -1e-3   % 之前 drying，现在变湿
            JJJH = 1;
            is_reversal = true;
        elseif JJJH == 1 && dh > 1e-3   % 之前 wetting，现在变干
            JJJH = 0;
            is_reversal = true;
        end

        if is_reversal
            % drying→wetting 反转: 新反向点是当前湿润起点，存在 IPSW1+1 的位置
            % wetting→drying 反转: 新反向点是当前干燥起点
            % 关键: drying 时 RHSW(2) 跟踪的是历史最干点，不能在反转时被覆盖
            new_idx = IPSW1 + 1;
            if new_idx > IPATH
                new_idx = IPATH;
            end
            SoilVariables.RHSW(i, new_idx) = h_now;
            SoilVariables.RASW(i, new_idx) = ASW_now;
            IPSW1 = new_idx;
        else
            % 无反转: 沿当前路径继续推进
            % drying (JJJH=0): 历史最干点存于 RHSW(2)，每次变干都更新
            % wetting (JJJH=1): 当前活跃反向点 IPSW1 存"湿润终点"
            if JJJH == 0
                % 主排水路径 IPSW1=1: 历史最干点放 RHSW(2)
                if h_now >= SoilVariables.RHSW(i, 2)
                    SoilVariables.RHSW(i, 2) = h_now;
                    SoilVariables.RASW(i, 2) = ASW_now;
                end
            elseif IPSW1 >= 2 && IPSW1 <= IPATH
                if JJJH == 1 && ASW_now > SoilVariables.RASW(i, IPSW1)
                    SoilVariables.RHSW(i, IPSW1) = h_now;
                    SoilVariables.RASW(i, IPSW1) = ASW_now;
                end
            end
        end

        SoilVariables.JJH(i) = JJJH;
        SoilVariables.IPSW(i) = IPSW1;
        SoilVariables.PHSW(i) = h_now;

        % 更新历史最低 ASW
        if ASW_now < SoilVariables.RSWAW(i)
            SoilVariables.RSWAW(i) = ASW_now;
        end
    end
end
