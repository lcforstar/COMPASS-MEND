function SoilVariables = hysterIni(SoilVariables, VanGenuchten, cfg)
% hysterIni - 初始化 HYDRUS-1D Lenhard (2004) 滞回模型状态
%             Initialize Lenhard hysteresis state arrays
%
% 假设初始条件在主排水曲线上 (iHyst=1)，参考 HYSTER.f line 188-199
%
% K-P 关系: ALPHAI = 2 * ALPHAD, n_W = n_D (Kool & Parker 1987)
%
% 输出 SoilVariables 新增字段:
%   RHSW(NL, IPATH)  反向点 |hh| 历史
%   RASW(NL, IPATH)  反向点 ASW (apparent saturation) 历史
%   IPSW(NL)         当前路径序号
%   JJH(NL)          0=干燥支, 1=湿润支
%   MPSW(NL)         路径模式 (0/1/3)
%   PHSW(NL)         上一步 |hh| (用于方向判断)
%   SARW(NL)         残余空气饱和度 (SARWI 初值)
%   RSWAW(NL)        实际剩余水
%   AlphaW(NL)       湿润曲线 alpha (= 2 * AlphaD)
%
% 修改者: Claude / 朱寄子星 (Zhu Jizixing)
% 修改日期: 2026-05-08

    if cfg.lHystrs ~= 1
        % 滞回未开启，留空字段防止下游报错
        return;
    end

    NL = cfg.Grid.NL;
    IPATH = 7;

    % 初始化反向点数组
    RHSW = zeros(NL, IPATH);
    RASW = zeros(NL, IPATH);

    % 第 1 个反向点恒为 (h=0, ASW=1) — 完全饱和参考点
    RASW(:, 1) = 1.0;
    RHSW(:, 1) = 0.0;

    % iHyst=1: 从主排水曲线起步
    % 偶/奇位反向点初始化（HYSTER.f line 141-167）
    for ij = 4:2:IPATH
        RASW(:, ij) = 1.0;
    end
    for ij = 3:2:(IPATH-1)
        RASW(:, ij) = 0.0;
    end
    if mod(IPATH, 2) == 1   % IPATH=7 奇数
        RASW(:, IPATH) = 0.0;
    end

    % iHyst=1 (main drainage) 起步状态
    RHSW(:, 2) = 0.0;
    RASW(:, 2) = 1.0;

    % K-P: AlphaW = 2 * AlphaD, n_W = n_D
    AlphaW = 2.0 * VanGenuchten.Alpha(:);

    % SARWI 残余空气饱和度（K-P 经验值: 主湿润饱和度 ~0.95 主排水饱和度）
    % 简化: SARWI = 0.05（5% 孔隙被气体捕获）
    SARWI = 0.05 * ones(NL, 1);

    SoilVariables.RHSW   = RHSW;
    SoilVariables.RASW   = RASW;
    SoilVariables.IPSW   = ones(NL, 1);     % 当前在反向点 1
    SoilVariables.JJH    = zeros(NL, 1);    % 0 = drying
    SoilVariables.MPSW   = 3 * ones(NL, 1); % 默认路径模式
    SoilVariables.PHSW   = zeros(NL, 1);
    SoilVariables.SARW   = SARWI;
    SoilVariables.RSWAW  = ones(NL, 1);     % 主排水起点 RSWAW=1
    SoilVariables.AlphaW = AlphaW;
    SoilVariables.SARWI  = SARWI;            % 常量副本，方便 path 调用
end
