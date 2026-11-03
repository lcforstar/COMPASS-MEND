function [SW, ESW, ASW, dSWdH, PERMW, JJJH_out, IPSW1_out, MPSW1_out, ILSW_out, RSW_out, SRAW_out] = ...
    pathDecide(HAW, AlphaD, AlphaW, n, m, SARWI, ...
               RHSW_row, RASW_row, JJJH_in, IPSW1_in, MPSW1_in, RSW_in, SRAW_in)
% pathDecide - 决定走哪条扫描线 (HYSTER.f HAWPATH 简化版)
%
% 简化策略: COMPASS 起步必在 main drainage (IPSW1=1, JJJH=0)
% 暴雨场景反向点最多 2-3 个，IPATH=7 足够
%
% 决策树:
%   if HAW >= RHSW(2): 走主排水曲线 (DRAIN)
%   elif HAW < PHSW (变湿): 走湿润扫描线 (WET)
%   else (变干): 走干燥扫描线 (DRY)
%
% 输入:
%   HAW       - 当前 max(0, h_air - h_water)
%   AlphaD/W  - VG alpha (drying/wetting)
%   n, m      - VG 形状参数
%   SARWI     - 残余气体饱和度
%   RHSW_row  - 该节点反向点 |hh| 数组 (1×IPATH)
%   RASW_row  - 该节点反向点 ASW 数组
%   JJJH_in, IPSW1_in, MPSW1_in - 当前路径状态
%   RSW_in    - 历史最低 ASW
%   SRAW_in   - 当前残余气体
%
% 输出:
%   SW/ESW/ASW/dSWdH/PERMW - 该节点 K-θ-h 关系值
%   JJJH_out, IPSW1_out, MPSW1_out, ILSW_out - 更新后路径状态
%   RSW_out, SRAW_out - 更新后历史最低/残余气

    ZERO = 0.0;
    ONE  = 1.0;
    IPATH = length(RHSW_row);

    % 残余气体因子
    if SARWI > 0
        RAW = (ONE / SARWI) - ONE;
    else
        RAW = ZERO;
    end

    % 默认输出
    JJJH_out = JJJH_in;
    IPSW1_out = IPSW1_in;
    MPSW1_out = MPSW1_in;
    ILSW_out = 0;
    RSW_out = RSW_in;
    SRAW_out = SRAW_in;

    % --- 1) 主排水曲线: HAW >= RHSW(2) 或 当前在 IPSW=1 (主曲线) ---
    % 修正: IPSW=1 + JJJH=0 时永远走主曲线，避免反向点未建立时退化到 RASW(1)=1
    on_main_drainage = (IPSW1_in <= 1) || (HAW >= RHSW_row(2));
    if on_main_drainage
        [SW, ESW, ASW, dSWdH, PERMW, RSW_out, SRAW_out] = ...
            soil.water.hysteresis.drainPath(HAW, AlphaD, n, m, SARWI, RSW_in, RAW);
        JJJH_out = 0;
        IPSW1_out = 1;
        MPSW1_out = 3;
        return;
    end

    % --- 2) 简化版（IPATH 用到 2-3 即足够）---
    % 当前位置在反向点 IPSW1_in 与 IPSW1_in-1 之间
    if IPSW1_in < 2
        IPSW1_in = 2;
    end
    if IPSW1_in > IPATH
        IPSW1_in = IPATH;
    end

    % 判断方向: 与上次位置 PHSW 比较 (但简化版用 JJJH 代理)
    %   JJJH=0 (drying): HAW 增大 → 继续 dry; HAW 减小 → 反转 wet
    %   JJJH=1 (wetting): HAW 减小 → 继续 wet; HAW 增大 → 反转 dry
    %
    % 走当前活跃扫描线: 若 IPSW1 是奇数（drying 反向点），走 dryScanline
    %                     若 IPSW1 是偶数（wetting 反向点），走 wetScanline
    %
    % 简化：根据 JJJH 决定
    if JJJH_in == 0
        % 干燥扫描线: 反向点 i = IPSW1, i-1 = IPSW1-1
        i  = IPSW1_in;
        im1 = IPSW1_in - 1;
        if im1 < 1
            im1 = 1;
        end
        [SW, ESW, ASW, dSWdH, PERMW] = ...
            soil.water.hysteresis.dryScanline(HAW, AlphaD, n, m, ...
                RHSW_row(i), RHSW_row(im1), RASW_row(i), RASW_row(im1), ...
                SRAW_in, RSW_in);
        IPSW1_out = i;
        MPSW1_out = 3;
    else
        % 湿润扫描线
        i  = IPSW1_in;
        im1 = IPSW1_in - 1;
        if im1 < 1
            im1 = 1;
        end
        [SW, ESW, ASW, dSWdH, PERMW] = ...
            soil.water.hysteresis.wetScanline(HAW, AlphaW, AlphaD, n, m, ...
                RHSW_row(i), RHSW_row(im1), RASW_row(i), RASW_row(im1), ...
                SRAW_in, RSW_in);
        IPSW1_out = i;
        MPSW1_out = 3;
    end
end
