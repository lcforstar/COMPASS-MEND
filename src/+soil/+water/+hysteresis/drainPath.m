function [SW, ESW, ASW, dSWdH, PERMW, RSW, SRAW] = drainPath(HAW, AlphaD, n, m, ...
    SARWI, RSW_in, RAW)
% drainPath - 主排水曲线 ASW/SW/K (HYDRUS HYSTER.f DRAIN subroutine)
%
% 输入:
%   HAW    - 当前空气-水水头差 max(0, h_air - h_water) [cm]
%            COMPASS 中 hh_air ~ 0 时 HAW = max(0, -hh_eff) = |hh_eff| 当 hh < 0
%   AlphaD - 主排水 alpha [1/cm]
%   n      - VG n
%   m      - 1 - 1/n
%   SARWI  - 主湿润最大捕获空气饱和度
%   RSW_in - 历史最低 ASW (residual ASW 跟踪)
%   RAW    - 残余空气因子 (1/SARWI - 1)
%
% 输出:
%   SW     - 实际水饱和度 SW = ESW*(1-Sm) + Sm; 但此函数中 ESW 是相对量，SW 直接是 saturation
%            注：HYDRUS 将 SW 解释为相对于 (theta_s - theta_r) 的饱和度
%            此处简化为直接返回 ASW（apparent saturation = θ/θs 的近似）
%   ESW    - 有效水饱和度 = ASW - ESAT
%   ASW    - 表观水饱和度 (= (1+(α·HAW)^n)^(-m))
%   dSWdH  - dθ/dh (用于容水率 Cap)
%   PERMW  - 相对水导率 (Mualem)
%   RSW    - 更新后的历史最低 ASW
%   SRAW   - 当前残余气体饱和度
%
% 简化版（IPATH > 1 hysteresis 路径，HYSTER.f line 506-514）

    ZERO = 0.0;
    ONE = 1.0;

    XXM = ONE / m;

    ASW = (ONE + (AlphaD * HAW)^n)^(-m);
    ASW = min(ASW, ONE);

    % 跟踪历史最低 ASW
    RSW = RSW_in;
    if ASW < RSW
        RSW = ASW;
    end

    % 残余气体饱和度
    if SARWI > 0
        SRAW = (ONE - ASW) / (ONE + RAW * (ONE - ASW));
    else
        SRAW = ZERO;
    end

    ESAT = ZERO;   % 主排水曲线无气体被困
    ESW  = ASW - ESAT;
    SW   = ESW;

    % Mualem K (HYSTER.f line 513-514)
    P0 = ONE - (ONE - ESW^XXM)^m;
    PERMW = sqrt(ESW) * (P0^2);

    % dθ/dh (HYSTER.f line 508-509)
    if ASW > ZERO && ASW < ONE
        X1 = ASW^XXM;
        dSWdH = AlphaD * (n - ONE) * X1 * (ONE - X1)^m;
    else
        dSWdH = ZERO;
    end
end
