function [SW, ESW, ASW, dSWdH, PERMW] = wetPath(HAW, AlphaW, AlphaD, n, m, RSW)
% wetPath - 主湿润曲线 (path 2, HAW=0 收敛端)
%
% K-P 关系: AlphaW = 2 * AlphaD, n_W = n_D
%
% 输入:
%   HAW    - 当前空气-水水头差 [cm]
%   AlphaW - 湿润 alpha = 2 * AlphaD
%   AlphaD - 排水 alpha
%   n, m   - VG 形状参数
%   RSW    - 历史最低 ASW
%
% 输出: 与 drainPath 一致

    ZERO = 0.0;
    ONE = 1.0;
    XXM = ONE / m;

    SWI = (ONE + (AlphaW * HAW)^n)^(-m);
    SWI = min(SWI, ONE);
    ASW = SWI;

    ESAT = ZERO;   % 简化
    ESW  = ASW - ESAT;
    ESW  = max(ESW, ZERO);
    SW   = ESW;

    % Mualem K
    P0 = ONE - (ONE - ESW^XXM)^m;
    PERMW = sqrt(max(ESW,0)) * (P0^2);

    % dθ/dh
    if ASW > ZERO && ASW < ONE
        X1 = ASW^XXM;
        dSWdH = AlphaW * (n - ONE) * X1 * (ONE - X1)^m;
    else
        dSWdH = ZERO;
    end
end
