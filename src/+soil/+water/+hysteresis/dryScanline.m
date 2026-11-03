function [SW, ESW, ASW, dSWdH, PERMW] = dryScanline(HAW, AlphaD, n, m, RHSW_i, RHSW_im1, RASW_i, RASW_im1, SRAW, RSW)
% dryScanline - 干燥扫描线插值 (HYSTER.f DRY subroutine)
%
% 在两个反向点 i-1 (湿润端, 较小 ASW) 和 i (干燥端, 较大 ASW) 之间插值
%
% 公式 (line 553-554):
%   ASW = ((SWD - RDISWD) * (RASW[i] - RASW[i-1])) / (RIDSWD - RDISWD) + RASW[i-1]
% 其中:
%   RIDSWD = (1 + (αD·R[i])^n)^(-m)        ← 反向点 i 在主排水曲线上的"如果值"
%   RDISWD = (1 + (αD·R[i-1])^n)^(-m)
%   SWD    = (1 + (αD·HAW)^n)^(-m)
%
% 这是关于 SWD 的线性插值，保证 HAW=R[i] 时 ASW=RASW[i]，连续过渡

    ZERO = 0.0;
    ONE = 1.0;
    XXM = ONE / m;

    RIDSWD = (ONE + (AlphaD * RHSW_i)^n)^(-m);
    RDISWD = (ONE + (AlphaD * RHSW_im1)^n)^(-m);
    SWD    = (ONE + (AlphaD * HAW)^n)^(-m);

    denom = RIDSWD - RDISWD;
    if abs(denom) < 1e-12
        ASW = RASW_im1;
    else
        ASW = ((SWD - RDISWD) * (RASW_i - RASW_im1)) / denom + RASW_im1;
    end
    ASW = min(ASW, ONE);
    ASW = max(ASW, 0);

    % HYSTER.f line 556-559
    if (ONE - RSW) > 1e-12
        ESAT = SRAW * ((ASW - RSW) / (ONE - RSW));
    else
        ESAT = 0;
    end
    ESW = ASW - ESAT;
    ESW = max(ESW, 0);
    SW  = ESW;

    % Mualem K (line 563-567)
    P0 = ONE - ((ONE - ASW^XXM)^m);
    if (ONE - RSW) > 1e-12
        CSATW = SRAW / (ONE - RSW);
    else
        CSATW = 0;
    end
    P3 = ((ONE - RSW^XXM)^m) - ((ONE - ASW^XXM)^m);
    P3 = CSATW * P3;
    PERMW = sqrt(max(ESW,0)) * ((P0 - P3)^2);

    % dθ/dh (line 571-574)
    if SWD > 0 && SWD < ONE
        X1 = SWD^XXM;
        X7 = max(0, ONE - CSATW);
        if abs(denom) < 1e-12
            X8 = 0;
        else
            X8 = (RASW_i - RASW_im1) / denom;
        end
        dSWdH = AlphaD * (n - ONE) * X1 * X7 * X8 * (ONE - X1)^m;
    else
        dSWdH = 0;
    end
end
