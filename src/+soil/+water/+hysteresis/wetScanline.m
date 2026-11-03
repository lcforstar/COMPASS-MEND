function [SW, ESW, ASW, dSWdH, PERMW] = wetScanline(HAW, AlphaW, AlphaD, n, m, ...
    RHSW_i, RHSW_im1, RASW_i, RASW_im1, SRAW, RSW)
% wetScanline - 湿润扫描线插值 (HYSTER.f WET subroutine)
%
% 输入:
%   AlphaW - 湿润 alpha = 2*AlphaD (K-P 关系)
%   AlphaD - 排水 alpha (用于 K 公式)
%   其他参数同 dryScanline
%
% 公式 (line 608-612):
%   ASW = ((SWI - RIDSWI) * (RASW[i] - RASW[i-1])) / (RDISWI - RIDSWI) + RASW[i-1]
% 其中:
%   RIDSWI = (1 + (αW·R[i-1])^n)^(-m)   ← 用 AlphaW
%   RDISWI = (1 + (αW·R[i])^n)^(-m)
%   SWI    = (1 + (αW·HAW)^n)^(-m)

    ZERO = 0.0;
    ONE = 1.0;
    XXM = ONE / m;

    % K-P: n_W = n, m_W = m
    RIDSWI = (ONE + (AlphaW * RHSW_im1)^n)^(-m);
    RDISWI = (ONE + (AlphaW * RHSW_i)^n)^(-m);
    SWI    = (ONE + (AlphaW * HAW)^n)^(-m);

    denom = RDISWI - RIDSWI;
    if abs(denom) < 1e-12
        ASW = RASW_im1;
    else
        ASW = ((SWI - RIDSWI) * (RASW_i - RASW_im1)) / denom + RASW_im1;
    end
    ASW = min(ASW, ONE);
    ASW = max(ASW, 0);

    if (ONE - RSW) > 1e-12
        ESAT = SRAW * ((ASW - RSW) / (ONE - RSW));
    else
        ESAT = 0;
    end
    ESW = ASW - ESAT;
    ESW = max(ESW, 0);
    SW  = ESW;

    % Mualem K (line 621-625, 与 dryScanline 同公式)
    P0 = ONE - ((ONE - ASW^XXM)^m);
    if (ONE - RSW) > 1e-12
        CSATW = SRAW / (ONE - RSW);
    else
        CSATW = 0;
    end
    P3 = ((ONE - RSW^XXM)^m) - ((ONE - ASW^XXM)^m);
    P3 = CSATW * P3;
    PERMW = sqrt(max(ESW,0)) * ((P0 - P3)^2);

    % dθ/dh (line 629-632)
    if SWI > 0 && SWI < ONE
        X1 = SWI^XXM;
        X7 = max(0, ONE - CSATW);
        if abs(denom) < 1e-12
            X8 = 0;
        else
            X8 = (RASW_i - RASW_im1) / denom;
        end
        dSWdH = AlphaW * (n - ONE) * X1 * X8 * X7 * (ONE - X1)^m;
    else
        dSWdH = 0;
    end
end
