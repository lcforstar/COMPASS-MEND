function theta = thetaFromHead(h, theta_s, theta_r, alpha, n, m)
%THETAFROMHEAD  van Genuchten water content from matric head.
%
%   theta = soil.water.thetaFromHead(h, theta_s, theta_r, alpha, n, m)
%   computes volumetric water content from matric head with saturated and dry
%   guards used by the COMPASS soil-water auxiliary refresh.
%
%   Inputs:
%     h        [in] [cm]   matric head
%     theta_s  [in] [-]    saturated water content
%     theta_r  [in] [-]    residual water content
%     alpha    [in] [1/cm] van Genuchten alpha
%     n        [in] [-]    van Genuchten n
%     m        [in] [-]    van Genuchten m
%
%   Output:
%     theta    [out] [-]   volumetric water content
%
%   Reference:
%     van Genuchten 1980 SSSAJ 44:892 Eq. 21.
%
%   Editor: Jizixing Zhu; Date: 2026-06-28.

    if h >= -1e-6
        theta = theta_s;
    elseif h <= -1e7
        theta = theta_r;
    else
        theta = theta_r + (theta_s - theta_r) / (1 + abs(alpha * h)^n)^m;
    end
end
