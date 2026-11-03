function dtheta = vgDerivative(h, theta_s, theta_r, alpha, n, m)
%VGDERIVATIVE  Derivative dtheta/dh for the van Genuchten retention curve.
%
%   dtheta = soil.water.vgDerivative(h, theta_s, theta_r, alpha, n, m)
%   computes the positive derivative of volumetric water content with respect
%   to matric head in COMPASS native units.
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
%     dtheta   [out] [1/cm] dtheta/dh
%
%   Reference:
%     van Genuchten 1980 SSSAJ 44:892 Eq. 21 differentiated analytically.
%
%   Editor: Jizixing Zhu; Date: 2026-06-28.

    if h >= -1e-6 || h <= -1e7
        dtheta = 0;
        return;
    end
    x = abs(alpha * h);
    dtheta = (theta_s - theta_r) * m * n * alpha * x^(n - 1) * ...
        (1 + x^n)^(-m - 1);
    dtheta = max(dtheta, 0);
end
