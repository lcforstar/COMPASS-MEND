function [r1, r2] = quadraticRoots(a, b, c)
%QUADRATICROOTS  Numerically stable quadratic solver (numerical recipes form).
%
%   Avoids cancellation when b² >> 4ac.
%
%   Inputs:
%     a, b, c  [-]  ax² + bx + c = 0 coefficients
%
%   Returns:
%     r1, r2   [-]  two roots (1e36 if degenerate)
%
%   Reference:
%     CLM5 LunaMod.F90 line 1325-1352 (Quadratic) — 1:1 port
%
%   Author: Zhu Jizixing 2026-05-20

    r1 = 1.0e36;
    r2 = 1.0e36;

    if a == 0
        return;
    end

    if b >= 0.0
        q = -0.5 * (b + sqrt(b*b - 4.0*a*c));
    else
        q = -0.5 * (b - sqrt(b*b - 4.0*a*c));
    end

    r1 = q / a;
    if q ~= 0.0
        r2 = c / q;
    end
end
