function [ev, eh] = compute_i2em_emissivity(fr, sig, L, theta_d, er, sp)
%I2EM_EMISSIVITY_MODEL  Improved Integral Equation Method rough-surface emissivity.
%
%   Returns:
%     ev, eh  [-]  Rough-surface emissivities for V and H polarizations
%
%   Inputs:
%     fr        [GHz]  Microwave frequency
%     sig       [m]    RMS surface height (function converts to cm internally)
%     L         [m]    Surface correlation length (function converts to cm)
%     theta_d   [deg]  Incidence angle
%     er        [-]    Complex dielectric constant of soil
%     sp        [-]    Spectrum: 1 = exponential, 2 = gaussian
%
%   Method:
%     I2EM is Fung et al. 2002's improved integral-equation rough surface
%     model. Numerically integrates the bistatic scattering coefficient over
%     the upper hemisphere and uses energy conservation to get emissivity.
%
%   References:
%     Fung, Liu, Chen, Tsang 2002 "An improved IEM model for bistatic
%       scattering from rough surfaces" J. Electromagn. Waves Appl. 16:689
%     Source: github.com/iganache/pvert (Public domain), with sig/L units
%       adapted to SI on input (matches PassiveRoughness signature).
%
%   Note: dblquad integration is slow (~50 ms per call). For full-year hourly
%     runs (8760 calls), expect ~7 min. Use I2EM only when σ > 1.5 cm where
%     Wigneron 2011 empirical fit degrades.

sig = sig * 100;   % m → cm
L   = L   * 100;   % m → cm

k = 2*pi *fr/30;                                  % wavenumber [cm^-1] (c=30 cm/ns)
theta = theta_d * pi/180;

ks  = k * sig;
kl  = k * L;
ks2 = ks * ks;

cs = cos(theta);
s  = sin(theta);
s2 = s*s;

sq  = sqrt(er - s2);
rv0 = (sqrt(er) - 1) / (sqrt(er) + 1);
rh0 = -rv0;                                        %#ok<NASGU>  reserved
rv  = (er*cs - sq) / (er*cs + sq);
rh  = (cs - sq) / (cs + sq);

pol = 1;
refv = dblquad(@(ths,phs) emsv_integralfunc(ths, phs, pol, sp, k, L, ks, cs, s, kl, er, sq, rv, rh), ...
               0.0, pi/2, 0, pi);

pol = 2;
refh = dblquad(@(ths,phs) emsv_integralfunc(ths, phs, pol, sp, k, L, ks, cs, s, kl, er, sq, rv, rh), ...
               0.0, pi/2, 0, pi);

ev = 1 - refv - exp(-ks2 * cs*cs) * (abs(rv))^2;
eh = 1 - refh - exp(-ks2 * cs*cs) * (abs(rh))^2;

ev = real(ev);
eh = real(eh);
end

% =========================================================================
function ref = emsv_integralfunc(ths, phs, pol, sp, k, L, ks, cs, s, kl, er, sq, rv, rh)
%EMSV_INTEGRALFUNC  Bistatic scattering integrand for I2EM emissivity (Fung 2002).

cs2 = cs^2;
ks2 = ks^2;
s2  = s*s;
kl2 = kl*kl;

css = cos(ths);
css2 = css.*css;
ss   = sin(ths);
ss2  = ss.*ss;
sf   = sin(phs);
csf  = cos(phs);

sqs = sqrt(er - ss2);
rc  = (rv - rh) / 2;
tv  = 1 + rv;
th  = 1 + rh;
tv2 = tv*tv;
th2 = th*th;

% surface correlation spectrum argument
wvnb = k .* sqrt(s2 - 2*s.*ss.*csf + ss2);

% number of spectral components (use smallest ths in vector)
costhsmin = cos(min(ths(:)));
n_spec = 1;
err_target = 1.0e-3;
while true
    n_spec = n_spec + 1;
    err_val = (ks2 * (cs + costhsmin)^2)^n_spec / factorial(n_spec);
    if err_val < err_target || n_spec > 30, break; end
end

nr = numel(ths);
wn = spectrm(sp, kl2, L, wvnb, n_spec, nr);

% Kirchhoff-like field coefficients
ff  = 2 .* (s.*ss - (1 + cs.*css).*csf) ./ (cs + css);
fvv =  rv .* ff;
fhh = -rh .* ff;
fvh = -2 * rc .* sf;
fhv =  2 * rc .* sf;

cm1 = s .* (ss - s.*csf) ./ (cs2 .* css);
T   = (sq .* (cs + sq) + cs .* (er*cs + sq)) ./ ...
      (er .* cs .* (cs + sq) + sq .* (er*cs + sq));
cm2 = css .* sq ./ cs ./ sqs - 1;
ex  = exp(-ks2 .* cs .* css);
de  = 0.5 .* exp(-ks2 .* (cs2 + css2));

if pol == 1
    Fvv = (er - 1) .* s2 .* tv2 .* cm1 ./ er.^2;
    Fhv = (T .* s.*s - 1 + cs./css + (er .* T .* cs .* css .* (er.*T - s.*s) - sq.*sq) ...
           ./ (T .* er .* sq .* css)) .* (1 - rc.*rc) .* sf;
    Fvvs = -cm2 .* sq .* tv2 .* (csf - s.*ss) ./ cs2 ./ er - cm2 .* sqs .* tv2 .* csf ./ er ...
           -(css .* sq ./ cs ./ sqs ./ er - 1) .* ss .* tv2 .* (s - ss.*csf) ./ cs;
    Fhvs = -(ss.*ss ./ T - 1 + css ./ cs + (cs .* css .* (1 - ss.*ss.*T) - T.*T.*sqs.*sqs) ...
           ./ (T .* sqs .* cs)) .* (1 - rc.*rc) .* sf;

    svv = zeros(n_spec, nr);
    for n = 1:n_spec
        Ivv = fvv .* ex .* (ks*(cs + css)).^n + (Fvv .* (ks*css).^n + Fvvs .* (ks*cs).^n) ./ 2;
        Ihv = fhv .* ex .* (ks*(cs + css)).^n + (Fhv .* (ks*css).^n + Fhvs .* (ks*cs).^n) ./ 2;
        wnn = wn(n,:) ./ factorial(n);
        vv  = wnn .* (abs(Ivv)).^2;
        hv  = wnn .* (abs(Ihv)).^2;
        svv(n,:) = de .* (vv + hv) .* ss ./ (4*pi*cs);
    end
    ref = sum(svv, 1);
else
    Fhh = -(er - 1) .* th2 .* cm1;
    Fvh = (s.*s ./ T - 1 + cs./css + (cs .* css .* (1 - s.*s.*T) - T.*T.*sq.*sq) ...
          ./ (T .* sq .* css)) .* (1 - rc.*rc) .* sf;
    Fhhs = cm2 .* sq .* th2 .* (csf - s.*ss) ./ cs2 + cm2 .* sqs .* th2 .* csf ...
           + cm2 .* ss .* th2 .* (s - ss.*csf) ./ cs;
    Fvhs = -(T .* ss.*ss - 1 + css ./ cs + (er .* T .* cs .* css .* (er*T - ss.*ss) - sqs.*sqs) ...
           ./ (T .* er .* sqs .* cs)) .* (1 - rc.*rc) .* sf;

    shh = zeros(n_spec, nr);
    for n = 1:n_spec
        Ihh = fhh .* ex .* (ks*(cs + css)).^n + (Fhh .* (ks*css).^n + Fhhs .* (ks*cs).^n) ./ 2;
        Ivh = fvh .* ex .* (ks*(cs + css)).^n + (Fvh .* (ks*css).^n + Fvhs .* (ks*cs).^n) ./ 2;
        wnn = wn(n,:) ./ factorial(n);
        hh  = wnn .* (abs(Ihh)).^2;
        vh  = wnn .* (abs(Ivh)).^2;
        shh(n,:) = de .* (hh + vh) .* ss ./ (4*pi*cs);
    end
    ref = sum(shh, 1);
end
end

% =========================================================================
function wn = spectrm(sp, kl2, L, wvnb, np, nr)                              %#ok<INUSL>
%SPECTRM  Surface correlation power spectrum (Fung 1994 §3.3).
wn = zeros(np, nr);
if sp == 1     % exponential
    for n = 1:np
        wn(n,:) = n * kl2 ./ (n^2 + (wvnb.*L).^2).^1.5;
    end
elseif sp == 2 % gaussian
    for n = 1:np
        wn(n,:) = 0.5 * kl2 ./ n .* exp(-(wvnb.*L).^2 ./ (4*n));
    end
end
end
