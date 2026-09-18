function [esv, esh] = compute_rough_surface_emissivity(freq, theta, sig, cl, rv, rh, PRflag)

switch PRflag
    case 1
        [esh,esv]=matrix_PassiveRoughnessWigneron(sig,cl,rh,rv);
    case 2
        [esh,esv]=matrix_PassiveRoughnessWang(sig,rh,rv,theta,freq);
    case 3
        [esh,esv]=matrix_PassiveRoughnessLawrence(sig,cl,rh,rv,theta);
    case 4
        [esh,esv]=matrix_PassiveRoughnessZhao(sig,cl,rh,rv,theta);
    case 5
        [esh,esv]=matrix_PassiveRoughnessQHN(rh,rv,theta,freq);
end

end
%---------------------------PassiveRoughnessWigneron-----------------------
function [esh,esv]=matrix_PassiveRoughnessWigneron(sig,cl,rh,rv)
%    correct for surface roughness effect
%    J-P. Wigneron, L. Laguerre, and Y. H. Kerr, 2001.
%    A simple parameterization of the L-band microwave emission from rough agricultural soils
%    IEEE Transactions on Geoscience and Remote Sensing
%    Vol.39, No.8, 1697-1707.
%
%    Variables used in this subroutine:
%
%    Input
%    sig       - RMS height, in cm
%    cl        - correlation length, in cm
%    r0        - reflectivity at flat surface
%
%    Output
%    es        - emissivity at rough surface

sl=sig/cl;
% hs=1.3972*(sl).^0.5879;
hs=((0.9437*sig)./(0.8865*sig+2.2913)).^6; % Wigneron 2011
% hs=0.094;
esh=1-rh*exp(-hs);
esv=1-rv*exp(-hs);
% H=exp(-hs);
end
%---------------------------PassiveRoughnessWang---------------------------
function [esh,esv]=matrix_PassiveRoughnessWang(sig,rh,rv,theta,freq)
%    correct for surface roughness effect
%    J. R. Wang, and B. J. Choudhury, 1981.
%    Remote sensing of soil moisture content over bare field at 1.4 GHz frequency
%    Journal of Geophysical Research
%    Vol.86, No.C6, 5277-5282.
%
%    Variables used in this subroutine:
%
%    Input
%    rad       - incidence angle, in radian
%    theta     - incidence angle, in degree
%    sig       - RMS height, in cm
%    r0        - reflectivity at flat surface
%
%    Output
%    es        - emissivity at rough surface
rad=theta/180*pi;
ks=(pi*freq/15)*sig;
hs=(2*ks).^2;
% hs=(2*ks);
% hs=ks;
nr=2;
H=hs*cos(rad)^nr;
esh=1-rh*exp(-H);
esv=1-rv*exp(-H);

end

%---------------------------PassiveRoughnessLawrence-----------------------
function [esh,esv] = matrix_PassiveRoughnessLawrence(s,l,rh,rv,theta)
%Description: Code computes polarized emissivity for a random surface using
%the parametric model equations developed by Lawrence et al. (2013). Model
%is applicable to L-band only.

%Input Variables:
%s: rms heigh (cm)
%l: correlation length (cm)
%theta: Incidence angle (degrees)

%Output Product:
%e_v and e_h: emissivity for v and h polarizations.

thetar = theta * pi/180;

Zs = s^2 /l; %cm

if Zs <= 1.235
    h_pr = 2.615*(1 - exp(-Zs/2.473));
else
    h_pr = 1.0279;
end
Q = 0.1771 * h_pr;

n_v = 1.615 *(1-exp(-h_pr/0.359)) - 0.238; %for v-pol
n_h = 0.767* h_pr - 0.099; % for h-pol

gamma_rs_v = ( (1-Q)*rv + Q* rh) * exp(-h_pr * (cos(thetar))^(n_v));

gamma_rs_h = ( (1-Q)*rh + Q* rv) * exp(-h_pr * (cos(thetar))^(n_h));

esv = 1- gamma_rs_v;
esh = 1- gamma_rs_h;

end


%---------------------------PassiveRoughnessZhao-----------------------
function [esh,esv]=matrix_PassiveRoughnessZhao(sig,cl,rh,rv,theta)
%    correct for surface roughness effect, developed by Tianjie Zhao (2011)
%    for Exponential correlated surface
%    Variables used in this subroutine:
%    Input
%    rad       - incidence angle, in radian
%    theta     - incidence angle, in degree
%    polar     - polarization mode
%    sig       - RMS height, in cm
%    cl        - correlation length, in cm
%    r0        - reflectivity at flat surface
%    Output
%    es        - emissivity at rough surface
rad=theta/180*pi;
s2l=sig^2/cl;

%     ; H polarization
Ah= 0.068502*rad^2 -0.058486*rad +0.976321;
Bh=-0.051377*rad^2 +0.014978*rad +0.045456;
Ch= 0.601618*rad^2 -0.151848*rad -0.607679;
Hh=Ah*exp(Bh*s2l^2+Ch*s2l);
%     H=Hh;
esh=1.0-rh*Hh;

%     ; V polarization
Av= 0.296412*rad^2 -0.169075*rad +0.993309;
Bv=-0.049491*rad^2 -0.012277*rad +0.048694;
Cv= 0.552790*rad^2 +0.086331*rad -0.636948;
Hv=Av*exp(Bv*s2l.^2+Cv*s2l);
%     H=Hv;
esv=1.0-rv*Hv;

end

%---------------------------PassiveRoughnessQHN---------------------------
function [esh,esv]=matrix_PassiveRoughnessQHN(rh,rv,theta,freq)
%    correct for surface roughness effect
%    J. R. Wang, and B. J. Choudhury, 1981.
%    Remote sensing of soil moisture content over bare field at 1.4 GHz frequency
%    Journal of Geophysical Research
%    Vol.86, No.C6, 5277-5282.
%
%    Variables used in this subroutine:
%
%    Input
%    rad       - incidence angle, in radian
%    theta     - incidence angle, in degree
%    sig       - RMS height, in cm
%    r0        - reflectivity at flat surface
%
%    Output
%    es        - emissivity at rough surface

if freq>(1.41-0.5) && freq<(1.41+0.5)
    Q=0.0327;
    H=0.0967;
elseif freq>(6.925-0.5) && freq<(6.925+0.5)
    Q=0.2783;
    H=0.1042;
elseif freq>(10.65-0.5) && freq<(10.65+0.5)
    Q=0.3143;
    H=0.2018;
end

nr=2;

RH=((1-Q)*rh+Q*rv)*exp(-H*(cosd(theta))^nr);
RV=((1-Q)*rv+Q*rh)*exp(-H*(cosd(theta))^nr);

esh=1-RH;
esv=1-RV;

end