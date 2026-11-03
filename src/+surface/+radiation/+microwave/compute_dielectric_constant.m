function dc = compute_dielectric_constant(freq,sand,clay,bd,temp,sm,DCflag)

switch DCflag
    case 1
        dc=matrix_diel_mironov(freq,sm,clay);
    case 2
        dc=matrix_diel_dobson(freq,sm,temp,sand,clay);
    case 3
        dc=matrix_diel_wang(freq,sm,temp,sand,clay);
    case 4
        dc=matrix_diel_zhang(freq,sm,temp,sand,clay,bd);
end

end


%-----------------------------------diel_mironov---------------------------
function dc=matrix_diel_mironov(freq,sm,clay)

%   ; Purpose :
%   ;   Calculate the dielectric constant of a wet soil
%   ;   Developed and validated from 1 to 10 GHz.
%   ;   adapted for a large range of soil moisture
%
%   ; Reference:
%   ;   Mironov et al: Generalized Refractive Mixing Dielectric Model for moist soil
%   ;    IEEE Trans. Geosc. Rem. Sens., vol 42 (4), 773-785. 2004.

eps_winf=4.9 ; %constant of high-frequency limit
eps_0=8.854e-12 ; %dielectric constant of free space
f=freq*1.0e9;
% Initializing the GRMDM spectroscopic parameters with clay (fraction)

%  RI & NAC of dry soils
znd=1.634-0.539*clay+0.2748*clay^2.0;
zkd=0.03952-0.04038*clay;

% Maximum bound water fraction
zxmvt=0.02863+0.30673*clay;

% Bound water parameters
zep0b=79.8-85.4*clay+32.7*clay^2.0;
ztaub=1.062e-11+3.450e-12*clay;
zsigmab=0.3112+0.467*clay;

% Unbound (free) water parameters
zep0u=100.0;
ztauu=8.5e-12;
zsigmau=0.3631+1.217*clay;

% Computation of epsilon water (bound & unbound)

zcxb=(zep0b-eps_winf)/(1.0+(2.0*pi*f*ztaub)^2.0);
zepwbx=eps_winf+zcxb;
zepwby=zcxb*(2.0*pi*f*ztaub)+zsigmab/(2.0*pi*eps_0*f);

zcxu=(zep0u-eps_winf)/(1.0+(2.0*pi*f*ztauu)^2.0);
zepwux=eps_winf+zcxu;
zepwuy=zcxu*(2.0*pi*f*ztauu)+zsigmau/(2.0*pi*eps_0*f);

% Computation of refractive index of water (bound & unbound)

znb=sqrt(sqrt(zepwbx^2.0+zepwby^2.0)+zepwbx)/sqrt(2.0);
zkb=sqrt(sqrt(zepwbx^2.0+zepwby^2.0)-zepwbx)/sqrt(2.0);

znu=sqrt(sqrt(zepwux^2.0+zepwuy^2.0)+zepwux)/sqrt(2.0);
zku=sqrt(sqrt(zepwux^2.0+zepwuy^2.0)-zepwux)/sqrt(2.0);

% Computation of soil refractive index (nm & km):
zxmvt2=zxmvt;
if (sm < zxmvt)
    zxmvt2=sm;
end
zflag=0.0;
if (sm >= zxmvt)
    zflag=1.0;
end

znm=znd+(znb-1.0)*zxmvt2+(znu-1.0)*(sm-zxmvt)*zflag;
zkm=zkd+zkb*zxmvt2+zku*(sm-zxmvt)*zflag;

% computation of soil dielectric constant:

zepmx=znm^2.0-zkm^2.0;
zepmy=znm*zkm*2.0;

dc=complex(zepmx,zepmy);

end


%-----------------------------------diel_dobson----------------------------
function dc=matrix_diel_dobson(freq,sm,temp,sand,clay)

%   The soil texture is in decimal
%   Bulk density could be calculated from soil texture
%   There is a correction when frequency is lower than 1.4 GHz

fi=double(freq) ; % frequency (GHz)
wc=double(sm) ; % soil moisture (m3/m3)
ts=double(temp) ; % soil temperature (Degrees Celsius)
sf=double(sand) ; % sand fraction (decimal)
cf=double(clay) ; % clay fraction (decimal)
% %   fi, frequency, GHz
%   fi=float(fi);
% %   sm, soil moisture, m3/m3
%   wc=float(sm);
% %   ts, soil temperature, Degrees Celsius
%   ts=float(ts);
% %   sf, sand fraction, decimal
%   sf=float(sf);
% %   cf, clay fraction, decimal
%   cf=float(cf);

%   Set following values
alpha=0.65;
epso=8.854e-12 ; %dielectric constant of free space
epswinf=4.9 ; %constant of high-frequency limit

%   Compute dielectric constant of soil solid particles
rhos=2.66 ; %soil specific density
epss=(1.01+0.44*rhos)^2-0.062 ; %dielectric constant of soil solid particles, 4.69
por=0.505-0.142*sf-0.037*cf;
fv=1.0-por;
rhob=fv*rhos ; %soil bulk density

%   Compute optimized coefficient values
betar=1.27480-0.519*sf-0.152*cf;
betai=1.33797-0.603*sf-0.166*cf;

%   Compute empirical expression for effective conductivity
sigmae=-1.645+1.939*rhob-2.25622*sf+1.594*cf ; %2.25622 > 2.013 ?
if (fi < 1.4)
    sigmae=0.0467+0.2204*rhob-0.4111*sf+0.6614*cf;
end

%   Compute dielectric constant of free water (expressions for
%   relaxation time and static dielectric constant obtained from
%   Ulaby et al., Vol. 3)
omtau=fi*(1.1109e-1-3.824e-3*ts+6.938e-5*ts^2-5.096e-7*ts^3) ; %relaxation time of pure water
epswo=88.045-0.4147*ts+6.295e-4*ts^2+1.075e-5*ts^3 ; %static dielectric constant of pure water
fac=(epswo-epswinf)/(1.0+omtau^2);
epsfwr=epswinf+fac ; %real part of dielectric constant of water
epsfwi=(omtau*fac)+(sigmae*(rhos-rhob)/(2.0*pi*fi*epso*1.0e9*rhos*wc)) ; %imag part of dielectric constant of water

%   Compute dielectric constant of soil
dr=(1.0+(epss^alpha-1.0)*rhob/rhos+(wc^betar)*(epsfwr^alpha)-wc)^(1.0/alpha);
di=((wc^betai)*(epsfwi^alpha))^(1.0/alpha);

%   Correction factor (Peplinski et al., 1995)
if (fi < 1.4)
    dr=1.15*dr-0.68;
end

dc=complex(dr,di);

end


%-----------------------------------diel_wang------------------------------
function dc=matrix_diel_wang(freq,sm,temp,sand,clay)

%    Purpose :
%      Calculate the dielectric constant of a wet soil
%      Developed and validated for 1.4 and 5 GHz.
%
%   Reference:
%      Wang and Schmugge, 1980: An empirical model for the
%      complex dielectric permittivity of soils as a function of water
%      content. IEEE Trans. Geosci. Rem. Sens., GE-18, No. 4, 288-295.
%
%      wt    : transition moisture point (m3/m3)
%      rhob  : soil bulk density (g/cm3)
%      gamma : fitting parameter
%      ecl   : conductivity loss
%      ei    : dielectric constant of ice, for initially absorbed water
%      ex    : dielectric constant of the initially absorbed water

fi=double(freq) ; % frequency (GHz)
wc=double(sm) ; % soil moisture (m3/m3)
ts=double(temp) ; % soil temperature (Degrees Celsius)
sf=double(sand) ; % sand fraction (decimal)
cf=double(clay) ; % clay fraction (decimal)

ea=complex(1.0,0.0) ; % dielectric constant of air
er=complex(5.5,0.2) ; % dielectric constant of rock
ei=complex(3.2,0.1) ; % dielectric constant of ice

rhos=2.66 ; % soil specific density (g/cm3)
epso=8.854e-12 ; % dielectric constant of free space
epswinf=4.9 ; % constant of high-frequency limit
p=0.505-0.142*sf-0.037*cf;
%    p=1-rhob/rhos ; soil porosity

%    Calculate dielectric constant of free water
omtau=fi*(1.1109e-1-3.824e-3*ts+6.938e-5*ts^2-5.096e-7*ts^3) ; % relaxation time of pure water
epswo=88.045-0.4147*ts+6.295e-4*ts^2+1.075e-5*ts^3 ; % static dielectric constant of pure water
fac=(epswo-epswinf)/(1.0+omtau^2);
epsfwr=epswinf+fac ; % real part of dielectric constant of water
epsfwi=omtau*fac ; % imag part of dielectric constant of water
ew=complex(epsfwr,epsfwi);

%   Calculate dielectric constant of soil-water mixture
wp=0.06774-0.064*sf+0.478*cf;
gamma=-0.57*wp+0.481;
wt=0.49*wp+0.165;
if (wc <= wt)
    ex=ei+(ew-ei)*(wc/wt)*gamma;
    eps=wc*ex+(p-wc)*ea+(1.0-p)*er;
else
    ex=ei+(ew-ei)*gamma;
    eps=wt*ex+(wc-wt)*ew+(p-wc)*ea+(1.0-p)*er;
end
%   Add conductivity loss
%    ecl=alpha*wc^2.0
%   eps=eps+complex(0.0,1.0)*ecl

dc=eps;

end

function  dc = matrix_diel_zhang (fre, vms, ts, sand, clay, bd)

% sandc: sand content(%)
% clayc: clay content(%)
% s: specific surface area of soil mixture with sand, silt and clay
% bs: specific density
% bd: bulk density
% vmsu: volumic unfrozen soil moisture

sandc=sand*100;
clayc=clay*100;
bs=2.60;
siltc=100-sandc-clayc;

s=0.042+4.23*clayc+1.12*siltc-1.16*sandc; % modified by Tianjie
volsolid=bd/bs;
ice_epson=3.15;
alpha=0.65;
dcsr=4.7;

% a=5.2752; b=0.5675
% The above two parameters are  related to the texure of soil and dependtent on "cv" and "sv"

if (ts < 0)    
    vmsu=exp( 0.2618 + 0.5519*log(s) - exp(-0.264*log(s)+0.3711)*log(abs(ts)) ) * bd /100.0;
    vi=(vms-vmsu)/0.9;    
    if (vi < 0)
        vi=0;
        vmsu=vms;
    end    
else    
    vi=0;
    vmsu=vms;    
end

% set constant of high-frequency limit (dimensionless)
ewin = 4.9;

% calculate texture dependent parameters
betar = 1.2748 - 0.00519 * sandc - 0.00152 * clayc;
betai = 1.33797 - 0.00603 * sandc - 0.00166 * clayc;
sigmae = -1.645 + 1.939 * bd - 0.02013 * sandc + 0.01594 * clayc;

% calculate temperature dependence of static dielectric constant of pure
% water (dimensionless)
ew0 = 88.045 - 0.4147 * ts + 6.295 * 0.0001 * ts * ts + 1.075 * 0.00001 * ts * ts * ts;

% calculate temperature dependence of relaxation time of pure water
% where fre in frequncy in GHz
tauf = fre * (0.11109 - 0.003824 * ts + 0.00006938 * ts * ts - 5.096 * 0.0000001 * ts * ts * ts);

% calculate real part of dielectric constant of water
rdcw = ewin + ( ew0 - ewin ) / ( 1.0 + tauf * tauf );

% calculate imag part of dielectric constant of water
idcw = tauf * ( ew0 - ewin ) / ( 1.0 + tauf * tauf ) ...
    + sigmae * ( 1.0 - volsolid ) / (8.0 * atan(1.0) * 8.854 * 0.001 * fre * vmsu);

tmp = 1.0 + volsolid * ((dcsr^alpha) - 1.0 ) ...
    + (vmsu^betar) * (rdcw^alpha) ...
    -vmsu+vi*(ice_epson^alpha)-vi;

rdcs = (tmp^(1.0/alpha));
idcs = (vmsu^(betai/alpha)) * idcw;

dc = complex(rdcs,idcs);

end

