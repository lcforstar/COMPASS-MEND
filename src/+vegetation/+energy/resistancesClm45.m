function [resist_out] = resistancesClm45(resist_in)
%{
 CLM4.5-style aerodynamic resistance branch for ebal_clm45.

 Scientific role:
   This function preserves the original SCOPE/Wallace-Verhoef within-canopy
   resistance network, but augments the above-canopy and boundary-layer
   terms with CLM4.5/CLM5-style parameterizations. It is a hybrid SCOPE plus
   CLM4.5 aerodynamic scheme, not a pure CLM4.5 big-leaf resistance module.

 CLM4.5-style improvement:
   Jizixing Zhu (朱寄子星), 18 Apr 2026.

 Main updates:
   - CLM4.5 four-regime Monin-Obukhov stability corrections.
   - Stanton-type scalar roughness lengths z0h and z0w.
   - Explicit heat and water-vapor aerodynamic resistances rah and raw.
   - Choudhury-Monteith style canopy boundary-layer resistance rbc.
   - Dynamic leaf-width scaling driven by SCOPE RTMo direct-beam k.

 Retained Wallace-Verhoef network:
   canopy path: canopy leaf -> rawc -> raa -> atmosphere
   soil path:   soil surface -> raws -> raa -> atmosphere

 Method document:
   docs/clm45_aerodynamic_resistance_method.md

 Editor: Jizixing Zhu; Date: 2026-07-20.
%}
%
%   function resistances calculates aerodynamic and boundary resistances
%   for soil and vegetation
%
%   Date:       01 Feb 2008
%   Authors:    Anne Verhoef            (a.verhoef@reading.ac.uk)
%               Christiaan van der Tol  (tol@itc.nl)
%               Joris Timmermans        (j_timmermans@itc.nl)
%   Source:     Wallace and Verhoef (2000) 'Modelling interactions in
%               mixed-plant communities: light, water and carbon dioxide', in: Bruce
%               Marshall, Jeremy A. Roberts (ed), 'Leaf Development and Canopy Growth',
%               Sheffield Academic Press, UK. ISBN 0849397693
%
%               ustar:  Tennekes, H. (1973) 'The logaritmic wind profile', J.
%               Atmospheric Science, 30, 234-238
%               Psih:   Paulson, C.A. (1970), The mathematical
%               representation of wind speed and temperature in the
%               unstable atmospheric surface layer. J. Applied Meteorol. 9,
%               857-861
%
% Note: Equation numbers refer to equation numbers in Wallace and Verhoef (2000)
%
% Usage:
%   [resist_out] = resistances(resist_in)
%
% The input and output are structures. These structures are further
% specified in a readme file.
%
% Input:
%   resist_in   aerodynamic resistance parameters and wind speed
%
% The strucutre resist_in contains the following elements:
% u         =   windspeed
% L         =   stability
% LAI       =   Leaf Area Index

% rbs       =   Boundary Resistance of soil                         [s m-1]
% rss       =   Surface resistance of soil for vapour transport     [s m-1]
% rwc       =   Within canopy Aerodynamic Resistance canopy         [s m-1]

% z0m       =   Roughness lenght for momentum for the vegetation    [m]
% d         =   Displacement height (Zero plane)                    [m]
% z         =   Measurement height                                  [m]
% h         =   Vegetation height                                   [m]

%
% Output:
%   resist_out  aeorodynamic resistances
%
% The strucutre resist_out contains the following elements:
% ustar     =   Friction velocity                                   [m s-1]
% raa       =   Aerodynamic resistance above the canopy             [s m-1]
% rawc      =   Total resistance within the canopy (canopy)         [s m-1]
% raws      =   Total resistance within the canopy (soil)           [s m-1]

% rai       =   Aerodynamic resistance in inertial sublayer         [s m-1]
% rar       =   Aerodynamic resistance in roughness sublayer        [s m-1]
% rac       =   Aerodynamic resistance in canopy layer (above z0+d) [s m-1]

% rbc       =   Boundary layer resistance (canopy)                  [s m-1]
% rwc       =   Aerodynamic Resistance within canopy(canopy)(Update)[s m-1]

% rbs       =   Boundary layer resistance (soil) (Update)           [s m-1]
% rws       =   Aerodynamic resistance within canopy(soil)          [s m-1]

% rss       =   Surface resistance vapour transport(soil)(Update)   [s m-1]

% uz0       =   windspeed at z0                                     [m s-1]
% Kh        =   Diffusivity for heat                                [m2s-1]

% load physical constants
phys = equations.physicalConstants();
kappa = phys.k_VK;

%% parameters

Cd        =  resist_in.Cd;

u         =  resist_in.u;
L         =  resist_in.L;
LAI       =  resist_in.LAI;

rbs       =  resist_in.rbs;
% rss       =  resist_in.rss;
rwc       =  resist_in.rwc;

z0m       =  resist_in.zo;
d         =  resist_in.d;
z         =  resist_in.z;
h         =  resist_in.hc;
w         =  resist_in.w;

LAI_max = 4.0;
w_leaf_max = w;
CropCode = 1;
Tao = 0.6;
if isfield(resist_in, 'LAI_max') && ~isempty(resist_in.LAI_max), LAI_max = resist_in.LAI_max; end
if isfield(resist_in, 'w_leaf_max') && ~isempty(resist_in.w_leaf_max), w_leaf_max = resist_in.w_leaf_max; end
if isfield(resist_in, 'CropCode') && ~isempty(resist_in.CropCode), CropCode = resist_in.CropCode; end
if isfield(resist_in, 'Tao') && ~isempty(resist_in.Tao), Tao = resist_in.Tao; end
resist_out.Tao = Tao;       % extinction coefficient used by model-calculated leaf-width scaling

if h <= 0 || z0m <= 0 || z <= d + z0m
    error('Aerodynamic:InvalidReferenceHeight', ...
        ['Invalid aerodynamic geometry: measurement height z must be above ', ...
        'displacement height plus roughness length (z > d + z0m). ', ...
        'Got z=%.3f m, d=%.3f m, z0m=%.3f m, hc=%.3f m.'], ...
        z, d, z0m, h);
end
% [改进 2026-04-18 朱寄子星] 冠层接近测风高度时的预警（只提醒一次）
% Warn once when z < 2·hc (geometry borderline, z > d+z0m still holds but margin thin)
persistent tallCanopyWarned
if isempty(tallCanopyWarned), tallCanopyWarned = false; end
if z < 2 * h && ~tallCanopyWarned
    warning('Aerodynamic:TallCanopyRisk', ...
        ['测风高度偏低: z=%.2f m 仅为 hc=%.2f m 的 %.1fx，' ...
         '推荐 z ≥ 2·hc_max。冠层继续长高会触发 InvalidReferenceHeight 错误。' ...
         '建议 config_environment.txt 里加大 Z_measure。'], ...
        z, h, z/max(h,1e-6));
    tallCanopyWarned = true;
end

% derived parameters
% zr: top of roughness sublayer, bottom of intertial sublayer
zr          = 2.5 * h;                   %                            [m]
% n: dimensionless wind extinction coefficient                       W&V Eq 33
n           = Cd * LAI / (2 * kappa^2);      %                            []

%% stability correction for non-neutral conditions
% neu        = find(L >= -.001 & L <= .001);
unst        = find(L <  -4);
st          = find(L >  4E3);

% stability correction functions, friction velocity and Kh=Km=Kv
pm_z        = psim(z - d, L, unst, st);
ph_z        = psih(z - d, L, unst, st);
pm_h        = psim(h - d, L, unst, st);
% ph_h       = psih(h -d,L,unst,st);
ph_zr       = psih(zr - d, L, unst, st) .* (z >= zr) + ph_z .* (z < zr);
phs_zr      = phstar(zr, zr, d, L, st, unst);
phs_h       = phstar(h, zr, d, L, st, unst);

ustar       = max(.001, kappa * u ./ (log((z - d) / z0m) - pm_z)); %          W&V Eq 30

zeta_kh     = safe_zeta_clm45((h - d) .* ones(size(L)), L);
phi_h_kh    = phi_h_clm45(zeta_kh);
Kh          = kappa * ustar * (zr - d) ./ phi_h_kh;                  % W&V Eq 35, CLM4.5 four-regime stability
resist_out.Kh = Kh;

%% wind speed at height h and z0m
uh          = max(ustar / kappa .* (log((h - d) / z0m) - pm_h), .01);
uz0         = uh * exp(n * ((z0m + d) / h - 1));                      %       W&V Eq 32

%% resistances

resist_out.uz0 = uz0;
resist_out.ustar = ustar;
resist_out.uh = uh;             % [m s-1] wind speed at canopy top
resist_out.n_canopy = n;        % [-] W&V canopy wind/K extinction coefficient
rai = (z > zr) .* (1 ./ (kappa * ustar) .* (log((z - d) / (zr - d))  - ph_z   + ph_zr)); % W&V Eq 41
rar = 1 ./ (kappa * ustar) .* ((zr - h) / (zr - d))      - phs_zr + phs_h; % W&V Eq 39
rac = h * sinh(n) ./ (n * Kh) * (log((exp(n) - 1) / (exp(n) + 1)) - log((exp(n * (z0m + d) / h) - 1) / (exp(n * (z0m + d) / h) + 1))); % W&V Eq 42
rws = h * sinh(n) ./ (n * Kh) * (log((exp(n * (z0m + d) / h) - 1) / (exp(n * (z0m + d) / h) + 1)) - log((exp(n * (.01) / h) - 1) / (exp(n * (.01) / h) + 1))); % W&V Eq 43
rbc_wv = 70 / LAI * sqrt(w ./ uz0);                        % previous W&V empirical canopy boundary resistance

n_bl = 2.5 .* ones(size(h));
n_bl(h > 1 & h < 10) = 2.306 + 0.194 .* h(h > 1 & h < 10);
n_bl(h >= 10) = 4.25;

U_h_bl = (ustar ./ kappa) .* log((h - d) ./ z0m);
U_h_bl = max(0.01, U_h_bl);

LAI_safe = max(LAI, 1e-6);
w_leaf_actual = zeros(size(LAI_safe));
if isscalar(w_leaf_max), w_leaf_max = repmat(w_leaf_max, size(LAI_safe)); end
if isscalar(CropCode), CropCode = repmat(CropCode, size(LAI_safe)); end

mask_static = (CropCode == 1);
w_leaf_actual(mask_static) = w_leaf_max(mask_static);

mask_dynamic = ~mask_static;
if any(mask_dynamic)
    % CropCode==1 uses user-provided w_leaf_max directly; all other values
    % use model-calculated leaf width driven by Tao, normally SCOPE RTMo gap.k.
    if isscalar(Tao)
        exponent_term = -Tao .* LAI_safe(mask_dynamic);
    else
        exponent_term = -Tao(mask_dynamic) .* LAI_safe(mask_dynamic);
    end
    w_leaf_actual(mask_dynamic) = w_leaf_max(mask_dynamic) .* (1 - exp(exponent_term));
end
w_leaf_actual = max(w_leaf_actual, 1e-6);

r_b = (100 ./ n_bl) .* (w_leaf_actual ./ U_h_bl).^0.5 ./ (1 - exp(-n_bl ./ 2));

LAI_max_val = LAI_max;
if isempty(LAI_max_val), LAI_max_val = 4.0; end
threshold_val = 0.5 .* LAI_max_val;
denom = 2 .* LAI_safe;
idx_high = LAI_safe > threshold_val;
if isscalar(LAI_max_val)
    denom(idx_high) = 0.5 .* LAI_max_val;
else
    denom(idx_high) = 0.5 .* LAI_max_val(idx_high);
end
denom = max(denom, 1e-6);

rbc = r_b ./ denom;
rbc(LAI < 0.05) = 0;

resist_out.rai = rai;
resist_out.rar = rar;
resist_out.rac = rac;
resist_out.rws = rws;
resist_out.rbc = rbc;
resist_out.rbc_wv = rbc_wv;
resist_out.r_b = r_b;
resist_out.rbc_denom = denom;
resist_out.w_leaf_actual = w_leaf_actual;
resist_out.n_bl = n_bl;
resist_out.U_h_bl = U_h_bl;

raa  = rai + rar + rac;
rawc = rwc + rbc;
raws = rws + rbs;

resist_out.raa  = raa;          % aerodynamic resistance above the canopy           W&V Figure 8.6
resist_out.rawc = rawc;         % aerodynamic resistance within the canopy (canopy)
resist_out.raws = raws;         % aerodynamic resistance within the canopy (soil)

resist_out.raa  = min(4E2, raa);         % to prevent unrealistically high resistances
resist_out.rawc = min(4E2, rawc);        % to prevent unrealistically high resistances
resist_out.raws = min(4E2, raws);        % to prevent unrealistically high resistances

% CLM4.5-style scalar aerodynamic resistance above the canopy.
% Original Wallace-Verhoef terms are retained for the within-canopy network.
nu_air = 1.5e-5;      % [m2 s-1]
z0m_clm = max(z0m, 1e-6);
zref_eff = max(z - d, 1.01 .* z0m_clm);
Va = u;     % u 已在 MeteoData 入口钳到 max(1.0, u)，此处不再二次钳位

zeta_ref = safe_zeta_clm45(zref_eff, L);
zeta_m0  = safe_zeta_clm45(z0m_clm, L);
psi_m_ref = psi_m_clm45(zeta_ref);
psi_m_0   = psi_m_clm45(zeta_m0);
fm = max(0.1, log(zref_eff ./ z0m_clm) - psi_m_ref + psi_m_0);
ustar_clm = max(0.001, kappa .* Va ./ fm);

z0h = z0m_clm .* exp(-0.13 .* max(ustar_clm .* z0m_clm ./ nu_air, 0).^0.45);
z0h = min(max(z0h, 1e-7), 0.99 .* z0m_clm);
z0w = z0h;

psi_h_ref = psi_h_clm45(zeta_ref);
zeta_h0 = safe_zeta_clm45(z0h, L);
psi_h_0 = psi_h_clm45(zeta_h0);
fh = max(0.1, log(zref_eff ./ z0h) - psi_h_ref + psi_h_0);

zeta_w0 = safe_zeta_clm45(z0w, L);
psi_w_0 = psi_h_clm45(zeta_w0);
fw = max(0.1, log(zref_eff ./ z0w) - psi_h_ref + psi_w_0);

resist_out.rah = min(4E2, max(1, fm .* fh ./ (kappa.^2 .* Va)));
resist_out.raw = min(4E2, max(1, fm .* fw ./ (kappa.^2 .* Va)));
resist_out.z0h = z0h;
resist_out.z0w = z0w;
resist_out.Va = Va;
resist_out.ustar_clm = ustar_clm;

return

%% subfunction pm for stability correction (CLM4.5 four-regime)
function pm = psim(z, L, ~, ~)
    zeta        = safe_zeta_clm45(z .* ones(size(L)), L);
    pm          = psi_m_clm45(zeta);
    return

    %% subfunction ph for stability correction (CLM4.5 four-regime)
function ph = psih(z, L, ~, ~)
    zeta        = safe_zeta_clm45(z .* ones(size(L)), L);
    ph          = psi_h_clm45(zeta);
    return

    %% subfunction ph for roughness sublayer correction (CLM4.5 four-regime)
function phs = phstar(z, zR, d, L, ~, ~)
    zeta        = safe_zeta_clm45((z - d) .* ones(size(L)), L);
    phs         = (z - d) ./ max(zR - d, 1e-6) .* psi_h_clm45(zeta);
    return

function zeta = safe_zeta_clm45(z_over_d, L)
    zeta = zeros(size(L));
    neutral = abs(L) < 1e-6 | ~isfinite(L);
    zeta(~neutral) = z_over_d(~neutral) ./ L(~neutral);
    zeta = max(-100, min(2, zeta));
    return

function psi = psi_m_clm45(zeta)
    psi = zeros(size(zeta));
    zc = -1.574;
    very_unstable = zeta < zc;
    unstable = zeta >= zc & zeta < 0;
    stable = zeta >= 0 & zeta <= 1;
    very_stable = zeta > 1;

    if any(very_unstable)
        z = zeta(very_unstable);
        xcrit = (1 - 16 * zc).^0.25;
        psicrit = 2*log((1+xcrit)/2) + log((1+xcrit.^2)/2) - 2*atan(xcrit) + pi/2;
        psi(very_unstable) = -log(zc ./ z) + psicrit - 1.14 .* ((-z).^(1/3) - (-zc).^(1/3));
    end
    if any(unstable)
        z = zeta(unstable);
        x = (1 - 16 .* z).^0.25;
        psi(unstable) = 2*log((1+x)/2) + log((1+x.^2)/2) - 2*atan(x) + pi/2;
    end
    if any(stable)
        psi(stable) = -5 .* zeta(stable);
    end
    if any(very_stable)
        z = zeta(very_stable);
        psi(very_stable) = -(4 .* log(z) + z + 4);
    end
    return

function psi = psi_h_clm45(zeta)
    psi = zeros(size(zeta));
    zc = -0.465;
    very_unstable = zeta < zc;
    unstable = zeta >= zc & zeta < 0;
    stable = zeta >= 0 & zeta <= 1;
    very_stable = zeta > 1;

    if any(very_unstable)
        z = zeta(very_unstable);
        xcrit = (1 - 16 * zc).^0.25;
        psicrit = 2*log((1+xcrit.^2)/2);
        psi(very_unstable) = -log(zc ./ z) + psicrit - 0.8 .* ((-zc).^(-1/3) - (-z).^(-1/3));
    end
    if any(unstable)
        z = zeta(unstable);
        x = (1 - 16 .* z).^0.25;
        psi(unstable) = 2*log((1+x.^2)/2);
    end
    if any(stable)
        psi(stable) = -5 .* zeta(stable);
    end
    if any(very_stable)
        z = zeta(very_stable);
        psi(very_stable) = -(4 .* log(z) + z + 4);
    end
    return

function phi = phi_h_clm45(zeta)
    phi = ones(size(zeta));
    zc = -0.465;
    very_unstable = zeta < zc;
    unstable = zeta >= zc & zeta < 0;
    stable = zeta >= 0 & zeta <= 1;
    very_stable = zeta > 1;

    if any(very_unstable)
        z = zeta(very_unstable);
        phi(very_unstable) = 0.9 .* max(-z, 1e-6).^(-1/3);
    end
    if any(unstable)
        z = zeta(unstable);
        phi(unstable) = (1 - 16 .* z).^(-0.5);
    end
    if any(stable)
        phi(stable) = 1 + 5 .* zeta(stable);
    end
    if any(very_stable)
        z = zeta(very_stable);
        phi(very_stable) = 5 + z;
    end
    phi = max(phi, 0.05);
    return
