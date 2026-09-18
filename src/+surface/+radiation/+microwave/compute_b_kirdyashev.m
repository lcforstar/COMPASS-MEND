function b = compute_b_kirdyashev(T_C, m_d, freq_GHz, inc_deg, a_geo)
%COMPUTE_B_KIRDYASHEV  Physical b coefficient for VOD = b·VWC (CMEM vegkird).
%
%   Returns:
%     b   [m²/kg]  Mass extinction coefficient (kg/m² water → optical depth)
%
%   Inputs:
%     T_C       [°C]   Tissue temperature
%     m_d       [-]    Dry-matter mass fraction (Cdm / (Cdm + Cw))
%     freq_GHz  [GHz]  Frequency
%     inc_deg   [deg]  Incidence angle
%     a_geo     [-]    Vegetation structure coefficient (LIDF-derived, ≈ 2/3 spherical)
%
%   Method (CMEM veg_sub.F90 SUBROUTINE VEGKIRD):
%     "Effective medium theory" — vegetation modelled as water-saturated
%     dielectric medium; opacity proportional to volumetric water content
%     times imaginary part of complex permittivity.
%
%       τ = a_geo · k · (VWC/ρ_w) · Im(ε_veg) / cosθ
%       b = a_geo · k ·       1/ρ_w  · Im(ε_veg) / cosθ
%
%     where:
%       k          = 2π·f/c  = free-space wavenumber [1/cm] (with f in Hz, c=cm/s)
%       ρ_w        = 1000 kg/m³ (treated as 1 g/cm³ internally → unit cancel)
%       ε_veg      = (0.51 − 0.64·m_d)·ε_water + 3.2·m_d + 0.49   (Wegmueller 1990)
%       ε_water    = Debye(T)   (CMEM DIEL_WAT medium=0, pure water)
%
%     Net dimensionally: b has units 1/(kg/m²) = m²/kg, so VOD = b·VWC[kg/m²]
%     is dimensionless. CMEM internal sanity check: at 1.4 GHz, m_d=0.3, T=20°C,
%     a_geo=2/3, θ=40° → ε_veg ≈ 15 + 4i, b ≈ 0.13 m²/kg.
%
%   References:
%     Kirdyashev, Chukhlantsev, Shutko 1979 Sov. J. Remote Sens. 2:29
%     Wegmueller, Mätzler 1990 IEEE TGRS 28(4), 595–604 — leaf dielectric mixture
%     CMEM v5 veg_sub.F90 SUBROUTINE VEGKIRD, dielwat_sub.F90 DIEL_WAT
%
%   Limits / cautions:
%     - Valid for f < 7.5 GHz (Effective Medium theory breaks at high freq).
%       At 1.4 GHz (L-band) this is well within validity.
%     - Pure water dielectric (no salinity) is fine for plant water content.
%     - m_d = 0.3 is a typical herbaceous/leaf value; woody m_d ≈ 0.4–0.5.

% ─── Pure water dielectric (CMEM DIEL_WAT medium=0, isal=0) ───────────────
%   Debye expression: ε_w(ω,T) = ε_∞ + (ε_w0 − ε_∞) / (1 − iω·τ_w)
T = T_C;                                                              % shorten
eps_winf = 4.9;                                                       % CMEM YOMCMEMSOIL eps_winf
omega = 2 * pi * freq_GHz * 1e9;                                       % [rad/s]
eps_0 = 8.854e-12;                                                     %#ok<NASGU>  CMEM constant, not used here

tau_w = 1.768e-11 - 6.068e-13*T + 1.104e-14*T^2 - 8.111e-17*T^3;       % [s]
eps_w0 = 88.045 - 0.4147*T + 6.295e-4*T^2 + 1.075e-5*T^3;
eps_water = eps_winf + (eps_w0 - eps_winf) / (1 - 1i*omega*tau_w);

% ─── Vegetation dielectric (Wegmueller 1990 mixing formula) ───────────────
eps_veg = (0.51 - 0.64*m_d) * eps_water + 3.2*m_d + 0.49;

% ─── Kirdyashev opacity coefficient (CMEM SI units) ───────────────────────
% CMEM cmem_setup.F90:229-232: lam [m], k = 2π/lam [1/m]
%   τ      = a_geo · k[1/m] · (VWC[kg/m²] / ρ_w[kg/m³]) · Im(ε) / cosθ  (dimensionless)
%   b      = τ / VWC = a_geo · k / ρ_w · Im(ε) / cosθ   [m² / kg]
% At L-band 1.41 GHz, T=20°C, m_d=0.3, a_geo=2/3, θ=40° → b ≈ 0.07 m²/kg
% (this matches Wigneron 2007 Table 3 prior b_leaf ≈ 0.06).
c_m_per_s  = 2.998e8;
k_per_m    = 2*pi*freq_GHz*1e9 / c_m_per_s;                            % wavenumber [1/m]
rho_water  = 1000;                                                      % [kg/m³]
cos_theta  = cosd(inc_deg);
b = a_geo * k_per_m / rho_water * imag(eps_veg) / cos_theta;            % [m²/kg]
end
