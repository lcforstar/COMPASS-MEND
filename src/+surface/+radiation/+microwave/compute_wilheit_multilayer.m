function [R_pol, Tb_eff, f_layers] = compute_wilheit_multilayer(freq_GHz, theta_deg, ...
                                                       n_prof, dz_cm, T_prof_K, pol)
%WILHEIT_MULTILAYER  Stratified-dielectric brightness temperature (Wilheit 1978).
%
%   Returns:
%     R_pol    [-]      Power reflectivity of the multilayer at the top interface
%     Tb_eff   [K]      Effective brightness temperature emerging upward from
%                       the multilayer = Σ f_j · t_j (layers 2..N-1)
%     f_layers [-, N×1] Per-layer thermal weighting function. f(1) = R_pol (the
%                       weight for the incident sky radiation), f(2..N-1) are
%                       layer absorption fractions, f(N) ≈ 0 (semi-infinite tail).
%
%   Inputs:
%     freq_GHz   [GHz]      Microwave frequency
%     theta_deg  [deg]      Incidence angle in the first (air) layer
%     n_prof     [-, N×1]   Complex refractive index per layer (top → bottom).
%                            Layer 1 = air = (1, 0); layer N = semi-infinite bottom.
%     dz_cm      [cm, N×1]  Layer thicknesses. dz(1) and dz(N) are unused
%                            (semi-infinite top/bottom).
%     T_prof_K   [K, N×1]   Layer physical temperatures.
%     pol        char       'H' or 'V'
%
%   Method:
%     Direct MATLAB port of Wilheit 1978 (NASA TM-X-71051) Appendix I
%     FORTRAN subroutines HPLD / VPLD. Algorithm:
%
%       1. Forward pass: build propagator P_j = exp(i · 2π/λ · n_j · cos θ_j · δ_j)
%          cumulatively from layer 2 downward. Truncate at layer NMAX where
%          |P| < 1e-4 (subsequent layers are EM-opaque, set P = 0/(i·1e50)).
%
%       2. Boundary at layer NMAX: set E+ = 1, E- = 0 (no upcoming wave).
%
%       3. Backward sweep (J = NMAX-1 down to 1): apply Fresnel boundary
%          conditions (continuity of E_t and H_t) at each j↔j+1 interface to
%          solve for E+_j, E-_j. For H-pol:
%              SA = 2 n_j cos θ_j / (n_j cos θ_j + n_{j+1} cos θ_{j+1})
%              SR = (n_j cos θ_j − n_{j+1} cos θ_{j+1}) / [(n_j cos θ_j + n_{j+1} cos θ_{j+1}) P_j]
%              E+_j = E+_{j+1}/SA  +  SR · E-_{j+1}/SA
%              E-_j = E-_{j+1} + (E+_{j+1} − E+_j) · P_j
%          For V-pol (Wilheit Eq Section "Vertical polarization"):
%              SA = n_j cos θ_{j+1} + n_{j+1} cos θ_j
%              SR = (n_j cos θ_{j+1} − n_{j+1} cos θ_j) / P_j
%              E+_j = SA · E+_{j+1}/[denom × P_j]  +  SR · E-_{j+1}/[denom × P_j]
%       4. Renormalize E_j /= E+_1 so the incident wave has unit amplitude.
%
%       5. Per-layer weighting function f_j = (Φ_{j-1} − Φ_j) / Φ_0, where
%          Φ_j is the Poynting flux at interface j (Wilheit pp.6-7 "Radiative
%          Transfer" section). The Re(n cos θ) terms give the attenuated
%          forward/backward waves, the Im(n cos θ) cross-term gives the
%          standing-wave correction.
%
%       6. R_pol = |E-_1|^2 (reflectivity at top interface).
%
%       7. Tb_eff = Σ_{j=2..NMAX-1} f_j · T_j  (downwelling sky T_sky implicit
%          in the (1 − R) weighting; here we return Tb upwelling from the soil
%          only — caller adds (1 − R)·T_canopy_or_sky via τ-ω synthesis).
%
%   References:
%     Wilheit, T. T. 1978, IEEE Trans. Geosci. Electron. GE-16(2), 138-143
%       (originally NASA-TM-X-71051, March 1975) — Appendix I lines 1-50.
%     Tsang, Kong, Shin 1985 "Theory of Microwave Remote Sensing" Ch.3
%       (alternative derivation of the same recursion).
%
%   Notes:
%     - First layer must be air with n_prof(1) = complex(1, 0). dz(1) is
%       ignored (semi-infinite).
%     - Last layer is semi-infinite homogeneous half-space; dz(N) is ignored.
%     - Function does NOT include rough-surface or vegetation effects;
%       caller (microwave_forward) must compose with PassiveRoughness or I2EM
%       on top of (1 − R_pol) as the smooth-surface emissivity proxy.

if nargin < 6, pol = 'H'; end
n_prof  = n_prof(:);
dz_cm   = dz_cm(:);
T_prof_K = T_prof_K(:);
N = numel(n_prof);
assert(numel(dz_cm) == N && numel(T_prof_K) == N, 'compute_wilheit_multilayer: profile length mismatch');

c_mps     = 2.998e8;
lambda_cm = (c_mps / (freq_GHz * 1e9)) * 100;
theta_rad = theta_deg * pi/180;
s = sin(theta_rad);
ctheta_top = cos(theta_rad);

% ─── 1. Forward pass: cumulative propagator P_j ────────────────────────────
P = complex(zeros(N, 1));
P(1) = complex(1, 0);                                  % air layer
NMAX = N;
for j = 2:(N-1)
    ss = n_prof(1) * s / n_prof(j);                    % Snell sin θ_j
    sc = sqrt(complex(1, 0) - ss*ss);                  % cos θ_j
    arg = dz_cm(j) * 2 * pi / lambda_cm;
    s_arg = 2 * arg * n_prof(j) * sc * complex(0, 1);
    P(j) = exp(s_arg) * P(j-1);
    if abs(P(j)) < 1e-4
        NMAX = j;
        break;
    end
end
if NMAX > N-1, NMAX = N-1; end

% Saturate deeper-than-NMAX layers as opaque (Wilheit FORTRAN: SP(J) = i/1e50)
for j = (NMAX+1):N
    P(j) = complex(0, 1) / 1e50;
end

% ─── 2 & 3. Backward sweep: Fresnel boundary conditions ───────────────────
SEP = complex(zeros(N, 1));     % E+ per layer
SEM = complex(zeros(N, 1));     % E-
SEP(NMAX) = complex(1, 0);
SEM(NMAX) = complex(0, 0);

for jj = 2:NMAX
    j = NMAX - jj + 1;
    ssj   = n_prof(1) * s / n_prof(j);
    scj   = sqrt(complex(1, 0) - ssj*ssj);
    ssj1  = n_prof(1) * s / n_prof(j+1);
    scj1  = sqrt(complex(1, 0) - ssj1*ssj1);

    % Interface Fresnel coefficients (Tsang 1985 §3.2)
    if upper(pol) == 'H'
        % H-pol: n_j·cos θ_j ≡ Y_j (admittance)
        a_j   = n_prof(j) * scj;
        a_j1  = n_prof(j+1) * scj1;
    else  % V-pol: n_{j+1}·cos θ_j for the down-side, n_j·cos θ_{j+1} for upper
        % Wilheit p.3 V boundary conditions give Fresnel r_v in this form
        a_j   = n_prof(j+1) * scj;
        a_j1  = n_prof(j)   * scj1;
    end
    t_F = 2 * a_j / (a_j + a_j1);
    r_F = (a_j - a_j1) / (a_j + a_j1);

    SA = t_F;
    SR = r_F / P(j);

    SEP(j) = SEP(j+1) / SA + SR * SEM(j+1) / SA;
    SEM(j) = SEM(j+1) + (SEP(j+1) - SEP(j)) * P(j);
end

% ─── 4. Renormalize so E+_1 = 1 ────────────────────────────────────────────
SX = SEP(1);
if abs(SX) < 1e-12
    R_pol = 0.5;  Tb_eff = mean(T_prof_K(2:end-1));
    f_layers = zeros(N, 1);
    return;
end
for j = 1:NMAX
    SEP(j) = SEP(j) / SX;
    SEM(j) = SEM(j) / SX;
end

% ─── 5. Per-layer weighting function f_j ───────────────────────────────────
% Wilheit pp.6-7: f_j = (Φ_{j-1} − Φ_j) / Φ_0
% with Φ_j ∝ Re(n cos θ){|E+|²(|P_{j-1}|⁻² − |P_j|⁻²) − ...} − Im cross-term.
% Implementation follows FORTRAN HPLD lines 50-60.
f_layers = zeros(N, 1);
LL = NMAX - 1;
for jj = 1:LL
    j  = NMAX - jj + 1;
    Rj = abs(P(j))^2;
    Sj = abs(P(j-1))^2;
    F2 = (Sj - Rj) * abs(SEP(j))^2 + (1/Rj - 1/Sj) * abs(SEM(j))^2;
    ssj = n_prof(1) * s / n_prof(j);
    scj = sqrt(complex(1, 0) - ssj*ssj);
    dP  = F2 * real(n_prof(j) * scj) / ctheta_top;
    SXP = SEP(j) * conj(SEM(j));
    cross = 2 * imag(n_prof(j) * scj) / ctheta_top * ...
            (imag(SXP * P(j-1) / abs(P(j-1))) - imag(SXP * P(j) / abs(P(j))));
    f_layers(j) = dP - cross;
end

% Top reflection coefficient
R_pol = abs(SEM(1))^2;
f_layers(1) = R_pol;        % bookkeeping (sky-reflection weight)

% ─── 6 & 7. Effective brightness temperature from soil layers ──────────────
% (Sum f_j · T_j for j = 2..NMAX; f(1) is sky reflection slot, NOT emission)
mask = (2:NMAX);
fsum = sum(f_layers(mask));
if fsum < 1e-9
    Tb_eff = T_prof_K(2);   % degenerate case
else
    Tb_eff = sum(f_layers(mask) .* T_prof_K(mask)) / fsum * (1 - R_pol);
end

% Clip numerical noise
f_layers = real(f_layers);
R_pol = max(0, min(real(R_pol), 1));
Tb_eff = real(Tb_eff);
end
