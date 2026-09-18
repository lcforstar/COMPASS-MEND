function I_x = daisySingleRootFlux(C_bulk, C_root, q_w, D_x, r_r, L_rld)
%DAISYSINGLEROOTFLUX  Steady-state single-root N flux (DAISY single-root concept).
%
%   Solves the cylindrical advection-dispersion equation for solute X around
%   a single root in steady state, returning the radial flux per unit root length.
%   Three α-case branches (α = q_w/(2π D_x) — radial Péclet number) per
%   Kløffel 2026 (DAISY 7.1.0 Appendix 7.2 Eq. 25).
%
%   Returns:
%     I_x       [kg N / (m-root · s)]  Single-root radial N flux.
%                                       Positive = uptake from soil to root surface.
%
%   Inputs:
%     C_bulk    [kg/m³]   Far-field (bulk soil-solution) N concentration.
%     C_root    [kg/m³]   N concentration at the root surface (membrane boundary).
%                          Determined by coupled DAISY+Nacry surface solution.
%     q_w       [m³/(m-root·s)]  Single-root line water flow rate = 2π r_r v_w.
%                                  Positive towards root, negative for hydraulic
%                                  redistribution.
%     D_x       [m²/s]    Effective dispersion coefficient (= D_x* × b_x),
%                          species-specific via buffer power b_x.
%     r_r       [m]       Root radius.
%     L_rld     [m-root/m³-soil]  Root length density (RLD).
%
%   Method:
%     α  = q_w / (2π D_x)            radial Péclet number
%     β  = 1 / (r_r × √(π × L_rld))  half-distance / r_r ratio (single-root cylinder)
%
%     Case 1 (|α| < 0.05, pure diffusion limit):
%        I_x = (C_bulk - C_root) × 2π × D_x × [β² ln(β)/(β² - 1) - 1/2]^(-1)
%
%     Case 3 (|α - 2| < 0.05, removable singularity at α = 2):
%        I_x = q_w × [C_bulk × (β² - 1) - C_root × ln(β²)] / [(β² - 1) - ln(β²)]
%
%     Case 2 (general):
%        I_x = q_w × [C_bulk × (β²-1)(1-α/2) - C_root × (β^(2-α) - 1)] /
%                    [(β²-1)(1-α/2) - (β^(2-α) - 1)]
%
%   ε-guards:
%     - β floor at 1.05 (avoids ln(β) → 0 when RLD enormous)
%     - L_rld floor 1e-6 (avoids β → ∞ when zero root density)
%
%   Reference:
%     Kløffel, T. (2026). Derivation of the flux of substance X towards roots
%     using the single-root concept. DAISY 7.1.0 documentation, Appendix 7.2.
%
%   Phase C M1 (2026-05-18): new module, called by solveRootSurface +
%   nitrogenUptake. Companion to nacryMembraneFlux.

% ---- ε-guards on geometry ----
if L_rld < 1e-6
    I_x = 0;
    return;
end

beta = 1.0 / (r_r * sqrt(pi * L_rld));
if beta < 1.05
    beta = 1.05;     % collapse beta to safe floor (extremely dense root)
end

% ---- ε-guards on dispersion ----
if D_x <= 0
    I_x = 0;
    return;
end

alpha = q_w / (2.0 * pi * D_x);

% ---- Three α-case branches ----
EPS_ALPHA = 0.05;

beta2 = beta * beta;

if abs(alpha) < EPS_ALPHA
    % Case 1: pure diffusion limit
    denom = beta2 * log(beta) / (beta2 - 1.0) - 0.5;
    I_x = (C_bulk - C_root) * 2.0 * pi * D_x / denom;

elseif abs(alpha - 2.0) < EPS_ALPHA
    % Case 3: removable singularity at α = 2
    ln_b2 = log(beta2);
    num = C_bulk * (beta2 - 1.0) - C_root * ln_b2;
    den = (beta2 - 1.0) - ln_b2;
    I_x = q_w * num / den;

else
    % Case 2: general
    one_minus_half_alpha = 1.0 - 0.5 * alpha;
    pow_term = beta^(2.0 - alpha) - 1.0;
    num = C_bulk * (beta2 - 1.0) * one_minus_half_alpha - C_root * pow_term;
    den = (beta2 - 1.0) * one_minus_half_alpha - pow_term;
    I_x = q_w * num / den;
end

end
