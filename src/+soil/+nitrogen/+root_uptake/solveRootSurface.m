function [C_root, I_flux, iter, status] = solveRootSurface( ...
    C_bulk, q_w, D_x, r_r, L_rld, Vmax_umol_g_h, Km, C_min, SRL, tol, max_iter, ...
    demand_per_root)
%SOLVEROOTSURFACE  Bisection solver for steady-state root-surface N concentration.
%
%   Two modes (selected by demand_per_root):
%     M1 (coupled DAISY+Nacry, demand_per_root = [] or missing):
%         I_DAISY(C_root) == I_Nacry(C_root)
%         Nacry permanently saturates root membrane (always at max).
%     M2-B.2 (DAISY = demand, demand_per_root scalar ≥ 0):
%         I_DAISY(C_root) == demand_per_root
%         Plant-side N_labile feedback throttles uptake via demand;
%         Nacry kept as upstream cap applied by caller (not in this solver).
%
%   Returns:
%     C_root    [kg/m³]               Root-surface N concentration.
%     I_flux    [kg N/(m-root·s)]     DAISY flux at C_root (=Nacry in M1, =demand in M2-B.2).
%     iter      [-]                   Bisection iterations used.
%     status    string                'converged' | 'cbulk_below_cmin' |
%                                     'membrane_limited' | 'soil_limited' |
%                                     'demand_saturated' | 'failed'.
%
%   Inputs:
%     C_bulk          [kg/m³]              Far-field N concentration.
%     q_w             [m³/(m-root·s)]      Single-root water line flow.
%     D_x             [m²/s]               Effective dispersion (D_x* × b_x).
%     r_r             [m]                  Root radius.
%     L_rld           [m-root/m³-soil]     Root length density.
%     Vmax_umol_g_h   [μmol N/(g-root·h)]  Membrane max uptake rate (M1 only; ignored if demand).
%     Km              [kg/m³]              Membrane half-saturation (M1 only).
%     C_min           [kg/m³]              Membrane lower bound (used in both modes as floor).
%     SRL             [m-root/g-root]      Specific root length (M1 only).
%     tol             [kg/m³]              Bisection termination tolerance on C_root.
%     max_iter        [-]                  Maximum bisection iterations.
%     demand_per_root [kg N/(m-root·s)]    [optional] target uptake per unit root length.
%                                          [] or missing → M1 mode (default).
%                                          ≥ 0 scalar     → M2-B.2 mode.
%
%   Method (M1):
%     g(C_root) = I_DAISY(C_root) − I_Nacry(C_root) is monotonically decreasing
%     on [C_min, C_bulk]; bisect to find zero.
%
%   Method (M2-B.2):
%     g(C_root) = I_DAISY(C_root) − demand
%     I_DAISY is monotonically decreasing in C_root (浓度差变小→流变小).
%     - At C_root = C_min:  I_DAISY at max (土壤极限供应)
%     - At C_root = C_bulk: I_DAISY = 0
%     If demand > I_DAISY(C_min): demand_saturated, return C_root=C_min, I=I_DAISY(C_min).
%     If demand ≤ 0:                              return C_root=C_bulk, I=0.
%     Else bisect on [C_min, C_bulk].
%
%   Reference:
%     Coupled DAISY + Nacry framework — see cn_coupling_path_and_literature_audit
%     §12.8.1 (M1) and §12.19 (M2-B.2 demand-driven extension, 2026-05-19).
%
%   Phase C M1 (2026-05-18): new module.
%   Phase C M2-B.2 (2026-05-19): added optional demand_per_root mode switch.

% ---- Mode dispatch ----
if nargin < 12 || isempty(demand_per_root)
    use_demand = false;
else
    use_demand = true;
end

% ---- Boundary: bulk below membrane lower bound ----
if C_bulk <= C_min
    C_root = C_bulk;
    I_flux = 0;
    iter = 0;
    status = 'cbulk_below_cmin';
    return;
end

% ---- M2-B.2 demand mode: no Nacry inside, just I_DAISY = demand ----
if use_demand
    if demand_per_root <= 0
        % Plant doesn't want anything → 不抽, C_root 拉到 bulk
        C_root = C_bulk;
        I_flux = 0;
        iter = 0;
        status = 'demand_zero';
        return;
    end

    I_max_at_cmin = soil.nitrogen.root_uptake.daisySingleRootFlux( ...
                        C_bulk, C_min, q_w, D_x, r_r, L_rld);   % [kg N/(m-root·s)]
    if demand_per_root >= I_max_at_cmin
        % Demand 超过土壤极限供应 → 全力供应, 实际 = 土壤极限
        C_root = C_min;
        I_flux = I_max_at_cmin;
        iter = 0;
        status = 'demand_saturated';
        return;
    end

    % Bisection on g(C_root) = I_DAISY(C_root) − demand, 单调递减
    lo = C_min;
    hi = C_bulk;
    for iter = 1:max_iter
        mid = 0.5 * (lo + hi);
        I_mid = soil.nitrogen.root_uptake.daisySingleRootFlux( ...
                    C_bulk, mid, q_w, D_x, r_r, L_rld);
        % [F3 fix 2026-06-03] Flux convergence must use a flux-scaled (relative)
        % tolerance. The previous check compared a flux difference [kg/(m·s)]
        % (~1e-12) against `tol`, which is a C_root CONCENTRATION tolerance
        % [kg/m³] (~1e-6) — a unit mismatch that passed on iteration 1 and
        % returned I_DAISY at the C_root midpoint (~half soil supply), ignoring
        % demand. That let root N uptake run ~10-50x the plant demand. Use a
        % relative flux tolerance for the flux criterion; keep the C_root
        % interval check (tol is in C_root units, correct) as the backstop.
        if abs(I_mid - demand_per_root) <= 1e-3 * demand_per_root || (hi - lo) * 0.5 < tol
            C_root = mid;
            I_flux = I_mid;
            status = 'converged';
            return;
        end
        if I_mid > demand_per_root
            lo = mid;     % DAISY 太大 → C_root 调高
        else
            hi = mid;     % DAISY 太小 → C_root 调低
        end
    end

    C_root = 0.5 * (lo + hi);
    I_flux = soil.nitrogen.root_uptake.daisySingleRootFlux( ...
                C_bulk, C_root, q_w, D_x, r_r, L_rld);
    status = 'failed';
    return;
end

% ---- M1 mode: coupled DAISY = Nacry (original) ----
lo = C_min;
hi = C_bulk;

g_lo = soil.nitrogen.root_uptake.daisySingleRootFlux( ...
            C_bulk, lo, q_w, D_x, r_r, L_rld) ...
       - soil.nitrogen.root_uptake.nacryMembraneFlux( ...
            lo, Vmax_umol_g_h, Km, C_min, SRL);

g_hi = soil.nitrogen.root_uptake.daisySingleRootFlux( ...
            C_bulk, hi, q_w, D_x, r_r, L_rld) ...
       - soil.nitrogen.root_uptake.nacryMembraneFlux( ...
            hi, Vmax_umol_g_h, Km, C_min, SRL);

% Endpoints same sign → two distinct physical regimes:
%   (a) g_lo > 0 ∧ g_hi > 0:
%       即使 C_root 推到 C_bulk, DAISY 仍 > Nacry → 根膜成为瓶颈
%       物理: 对流送来的 N 比根膜能吸的多, 真实通量 = Nacry 在 C_root=C_bulk 处的值,
%       多余 N 累积在根周围 (这里抽象为 C_root=C_bulk).
%   (b) g_lo < 0 ∧ g_hi < 0:
%       即使 C_root = C_min, Nacry 仍 > DAISY → DAISY 端供应不足 (土壤瓶颈)
%       物理: 根膜想吸但土里供应跟不上, 真实通量 = DAISY 在 C_root=C_min 处的值.
% 两种 case 都不是数值"failed", 是物理上真实的瓶颈状态.
if g_lo * g_hi > 0
    if g_lo > 0
        % Membrane bottleneck (a)
        C_root = C_bulk;
        I_flux = soil.nitrogen.root_uptake.nacryMembraneFlux( ...
                    C_root, Vmax_umol_g_h, Km, C_min, SRL);
        status = 'membrane_limited';
    else
        % Soil bottleneck (b)
        C_root = C_min;
        I_flux = soil.nitrogen.root_uptake.daisySingleRootFlux( ...
                    C_bulk, C_root, q_w, D_x, r_r, L_rld);
        status = 'soil_limited';
    end
    iter = 0;
    return;
end

% ---- Bisection loop ----
for iter = 1:max_iter
    mid = 0.5 * (lo + hi);
    g_mid = soil.nitrogen.root_uptake.daisySingleRootFlux( ...
                C_bulk, mid, q_w, D_x, r_r, L_rld) ...
           - soil.nitrogen.root_uptake.nacryMembraneFlux( ...
                mid, Vmax_umol_g_h, Km, C_min, SRL);

    if abs(g_mid) < tol || (hi - lo) * 0.5 < tol
        C_root = mid;
        I_flux = soil.nitrogen.root_uptake.daisySingleRootFlux( ...
                    C_bulk, C_root, q_w, D_x, r_r, L_rld);
        status = 'converged';
        return;
    end

    if g_lo * g_mid < 0
        hi = mid;
    else
        lo = mid;
        g_lo = g_mid;
    end
end

% Didn't converge — return best guess + flag
C_root = 0.5 * (lo + hi);
I_flux = soil.nitrogen.root_uptake.daisySingleRootFlux( ...
            C_bulk, C_root, q_w, D_x, r_r, L_rld);
status = 'failed';

end
