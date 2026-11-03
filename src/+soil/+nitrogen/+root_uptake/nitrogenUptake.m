function [SinkS, SPUptake, uptakeDiag] = nitrogenUptake( ...
        Srt, NL, DeltZ, C_NH4, C_NO3, SoilState, RLD, daisyParams, nacryParams)
%NITROGENUPTAKE  Root N uptake via coupled DAISY single-root + Nacry membrane (Phase C M1).
%
% Replaces the legacy passive (max(Srt,0)×C) + active (MM-only-NO3) two-branch
% formulation with a unified single mechanism: at each layer, both NH4 and NO3
% radial fluxes are determined by enforcing steady state between the cylindrical
% advection-dispersion equation (DAISY, soil→root surface, Kløffel 2026) and the
% membrane Michaelis-Menten kinetics (Nacry 2013 AMT/NRT2). Bisection finds the
% root-surface concentration C_root that closes the coupling, giving a single
% per-layer uptake I_x [kg/m-root/s] which is then aggregated to layer sinks.
%
% Inputs:
%   Srt          [in]  [1/s, NL×2]      Root water uptake at layer-node pairs (water side).
%                                        Positive = soil→root, negative = HR back-flow.
%   NL           [in]  [-]               Number of layers.
%   DeltZ        [in]  [cm, NL×1]       Layer thickness.
%   C_NH4        [in]  [μg/cm³, NN×1]   NH4 concentration at nodes (= mg/L = g/m³).
%   C_NO3        [in]  [μg/cm³, NN×1]   NO3 concentration at nodes.
%   SoilState    [in]  struct           Required fields:
%                                        .Theta_LL   [m³/m³, NL×2]  current water content
%                                        .Theta_s    [m³/m³, NL×1] saturated water content
%                                        .RHO_bulk   [g/cm³, NL×1] soil bulk density
%                                        .RHOKG      [-, NL×1]     ρ × K_d (= b_NH4 - θ)
%                                        .DIF1111    [cm²/s, 1×NL×nS] free-water diffusion
%                                        .TT         [°C, NN×1]   soil temperature (reserved)
%   RLD          [in]  [m-root/m³-soil, NL×1]  Root length density profile.
%   daisyParams  [in]  struct           DAISY-side parameters:
%                                        .r_root_m            [m]  root radius
%                                        .SRL                 [m-root/g-root]  specific root length
%   nacryParams  [in]  struct           Nacry-side parameters (and solver):
%                                        .Vmax_NH4_umol_g_h, .Vmax_NO3_umol_g_h
%                                        .Km_NH4, .Km_NO3                     [mg/L]
%                                        .Cmin_NH4, .Cmin_NO3                 [mg/L]
%                                        .solver_tol, .solver_max_iter
%                                        Phase C M2-B.2 (2026-05-19) optional:
%                                        .demand_NH4_g_m2_s, .demand_NO3_g_m2_s [g N/(m²·s)]
%                                          If present: switch to demand-driven mode
%                                          (DAISY = demand per root), Nacry as upstream cap.
%                                          Missing/empty: M1 mode (DAISY = Nacry).
%
% Returns:
%   SinkS        [μg/(cm³·s), NL×2×2]  Layer sink terms; (:,:,1)=NH4, (:,:,2)=NO3.
%                                        Both nodes of each layer get equal half-share
%                                        of the depth-integrated uptake.
%   SPUptake     [μg/(cm²·s), 1×2]    Depth-integrated uptake [NH4, NO3].
%   uptakeDiag   struct                Per-species per-layer diagnostics:
%                                        .C_root_NH4 / .C_root_NO3   [μg/cm³ = mg/L]
%                                        .alpha_NH4 / .alpha_NO3     [-]
%                                        .beta                       [-]
%                                        .I_x_NH4 / .I_x_NO3         [kg/m-root/s]
%                                        .solver_iter_NH4/.solver_iter_NO3
%                                        .solver_failed_NH4/.solver_failed_NO3
%
% Method (per layer L):
%   1. Compute q_w per single root from Srt (water side already includes
%      Sperry water-stress reduction, so no double-counting).
%   2. For each species s ∈ {NH4, NO3}:
%      a. b_x = θ      (NO3)         or θ + RHOKG    (NH4)
%      b. D_x* = DIF1111 × Millington-Quirk (θ^(7/3)/θ_s²)
%      c. D_x = D_x* × b_x
%      d. Convert all to SI for DAISY:
%           C_bulk_SI = C_NH4 [μg/cm³] × 1e-3 → [kg/m³]
%           D_x_SI    = D_x   [cm²/s]  × 1e-4 → [m²/s]
%           q_w in [m³/(m-root·s)] already
%      e. solveRootSurface returns C_root [kg/m³] and I_x [kg/m-root/s].
%      f. Convert back to model units:
%           SinkS [μg/(cm³·s)] = I_x [kg/m-root/s] × RLD [m-root/m³] × 1000
%                                                          ↑↑↑↑
%                                                          (kg→μg)×1e9 / (m³→cm³)×1e6 = 1000
%
% Reference:
%   Kløffel 2026 (DAISY 7.1.0 Appendix 7.2) — single-root cylindrical ADE
%   Nacry et al. 2013 Plant and Soil 370:1-29 — root membrane M-M
%   cn_coupling_path_and_literature_audit.md §12.8 — full coupling derivation
%
% Phase C M1 (2026-05-18): full rewrite. Replaces v33.2 passive+active two-branch
% formulation. lActRSU switch removed (active uptake is now intrinsic). N_demand
% removed from signature (M2 will reintroduce as C_min(t) modulation, not separate arg).

% ---- 0. Init ----
SinkS    = zeros(NL, 2, 2);
SPUptake = zeros(1, 2);

uptakeDiag = struct( ...
    'C_root_NH4', zeros(NL, 1), ...
    'C_root_NO3', zeros(NL, 1), ...
    'alpha_NH4',  zeros(NL, 1), ...
    'alpha_NO3',  zeros(NL, 1), ...
    'beta',       zeros(NL, 1), ...
    'I_x_NH4',    zeros(NL, 1), ...
    'I_x_NO3',    zeros(NL, 1), ...
    'solver_iter_NH4',   zeros(NL, 1), ...
    'solver_iter_NO3',   zeros(NL, 1), ...
    'solver_failed_NH4', zeros(NL, 1), ...
    'solver_failed_NO3', zeros(NL, 1));

% ---- Constants ----
r_r = daisyParams.r_root_m;            % [m]
SRL = daisyParams.SRL;                 % [m-root/g-root]
tol      = nacryParams.solver_tol;     % [mg/L]
max_iter = nacryParams.solver_max_iter;

% Unit conversions
%   C [μg/cm³] → [kg/m³]:  × 1e-3
%   D [cm²/s]  → [m²/s]:   × 1e-4
%   I_x [kg/m-root/s] × RLD [m-root/m³] → [kg/m³/s] → [μg/cm³/s]: × 1000
%   tol [mg/L] = [g/m³] = [kg/m³] × 1e-3:  tol_SI = tol × 1e-3
C_to_SI       = 1e-3;
D_to_SI       = 1e-4;
Ix_to_SinkS   = 1e3;   % kg/m³/s → μg/cm³/s
tol_SI        = tol * 1e-3;

% Pre-extract soil state
Theta_LL  = SoilState.Theta_LL;          % [m³/m³, NL×2]
Theta_s   = SoilState.Theta_s(:);        % [m³/m³, NL×1]
RHOKG     = SoilState.RHOKG(:);          % [-, NL×1]  (= ρ × K_d, dimensionless when K_d in cm³/g)
DIF1111   = SoilState.DIF1111;           % [cm²/s, 1×NL×nS]; (:,:,1)=NH4, (:,:,2)=NO3

% Vmax in physiological units (μmol/(g-root·h))
Vmax_NH4 = nacryParams.Vmax_NH4_umol_g_h;
Vmax_NO3 = nacryParams.Vmax_NO3_umol_g_h;
Km_NH4_SI   = nacryParams.Km_NH4   * C_to_SI;     % [kg/m³]
Km_NO3_SI   = nacryParams.Km_NO3   * C_to_SI;
Cmin_NH4_SI = nacryParams.Cmin_NH4 * C_to_SI;
Cmin_NO3_SI = nacryParams.Cmin_NO3 * C_to_SI;

% ---- M2-B.2 demand mode: pre-allocate per-layer per-root demand ----
% demand_NH4/NO3 [g N/(m²·s)] from整株调度 → 按 RLD×DeltZ 加权分到 layer
% → demand_per_root [kg N/(m-root·s)] = demand_layer [g/m³/s] × 1e-3 / RLD [m/m³]
%                                     = demand_total × weight / DeltZ_m / 1000 / RLD
use_demand = isfield(nacryParams, 'demand_NH4_g_m2_s') && ...
             isfield(nacryParams, 'demand_NO3_g_m2_s') && ...
             ~isempty(nacryParams.demand_NH4_g_m2_s) && ...
             ~isempty(nacryParams.demand_NO3_g_m2_s);
if use_demand
    % RLD × DeltZ_m gives root length [m] per m² ground → weight
    DeltZ_m = DeltZ * 1e-2;                       % [m, NL×1]
    root_len_per_layer = RLD(:) .* DeltZ_m(:);    % [m-root/m²-ground, NL×1]
    total_root_len = sum(root_len_per_layer);     % [m-root/m²-ground]
    if total_root_len < 1e-6
        use_demand = false;                        % no roots, skip demand
    end
end

% ---- 1. Loop over layers ----
for ML = 1:NL
    L_rld = RLD(ML);                      % [m-root/m³-soil]

    if L_rld < 1e-6
        continue;                          % no roots → no uptake
    end

    % Layer-averaged θ (mid of top/bottom node)
    theta_mean = 0.5 * (Theta_LL(ML, 1) + Theta_LL(ML, 2));
    if theta_mean < 1e-4
        continue;                          % frozen / dry layer → skip
    end

    % Millington-Quirk tortuosity (matches soluteDispersion.m line 65)
    fw = theta_mean^(7.0/3.0) / Theta_s(ML)^2;

    % q_w from Srt: layer-average uptake rate [1/s] × DeltZ_m → m/s
    % Single-root line flow: q_w = 2π × r_r × v_w
    % v_w = volumetric flow per unit root area = (layer flux) / (RLD × layer thick)
    %     ≈ Srt_layer × DeltZ_m / (2π × r_r × RLD × DeltZ_m)
    %     = Srt_layer / (2π × r_r × RLD)
    % So q_w = 2π × r_r × v_w = Srt_layer / RLD     [m³/m-root/s]
    Srt_layer = 0.5 * (max(Srt(ML, 1), 0) + max(Srt(ML, 2), 0)) * (DeltZ(ML) * 1e-2);
    %                                                              ↑↑↑↑↑↑↑↑↑↑↑↑↑↑
    %                                                              cm → m
    q_w = Srt_layer / L_rld;               % [m³/(m-root·s)]

    % Node-average bulk concentrations [μg/cm³] → [kg/m³]
    C_NH4_bulk_SI = 0.5 * (C_NH4(ML) + C_NH4(ML + 1)) * C_to_SI;
    C_NO3_bulk_SI = 0.5 * (C_NO3(ML) + C_NO3(ML + 1)) * C_to_SI;

    if C_NH4_bulk_SI < 0, C_NH4_bulk_SI = 0; end
    if C_NO3_bulk_SI < 0, C_NO3_bulk_SI = 0; end

    % --- demand per root (M2-B.2 mode) [kg N/(m-root·s)] ---
    % 总 demand [g/m²/s] × 该层根长占比 / 该层总根长 = 单根 demand
    % weight_ML = root_len_per_layer(ML) / total_root_len
    % demand_layer_per_m2 = demand_total × weight_ML  [g N/(m²·s)]
    % demand_per_root = demand_layer_per_m2 / root_len_per_layer × 1e-3 (g→kg)
    %                 = demand_total / total_root_len × 1e-3
    % 即: 均匀按 root length 分摊 (Apple 假设全根等效)
    if use_demand
        demand_NH4_per_root = nacryParams.demand_NH4_g_m2_s / total_root_len * 1e-3;
        demand_NO3_per_root = nacryParams.demand_NO3_g_m2_s / total_root_len * 1e-3;
    else
        demand_NH4_per_root = [];
        demand_NO3_per_root = [];
    end

    % --- Species 1: NH4 ---
    b_NH4 = theta_mean + RHOKG(ML);        % buffer power (DAISY framework)
    D_NH4_star = DIF1111(1, ML, 1) * fw * D_to_SI;       % [m²/s]
    D_NH4_eff  = D_NH4_star * b_NH4;                      % [m²/s]

    [C_root_NH4_SI, I_x_NH4, iter_NH4, status_NH4] = ...
        soil.nitrogen.root_uptake.solveRootSurface( ...
            C_NH4_bulk_SI, q_w, D_NH4_eff, r_r, L_rld, ...
            Vmax_NH4, Km_NH4_SI, Cmin_NH4_SI, SRL, tol_SI, max_iter, ...
            demand_NH4_per_root);

    % M2-B.2 Nacry upstream cap: 在 demand 模式下, 实际吸 = min(DAISY=demand, Nacry@C_root)
    if use_demand && I_x_NH4 > 0
        I_x_nacry = soil.nitrogen.root_uptake.nacryMembraneFlux( ...
                        C_root_NH4_SI, Vmax_NH4, Km_NH4_SI, Cmin_NH4_SI, SRL);
        I_x_NH4 = min(I_x_NH4, I_x_nacry);
    end

    % Convert single-root flux to layer sink: kg/(m-root·s) × m-root/m³ × 1000 = μg/(cm³·s)
    sink_NH4 = max(I_x_NH4, 0) * L_rld * Ix_to_SinkS;
    SinkS(ML, 1, 1) = sink_NH4;
    SinkS(ML, 2, 1) = sink_NH4;

    % --- Species 2: NO3 ---
    b_NO3 = theta_mean;                    % no sorption
    D_NO3_star = DIF1111(1, ML, 2) * fw * D_to_SI;
    D_NO3_eff  = D_NO3_star * b_NO3;

    [C_root_NO3_SI, I_x_NO3, iter_NO3, status_NO3] = ...
        soil.nitrogen.root_uptake.solveRootSurface( ...
            C_NO3_bulk_SI, q_w, D_NO3_eff, r_r, L_rld, ...
            Vmax_NO3, Km_NO3_SI, Cmin_NO3_SI, SRL, tol_SI, max_iter, ...
            demand_NO3_per_root);

    if use_demand && I_x_NO3 > 0
        I_x_nacry = soil.nitrogen.root_uptake.nacryMembraneFlux( ...
                        C_root_NO3_SI, Vmax_NO3, Km_NO3_SI, Cmin_NO3_SI, SRL);
        I_x_NO3 = min(I_x_NO3, I_x_nacry);
    end

    sink_NO3 = max(I_x_NO3, 0) * L_rld * Ix_to_SinkS;
    SinkS(ML, 1, 2) = sink_NO3;
    SinkS(ML, 2, 2) = sink_NO3;

    % --- Diagnostics ---
    uptakeDiag.C_root_NH4(ML) = C_root_NH4_SI / C_to_SI;  % back to μg/cm³ = mg/L
    uptakeDiag.C_root_NO3(ML) = C_root_NO3_SI / C_to_SI;
    if D_NH4_eff > 0
        uptakeDiag.alpha_NH4(ML) = q_w / (2.0 * pi * D_NH4_eff);
    end
    if D_NO3_eff > 0
        uptakeDiag.alpha_NO3(ML) = q_w / (2.0 * pi * D_NO3_eff);
    end
    uptakeDiag.beta(ML) = 1.0 / (r_r * sqrt(pi * L_rld));
    uptakeDiag.I_x_NH4(ML) = I_x_NH4;
    uptakeDiag.I_x_NO3(ML) = I_x_NO3;
    uptakeDiag.solver_iter_NH4(ML) = iter_NH4;
    uptakeDiag.solver_iter_NO3(ML) = iter_NO3;
    uptakeDiag.solver_failed_NH4(ML) = double(strcmp(status_NH4, 'failed'));
    uptakeDiag.solver_failed_NO3(ML) = double(strcmp(status_NO3, 'failed'));

    % Depth-integrated uptake
    SPUptake(1) = SPUptake(1) + sink_NH4 * DeltZ(ML);
    SPUptake(2) = SPUptake(2) + sink_NO3 * DeltZ(ML);
end

end
