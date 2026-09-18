function hydraulic = solveEbalHydraulicsLayered( ...
    plantHydro, soil_state, soil_conductivity, KT, H, Tp_sun_layer_cm_s, Tp_sha_layer_cm_s, layer_weight, mode, soil_params)
%SOLVEEBALHYDRAULICSLAYERED Layer-resolved hydraulic stress for vegetation.energy.ebal.
%
% Phase 3: LAI-proportional per-layer psi_leaf and wt.
%
%   The aggregate two-leaf solve (unchanged) determines psi_stem and actual
%   total flows.  Each canopy layer then receives a fraction of the stem
%   conductance proportional to its LAI.  The resulting per-layer leaf water
%   potential and hydraulic stress factor follow from the same vulnerability-
%   curve lookup used in the aggregate solve.
%
%   Inputs:
%     soil_state         target soil slice (.h, .theta_liquid)
%     soil_conductivity  NL x 2 conductivity for root-soil path
%     soil_params        optional retention (.theta_s)
%
%   Physics:
%     K(psi_stem) - K(psi_leaf_il) = flow_il / lai_frac_il
%     flow_il   = actual_sun/sha * (demand_il / demand_total)   [m s-1]
%     lai_frac_il = LAI_layer_il / LAI_total                    [-]
%
%   For uniform LAI (lai_frac = 1/nl for all layers) this reduces exactly
%   to the aggregate result when demand per layer is also uniform.
%
% Editor: Jizixing Zhu; Date: 2026-04-20.
% Change: 2026-07-16 Jizixing Zhu — accept soil_state / soil_params target names (no legacy bridge).

    Tp_sun_layer_cm_s = max(0, Tp_sun_layer_cm_s(:));
    Tp_sha_layer_cm_s = max(0, Tp_sha_layer_cm_s(:));
    nl = numel(Tp_sun_layer_cm_s);

    if nargin < 8 || isempty(layer_weight)
        layer_weight = ones(nl, 1) / nl;
    else
        layer_weight = layer_weight(:);
    end
    use_layer_stress = nargin >= 9 && isequal(mode, 'LAI_proportional');

    if numel(Tp_sha_layer_cm_s) ~= nl || numel(layer_weight) ~= nl
        error('PlantHydraulic:BadLayerDemandSize', ...
            'Layer demand and layer_weight must have the same number of canopy layers.');
    end

    % --- Aggregate solve (root-soil-stem, unchanged) ---
    Tp_sun_total = sum(Tp_sun_layer_cm_s);
    Tp_sha_total = sum(Tp_sha_layer_cm_s);
    if nargin >= 10
        hydraulic = vegetation.hydraulics.solveEbalHydraulics( ...
            plantHydro, soil_state, soil_conductivity, KT, H, Tp_sun_total, Tp_sha_total, soil_params);
    else
        hydraulic = vegetation.hydraulics.solveEbalHydraulics( ...
            plantHydro, soil_state, soil_conductivity, KT, H, Tp_sun_total, Tp_sha_total);
    end

    % Demand and actual flow per layer
    hydraulic.Tp_sun_layer_demand = Tp_sun_layer_cm_s;
    hydraulic.Tp_sha_layer_demand = Tp_sha_layer_cm_s;
    hydraulic.Tp_layer_demand = Tp_sun_layer_cm_s + Tp_sha_layer_cm_s;

    hydraulic.Tp_sun_layer = zeros(nl, 1);
    hydraulic.Tp_sha_layer = zeros(nl, 1);
    if Tp_sun_total > 0
        hydraulic.Tp_sun_layer = hydraulic.Tp_sun * Tp_sun_layer_cm_s / Tp_sun_total;
    end
    if Tp_sha_total > 0
        hydraulic.Tp_sha_layer = hydraulic.Tp_sha * Tp_sha_layer_cm_s / Tp_sha_total;
    end
    hydraulic.Tp_layer = hydraulic.Tp_sun_layer + hydraulic.Tp_sha_layer;

    % Default: aggregate scalar (fallback when hydraulics not active)
    hydraulic.psil_sun_layer = hydraulic.psil_sun * ones(nl, 1);
    hydraulic.psil_sha_layer = hydraulic.psil_sha * ones(nl, 1);
    hydraulic.psir_layer     = hydraulic.psir     * ones(nl, 1);
    hydraulic.psis_layer     = hydraulic.psis     * ones(nl, 1);
    hydraulic.wt_sun_layer   = hydraulic.wt_sun   * ones(nl, 1);
    hydraulic.wt_sha_layer   = hydraulic.wt_sha   * ones(nl, 1);
    hydraulic.layer_hydraulic_mode = 'shared_path_scalar_repeated';

    % --- Phase 3: LAI-proportional per-layer ψ_leaf ---
    if ~use_layer_stress || Tp_sun_total + Tp_sha_total < 1e-12
        return;
    end

    % Stem water potential from aggregate result
    gravity_MPa = 1000 * 9.80665 * H / 1e6;   % [MPa]
    pstem = hydraulic.psir - gravity_MPa;

    % Build K(ψ) lookup table (same as calSPAC)
    psi_table = plantHydro.PSIL_RANGE(:);
    K_table   = plantHydro.KS(:);
    [psi_table, order] = sort(psi_table);
    K_table   = K_table(order);
    psi_min   = min(psi_table);
    psi_max   = max(psi_table);
    [K_unique, ia] = unique(K_table, 'stable');
    psi_for_K = psi_table(ia);
    K_min = min(K_unique);
    K_max = max(K_unique);

    pstem_cl = min(max(pstem, psi_min), psi_max);
    K_pstem  = interp1(psi_table, K_table, pstem_cl, 'linear');

    % LAI fraction per layer (normalised, sums to 1)
    lai_frac = layer_weight(:);
    lai_frac = max(lai_frac, 1e-6);   % guard against zero-LAI layers

    % Per-layer leaf water potential and wt
    actual_sun_m_s = hydraulic.Tp_sun / 100;   % cm/s → m/s
    actual_sha_m_s = hydraulic.Tp_sha / 100;

    for il = 1:nl
        % --- sunlit ---
        if Tp_sun_total > 0 && actual_sun_m_s > 0
            flow_sun_il = actual_sun_m_s * Tp_sun_layer_cm_s(il) / Tp_sun_total;
            K_tgt = K_pstem - flow_sun_il / lai_frac(il);
            K_tgt = min(max(K_tgt, K_min), K_max);
            hydraulic.psil_sun_layer(il) = interp1(K_unique, psi_for_K, K_tgt, 'linear');
            wt = vegetation.hydraulics.VulnerabilityCurve( ...
                    hydraulic.psil_sun_layer(il), plantHydro.gpmax, plantHydro.p50, plantHydro.aa) ...
                 / plantHydro.gpmax;
            hydraulic.wt_sun_layer(il) = max(0, min(1, wt));
        end

        % --- shaded ---
        if Tp_sha_total > 0 && actual_sha_m_s > 0
            flow_sha_il = actual_sha_m_s * Tp_sha_layer_cm_s(il) / Tp_sha_total;
            K_tgt = K_pstem - flow_sha_il / lai_frac(il);
            K_tgt = min(max(K_tgt, K_min), K_max);
            hydraulic.psil_sha_layer(il) = interp1(K_unique, psi_for_K, K_tgt, 'linear');
            wt = vegetation.hydraulics.VulnerabilityCurve( ...
                    hydraulic.psil_sha_layer(il), plantHydro.gpmax, plantHydro.p50, plantHydro.aa) ...
                 / plantHydro.gpmax;
            hydraulic.wt_sha_layer(il) = max(0, min(1, wt));
        end
    end

    hydraulic.psil_sun_layer = max(psi_min, hydraulic.psil_sun_layer);
    hydraulic.psil_sha_layer = max(psi_min, hydraulic.psil_sha_layer);
    hydraulic.layer_hydraulic_mode = 'LAI_proportional_v1';
end
