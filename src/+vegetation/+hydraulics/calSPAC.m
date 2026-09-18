function [Tp_sun, Tp_sha, psil_sun, psil_sha, psir, psis, gsr] = ...
    calSPAC(plantHydro, Khh, hh, Tp_sun, Tp_sha, KT, H)
%CALSPAC Couple soil-root-stem-leaf hydraulics for sunlit and shaded leaves.
%   Units:
%     Khh/KL_h : NL x 2 soil hydraulic conductivity, cm/s.
%     hh       : NN x 1 soil pressure head, cmH2O.
%     Tp       : transpiration demand, cm/s; converted to m/s internally.
%     H        : canopy height, m.
%     psis/psir/psil : MPa.
%
%   The sunlit and shaded leaves share one upstream soil-root-stem supply.
%   Their ebal transpiration inputs are treated as demand. The actual total
%   flow is solved by a scalar bisection so that root-stem potential drop and
%   maximum leaf-side supply are mutually consistent.

    rootContact = [];
    if isfield(plantHydro, 'rootContact')
        rootContact = plantHydro.rootContact;
    end

    [psis, gsr, sumgsr] = vegetation.hydraulics.Plantrootconductivity( ...
        plantHydro.RL, plantHydro.Lmax, plantHydro.NN, plantHydro.DeltZ, ...
        Khh, plantHydro.RAI, hh, KT, rootContact, plantHydro);

    if sumgsr <= 0
        % Zero soil-root conductance: stop actual transpiration and use the
        % configured lower leaf-water-potential bound.
        Tp_sun(KT) = 0;
        Tp_sha(KT) = 0;
        psil_sun(KT) = plantHydro.psil_min;
        psil_sha(KT) = plantHydro.psil_min;
        psir(KT)    = psis(KT);
        psis(KT)    = psis(KT);
        return;
    end

    if ~isreal(Tp_sun(KT)) || ~isreal(Tp_sha(KT))
        error('PlantHydraulic:ComplexTranspiration', ...
            'Transpiration demand is complex before SPAC. Check ebal latent heat and unit conversion.');
    end

    demand_sun_m_s = max(0, Tp_sun(KT)) / 100;       % cm/s -> m/s
    demand_sha_m_s = max(0, Tp_sha(KT)) / 100;       % cm/s -> m/s
    demand_total_m_s = demand_sun_m_s + demand_sha_m_s;

    gravity_MPa = 1000 * 9.80665 * H / 1e6;          % rho_w*g*H: Pa -> MPa

    psi_table = plantHydro.PSIL_RANGE(:);
    K_table = plantHydro.KS(:);
    validateLookupTable(psi_table, K_table);

    [psi_table, order] = sort(psi_table);
    K_table = K_table(order);
    [K_unique, ia] = unique(K_table, 'stable');
    psi_for_K = psi_table(ia);
    if numel(K_unique) < 2
        error('PlantHydraulic:BadLookupTable', ...
            'plantHydro.KS must contain at least two unique conductance values.');
    end

    psi_min = min(psi_table);
    psi_max = max(psi_table);
    K_min = min(K_unique);
    K_max = max(K_unique);

    psir = zeros(size(Tp_sun));
    psil_sun = zeros(size(Tp_sun));
    psil_sha = zeros(size(Tp_sha));

    if demand_total_m_s <= 0
        actual_total_m_s = 0;
        actual_sun_m_s = 0;
        actual_sha_m_s = 0;
    else
        actual_total_m_s = solveSharedSupply( ...
            demand_total_m_s, psis(KT), sumgsr, gravity_MPa, ...
            psi_table, K_table, K_min, psi_min, psi_max);

        actual_total_m_s = min(max(actual_total_m_s, 0), demand_total_m_s);
        actual_sun_m_s = actual_total_m_s * demand_sun_m_s / demand_total_m_s;
        actual_sha_m_s = actual_total_m_s * demand_sha_m_s / demand_total_m_s;
    end

    psir(KT) = psis(KT) - actual_total_m_s / sumgsr;
    pstem = psir(KT) - gravity_MPa;
    if ~isreal(pstem) || ~isfinite(pstem)
        error('PlantHydraulic:BadStemPotential', ...
            'Stem water potential is not real/finite. Check KL_h, hh, and transpiration demand.');
    end

    pstem_lookup = min(max(pstem, psi_min), psi_max);
    K_pstem = interp1(psi_table, K_table, pstem_lookup, 'linear');

    Tp_sun(KT) = actual_sun_m_s * 100;               % m/s -> cm/s
    Tp_sha(KT) = actual_sha_m_s * 100;               % m/s -> cm/s
    psil_sun(KT) = leafPotentialFromFlow(actual_sun_m_s, K_pstem, K_unique, psi_for_K, K_min, K_max);
    psil_sha(KT) = leafPotentialFromFlow(actual_sha_m_s, K_pstem, K_unique, psi_for_K, K_min, K_max);
end

function validateLookupTable(psi_table, K_table)
    if numel(psi_table) ~= numel(K_table)
        error('PlantHydraulic:BadLookupTable', ...
            'PSIL_RANGE and KS must have the same number of elements.');
    end
    if ~isreal(psi_table) || ~isreal(K_table) || ...
            any(~isfinite(psi_table)) || any(~isfinite(K_table))
        error('PlantHydraulic:BadLookupTable', ...
            'PSIL_RANGE and KS must be real and finite.');
    end
    if any(diff(psi_table(:)) == 0)
        error('PlantHydraulic:BadLookupTable', ...
            'PSIL_RANGE must not contain duplicate water potentials.');
    end
    if any(K_table < 0)
        error('PlantHydraulic:BadLookupTable', ...
            'KS must be non-negative.');
    end
end

function actual_total_m_s = solveSharedSupply( ...
    demand_total_m_s, psis, sumgsr, gravity_MPa, psi_table, K_table, K_min, psi_min, psi_max)

    supply_at_demand = sharedSupply( ...
        demand_total_m_s, psis, sumgsr, gravity_MPa, psi_table, K_table, K_min, psi_min, psi_max);

    if demand_total_m_s <= supply_at_demand
        actual_total_m_s = demand_total_m_s;
        return;
    end

    lo = 0;
    hi = demand_total_m_s;
    flow_tol = max(1e-14, demand_total_m_s * 1e-6);
    for iter = 1:80
        mid = 0.5 * (lo + hi);
        supply_mid = sharedSupply( ...
            mid, psis, sumgsr, gravity_MPa, psi_table, K_table, K_min, psi_min, psi_max);

        if mid <= supply_mid
            lo = mid;
        else
            hi = mid;
        end

        if (hi - lo) <= flow_tol
            break;
        end
    end

    actual_total_m_s = lo;
end

function supply_m_s = sharedSupply(flow_m_s, psis, sumgsr, gravity_MPa, psi_table, K_table, K_min, psi_min, psi_max)
    pstem = psis - flow_m_s / sumgsr - gravity_MPa;
    if ~isreal(pstem) || ~isfinite(pstem)
        error('PlantHydraulic:BadStemPotential', ...
            'Stem water potential is not real/finite during shared supply solve.');
    end

    pstem_lookup = min(max(pstem, psi_min), psi_max);
    K_pstem = interp1(psi_table, K_table, pstem_lookup, 'linear');
    supply_m_s = max(K_pstem - K_min, 0);
end

function psil = leafPotentialFromFlow(flow_m_s, K_pstem, K_unique, psi_for_K, K_min, K_max)
    K_target = K_pstem - max(0, flow_m_s);
    if ~isreal(K_target) || ~isfinite(K_target)
        error('PlantHydraulic:BadKTarget', ...
            'Leaf conductance target is not real/finite before inverse lookup.');
    end

    K_target = min(max(K_target, K_min), K_max);
    psil = interp1(K_unique, psi_for_K, K_target, 'linear');
end
