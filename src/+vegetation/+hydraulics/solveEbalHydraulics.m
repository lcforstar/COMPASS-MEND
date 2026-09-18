function hydraulic = solveEbalHydraulics(plantHydro, soil_state, soil_conductivity, KT, H, Tp_sun_cm_s, Tp_sha_cm_s, soil_params)
%SOLVEEBALHYDRAULICS Solve plant hydraulics from ebal sunlit/shaded transpiration.
%
%   Inside ebal's Picard loop this function ALWAYS runs the v1 steady SPAC
%   scheme (calSPAC bisection). The hydraulic_scheme switch only controls
%   what runs OUTSIDE this function:
%     'steady_spac'      — v1 only (no organ-storage wrapper)
%     'storage_from_psi' — v1 inside ebal + W diagnosed from psi after ebal
%     'storage_dynamic'  — same off-loop W path today (true dW/dt later)
%   Legacy aliases: equilibrium -> steady_spac; dynamic -> storage_from_psi.
%   Organ-storage W does NOT feed back into the current KT's ebal.
%
%   Inputs:
%     soil_state         .h [cmH2O], .theta_liquid [m3 m-3] (target names)
%     soil_conductivity  KL/KfL field used as Khh proxy [NL x 2]
%     soil_params        optional .theta_s for contact / stress helpers
%
%   This separation keeps ebal Picard numerically stable (no W state inside
%   the iteration) while still providing W storage as a prognostic output.
% Editor: Jizixing Zhu; Date: 2026-07-16.
% Change: 2026-07-16 Jizixing Zhu — accept soil_state.h / soil_params.theta_s (drop SoilVariables.hh bridge).

    hydraulic = struct( ...
        'Tp_sun', 0, ...
        'Tp_sha', 0, ...
        'Tp_t', 0, ...
        'Tp_sun_demand', 0, ...
        'Tp_sha_demand', 0, ...
        'Tp_demand', 0, ...
        'supply_ratio', 1, ...
        'psil_sun', plantHydro.init.psil_sun, ...
        'psil_sha', plantHydro.init.psil_sha, ...
        'psil', 0.5 * (plantHydro.init.psil_sun + plantHydro.init.psil_sha), ...
        'psir', plantHydro.init.psir, ...
        'psis', plantHydro.init.psir, ...
        'wt_sun', 1, ...
        'wt_sha', 1, ...
        'gsr_total', 0, ...
        'Kplant_diff_sun', 0, ...
        'Kplant_diff_sha', 0, ...
        'Kplant_diff_flux_weighted', 0, ...
        'Kplant_app_sun', NaN, ...
        'Kplant_app_sha', NaN, ...
        'Kplant_app_total', NaN, ...
        'Kwhole_app_total', NaN, ...
        'Ksoilroot_app_total', NaN, ...
        'gsr', zeros(plantHydro.NN, 1));

    if nargin < 7 || KT < 1 || ~isstruct(soil_state) || ~isfield(soil_state, 'h')
        return;
    end

    if isempty(soil_state.h)
        return;
    end

    Tp_sun = zeros(KT, 1);
    Tp_sha = zeros(KT, 1);
    Tp_sun(KT) = max(0, Tp_sun_cm_s);
    Tp_sha(KT) = max(0, Tp_sha_cm_s);

    if nargin < 3 || isempty(soil_conductivity)
        error('PlantHydraulic:MissingKLh', ...
            'soil_conductivity is required for plant hydraulics.');
    end

    NL = plantHydro.NN - 1;
    if ~isequal(size(soil_conductivity), [NL, 2])
        error('PlantHydraulic:BadKLhSize', ...
            'soil_conductivity must be NL x 2. Got %d x %d, expected %d x 2.', ...
            size(soil_conductivity, 1), size(soil_conductivity, 2), NL);
    end
    if ~isreal(soil_conductivity) || any(~isfinite(soil_conductivity(:))) || any(soil_conductivity(:) < 0)
        error('PlantHydraulic:BadKLhValue', ...
            'soil_conductivity must be real, finite, and non-negative.');
    end

    Khh_proxy = soil_conductivity;
    plantHydro_step = plantHydro;
    if nargin >= 8 && isstruct(soil_params) && ...
            isfield(soil_state, 'theta_liquid') && isfield(soil_params, 'theta_s')
        % Bottom-up soil order (layer 1 deepest, NL surface).
        theta_layer = soil_state.theta_liquid(1:NL, min(2, size(soil_state.theta_liquid, 2)));
        plantHydro_step.soil_theta_layer = theta_layer(:);
        plantHydro_step.soil_theta_s_layer = soil_params.theta_s(:);
    end

    [Tp_sun, Tp_sha, psil_sun, psil_sha, psir, psis, wt_sun, wt_sha, gsr] = ...
        vegetation.hydraulics.Plant_hydraulics( ...
            plantHydro_step, Khh_proxy, soil_state.h(:), Tp_sun, Tp_sha, KT, H);

    hydraulic.Tp_sun = Tp_sun(KT);
    hydraulic.Tp_sha = Tp_sha(KT);
    hydraulic.Tp_t = Tp_sun(KT) + Tp_sha(KT);
    hydraulic.Tp_sun_demand = max(0, Tp_sun_cm_s);
    hydraulic.Tp_sha_demand = max(0, Tp_sha_cm_s);
    hydraulic.Tp_demand = max(0, Tp_sun_cm_s) + max(0, Tp_sha_cm_s);
    if hydraulic.Tp_demand > 0
        hydraulic.supply_ratio = hydraulic.Tp_t / hydraulic.Tp_demand;
    else
        hydraulic.supply_ratio = 1;
    end
    hydraulic.psil_sun = psil_sun(KT);
    hydraulic.psil_sha = psil_sha(KT);
    hydraulic.psir = psir(KT);
    hydraulic.psis = psis(KT);
    hydraulic.gsr  = gsr;        % NN×1 [m/s/MPa]，各节点土壤-根系导度
    hydraulic.wt_sun = wt_sun;
    hydraulic.wt_sha = wt_sha;

    if hydraulic.Tp_t > 0
        hydraulic.psil = (hydraulic.Tp_sun * hydraulic.psil_sun + hydraulic.Tp_sha * hydraulic.psil_sha) / hydraulic.Tp_t;
    else
        hydraulic.psil = 0.5 * (hydraulic.psil_sun + hydraulic.psil_sha);
    end

    hydraulic.gsr_total = sum(gsr);

    % Diagnostic conductances.  The vulnerability curve gives the local
    % differential plant conductance, while the apparent values are inferred
    % from the realized flow and realized potential drops at this time step.
    gravity_MPa = 1000 * 9.80665 * H / 1e6;
    pstem = hydraulic.psir - gravity_MPa;
    Tp_sun_m_s = hydraulic.Tp_sun / 100;
    Tp_sha_m_s = hydraulic.Tp_sha / 100;
    Tp_total_m_s = hydraulic.Tp_t / 100;

    hydraulic.Kplant_diff_sun = vegetation.hydraulics.VulnerabilityCurve( ...
        hydraulic.psil_sun, plantHydro.gpmax, plantHydro.p50, plantHydro.aa);
    hydraulic.Kplant_diff_sha = vegetation.hydraulics.VulnerabilityCurve( ...
        hydraulic.psil_sha, plantHydro.gpmax, plantHydro.p50, plantHydro.aa);
    if Tp_total_m_s > 0
        hydraulic.Kplant_diff_flux_weighted = ...
            (Tp_sun_m_s * hydraulic.Kplant_diff_sun + ...
             Tp_sha_m_s * hydraulic.Kplant_diff_sha) / Tp_total_m_s;
    else
        hydraulic.Kplant_diff_flux_weighted = ...
            0.5 * (hydraulic.Kplant_diff_sun + hydraulic.Kplant_diff_sha);
    end

    hydraulic.Kplant_app_sun = conductanceFromDrop(Tp_sun_m_s, pstem - hydraulic.psil_sun);
    hydraulic.Kplant_app_sha = conductanceFromDrop(Tp_sha_m_s, pstem - hydraulic.psil_sha);
    hydraulic.Kplant_app_total = conductanceFromDrop(Tp_total_m_s, pstem - hydraulic.psil);
    hydraulic.Kwhole_app_total = conductanceFromDrop( ...
        Tp_total_m_s, hydraulic.psis - hydraulic.psil - gravity_MPa);
    hydraulic.Ksoilroot_app_total = conductanceFromDrop( ...
        Tp_total_m_s, hydraulic.psis - hydraulic.psir);
end

function K = conductanceFromDrop(flow_m_s, delta_psi_MPa)
    if flow_m_s > 0 && isfinite(delta_psi_MPa) && delta_psi_MPa > 0
        K = flow_m_s / delta_psi_MPa;
    else
        K = NaN;
    end
end
