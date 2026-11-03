function plantHydro = assemblePlantHydro(built)
%ASSEMBLEPLANTHYDRO  Reassemble the legacy-shaped plantHydro view for the ebal cluster.
%
%   plantHydro = init.assemblePlantHydro(built) gathers the STATIC plant-hydraulic
%   bundle that vegetation.energy.ebal and its sub-functions (solveEbalHydraulics /
%   calSPAC / Plantrootconductivity / assembleRWU) still consume as one monolithic
%   struct, out of the new chain's already-split containers. It is the faithful
%   equivalent of legacy Init_Plant.m section [F] (src/Init_Plant.m:108-161): in
%   COMPASS.m the plantHydro struct is built once before the loop and never mutated,
%   so this is a constant view — build it once at init and cache it (built.p.plantHydro).
%
%   Because it runs at init time (Y.plant_hydraulic still holds the config init values),
%   the .init.* seeds are frozen from Y.plant_hydraulic here and stay constant even after
%   the prognostic potentials evolve, matching the legacy constant-seed behaviour (a warm
%   start from the evolving state was tested and rejected -- no iteration-count benefit).
%   When hydraulic_scheme is storage_from_psi or storage_dynamic (legacy: dynamic)
%   the V2 plant-water-storage params (dyn_a_retention/W_capacity/
%   dz_rs/dz_sl/dyn_tau_stem_s) are ALSO built here (matching legacy Init_Plant), so the per-step
%   PR-A stepPlantHydraulicsV2 reads them straight off this static cache -- it never re-assembles the
%   legacy plantHydro shape inside the loop. (dyn_K_sat/psi_63/c_weibull of the Weibull demo path
%   stay absent -- that path is not wired.)
%
%   Inputs:
%     built  struct  needs .domain (tot_depth/NN/dz), .parameters.veg.hydraulic
%                    (hydraulic_p50/aa/psil_min/scheme), .parameters.veg.root
%                    (root_max_depth/RAI/rootProfile/root_* params), .p.plant_hydraulic
%                    (PSIL_RANGE/KS/gpmax cache), .Y.plant_hydraulic (init psi seeds).
%
%   Returns:
%     plantHydro  struct  legacy-shaped bundle. Geometry: .RL [cm], .Lmax [cm], .NN [-],
%                 .DeltZ [cm, NL x 1, bottom-up], .RAI [NL x 1], .rootProfile. Vulnerability:
%                 .PSIL_RANGE [MPa], .psil_min [MPa], .KS, .gpmax [m/s/MPa], .p50 [MPa], .aa.
%                 Root hydraulics: .root_krmax, .root_side_mode, .root_radial_impedance_Pr,
%                 .rootContact{enabled,psi_start,psi_end}, .gsr_formula, .r_root_m [m],
%                 optional .root_length_density_m_m3. Seeds: .init{psir,pstem,psil_sun,psil_sha}
%                 [MPa]. Scheme: .hydraulic_scheme.
%
%   Reference: src/Init_Plant.m:108-161 (section [F], the struct this reproduces).
% Editor: Jizixing Zhu; Date: 2026-07-03.
% Change: 2026-07-04 Jizixing Zhu -- build the dynamic-scheme V2 storage params (dyn_a_retention/
%   W_capacity/dz_rs/dz_sl/dyn_tau_stem_s) here when organ storage is on, so PR-A stepPlantHydraulicsV2
%   passes built.p.plantHydro directly instead of re-assembling the legacy shape each step.
    dom = built.domain;
    hy  = built.parameters.veg.hydraulic;
    rt  = built.parameters.veg.root;
    pph = built.p.plant_hydraulic;
    yph = built.Y.plant_hydraulic;

    % --- geometry (Init_Plant [F] :108-111,122-123) ---
    plantHydro.RL         = dom.tot_depth;          % [cm] root-zone bottom depth (= Tot_Depth)
    plantHydro.Lmax       = rt.root_max_depth;      % [cm] maximum rooting depth
    plantHydro.NN         = dom.NN;                 % [-]  node count (NL+1)
    plantHydro.DeltZ      = dom.dz;                 % [cm, NL x 1] layer thickness, bottom-up
    plantHydro.RAI        = rt.RAI;                 % [m2/m2, NL x 1] root area index per layer
    plantHydro.rootProfile = rt.rootProfile;

    % --- vulnerability lookup + params (Init_Plant [F] :124-129) ---
    plantHydro.PSIL_RANGE = pph.PSIL_RANGE;         % [MPa] lowest leaf psi .. 0
    plantHydro.psil_min   = hy.hydraulic_psil_min;  % [MPa]
    plantHydro.KS         = pph.KS;                 % leaf-side cumulative supply table
    plantHydro.gpmax      = pph.gpmax;              % [m/s/MPa] max plant conductance
    plantHydro.p50        = hy.hydraulic_p50;       % [MPa] xylem vulnerability P50
    plantHydro.aa         = hy.hydraulic_aa;        % [-]  vulnerability curve shape

    % --- root hydraulics (Init_Plant [F] :130-148) ---
    plantHydro.root_krmax = rt.root_hydraulic_krmax;
    plantHydro.root_side_mode = lower(strtrim(string(rt.root_side_mode)));
    plantHydro.root_radial_impedance_Pr = rt.root_radial_impedance_Pr;
    plantHydro.rootContact = struct( ...
        'enabled',   rt.root_contact_enabled ~= 0, ...
        'psi_start', rt.root_contact_psi_start, ...   % [MPa]
        'psi_end',   rt.root_contact_psi_end);        % [MPa]
    plantHydro.gsr_formula = lower(strtrim(string(rt.gsr_formula)));
    plantHydro.r_root_m   = max(rt.root_rroot_mm * 1e-3, 1e-5);  % [m] fine-root radius, floored
    if isfield(rt, 'rootProfile') && isstruct(rt.rootProfile) && ...
            isfield(rt.rootProfile, 'root_length_density_m_m3') && ...
            ~isempty(rt.rootProfile.root_length_density_m_m3)
        plantHydro.root_length_density_m_m3 = rt.rootProfile.root_length_density_m_m3;
    end

    % --- warm-start seeds: constant, frozen from Y at init (= config hydraulic_init_*) ---
    plantHydro.init.psir     = yph.psi_root;      % [MPa]
    plantHydro.init.pstem    = yph.psi_stem;      % [MPa]
    plantHydro.init.psil_sun = yph.psi_leaf_sun;  % [MPa]
    plantHydro.init.psil_sha = yph.psi_leaf_sha;  % [MPa]

    % --- scheme selector (canonical names; legacy equilibrium/dynamic accepted) ---
    plantHydro.hydraulic_scheme = vegetation.hydraulics.normalizeHydraulicScheme( ...
        hy.hydraulic_scheme);

    % --- V2 plant-water-storage params (only when organ storage is active) ---
    % Consumed by vegetation.hydraulics.updateOrganStorageFromPsi. Built once here so the

    % per-step wrapper reads this static cache -- no legacy-shape reassembly in the loop.
    if vegetation.hydraulics.usesOrganWaterStorage(plantHydro.hydraulic_scheme)
        plantHydro.dyn_a_retention = [hy.dyn_a_retention_root; hy.dyn_a_retention_stem; hy.dyn_a_retention_leaf];  % [m^-1]
        plantHydro.W_capacity      = [hy.dyn_W_capacity_root;  hy.dyn_W_capacity_stem;  hy.dyn_W_capacity_leaf];   % [m]
        plantHydro.dz_rs           = hy.dyn_dz_rs;        % [m] root->stem path length
        plantHydro.dz_sl           = hy.dyn_dz_sl;        % [m] stem->leaf path length
        plantHydro.dyn_tau_stem_s  = hy.dyn_tau_stem_s;   % [s] stem capacitance lag (0 = instantaneous)
    end
end
