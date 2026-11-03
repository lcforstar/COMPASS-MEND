function [yph, pph] = buildPlantHydraulicState(raw, built)
%BUILDPLANTHYDRAULICSTATE  Plant hydraulic init state + vulnerability cache (P1).
%
%   [yph, pph] = init.buildPlantHydraulicState(raw, built) builds the plant
%   hydraulic prognostic init state and the recomputable vulnerability lookup
%   cache, faithfully reproducing the hydraulic part of the legacy Init_Plant
%   script. Caller stores yph at built.Y.plant_hydraulic and pph at
%   built.p.plant_hydraulic.
%
%   Inputs:
%     raw    struct  raw.config.vegetation (hydraulic_init_psi* initial values)
%     built  struct  built.parameters.veg.hydraulic (gpmax/p50/aa/psil_min/scheme/
%                    dyn_W_capacity_*/tree_*) + .canopy.canopy_height (gpmax fallback)
%
%   Returns:
%     yph  struct  STATE: psi_root/psi_stem/psi_leaf_sun/psi_leaf_sha [MPa] +
%                  W_root/W_stem/W_leaf [m] (when storage_from_psi/storage_dynamic)
%     pph  struct  CACHE (recomputable): PSIL_RANGE [MPa], VC, KS, gpmax
%
%   Reuses vegetation.hydraulics.VulnerabilityCurve / initWaterStorage / calculate_gpmax
%   (== Init_Plant). RAI is built separately (init.buildVeg V2). Ref: src/Init_Plant.m.
% Editor: Jizixing Zhu; Date: 2026-06-27.
% Change: 2026-07-01 Jizixing Zhu — restore Init_Plant BadInitialWaterPotential gate (init psi finite and <= 0 MPa).

    hy = built.parameters.veg.hydraulic;
    vc = raw.config.vegetation;

    % --- effective gpmax: config value if > 0, else estimate from tree structure ---
    gpmax = hy.hydraulic_gpmax;
    if gpmax <= 0
        H_ref = built.parameters.veg.canopy.canopy_height;       % mean(hc) = canopy_height (constant)
        CA    = pi * (hy.tree_CAH / 2)^2;                         % crown area [m2]
        SCA   = pi * (hy.tree_DBH / 100)^2 / 4;                   % stem cross-section [m2]
        gpmax = vegetation.hydraulics.calculate_gpmax(H_ref, hy.tree_Ksmax, CA, SCA, 1000);
    end

    % --- vulnerability lookup cache (recomputable) ---
    DPSI = 0.01;
    PSIL_RANGE = hy.hydraulic_psil_min : DPSI : 0;                % [MPa] lowest leaf psi .. 0
    VC = vegetation.hydraulics.VulnerabilityCurve(PSIL_RANGE, gpmax, hy.hydraulic_p50, hy.hydraulic_aa);
    KS = cumsum(VC * DPSI);
    pph = struct('PSIL_RANGE', PSIL_RANGE, 'VC', VC, 'KS', KS, 'gpmax', gpmax);

    % --- prognostic init state: water potentials (MPa) ---
    yph = struct();
    yph.psi_root     = getf(vc, 'hydraulic_init_psir',     -0.2);
    yph.psi_stem     = getf(vc, 'hydraulic_init_pstem',    -0.25);
    yph.psi_leaf_sun = getf(vc, 'hydraulic_init_psil_sun', -0.3);
    yph.psi_leaf_sha = getf(vc, 'hydraulic_init_psil_sha', -0.3);

    % Validation (Init_Plant.m:50-54): initial plant water potentials must be
    % finite and <= 0 MPa (getf already floors non-finite config to the default,
    % so this mainly rejects a positive hydraulic_init_psi* config value).
    init_psi = [yph.psi_root, yph.psi_stem, yph.psi_leaf_sun, yph.psi_leaf_sha];
    if any(~isfinite(init_psi)) || any(init_psi > 0)
        error('buildPlantHydraulicState:badInitialWaterPotential', ...
            'hydraulic_init_psi* must be finite and <= 0 MPa.');
    end

    % --- organ water storage W [m] (90% of capacity) when storage scheme on ---
    if vegetation.hydraulics.usesOrganWaterStorage(hy.hydraulic_scheme)
        Wcap = [hy.dyn_W_capacity_root; hy.dyn_W_capacity_stem; hy.dyn_W_capacity_leaf];
        W = vegetation.hydraulics.initWaterStorage(struct('W_capacity', Wcap), []);
        yph.W_root = W(1);  yph.W_stem = W(2);  yph.W_leaf = W(3);
    end
end


% ---------------------------------------------------------------------------
function val = getf(s, name, default)
% Field by name with fallback; treat empty/non-finite as unset (legacy convention).
% Editor: Jizixing Zhu; Date: 2026-06-27.
    val = default;
    if isfield(s, name)
        v = s.(name);
        if ~isempty(v) && isnumeric(v) && isfinite(v); val = v; end
    end
end
