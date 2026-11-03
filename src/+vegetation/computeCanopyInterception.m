function intc = computeCanopyInterception(scheme, precip, dt, LAI, kV, S_L, F0, cw)
%COMPUTECANOPYINTERCEPTION  Pre-evaporation canopy interception, storage, and fwet.
%
%   intc = vegetation.computeCanopyInterception(scheme, precip, dt, LAI, kV, S_L, F0, cw)
%   is the PRE-evaporation half of the canopy-water chain. It (1) splits rainfall
%   with vegetation.calculateRainfallPartitioning (Gash), (2) adds intercepted water to
%   canopy reservoir(s) and spills excess as drip, and (3) returns wet-canopy
%   fractions fwet = (h2ocan/max_storage)^(2/3). Pure: no Y/runtime mutation; the
%   caller stores intc (e.g. workspace.scope_step.interception). Reservoir drawdown
%   and Y.canopy_water commit are post-ebal (vegetation.computeCanopyWaterEvap + accept).
%
%   Returns:
%     intc [struct] uniform fields for both schemes:
%       .scheme            [char]   'single' | 'four_quadrant'
%       .fwet              [-]      bulk wet fraction (single value or Ps-weighted)
%       .fwet_sun          [-]      sunlit wet fraction
%       .fwet_sha          [-]      shaded wet fraction
%       .h2ocan            [mm]     total pre-evaporation storage
%       .h2ocan_sun        [mm]     sunlit reservoir (0 in single)
%       .h2ocan_sha        [mm]     shaded reservoir (0 in single)
%       .max_storage       [mm]     total capacity
%       .max_storage_sun   [mm]     sunlit capacity
%       .max_storage_sha   [mm]     shaded capacity
%       .throughfall_cms   [cm/s]   free throughfall + drip to ground
%       .CI_potential      [cm/s]   Gash interception rate
%       .drip_mm           [mm]     ground drip (single) or drip_to_ground (alias)
%       .drip_sun_to_sha   [mm]     sunlit overflow to shaded (0 in single)
%       .drip_to_ground    [mm]     shaded overflow to ground (0 in single if none)
%
%   Inputs:
%     scheme  [in]  char       'single' | 'four_quadrant'
%     precip  [in]  [cm/s]     rainfall rate this outer step
%     dt      [in]  [s]        outer step length
%     LAI     [in]  [m2/m2]    leaf area index
%     kV      [in]  [-]        canopy extinction coefficient (Gash)
%     S_L     [in]  [mm/LAI]   storage capacity per unit LAI
%     F0      [in]  [-]        free-throughfall / wet-canopy ratio param for Gash
%     cw      [in]  struct     caller-assembled param-view (no runtime read):
%                              .state_h2ocan [mm], .state_h2ocan_sun [mm],
%                              .state_h2ocan_sha [mm] from Y.canopy_water;
%                              .Ps_max_diurnal [-], .Ps_min_diurnal [-] from
%                              params.veg.canopy (fixed capacity split)
%
%   Method:
%     Call Gash for CI and free throughfall rates. Then:
%       single — one reservoir max = S_L*LAI; add CI*dt*10 [mm]; spill to ground;
%         fwet_sun = fwet_sha = fwet = wetFraction(h2o, max).
%       four_quadrant — capacities max_sun = S_L*LAI*Ps_max,
%         max_sha = S_L*LAI*(1-Ps_min); split CI mass by Ps_max into sun/sha;
%         Rutter cascade: sun spill -> sha, sha spill -> ground throughfall;
%         independent fwet_sun/fwet_sha; bulk fwet = Ps_max-weighted average.
%     Capacities use fixed Ps_max/min (not instantaneous sunlit fraction) to avoid
%     artificial mass transfer when Ps changes.
%
%     Four latent-heat pathways (dry-sun / wet-sun / dry-shade / wet-shade) are
%     NOT assembled here. This function only supplies fwet_sun and fwet_sha (and
%     storages). vegetation.energy.ebal / ebalProfile pass those into
%     heatfluxes_clm45 separately for sunlit and shaded leaves; heatfluxes_clm45
%     then splits each call into dry (1-fwet) stomatal and wet (fwet) free-water
%     latent heat. Post-ebal drawdown: vegetation.computeCanopyWaterEvap.
%
%   Reference:
%     Legacy port of COMPASS.m interception block (approx. lines 733-837).
%     Gash bulk partitioning: vegetation.calculateRainfallPartitioning (Gash 1979 form;
%       see that file for SOURCE UNVERIFIED note if local PDF absent).
%     fwet exponent 2/3: CLM5/NoahMP-style wet fraction; form commonly linked to
%       Deardorff (1978) — [SOURCE UNVERIFIED against local PDF 2026-07-16].
%     Dry/wet leaf LE split: heatfluxes_clm45.m (CLM5 Tech Note-style Eq 5.97 path
%       as documented in that file); sun/sha fwet wiring: ebal.m /
%       ebalProfile.m interception parse + heatfluxes calls.
%     Cascaded sun->sha storage: Rutter-style multi-layer interception concept —
%       [SOURCE UNVERIFIED against local PDF 2026-07-16]; COMPASS implementation
%       is the v30 four_quadrant scheme (Jizixing Zhu).
%
%   Editor: Jizixing Zhu; Date: 2026-07-02.
%   Change: 2026-07-05 Jizixing Zhu — doc: cw is a caller-assembled param-view
%           (state from Y.canopy_water, config from params.veg.canopy); no runtime read.
%   Change: 2026-07-16 Jizixing Zhu — compass-coding-standard English doc: Method
%           clarifies pre-evap scope; four-quadrant fwet feeds dry/wet sun/shade LE
%           in ebal/heatfluxes, not inside this function; references + SOURCE UNVERIFIED.

    [CI_cms, TF_cms] = vegetation.calculateRainfallPartitioning(precip, dt, LAI, kV, S_L, F0);

    if strcmpi(scheme, 'four_quadrant')
        % Sunlit/shaded reservoirs + Rutter cascade (pre-evap state only).
        maxSun = S_L * LAI * cw.Ps_max_diurnal;
        maxSha = S_L * LAI * (1 - cw.Ps_min_diurnal);

        CI_total_mm = CI_cms * dt * 10;                       % [mm] potential interception
        CI_to_sun = CI_total_mm * cw.Ps_max_diurnal;
        CI_to_sha = CI_total_mm * (1 - cw.Ps_max_diurnal);

        h2oSun = cw.state_h2ocan_sun + CI_to_sun;
        [h2oSun, dripSunToSha] = spill(h2oSun, maxSun);       % sunlit overflow -> shaded

        h2oSha = cw.state_h2ocan_sha + CI_to_sha + dripSunToSha;
        [h2oSha, dripToGround] = spill(h2oSha, maxSha);       % shaded overflow -> ground

        TF_cms = TF_cms + dripToGround / dt / 10;             % drip joins throughfall
        fwetSun = wetFraction(h2oSun, maxSun);
        fwetSha = wetFraction(h2oSha, maxSha);
        fwetTot = fwetSun * cw.Ps_max_diurnal + fwetSha * (1 - cw.Ps_max_diurnal);

        intc = packInterception('four_quadrant', fwetTot, fwetSun, fwetSha, ...
            h2oSun, h2oSha, maxSun, maxSha, TF_cms, CI_cms, dripToGround, dripSunToSha, dripToGround);
    else
        % Single reservoir + single fwet (CLM5/NoahMP-style baseline).
        maxStore = S_L * LAI;
        h2o = cw.state_h2ocan + CI_cms * dt * 10;               % [mm]
        [h2o, dripMm] = spill(h2o, maxStore);
        TF_cms = TF_cms + dripMm / dt / 10;
        fwet = wetFraction(h2o, maxStore);
        % Mirror sun/sha fields so the output struct is scheme-uniform for ebal.
        intc = packInterception('single', fwet, fwet, fwet, ...
            h2o, 0, maxStore, 0, TF_cms, CI_cms, dripMm, 0, dripMm);
    end
end

% ---------------------------------------------------------------------------
function [h2o, drip] = spill(h2o, maxStore)
%SPILL  Cap reservoir at maxStore; excess is drip [mm] (0 if maxStore <= 0).
% Editor: Jizixing Zhu; Date: 2026-07-02.
    if h2o > maxStore && maxStore > 0
        drip = h2o - maxStore;
        h2o  = maxStore;
    else
        drip = 0;
    end
end

% ---------------------------------------------------------------------------
function f = wetFraction(h2o, maxStore)
%WETFRACTION  Wet-canopy fraction (h2o/max)^(2/3) in [0,1]; 0 if no capacity.
% Editor: Jizixing Zhu; Date: 2026-07-02.
    if maxStore > 1e-9
        f = min(1.0, (h2o / maxStore) ^ (2/3));
    else
        f = 0;
    end
end

% ---------------------------------------------------------------------------
function intc = packInterception(scheme, fwet, fwetSun, fwetSha, h2oSun, h2oSha, ...
        maxSun, maxSha, throughfallCms, CIpotential, dripMm, dripSunToSha, dripToGround)
%PACKINTERCEPTION  Fixed-field interception struct shared by both schemes.
% Editor: Jizixing Zhu; Date: 2026-07-02.
    intc = struct( ...
        'scheme',          scheme, ...
        'fwet',            fwet, ...
        'fwet_sun',        fwetSun, ...
        'fwet_sha',        fwetSha, ...
        'h2ocan',          h2oSun + h2oSha, ...
        'h2ocan_sun',      h2oSun, ...
        'h2ocan_sha',      h2oSha, ...
        'max_storage',     maxSun + maxSha, ...
        'max_storage_sun', maxSun, ...
        'max_storage_sha', maxSha, ...
        'throughfall_cms', throughfallCms, ...
        'CI_potential',    CIpotential, ...
        'drip_mm',         dripMm, ...
        'drip_sun_to_sha', dripSunToSha, ...
        'drip_to_ground',  dripToGround);
end
