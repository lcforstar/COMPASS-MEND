function soil = buildSoilThermal(soil)
%BUILDSOILTHERMAL  Stage B3: static soil thermal parameters (de Vries / Tian /
%                  Peters-Lidard / Farouki), the dry-state inputs the runtime
%                  thermal kernel needs each step.
%
%   soil = init.buildSoilThermal(soil) ports legacy SoilProperties.
%   compute_thermal_parameters verbatim, operating per layer on the already
%   bottom-up soil parameter struct. Because every formula is LOCAL to a layer,
%   computing directly on the bottom-up arrays equals the legacy "compute then
%   flip" (flip(f(x)) == f(flip(x))). Depends on stage B2 (wilting_point, used
%   in GB2) and on the texture->quartz lookup done here.
%
%   Inputs:
%     soil   struct  built.parameters.soil with [NL x 1] porosity, msoc_volume,
%                    wilting_point, bulk_density, sand, silt, clay
%
%   Returns soil augmented with (bottom-up):
%     quartz             [-]            quartz fraction of solids (texture lookup)
%     hcap_phase [1x6]   [J/cm3/C]      phase volumetric heat capacities
%     sf_phase   [1x6]   [-]            phase shape factors
%     tca,ga1,ga2        [J/cm/s/C],[-] dry-air conductivity + gas shape slopes
%     heat_capacity_dry  [J/cm3/C]      dry-soil volumetric heat capacity (de Vries)
%     s1,s2,zeta0,con0,ps1,ps2,gb1,gb2  de Vries (1963) coefficients
%     tps1,tps2,tcon0    Tian (2016) particle-size weighted mineral conductivity
%     tcon_dry,tcon_s    Peters-Lidard (1998)/Johansen dry + solid conductivity [W/m/K]
%     fehcap,tcon_s_farouki  Farouki dry heat capacity [J/m3/K] + solid cond. [W/m/K]
%
%   Reference: legacy SoilProperties.compute_thermal_parameters (de Vries 1963;
%   Tian 2016; Peters-Lidard et al. 1998; Farouki 1981) + classify_soil_texture
%   (USDA triangle) + calculate_quartz_fraction (Peters-Lidard quartz table).
% Editor: Jizixing Zhu; Date: 2026-06-27.

    NL = numel(soil.theta_s);

    % --- de Vries phase material constants (CGS, J cm^-1 s^-1 C^-1) ---
    % index 1=liquid water 2=dry air 3=quartz 4=other mineral 5=organic 6=ice
    hcap = [0.998 * 4.182, 0.0003 * 4.182, 0.46 * 4.182, 0.46 * 4.182, 0.6 * 4.182, 0.45 * 4.182];
    tcon = [1.37e-3 * 4.182, 6e-5 * 4.182, 2.1e-2 * 4.182, 7e-3 * 4.182, 6e-4 * 4.182, 5.2e-3 * 4.182];
    sf   = [0, 0, 0.125, 0.125, 0.5, 0.125];
    tca  = 6e-5 * 4.182;     % dry-air thermal conductivity
    ga1  = 0.035;
    ga2  = 0.013;

    % --- Johansen/Tian mineral conductivities (SI, W/(m.K)) ---
    TCON_qtz = 7.7;    TCON_o   = 2.0;
    TCON_sa  = 7.7;    SF_sa  = 0.182;
    TCON_slt = 2.74;   SF_slt = 0.0534;
    TCON_cl  = 1.93;   SF_cl  = 0.00775;
    TCON_Soc = 0.05;   % organic-matter conductivity

    % phase constants (runtime kernel reads these directly)
    soil.hcap_phase = hcap;
    soil.sf_phase   = sf;
    soil.tca        = tca;
    soil.ga1        = ga1;
    soil.ga2        = ga2;

    % --- texture class -> quartz fraction (Peters-Lidard 1998 table) ---
    quartz_db = [0.92, 0.82, 0.60, 0.40, 0.25, 0.60, 0.35, 0.10, 0.52, 0.10, 0.25];
    quartz = zeros(NL, 1);
    for j = 1:NL
        tc = classifyTexture(soil.sand(j), soil.silt(j), soil.clay(j));
        quartz(j) = quartz_db(tc + 1);
    end
    soil.quartz = quartz;

    [heat_capacity_dry, s1, s2, zeta0, con0, ps1, ps2, gb1, gb2, ...
     tps1, tps2, tcon0, tcon_dry, tcon_s, fehcap, tcon_s_farouki] = deal(zeros(NL, 1));

    for j = 1:NL
        por_j   = soil.porosity(j);
        xsoc_j  = soil.msoc_volume(j);     % organic volume fraction
        xwilt_j = soil.wilting_point(j);   % from stage B3
        qtz_j   = quartz(j);
        rho_b_j = soil.bulk_density(j);    % g/cm3
        sand_j  = soil.sand(j);
        silt_j  = soil.silt(j);
        clay_j  = soil.clay(j);
        solid_j = 1 - por_j;

        % component volumes (de Vries / Tian shared)
        min_solid_j  = (1 - xsoc_j) * solid_j;
        vper = [qtz_j * min_solid_j, (1 - qtz_j) * min_solid_j, xsoc_j * solid_j];

        % [1] de Vries (1963): ZETA0, CON0, PS1, PS2, HCD, GB1, GB2
        s1v = por_j * tca; s2v = por_j; hcd_j = 0; ps1_j = 0; ps2_j = 0;
        for i = 1:3
            idx = i + 2;   % 3=quartz, 4=mineral, 5=organic
            v = vper(i);
            targ1 = tcon(idx) / tca - 1;
            g1 = 2/3 / (1 + targ1 * sf(idx)) + 1/3 / (1 + targ1 * (1 - 2 * sf(idx)));
            s1v   = s1v + g1 * tcon(idx) * v;
            s2v   = s2v + g1 * v;
            hcd_j = hcd_j + hcap(idx) * v;
            targ2 = tcon(idx) / tcon(1) - 1;
            g2 = 2/3 / (1 + targ2 * sf(idx)) + 1/3 / (1 + targ2 * (1 - 2 * sf(idx)));
            ps1_j = ps1_j + g2 * v * tcon(idx);
            ps2_j = ps2_j + g2 * v;
        end
        heat_capacity_dry(j) = hcd_j;
        s1(j) = s1v; s2(j) = s2v;
        zeta0(j) = 1 / s2v;
        con0(j)  = 1.25 * s1v / s2v;
        ps1(j) = ps1_j; ps2(j) = ps2_j;
        slope_term = 1/3 - ga1;
        gb1(j) = slope_term / por_j;
        gb2(j) = (ga1 - ga2) / xwilt_j + slope_term / por_j;

        % [2] Tian (2016): particle-size weighted mineral conductivity
        tcon_min_cgs = (TCON_sa^sand_j * TCON_slt^silt_j * TCON_cl^clay_j) / 100;
        SF_min_j     = SF_sa * sand_j + SF_slt * silt_j + SF_cl * clay_j;
        ttarg1 = tcon_min_cgs / tca - 1;
        tgrat  = 2/3 / (1 + ttarg1 * SF_min_j) + 1/3 / (1 + ttarg1 * (1 - 2 * SF_min_j));
        ts1 = por_j * tca + tgrat * tcon_min_cgs * solid_j;
        ts2 = por_j + tgrat * solid_j;
        tcon0(j) = 1.25 * ts1 / ts2;
        ttarg2 = tcon_min_cgs / tcon(1) - 1;
        tgrat2 = 2/3 / (1 + ttarg2 * SF_min_j) + 1/3 / (1 + ttarg2 * (1 - 2 * SF_min_j));
        tterm  = tgrat2 * solid_j;
        tps1(j) = tterm * tcon_min_cgs;
        tps2(j) = tterm;

        % [3] Peters-Lidard (1998) / Johansen: TCON_dry, TCON_s (SI)
        RHO_bulk_SI = rho_b_j * 1000;   % g/cm3 -> kg/m3
        tcon_dry_j = (0.135 * RHO_bulk_SI + 64.7) / (2700 - 0.947 * RHO_bulk_SI);
        tcon_dry_j = tcon_dry_j * (1 - xsoc_j) + xsoc_j * TCON_Soc;   % organic correction
        tcon_dry(j) = tcon_dry_j;
        tcon_s(j) = TCON_qtz^qtz_j * TCON_o^(1 - qtz_j);   % Johansen geometric mean

        % [4] Farouki: FEHCAP dry heat capacity + TCON_s_Farouki solid conductivity
        if (sand_j + clay_j) > 0
            fehcap_j = (2.128 * sand_j + 2.385 * clay_j) / (sand_j + clay_j) * 1e6;
            fehcap_j = fehcap_j * (1 - xsoc_j) + xsoc_j * 2.5e6;
            tcon_s_far = (8.8 * sand_j + 2.92 * clay_j) / (sand_j + clay_j);   % linear weight (NOT geometric)
            tcon_s_far = tcon_s_far * (1 - xsoc_j) + xsoc_j * 0.25;
        else
            fehcap_j   = 2.0e6;
            tcon_s_far = tcon_s(j);   % degenerate -> Johansen value
        end
        fehcap(j)         = fehcap_j;
        tcon_s_farouki(j) = tcon_s_far;
    end

    soil.heat_capacity_dry = heat_capacity_dry;
    soil.s1 = s1; soil.s2 = s2; soil.zeta0 = zeta0; soil.con0 = con0;
    soil.ps1 = ps1; soil.ps2 = ps2; soil.gb1 = gb1; soil.gb2 = gb2;
    soil.tps1 = tps1; soil.tps2 = tps2; soil.tcon0 = tcon0;
    soil.tcon_dry = tcon_dry; soil.tcon_s = tcon_s;
    soil.fehcap = fehcap; soil.tcon_s_farouki = tcon_s_farouki;
end

% ---------------------------------------------------------------------------
function tc = classifyTexture(sand, silt, clay)
%CLASSIFYTEXTURE  USDA texture-triangle class index 0..10 from fractions (0-1).
%   Verbatim from legacy SoilProperties.classify_soil_texture.
% Editor: Jizixing Zhu; Date: 2026-06-27.
    sand = sand * 100; silt = silt * 100; clay = clay * 100;
    if silt + 1.5 * clay < 15
        tc = 0;   % Sand
    elseif silt + 1.5 * clay >= 15 && silt + 2.0 * clay < 30
        tc = 1;   % LoamySand
    elseif (clay >= 7 && clay < 20 && sand > 52 && silt + 2.0 * clay >= 30) || ...
            (clay < 7 && silt < 50 && silt + 2.0 * clay >= 30)
        tc = 2;   % SandyLoam
    elseif clay >= 7 && clay < 27 && silt >= 28 && silt < 50 && sand <= 52
        tc = 3;   % Loam
    elseif (silt >= 50 && clay >= 12 && clay < 27) || (silt >= 50 && silt < 80 && clay < 12)
        tc = 4;   % SiltLoam
    elseif clay >= 20 && clay < 35 && silt < 28 && sand > 45
        tc = 5;   % SandyClayLoam
    elseif clay >= 27 && clay < 40 && sand > 20 && sand <= 45
        tc = 6;   % ClayLoam
    elseif clay >= 27 && clay < 40 && sand <= 20
        tc = 7;   % SiltyClayLoam
    elseif clay >= 35 && sand > 45
        tc = 8;   % SandyClay
    elseif clay >= 40 && silt >= 40
        tc = 9;   % SiltyClay
    elseif clay >= 40 && sand <= 45 && silt < 40
        tc = 10;  % Clay
    else
        tc = 3;   % default Loam
    end
end
