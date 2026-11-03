function [psis2, gsr, sumgsr, contactStats] = Plantrootconductivity(RL, Lmax, NN, DeltZ, Khh, RAI, hh, KT, rootContact, plantHydro)
%PLANTROOTCONDUCTIVITY Root-weighted soil water potential and conductance.
%   Khh is soil hydraulic conductivity at each layer endpoint (NL x 2,
%   cm/s). It is converted to m/s inside this function.
%
%   RAI is root area index at the same layer endpoints (NL x 2, m2/m2).
%   For compatibility, an NL x 1 layer-total RAI profile is split
%   conservatively over the two layer endpoints.
%
%   hh is nodal soil pressure head (NN x 1, cmH2O), converted to MPa by
%   1 cmH2O = 9.80665e-5 MPa.
%
%   In clm5/radial_darcy_series modes, soil-side conductance can be combined
%   with an optional CLM5-style root tissue-side conductance. This root-side
%   term is controlled by plantHydro.root_krmax [m/s], interpreted as root
%   conductance per root area per path length. If root_krmax is absent or
%   <= 0, root tissue is treated as non-limiting here. Whole-plant
%   conductance plantHydro.gpmax is deliberately not used in this soil-root
%   function because calSPAC already applies gpmax*vulnerability(psil) along
%   the plant hydraulic path.
%
%   Optional rootContact applies a piecewise-linear root-soil contact
%   reduction factor to endpoint conductance:
%       f_contact = 1,                                  psi >= psi_start
%                 = (psi - psi_end)/(psi_start-psi_end), psi_end < psi < psi_start
%                 = 0,                                  psi <= psi_end
%   Formula source: CONTACT mechanism in Lei et al. (2023, Journal of
%   Hydrology 617:128879), following the Couvreur et al. style of reducing
%   rhizosphere conductivity under root-soil contact loosening.
%   Default parameter source in config_vegetation.txt:
%       psi_start = -0.05 MPa, psi_end = -1.5 MPa
%   These are adapted from the sunflower case discussed by Lei et al. (2023)
%   and must be calibrated for the target species, soil, and experiment.

    rhow = 1000;          % water density (kg/m3)
    g = 9.80665;          % gravitational acceleration (m/s2)
    pa_per_mpa = 1e6;     % Pa/MPa
    if nargin < 9
        rootContact = [];
    end
    if nargin < 10
        plantHydro = struct();
    end
    rootContact = normalizeRootContact(rootContact);

    gsr_mode = 'katul';
    r_root_m = 5e-4;      % default fine-root radius: 0.5 mm
    if nargin >= 10 && isstruct(plantHydro)
        if isfield(plantHydro, 'gsr_formula') && ~isempty(plantHydro.gsr_formula)
            gsr_mode = char(plantHydro.gsr_formula);
        end
        if isfield(plantHydro, 'r_root_m') && plantHydro.r_root_m > 0
            r_root_m = plantHydro.r_root_m;
        end
    end

    NL = NN - 1;
    DeltZ = DeltZ(:);
    hh = hh(:);
    hhMPA = hh * 9.80665e-5;               % cmH2O -> MPa

    if numel(DeltZ) ~= NL
        error('PlantHydraulic:BadDeltZ', ...
            'DeltZ must have NN-1 elements. Got %d, expected %d.', numel(DeltZ), NL);
    end
    if numel(hh) ~= NN
        error('PlantHydraulic:BadSoilHead', ...
            'hh must have NN elements. Got %d, expected %d.', numel(hh), NN);
    end
    if ~isequal(size(Khh), [NL, 2])
        error('PlantHydraulic:BadKLhSize', ...
            'Khh/KL_h must be NL x 2. Got %d x %d, expected %d x 2.', ...
            size(Khh, 1), size(Khh, 2), NL);
    end
    if ~isreal(Khh) || any(~isfinite(Khh(:))) || any(Khh(:) < 0)
        error('PlantHydraulic:BadKLhValue', ...
            'Khh/KL_h must be real, finite, and non-negative before plant hydraulics.');
    end

    if isvector(RAI) && numel(RAI) == NL
        % Preserve layer-total root area: each endpoint carries half of the
        % layer support instead of smoothing root area onto adjacent layers.
        RAI_layer = RAI(:);
        RAI = [0.5 * RAI_layer, 0.5 * RAI_layer];
    end
    if ~isequal(size(RAI), [NL, 2])
        error('PlantHydraulic:BadRAISize', ...
            'RAI must be NL x 2, or NL x 1 for automatic endpoint expansion.');
    end
    if ~isreal(RAI) || any(~isfinite(RAI(:))) || any(RAI(:) < 0)
        error('PlantHydraulic:BadRAIValue', ...
            'RAI must be real, finite, and non-negative before plant hydraulics.');
    end
    [root_density_endpoint, root_area_endpoint] = rootGeometryEndpointFromPlantHydro( ...
        plantHydro, RAI, DeltZ, r_root_m, NL);

    psis2 = zeros(KT, 1);
    gsr = zeros(NN, 1);
    contact_factor_sum = 0;
    contact_weight_sum = 0;
    contact_factor_min = inf;
    contact_factor_max = -inf;
    contact_endpoint_count = 0;

    root_boundary_from_bottom = max(RL - Lmax, 0);
    Elmn_Lnth = 0;
    for ML = 1:NL
        Elmn_Lnth = Elmn_Lnth + DeltZ(ML);
        layer_top_from_bottom = Elmn_Lnth;

        % This model stores soil layers bottom-up. KL_h is passed as the
        % full soil profile, but only the active root zone can contribute to
        % soil-root conductance. root_area_endpoint is the root-presence mask;
        % the geometric test is retained as a guard against inconsistent input.
        has_root_overlap = layer_top_from_bottom > root_boundary_from_bottom;
        if ~has_root_overlap || all(root_area_endpoint(ML, :) <= 0)
            continue;
        end

        % Soil-root conductance options:
        %   katul: bulk root-zone length scale from Katul-style conductance.
        %   clm5 : CTSM-like soil-side and root-side conductances in series.
        %   radial_darcy_series: radial Darcy soil conductance plus root side.
        Zrj = Lmax / 100;                       % cm -> m
        depth_factor = Zrj * pi * rhow * g;     % Katul denominator
        if depth_factor <= 0 || ~isfinite(depth_factor)
            error('PlantHydraulic:BadLayerThickness', ...
                'Layer %d has invalid Lmax for soil-root conductance.', ML);
        end
        for ND = 1:2
            if root_area_endpoint(ML, ND) <= 0
                continue;
            end
            MN = ML + ND - 1;
            K_m_s = Khh(ML, ND) / 100;          % cm/s -> m/s
            f_contact = rootSoilContactFactor(hhMPA(MN), rootContact);
            K_eff_m_s = K_m_s * f_contact;
            m_water_per_MPa = pa_per_mpa / (rhow * g);

            switch gsr_mode
                case 'clm5'
                    % CLM5/CTSM structure:
                    %   r_soil = sqrt(1 / (pi * L_root))             [m]
                    %   k_s = K_soil / r_soil                        [1/s]
                    %   k_r includes root amount once through root area [1/s]
                    %   k_soil_root = k_s * k_r / (k_s + k_r)         [1/s]
                    %   gsr = k_soil_root * m_water_per_MPa           [m/s/MPa]
                    % CTSM puts RAI inside root_conductance, not as an
                    % additional multiplier after the series conductance.
                    L_root_density = max(root_density_endpoint(ML, ND), 1.0);
                    r_soil_m = sqrt(1.0 / (pi * L_root_density));
                    k_s = K_eff_m_s / r_soil_m;

                    [k_r, has_root_side] = rootSideConductance(plantHydro, root_area_endpoint, ...
                        root_density_endpoint, DeltZ, ML, ND, k_s, hhMPA(MN));
                    if has_root_side
                        k_combined = (k_s * k_r) / (k_s + k_r + 1e-30);
                    else
                        k_combined = k_s;
                    end
                    layer_endpoint_gsr = k_combined * m_water_per_MPa;

                case {'radial_darcy', 'radial_darcy_series'}
                    % Full radial Darcy soil-side conductance:
                    %   B       = measured root length density when available;
                    %             otherwise equivalent density from root area [m/m3]
                    %   r_outer = 1 / sqrt(pi*B)                      [m]
                    %   k_s     = K*A_root / (r_root*log(r_outer/r_root)) [1/s]
                    % Then use the same CLM5-style root-side series term.
                    B_m_per_m3 = root_density_endpoint(ML, ND);
                    B_m_per_m3 = max(B_m_per_m3, 1e-3);
                    r_outer_m = 1.0 / sqrt(pi * B_m_per_m3);
                    r_outer_m = max(r_outer_m, r_root_m * 1.01);
                    log_term = max(log(r_outer_m / r_root_m), 1e-6);
                    k_s = K_eff_m_s * root_area_endpoint(ML, ND) / (r_root_m * log_term);

                    [k_r, has_root_side] = rootSideConductance(plantHydro, root_area_endpoint, ...
                        root_density_endpoint, DeltZ, ML, ND, k_s, hhMPA(MN));
                    if has_root_side
                        k_combined = (k_s * k_r) / (k_s + k_r + 1e-30);
                    else
                        k_combined = k_s;
                    end
                    layer_endpoint_gsr = k_combined * m_water_per_MPa;

                otherwise   % 'katul'
                    layer_endpoint_gsr = K_eff_m_s * sqrt(root_area_endpoint(ML, ND)) / depth_factor * pa_per_mpa;
            end

            gsr(MN) = gsr(MN) + layer_endpoint_gsr;

            % Diagnostic only: the CONTACT factor itself is dimensionless.
            % Weighting by the uncorrected endpoint conductance contribution
            % reports how much the active root-soil hydraulic interface is
            % reduced by the Lei et al. (2023) CONTACT mechanism.
            uncorrected_gsr = K_m_s * sqrt(root_area_endpoint(ML, ND)) / depth_factor * pa_per_mpa;
            contact_factor_sum = contact_factor_sum + f_contact * uncorrected_gsr;
            contact_weight_sum = contact_weight_sum + uncorrected_gsr;
            contact_factor_min = min(contact_factor_min, f_contact);
            contact_factor_max = max(contact_factor_max, f_contact);
            contact_endpoint_count = contact_endpoint_count + 1;
        end
    end

    soil_gravity_MPa = soilNodeGravityMPa(DeltZ, rhow, g);
    hhMPA_effective = hhMPA - soil_gravity_MPa;

    sumgsr = sum(gsr);
    if ~isfinite(sumgsr)
        error('PlantHydraulic:NoRootConductance', ...
            'Total soil-root conductance is not finite. Check KL_h and RAI.');
    end
    if sumgsr <= 0
        % Extremely dry soil or complete root-soil contact loss. calSPAC
        % handles sumgsr=0 by setting actual transpiration to zero and using
        % the lower leaf-water-potential limit.
        psis2(KT) = mean(hhMPA_effective);
        contactStats = struct('enabled', rootContact.enabled, ...
            'weighted_mean', NaN, 'min', NaN, 'max', NaN, 'active_endpoint_count', 0);
        return;
    end

    psis2(KT) = sum(gsr .* hhMPA_effective) / sumgsr;

    if contact_endpoint_count > 0 && contact_weight_sum > 0
        contact_mean = contact_factor_sum / contact_weight_sum;
    else
        contact_mean = NaN;
        contact_factor_min = NaN;
        contact_factor_max = NaN;
    end
    contactStats = struct( ...
        'enabled', rootContact.enabled, ...
        'weighted_mean', contact_mean, ...
        'min', contact_factor_min, ...
        'max', contact_factor_max, ...
        'active_endpoint_count', contact_endpoint_count);
end

function [k_r, has_root_side] = rootSideConductance(plantHydro, root_area_endpoint, root_density_endpoint, DeltZ, ML, ND, k_s, psi_mpa)
%ROOTSIDECONDUCTANCE Root-side conductance for soil/root series coupling.
%   CLM5/CTSM computes root tissue conductance as
%       root_conductance = fs * RAI * krmax / (croot_length + z)
%   where fs is a root vulnerability factor, RAI is root area index, krmax
%   is root conductance per root area per path length, and z is soil depth.
%
%   STEMMUS-SCOPE legacy mode uses the empirical root radial resistance rrr:
%       rrr = Pr * (theta_s/theta) / (Rl * dz)
%   where Pr is the root radial impedance parameter [s m-1].
%   in the original water-head equation Trans=(PSIs-PSI)/(rsss+rrr+rxx).
%   Since PSIs and PSI are in meters and Trans is m/s, the resistance used
%   by that code acts as seconds in this layer-integrated relation, so its
%   inverse is a water-head conductance [1/s].
%   plantHydro.root_krmax follows that role and is separate from whole-plant
%   gpmax, which is applied later in calSPAC through the plant vulnerability
%   curve. This avoids counting whole-plant conductance once in gsr and
%   again in the plant hydraulic path.
    root_krmax = 0;
    root_path_offset_m = 0.25;
    root_side_mode = 'none';
    if isstruct(plantHydro)
        if isfield(plantHydro, 'root_side_mode') && ~isempty(plantHydro.root_side_mode)
            root_side_mode = lower(strtrim(char(plantHydro.root_side_mode)));
        end
        if isfield(plantHydro, 'root_krmax') && isfinite(plantHydro.root_krmax)
            root_krmax = plantHydro.root_krmax;
        elseif isfield(plantHydro, 'soil_root_krmax') && isfinite(plantHydro.soil_root_krmax)
            root_krmax = plantHydro.soil_root_krmax;
        end
        if isfield(plantHydro, 'root_path_offset_m') && isfinite(plantHydro.root_path_offset_m) ...
                && plantHydro.root_path_offset_m >= 0
            root_path_offset_m = plantHydro.root_path_offset_m;
        end
    end

    if any(strcmp(root_side_mode, {'stemmus_scope_rrr', 'stemmus_rrr', 'rrr'}))
        [k_r, has_root_side] = stemmusScopeRootRadialConductance( ...
            plantHydro, root_density_endpoint, DeltZ, ML, ND);
    elseif root_krmax > 0
        has_root_side = true;
        fs_root = rootVulnerabilityFactor(plantHydro, psi_mpa);
        MN = ML + ND - 1;
        depth_node_m = soilNodeDepthFromSurfaceM(DeltZ);
        path_length_m = root_path_offset_m + depth_node_m(MN);
        path_length_m = max(path_length_m, 1e-6);
        k_r = fs_root * root_area_endpoint(ML, ND) * root_krmax / path_length_m;
    else
        has_root_side = false;
        k_r = NaN;
    end
end

function [k_r, has_root_side] = stemmusScopeRootRadialConductance( ...
        plantHydro, root_density_endpoint, DeltZ, ML, ND)
%STEMMUSSCOPEROOTRADIALCONDUCTANCE Legacy STEMMUS-SCOPE rrr as root side.
%   Original STEMMUS-SCOPE calc_rsoil.m:
%       rrr = Pr * (theta_s/theta) / Rl / dz_m
%   where Pr is root radial impedance [s m-1], Rl is root length density
%   [m/m3], and dz_m is layer thickness.
%   Here each finite-element endpoint represents half a layer, so the
%   denominator uses B_endpoint * (0.5*dz_m). The inverse is a water-head
%   conductance [1/s], ready for series coupling with k_s.
    has_root_side = false;
    k_r = NaN;
    if ~isstruct(plantHydro) || ...
            ~isfield(plantHydro, 'soil_theta_layer') || ...
            ~isfield(plantHydro, 'soil_theta_s_layer')
        return;
    end
    Pr = 4e11;
    if isfield(plantHydro, 'root_radial_impedance_Pr') && ...
            isfinite(plantHydro.root_radial_impedance_Pr) && ...
            plantHydro.root_radial_impedance_Pr > 0
        Pr = plantHydro.root_radial_impedance_Pr;
    end

    theta = plantHydro.soil_theta_layer(:);
    theta_s = plantHydro.soil_theta_s_layer(:);
    NL = numel(DeltZ);
    if numel(theta) < NL || numel(theta_s) < NL
        return;
    end

    B_m_per_m3 = max(root_density_endpoint(ML, ND), 1e-12);
    z_support_m = 0.5 * DeltZ(ML) / 100;
    theta_eff = max(theta(ML), 1e-6);
    theta_s_eff = max(theta_s(ML), theta_eff);

    rrr_seconds = Pr * (theta_s_eff / theta_eff) / ...
        max(B_m_per_m3 * z_support_m, 1e-30);
    if isfinite(rrr_seconds) && rrr_seconds > 0
        k_r = 1 / rrr_seconds;
        has_root_side = true;
    end
end

function fs_root = rootVulnerabilityFactor(plantHydro, psi_mpa)
%ROOTVULNERABILITYFACTOR Optional CLM5-style root PLC factor.
%   If root-specific vulnerability parameters are absent, root tissue
%   conductance is not additionally downregulated here. Whole-plant/leaf
%   vulnerability is handled by calSPAC.
    fs_root = 1;
    if ~isstruct(plantHydro)
        return;
    end

    if isfield(plantHydro, 'root_p50') && isfield(plantHydro, 'root_aa') && ...
            isfinite(plantHydro.root_p50) && isfinite(plantHydro.root_aa) && ...
            plantHydro.root_p50 < 0 && plantHydro.root_aa > 0
        fs_root = vegetation.hydraulics.VulnerabilityCurve( ...
            psi_mpa, 1, plantHydro.root_p50, plantHydro.root_aa);
        fs_root = max(0, min(1, fs_root));
    end
end

function [root_density_endpoint, root_area_endpoint] = rootGeometryEndpointFromPlantHydro( ...
        plantHydro, RAI, DeltZ, r_root_m, NL)
%ROOTGEOMETRYENDPOINTFROMPLANTHYDRO Source-consistent root endpoint geometry.
%   Measured-mode runs carry plantHydro.root_length_density_m_m3 from the
%   root CSV and use that B directly, deriving endpoint root area from the
%   same B. Empirical/Jackson modes do not have measured B, so they use RAI
%   as endpoint root area and derive equivalent B from it.
    measured_B = [];
    if isstruct(plantHydro) && isfield(plantHydro, 'root_length_density_m_m3')
        measured_B = plantHydro.root_length_density_m_m3;
    end

    DeltZ_m = DeltZ(:) / 100;
    endpoint_support_m = 0.5 * DeltZ_m;
    endpoint_support = [endpoint_support_m, endpoint_support_m];

    if ~isempty(measured_B)
        measured_B = measured_B(:);
        if numel(measured_B) ~= NL || any(~isfinite(measured_B)) || any(measured_B < 0)
            error('PlantHydraulic:BadMeasuredRootDensity', ...
                'plantHydro.root_length_density_m_m3 must be finite, non-negative, and NL x 1.');
        end
        root_density_endpoint = [measured_B, measured_B];
        root_area_endpoint = 2 * pi * r_root_m * root_density_endpoint .* endpoint_support;
        return;
    end

    root_area_endpoint = RAI;
    root_density_endpoint = root_area_endpoint ./ ...
        max(2 * pi * r_root_m * endpoint_support, 1e-30);
end

function soil_gravity_MPa = soilNodeGravityMPa(DeltZ, rhow, g)
%SOILNODEGRAVITYMPA Gravitational potential from each soil node to surface.
%   Grid convention is bottom-up: node 1 is deepest, node NN is the surface.
%   The returned positive term is subtracted from soil matric potential when
%   water must move upward to the root collar near the surface.
    depth_node_m = soilNodeDepthFromSurfaceM(DeltZ);
    soil_gravity_MPa = rhow * g * depth_node_m / 1e6;
end

function depth_node_m = soilNodeDepthFromSurfaceM(DeltZ)
%SOILNODEDEPTHFROMSURFACEM Node depth for the bottom-up soil grid.
%   Node 1 is deepest and node NN is the soil surface.
    DeltZ = DeltZ(:);
    NL = numel(DeltZ);
    NN = NL + 1;
    depth_node_m = zeros(NN, 1);
    depth_from_surface_cm = flip(cumsum(flip(DeltZ)));
    depth_node_m(1:NL) = depth_from_surface_cm / 100;
end

function rootContact = normalizeRootContact(rootContact)
    if nargin < 1 || isempty(rootContact) || ~isstruct(rootContact)
        rootContact.enabled = false;
        rootContact.psi_start = -0.05;
        rootContact.psi_end = -1.5;
        return;
    end
    if ~isfield(rootContact, 'enabled'), rootContact.enabled = false; end
    if ~isfield(rootContact, 'psi_start'), rootContact.psi_start = -0.05; end
    if ~isfield(rootContact, 'psi_end'), rootContact.psi_end = -1.5; end

    if rootContact.enabled && ...
            (~isfinite(rootContact.psi_start) || ~isfinite(rootContact.psi_end) || ...
             rootContact.psi_start >= 0 || rootContact.psi_end >= 0 || ...
             rootContact.psi_start <= rootContact.psi_end)
        error('PlantHydraulic:BadRootContact', ...
            ['Root contact thresholds must be finite negative MPa values with ', ...
             'psi_start > psi_end.']);
    end
end

function f_contact = rootSoilContactFactor(psi_mpa, rootContact)
%ROOTSOILCONTACTFACTOR Piecewise root-soil contact reduction factor.
%   Formula source: Lei et al. (2023, Journal of Hydrology 617:128879),
%   CONTACT mechanism; rhizosphere conductivity decreases linearly after
%   root-soil contact loosening begins.
%   Parameter source: default thresholds in config_vegetation.txt follow the
%   sunflower scenario in Lei et al. (2023): -0.05 MPa onset and -1.5 MPa
%   complete contact loss, converted to the negative water-potential sign
%   convention used by this code.
    if ~rootContact.enabled
        f_contact = 1;
        return;
    end
    if ~isreal(psi_mpa) || ~isfinite(psi_mpa)
        error('PlantHydraulic:BadContactPotential', ...
            'Soil/root contact potential must be real and finite.');
    end

    if psi_mpa >= rootContact.psi_start
        f_contact = 1;
    elseif psi_mpa <= rootContact.psi_end
        f_contact = 0;
    else
        f_contact = (psi_mpa - rootContact.psi_end) / ...
                    (rootContact.psi_start - rootContact.psi_end);
    end
    f_contact = max(0, min(1, f_contact));
end
