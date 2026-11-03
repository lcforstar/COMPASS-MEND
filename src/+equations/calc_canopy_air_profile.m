function air = calc_canopy_air_profile(canopy, gap, resist_out, Ta_ref, ea_ref, CO2_ref, ...
        Tch, Tcu1d, ech, ecu1d, Cch, Ccu1d, flux_gradient)
%CALC_CANOPY_AIR_PROFILE Diagnostic layer-resolved canopy air profile.
%
% Phase 4 canopy-air layer:
%   This function estimates Ta_layer, ea_layer, CO2_layer, and Kh_layer.
%   If flux_gradient is provided, Ta_layer and ea_layer are solved from a
%   steady flux-gradient closure using layer heat/water sources and Kh_layer.
%   CO2_layer remains diagnostic until the biochemical feedback is explicitly
%   connected.
%
% Layer order follows SCOPE:
%   layer 1  = canopy top
%   layer nl = canopy bottom
%
% 改进者 / Modified by:   朱寄子星 (Jizixing Zhu)
% 改进日期 / Modified date: 2026-04-19

nl = canopy.nlayers;
if nargin < 13
    flux_gradient = [];
end

if isfield(canopy, 'layer_weight') && ~isempty(canopy.layer_weight)
    layer_weight = canopy.layer_weight(:);
else
    layer_weight = ones(nl, 1) / nl;
end
layer_weight = layer_weight ./ max(sum(layer_weight), eps);

depth_geom = ((1:nl)' - 0.5) ./ nl;                 % [-] 0=top, 1=bottom
depth_lai  = cumsum(layer_weight) - 0.5 .* layer_weight; % [-] cumulative LAI depth
if isfield(canopy, 'hc') && ~isempty(canopy.hc)
    hc = max(canopy.hc, 1e-6);
    z_layer_m = hc .* (1 - depth_geom);             % [m] height above ground
else
    hc = 1;
    z_layer_m = 1 - depth_geom;
end

Ps = gap.Ps(1:nl);
Tleaf_layer = (1 - Ps) .* Tch(:)   + Ps .* Tcu1d(:);
eleaf_layer = (1 - Ps) .* ech(:)   + Ps .* ecu1d(:);
Cleaf_layer = (1 - Ps) .* Cch(:)   + Ps .* Ccu1d(:);

Kh_top = 0.05;
if isfield(resist_out, 'Kh') && ~isempty(resist_out.Kh) && isfinite(resist_out.Kh)
    Kh_top = max(resist_out.Kh, 1e-5);
end
n_canopy = 2.0;
if isfield(resist_out, 'n_canopy') && ~isempty(resist_out.n_canopy) && isfinite(resist_out.n_canopy)
    n_canopy = max(resist_out.n_canopy, 1e-6);
end

% W&V-style within-canopy turbulent diffusivity decay. Kh_top is the
% top/roughness-layer diffusivity already calculated by resistances_clm45.
Kh_layer = max(1e-5, Kh_top .* exp(-n_canopy .* depth_lai));

relax_strength = 1.5;
mix = 1 - exp(-relax_strength .* depth_lai);

if isstruct(flux_gradient) && isfield(flux_gradient, 'H_layer_ground')
    phys = equations.physicalConstants();
    rhoa = phys.Pstd ./ (phys.R_da_SI .* (phys.Tr + phys.C_to_K));
    cp   = phys.Cp_air_SI;

    dz_layer = hc ./ nl .* ones(nl, 1);  % geometric layer thickness for diffusion

    H_layer = flux_gradient.H_layer_ground(:);
    LE_layer = flux_gradient.LE_layer_ground(:);
    H_soil = 0;
    LE_soil = 0;
    if isfield(flux_gradient, 'H_soil_ground'), H_soil = flux_gradient.H_soil_ground; end
    if isfield(flux_gradient, 'LE_soil_ground'), LE_soil = flux_gradient.LE_soil_ground; end

    e_to_q = [];
    if isfield(flux_gradient, 'e_to_q'), e_to_q = flux_gradient.e_to_q; end
    lambda_ref = [];
    if isfield(flux_gradient, 'lambda'), lambda_ref = flux_gradient.lambda; end
    if isempty(e_to_q) || ~isfinite(e_to_q) || e_to_q <= 0
        e_to_q = 0.0062;
    end
    if isempty(lambda_ref) || ~isfinite(lambda_ref) || lambda_ref <= 0
        lambda_ref = (2.501 - 0.002361 .* Ta_ref) .* 1e6;
    end

    E_layer = LE_layer ./ lambda_ref;    % [kg m-2 s-1] layer water source
    E_soil  = LE_soil ./ lambda_ref;     % [kg m-2 s-1] soil water source

    Ta_layer = Ta_ref .* ones(nl, 1);
    q_layer  = ea_ref .* e_to_q .* ones(nl, 1);

    int_T = 0;
    int_q = 0;
    for il = 1:nl
        H_mid = H_soil + sum(H_layer(il + 1:end)) + 0.5 .* H_layer(il);
        E_mid = E_soil + sum(E_layer(il + 1:end)) + 0.5 .* E_layer(il);
        if il > 1
            H_prev = H_soil + sum(H_layer(il:end)) + 0.5 .* H_layer(il - 1);
            E_prev = E_soil + sum(E_layer(il:end)) + 0.5 .* E_layer(il - 1);
            int_T = int_T + dz_layer(il - 1) .* H_prev ./ (rhoa .* cp .* Kh_layer(il - 1));
            int_q = int_q + dz_layer(il - 1) .* E_prev ./ (rhoa .* Kh_layer(il - 1));
        end
        Ta_layer(il) = Ta_ref + int_T + 0.5 .* dz_layer(il) .* H_mid ./ (rhoa .* cp .* Kh_layer(il));
        q_layer(il)  = ea_ref .* e_to_q + int_q + 0.5 .* dz_layer(il) .* E_mid ./ (rhoa .* Kh_layer(il));
    end

    ea_layer = q_layer ./ e_to_q;

    % CO2 flux-gradient: integrate from canopy top using layer assimilation and Kh_layer.
    % [改进 2026-04-19 朱寄子星] 用逐层净光合通量代替诊断松弛公式
    M_air = 0.02897;               % [kg mol-1] dry air molar mass
    rho_mol = rhoa ./ M_air;       % [mol m-3] molar air density
    A_layer_g = zeros(nl, 1);
    if isfield(flux_gradient, 'A_layer_ground')
        A_layer_g = max(0, flux_gradient.A_layer_ground(:));
    end
    A_soil_CO2 = 0;
    if isfield(flux_gradient, 'A_soil_ground')
        A_soil_CO2 = max(0, flux_gradient.A_soil_ground);
    end
    CO2_layer = CO2_ref .* ones(nl, 1);
    int_C = 0;
    for il = 1:nl
        A_mid = A_soil_CO2 + sum(A_layer_g(il+1:end)) + 0.5 .* A_layer_g(il);
        if il > 1
            A_prev = A_soil_CO2 + sum(A_layer_g(il:end)) + 0.5 .* A_layer_g(il-1);
            int_C = int_C + dz_layer(il-1) .* A_prev ./ (Kh_layer(il-1) .* rho_mol);
        end
        CO2_layer(il) = CO2_ref - int_C - 0.5 .* dz_layer(il) .* A_mid ./ (Kh_layer(il) .* rho_mol);
    end
    CO2_layer = max(0, CO2_layer);
else
    % Fallback diagnostic relaxation when no flux sources are available.
    Ta_layer  = Ta_ref + mix .* (Tleaf_layer - Ta_ref);
    ea_layer  = ea_ref + mix .* (eleaf_layer - ea_ref);
    CO2_layer = CO2_ref + mix .* (Cleaf_layer - CO2_ref);
end

air.Ta_layer       = Ta_layer;      % [degC] layer diagnostic canopy-air temperature
air.ea_layer       = max(0, ea_layer); % [hPa] layer diagnostic vapour pressure
air.CO2_layer      = max(0, CO2_layer); % [umol mol-1] layer diagnostic CO2
air.Kh_layer       = Kh_layer;      % [m2 s-1] layer turbulent diffusivity
air.z_layer_m      = z_layer_m;     % [m] layer center height above ground
air.depth_geom     = depth_geom;    % [-] geometric depth, 0=top, 1=bottom
air.depth_lai      = depth_lai;     % [-] cumulative LAI-weighted depth
air.layer_weight   = layer_weight;  % [-] layer LAI fraction
air.mix_air_source = mix;           % [-] diagnostic relaxation strength per layer
air.Tleaf_layer    = Tleaf_layer;   % [degC] two-leaf weighted leaf source
air.eleaf_layer    = eleaf_layer;   % [hPa] two-leaf weighted leaf-source vapour pressure
air.Cleaf_layer    = Cleaf_layer;   % [umol mol-1] two-leaf weighted leaf-source CO2
if isstruct(flux_gradient) && isfield(flux_gradient, 'H_layer_ground')
    air.method = 'flux_gradient_CO2_v1';
else
    air.method = 'diagnostic_relaxation_v1';
end
end
