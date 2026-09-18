function air = solveCanopyAirCO2Profile(canopy, resist_out, CO2_ref, A_layer_ground, opts)
%SOLVECANOPYAIRCO2PROFILE Implicit layered canopy-air CO2 balance.
%
% CO2 concentration is carried as umol mol-1 (ppm). The matrix conductances
% are mol_air m-2 s-1, so multiplying by ppm gives umol_CO2 m-2 s-1.
% A_layer_ground is therefore used directly as a canopy uptake sink
% [umol_CO2 m-2 ground s-1], without an additional 1e-6 factor.

arguments
    canopy struct
    resist_out struct
    CO2_ref double
    A_layer_ground double
    opts.ra_atm_c double = NaN
    opts.CO2_layer_old double = []
    opts.A_soil_ground double = 0
    opts.Delt_t double {mustBePositive} = 3600
end

phys = equations.physicalConstants();
rhoa = phys.Pstd ./ (phys.R_da_SI .* (phys.Tr + phys.C_to_K));
rho_mol = rhoa ./ (phys.MW_da ./ 1000); % [mol m-3], MW_da is kg kmol-1 numerically

nl = canopy.nlayers;
if isfield(canopy, 'layer_weight') && ~isempty(canopy.layer_weight)
    layer_weight = canopy.layer_weight(:);
else
    layer_weight = ones(nl, 1) / nl;
end
layer_weight = layer_weight ./ max(sum(layer_weight), eps);

if isfield(canopy, 'hc') && ~isempty(canopy.hc)
    hc = max(canopy.hc, 1e-6);
else
    hc = 1;
end
z_face = linspace(hc, 0, nl + 1).';
z_layer = 0.5 .* (z_face(1:end - 1) + z_face(2:end));
dz = abs(diff(z_face));
depth_lai = cumsum(layer_weight) - 0.5 .* layer_weight;

Kh_top = 0.05;
if isfield(resist_out, 'Kh') && ~isempty(resist_out.Kh) && isfinite(resist_out.Kh)
    Kh_top = max(resist_out.Kh, 1e-5);
end
n_canopy = 2.0;
if isfield(resist_out, 'n_canopy') && ~isempty(resist_out.n_canopy) && isfinite(resist_out.n_canopy)
    n_canopy = max(resist_out.n_canopy, 1e-6);
end
Kh_layer = max(1e-5, Kh_top .* exp(-n_canopy .* depth_lai));

A_layer = A_layer_ground(:);
if numel(A_layer) ~= nl || any(~isfinite(A_layer))
    A_layer = zeros(nl, 1);
end

storage = rho_mol .* dz ./ opts.Delt_t;
CO2_old = CO2_ref .* ones(nl, 1);
if ~isempty(opts.CO2_layer_old)
    old_candidate = opts.CO2_layer_old(:);
    if numel(old_candidate) == nl && all(isfinite(old_candidate))
        CO2_old = old_candidate;
    end
end

A = zeros(nl, nl);
b = storage .* CO2_old - A_layer;

for i = 1:nl
    A(i, i) = A(i, i) + storage(i);

    if i == 1
        if isfinite(opts.ra_atm_c) && opts.ra_atm_c > 0
            g_top = rho_mol ./ opts.ra_atm_c;
        else
            g_top = rho_mol .* Kh_layer(i) ./ max(0.5 .* dz(i), 1e-6);
        end
        A(i, i) = A(i, i) + g_top;
        b(i) = b(i) + g_top .* CO2_ref;
    else
        K_face = 0.5 .* (Kh_layer(i - 1) + Kh_layer(i));
        dz_face = 0.5 .* (dz(i - 1) + dz(i));
        g_up = rho_mol .* K_face ./ max(dz_face, 1e-6);
        A(i, i) = A(i, i) + g_up;
        A(i, i - 1) = A(i, i - 1) - g_up;
    end

    if i < nl
        K_face = 0.5 .* (Kh_layer(i) + Kh_layer(i + 1));
        dz_face = 0.5 .* (dz(i) + dz(i + 1));
        g_down = rho_mol .* K_face ./ max(dz_face, 1e-6);
        A(i, i) = A(i, i) + g_down;
        A(i, i + 1) = A(i, i + 1) - g_down;
    else
        b(i) = b(i) + opts.A_soil_ground;
    end
end

CO2_layer = max(0, A \ b);

air.CO2_layer = CO2_layer;
air.CO2_layer_old = CO2_old;
air.Kh_layer = Kh_layer;
air.z_layer = z_layer;
air.dz = dz;
air.g_top = g_top;
air.A_layer_ground = A_layer;
air.A_soil_ground = opts.A_soil_ground;
air.method = 'implicit_layer_co2_v1';
end
