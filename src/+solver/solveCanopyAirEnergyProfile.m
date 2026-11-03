function air = solveCanopyAirEnergyProfile(canopy, gap, resist_out, Ta_ref, Tch, Tcu1d, Ts, opts)
%SOLVECANOPYAIRENERGYPROFILE Implicit layered canopy-air heat balance.
%
% This solves a linear, layer-resolved sensible-heat balance for canopy air:
%
%   storage_i * (Ta_i - Ta_i_old)
%     = diffusion_i + leaf_heat_i + soil_heat_i
%
% The returned Ta_layer can be used inside an EB iteration without the
% diagnostic H/Kh feedback loop used by calc_canopy_air_profile.

arguments
    canopy struct
    gap struct
    resist_out struct
    Ta_ref double
    Tch double
    Tcu1d double
    Ts double
    opts.ra_leaf_h double {mustBePositive}
    opts.ra_soil_h double {mustBePositive}
    opts.ra_atm_h double = NaN
    opts.Ta_layer_old double = []
    opts.Fs double = [1 0]
    opts.Delt_t double {mustBePositive} = 3600
end

phys = equations.physicalConstants();
rhoa = phys.Pstd ./ (phys.R_da_SI .* (phys.Tr + phys.C_to_K));
cp = phys.Cp_air_SI;

nl = canopy.nlayers;
if isfield(canopy, 'layer_weight') && ~isempty(canopy.layer_weight)
    layer_weight = canopy.layer_weight(:);
else
    layer_weight = ones(nl, 1) / nl;
end
layer_weight = layer_weight ./ max(sum(layer_weight), eps);

if isfield(canopy, 'LAI_layer') && ~isempty(canopy.LAI_layer)
    LAI_layer = canopy.LAI_layer(:);
else
    LAI_layer = ones(nl, 1) * canopy.LAI / nl;
end

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

Ps = gap.Ps(1:nl);
Tleaf = (1 - Ps(:)) .* Tch(:) + Ps(:) .* Tcu1d(:);

g_leaf = LAI_layer .* rhoa .* cp ./ opts.ra_leaf_h;
Tsoil = opts.Fs(:)' * Ts(:);
g_soil = rhoa .* cp ./ opts.ra_soil_h;

storage = rhoa .* cp .* dz ./ opts.Delt_t;
Ta_layer_old = Ta_ref .* ones(nl, 1);
if ~isempty(opts.Ta_layer_old)
    old_candidate = opts.Ta_layer_old(:);
    if numel(old_candidate) == nl && all(isfinite(old_candidate))
        Ta_layer_old = old_candidate;
    end
end

A = zeros(nl, nl);
b = storage .* Ta_layer_old + g_leaf .* Tleaf;

for i = 1:nl
    A(i, i) = A(i, i) + storage(i) + g_leaf(i);

    if i == 1
        if isfinite(opts.ra_atm_h) && opts.ra_atm_h > 0
            g_top = rhoa .* cp ./ opts.ra_atm_h;
        else
            g_top = rhoa .* cp .* Kh_layer(i) ./ max(0.5 .* dz(i), 1e-6);
        end
        A(i, i) = A(i, i) + g_top;
        b(i) = b(i) + g_top .* Ta_ref;
    else
        K_face = 0.5 .* (Kh_layer(i - 1) + Kh_layer(i));
        dz_face = 0.5 .* (dz(i - 1) + dz(i));
        g_up = rhoa .* cp .* K_face ./ max(dz_face, 1e-6);
        A(i, i) = A(i, i) + g_up;
        A(i, i - 1) = A(i, i - 1) - g_up;
    end

    if i < nl
        K_face = 0.5 .* (Kh_layer(i) + Kh_layer(i + 1));
        dz_face = 0.5 .* (dz(i) + dz(i + 1));
        g_down = rhoa .* cp .* K_face ./ max(dz_face, 1e-6);
        A(i, i) = A(i, i) + g_down;
        A(i, i + 1) = A(i, i + 1) - g_down;
    else
        A(i, i) = A(i, i) + g_soil;
        b(i) = b(i) + g_soil .* Tsoil;
    end
end

Ta_layer = A \ b;

air.Ta_layer = Ta_layer;
air.Ta_layer_old = Ta_layer_old;
air.Kh_layer = Kh_layer;
air.z_layer = z_layer;
air.dz = dz;
air.Tleaf_layer = Tleaf;
air.g_leaf = g_leaf;
air.g_soil = g_soil;
air.g_top = g_top;
air.ra_atm_h = opts.ra_atm_h;
air.Tsoil = Tsoil;
air.method = 'implicit_layer_energy_v1';
end
