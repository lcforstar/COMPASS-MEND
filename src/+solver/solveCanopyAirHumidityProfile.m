function air = solveCanopyAirHumidityProfile(canopy, gap, resist_out, ea_ref, Tch, Tcu1d, Ts, opts)
%SOLVECANOPYAIRHUMIDITYPROFILE Implicit layered canopy-air humidity balance.
%
% The unknown is specific humidity q; ea_layer is returned for compatibility
% with the existing SCOPE heat-flux routines. PSI_leaf and PSI_soil are
% water potentials in J kg-1, matching heatfluxes_clm45.

arguments
    canopy struct
    gap struct
    resist_out struct
    ea_ref double
    Tch double
    Tcu1d double
    Ts double
    opts.ra_leaf_w double {mustBePositive}
    opts.ra_soil_w double {mustBePositive}
    opts.ra_atm_w double = NaN
    opts.rcw_sha double = []
    opts.rcw_sun double = []
    opts.rss double = 0
    opts.ea_layer_old double = []
    opts.e_to_q double {mustBePositive} = 0.0062
    opts.PSI_leaf double = 0
    opts.PSI_soil double = 0
    opts.Fs double = [1 0]
    opts.Delt_t double {mustBePositive} = 3600
end

phys = equations.physicalConstants();
rhoa = phys.Pstd ./ (phys.R_da_SI .* (phys.Tr + phys.C_to_K));

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
rcw_sha = normalize_layer_resistance(opts.rcw_sha, nl, 0);
rcw_sun = normalize_layer_resistance(opts.rcw_sun, nl, 0);

e_leaf_sha = vapor_pressure_with_water_potential(Tch(:), opts.PSI_leaf);
e_leaf_sun = vapor_pressure_with_water_potential(Tcu1d(:), opts.PSI_leaf);
q_leaf_sha = e_leaf_sha .* opts.e_to_q;
q_leaf_sun = e_leaf_sun .* opts.e_to_q;

g_leaf_sha = LAI_layer .* (1 - Ps(:)) .* rhoa ./ max(opts.ra_leaf_w + rcw_sha, 1e-6);
g_leaf_sun = LAI_layer .* Ps(:) .* rhoa ./ max(opts.ra_leaf_w + rcw_sun, 1e-6);

Tsoil = opts.Fs(:)' * Ts(:);
e_soil = vapor_pressure_with_water_potential(Tsoil, opts.PSI_soil);
q_soil = e_soil .* opts.e_to_q;

storage = rhoa .* dz ./ opts.Delt_t;
q_old = ea_ref .* opts.e_to_q .* ones(nl, 1);
if ~isempty(opts.ea_layer_old)
    old_candidate = opts.ea_layer_old(:);
    if numel(old_candidate) == nl && all(isfinite(old_candidate))
        q_old = max(0, old_candidate) .* opts.e_to_q;
    end
end

if q_old(nl) > q_soil
    % Condensation: vapor moves from canopy air to the soil surface, so the
    % internal soil surface resistance should not limit the air-side sink.
    g_soil = rhoa ./ max(opts.ra_soil_w, 1e-6);
    soil_vapor_mode = 'condensation_air_to_soil';
else
    % Evaporation: vapor moves from soil pores/surface to canopy air and is
    % limited by both aerodynamic and soil surface resistance.
    g_soil = rhoa ./ max(opts.ra_soil_w + opts.rss, 1e-6);
    soil_vapor_mode = 'evaporation_soil_to_air';
end

A = zeros(nl, nl);
b = storage .* q_old + g_leaf_sha .* q_leaf_sha + g_leaf_sun .* q_leaf_sun;

for i = 1:nl
    A(i, i) = A(i, i) + storage(i) + g_leaf_sha(i) + g_leaf_sun(i);

    if i == 1
        if isfinite(opts.ra_atm_w) && opts.ra_atm_w > 0
            g_top = rhoa ./ opts.ra_atm_w;
        else
            g_top = rhoa .* Kh_layer(i) ./ max(0.5 .* dz(i), 1e-6);
        end
        A(i, i) = A(i, i) + g_top;
        b(i) = b(i) + g_top .* ea_ref .* opts.e_to_q;
    else
        K_face = 0.5 .* (Kh_layer(i - 1) + Kh_layer(i));
        dz_face = 0.5 .* (dz(i - 1) + dz(i));
        g_up = rhoa .* K_face ./ max(dz_face, 1e-6);
        A(i, i) = A(i, i) + g_up;
        A(i, i - 1) = A(i, i - 1) - g_up;
    end

    if i < nl
        K_face = 0.5 .* (Kh_layer(i) + Kh_layer(i + 1));
        dz_face = 0.5 .* (dz(i) + dz(i + 1));
        g_down = rhoa .* K_face ./ max(dz_face, 1e-6);
        A(i, i) = A(i, i) + g_down;
        A(i, i + 1) = A(i, i + 1) - g_down;
    else
        A(i, i) = A(i, i) + g_soil;
        b(i) = b(i) + g_soil .* q_soil;
    end
end

q_layer = max(0, A \ b);
ea_layer = q_layer ./ opts.e_to_q;

air.ea_layer = ea_layer;
air.q_layer = q_layer;
air.ea_layer_old = q_old ./ opts.e_to_q;
air.Kh_layer = Kh_layer;
air.z_layer = z_layer;
air.dz = dz;
air.g_leaf_sha = g_leaf_sha;
air.g_leaf_sun = g_leaf_sun;
air.g_soil = g_soil;
air.g_top = g_top;
air.eleaf_sha = e_leaf_sha;
air.eleaf_sun = e_leaf_sun;
air.esoil = e_soil;
air.soil_vapor_mode = soil_vapor_mode;
air.method = 'implicit_layer_humidity_v1';
end

function r = normalize_layer_resistance(r_in, nl, fallback)
    if isempty(r_in)
        r = fallback .* ones(nl, 1);
        return;
    end
    r = r_in(:);
    if numel(r) ~= nl || any(~isfinite(r))
        r = fallback .* ones(nl, 1);
    end
    r = max(r, 0);
end

function e_hPa = vapor_pressure_with_water_potential(T_C, PSI_Jkg)
    phys = equations.physicalConstants();
    es_Pa = equations.Calc_Esat_Flatau(T_C, 101325);   % 只需 es
    es_hPa = reshape(es_Pa, size(T_C)) .* 0.01;
    mw_wv_kg_mol = phys.MW_wv ./ 1000;
    e_hPa = es_hPa .* exp(PSI_Jkg .* mw_wv_kg_mol ./ phys.R ./ (T_C + 273.15));
    e_hPa = max(0, e_hPa);
end
