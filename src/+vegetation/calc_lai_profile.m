function [LAI_layer, layer_weight, shape, z] = calc_lai_profile(LAI, nl, mode, peak, sigma, skew, k)
%calc_lai_profile Build a canopy-top to canopy-bottom LAI profile.
%
% Layer order follows SCOPE canopy convention:
%   layer 1  = canopy top
%   layer nl = canopy bottom
%
% LAI_layer always sums to LAI. layer_weight is independent of total LAI
% and can be reused when LAI changes in the time series.

if nargin < 3 || isempty(mode), mode = 'uniform'; end
if nargin < 4 || isempty(peak), peak = 0.50; end
if nargin < 5 || isempty(sigma), sigma = 0.24; end
if nargin < 6 || isempty(skew), skew = 0.8; end
if nargin < 7 || isempty(k), k = 1.5; end

if ~isscalar(nl) || nl < 1 || ~isfinite(nl)
    error('Plant:BadLAIProfileLayers', 'nl must be a positive finite scalar.');
end
nl = round(nl);

if ~isscalar(LAI) || ~isfinite(LAI) || LAI < 0
    error('Plant:BadLAIProfileLAI', 'LAI must be a finite non-negative scalar.');
end

peak = min(max(peak, 0.01), 0.99);
sigma = max(sigma, 1e-3);
k = max(k, 0);

z = ((1:nl)' - 0.5) / nl;       % 0 = canopy top, 1 = canopy bottom
zz = min(max(z, 1e-3), 1 - 1e-3);

switch lower(strtrim(mode))
    case 'uniform'
        shape = ones(nl, 1);

    case 'pine_cone'
        % Conifer-like conical crown: small LAI at the pointed top,
        % monotonically increasing toward lower crown. Peak of z^k*exp(-0.5z)
        % lies at z=2k (outside [0,1] for k>=0.5), so no bottom taper occurs.
        shape = z.^k .* exp(-0.5 * z);

    case 'fruit_ellipsoid'
        % Fruit-tree / ellipsoid crown: foliage concentrated around the
        % middle crown, less at top and bottom.
        shape = exp(-((z - peak).^2) ./ (2 * sigma^2));

    case 'willow_drooping'
        % Drooping willow-like crown: middle-upper peak with a downward tail.
        base = exp(-((z - peak).^2) ./ (2 * sigma^2));
        shape = base .* max(1 + skew * z, 1e-3);

    case 'broadleaf_tree'
        % Generic broadleaf tree: Beta-distribution crown, naturally tapers
        % to zero at both canopy top and bottom.
        % peak sets the depth of maximum LAI (0=top, 1=bottom).
        % k controls concentration: higher k = narrower peak.
        a = 1 + peak * k;
        b = 1 + (1 - peak) * k;
        shape = zz.^(a - 1) .* (1 - zz).^(b - 1);

    case 'crop_upper'
        % Dense upper canopy, common first approximation for crops/grass.
        shape = exp(-k * z);

    case 'maize'
        % Tall crop with Gaussian LAI peak in the upper-middle canopy.
        % Defaults (peak=0.40, sigma=0.22) represent typical maize structure
        % and intentionally differ from fruit_ellipsoid (peak=0.50).
        % To override, use fruit_ellipsoid mode with custom peak/sigma.
        shape = exp(-((z - 0.40).^2) ./ (2 * 0.22^2));

    case 'shrub_low'
        % Shrub or low crown with more foliage near lower canopy.
        shape = exp(k * z);

    case 'orchard_open_center'
        % Managed orchard/open-center crown: ellipsoid-like but slightly
        % weighted toward middle-lower layers.
        base = exp(-((z - peak).^2) ./ (2 * sigma^2));
        shape = base .* max(0.7 + skew * z, 1e-3);

    case 'tropical_rainforest'
        % Tropical evergreen rainforest: bimodal distribution with a dominant
        % main-canopy peak (~30% depth) and a weaker understory peak (~70%).
        % Ref: Asner et al. (2002), Weiss et al. (2004).
        shape = exp(-((z - 0.30).^2) ./ (2 * 0.15^2)) ...
              + 0.40 .* exp(-((z - 0.70).^2) ./ (2 * 0.12^2));

    case 'boreal_needleleaf'
        % Boreal conifer (spruce, fir, larch): cylindrical crown, near-uniform
        % with a broad upper-middle concentration. Distinct from pine_cone
        % (conical, monotonically increasing). Uses broad Gaussian.
        shape = exp(-((z - 0.35).^2) ./ (2 * 0.32^2));

    case 'c4_grass'
        % C4 grass / savanna grass: very top-concentrated thin canopy.
        % Steeper exponential than crop_upper (hardcoded k=3.5).
        % For custom k, use crop_upper mode with LAI_profile_k.
        shape = exp(-3.5 .* z);

    case 'tundra'
        % Arctic tundra / dwarf shrub: very short canopy, near-uniform with
        % a gentle bottom-weight (older stems near ground bear more biomass).
        shape = 0.85 + 0.30 .* z;

    otherwise
        error('Plant:UnknownLAIProfileMode', ...
            ['Unknown LAI_profile_mode "%s". Valid modes: uniform, ', ...
             'pine_cone, boreal_needleleaf, tropical_rainforest, ', ...
             'broadleaf_tree, fruit_ellipsoid, willow_drooping, ', ...
             'orchard_open_center, maize, crop_upper, c4_grass, ', ...
             'shrub_low, tundra.'], mode);
end

shape(~isfinite(shape)) = 0;
shape = max(shape(:), 1e-9);
layer_weight = shape / sum(shape);
LAI_layer = LAI * layer_weight;
end
