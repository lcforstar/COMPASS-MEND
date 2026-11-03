function [dz, dzTopDown, NL, ML] = buildVerticalGrid( ...
        legacyRequestedNL, totDepth, maxElements) %#ok<INUSD>
%BUILDVERTICALGRID  Build a depth-banded 1-D soil-column grid.
%
%   [dz, dzTopDown, NL, ML] = init.buildVerticalGrid( ...
%       legacyRequestedNL, totDepth, maxElements)
%   builds layer thicknesses from the surface downward by a deterministic
%   depth schedule (not an adaptive mesh). The first argument is retained
%   only for the legacy Dtrmn_Z call signature; it does not select the
%   number of layers. Layer count is fixed by accumulating thicknesses
%   until totDepth.
%
%   Inputs:
%     legacyRequestedNL  [-]  unused legacy layer-count request (was NL in
%                            site.toml / Dtrmn_Z). Depth banding owns NL.
%     totDepth           [cm] total soil-column depth, positive downward
%     maxElements        [-]  maximum number of layers allowed
%
%   Returns:
%     dz           [cm, NL x 1]  layer thickness, bottom-up; dz(1) is deepest
%                               and dz(NL) is the surface layer
%     dzTopDown    [cm, 1 x maxElements] top-down thicknesses (zeros after NL)
%     NL           [-] actual number of active layers
%     ML           [-] index of the final active top-down layer (= NL)
%
%   Depth schedule (from surface, z = top of the new layer) [cm]:
%     0-20          -> 1      surface coupling (stable cm-scale; was 0.25)
%     20-100        -> 5      shallow hydraulic response
%     100-200       -> 10     upper root zone
%     200-500       -> 25     0.2-5 m (covers a 5 m production column)
%     500-1000      -> 50     5-10 m
%     1000-2500     -> 75     10-25 m deep-root band (~0.75 m)
%     2500-5000     -> 125    25-50 m
%     5000-10000    -> 250    50-100 m
%     10000-50000   -> 500    100-500 m fixed deep-vadose resolution (5 m)
%     z >= 50000    -> 1000   >500 m fixed deep-vadose resolution (10 m)
%
%   Resolution follows multiphysics time scales (water, vapor, energy, gas,
%   solutes, roots), not pure aquifer storage. Deep layers use a fixed
%   deep-vadose resolution strategy (not a runtime adaptive algorithm).
%
%   The last layer is shortened so sum(dz) == totDepth exactly.
%
%   Design notes:
%     - Static depth banding only. Wetting / freeze fronts and water-table
%       motion are NOT remeshed. Adaptive mesh would force checkpoint remap
%       and conservation re-projection; deferred past the first production
%       grid contract. Prefer static grids for future global runs.
%     - Depth, not legacyRequestedNL, controls resolution.
%
%   Legacy provenance:
%     Replaces the former fixed-index STEMMUS/COMPASS Dtrmn_Z table
%     (0.25/0.5/1/2/.../100 cm bands).
%
%   Editor: Jizixing Zhu; Date: 2026-06-26.
%   Change: 2026-08-10 Jizixing Zhu -- depth-banded grid: 1 cm surface,
%           mid-resolution through 25 m roots, deep fixed resolution below
%           100 m; rename unused NL arg to legacyRequestedNL; static mesh.
%   Change: 2026-08-10 Jizixing Zhu -- deep band 100-500 m at 5 m and
%           >500 m at 10 m (not a single 15 m storage-only spacing).

    if nargin < 3 || isempty(maxElements)
        maxElements = 5000;
    end
    if totDepth <= 0
        error('buildVerticalGrid:badDepth', 'totDepth must be positive [cm].');
    end
    if maxElements < 1
        error('buildVerticalGrid:badMaxElements', ...
            'maxElements must be positive.');
    end

    dzTopDown = zeros(1, maxElements);
    cumulativeDepth = 0;     % [cm] depth of the bottom of the last filled layer
    ML = 0;

    while cumulativeDepth < totDepth - max(1e-12, 10 * eps(totDepth))
        ML = ML + 1;
        if ML > maxElements
            error('buildVerticalGrid:depthNotCovered', ...
                ['Total depth %.3f cm was not fully covered by %d elements. ', ...
                 'Increase maxElements or coarsen the deep fixed resolution.'], ...
                totDepth, maxElements);
        end
        layerTop_cm = cumulativeDepth;
        dzLayer = localLayerThickness_cm(layerTop_cm);
        remaining = totDepth - cumulativeDepth;
        if dzLayer > remaining
            dzLayer = remaining;
        end
        if dzLayer <= 0
            error('buildVerticalGrid:nonPositiveLayer', ...
                'Encountered a non-positive layer thickness at depth %.12g cm.', ...
                layerTop_cm);
        end
        dzTopDown(ML) = dzLayer;
        cumulativeDepth = cumulativeDepth + dzLayer;
    end

    NL = ML;
    dzFlipped = flip(dzTopDown(1:NL));
    dz = dzFlipped(dzFlipped > 0).';
    NL = numel(dz);
    ML = NL;

    if abs(sum(dz) - totDepth) > max(1e-9, eps(totDepth) * 10)
        error('buildVerticalGrid:depthMismatch', ...
            'Layer thickness sum %.12g cm does not match totDepth %.12g cm.', ...
            sum(dz), totDepth);
    end
end

% -----------------------------------------------------------------------------
function d = localLayerThickness_cm(zTop_cm)
%LOCALAYERTHICKNESSCM  Layer thickness for a new layer whose top is zTop_cm.
%
%   zTop_cm is the depth below the surface of the top face of the new layer.
%   Thickness increases with depth. Below 100 m a fixed deep-vadose resolution
%   is used (5 m to 500 m depth, 10 m below 500 m) so multiphysics gradients
%   remain representable without an adaptive mesh.

    if zTop_cm < 20
        d = 1;           % [cm] surface coupling zone
    elseif zTop_cm < 100
        d = 5;           % [cm] shallow hydraulic response
    elseif zTop_cm < 200
        d = 10;          % [cm] upper root zone
    elseif zTop_cm < 500
        d = 25;          % [cm] 0.2-5 m
    elseif zTop_cm < 1000
        d = 50;          % [cm] 5-10 m
    elseif zTop_cm < 2500
        d = 75;          % [cm] 10-25 m deep-root band
    elseif zTop_cm < 5000
        d = 125;         % [cm] 25-50 m
    elseif zTop_cm < 10000
        d = 250;         % [cm] 50-100 m
    elseif zTop_cm < 50000
        d = 500;         % [cm] 100-500 m fixed deep-vadose resolution (5 m)
    else
        d = 1000;        % [cm] >500 m fixed deep-vadose resolution (10 m)
    end
end
