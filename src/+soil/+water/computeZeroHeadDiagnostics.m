function diag = computeZeroHeadDiagnostics(h, nodeDepthCm)
%COMPUTEZEROHEADDIAGNOSTICS  Locate all saturated-unsaturated interfaces.
%
%   diag = soil.water.computeZeroHeadDiagnostics(h, nodeDepthCm) derives
%   write-only diagnostics from an accepted pressure-head profile. It never
%   feeds the Richards solve or changes the active set.
%
%   Inputs:
%     h           [in] [cm, NN x nCell] signed pressure head. h >= 0 is
%                 saturated and h < 0 is unsaturated.
%     nodeDepthCm [in] [cm, NN x 1] node depth below the surface in the same
%                 bottom-up order as h.
%
%   Returns:
%     diag [out] struct with:
%       node_saturated [logical, NN x nCell]
%       crossing_depth_cm [cm, NL x nCell], shallowest first, NaN padded
%       crossing_count [-, 1 x nCell]
%       shallowest_crossing_depth_cm [cm, 1 x nCell]
%       saturated_thickness_cm [cm, 1 x nCell]
%       saturated_fraction [-, 1 x nCell]
%
%   Method:
%     Within a sign-changing linear element, solve h(z)=0 by linear
%     interpolation. Exact zero-head nodes are retained once. Multiple
%     crossings are allowed so perched saturated zones are not collapsed into
%     one assumed groundwater table.
%
%   Editor: Jizixing Zhu; Date: 2026-07-12.

    if isvector(h), h = h(:); end
    depth = nodeDepthCm(:);
    NN = size(h, 1);
    nCell = size(h, 2);
    NL = NN - 1;
    assert(numel(depth) == NN && all(isfinite(depth)), ...
        'computeZeroHeadDiagnostics:depth', ...
        'nodeDepthCm must contain NN finite values.');
    assert(all(isfinite(h(:))), ...
        'computeZeroHeadDiagnostics:head', 'h must be finite.');
    assert(all(diff(depth) ~= 0) && ...
        (all(diff(depth) > 0) || all(diff(depth) < 0)), ...
        'computeZeroHeadDiagnostics:depthOrder', ...
        'nodeDepthCm must be strictly monotonic.');

    crossingDepth = NaN(NL, nCell);
    crossingCount = zeros(1, nCell);
    shallowest = NaN(1, nCell);
    saturatedThickness = zeros(1, nCell);
    profileThickness = abs(depth(end) - depth(1));

    for c = 1:nCell
        hc = h(:, c);
        depths = zeros(0, 1);
        for i = 1:NL
            h1 = hc(i); h2 = hc(i + 1);
            d1 = depth(i); d2 = depth(i + 1);
            if h1 == 0
                depths(end + 1, 1) = d1; %#ok<AGROW>
            end
            if (h1 < 0 && h2 > 0) || (h1 > 0 && h2 < 0)
                d0 = d1 + (-h1) / (h2 - h1) * (d2 - d1);
                depths(end + 1, 1) = d0; %#ok<AGROW>
            end
            saturatedThickness(c) = saturatedThickness(c) ...
                + saturatedSegmentLength(h1, h2, d1, d2);
        end
        if hc(end) == 0
            depths(end + 1, 1) = depth(end); %#ok<AGROW>
        end
        depths = sort(unique(depths), 'ascend');
        crossingCount(c) = numel(depths);
        if ~isempty(depths)
            crossingDepth(1:numel(depths), c) = depths;
            shallowest(c) = depths(1);
        end
    end

    diag = struct();
    diag.node_saturated = h >= 0;
    diag.crossing_depth_cm = crossingDepth;
    diag.crossing_count = crossingCount;
    diag.shallowest_crossing_depth_cm = shallowest;
    diag.saturated_thickness_cm = saturatedThickness;
    diag.saturated_fraction = saturatedThickness / profileThickness;
end

function lengthCm = saturatedSegmentLength(h1, h2, d1, d2)
% Integrate the h >= 0 portion of one linearly interpolated element.
    fullLength = abs(d2 - d1);
    if h1 >= 0 && h2 >= 0
        lengthCm = fullLength;
    elseif h1 < 0 && h2 < 0
        lengthCm = 0;
    else
        d0 = d1 + (-h1) / (h2 - h1) * (d2 - d1);
        if h1 >= 0
            lengthCm = abs(d0 - d1);
        else
            lengthCm = abs(d2 - d0);
        end
    end
end
