function k = findNearestNode(depth, depths_target)
% Index of the depths_target node nearest `depth` (exact match preferred),
% mapped to the VG-parameter index (params hold NL = numel(depths_target)-1).
% Editor: Jizixing Zhu; Date: 2026-06-27.
    tol = 1e-10;
    hit = abs(depths_target - depth) < tol;
    if any(hit)
        k = find(hit, 1);
    else
        [~, k] = min(abs(depths_target - depth));
    end
    k = min(max(k - 1, 1), numel(depths_target) - 2);
end
