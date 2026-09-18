function v = getNodeColumn(A, cellIndex)
%GETNODECOLUMN  Extract one vertical node column from column or regional arrays.
%
%   v = soil.water.getNodeColumn(A, cellIndex) returns an [NN x 1] vector from
%   either an [NN x 1] single-column array or an [NN x nCell] regional array.
%
%   Inputs:
%     A          [in] numeric  node field, [NN x 1] or [NN x nCell]
%     cellIndex  [in] scalar   horizontal cell index
%
%   Output:
%     v          [out] numeric [NN x 1] selected node column
%
%   Editor: Jizixing Zhu; Date: 2026-06-28.

    if isvector(A)
        v = A(:);
    else
        v = A(:, cellIndex);
    end
end
