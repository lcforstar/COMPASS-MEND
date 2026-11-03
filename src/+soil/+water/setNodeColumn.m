function A = setNodeColumn(A, cellIndex, v)
%SETNODECOLUMN  Insert one vertical node column into column or regional arrays.
%
%   A = soil.water.setNodeColumn(A, cellIndex, v) writes v into either an
%   [NN x 1] single-column field or an [NN x nCell] regional field.
%
%   Inputs:
%     A          [inout] numeric node field, [NN x 1] or [NN x nCell]
%     cellIndex  [in]    scalar  horizontal cell index
%     v          [in]    numeric [NN x 1] node column
%
%   Output:
%     A          [out] numeric updated node field
%
%   Editor: Jizixing Zhu; Date: 2026-06-28.

    if isvector(A)
        A = v(:);
    else
        A(:, cellIndex) = v(:);
    end
end
