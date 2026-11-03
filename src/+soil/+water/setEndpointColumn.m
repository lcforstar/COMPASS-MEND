function A = setEndpointColumn(A, cellIndex, E)
%SETENDPOINTCOLUMN  Insert one layer-endpoint matrix into column/regional arrays.
%
%   A = soil.water.setEndpointColumn(A, cellIndex, E) writes E into either an
%   [NL x nD] single-column endpoint field or an [NL x nD x nCell] regional
%   endpoint field.
%
%   Inputs:
%     A          [inout] numeric endpoint field
%     cellIndex  [in]    scalar  horizontal cell index
%     E          [in]    numeric [NL x nD] endpoint matrix
%
%   Output:
%     A          [out] numeric updated endpoint field
%
%   Editor: Jizixing Zhu; Date: 2026-06-28.

    if ismatrix(A)
        A = E;
    else
        A(:, :, cellIndex) = E;
    end
end
