function SoilVariables = unpackSoilPicardState(SoilVariables, x, meta)
%UNPACKSOILPICARDSTATE Write primary Picard unknowns back to SoilVariables.
%
%   SoilVariables = solver.unpackSoilPicardState(SoilVariables, x, meta)
%
%   The meta struct should come from solver.packSoilPicardState. This
%   function updates only hh, TT or UU, and optionally P_gg. All derived physical
%   fields must be refreshed by the normal model update sequence.
%
%   Editor: Jizixing Zhu; Date: 2026-07-14.
%   Change: 2026-07-14 Jizixing Zhu -- unpack U-primary Anderson states and
%           support the packer's frozen nodal energy scale.

x = x(:);
expectedLength = sum(meta.lengths);
if numel(x) ~= expectedLength
    error('solver:unpackSoilPicardState:SizeMismatch', ...
        'State vector length does not match meta information.');
end

idx0 = 0;
for k = 1:numel(meta.names)
    len = meta.lengths(k);
    idx = idx0 + (1:len);
    if isfield(meta, 'scaleVectors') && numel(meta.scaleVectors) >= k
        scale = meta.scaleVectors{k};
    else
        scale = meta.scales(k);
    end
    values = x(idx) .* scale;

    switch meta.names{k}
        case 'hh'
            SoilVariables.hh = values(:);
        case 'TT'
            SoilVariables.TT = values(:);
        case 'UU'
            SoilVariables.UU = values(:);
        case 'P_gg'
            SoilVariables.P_gg = values(:);
        otherwise
            error('solver:unpackSoilPicardState:UnknownField', ...
                'Unknown packed state field: %s', meta.names{k});
    end

    idx0 = idx0 + len;
end
end
