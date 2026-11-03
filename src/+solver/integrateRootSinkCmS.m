function qSink = integrateRootSinkCmS(Srt, DeltZ)
%INTEGRATEROOTSINKCMS  Integrate the Richards root sink over soil depth.
%
%   qSink = solver.integrateRootSinkCmS(Srt, DeltZ) applies the same linear
%   finite-element endpoint quadrature as the Richards sink assembly.
%
%   Inputs:
%     Srt   [in] [1/s, NL x 2] root-water sink at element endpoints.
%     DeltZ [in] [cm, NL x 1] element thickness, bottom-up.
%
%   Returns:
%     qSink [out] [cm/s] depth-integrated water removal. Positive removes
%           water; negative values retain hydraulic redistribution.
%
%   Method:
%     q_sink = sum_e dz_e*(Srt_left + Srt_right)/2.
%
%   Editor: Jizixing Zhu; Date: 2026-07-12.

    NL = numel(DeltZ);
    assert(isequal(size(Srt), [NL, 2]), ...
        'integrateRootSinkCmS:shape', 'Srt must be NL x 2.');
    assert(all(isfinite(Srt(:))) && all(isfinite(DeltZ(:))) ...
        && all(DeltZ(:) > 0), ...
        'integrateRootSinkCmS:finite', ...
        'Srt and positive DeltZ must be finite.');
    qSink = sum(0.5 * (Srt(:, 1) + Srt(:, 2)) .* DeltZ(:));
end
