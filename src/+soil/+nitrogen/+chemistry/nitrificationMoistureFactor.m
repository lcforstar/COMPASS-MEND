function factor = nitrificationMoistureFactor(matric_head_cm)
%NITRIFICATIONMOISTUREFACTOR Dry-side water-potential limitation of nitrification.
%
%   factor = soil.nitrogen.chemistry.nitrificationMoistureFactor(h)
%   accepts a nonempty real finite numeric node vector h [cm of water],
%   bottom-to-surface, and returns a column vector factor [-] in [0,1].
%   The caller supplies the converged carrier-tile matric head; this pure
%   function never reads or changes accepted state or reaction ledgers.
%
%   For unsaturated h<0, pF=log10(-h [cm]). The retained DAISY dry branch is:
%     h >= -10^2.5 cm          : factor = 1;
%     -10^5 < h < -10^2.5 cm   : factor = (5-log10(-h))/2.5;
%     h <= -10^5 cm            : factor = 0.
%   Positive/saturated heads also return one. This intentionally excludes
%   DAISY's wet branch: COMPASS already resolves oxygen transport and its
%   f_aerobic response. It is a dry-side adaptation, not full DAISY parity.
%
%   Source checked 2026-10-28: DAISY 7.1.0
%     src/daisy/chemicals/nitrification.C, Nitrification::f_h, lines 44-60;
%     src/daisy/chemicals/nitrification_soil.C, NitrificationSoil::tick,
%       lines 65-69, multiplies the water factor into the total-NH4 MM rate;
%     include/util/mathlib.h, h2pF, defines log10(-h).
%   Source parameters pF=2.5 and pF=5 are not locally calibrated constants.
%   Invalid/nonfinite/complex/nonvector input raises
%   nitrificationMoistureFactor:head instead of silently selecting no stress.
%
%   Editor: Wenjie Liu; Date: 2026-10-28.

    assert(isnumeric(matric_head_cm)&&isreal(matric_head_cm) ...
        &&isvector(matric_head_cm)&&~isempty(matric_head_cm) ...
        &&all(isfinite(matric_head_cm(:))), ...
        'nitrificationMoistureFactor:head', ...
        'Matric head must be a nonempty real finite numeric vector in cm.');
    head = double(matric_head_cm(:));
    factor = ones(size(head));
    dryLimit = -1e5;
    transition = head < -10^2.5 & head > dryLimit;
    factor(transition) = (5-log10(-head(transition)))/2.5;
    factor(head <= dryLimit) = 0;
end
