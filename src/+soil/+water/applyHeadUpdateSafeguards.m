function [hUpdated, rawIncrement, dryHeadProjected] = applyHeadUpdateSafeguards( ...
        hCurrent, hRaw, hcrit, isRaining, allowPositiveHead, firstActiveNode, ...
        enforceDryHead)
%APPLYHEADUPDATESAFEGUARDS  Apply the shared Richards head-update policy.
%
%   [hUpdated, rawIncrement] = soil.water.applyHeadUpdateSafeguards(...)
%   applies the existing COMPASS large-jump damping and physical head guards
%   to either an absolute Thomas map or a correction-coordinate map.
%
%   Inputs:
%     hCurrent/hRaw [cm, NN x 1] current and raw proposed matric head
%     hcrit         [cm] dry surface lower bound
%     isRaining     [-] true releases the dry-surface hcrit guard
%     allowPositiveHead [-] true permits variably saturated positive heads
%     firstActiveNode [-] first node solved by the water equation
%     enforceDryHead  [-] optional legacy post-solve hcrit projection;
%                         default true for compatibility. Typed boundaries
%                         set false because their active set owns this choice.
%
%   rawIncrement is abs(hRaw-hCurrent), before damping.  When its maximum is
%   above 200 cm, negative heads are blended in log(abs(h)) space with
%   omega=0.5, exactly as in the former local implementation.
%
%   Editor: Jizixing Zhu; Date: 2026-07-14.
%   Change: 2026-07-17 Jizixing Zhu -- report whether the post-damping top
%           candidate was actually projected to hcrit for active-set logic.
%   Change: 2026-07-17 Jizixing Zhu -- make the legacy dry projection optional so it
%           cannot overwrite a typed Neumann/Dirichlet active-set row.

    if nargin < 7 || isempty(enforceDryHead)
        enforceDryHead = true;
    end

    hCurrent = hCurrent(:);
    hRaw = hRaw(:);
    NN = numel(hCurrent);
    assert(numel(hRaw) == NN && all(isfinite(hCurrent)) && all(isfinite(hRaw)) ...
        && isscalar(hcrit) && isfinite(hcrit) ...
        && isscalar(firstActiveNode) && firstActiveNode >= 1 ...
        && firstActiveNode <= NN && firstActiveNode == floor(firstActiveNode) ...
        && isscalar(enforceDryHead), ...
        'applyHeadUpdateSafeguards:badInput', ...
        'Head vectors and safeguard controls must be finite and consistent.');

    rawIncrement = abs(hRaw - hCurrent);
    if max(rawIncrement) > 200
        omega = 0.5;
        absRaw = max(abs(hRaw), 1e-6);
        absCurrent = max(abs(hCurrent), 1e-6);
        hUpdated = -exp(omega * log(absRaw) ...
            + (1 - omega) * log(absCurrent));
        positive = hRaw > 0;
        hUpdated(positive) = hRaw(positive);
    else
        hUpdated = hRaw;
    end

    dryHeadProjected = false;
    for i = firstActiveNode:NN
        if ~allowPositiveHead && hUpdated(i) > -1e-6
            hUpdated(i) = -1e-6;
        elseif i == NN && logical(enforceDryHead) ...
                && hUpdated(i) < hcrit && ~isRaining
            hUpdated(i) = hcrit;
            dryHeadProjected = true;
        elseif hUpdated(i) < -1e7
            hUpdated(i) = -1e7;
        end
    end
end
