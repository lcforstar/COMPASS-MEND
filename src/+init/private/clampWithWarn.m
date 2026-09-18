function v = clampWithWarn(v, lo, hi, name)
% Clamp v to [lo, hi], warning once if any element was out of range.
% Editor: Jizixing Zhu; Date: 2026-06-27.
    if any(v < lo) || any(v > hi)
        warning('init:clamped', '%s out of [%.2f, %.2f] (range [%.3f, %.3f]); clamping.', ...
            name, lo, hi, min(v), max(v));
    end
    v = min(max(v, lo), hi);
end
