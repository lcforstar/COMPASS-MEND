function x = getEnvDouble(suffix, defaultValue)
%GETENVDOUBLE  Read one numeric COMPASS_<suffix> env override, with STEMMUS_<suffix> fallback.
%
%   x = io.getEnvDouble(suffix, defaultValue) returns str2double(io.getEnv(suffix)); when the env is
%   unset OR does not parse to a finite number it returns defaultValue. Shared by prepareRuntime and
%   buildVeg so numeric env overrides keep ONE parse (COMPASS_ priority, STEMMUS_ fallback, invalid
%   -> default). Centralizes what was a per-caller readEnvDouble local.
%
%   Inputs:
%     suffix        [char]   variable name without prefix (e.g. 'PS_MAX' for COMPASS_PS_MAX)
%     defaultValue  [double] value used when unset or non-finite
%
%   Returns:
%     x             [double] the parsed override, or defaultValue
%
%   Reference: wraps io.getEnv (COMPASS_/STEMMUS_ dual-prefix lookup).
%   Editor: Jizixing Zhu; Date: 2026-07-05.
    raw = io.getEnv(suffix);
    if isempty(raw)
        x = defaultValue;
        return;
    end
    x = str2double(raw);
    if ~isfinite(x)
        x = defaultValue;
    end
end
