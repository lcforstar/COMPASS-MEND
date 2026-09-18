function [value, name] = getEnv(suffix)
%GETENV  Read a COMPASS_<suffix> env var, with legacy STEMMUS_<suffix> fallback.
%
%   [value, name] = io.getEnv(suffix) returns getenv('COMPASS_<suffix>'); when that
%   is unset it falls back to getenv('STEMMUS_<suffix>'), so run scripts written for
%   the legacy STEMMUS_* prefix keep working during the prefix migration.
%
%   Inputs:
%     suffix  [char]  variable name without prefix, e.g. 'OUTPUT' for COMPASS_OUTPUT.
%
%   Outputs:
%     value  [char]  the env value; '' when neither prefix is set (matches getenv).
%     name   [char]  the variable actually read (COMPASS_* when set, else the
%                    STEMMUS_* legacy name), useful for provenance/logging.
%
%   Editor: Jizixing Zhu; Date: 2026-06-28.

    name  = ['COMPASS_' suffix];
    value = getenv(name);
    if isempty(value)
        name  = ['STEMMUS_' suffix];
        value = getenv(name);
    end
end
