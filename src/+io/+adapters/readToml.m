function cfg = readToml(filepath)
%READTOML  Minimal TOML reader for the COMPASS site configuration (subset).
%
%   cfg = io.adapters.readToml(filepath) parses the constrained TOML subset
%   used by site.toml and returns a nested struct: cfg.<section>.<key> = value.
%   This is the format adapter for stage 1 (pure read); it performs NO grid
%   construction, interpolation, or unit conversion (those live in buildState).
%
%   Inputs:
%     filepath  [char]  path to a .toml file
%
%   Returns:
%     cfg       struct  sections become sub-structs; values typed as below.
%
%   Supported subset (sufficient for site.toml):
%     - [section] and [section.sub] headers (dotted name -> nested struct)
%     - key = value, value in { "string", true|false, number, [array] }
%     - array of numbers -> [1xN double]; array of "strings" -> cellstr
%     - full-line '#' comments (incl. the '#:schema ...' directive), blanks,
%       and inline '# ...' comments (not stripped inside quoted strings)
%   Datetimes are kept as quoted strings (site.toml quotes them); the caller
%   converts. NOT a general TOML parser (no multiline strings / inline tables /
%   nested arrays / dotted keys outside section headers).

    fid = fopen(filepath, 'r');
    if fid < 0
        error('readToml:cannotOpen', 'Cannot open TOML file: %s', filepath);
    end
    closer = onCleanup(@() fclose(fid));   % guarantee close on any exit

    cfg = struct();
    sectionPath = {};                       % current [section] split on '.'

    while true
        raw = fgetl(fid);
        if ~ischar(raw)                     % EOF
            break;
        end
        line = strtrim(raw);
        if isempty(line) || line(1) == '#'  % blank / comment / '#:schema'
            continue;
        end
        if line(1) == '['                   % section header
            closeIdx = find(line == ']', 1);
            if isempty(closeIdx)
                error('readToml:badSection', 'Malformed section header: %s', raw);
            end
            sectionPath = strsplit(strtrim(line(2:closeIdx - 1)), '.');
            continue;
        end
        eqIdx = find(line == '=', 1);       % key = value
        if isempty(eqIdx)
            error('readToml:badLine', 'Line is not [section] or key=value: %s', raw);
        end
        key    = strtrim(line(1:eqIdx - 1));
        valStr = strtrim(stripInlineComment(line(eqIdx + 1:end)));
        cfg    = assignNested(cfg, [sectionPath, {key}], parseTomlValue(valStr));
    end
end

% ---------------------------------------------------------------------------
function s = stripInlineComment(s)
% Remove a trailing '# ...' inline comment, ignoring '#' inside quoted strings.
    inStr = false;
    for i = 1:numel(s)
        if s(i) == '"'
            inStr = ~inStr;
        elseif s(i) == '#' && ~inStr
            s = s(1:i - 1);
            return;
        end
    end
end

% ---------------------------------------------------------------------------
function v = parseTomlValue(s)
% Type a trimmed TOML scalar/array value string.
    s = strtrim(s);
    if isempty(s)
        v = '';
        return;
    end
    if s(1) == '"'                          % quoted string
        v = stripQuotes(s);
        return;
    end
    if s(1) == '['                          % array
        inner = strtrim(s(2:end));
        if ~isempty(inner) && inner(end) == ']'
            inner = strtrim(inner(1:end - 1));
        end
        parts = splitTopLevel(inner);
        if isempty(parts)
            v = [];
            return;
        end
        first = strtrim(parts{1});
        if ~isempty(first) && first(1) == '"'   % cellstr array
            v = cellfun(@(p) stripQuotes(strtrim(p)), parts, 'UniformOutput', false);
        else                                    % numeric array
            v = cellfun(@(p) str2double(strtrim(p)), parts);
        end
        return;
    end
    if strcmp(s, 'true')                    % boolean
        v = true;
        return;
    end
    if strcmp(s, 'false')
        v = false;
        return;
    end
    num = str2double(s);                     % number
    if ~isnan(num)
        v = num;
        return;
    end
    v = s;                                   % fallback: bare (unquoted) string
end

% ---------------------------------------------------------------------------
function s = stripQuotes(s)
    s = strtrim(s);
    if numel(s) >= 2 && s(1) == '"' && s(end) == '"'
        s = s(2:end - 1);
    end
end

% ---------------------------------------------------------------------------
function parts = splitTopLevel(s)
% Split on commas that are not inside quoted strings.
    parts = {};
    inStr = false;
    startIdx = 1;
    for i = 1:numel(s)
        if s(i) == '"'
            inStr = ~inStr;
        elseif s(i) == ',' && ~inStr
            parts{end + 1} = s(startIdx:i - 1); %#ok<AGROW>
            startIdx = i + 1;
        end
    end
    tailStr = strtrim(s(startIdx:end));
    if ~isempty(tailStr)
        parts{end + 1} = tailStr; 
    end
end

% ---------------------------------------------------------------------------
function root = assignNested(root, path, val)
% Assign val into root.(path{1}).(path{2})... creating sub-structs as needed.
    if isscalar(path)
        root.(path{1}) = val;
    else
        head = path{1};
        if ~isfield(root, head) || ~isstruct(root.(head))
            root.(head) = struct();
        end
        root.(head) = assignNested(root.(head), path(2:end), val);
    end
end
