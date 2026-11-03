function value = getVegetationProperty(veg, primary_name, aliases, default_value)
%GETVEGETATIONPROPERTY  Read a canonical vegetation field or legacy alias.
% Editor: Jizixing Zhu; Date: 2026-07-20.
if nargin < 3 || isempty(aliases)
    aliases = {};
end
if nargin < 4
    default_value = [];
end

names = [{primary_name}, aliases(:)'];
for i = 1:numel(names)
    name = names{i};
    if isstruct(veg)
        if isfield(veg, name)
            candidate = veg.(name);
        else
            continue;
        end
    else
        if isprop(veg, name)
            candidate = veg.(name);
        else
            continue;
        end
    end
    if ~isempty(candidate)
        value = candidate;
        return;
    end
end

value = default_value;
end
