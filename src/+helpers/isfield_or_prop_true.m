function tf = isfield_or_prop_true(s, fieldName)
% isfield_or_prop_true True when a struct field or object property is truthy.

tf = false;
if isstruct(s) && isfield(s, fieldName) && ~isempty(s.(fieldName))
    tf = logical(s.(fieldName));
elseif isobject(s) && isprop(s, fieldName) && ~isempty(s.(fieldName))
    tf = logical(s.(fieldName));
end
end
