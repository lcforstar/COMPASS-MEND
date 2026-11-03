function value = getfield_safe_numeric(s, fieldName, defaultValue)
% getfield_safe_numeric Read the first numeric field value or return a default.

value = defaultValue;
if isstruct(s) && isfield(s, fieldName) && isnumeric(s.(fieldName)) && ~isempty(s.(fieldName))
    value = double(s.(fieldName)(1));
end
end
