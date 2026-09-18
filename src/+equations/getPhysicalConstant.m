function value = getPhysicalConstant(name)
%GETPHYSICALCONSTANT  Return one named value from the immutable constants set.
%
%   This accessor exists because MATLAB does not support direct field access
%   on a function return value. Process kernels that need several constants
%   should call equations.physicalConstants once and retain the local struct.
%
%   Editor: Jizixing Zhu; Date: 2026-07-20.

    constants = equations.physicalConstants();
    assert((ischar(name) || (isstring(name) && isscalar(name))), ...
        'equations:getPhysicalConstant:invalidName', ...
        'Physical-constant name must be scalar text.');
    name = char(name);
    assert(isfield(constants, name), ...
        'equations:getPhysicalConstant:unknownName', ...
        'Unknown physical constant: %s.', name);
    value = constants.(name);
end
