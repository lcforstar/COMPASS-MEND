function tf = canopyLayerProfilesEnabled(options)
%CANOPYLAYERPROFILESENABLED  Return the unified canopy-layer profile switch.
%   New site configs use options.canopy_layer_profiles. Legacy SCOPE input
%   files use options.calc_vert_profiles. This adapter keeps old inputs
%   readable while giving COMPASS one semantic switch for canopy layer
%   structure, layer diagnostics, and canopy_profiles output.

    tf = false;
    if isstruct(options) && isfield(options, 'canopy_layer_profiles')
        tf = logical(options.canopy_layer_profiles);
    elseif isstruct(options) && isfield(options, 'calc_vert_profiles')
        tf = logical(options.calc_vert_profiles);
    end
end
