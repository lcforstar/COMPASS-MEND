function model = applyCanopyInterception(model, kt, dt_outer_s)
%APPLYCANOPYINTERCEPTION  Stage canopy rainfall storage and wet fraction.
%
%   Inputs:
%     model      [inout] Model shell with finalized current canopy geometry.
%     kt         [-] One-based outer-step index.
%     dt_outer_s [s] Actual outer-step duration.
%
%   Return:
%     model      [out] Model with interception auxiliary state staged for
%                canopy energy and later water-boundary producers.
%
%   The previous accepted canopy-water reservoir is read from prognostic
%   state. Reservoir evaporation and accept-only state advancement occur in
%   the post-energy stage.
%   Editor: Jizixing Zhu; Date: 2026-07-16.
    built = model.built;
    canopy_config = built.parameters.veg.canopy;
    canopy = built.workspace.scope_step.canopy;
    canopy_water_view = struct( ...
        'state_h2ocan', built.Y.canopy_water.state_h2ocan, ...
        'state_h2ocan_sun', built.Y.canopy_water.state_h2ocan_sun, ...
        'state_h2ocan_sha', built.Y.canopy_water.state_h2ocan_sha, ...
        'Ps_max_diurnal', readField( ...
            canopy_config, 'Ps_max_diurnal', 0.7), ...
        'Ps_min_diurnal', readField( ...
            canopy_config, 'Ps_min_diurnal', 0.0));
    interception = vegetation.computeCanopyInterception( ...
        readField(canopy_config, 'canopy_water_scheme', 'single'), ...
        surfacePrecipitation(built, kt), dt_outer_s, canopy.LAI, ...
        readField(canopy, 'kV', 0.6), ...
        readField(canopy_config, 'canopy_storage_per_LAI', 0.1), ...
        readField(canopy_config, 'interc_F0', 0.15), ...
        canopy_water_view);
    model.built.workspace.scope_step.interception = interception;
end

function precipitation_cm_s = surfacePrecipitation(built, kt)
% Read the current liquid precipitation forcing in the surface-water units.
    precipitation_cm_s = 0;
    if isfield(built, 'forcing') && isfield(built.forcing, 'hourly') ...
            && isfield(built.forcing.hourly, 'Precip') ...
            && numel(built.forcing.hourly.Precip) >= kt
        precipitation_cm_s = built.forcing.hourly.Precip(kt);
    end
end

function value = readField(container, name, fallback)
% Read a nonempty struct field or return a local fallback.
    value = fallback;
    if isfield(container, name) && ~isempty(container.(name))
        value = container.(name);
    end
end
