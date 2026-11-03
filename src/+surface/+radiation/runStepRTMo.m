function model = runStepRTMo(model)
%RUNSTEPRTMO  Run optical canopy transfer on finalized step canopy inputs.
%
%   Input:
%     model [inout] Model shell; reads scope_static and scope_step inputs and
%           replaces the current step rad, gap, and profiles auxiliaries.
%
%   Return:
%     model [out] Model with updated optical radiation auxiliaries.
%
%   This producer must run after the final LAI and canopy height are staged
%   and before the canopy energy balance. It updates auxiliary state only.
%   Editor: Jizixing Zhu; Date: 2026-07-16.
    scope_step = model.built.workspace.scope_step;
    scope_static = model.built.p.scope_static;
    [rad, gap, profiles] = surface.radiation.RTMo( ...
        scope_static.spectral, scope_step.atmo, scope_step.Soil, ...
        scope_step.leafopt, scope_step.canopy, scope_step.angles, ...
        scope_step.meteo, scope_step.rad, scope_static.options);
    scope_step.rad = rad;
    scope_step.gap = gap;
    scope_step.profiles = profiles;

    state = model.built.Y;
    scope_step.profiles.Ta_layerOLD = state.canopy_air.Ta_layer;
    scope_step.profiles.ea_layerOLD = state.canopy_air.ea_layer;
    scope_step.profiles.CO2_layerOLD = state.canopy_air.CO2_layer;
    scope_step.Soil.Tsold = state.scope_soil.Tsold;
    model.built.workspace.scope_step = scope_step;
end

