function model = applyEbalOutput(model, out)
%APPLYEBALOUTPUT  Stage energy-balance results into workspace.scope_step.
%
%   model = vegetation.energy.applyEbalOutput(model, out) writes target-name output
%   fields onto the legacy scope_step keys expected by post-ebal consumers.
%
%   Returns:
%     model [out]  Model with scope_step.iter/fluxes/rad/gap/thermal/... staged
%
%   Inputs:
%     model [inout] struct  COMPASS model shell
%     out   [in]    struct  From vegetation.energy.solveCanopySoilEnergy
%
%   Method:
%     Maps out.* target names to existing scope_step consumer keys (contract
%     §5.3). Does not commit Y; accept path remains step.commitPending.
%
%   Editor: Jizixing Zhu; Date: 2026-07-16.

    assert(isstruct(out) && isfield(out, 'fluxes') && isfield(out, 'thermal'), ...
        'applyEbalOutput:badOutput', ...
        'out must be a target-name struct from solveCanopySoilEnergy.');

    scope_step = model.built.workspace.scope_step;
    scope_step.iter = out.iteration;
    scope_step.fluxes = out.fluxes;
    scope_step.rad = out.radiation;
    scope_step.gap = out.gap;
    scope_step.thermal = out.thermal;
    scope_step.profiles = out.profiles;
    scope_step.Soil = out.soil_surface;
    scope_step.RWU = out.root_water_uptake;
    scope_step.frac = out.frac;
    scope_step.hydraulic = out.hydraulic;
    model.built.workspace.scope_step = scope_step;
end
