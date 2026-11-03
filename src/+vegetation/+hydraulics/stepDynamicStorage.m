function [water_storage_next, diagnostics_out] = stepDynamicStorage(model, dt_outer_s)
%STEPDYNAMICSTORAGE  Stage dynamic organ-water storage from ebal potentials.
%
%   Inputs:
%     model      [in] Model shell after the canopy energy balance.
%     dt_outer_s [s] Actual outer-step duration.
%
%   Returns:
%     water_storage_next [m, 3x1] Root, stem, and leaf water storage; empty
%                        when the dynamic hydraulic scheme is disabled.
%     diagnostics_out    [struct] Dynamic hydraulic diagnostics.
%
%   Maps accepted storage and current algebraic water potentials through
%   updateOrganStorageFromPsi. It commits no state.
%   Active when hydraulic_scheme is storage_from_psi or storage_dynamic
%   (legacy alias: dynamic).
%   Editor: Jizixing Zhu; Date: 2026-07-16.
%   Change: 2026-07-17 Jizixing Zhu — rename schemes to steady_spac /
%     storage_from_psi / storage_dynamic with legacy aliases.
%   Change: 2026-07-17 Jizixing Zhu — call updateOrganStorageFromPsi
%     (was Plant_hydraulics_v2_wrapper).

    water_storage_next = [];
    diagnostics_out = struct();
    built = model.built;
    scheme = readField( ...
        built.parameters.veg.hydraulic, ...
        'hydraulic_scheme', 'steady_spac');
    if ~vegetation.hydraulics.usesOrganWaterStorage(scheme)
        return;
    end

    scope_step = built.workspace.scope_step;
    if ~isfield(scope_step, 'hydraulic') ...
            || ~isfield(scope_step.hydraulic, 'psil') ...
            || ~isfield(scope_step.hydraulic, 'psir')
        return;
    end
    plant_state = built.Y.plant_hydraulic;
    water_storage_old = [ ...
        plant_state.W_root; plant_state.W_stem; plant_state.W_leaf];
    lai = readField(scope_step.canopy, 'LAI', 0);
    [water_storage_next, ~, diagnostics_out] = ...
        vegetation.hydraulics.updateOrganStorageFromPsi( ...
        water_storage_old, scope_step.hydraulic, struct(), ...
        built.p.plantHydro, lai, dt_outer_s);

end

function value = readField(container, name, fallback)
% Read a nonempty struct field or return a local fallback.
    value = fallback;
    if isfield(container, name) && ~isempty(container.(name))
        value = container.(name);
    end
end

