function [luna_state_next, buffer_struct_next] = runDailyLuna(model, kt, buffer)
%RUNDAILYLUNA  Advance the lagged daily LUNA state after an accepted step.
%
%   Inputs:
%     model  [in] Model state after accepted prognostic candidates commit.
%     kt     [-] One-based outer-step index.
%     buffer [in] Climate10dBuffer handle after the current entry is pushed.
%
%   Returns:
%     luna_state_next   [struct] Lagged LUNA state consumed on the next day.
%     buffer_struct_next [struct] Serialized ten-day climate buffer.
%
%   The daily optimization runs only at the configured day boundary after
%   the buffer is warm. This function never pushes the buffer and is called
%   only from the accepted outer-step commit path.
%   Editor: Jizixing Zhu; Date: 2026-07-16.
% Change: 2026-09-30 Wenjie Liu -- actual accepted green-leaf nitrogen,
% astronomical daylength, annual season reset and functional-N bounds.
% Incident RTMo PAR no longer divides by thermal-band emissivity in the
% physical-leaf-N mode; the explicit compatibility mode retains its adapter.
% Change: 2026-09-30 Wenjie Liu -- roll back this round's leaf-N LUNA
%   extension to the frozen COMPASS implementation; historical records above
%   are retained. Winter-wheat C/N allocation is handled separately.
    luna_state = model.built.Y.plant_acclimation.luna_state;
    luna_state_next = luna_state;
    buffer_struct_next = buffer.toStruct();
    acclimation = model.built.p.plant_acclimation;
    if ~acclimation.luna_enabled
        return;
    end
    if mod(kt, acclimation.steps_per_day) ~= 0 || ~buffer.isWarm()
        return;
    end

    built = model.built;
    pft_type = built.parameters.veg.structure.PFT_type;
    luna_parameters = vegetation.integrated.LUNA.lunaParams( ...
        vegetation.getPFTValue(pft_type, 'Jmaxb1'));

    leaf_bio = built.workspace.scope_step.leafbio;
    leaf_canopy = struct();
    leaf_canopy.LAI = built.workspace.scope_step.canopy.LAI;
    leaf_canopy.slatop = vegetation.getPFTValue(pft_type, 'SLA0');
    leaf_canopy.leafcn = vegetation.getPFTValue(pft_type, 'CN_L');
    leaf_canopy.rabsorb = ...
        1 - leaf_bio.rho_thermal - leaf_bio.tau_thermal;
    leaf_canopy.nrad = 1;

    c3_photosynthesis = 1;
    photo_type = readField( ...
        built.parameters.veg.photosynthesis, 'Photo_type', 'C3');
    if strcmpi(photo_type, 'C4')
        c3_photosynthesis = 0;
    end
    day_length_s = 12 * 3600;

    luna_input = buffer.assembleLunaInput( ...
        0, leaf_canopy, day_length_s, c3_photosynthesis);
    luna_output = vegetation.integrated.LUNA.computeLuna( ...
        luna_input, luna_state, luna_parameters);

    luna_state.Vcmax25_z = luna_output.Vcmax25_z;
    luna_state.Jmax25_z = luna_output.Jmax25_z;
    luna_state.PNlc_z = luna_output.PNlc_z;
    luna_state.enzs_z = luna_output.enzs_z;
    luna_state.first_day = false;
    luna_state.last_Vcmax25 = luna_output.Vcmax25_z(1);
    luna_state.last_Jmax25 = luna_output.Jmax25_z(1);
    luna_state.last_FNCa = luna_output.FNCa_z(1);
    luna_state.last_PNlc = luna_output.PNlc_z(1);
    luna_state.last_PNcb = luna_output.PNcb_z(1);
    light_capture_n = luna_output.PNlc_z(1) * luna_output.FNCa_z(1);
    luna_state.last_Cab_nlc = ...
        vegetation.integrated.LUNA.nLeafToCab(light_capture_n);

    luna_state_next = luna_state;
end

function value = readField(container, name, fallback)
% Read a nonempty struct field or return a local fallback.
    value = fallback;
    if isfield(container, name) && ~isempty(container.(name))
        value = container.(name);
    end
end

