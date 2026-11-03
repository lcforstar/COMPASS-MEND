function pa = buildPlantAcclimationConfig(built)
%BUILDPLANTACCLIMATIONCONFIG  Static derived config for LUNA + leaf-respiration acclimation.
%
%   pa = init.buildPlantAcclimationConfig(built) returns the static, recomputable config
%   scalars that gate the LUNA daily optimization and size the leaf-respiration growth-
%   temperature window. These are derived from the vegetation parameters and the model
%   timestep, so they live in the auxiliary cache built.p (NOT prognostic, NOT checkpointed);
%   the prognostic acclimation STATE lives in Y.plant_acclimation (buildInitialState).
%
%   Inputs:
%     built  [in]  built model struct; reads built.parameters.veg.photosynthesis
%                  (lLUNA, Photo_LeafResp_Method) and built.control.time.Delt_t [s]
%
%   Returns:
%     pa     struct, fields:
%       .luna_enabled              [-]  logical, LUNA on/off (lLUNA)
%       .steps_per_day             [-]  round(86400/Delt_t) (shared LUNA + leaf_resp)
%       .leaf_resp_atkin_enabled   [-]  1 when Photo_LeafResp_Method == 2 (CLM5 Atkin2015), else 0
%       .leaf_resp_Tgrowth_window  [-]  10 * steps_per_day (10-day growth-temp window)
%
%   Method: mirrors the derived fields of the deleted prepareRuntime.prepareLuna /
%   prepareLeafResp (state-ownership migration increment 6a); the state those functions
%   also built now initializes in Y.plant_acclimation.
%
%   Reference: Ali et al. 2016 GMD 9:587 (LUNA 10-day mean climate driver);
%   Atkin et al. 2015 New Phytol 206:614 (leaf Rd25 growth-T acclimation).
%   Editor: Jizixing Zhu; Date: 2026-07-05.
% Change: 2026-09-30 Wenjie Liu -- make leaf-N source and coupled annual
% crop applicability explicit; fixed PFT nitrogen is a compatibility mode.
% Change: 2026-09-30 Wenjie Liu -- roll back this round's leaf-N LUNA
%   extension to the frozen COMPASS implementation; historical records above
%   are retained. Winter-wheat C/N allocation is handled separately.
    dt = built.control.time.Delt_t;
    stepsPerDay = max(1, round(86400 / dt));

    ps = struct();
    if isfield(built.parameters.veg, 'photosynthesis')
        ps = built.parameters.veg.photosynthesis;
    end
    lunaEnabled = isfield(ps, 'lLUNA') && logical(ps.lLUNA);
    method = 0;
    if isfield(ps, 'Photo_LeafResp_Method')
        method = ps.Photo_LeafResp_Method;
    end

    pa = struct( ...
        'luna_enabled',             lunaEnabled, ...
        'steps_per_day',            stepsPerDay, ...
        'leaf_resp_atkin_enabled',  double(method == 2), ...
        'leaf_resp_Tgrowth_window', 10 * stepsPerDay);
end
