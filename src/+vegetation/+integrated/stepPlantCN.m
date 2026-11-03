function [plantCnNext, neeOut, dbg] = stepPlantCN( ...
        model, kt, dtOuterActual, phenNext, demandContext, rootNitrogenResult)
%STEPPLANTCN  Advance PlantCN using independently calculated root-N uptake.
%
%   [plantCnNext,neeOut,dbg] = vegetation.integrated.stepPlantCN(...)
%   advances the plant carbon/nitrogen pools and NEE diagnostics. Root water
%   and nitrogen uptake are upstream inputs: this function neither calculates
%   nor commits Srt or SinkS. The caller commits plantCnNext only after the
%   outer step is accepted.
%
%   Inputs:
%     model               [in] struct  current staged model; read-only
%     kt                  [in] [-]     real outer-step index
%     dtOuterActual       [in] [s]     current outer-step duration
%     phenNext            [in] struct  pending phenology state
%     demandContext       [in] struct  output of preparePlantNitrogenDemand
%     rootNitrogenResult  [in] struct  independent root-N flux result
%
%   Returns:
%     plantCnNext  [out] struct  pending two-slot PlantCN working frame
%     neeOut       [out] struct  NEE and respiration diagnostics
%     dbg          [out] struct  golden-test intermediates
%
%   Physical closure:
%     PlantCN may constrain uptake demand upstream and consumes accepted-step
%     U_N_total here. It cannot alter the independently computed soil SinkS.
%
%   Reference: COMPASS.m:1200-1421 plant C/N and NEE sequence.
% Editor: Jizixing Zhu; Date: 2026-07-04.
% Change: 2026-07-05 Jizixing Zhu -- read lagged committed phenology from Y.vegetation.
% Change: 2026-07-11 Jizixing Zhu -- remove root water/N uptake ownership; consume independent U_N_total only.
% Change: 2026-08-21 Jizixing Zhu -- tolerate sites whose non-crop phenology
%   driver has no phen_HUI field when live litter coupling enables PlantCN.
% Change: 2026-09-30 Wenjie Liu -- retain the immutable accepted PlantCN
% input and frozen forcing for winter-wheat replay with actual soil N uptake.

    if nargin < 5 || isempty(demandContext), demandContext = struct(); end
    if nargin < 6 || isempty(rootNitrogenResult)
        rootNitrogenResult = struct('U_N_total', 0);
    end

    built = model.built;
    sw = built.control.switches;
    Delt_t = built.control.time.Delt_t;
    sc = built.workspace.scope_step;
    canopy = sc.canopy;
    fluxes = struct();
    if isfield(sc, 'fluxes') && isstruct(sc.fluxes), fluxes = sc.fluxes; end

    elapsed = 0;
    if isfield(model.runtime, 'step') && isfield(model.runtime.step, 'TIME')
        elapsed = model.runtime.step.TIME;
    end
    curDate = built.control.time.start_time + seconds(elapsed + dtOuterActual);
    NV = init.assembleNitrogenVariables(built, kt);
    neeOut = [];
    dbg = struct();

    if sw.lNefc == 1
        assert(isfield(demandContext, 'demand_ctx') && ~isempty(demandContext.demand_ctx), ...
            'stepPlantCN:missingDemandContext', ...
            'PlantCN with soil nitrogen enabled requires prepared plant N demand.');
        demandCtx = demandContext.demand_ctx;
        acceptedNV = NV;
        U_N_total = getField(rootNitrogenResult, 'U_N_total', 0);
        NV.cum_U_N = NV.cum_U_N + U_N_total * Delt_t;

        envTurn = buildEnvTurn(built, phenNext, canopy, curDate, kt, Delt_t);
        [NV, R_a_plant_kg] = vegetation.integrated.runPlantCNStep( ...
            NV, demandCtx, U_N_total, envTurn, Delt_t, 2);
        % The final AOM/MEND CO2 and root-N reconciliation are downstream.
        % Exchange is materialized from final owner increments in runOneStep.
        % Change: 2026-09-30 Wenjie Liu -- retire empirical soil Rh from PlantCN.
        neeOut = [];
        dbg = struct('demand_ctx', demandCtx, 'U_N_total', U_N_total, ...
            'R_a_plant_kg', R_a_plant_kg, 'env_turn', envTurn);
        if demandCtx.plant_cn_params.resolved.winter_wheat_cn == 1
            % Ephemeral request material: this is not a second accepted
            % physical state and is never installed in Y or checkpoints.
            dbg.plant_cn_replay = struct('schema_version',1, ...
                'accepted_nv',acceptedNV,'demand_ctx',demandCtx, ...
                'env_turn',envTurn,'dt_s',Delt_t);
        end
    elseif isfield(fluxes, 'GPP')
        R_h = 0;
        if isfield(fluxes, 'Resp'), R_h = fluxes.Resp * 12 / 1e9; end
        [nee, ~] = vegetation.integrated.computeNEE(fluxes.GPP, R_h, 0);
        neeOut = struct('NEE', nee, 'Resp', getField(fluxes, 'Resp', NaN));
    end

    plantCnNext = NV;
end

% -------------------------------------------------------------------------
function envTurn = buildEnvTurn(built, phenNext, canopy, curDate, kt, Delt_t)
% Build turnover forcing from pending phenology and current forcing.
% Editor: Jizixing Zhu; Date: 2026-07-04.
% Change: 2026-07-11 Jizixing Zhu -- retain PlantCN-only turnover assembly after uptake separation.
    stepsPerDay = round(86400 / Delt_t);
    dayIndex = ceil(kt / stepsPerDay);
    frPHU = 0;
    phenHui = [];
    if isfield(built.vegetation_state.phen_driver,'phen_HUI')
        phenHui = built.vegetation_state.phen_driver.phen_HUI;
    end
    if ~isempty(phenHui) && dayIndex <= numel(phenHui), frPHU = phenHui(dayIndex); end
    veg = built.parameters.veg;
    envTurn = struct( ...
        'frPHU', frPHU, ...
        'phen', phenNext.now, ...
        'phen_state', phenNext.state, ...
        'doy_now', day(curDate, 'dayofyear'), ...
        'LAI', canopy.LAI, ...
        'DLAI', veg.phenology.DLAI, ...
        'use_leafc_sla', strcmpi(veg.canopy.LAI_engine, 'leafc_sla_phen'), ...
        'year', year(curDate));
end

% -------------------------------------------------------------------------
function value = getField(s, name, fallback)
% Read one optional field.
% Editor: Jizixing Zhu; Date: 2026-07-11.
    value = fallback;
    if isstruct(s) && isfield(s, name), value = s.(name); end
end
