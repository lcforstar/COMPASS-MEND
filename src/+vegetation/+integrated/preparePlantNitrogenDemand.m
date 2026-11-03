function context = preparePlantNitrogenDemand(model, kt, dtOuterActual)
%PREPAREPLANTNITROGENDEMAND  Build optional PlantCN demand for root N uptake.
%
%   context = vegetation.integrated.preparePlantNitrogenDemand(model,kt,dt)
%   computes plant nitrogen demand and the corresponding DAISY/Nacry demand
%   cap without calculating or committing root uptake. When PlantCN or soil N
%   is disabled, it returns an empty struct; the independent root-N module can
%   then run its kinetic, supply-limited formulation without plant demand.
%
%   Inputs:
%     model          [in] struct  initialized model; state is read-only
%     kt             [in] [-]     real outer-step index
%     dtOuterActual  [in] [s]     current outer-step duration
%
%   Returns:
%     context.demand_ctx  plant C/N demand inputs consumed by stepPlantCN
%     context.root_ctx    optional demand caps and root-uptake parameters
%
%   Physical closure:
%     PlantCN supplies an optional demand cap only. It neither calculates nor
%     owns the soil NH4/NO3 sink. Root uptake remains valid without context.
%
%   Reference: COMPASS.m:1200-1285 legacy demand assembly, reused through
%   vegetation.integrated.runPlantNDemand.
%
%   Editor: Jizixing Zhu; Date: 2026-07-11.
% Change: 2026-09-30 Wenjie Liu -- pass the immutable compiled PlantCN view.

    context = struct();
    built = model.built;
    sw = built.control.switches;
    if getFlag(sw, 'lPlantCN', 0) ~= 1 || getFlag(sw, 'lNefc', 0) ~= 1
        return;
    end

    sc = built.workspace.scope_step;
    fluxes = struct();
    if isfield(sc, 'fluxes') && isstruct(sc.fluxes), fluxes = sc.fluxes; end
    elapsed = 0;
    if isfield(model.runtime, 'step') && isfield(model.runtime.step, 'TIME')
        elapsed = model.runtime.step.TIME;
    end
    curDate = built.control.time.start_time + seconds(elapsed + dtOuterActual);
    NV = init.assembleNitrogenVariables(built, kt);
    cfg = assembleDemandConfig(built);
    meteo = vegetation.integrated.buildDemandMeteoView(built, kt);
    soilState = struct('TT', built.Y.T(:));
    acceptedPhenAlloc = buildPhenologyAllocation(built.Y.vegetation.phenology, ...
        built.parameters.veg.structure.PFT_type);
    currentPhenology=built.Y.vegetation.phenology;
    if isfield(sc,'phenology_candidate') && ~isempty(sc.phenology_candidate)
        currentPhenology=sc.phenology_candidate;
    end
    phenAlloc=buildPhenologyAllocation(currentPhenology, ...
        built.parameters.veg.structure.PFT_type);
    [demandCtx, rootCtx] = vegetation.integrated.runPlantNDemand( ...
        2, cfg, meteo, soilState, NV, built.p.plantHydro, fluxes, ...
        sc.canopy, curDate, built.control.time.Delt_t, phenAlloc);
    demandCtx.accepted_phen_alloc=acceptedPhenAlloc;
    demandCtx.previous_crop_phase=acceptedPhenAlloc.crop_phase;
    context = struct('demand_ctx', demandCtx, 'root_ctx', rootCtx, ...
        'phen_alloc', phenAlloc,'accepted_phen_alloc',acceptedPhenAlloc);
end

% -------------------------------------------------------------------------
function cfg = assembleDemandConfig(built)
% Build the legacy-shaped parameter view consumed by runPlantNDemand.
% Editor: Jizixing Zhu; Date: 2026-07-11.
    veg = built.parameters.veg;
    pcn = built.parameters.plant_cn;
    cfg.Veg = struct( ...
        'plant_cn', pcn, ...
        'PFT_type', veg.structure.PFT_type, ...
        'cn_flex_offset', pcn.coupling.cn_flex_offset, ...
        'cn_coupling_mode', pcn.coupling.cn_coupling_mode, ...
        'fcur', pcn.allocation.fcur, ...
        'storage_cushion', pcn.storage.storage_cushion, ...
        'root_solver_tol', veg.root.root_solver_tol, ...
        'root_solver_max_iter', veg.root.root_solver_max_iter, ...
        'laimx', veg.canopy.laimx);
end

% -------------------------------------------------------------------------
function pa = buildPhenologyAllocation(phenology, pftType)
% Convert the explicit accepted/current candidate phenology to allocation.
% Editor: Jizixing Zhu; Date: 2026-07-11.
    state = phenology.state;
    gddBudburst = NaN;
    fruitStage = 0;
    harvested = 0;
    if isstruct(state)
        if isfield(state, 'gdd5') && isfield(state, 'gdd5_at_budburst')
            gddBudburst = state.gdd5 - state.gdd5_at_budburst;
        end
        if isfield(state, 'DS_fruit') && ~isnan(state.DS_fruit)
            fruitStage = state.DS_fruit;
        end
        if isfield(state, 'harvested_flag')
            harvested = double(logical(state.harvested_flag));
        end
    end
    pa = struct('gdd_budburst', gddBudburst, 'DS_fruit', fruitStage, ...
        'phen', phenology.now, 'harvested', harvested, ...
        'is_crop', vegetation.getPFTValue(pftType, 'winter_wheat_cn', 0) == 1, ...
        'crop_live', fieldOr(state, 'crop_live', false), ...
        'crop_phase', fieldOr(state, 'crop_phase', 0), ...
        'crop_mature_flag',fieldOr(state,'crop_mature_flag',false), ...
        'crop_a_leaf', fieldOr(state, 'crop_a_leaf', 0), ...
        'crop_a_root', fieldOr(state, 'crop_a_root', 0), ...
        'crop_a_stem', fieldOr(state, 'crop_a_stem', 0), ...
        'crop_a_repr', fieldOr(state, 'crop_a_repr', 0));
end

function value = fieldOr(s, name, fallback)
value = fallback;
if isstruct(s) && isfield(s, name) && ~isempty(s.(name)), value = s.(name); end
end

% -------------------------------------------------------------------------
function value = getFlag(s, name, fallback)
% Read one normalized switch.
% Editor: Jizixing Zhu; Date: 2026-07-11.
    value = fallback;
    if isfield(s, name) && ~isempty(s.(name)), value = s.(name); end
end
