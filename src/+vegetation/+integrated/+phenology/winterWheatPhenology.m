function [state, diag] = winterWheatPhenology(state, p, forcing, doy)
%WINTERWHEATPHENOLOGY  CoLM/CLM5 winter-wheat crop state machine.
%
% Inputs are one accepted model-day forcing slice: T_avg [degC], configured
% sowing_day [DOY], and immutable crop parameters. State is a pending
% PhenologyState candidate and is committed only with the outer transaction.
%
% Equations and branch order are retained from CoLM
% main/BGC/MOD_BGC_Veg_CNPhenology.F90:293-313, 1044-1189, 1591-1610 and
% MOD_BGC_Veg_NutrientCompetition.F90:389-454. COMPASS only adapts the
% time representation from CoLM subdaily rates to one daily update.
%
% Editor: Jizixing Zhu; Date: 2026-08-24.
% Change: 2026-09-30 Wenjie Liu -- observed_lai_zero management separates
%   the unchanged CoLM thermal maturity criterion from physical harvest.
%   The awaiting-harvest crop remains phase 3/live, stops organ allocation,
%   and retains its physical pools. thermal_maturity preserves the reference
%   branch. Observed harvest is an explicit COMPASS management adaptation.

state.crop_planting_event = false;
state.crop_emergence_event = false;
state.crop_harvest_event = false;
state.crop_grainfill_event = false;
state.crop_maturity_event = false;
if ~isfield(state,'crop_mature_flag'), state.crop_mature_flag = false; end
if ~isfield(state,'crop_season_id'), state.crop_season_id = 0; end
observedHarvest = isfield(p,'crop_harvest_rule') ...
    && strcmpi(p.crop_harvest_rule,'observed_lai_zero');
state.onset_flag = false;
state.offset_flag = false;

if ~state.crop_live && doy == p.sowing_day
    state.crop_live = true;
    state.crop_phase = 1;
    state.crop_planting_event = true;
    state.crop_days_after_planting = 0;
    state.crop_gddplant = 0;
    state.crop_cumvd = 0;
    state.crop_vf = 0;
    state.crop_a_leaf = 1e-5;
    state.crop_a_root = 0;
    state.crop_a_stem = 0;
    state.crop_a_repr = 0;
    state.budburst_doy = doy;
    state.harvested_flag = false;
    state.crop_mature_flag = false;
    state.crop_season_id = state.crop_season_id + 1;
end

if state.crop_live
    phaseBefore = state.crop_phase;
    if state.crop_phase == 2
        gddFactor = state.crop_vf;
    else
        gddFactor = 1;
    end
    state.crop_gddplant = state.crop_gddplant + gddFactor * ...
        max(0, forcing.T_avg - p.baset);
    state.crop_gddmaturity = 0.42 * p.gdd1020 + 440;
    state.crop_hui = state.crop_gddplant / state.crop_gddmaturity;
    state.crop_days_after_planting = state.crop_days_after_planting + 1;

    hui = state.crop_hui;
    idpp = state.crop_days_after_planting;
    reachesMaturity = hui >= 1 || idpp >= p.mxmat;
    if reachesMaturity && ~state.crop_mature_flag
        state.crop_mature_flag = true;
        state.crop_maturity_event = true;
    end
    if reachesMaturity && ~observedHarvest
        state.crop_live = false;
        state.crop_phase = 4;
        state.crop_harvest_event = true;
        state.offset_flag = true;
        state.harvested_flag = true;
        state.harvest_doy = doy;
        state.crop_hui = 0;
        state.crop_a_leaf = 0;
        state.crop_a_root = 0;
        state.crop_a_stem = 0;
        state.crop_a_repr = 0;
    elseif observedHarvest && state.crop_mature_flag
        state.crop_phase = 3;
        state.crop_a_leaf = 0;
        state.crop_a_root = 0;
        state.crop_a_stem = 0;
        state.crop_a_repr = 0;
    elseif hui >= p.grnfill
        state.crop_phase = 3;
        state.crop_grainfill_event = phaseBefore < 3;
        [state.crop_a_leaf, state.crop_a_root, state.crop_a_stem, ...
            state.crop_a_repr] = grainFillAllocation(state, p, hui);
    elseif hui >= p.lfemerg
        state.crop_phase = 2;
        state.crop_emergence_event = phaseBefore ~= 2;
        state.onset_flag = state.crop_emergence_event;
        if state.crop_vf ~= 1 && hui < 0.8 * p.grnfill
            [state.crop_cumvd, state.crop_vf] = vernalization( ...
                state.crop_cumvd, forcing.T_avg);
        end
        [state.crop_a_leaf, state.crop_a_root, state.crop_a_stem, ...
            state.crop_a_repr] = vegetativeAllocation(p, hui);
    else
        state.crop_phase = 1;
        state.crop_a_leaf = 1e-5;
        state.crop_a_root = 0;
        state.crop_a_stem = 0;
        state.crop_a_repr = 0;
    end
end

state.dormant_flag = ~state.crop_live;
state.phen = double(state.crop_live && state.crop_phase >= 2 && state.crop_phase <= 3);
state.aphen = state.aphen + state.phen;
diag = struct('phen', state.phen, 'crop_phase', state.crop_phase, ...
    'hui', state.crop_hui, 'gddplant', state.crop_gddplant, ...
    'gddmaturity', state.crop_gddmaturity, 'vf', state.crop_vf, ...
    'planting_event', state.crop_planting_event, ...
    'emergence_event', state.crop_emergence_event, ...
    'harvest_event', state.crop_harvest_event, ...
    'grainfill_event', state.crop_grainfill_event, ...
    'maturity_event', state.crop_maturity_event, ...
    'mature', state.crop_mature_flag);
end

function [aleaf, aroot, astem, arepr] = vegetativeAllocation(p, hui)
% CoLM MOD_BGC_Veg_NutrientCompetition.F90:398-404.
aroot = p.arooti - (p.arooti - p.arootf) * hui;
fleaf = p.fleafi * (exp(-p.bfact) - exp(-p.bfact * hui / p.grnfill)) / ...
    (exp(-p.bfact) - 1);
arepr = 0;
aleaf = max(1e-5, (1 - aroot) * fleaf);
astem = 1 - arepr - aleaf - aroot;
end

function [aleaf, aroot, astem, arepr] = grainFillAllocation(state, p, hui)
% CoLM MOD_BGC_Veg_NutrientCompetition.F90:420-454.
aroot = p.arooti - (p.arooti - p.arootf) * min(1, hui);
decline = max(0, (1 - hui) / (1 - p.grnfill));
astem = max(p.astemf, state.crop_a_stem * decline^p.allconss);
aleaf = 1e-5;
arepr = 1 - aroot - astem - aleaf;
arepr = arepr * state.crop_vf;
aroot = 1 - aleaf - astem - arepr;
end

function [cumvd, vf] = vernalization(cumvd, tc)
% CoLM MOD_BGC_Veg_CNPhenology.F90:1599-1610, daily dt/24 = 1.
vtmin = -1.3;
vtopt = 4.9;
vtmax = 15.7;
alpha = log(2) / log((vtmax - vtmin) / (vtopt - vtmin));
if tc >= vtmin && tc <= vtmax
    cumvd = cumvd + (2 * (tc - vtmin)^alpha * (vtopt - vtmin)^alpha - ...
        (tc - vtmin)^(2 * alpha)) / (vtopt - vtmin)^(2 * alpha);
end
vf = cumvd^5 / (22.5^5 + cumvd^5);
end
