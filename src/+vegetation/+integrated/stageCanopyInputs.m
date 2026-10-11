function [model, candidates] = stageCanopyInputs(model, kt, dt_outer_s)
% Change: 2026-09-30 Wenjie Liu -- selectively restore COMPASS canopy
% behavior; retain winter-wheat phenology and CoLM C/N allocation.
%STAGECANOPYINPUTS  Finalize LAI, growth temperature, and leaf optics.
%
%   Inputs:
%     model      [inout] Model shell before the canopy producers run.
%     kt         [-] One-based outer-step index.
%     dt_outer_s [s] Actual outer-step duration.
%
%   Returns:
%     model      [out] Model with the current SCOPE input bundle staged.
%     candidates [struct] Accept-only LAI and leaf-temperature-history state.
%
%   The stage preserves the load-bearing order prepareScopeStep -> daily LAI
%   advance -> final canopy LAI/height -> growth temperature -> dynamic leaf
%   optics. prepareScopeStep retains Soil_bare plus an accepted-snow seed;
%   stageCanopyEnergy replaces that seed with the current precipitation-only
%   snow prediction before radiation transfer.
%   Editor: Jizixing Zhu; Date: 2026-07-16.
%   Change: 2026-07-17 Jizixing Zhu -- document ownership of the bare-surface
%     template and the accepted-snow seed used by same-step snow prediction.
% Change: 2026-09-29 Wenjie Liu -- optional BGC replay context supplies the
% accepted preceding-cycle air-temperature tail; ordinary runs unchanged.
% Change: 2026-09-30 Wenjie Liu -- freeze daily/crop events before canopy
% exchange and expose that candidate to all same-step downstream consumers.
    model.built.workspace.scope_step = init.prepareScopeStep(model, kt);
    phenology_next = vegetation.integrated.phenology.stageCandidate(model,kt,dt_outer_s);
    model.built.workspace.scope_step.phenology_candidate = phenology_next;
    lai_state_next = advanceHeatAccumulationLai( ...
        model, kt, dt_outer_s);
    model = applyCanopyLAI(model, kt, lai_state_next);
    [model, air_temperature] = applyGrowthTemperature(model, kt);
    model = applyDynamicCabLeafOptics(model);
    candidates = struct( ...
        'phenology', phenology_next, ...
        'lai_state', lai_state_next, ...
        'leaf_resp_Ta', air_temperature);
end

function lai_state_next = advanceHeatAccumulationLai(model, kt, ~)
% Advance the EPIC heat-unit LAI candidate at the first step of a new day.
    lai_state_next = model.built.Y.vegetation.lai_state;
    built = model.built;
    canopy_config = built.parameters.veg.canopy;
    if ~strcmpi(readField(canopy_config, 'LAI_mode', 'constant'), 'model')
        return;
    end
    engine = char(string(readField(canopy_config,'LAI_engine','epic_phen_mix')));
    if strcmpi(engine,'leafc_sla_phen') && readField(built.control.switches,'lPlantCN',0)==1 ...
            && readField(built.control.switches,'lNefc',0)==1
        return; % leaf-C engine does not advance a second empirical LAI state.
    end
    time_control = built.control.time;
    steps_per_day = max(1, round(86400 / time_control.Delt_t));
    timestamp = time_control.start_time + seconds(step.readAcceptedTime(model.runtime));
    phenology = built.parameters.veg.phenology;
    driver_state = built.vegetation_state.phen_driver;
    forcing_start = readField(driver_state,'forcing_start_time',time_control.start_time);
    day_index = floor(days(dateshift(timestamp,'start','day') ...
        - dateshift(forcing_start,'start','day'))) + 1;
    is_crop = vegetation.getPFTValue(built.parameters.veg.structure.PFT_type,'crop_phenology',0)==1;
    if is_crop
        if readField(lai_state_next,'last_advance_day_index',0)==day_index,return;end
    elseif kt<=1 || mod(kt-1,steps_per_day)~=0
        return;
    end
    parameters = struct( ...
        'DLAI', readField(phenology, 'DLAI', 0), ...
        'yr_cur', readField(phenology, 'yr_cur', 1), ...
        'yr_full', readField(phenology, 'yr_full', 1), ...
        'LAI_max', readField(canopy_config, 'LAI_max', 6), ...
        'T_BASE', readField(phenology, 'T_BASE', 5), ...
        'T_OPT', readField(phenology, 'T_OPT', 24), ...
        'ALAI_MIN', readField(phenology, 'ALAI_MIN', 0.1), ...
        'GMHU', readField(phenology, 'GMHU', 0), ...
        'is_deciduous', ~is_crop && ~vegetation.getPFTValue( ...
            built.parameters.veg.structure.PFT_type, 'evergreen'), ...
        'is_crop', is_crop, ...
        'calendar_doy', day(timestamp,'dayofyear'), ...
        'calendar_year', year(timestamp), ...
        'crop_planting_event', model.built.workspace.scope_step.phenology_candidate.state.crop_planting_event, ...
        'crop_live', model.built.workspace.scope_step.phenology_candidate.state.crop_live, ...
        'crop_phase', model.built.workspace.scope_step.phenology_candidate.state.crop_phase, ...
        'sowing_day', readField(phenology, 'sowing_day', 1), ...
        'beta_decay', readField(phenology, 'beta_decay', 1));
    driver = struct( ...
        'phen_CHU_model', driver_state.phen_CHU, ...
        'phen_HUI_model', driver_state.phen_HUI, ...
        'phen_dHUF_model', driver_state.phen_dHUF);
    max_temperature = readSeries( ...
        driver_state, 'phen_Tmax_daily', day_index, 15);
    min_temperature = readSeries( ...
        driver_state, 'phen_Tmin_daily', day_index, 15);
    previous_stress = readField(lai_state_next, 'wt_day_min', 1.0);
    lai_state_next = vegetation.integrated.phenology.updateLaiHydraulic( ...
        lai_state_next, day_index, previous_stress, max_temperature, ...
        min_temperature, parameters, driver);
    lai_state_next.wt_day_min = 1.0;
    if is_crop,lai_state_next.last_advance_day_index=day_index;end
end

function model = applyCanopyLAI(model, kt, lai_state)
% Apply the configured canopy height and LAI consumer to this step bundle.
    built = model.built;
    vegParameters = built.parameters.veg;
    canopy_config = vegParameters.canopy;
    canopy = built.workspace.scope_step.canopy;

    if strcmpi(readField(canopy_config, 'hc_mode', 'constant'), 'constant')
        canopy.hc = readField( ...
            canopy_config, 'canopy_height', canopy.hc);
    else
        canopy.hc = readSeries(lai_state, 'hc', kt, canopy.hc);
    end
    canopy.hot = canopy.leafwidth / max(canopy.hc, eps);

    switch lower(readField(canopy_config, 'LAI_mode', 'constant'))
        case 'constant'
            canopy.LAI = readField( ...
                canopy_config, 'LAI_constant', canopy.LAI);
        case 'model'
            switches = built.control.switches;
            if readField(switches, 'lPlantCN', 0) ~= 1 ...
                    || readField(switches, 'lNefc', 0) ~= 1
                canopy.LAI = boundedHeatAccumulationLai( ...
                    lai_state, vegParameters);
            else
                engine = lower(char(string(readField( ...
                    canopy_config, 'LAI_engine', 'epic_phen_mix'))));
                switch engine
                    case 'leafc_sla_phen'
                        phenology = model.built.Y.vegetation.phenology.state.phen;
                        leaf_carbon = built.Y.plant_cn.C.leaf;
                        canopy.LAI = ...
                            vegetation.integrated.phenology.computeLAIfromPhen( ...
                            laiEngineConfig(vegParameters), ...
                            vegParameters.structure.PFT_type, phenology, ...
                            readSeries(lai_state, 'LAI', kt, canopy.LAI), ...
                            leaf_carbon);
                    case 'epic_phen_mix'
                        canopy.LAI = boundedHeatAccumulationLai( ...
                            lai_state, vegParameters);
                    otherwise
                        error('runOneStep:badLaiEngine', ...
                            ['Unknown LAI_engine "%s". Expected ' ...
                             'epic_phen_mix or leafc_sla_phen.'], engine);
                end
            end
        otherwise
            canopy.LAI = readSeries(lai_state, 'LAI', kt, canopy.LAI);
    end

    if strcmpi(readField(canopy_config,'LAI_mode','constant'),'model') ...
            && vegetation.getPFTValue(vegParameters.structure.PFT_type,'winter_wheat_cn',0)==1
        ps = built.workspace.scope_step.phenology_candidate.state;
        if ~ps.crop_live || ps.crop_phase<2
            canopy.LAI = 0; % explicit pre-emergence/harvest green-canopy state.
        end
    end
    canopy.green_LAI = canopy.LAI; % this bundle's LAI is photosynthetic green area.
    if isfield(canopy, 'layer_weight') && ~isempty(canopy.layer_weight)
        canopy.LAI_layer = canopy.layer_weight * canopy.LAI;
    end
    model.built.workspace.scope_step.canopy = canopy;
end

function lai = boundedHeatAccumulationLai(lai_state, vegParameters)
% Read the heat-unit LAI state and enforce configured physical bounds.
    lai_min = readField(vegParameters.phenology, 'ALAI_MIN', 0);
    lai_max = readField(vegParameters.canopy, 'LAI_max', Inf);
    lai_raw = lai_state.lai_model_current;
    assert(isscalar(lai_raw) && isfinite(lai_raw), ...
        'runOneStep:badHeatAccumulationLai', ...
        'lai_model_current must be a finite scalar, got size %s.', ...
        mat2str(size(lai_raw)));
    lai = min(lai_max, max(lai_min, lai_raw));
end

function config = laiEngineConfig(vegParameters)
% Build the transitional cfg.Veg view used by computeLAIfromPhen.
    config = struct('Veg', struct( ...
        'LAI_engine', readField( ...
            vegParameters.canopy, 'LAI_engine', 'epic_phen_mix'), ...
        'ALAI_MIN', readField(vegParameters.phenology, 'ALAI_MIN', 0), ...
        'LAI_max', readField(vegParameters.canopy, 'LAI_max', 6)));
end

function [model, air_temperature] = applyGrowthTemperature(model, kt)
% Stage the ten-day growth-temperature mean and current accepted candidate.
    history = model.built.Y.plant_acclimation.leaf_resp.Ta_history;
    window = ...
        model.built.p.plant_acclimation.leaf_resp_Tgrowth_window;
    air_temperature = model.built.workspace.scope_step.meteo.Ta;
    history = history(:);
    if numel(history) < kt - 1
        history(end + 1:kt - 1, 1) = NaN;
    end
    first_index = max(1, kt - window + 1);
    values = [history(first_index:kt - 1); air_temperature];
    if kt<window&&isfield(model.runtime,'biogeochemistry_cycle') ...
            &&model.runtime.biogeochemistry_cycle.active
        tail=model.runtime.biogeochemistry_cycle.previous_air_temperature_tail_C(:);
        need=window-kt;
        values=[tail(max(1,numel(tail)-need+1):end);values];
    end
    model.built.workspace.scope_step.canopy.Tgrowth = ...
        mean(values, 'omitnan');
end

function model = applyDynamicCabLeafOptics(model)
% Apply lagged LUNA chlorophyll to the current optical-band leaf properties.
    built = model.built;
    if ~logical(readField( ...
            built.parameters.veg.photosynthesis, 'cab_from_luna', 0))
        return;
    end
    luna_state = built.Y.plant_acclimation.luna_state;
    if ~built.p.plant_acclimation.luna_enabled || luna_state.first_day ...
            || ~isfinite(luna_state.last_Cab_nlc)
        return;
    end
    dynamic_cab = luna_state.last_Cab_nlc;
    if ~(isfinite(dynamic_cab) && dynamic_cab > 0 && dynamic_cab < 120)
        return;
    end

    scope_static = built.p.scope_static;
    leaf_bio = built.workspace.scope_step.leafbio;
    leaf_optical = built.workspace.scope_step.leafopt;
    if leaf_bio.Cab > 0
        leaf_bio.Cca = leaf_bio.Cca * (dynamic_cab / leaf_bio.Cab);
    end
    leaf_bio.Cab = dynamic_cab;

    if readField(scope_static.options, 'calc_PSI', 0)
        fluspect = @surface.radiation.fluspect_B_CX;
    else
        fluspect = @surface.radiation.fluspect_B_CX_PSI_PSII_combined;
    end
    if ~isfield(leaf_bio, 'V2Z')
        leaf_bio.V2Z = 0;
    end
    dynamic_optical = fluspect( ...
        scope_static.spectral, leaf_bio, scope_static.library.optipar);
    optical_index = scope_static.spectral.IwlP;
    leaf_optical.refl(optical_index) = dynamic_optical.refl;
    leaf_optical.tran(optical_index) = dynamic_optical.tran;

    built.workspace.scope_step.leafbio = leaf_bio;
    built.workspace.scope_step.leafopt = leaf_optical;
    model.built = built;
end

function value = readField(container, name, fallback)
% Read a nonempty struct field or return a local fallback.
    value = fallback;
    if isfield(container, name) && ~isempty(container.(name))
        value = container.(name);
    end
end

function value = readSeries(container, name, index, fallback)
% Read one indexed series value or return a local fallback.
    value = fallback;
    if isfield(container, name) && numel(container.(name)) >= index
        value = container.(name)(index);
    end
end
