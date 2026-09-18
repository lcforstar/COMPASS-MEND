function [model, candidates] = stageCanopyInputs(model, kt, dt_outer_s)
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
    model.built.workspace.scope_step = init.prepareScopeStep(model, kt);
    lai_state_next = advanceHeatAccumulationLai( ...
        model, kt, dt_outer_s);
    model = applyCanopyLAI(model, kt, lai_state_next);
    [model, air_temperature] = applyGrowthTemperature(model, kt);
    model = applyDynamicCabLeafOptics(model);
    candidates = struct( ...
        'lai_state', lai_state_next, ...
        'leaf_resp_Ta', air_temperature);
end

function lai_state_next = advanceHeatAccumulationLai(model, kt, dt_outer_s)
% Advance the EPIC heat-unit LAI candidate at the first step of a new day.
    lai_state_next = model.built.Y.vegetation.lai_state;
    built = model.built;
    canopy_config = built.parameters.veg.canopy;
    if ~strcmpi(readField(canopy_config, 'LAI_mode', 'constant'), 'model')
        return;
    end
    time_control = built.control.time;
    steps_per_day = max(1, round(86400 / time_control.Delt_t));
    if kt <= 1 || mod(kt - 1, steps_per_day) ~= 0
        return;
    end
    timestamp = time_control.start_time + seconds( ...
        step.readAcceptedTime(model.runtime) + dt_outer_s);
    day_index = (year(timestamp) - year(time_control.start_time)) * 365 ...
        + day(timestamp, 'dayofyear');
    phenology = built.parameters.veg.phenology;
    driver_state = built.vegetation_state.phen_driver;
    parameters = struct( ...
        'DLAI', readField(phenology, 'DLAI', 0), ...
        'yr_cur', readField(phenology, 'yr_cur', 1), ...
        'yr_full', readField(phenology, 'yr_full', 1), ...
        'LAI_max', readField(canopy_config, 'LAI_max', 6), ...
        'T_BASE', readField(phenology, 'T_BASE', 5), ...
        'T_OPT', readField(phenology, 'T_OPT', 24), ...
        'ALAI_MIN', readField(phenology, 'ALAI_MIN', 0.1), ...
        'GMHU', readField(phenology, 'GMHU', 0), ...
        'is_deciduous', ~vegetation.getPFTValue( ...
            built.parameters.veg.structure.PFT_type, 'evergreen'), ...
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
