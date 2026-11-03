function laiState = updateLaiHydraulic(laiState, dayIdx, wt, Tmax_d, Tmin_d, params, driver)
%UPDATELAIHYDRAULIC  Pure daily heat-unit LAI update with hydraulic stress.
%
%   laiState = vegetation.integrated.phenology.updateLaiHydraulic(laiState, dayIdx, wt, ...
%       Tmax_d, Tmin_d, params, driver) advances the per-day EPIC/frPHU LAI state one
%   day: a dormancy -> pre-germination -> growth -> senescence state machine whose
%   growth is gated by temperature stress f_T and a hydraulic stress factor wt in [0,1]
%   (replacing the Feddes f_theta of the original scheme).
%
%   This is a PURE, de-handled port of legacy VegetationData.update_lai_with_hydraulic
%   (that handle method is left UNCHANGED so the legacy COMPASS.m run is unaffected).
%   Value in / value out: it mutates no object -- the caller commits laiState on accept.
%
%   Inputs:
%     laiState  [in]  struct, mutable per-day state: .lai_day_dormant [bool],
%                     .lai_day_LAI_D / .lai_day_HUI_D [-], .lai_day_yr [-],
%                     .lai_day_prev [m2/m2], .lai_model_current [m2/m2]
%     dayIdx    [in]  [-]        absolute model day from simulation start (1-based)
%     wt        [in]  [0,1]      daily min hydraulic stress min(wt_sun,wt_sha)
%     Tmax_d    [in]  [degC]     daily max temperature
%     Tmin_d    [in]  [degC]     daily min temperature
%     params    [in]  struct     static EPIC/frPHU params: DLAI, yr_cur, yr_full,
%                     LAI_max, T_BASE, T_OPT, ALAI_MIN, GMHU, is_deciduous,
%                     sowing_day, beta_decay
%     driver    [in]  struct     per-day heat-unit arrays: phen_CHU_model,
%                     phen_HUI_model, phen_dHUF_model (indexed by dayIdx)
%
%   Returns:
%     laiState  updated struct (.lai_day_* + .lai_model_current)
%
%   Reference: legacy VegetationData.m:748-831 (Method B); Williams 1989 EPIC LAI.
%   Editor: Jizixing Zhu; Date: 2026-07-02.
%   Change: 2026-09-30 Wenjie Liu -- annual crops use explicit planting and
%   harvest state, no tree-age multiplier, autumn leaf-off, or January reset.
    if dayIdx < 1 || dayIdx > numel(driver.phen_CHU_model)
        return;   % outside the precomputed phenology window -> no update
    end

    % Year boundary: reset dormancy + decay state every 365 days (VegetationData.m:766-772)
    isCrop = isfield(params,'is_crop') && params.is_crop;
    withinYearDay = mod(dayIdx - 1, 365) + 1;
    if isCrop && isfield(params,'calendar_doy'),withinYearDay=params.calendar_doy;end
    if isCrop && params.crop_planting_event
        laiState.lai_day_dormant = false;
        laiState.lai_day_LAI_D = 0;
        laiState.lai_day_HUI_D = params.DLAI;
        laiState.lai_day_prev = 0;
    elseif ~isCrop && withinYearDay == 1 && dayIdx > 1
        laiState.lai_day_dormant = false;
        laiState.lai_day_LAI_D   = 0;
        laiState.lai_day_HUI_D   = params.DLAI;
        laiState.lai_day_yr      = laiState.lai_day_yr + 1;
    end

    if isCrop
        XLAI = params.LAI_max;
    else
        yrActual = params.yr_cur + laiState.lai_day_yr;
        XLAI = min(yrActual / params.yr_full, 1) * params.LAI_max;
    end

    % Temperature stress f_T (VegetationData.m:778-784)
    Tav = (Tmax_d + Tmin_d) / 2;
    R_T = (Tav - params.T_BASE) / (params.T_OPT - params.T_BASE);
    if R_T >= 0 && R_T <= 2
        f_T = sin(pi/2 * R_T);
    else
        f_T = 0;
    end

    % Hydraulic stress replaces Feddes f_theta (VegetationData.m:787-788)
    wtSafe = max(0, min(1, wt));
    fEff   = min(f_T, wtSafe);

    CHU_d  = driver.phen_CHU_model(dayIdx);
    HUI_d  = driver.phen_HUI_model(dayIdx);
    dHUF_d = driver.phen_dHUF_model(dayIdx);

    % Deciduous dormancy lock (irreversible within-year) (VegetationData.m:795-799)
    if ~isCrop && params.is_deciduous && ~laiState.lai_day_dormant
        if withinYearDay > 300 && (Tmin_d < 0 || withinYearDay >= 320)
            laiState.lai_day_dormant = true;
        end
    end

    % LAI state machine: dormancy -> pre-germination -> growth -> senescence
    cropBare = isCrop && (~params.crop_live || params.crop_phase<2);
    if cropBare
        LAI_new = 0;
    elseif laiState.lai_day_dormant
        LAI_new = params.ALAI_MIN;
    elseif CHU_d <= params.GMHU
        LAI_new = params.ALAI_MIN;
    elseif HUI_d > params.DLAI
        if laiState.lai_day_LAI_D == 0
            laiState.lai_day_LAI_D = laiState.lai_day_prev;
            laiState.lai_day_HUI_D = HUI_d;
        end
        decay   = max(0, (1.0 - HUI_d) / (1.0 - laiState.lai_day_HUI_D)) ^ params.beta_decay;
        LAI_new = max(params.ALAI_MIN, laiState.lai_day_LAI_D * decay);
    else
        if withinYearDay == params.sowing_day
            LAI_new = laiState.lai_day_prev;
        else
            growthRes = max(0, 1 - exp(5 * (laiState.lai_day_prev - XLAI)));
            LAI_new   = laiState.lai_day_prev + dHUF_d * XLAI * growthRes * fEff;
        end
    end

    if ~cropBare,LAI_new = max(params.ALAI_MIN, LAI_new);end
    if ~isfinite(LAI_new), LAI_new = laiState.lai_day_prev; end

    laiState.lai_day_prev      = LAI_new;
    laiState.lai_model_current = LAI_new;
end
