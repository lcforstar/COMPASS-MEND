function scopeStep = guardSurfaceTemperature(scopeStep, prevTss, kt)
%GUARDSURFACETEMPERATURE  Sanity-check the ebal surface temperature; fall back on reject.
%
%   Accept the ebal surface temperature th.Tsave when it is real, finite, in [-60,90] degC and
%   within 60 degC of Ta. On REJECT, fall back to the previous accepted Tss
%   (prevTss; or Ta when prevTss is NaN, i.e. kt=1 / no accepted step yet), zero the latent-heat
%   fluxes, and reset thermal/profiles/Soil to Tss/Ta. Pure: no model / owner access -- the caller
%   passes the fallback source in and stores the guarded Tss out (Julia-friendly (state, prev, t)).
%
%   Inputs:
%     scopeStep [in]  struct  SCOPE step bundle: .thermal (Tsave/Ts [degC][/Tcu/Tch/Tcave]),
%                             .fluxes (Actot/lEstot/lEctot [W m-2]), .meteo.Ta [degC], .profiles, .Soil
%     prevTss   [in]  [degC]  previous ACCEPTED surface temperature (fallback source); NaN if none yet
%     kt        [in]  [-]     outer step index (diagnostic print only)
%
%   Returns:
%     scopeStep [out] struct  guarded bundle: .thermal.Tss set; on reject also lEstot/lEctot zeroed
%                             and thermal/profiles/Soil reset to Tss/Ta
%
%   Method:
%     ok = real & finite & in [-60,90] degC & |Tsave - Ta| <= 60 degC. The 90 degC ceiling retains
%     a buffer above the observed 80.8 degC global land-skin record and the physically modelled
%     >85 degC hotspot range, while remaining far below the broad 400 K WRF input sanity ceiling.
%     On !ok, Tss = isfinite(prevTss)
%     ? prevTss : Ta, neutralize latent heat, and overwrite the step's canopy/soil temperature scratch
%     so a rejected ebal cannot poison the downstream snow/EVAP/soil-thermal top BC.
%
%   Reference:
%     COMPASS.m:1117-1146 (legacy SCOPE-GUARD, v30 thresholds); extracted verbatim from the
%     runOneStep.guardSurfaceTemperature local (state-ownership migration increment 2).
%
%   Editor: Jizixing Zhu; Date: 2026-07-05.
%   Change: 2026-07-14 Jizixing Zhu -- raise the hot-surface gates to 90
%           degC absolute and 60 degC above air using global LST evidence.
%   Change: 2026-07-20 Jizixing Zhu -- zero the sunlit/shaded wet/dry diagnostic
%           partition together with rejected canopy latent heat.
    sc = scopeStep;
    th = sc.thermal;
    fl = sc.fluxes;
    Ta = sc.meteo.Ta;

    ok = isreal(fl.Actot) && isreal(th.Tsave) && isreal(fl.lEstot) && isreal(fl.lEctot) && ...
         isfinite(th.Tsave) && all(isfinite(th.Ts(:)));
    if ok
        T_hi = max([th.Tsave; th.Ts(:)]);
        T_lo = min([th.Tsave; th.Ts(:)]);
        ok = T_lo >= -60 && T_hi <= 90 && abs(th.Tsave - Ta) <= 60;
    end

    if ok
        Tss = th.Tsave;
    else
        if isfinite(prevTss), Tss = prevTss; else, Tss = Ta; end   % kt=1 / no previous accepted Tss
        fl.lEstot = 0;
        fl.lEctot = 0;
        if isfield(fl, 'lE_wet_canopy'), fl.lE_wet_canopy = 0; end
        if isfield(fl, 'lE_dry_canopy'), fl.lE_dry_canopy = 0; end
        for name = {'lE_wet_sun','lE_wet_sha','lE_dry_sun','lE_dry_sha'}
            if isfield(fl, name{1}), fl.(name{1}) = 0; end
        end
        fprintf(['[SCOPE-GUARD] kt=%d rejected ebal surface temperature. Ta=%.3f Tss_fallback=%.3f ' ...
            'thermal.Tsave=%.3f thermal.Ts=[%.3f %.3f]\n'], kt, Ta, Tss, th.Tsave, th.Ts(1), th.Ts(2));
        th.Ts    = Tss * ones(2, 1);
        th.Tsave = Tss;
        if isfield(th, 'Tcu'),   th.Tcu(:)   = Ta; end
        if isfield(th, 'Tch'),   th.Tch(:)   = Ta; end
        if isfield(th, 'Tcave'), th.Tcave    = Ta; end
        prof = sc.profiles;
        for f = {'Ta_layer','Tleaf_layer','Tch','Tcu1d','Tc1d'}
            if isfield(prof, f{1}), prof.(f{1})(:) = Ta; end
        end
        sc.profiles = prof;
        if isfield(sc, 'Soil')
            if isfield(sc.Soil, 'Ts'), sc.Soil.Ts = Tss * ones(2, 1); end
            if isfield(sc.Soil, 'Tsold') && isnumeric(sc.Soil.Tsold) && ~isempty(sc.Soil.Tsold)
                sc.Soil.Tsold(1, :) = Tss;
            end
        end
    end

    th.Tss     = Tss;       % [degC] guarded surface temperature (caller commits accept-only to its owner)
    sc.thermal = th;
    sc.fluxes  = fl;
    scopeStep  = sc;
end
